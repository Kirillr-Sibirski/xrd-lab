#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from amcsd_cache import ensure_amcsd_library, select_amcsd_candidate_files
from cif_matcher import collect_reference_files, identify_likely_crystals


DEFAULT_WAVELENGTH = 1.5406  # Cu K-alpha, in Angstrom


@dataclass(frozen=True)
class Peak:
    two_theta_deg: float
    intensity: float
    prominence: float
    left_index: int
    right_index: int


@dataclass(frozen=True)
class ReflectionFamily:
    n_value: int
    families: tuple[str, ...]


@dataclass(frozen=True)
class LatticeFit:
    lattice_type: str
    first_n_value: int
    assigned_n_values: tuple[int, ...]
    a_values_angstrom: tuple[float, ...]
    a_mean_angstrom: float
    a_std_angstrom: float
    score: float


def parse_int_file(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text().splitlines()
    data_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        try:
            two_theta = float(parts[0])
            intensity = float(parts[1])
        except ValueError:
            continue
        data_lines.append((two_theta, intensity))

    if not data_lines:
        raise ValueError(f"No diffraction points found in {path}")

    data = np.asarray(data_lines, dtype=float)
    return data[:, 0], data[:, 1]


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return values.copy()
    if window % 2 == 0:
        window += 1
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(values, kernel, mode="same")


def detect_peaks(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    smooth_window: int = 9,
    prominence_window: int = 50,
    relative_threshold: float = 0.015,
    merge_tolerance_deg: float = 0.35,
) -> list[Peak]:
    smoothed = moving_average(intensity, smooth_window)
    max_signal = float(smoothed.max())
    candidates: list[Peak] = []

    for idx in range(1, len(smoothed) - 1):
        if smoothed[idx] <= smoothed[idx - 1] or smoothed[idx] < smoothed[idx + 1]:
            continue

        left = max(0, idx - prominence_window)
        right = min(len(smoothed), idx + prominence_window + 1)
        baseline = max(float(smoothed[left:idx].min()), float(smoothed[idx + 1 : right].min()))
        prominence = float(smoothed[idx] - baseline)
        if prominence < relative_threshold * max_signal:
            continue

        candidates.append(
            Peak(
                two_theta_deg=float(two_theta[idx]),
                intensity=float(intensity[idx]),
                prominence=prominence,
                left_index=idx,
                right_index=idx,
            )
        )

    if not candidates:
        return []

    candidates.sort(key=lambda peak: peak.two_theta_deg)
    merged: list[Peak] = []
    cluster: list[Peak] = [candidates[0]]

    def flush_cluster() -> None:
        if not cluster:
            return
        weights = np.asarray(
            [max(member.prominence, 1e-9) * max(member.intensity, 1e-9) for member in cluster],
            dtype=float,
        )
        positions = np.asarray([member.two_theta_deg for member in cluster], dtype=float)
        representative = max(cluster, key=lambda member: member.intensity)
        merged.append(
            Peak(
                two_theta_deg=float(np.average(positions, weights=weights)),
                intensity=representative.intensity,
                prominence=max(member.prominence for member in cluster),
                left_index=min(member.left_index for member in cluster),
                right_index=max(member.right_index for member in cluster),
            )
        )

    for peak in candidates[1:]:
        if peak.two_theta_deg - cluster[-1].two_theta_deg <= merge_tolerance_deg:
            cluster.append(peak)
            continue
        flush_cluster()
        cluster = [peak]
    flush_cluster()

    return merged


def cubic_reflection_families(max_index: int = 10) -> dict[str, list[ReflectionFamily]]:
    by_lattice: dict[str, dict[int, set[str]]] = {
        "primitive": {},
        "bcc": {},
        "fcc": {},
        "diamond": {},
    }

    for h in range(max_index + 1):
        for k in range(h + 1):
            for l in range(k + 1):
                if h == k == l == 0:
                    continue

                n_value = h * h + k * k + l * l
                family = f"({h}{k}{l})"
                parity = (h % 2, k % 2, l % 2)
                same_parity = parity[0] == parity[1] == parity[2]
                is_bcc = (h + k + l) % 2 == 0
                is_fcc = same_parity
                is_diamond = is_fcc and ((h % 2 == 1) or ((h + k + l) % 4 == 0))

                by_lattice["primitive"].setdefault(n_value, set()).add(family)
                if is_bcc:
                    by_lattice["bcc"].setdefault(n_value, set()).add(family)
                if is_fcc:
                    by_lattice["fcc"].setdefault(n_value, set()).add(family)
                if is_diamond:
                    by_lattice["diamond"].setdefault(n_value, set()).add(family)

    families: dict[str, list[ReflectionFamily]] = {}
    for lattice_type, mapping in by_lattice.items():
        families[lattice_type] = [
            ReflectionFamily(n_value=n_value, families=tuple(sorted(mapping[n_value])))
            for n_value in sorted(mapping)
        ]
    return families


def fit_cubic_lattice(
    peaks: list[Peak],
    wavelength_angstrom: float,
    max_candidates: int = 12,
) -> LatticeFit:
    if len(peaks) < 3:
        raise ValueError("Need at least three peaks for a cubic fit")

    theta = np.radians(np.asarray([peak.two_theta_deg for peak in peaks], dtype=float) / 2.0)
    q_values = np.sin(theta) ** 2
    ratios = q_values / q_values[0]

    reflection_db = cubic_reflection_families()
    best_fit: LatticeFit | None = None
    lattice_penalty = {
        "primitive": 0.0,
        "bcc": 0.0002,
        "fcc": 0.0004,
        "diamond": 0.0006,
    }

    for lattice_type, families in reflection_db.items():
        n_values = [entry.n_value for entry in families]
        for start_idx, first_n_value in enumerate(n_values[:max_candidates]):
            assigned = [first_n_value]
            cursor = start_idx
            ratio_error_sum = 0.0

            for ratio in ratios[1:]:
                best_j = None
                best_error = None
                for j in range(cursor + 1, len(n_values)):
                    expected = n_values[j] / first_n_value
                    error = abs(expected - ratio)
                    if best_error is None or error < best_error:
                        best_error = error
                        best_j = j
                    if expected > ratio and best_error is not None and error > best_error:
                        break

                if best_j is None or best_error is None:
                    assigned = []
                    break

                assigned.append(n_values[best_j])
                cursor = best_j
                ratio_error_sum += best_error

            if len(assigned) != len(peaks):
                continue

            a_values = tuple(
                wavelength_angstrom * math.sqrt(n_value) / (2.0 * math.sin(theta_i))
                for n_value, theta_i in zip(assigned, theta)
            )
            a_mean = float(np.mean(a_values))
            a_std = float(np.std(a_values, ddof=0))
            score = (
                ratio_error_sum / max(len(peaks) - 1, 1)
                + (a_std / a_mean)
                + 0.0002 * start_idx
                + lattice_penalty[lattice_type]
            )

            candidate = LatticeFit(
                lattice_type=lattice_type,
                first_n_value=first_n_value,
                assigned_n_values=tuple(assigned),
                a_values_angstrom=a_values,
                a_mean_angstrom=a_mean,
                a_std_angstrom=a_std,
                score=score,
            )

            if best_fit is None or candidate.score < best_fit.score:
                best_fit = candidate

    if best_fit is None:
        raise ValueError("Unable to fit a cubic lattice to the detected peaks")

    return best_fit


def gcd_of_values(values: Iterable[int]) -> int:
    iterator = iter(values)
    try:
        result = next(iterator)
    except StopIteration:
        return 0
    for value in iterator:
        result = math.gcd(result, value)
    return result


def build_report(
    source_path: Path,
    peaks: list[Peak],
    fit: LatticeFit,
    wavelength_angstrom: float,
    top_n: int,
) -> dict[str, object]:
    theta = np.radians(np.asarray([peak.two_theta_deg for peak in peaks], dtype=float) / 2.0)
    d_spacings = wavelength_angstrom / (2.0 * np.sin(theta))
    q_values = np.sin(theta) ** 2
    ratios = q_values / q_values[0]
    common_factor = gcd_of_values(fit.assigned_n_values)
    reflection_lookup = {
        entry.n_value: entry.families for entry in cubic_reflection_families()[fit.lattice_type]
    }

    peak_rows = []
    for index, (peak, d_spacing, ratio, n_value, a_value) in enumerate(
        zip(peaks[:top_n], d_spacings[:top_n], ratios[:top_n], fit.assigned_n_values[:top_n], fit.a_values_angstrom[:top_n]),
        start=1,
    ):
        peak_rows.append(
            {
                "peak_index": index,
                "two_theta_deg": round(peak.two_theta_deg, 4),
                "intensity": round(peak.intensity, 4),
                "d_spacing_angstrom": round(float(d_spacing), 6),
                "sin2_ratio": round(float(ratio), 6),
                "assigned_n": n_value,
                "possible_hkl_families": list(reflection_lookup.get(n_value, ())),
                "a_from_peak_angstrom": round(a_value, 6),
            }
        )

    return {
        "source_file": str(source_path),
        "wavelength_angstrom": wavelength_angstrom,
        "candidate_search": {
            "searched_reference_files": [],
            "reference_file_count": 0,
            "note": "No CIF reference library scanned yet.",
        },
        "slide_workflow_summary": {
            "assumption": "cubic unit cell",
            "best_fit_lattice_type": fit.lattice_type,
            "first_peak_n_value": fit.first_n_value,
            "common_factor_in_assigned_n_values": common_factor,
            "estimated_a_angstrom": round(fit.a_mean_angstrom, 6),
            "a_std_angstrom": round(fit.a_std_angstrom, 6),
            "lowest_a_angstrom": round(min(fit.a_values_angstrom), 6),
            "highest_a_angstrom": round(max(fit.a_values_angstrom), 6),
        },
        "detected_peaks": peak_rows,
    }


def add_candidate_matches(
    report: dict[str, object],
    reference_files: list[Path],
    matches: list[object],
    used_amcsd_cache: bool = False,
) -> None:
    if used_amcsd_cache:
        note = (
            "Probabilities are relative to the scanned AMCSD candidate subset filtered from the cached CIF archive."
        )
    elif reference_files:
        note = "Probabilities are relative to the currently scanned CIF library, not absolute confidence."
    else:
        note = "No CIF files found. Add references with --references or let the tool fetch AMCSD automatically."
    report["candidate_search"] = {
        "searched_reference_files": [str(path) for path in reference_files],
        "reference_file_count": len(reference_files),
        "note": note,
    }
    report["likely_crystals"] = [
        {
            "rank": index,
            "name": match.name,
            "formula": match.formula,
            "space_group": match.space_group,
            "lattice_a_angstrom": round(match.lattice_a_angstrom, 6),
            "probability_score": round(match.probability, 6),
            "match_score": round(match.score, 6),
            "matched_peak_count": match.matched_peak_count,
            "reference_peak_count": match.total_reference_peaks,
            "reference_file": str(match.source_path),
        }
        for index, match in enumerate(matches, start=1)
    ]


def format_text_report(report: dict[str, object]) -> str:
    summary = report["slide_workflow_summary"]
    peaks = report["detected_peaks"]
    candidate_search = report.get("candidate_search", {})
    likely_crystals = report.get("likely_crystals", [])

    lines = [
        f"File: {report['source_file']}",
        f"Wavelength: {report['wavelength_angstrom']:.4f} A",
        "",
        "Cubic analysis summary",
        f"  Best-fit lattice type: {summary['best_fit_lattice_type']}",
        f"  First assigned h^2+k^2+l^2 value: {summary['first_peak_n_value']}",
        f"  Common factor in assigned h^2+k^2+l^2 values: {summary['common_factor_in_assigned_n_values']}",
        f"  Estimated lattice parameter a: {summary['estimated_a_angstrom']:.6f} A",
        f"  Spread across fitted peaks: {summary['a_std_angstrom']:.6f} A",
        f"  Lowest per-peak a: {summary['lowest_a_angstrom']:.6f} A",
        f"  Highest per-peak a: {summary['highest_a_angstrom']:.6f} A",
        "",
        "Peak table",
        "idx  2theta(deg)  intensity   d(A)      sin^2 ratio  N=h^2+k^2+l^2  families        a_from_peak(A)",
    ]

    for row in peaks:
        families = ",".join(row["possible_hkl_families"]) or "-"
        lines.append(
            f"{row['peak_index']:>3}  "
            f"{row['two_theta_deg']:>10.4f}  "
            f"{row['intensity']:>9.4f}  "
            f"{row['d_spacing_angstrom']:>8.4f}  "
            f"{row['sin2_ratio']:>11.4f}  "
            f"{row['assigned_n']:>14}  "
            f"{families:<14}  "
            f"{row['a_from_peak_angstrom']:>13.6f}"
        )

    lines.extend(
        [
            "",
            "Likely crystals",
            f"  Reference CIF files scanned: {candidate_search.get('reference_file_count', 0)}",
            f"  Note: {candidate_search.get('note', '')}",
        ]
    )
    if not likely_crystals:
        lines.append("  No candidate matches yet.")
        return "\n".join(lines)

    lines.append("rank  probability  score    a(A)      name / formula / space group")
    for candidate in likely_crystals:
        descriptor = " / ".join(
            part
            for part in [candidate["name"], candidate["formula"], candidate["space_group"]]
            if part
        )
        lines.append(
            f"{candidate['rank']:>4}  "
            f"{candidate['probability_score']:>10.4f}  "
            f"{candidate['match_score']:>7.4f}  "
            f"{candidate['lattice_a_angstrom']:>8.4f}  "
            f"{descriptor}"
        )

    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze cubic powder XRD .int files using the workflow from the Intro XRD slides."
    )
    parser.add_argument("input_file", type=Path, help="Path to the .int diffraction file")
    parser.add_argument(
        "--wavelength",
        type=float,
        default=DEFAULT_WAVELENGTH,
        help="X-ray wavelength in Angstrom (default: Cu K-alpha 1.5406 A)",
    )
    parser.add_argument("--smooth-window", type=int, default=9, help="Moving-average window for noise reduction")
    parser.add_argument(
        "--relative-threshold",
        type=float,
        default=0.015,
        help="Minimum prominence as a fraction of the strongest smoothed peak",
    )
    parser.add_argument(
        "--prominence-window",
        type=int,
        default=50,
        help="Number of points to each side used when estimating local prominence",
    )
    parser.add_argument(
        "--merge-tolerance",
        type=float,
        default=0.35,
        help="Merge nearby maxima within this 2theta distance, useful for K-alpha splitting",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=12,
        help="Maximum number of peaks to include in the cubic fit and the report",
    )
    parser.add_argument(
        "--references",
        type=Path,
        nargs="*",
        help="Optional CIF files or directories used to identify likely crystals",
    )
    parser.add_argument(
        "--top-candidates",
        type=int,
        default=5,
        help="How many likely crystal matches to show when CIF references are available",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Optional cache directory for the downloaded AMCSD archive and its metadata index",
    )
    parser.add_argument(
        "--max-elements",
        type=int,
        default=3,
        help="Maximum number of unique chemical elements allowed in candidate CIFs",
    )
    parser.add_argument(
        "--amcsd-limit",
        type=int,
        default=400,
        help="How many AMCSD candidates to scan after metadata filtering; use 0 for the full filtered set",
    )
    parser.add_argument(
        "--no-fetch-amcsd",
        action="store_true",
        help="Do not auto-download and cache the AMCSD CIF archive when no local references are provided",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON instead of a text report")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    two_theta, intensity = parse_int_file(args.input_file)
    peaks = detect_peaks(
        two_theta,
        intensity,
        smooth_window=args.smooth_window,
        prominence_window=args.prominence_window,
        relative_threshold=args.relative_threshold,
        merge_tolerance_deg=args.merge_tolerance,
    )

    if not peaks:
        raise SystemExit("No peaks detected. Lower --relative-threshold or inspect the input data.")

    selected_peaks = peaks[: args.top_n]
    fit = fit_cubic_lattice(selected_peaks, wavelength_angstrom=args.wavelength)
    report = build_report(args.input_file, selected_peaks, fit, args.wavelength, args.top_n)
    reference_files = collect_reference_files(args.input_file, args.references)
    used_amcsd_cache = False
    if args.references is None and not args.no_fetch_amcsd:
        cif_dir, index_path = ensure_amcsd_library(args.cache_dir)
        amcsd_reference_files = select_amcsd_candidate_files(
            cif_dir=cif_dir,
            index_path=index_path,
            observed_a_angstrom=fit.a_mean_angstrom,
            max_elements=args.max_elements,
            limit=args.amcsd_limit,
        )
        reference_files = sorted({*reference_files, *amcsd_reference_files})
        used_amcsd_cache = bool(amcsd_reference_files)
    observed_peak_pairs = [(peak.two_theta_deg, peak.intensity) for peak in selected_peaks]
    candidate_matches = identify_likely_crystals(
        observed_peaks=observed_peak_pairs,
        observed_a_angstrom=fit.a_mean_angstrom,
        observed_lattice_type=fit.lattice_type,
        wavelength_angstrom=args.wavelength,
        max_two_theta_deg=float(selected_peaks[-1].two_theta_deg),
        reference_files=reference_files,
        top_k=args.top_candidates,
    )
    add_candidate_matches(report, reference_files, candidate_matches, used_amcsd_cache=used_amcsd_cache)

    if args.json:
        print(json.dumps(report, indent=2))
        return

    print(format_text_report(report))


if __name__ == "__main__":
    main()

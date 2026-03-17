from __future__ import annotations

import math
import re
import shlex
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np


PERIODIC_TABLE = {
    symbol: index
    for index, symbol in enumerate(
        [
            "H",
            "He",
            "Li",
            "Be",
            "B",
            "C",
            "N",
            "O",
            "F",
            "Ne",
            "Na",
            "Mg",
            "Al",
            "Si",
            "P",
            "S",
            "Cl",
            "Ar",
            "K",
            "Ca",
            "Sc",
            "Ti",
            "V",
            "Cr",
            "Mn",
            "Fe",
            "Co",
            "Ni",
            "Cu",
            "Zn",
            "Ga",
            "Ge",
            "As",
            "Se",
            "Br",
            "Kr",
            "Rb",
            "Sr",
            "Y",
            "Zr",
            "Nb",
            "Mo",
            "Tc",
            "Ru",
            "Rh",
            "Pd",
            "Ag",
            "Cd",
            "In",
            "Sn",
            "Sb",
            "Te",
            "I",
            "Xe",
            "Cs",
            "Ba",
            "La",
            "Ce",
            "Pr",
            "Nd",
            "Pm",
            "Sm",
            "Eu",
            "Gd",
            "Tb",
            "Dy",
            "Ho",
            "Er",
            "Tm",
            "Yb",
            "Lu",
            "Hf",
            "Ta",
            "W",
            "Re",
            "Os",
            "Ir",
            "Pt",
            "Au",
            "Hg",
            "Tl",
            "Pb",
            "Bi",
            "Po",
            "At",
            "Rn",
            "Fr",
            "Ra",
            "Ac",
            "Th",
            "Pa",
            "U",
            "Np",
            "Pu",
            "Am",
            "Cm",
            "Bk",
            "Cf",
            "Es",
            "Fm",
            "Md",
            "No",
            "Lr",
            "Rf",
            "Db",
            "Sg",
            "Bh",
            "Hs",
            "Mt",
            "Ds",
            "Rg",
            "Cn",
            "Nh",
            "Fl",
            "Mc",
            "Lv",
            "Ts",
            "Og",
        ],
        start=1,
    )
}


@dataclass(frozen=True)
class AtomSite:
    element: str
    x: float
    y: float
    z: float
    occupancy: float


@dataclass(frozen=True)
class CifStructure:
    source_path: Path
    name: str
    formula: str
    space_group: str
    cell_a: float
    cell_b: float
    cell_c: float
    alpha: float
    beta: float
    gamma: float
    symmetry_ops: tuple[str, ...]
    atom_sites: tuple[AtomSite, ...]

    @property
    def is_cubic(self) -> bool:
        return (
            abs(self.cell_a - self.cell_b) < 1e-3
            and abs(self.cell_a - self.cell_c) < 1e-3
            and abs(self.alpha - 90.0) < 1e-3
            and abs(self.beta - 90.0) < 1e-3
            and abs(self.gamma - 90.0) < 1e-3
        )


@dataclass(frozen=True)
class SimulatedPeak:
    two_theta_deg: float
    intensity: float
    n_value: int


@dataclass(frozen=True)
class CandidateMatch:
    source_path: Path
    name: str
    formula: str
    space_group: str
    lattice_a_angstrom: float
    score: float
    probability: float
    matched_peak_count: int
    total_reference_peaks: int


def parse_numeric_token(token: str | None) -> float | None:
    if token is None:
        return None
    cleaned = token.strip().strip("'").strip('"')
    if cleaned in {"?", ".", ""}:
        return None
    cleaned = re.sub(r"\([0-9]+\)$", "", cleaned)
    try:
        return float(cleaned)
    except ValueError:
        return None


def tokenize_cif(text: str) -> list[str]:
    tokens: list[str] = []
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        raw_line = lines[index].rstrip("\n")
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            index += 1
            continue
        if stripped == ";":
            block: list[str] = []
            index += 1
            while index < len(lines) and lines[index].strip() != ";":
                block.append(lines[index].rstrip("\n"))
                index += 1
            tokens.append("\n".join(block).strip())
            index += 1
            continue
        tokens.extend(shlex.split(raw_line, posix=True))
        index += 1
    return tokens


def infer_element_symbol(raw_symbol: str) -> str:
    cleaned = raw_symbol.strip().strip("'").strip('"')
    if cleaned in PERIODIC_TABLE:
        return cleaned
    match = re.match(r"([A-Z][a-z]?)", cleaned)
    if not match:
        raise ValueError(f"Unable to infer element symbol from {raw_symbol!r}")
    symbol = match.group(1)
    if symbol not in PERIODIC_TABLE:
        raise ValueError(f"Unsupported element symbol {symbol!r}")
    return symbol


def parse_cif_file(path: Path) -> CifStructure:
    tokens = tokenize_cif(path.read_text(errors="ignore"))
    scalar_values: dict[str, str] = {}
    loops: list[tuple[list[str], list[list[str]]]] = []

    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token == "loop_":
            index += 1
            tags: list[str] = []
            while index < len(tokens) and tokens[index].startswith("_"):
                tags.append(tokens[index])
                index += 1
            rows: list[list[str]] = []
            width = len(tags)
            while index < len(tokens):
                token = tokens[index]
                if token == "loop_" or token.startswith("_") or token.startswith("data_"):
                    break
                row = tokens[index : index + width]
                if len(row) < width:
                    break
                rows.append(row)
                index += width
            loops.append((tags, rows))
            continue

        if token.startswith("_"):
            scalar_values[token] = tokens[index + 1] if index + 1 < len(tokens) else ""
            index += 2
            continue

        index += 1

    a = parse_numeric_token(scalar_values.get("_cell_length_a"))
    b = parse_numeric_token(scalar_values.get("_cell_length_b"))
    c = parse_numeric_token(scalar_values.get("_cell_length_c"))
    alpha = parse_numeric_token(scalar_values.get("_cell_angle_alpha"))
    beta = parse_numeric_token(scalar_values.get("_cell_angle_beta"))
    gamma = parse_numeric_token(scalar_values.get("_cell_angle_gamma"))
    if None in {a, b, c, alpha, beta, gamma}:
        raise ValueError(f"Incomplete unit cell metadata in {path}")

    name = (
        scalar_values.get("_chemical_name_mineral")
        or scalar_values.get("_chemical_name_common")
        or scalar_values.get("_chemical_formula_sum")
        or path.stem
    )
    formula = scalar_values.get("_chemical_formula_sum", "").strip()
    space_group = (
        scalar_values.get("_space_group_name_H-M_alt")
        or scalar_values.get("_symmetry_space_group_name_H-M")
        or ""
    ).strip()

    symmetry_ops: list[str] = []
    atom_sites: list[AtomSite] = []
    for tags, rows in loops:
        tag_set = set(tags)
        if "_space_group_symop_operation_xyz" in tag_set or "_symmetry_equiv_pos_as_xyz" in tag_set:
            symop_tag = "_space_group_symop_operation_xyz"
            if symop_tag not in tag_set:
                symop_tag = "_symmetry_equiv_pos_as_xyz"
            sym_index = tags.index(symop_tag)
            symmetry_ops.extend(row[sym_index] for row in rows)
            continue

        required = {"_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"}
        if not required.issubset(tag_set):
            continue

        x_index = tags.index("_atom_site_fract_x")
        y_index = tags.index("_atom_site_fract_y")
        z_index = tags.index("_atom_site_fract_z")
        occ_index = tags.index("_atom_site_occupancy") if "_atom_site_occupancy" in tag_set else None
        if "_atom_site_type_symbol" in tag_set:
            symbol_index = tags.index("_atom_site_type_symbol")
        elif "_atom_site_label" in tag_set:
            symbol_index = tags.index("_atom_site_label")
        else:
            raise ValueError(f"No atom type information found in {path}")

        for row in rows:
            x = parse_numeric_token(row[x_index])
            y = parse_numeric_token(row[y_index])
            z = parse_numeric_token(row[z_index])
            if None in {x, y, z}:
                continue
            occupancy = 1.0
            if occ_index is not None:
                occupancy = parse_numeric_token(row[occ_index]) or 1.0
            atom_sites.append(
                AtomSite(
                    element=infer_element_symbol(row[symbol_index]),
                    x=x,
                    y=y,
                    z=z,
                    occupancy=occupancy,
                )
            )

    if not atom_sites:
        raise ValueError(f"No atom sites found in {path}")

    if not symmetry_ops:
        symmetry_ops = ["x,y,z"]

    return CifStructure(
        source_path=path,
        name=name.strip("'").strip('"'),
        formula=formula.strip("'").strip('"'),
        space_group=space_group.strip("'").strip('"'),
        cell_a=float(a),
        cell_b=float(b),
        cell_c=float(c),
        alpha=float(alpha),
        beta=float(beta),
        gamma=float(gamma),
        symmetry_ops=tuple(symmetry_ops),
        atom_sites=tuple(atom_sites),
    )


def eval_symmetry_component(component: str, x: float, y: float, z: float) -> float:
    cleaned = component.replace(" ", "")
    if not cleaned:
        return 0.0
    if cleaned[0] not in "+-":
        cleaned = "+" + cleaned

    coords = {"x": x, "y": y, "z": z}
    total = 0.0
    for term in re.findall(r"[+-][^+-]+", cleaned):
        sign = -1.0 if term[0] == "-" else 1.0
        body = term[1:]
        if body in coords:
            total += sign * coords[body]
            continue
        total += sign * float(Fraction(body))
    return total


def expand_atom_sites(structure: CifStructure) -> list[AtomSite]:
    expanded: dict[tuple[str, int, int, int], float] = {}
    for site in structure.atom_sites:
        for operation in structure.symmetry_ops:
            components = [part.strip() for part in operation.split(",")]
            if len(components) != 3:
                continue
            x = eval_symmetry_component(components[0], site.x, site.y, site.z) % 1.0
            y = eval_symmetry_component(components[1], site.x, site.y, site.z) % 1.0
            z = eval_symmetry_component(components[2], site.x, site.y, site.z) % 1.0
            key = (site.element, round(x * 1_000_000), round(y * 1_000_000), round(z * 1_000_000))
            expanded[key] = site.occupancy

    sites = []
    for (element, x, y, z), occupancy in expanded.items():
        sites.append(
            AtomSite(
                element=element,
                x=x / 1_000_000.0,
                y=y / 1_000_000.0,
                z=z / 1_000_000.0,
                occupancy=occupancy,
            )
        )
    return sites


def simulate_powder_pattern(
    structure: CifStructure,
    wavelength_angstrom: float,
    max_two_theta_deg: float,
    intensity_threshold_percent: float = 0.5,
) -> list[SimulatedPeak]:
    if not structure.is_cubic:
        return []

    sites = expand_atom_sites(structure)
    positions = np.asarray([[site.x, site.y, site.z] for site in sites], dtype=float)
    occupancies = np.asarray([site.occupancy for site in sites], dtype=float)
    atomic_numbers = np.asarray([PERIODIC_TABLE[site.element] for site in sites], dtype=float)

    theta_max = math.radians(max_two_theta_deg / 2.0)
    n_max = int((2.0 * structure.cell_a * math.sin(theta_max) / wavelength_angstrom) ** 2) + 1
    h_limit = int(math.sqrt(max(n_max, 1))) + 1
    grouped_intensity: dict[int, float] = {}

    for h in range(-h_limit, h_limit + 1):
        for k in range(-h_limit, h_limit + 1):
            for l in range(-h_limit, h_limit + 1):
                n_value = h * h + k * k + l * l
                if n_value == 0 or n_value > n_max:
                    continue
                argument = wavelength_angstrom * math.sqrt(n_value) / (2.0 * structure.cell_a)
                if argument >= 1.0:
                    continue
                theta = math.asin(argument)
                two_theta_deg = math.degrees(2.0 * theta)
                if two_theta_deg > max_two_theta_deg + 1e-9:
                    continue

                sin_theta_over_lambda = math.sin(theta) / wavelength_angstrom
                scattering = atomic_numbers * np.exp(-8.0 * sin_theta_over_lambda * sin_theta_over_lambda)
                phase = np.exp(2j * math.pi * (h * positions[:, 0] + k * positions[:, 1] + l * positions[:, 2]))
                structure_factor = np.sum(occupancies * scattering * phase)
                intensity = float(abs(structure_factor) ** 2)
                if intensity < 1e-10:
                    continue
                grouped_intensity[n_value] = grouped_intensity.get(n_value, 0.0) + intensity

    if not grouped_intensity:
        return []

    max_intensity = max(grouped_intensity.values())
    peaks = []
    for n_value, intensity in sorted(grouped_intensity.items()):
        scaled_intensity = 100.0 * intensity / max_intensity
        if scaled_intensity < intensity_threshold_percent:
            continue
        theta = math.asin(wavelength_angstrom * math.sqrt(n_value) / (2.0 * structure.cell_a))
        peaks.append(
            SimulatedPeak(
                two_theta_deg=math.degrees(2.0 * theta),
                intensity=scaled_intensity,
                n_value=n_value,
            )
        )

    return peaks


def infer_lattice_letter(space_group: str) -> str:
    compact = space_group.replace(" ", "")
    for char in compact:
        if char in {"P", "I", "F", "R", "A", "B", "C"}:
            return char
    return "P"


def collect_reference_files(input_file: Path, reference_inputs: Sequence[Path] | None) -> list[Path]:
    search_roots: list[Path] = []
    if reference_inputs:
        search_roots.extend(reference_inputs)
    else:
        search_roots.extend([Path.cwd() / "references", input_file.parent])

    found: set[Path] = set()
    for root in search_roots:
        if root.is_file() and root.suffix.lower() == ".cif":
            found.add(root.resolve())
            continue
        if not root.exists() or not root.is_dir():
            continue
        for path in root.rglob("*.cif"):
            found.add(path.resolve())

    return sorted(found)


def score_candidate_pattern(
    observed_peaks: Sequence[tuple[float, float]],
    observed_a_angstrom: float,
    observed_lattice_type: str,
    candidate: CifStructure,
    predicted_peaks: Sequence[SimulatedPeak],
    position_tolerance_deg: float = 0.4,
) -> tuple[float, int]:
    if not predicted_peaks:
        return float("inf"), 0

    obs_max = max(intensity for _, intensity in observed_peaks)
    pred_max = max(peak.intensity for peak in predicted_peaks)
    obs_norm = [(two_theta, 100.0 * intensity / obs_max) for two_theta, intensity in observed_peaks]
    pred_norm = [(peak.two_theta_deg, 100.0 * peak.intensity / pred_max) for peak in predicted_peaks]

    matched_pred: set[int] = set()
    weighted_position_error = 0.0
    weighted_intensity_error = 0.0
    total_weight = 0.0
    matched_peak_count = 0

    for obs_two_theta, obs_intensity in obs_norm:
        weight = max(obs_intensity / 100.0, 0.02)
        best_index = None
        best_distance = None
        for index, (pred_two_theta, _pred_intensity) in enumerate(pred_norm):
            distance = abs(pred_two_theta - obs_two_theta)
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_index = index

        total_weight += weight
        if best_index is None or best_distance is None or best_distance > position_tolerance_deg:
            weighted_position_error += weight * 1.5
            weighted_intensity_error += weight * 1.0
            continue

        matched_pred.add(best_index)
        matched_peak_count += 1
        pred_two_theta, pred_intensity = pred_norm[best_index]
        weighted_position_error += weight * (best_distance / position_tolerance_deg)
        weighted_intensity_error += weight * (abs(pred_intensity - obs_intensity) / 100.0)

    strong_pred_total = 0.0
    strong_pred_unmatched = 0.0
    for index, (_pred_two_theta, pred_intensity) in enumerate(pred_norm[: max(len(obs_norm), 12)]):
        if pred_intensity < 5.0:
            continue
        strong_pred_total += pred_intensity
        if index not in matched_pred:
            strong_pred_unmatched += pred_intensity

    coverage_penalty = strong_pred_unmatched / strong_pred_total if strong_pred_total else 0.0
    a_penalty = abs(candidate.cell_a - observed_a_angstrom) / max(0.02, observed_a_angstrom * 0.01)
    lattice_letter = infer_lattice_letter(candidate.space_group)
    observed_letter = {
        "primitive": "P",
        "bcc": "I",
        "fcc": "F",
        "diamond": "F",
    }.get(observed_lattice_type, "P")
    lattice_penalty = 0.0 if lattice_letter == observed_letter else 0.35

    score = (
        1.2 * (weighted_position_error / max(total_weight, 1e-9))
        + 0.35 * (weighted_intensity_error / max(total_weight, 1e-9))
        + 0.45 * coverage_penalty
        + 0.9 * a_penalty
        + lattice_penalty
    )
    return score, matched_peak_count


def identify_likely_crystals(
    observed_peaks: Sequence[tuple[float, float]],
    observed_a_angstrom: float,
    observed_lattice_type: str,
    wavelength_angstrom: float,
    max_two_theta_deg: float,
    reference_files: Sequence[Path],
    top_k: int = 5,
) -> list[CandidateMatch]:
    scored: list[tuple[CifStructure, float, int, int]] = []

    for path in reference_files:
        try:
            structure = parse_cif_file(path)
        except Exception:
            continue
        if not structure.is_cubic:
            continue
        if abs(structure.cell_a - observed_a_angstrom) / observed_a_angstrom > 0.25:
            continue

        predicted_peaks = simulate_powder_pattern(
            structure,
            wavelength_angstrom=wavelength_angstrom,
            max_two_theta_deg=max_two_theta_deg,
        )
        if not predicted_peaks:
            continue

        score, matched_peak_count = score_candidate_pattern(
            observed_peaks=observed_peaks,
            observed_a_angstrom=observed_a_angstrom,
            observed_lattice_type=observed_lattice_type,
            candidate=structure,
            predicted_peaks=predicted_peaks,
        )
        if not math.isfinite(score):
            continue
        scored.append((structure, score, matched_peak_count, len(predicted_peaks)))

    if not scored:
        return []

    scored.sort(key=lambda item: item[1])
    deduplicated: list[tuple[CifStructure, float, int, int]] = []
    seen_keys: set[tuple[str, str, str]] = set()
    for structure, score, matched_peak_count, ref_peak_count in scored:
        key = (
            structure.name.strip().lower(),
            structure.formula.strip().lower(),
            structure.space_group.strip().lower(),
        )
        if key in seen_keys:
            continue
        seen_keys.add(key)
        deduplicated.append((structure, score, matched_peak_count, ref_peak_count))
        if len(deduplicated) >= top_k:
            break

    kept = deduplicated
    logits = np.asarray([-item[1] for item in kept], dtype=float)
    logits -= float(logits.max())
    probabilities = np.exp(logits)
    probabilities /= float(probabilities.sum())

    matches = []
    for (structure, score, matched_peak_count, ref_peak_count), probability in zip(kept, probabilities, strict=True):
        matches.append(
            CandidateMatch(
                source_path=structure.source_path,
                name=structure.name,
                formula=structure.formula,
                space_group=structure.space_group,
                lattice_a_angstrom=structure.cell_a,
                score=float(score),
                probability=float(probability),
                matched_peak_count=matched_peak_count,
                total_reference_peaks=ref_peak_count,
            )
        )

    return matches

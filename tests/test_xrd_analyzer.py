from pathlib import Path
import tempfile
import unittest

import numpy as np

from cif_matcher import collect_reference_files, identify_likely_crystals, parse_cif_file, simulate_powder_pattern
from xrd_analyzer import detect_peaks, fit_cubic_lattice, parse_int_file


class XrdAnalyzerTests(unittest.TestCase):
    SIMPLE_CUBIC_TEMPLATE = """data_test
_chemical_name_mineral '{name}'
_chemical_formula_sum '{formula}'
_cell_length_a {a}
_cell_length_b {a}
_cell_length_c {a}
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P m 3 m'
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
{element} 0 0 0
"""

    def test_parse_int_file_skips_header(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample.int"
            path.write_text(
                "GENERAL$\n"
                "3\n"
                "10.000000 1.0 0.0\n"
                "10.010000 2.0 0.0\n"
                "10.020000 3.0 0.0\n"
            )
            two_theta, intensity = parse_int_file(path)

        self.assertEqual(two_theta.tolist(), [10.0, 10.01, 10.02])
        self.assertEqual(intensity.tolist(), [1.0, 2.0, 3.0])

    def test_detect_peaks_and_fit_simple_cubic_pattern(self) -> None:
        two_theta = np.linspace(10.0, 90.0, 8001)
        intensity = np.zeros_like(two_theta)
        wavelength = 1.0
        a_value = 4.0
        n_values = [1, 2, 3, 4]
        two_theta_positions = []
        for n_value in n_values:
            d_spacing = a_value / np.sqrt(n_value)
            theta = np.arcsin(wavelength / (2.0 * d_spacing))
            two_theta_positions.append(np.degrees(2.0 * theta))

        for center, height in zip(two_theta_positions, [30.0, 100.0, 40.0, 60.0]):
            intensity += height * np.exp(-0.5 * ((two_theta - center) / 0.08) ** 2)

        peaks = detect_peaks(two_theta, intensity, relative_threshold=0.02, merge_tolerance_deg=0.2)
        fit = fit_cubic_lattice(peaks[:4], wavelength_angstrom=wavelength)

        self.assertEqual(
            [round(peak.two_theta_deg, 1) for peak in peaks[:4]],
            [round(value, 1) for value in two_theta_positions],
        )
        self.assertEqual(fit.lattice_type, "primitive")
        self.assertEqual(fit.assigned_n_values[:4], (1, 2, 3, 4))

    def test_identify_likely_crystal_prefers_best_matching_cif(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            matching_cif = tmp_path / "match.cif"
            mismatch_cif = tmp_path / "mismatch.cif"
            matching_cif.write_text(
                self.SIMPLE_CUBIC_TEMPLATE.format(name="Matchite", formula="Na", a="4.0", element="Na")
            )
            mismatch_cif.write_text(
                self.SIMPLE_CUBIC_TEMPLATE.format(name="Mismatchite", formula="Na", a="4.8", element="Na")
            )

            structure = parse_cif_file(matching_cif)
            predicted = simulate_powder_pattern(structure, wavelength_angstrom=1.0, max_two_theta_deg=90.0)
            two_theta = np.linspace(10.0, 90.0, 8001)
            intensity = np.zeros_like(two_theta)
            for peak in predicted[:6]:
                intensity += peak.intensity * np.exp(-0.5 * ((two_theta - peak.two_theta_deg) / 0.08) ** 2)

            observed_peaks = detect_peaks(two_theta, intensity, relative_threshold=0.02, merge_tolerance_deg=0.2)
            fit = fit_cubic_lattice(observed_peaks[:6], wavelength_angstrom=1.0)
            reference_files = collect_reference_files(tmp_path / "pattern.int", [tmp_path])
            matches = identify_likely_crystals(
                observed_peaks=[(peak.two_theta_deg, peak.intensity) for peak in observed_peaks[:6]],
                observed_a_angstrom=fit.a_mean_angstrom,
                observed_lattice_type=fit.lattice_type,
                wavelength_angstrom=1.0,
                max_two_theta_deg=90.0,
                reference_files=reference_files,
                top_k=2,
            )

        self.assertGreaterEqual(len(matches), 2)
        self.assertEqual(matches[0].name, "Matchite")
        self.assertGreater(matches[0].probability, matches[1].probability)


if __name__ == "__main__":
    unittest.main()

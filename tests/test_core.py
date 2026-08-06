import numpy as np

from ti_vpl.analytic import PreviewParameters, preview_metrics
from ti_vpl.coordinates import DEFAULT_VPL, crossed_montage
from ti_vpl.envelope import directional_envelope, maximum_orientation_envelope
from ti_vpl.fem import TetMesh, assemble_stiffness, solve_current_pair
from ti_vpl.registration import apply_affine, fit_affine
from scripts.convert_efmouse_mesh import normalize_node_elem


def test_updated_default_vpl_coordinate():
    assert DEFAULT_VPL.AP == -1.94
    assert DEFAULT_VPL.ML == -1.75
    assert DEFAULT_VPL.DV == -3.10


def test_xu_lab_crossed_montage_has_exact_default_coordinates_and_geometry():
    montage = crossed_montage(DEFAULT_VPL)
    expected = np.asarray([
        [0.14, -1.75, 0.50],
        [-4.02, -1.75, 0.50],
        [-1.94, 0.33, 0.50],
        [-1.94, -3.83, 0.50],
    ])
    actual = np.vstack([p.as_array() for pair in montage.values() for p in pair])
    assert np.allclose(actual, expected)
    for pair in montage.values():
        midpoint = (pair[0].as_array() + pair[1].as_array()) / 2
        assert np.allclose(midpoint[:2], DEFAULT_VPL.as_array()[:2])
        trajectory = np.linalg.norm(pair[0].as_array() - DEFAULT_VPL.as_array())
        angle = np.degrees(np.arctan2(3.60, 2.08))
        assert np.isclose(trajectory, 4.16, atol=0.01)
        assert np.isclose(angle, 60.0, atol=0.1)


def test_envelope_aligned_equal():
    e = np.array([[1.0, 0.0, 0.0]])
    assert np.allclose(maximum_orientation_envelope(e, e), 2.0)


def test_envelope_aligned_unequal():
    e1 = np.array([[1.0, 0.0, 0.0]])
    e2 = np.array([[3.0, 0.0, 0.0]])
    assert np.allclose(maximum_orientation_envelope(e1, e2), 2.0)


def test_envelope_equal_orthogonal():
    e1 = np.array([[1.0, 0.0, 0.0]])
    e2 = np.array([[0.0, 1.0, 0.0]])
    assert np.allclose(maximum_orientation_envelope(e1, e2), np.sqrt(2.0))
    assert np.allclose(directional_envelope(e1, e2, np.array([1.0, 0.0, 0.0])), 0.0)


def test_envelope_zero_channel():
    assert np.allclose(maximum_orientation_envelope(np.zeros((2, 3)), np.ones((2, 3))), 0.0)


def test_preview_scales_linearly_with_current():
    one = preview_metrics(PreviewParameters(channel_1_current_mA=0.25, channel_2_current_mA=0.25), shape=(26, 27, 22))
    two = preview_metrics(PreviewParameters(channel_1_current_mA=0.50, channel_2_current_mA=0.50), shape=(26, 27, 22))
    assert np.isclose(two["vpl_median_V_per_m"] / one["vpl_median_V_per_m"], 2.0, rtol=1e-10)


def test_equal_frequencies_have_no_ti_beat_envelope():
    metrics = preview_metrics(PreviewParameters(f1_Hz=2000.0, f2_Hz=2000.0), shape=(26, 27, 22))
    assert metrics["vpl_median_V_per_m"] == 0.0


def test_single_tetra_stiffness_is_symmetric_with_zero_row_sum():
    points = np.array([[0, 0, 0], [1e-3, 0, 0], [0, 1e-3, 0], [0, 0, 1e-3]], dtype=float)
    tets = np.array([[0, 1, 2, 3]], dtype=int)
    matrix, volumes, _ = assemble_stiffness(points, tets, np.array([0.3]))
    dense = matrix.toarray()
    assert np.allclose(dense, dense.T)
    assert np.allclose(dense.sum(axis=1), 0.0, atol=1e-14)
    assert np.isclose(volumes[0], 1e-9 / 6.0)


def test_registration_recovers_known_affine_with_submicron_residual():
    source = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, -3], [2, -1, -2]], dtype=float)
    affine = np.array([[1.1, .02, 0, 4], [-.01, .9, .03, -2], [0, .04, 1.05, 7], [0, 0, 0, 1]], dtype=float)
    target = apply_affine(source, affine)
    result = fit_affine(source, target)
    assert np.allclose(result.affine, affine, atol=1e-12)
    assert result.rms_mm < 1e-12


def test_efmouse_converter_accepts_matlab_transposed_one_based_arrays():
    node = np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=float)
    elem = np.array([[1], [2], [3], [4], [1]], dtype=float)
    points, tets, labels = normalize_node_elem(node, elem)
    assert points.shape == (4, 3)
    assert np.array_equal(tets, [[0, 1, 2, 3]])
    assert labels.tolist() == [1]


def test_single_tetra_neumann_pair_is_balanced_and_solved():
    points = np.array([[0, 0, 0], [1e-3, 0, 0], [0, 1e-3, 0], [0, 0, 1e-3]], dtype=float)
    mesh = TetMesh(points, np.array([[0, 1, 2, 3]]), np.array([1]))
    result = solve_current_pair(
        mesh, {1: .275}, np.array([[0, 1, 2]]), np.array([[0, 1, 3]]),
        current_A=1e-4, direct_max_nodes=20,
    )
    assert float(result["boundary_current_imbalance_A"]) < 1e-15
    assert float(result["linear_relative_residual"]) < 1e-10
    assert np.all(np.isfinite(result["E_V_per_m"]))

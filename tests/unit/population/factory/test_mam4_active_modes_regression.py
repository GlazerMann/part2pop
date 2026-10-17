import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")
from part2pop.population.factory import mam4
from part2pop.utilities import power_moments_from_lognormal
from part2pop.constants import DENSITY_LIQUID_WATER, MOLAR_MASS_DRY_AIR, R


def _output(tmp_path, numbers):
    with netCDF4.Dataset(tmp_path / "mam_output.nc", "w") as ds:
        ds.createDimension("mode", 4)
        ds.createDimension("time", 1)
        values = {
            "num_aer": numbers,
            "so4_aer": [1e-9] * 4,
            "soa_aer": [2e-9] * 4,
            "dgn_a": [1e-7, 2e-7, 3e-7, 4e-7],
            "dgn_awet": [1.2e-7, 2.4e-7, 3.6e-7, 4.8e-7],
        }
        for name, data in values.items():
            ds.createVariable(name, "f8", ("mode", "time"))[:, 0] = data
    return {"mam4_dir": str(tmp_path), "timestep": 2,
            "GSD": [1.2, 1.4, 1.6, 1.8], "N_bins": [2, 3, 4, 5],
            "p": 90000.0, "T": 298.0}


@pytest.mark.parametrize("wet", [True, False])
@pytest.mark.parametrize("bins", [[2, 3, 4, 5], 3])
def test_middle_inactive_mode_preserves_mode_indices(tmp_path, monkeypatch, wet, bins):
    cfg = _output(tmp_path, [1e6, 0.0, 3e6, 4e6])
    cfg.update(D_is_wet=wet, N_bins=bins)
    captured = {}

    def capture(config):
        captured.update(config)
        return "built"

    monkeypatch.setattr(mam4, "build_binned_lognormals", capture)
    assert mam4.build(cfg) == "built"
    assert captured["GSD"] == [1.2, 1.6, 1.8]
    assert captured["N_bins"] == ([2, 4, 5] if isinstance(bins, list) else 3)
    expected_names = ["SO4", "OC", "H2O"] if wet else ["SO4", "OC"]
    assert captured["aero_spec_names"] == [expected_names] * 3
    np.testing.assert_allclose(captured["GMD"],
                               np.array([1e-7, 3e-7, 4e-7]) * (1.2 if wet else 1.0),
                               rtol=1e-12, atol=0.0)
    if wet:
        density = MOLAR_MASS_DRY_AIR * cfg["p"] / (R * cfg["T"])
        for row, n, dry_d, gsd in zip(captured["aero_spec_fracs"],
                                    np.array([1e6, 3e6, 4e6]) * density,
                                    [1e-7, 3e-7, 4e-7], [1.2, 1.6, 1.8]):
            water = DENSITY_LIQUID_WATER * np.pi / 6.0 * (
                power_moments_from_lognormal(3, n, dry_d * 1.2, gsd)
                - power_moments_from_lognormal(3, n, dry_d, gsd))
            expected = np.array([1e-9 * density, 2e-9 * density, water])
            np.testing.assert_allclose(row, expected / expected.sum(), rtol=1e-12, atol=0.0)


def test_real_builder_accepts_middle_inactive_mode(tmp_path):
    cfg = _output(tmp_path, [1e6, 0.0, 3e6, 4e6])
    pop = mam4.build(cfg)
    assert len(pop.ids) == 2 + 4 + 5
    expected = 8e6 * MOLAR_MASS_DRY_AIR * cfg["p"] / (R * cfg["T"])
    assert pop.get_Ntot() == pytest.approx(expected)


def test_all_inactive_modes_return_empty_population(tmp_path):
    cfg = _output(tmp_path, [0.0] * 4)
    pop = mam4.build(cfg)
    assert len(pop.ids) == 0
    assert pop.get_Ntot() == 0.0


@pytest.mark.parametrize("field,value", [("GSD", [1.2, 1.4, 1.6]),
                                         ("N_bins", [2, 3, 4])])
def test_mode_parameters_must_describe_original_modes(tmp_path, field, value):
    cfg = _output(tmp_path, [1e6, 0.0, 3e6, 4e6])
    cfg[field] = value
    with pytest.raises(ValueError, match=field + " must have one entry per original MAM mode"):
        mam4.build(cfg)

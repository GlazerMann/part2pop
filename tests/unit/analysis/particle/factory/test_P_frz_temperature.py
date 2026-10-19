import numpy as np
import pytest
from part2pop.analysis.builder import build_variable

from part2pop import ParticlePopulation, make_particle


def population():
    p = make_particle(100e-9, ["SO4"], [1.0])
    return ParticlePopulation(p.species, np.array([p.masses]), np.array([1.0]), [0])


@pytest.mark.parametrize("celsius", [0.0, -20.0])
def test_celsius_probability_matches_kelvin(celsius):
    c = build_variable(
        "P_frz", scope="particle", var_cfg={"T": celsius, "T_units": "C"}
    ).compute_all(population())
    k = build_variable(
        "P_frz", scope="particle", var_cfg={"T": celsius + 273.15, "T_units": "K"}
    ).compute_all(population())
    assert np.isfinite(c).all()
    np.testing.assert_allclose(c, k, rtol=1e-12, atol=0.0)


@pytest.mark.parametrize(
    "temperature,units",
    [
        (None, "K"),
        (0.0, "K"),
        (-1.0, "K"),
        (np.nan, "K"),
        (np.inf, "K"),
        (-273.15, "C"),
        (250.0, "invalid"),
    ],
)
def test_invalid_temperature_is_rejected_before_water_mutation(temperature, units):
    pop = population()
    before = pop.spec_masses.copy()
    with pytest.raises(ValueError):
        build_variable(
            "P_frz", scope="particle", var_cfg={"T": temperature, "T_units": units}
        ).compute_all(pop)
    np.testing.assert_array_equal(pop.spec_masses, before)


@pytest.mark.parametrize("inherited", [False, True])
def test_zero_celsius_preserves_rh_equilibration_and_coefficient_overrides(inherited):
    from part2pop.freezing.builder import build_freezing_population

    p = make_particle(100e-9, ["BC", "SO4"], [0.5, 0.5])

    def make_pop():
        return ParticlePopulation(p.species, np.array([p.masses]), np.array([1.0]), [0])

    pop, control = make_pop(), make_pop()
    modifications = {"BC": {"m_log10_Jhet": 0.0, "b_log10_Jhet": 14.0}}
    cfg = {"T": 0.0, "T_units": "C", "RH": 0.7}
    if inherited:
        pop.species_modifications = modifications
    else:
        cfg["species_modifications"] = modifications
    actual = build_variable("P_frz", scope="particle", var_cfg=cfg).compute_all(pop)
    control._equilibrate_h2o(0.7, 273.15)
    freezing_cfg = {"morphology": "homogeneous", "species_modifications": modifications}
    frozen = build_freezing_population(control, freezing_cfg)
    expected = frozen.get_freezing_probs(273.15, freezing_cfg)
    assert 0.0 < actual[0] < 1.0
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(
        pop.spec_masses, control.spec_masses, rtol=1e-12, atol=0.0
    )

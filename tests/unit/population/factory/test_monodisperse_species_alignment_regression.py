import numpy as np
import pytest
from part2pop.population.builder import build_population


def test_permuted_species_order_is_aligned_before_particle_construction():
    cfg = {
        "type": "monodisperse", "N": [1., 1.], "D": [1e-7, 1e-7],
        "aero_spec_names": [["SO4", "OC"], ["OC", "SO4"]],
        "aero_spec_fracs": [[.9, .1], [.8, .2]],
    }
    population = build_population(cfg)
    names = [s.name for s in population.species]
    masses = population.spec_masses
    for i, expected in enumerate((.9, .2)):
        dry = masses[i, names.index("SO4")] + masses[i, names.index("OC")]
        assert masses[i, names.index("SO4")] / dry == pytest.approx(expected)


def test_shared_species_row_remains_compatible():
    pop = build_population({"type":"monodisperse", "N":[1., 2.],
                            "D":[1e-7, 2e-7], "aero_spec_names":[["SO4"]],
                            "aero_spec_fracs":[[1.], [1.]]})
    assert len(pop.ids) == 2


def test_species_row_length_mismatch_fails_early():
    with pytest.raises(ValueError, match="mismatched"):
        build_population({"type":"monodisperse", "N":[1.], "D":[1e-7],
                          "aero_spec_names":[["SO4", "OC"]],
                          "aero_spec_fracs":[[1.]]})

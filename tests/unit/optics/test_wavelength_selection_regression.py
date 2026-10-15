from part2pop.optics.base import OpticalPopulation
from part2pop.population.base import ParticlePopulation
import numpy as np
import pytest


def test_select_505nm_does_not_match_500nm():
    base = ParticlePopulation(species=(), spec_masses=np.empty((1, 0)),
                              num_concs=np.array([1.]), ids=[1])
    pop = OpticalPopulation(base, [0.], [500e-9, 505e-9])
    pop.Cext[0, 0, :] = [7., 13.]
    assert pop.get_optical_coeff('b_ext', rh=0, wvl=505e-9) == 13.
    with pytest.raises(ValueError, match="Wavelength"):
        pop.get_optical_coeff('b_ext', rh=0, wvl=502e-9)

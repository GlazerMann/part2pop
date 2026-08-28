# optics/factory/homogeneous.py
import numpy as np
import math
import warnings

from .registry import register
from ..base import OpticalParticle
from ..refractive_index import build_refractive_index

try:
    from part2pop._patch import patch_pymiescatt
    patch_pymiescatt()
    from PyMieScatt import MieQ
    _PMS_ERR = None
except Exception as e:
    MieQ = None
    _PMS_ERR = e



@register("homogeneous")
class HomogeneousParticle(OpticalParticle):
    """
    Homogeneous sphere morphology optical particle model with RH and wavelength dependence.

    Constructor expects (base_particle, config) to align with the factory builder.

    Optional config (read by OpticalParticle or here):
      - rh_grid, wvl_grid, temp (K), specdata_path, species_modifications

    If PyMieScatt is unavailable, a Rayleigh-sphere approximation is used.
    """

    def __init__(self, base_particle, config):
        super().__init__(base_particle, config)

        # Refractive indices are attached at the population level by the
        # optics builder; the base class's _attach_refractive_indices is
        # guarded and will no-op if the species already have wavelength-aware
        # RIs. Keep the call to the base preparation intact.

        # Precompute geometry & per-wavelength dry/water RIs
        self._prepare_geometry_and_ris()

        # Do the optics
        self.compute_optics()

    def _prepare_geometry_and_ris(self):
        """
        Precompute:
          - dry particle volume
          - water volumes vs RH (from Dwet - Ddry)
          - wavelength-dependent RIs for dry mix and water
        """
        # Water volumes vs RH from Dwet(RH) and Ddry
        Ddry = float(self.get_Ddry())

        self.h2o_vols = np.zeros(len(self.rh_grid))
        for rr, rh in enumerate(self.rh_grid):
            Dw = float(self.get_Dwet(RH=float(rh), T=self.temp))
            self.h2o_vols[rr] = (math.pi / 6.0) * (Dw ** 3 - Ddry ** 3) if rh > 0.0 else 0.0

        # Per-wavelength complex RIs
        Nw = len(self.wvl_grid)
        self.dry_ris = np.zeros(Nw, dtype=complex)
        self.h2o_ris = np.zeros(Nw, dtype=complex)

        # Species partitioning (use Particle-provided helpers if available)
        vks = self.get_vks()  # species "dry" volumes for partitioning
        h2o_idx = self.idx_h2o()

        # Build dry mixture RI by volume-weighted averaging of all non-water species
        dry_indices = [ii for ii in range(len(self.species)) if ii != h2o_idx]
        self.dry_vol = np.sum([vks[ii] for ii in dry_indices])
        for ww in range(Nw):
            # fixme: should this go in some effective_ri helper?
            n_dry = 0.0
            k_dry = 0.0
            for ii in dry_indices:
                f = float(vks[ii] / self.dry_vol)
                n_dry += self.species[ii].refractive_index.real_ri_fun(self.wvl_grid[ww]) * f
                k_dry += self.species[ii].refractive_index.imag_ri_fun(self.wvl_grid[ww]) * f
            self.dry_ris[ww] = complex(n_dry, k_dry)
            # Water RI (always set)
            n_w = self.species[h2o_idx].refractive_index.real_ri_fun(self.wvl_grid[ww])
            k_w = self.species[h2o_idx].refractive_index.imag_ri_fun(self.wvl_grid[ww])
            self.h2o_ris[ww] = complex(n_w, k_w)
    
    def _mixture_ri(self, rr: int, ww: int) -> complex:
        """
        Volume-weighted homogeneous mixture of dry material and water at a given RH and wavelength.
        """
        v_h2o = self.h2o_vols[rr]
        v_dry = self.dry_vol
        if (v_h2o + v_dry) <= 0.0:
            return complex(1.0, 0.0)
        return (self.h2o_ris[ww] * v_h2o + self.dry_ris[ww] * v_dry) / (v_h2o + v_dry)

    @staticmethod
    def _rayleigh_cross_sections(m: complex, wavelength_m: float, radius_m: float):
        """Return Rayleigh Cext, Csca, Cabs, and g for a homogeneous sphere.

        This is the small-particle approximation and is most accurate when
        x = 2*pi*r/lambda << 1.
        """
        k = 2.0 * math.pi / wavelength_m
        m2 = m * m
        polarizability_factor = (m2 - 1.0) / (m2 + 2.0)

        csca = (
            (8.0 * math.pi / 3.0)
            * (k ** 4)
            * (radius_m ** 6)
            * abs(polarizability_factor) ** 2
        )
        cabs = (
            4.0
            * math.pi
            * k
            * (radius_m ** 3)
            * max(0.0, float(polarizability_factor.imag))
        )
        cext = csca + cabs

        # Rayleigh scattering is symmetric fore/aft, so <cos(theta)> = 0.
        return cext, csca, cabs, 0.0

    def compute_optics(self):
        """
        Compute cross-sections and asymmetry parameter per (RH, wavelength).
        Prefer PyMieScatt if available; otherwise use the Rayleigh approximation.
        """
        if MieQ is None:
            warnings.warn(
                "PyMieScatt is unavailable; homogeneous optics are using the "
                "Rayleigh-sphere approximation. This fallback is only reliable "
                "when the particle size parameter x = 2*pi*r/lambda is much "
                "smaller than 1.",
                RuntimeWarning,
                stacklevel=2,
            )

        for rr, rh in enumerate(self.rh_grid):
            D_m = float(self.get_Dwet(RH=float(rh), T=self.temp, sigma_sa=self.get_surface_tension()))
            r_m = 0.5 * D_m
            area = math.pi * r_m * r_m  # geometric cross-section

            D_nm = D_m * 1e9
            for ww, lam_m in enumerate(self.wvl_grid):
                lam_nm = float(lam_m * 1e9)
                m = complex(self._mixture_ri(rr, ww))

                if MieQ is not None:
                    out = MieQ(m, lam_nm, D_nm, asDict=True, asCrossSection=False)
                    # Convert efficiencies to absolute cross sections.
                    geom_area = np.pi / 4.0 * D_nm ** 2 * 1e-18
                    self.Cext[rr, ww] = out["Qext"] * geom_area
                    self.Csca[rr, ww] = out["Qsca"] * geom_area
                    self.Cabs[rr, ww] = out["Qabs"] * geom_area
                    self.g[rr, ww] = out["g"]
                else:
                    cext, csca, cabs, g = self._rayleigh_cross_sections(
                        m=m,
                        wavelength_m=float(lam_m),
                        radius_m=r_m,
                    )
                    self.Cext[rr, ww] = cext
                    self.Csca[rr, ww] = csca
                    self.Cabs[rr, ww] = cabs
                    self.g[rr, ww] = g
    
    # Convenience getters (unchanged)
    def get_cross_sections(self):
        return {
            "Cabs": self.Cabs,
            "Csca": self.Csca,
            "Cext": self.Cext,
            "g": self.g,
        }

    def get_refractive_indices(self):
        return {"dry_ri": self.dry_ris, "h2o_ri": self.h2o_ris}

    def get_cross_section(self, optics_type, rh_idx=None, wvl_idx=None):
        key = str(optics_type).lower()
        if key in ("b_abs", "absorption", "abs"):
            arr = self.Cabs
        elif key in ("b_scat", "scattering", "scat"):
            arr = self.Csca
        elif key in ("b_ext", "extinction", "ext"):
            arr = self.Cext
        elif key in ("g", "asymmetry"):
            arr = self.g
        else:
            raise ValueError(f"Unknown optics_type: {optics_type}")
        if rh_idx is not None and wvl_idx is not None:
            return arr[rh_idx, wvl_idx]
        return arr


def build(base_particle, config):
    """Optional fallback factory callable for discovery."""
    return HomogeneousParticle(base_particle, config)

from __future__ import annotations
from ..base import PopulationVariable, VariableMeta
from .registry import register_variable
from part2pop.freezing.builder import build_freezing_population
import numpy as np

@register_variable("nucleating_sites")
class NucleatingSites(PopulationVariable):
    """Integrate heterogeneous rates using the supplied particle water content.

    No RH equilibration is performed here. Omitted species_modifications
    inherit population overrides; an explicit mapping, including {}, replaces
    them. The cooling-integral definition still requires scientific review.
    """
    meta = VariableMeta(
        name="nucleating_sites",
        description='Cumulative heterogeneous ice-nucleation active sites per surface area during cooling.',
        units = r'm$^{-2}$',
        axis_names=("T_grid",),
        default_cfg={},
        aliases = ('n_s',),
        scale = 'log',
        short_label = 'n_s',
        long_label = 'active site density',
    )

    def compute(self, population, as_dict=False):
        cfg = self.cfg
        units = cfg.get("T_units", "K")
        freezing_cfg = {
            "T_grid": list(cfg["T_grid"]),
            "morphology": cfg.get("morphology", "homogeneous"),
            "species_modifications": cfg.get(
                "species_modifications", getattr(population, "species_modifications", {})
            ) or {},
            "T_units": units
        }
        freezing_pop = build_freezing_population(population, freezing_cfg)
        arr = freezing_pop.get_nucleating_sites(
            cfg["T_grid"], cfg["cooling_rate"], freezing_cfg, T_units=units
        )
        if as_dict:
            return {"T_grid": np.asarray(cfg["T_grid"]), "cooling_rate": cfg["cooling_rate"],
                    "T_units": cfg.get("T_units", "K"), "nucleating_sites": arr,}
        return arr


def build(cfg=None):
    cfg = cfg or {}
    return NucleatingSites(cfg)

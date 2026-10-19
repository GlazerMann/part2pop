# Freezing analysis species overrides

`P_frz`, `avg_Jhet`, `frozen_frac`, `unfrozen_frac`, and `nucleating_sites`
use the same authoritative species-modification policy:

| Analysis configuration | Effective mapping |
| --- | --- |
| Key omitted | Population `species_modifications` |
| `species_modifications=None` | Population `species_modifications` |
| `species_modifications={}` | Empty mapping; tabulated/default coefficients |
| Explicit nonempty mapping | Exactly that mapping; no per-species merge |

Overrides apply to one calculation. These analysis entry points do not replace
the population's modification mapping or edit either input mapping. The shared
resolver creates an analysis-local configuration and borrows the selected mapping;
freezing construction treats that mapping as read-only. Direct freezing-builder
APIs retain their existing behavior.

Humidity equilibration retains its existing water-mass side effects. This patch
does not make all diagnostic calculations pure or change temperature conversion,
freezing geometry, rate formulas, stochastic treatment, or time stepping.

The separate P05 custom-area patch gives `INSA_distribution` this same policy;
apply both to cover all six analysis entry points. P11 is independently applicable
to the pinned upstream base and supplies the other five entry points. The missing
nucleating-site numerical API is repaired by existing PR #117, not by this patch.

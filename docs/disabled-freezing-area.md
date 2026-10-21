# Disabled freezing material and active area

A species contributes its positive, finite volume to the equivalent-sphere
active surface area only when both effective heterogeneous coefficients
(`m_log10_Jhet` and `b_log10_Jhet`) are finite. Both tabulated species and
custom species with complete override pairs follow this rule. A zero/zero pair
is valid, whereas NaN or infinite coefficients disable the material.

Volumes are selected by species position, not substring name matching, so
similarly named materials do not contaminate each other's active area.
The area is calculated from the pooled active volume, not by summing
individual particle surface areas.

PR #127 implements matching custom-material and volume-alignment behavior.
PR #135 independently fixes exact freezing-table lookup. Merge #135 first;
then reconcile #127 and #136 against their shared implementation before the
second of those PRs merges. PR #123 independently handles heterogeneous rate
calculation for absent or disabled active species.

The physical coefficient fits, morphology and water-activity relation are
unchanged. The mapping between disabled materials and solutes is not changed
by this active-area calculation.

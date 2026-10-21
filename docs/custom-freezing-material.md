# Custom heterogeneous active material

A complete custom `m_log10_Jhet`/`b_log10_Jhet` pair contributes material volume
to the same equivalent-sphere ice-active area as tabulated material. Zero-valued
coefficients are valid. Missing or partial custom pairs retain their existing
warning/fallback behavior. This changes active geometry, not water activity or
the morphology's area convention.
Volumes are taken from the species-aligned mass and density arrays, so material
names such as `Na` and `NaCl` remain distinct.

A tabulated material whose effective coefficients include NaN is disabled and
adds no area, as its rate is already excluded. #136 uses the same active-area
code and also keeps disabled materials out of the water-activity solutes.
Without #135, the table lookup matches name prefixes, so a custom species named
`O` or `B` picks up the OC or biological coefficients; merge #135 first.

`INSA_distribution` forwards the authoritative species-modification mapping to
freezing construction. An absent or `None` mapping inherits population defaults;
an explicit mapping replaces them, and an explicit empty mapping disables them.
The calculation does not mutate the population or either mapping.

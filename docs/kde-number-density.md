# KDE number density and single-valued populations

With `normalize=False`, KDE output retains the sum of positive number-concentration
weights; a restricted query interval retains only the density supported there.
With `normalize=True`, the distribution helper normalizes a positive integral.
Normalization requires at least two strictly increasing query nodes with finite
spacing in the integration coordinate. A regular, unnormalized KDE may still be
evaluated at a single point or in arbitrary query order. Zero supported density
remains zero, including when all weights are zero or the evaluated tails underflow.
The INSA distribution additionally retains its existing peak-normalization rule;
an empty supported distribution stays zero.

A single-valued population is represented as a point mass deposited onto adjacent
query nodes in the declared linear or logarithmic integration coordinate. The
trapezoidal integral conserves its supported number concentration, including at
endpoints. This fallback requires at least two finite increasing query nodes;
points outside the interval contribute zero. It introduces no physical spread.

The diameter and INSA variables pass their outer bin edges as point-mass support.
A point between an outer edge and the nearest center is deposited at that center
with its count conserved under the same trapezoidal quadrature. This boundary
representation preserves count, but moves the first moment to the nearest node.
Points outside the outer edges still contribute zero. Generic helper calls retain
query-interval support unless `point_mass_bounds` is supplied.

Zero-area particles are excluded from the INSA log distribution, since ln(0) is
undefined. Absolute densities retain only the concentration of positive-area
particles; normalization applies to that positive-area subset. An all-inactive
population returns zeros for both normalized and absolute KDE output.
# Normalized distributions at extreme concentrations

Normalized KDEs use SciPy's unscaled probability density, and normalized point
masses start with unit weight. They are then normalized on the query grid as
before. Multiplying by total concentration and subsequently dividing it out
can underflow or overflow unnecessarily; raw number-density output still uses
the total concentration.

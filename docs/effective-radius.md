# Population effective radius

`ParticlePopulation.get_effective_radius()` uses wet radii in meters and returns
`sum(N*r**3) / sum(N*r**2)`, the surface-area-weighted effective radius.
This differs from the number-mean radius for polydisperse populations. With no
total particle area (including an empty or zero-concentration population), the
result is undefined and returned as NaN.

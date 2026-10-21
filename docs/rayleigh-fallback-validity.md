# Rayleigh fallback validity and unresolved acceptance policy

The opt-in fallback currently enforces only x = 2*pi*r/lambda < 0.1. This is
a necessary size restriction, not a universal error bound for arbitrary
complex refractive index m. Internal optical size |m|*x and refractive-index
contrast also matter. Custom high-index materials have not been validated.

At x=0.09 and wavelength 550 nm, Rayleigh/Mie absorption ratios measured with
the implemented formula and PyMieScatt 1.8.1.1 were about 0.9955 for
m=1.5+0.02i, 0.4734 for 10+0.1i, and 0.000994 for 30+0.1i. The current
constructor accepts these indices when the fallback is enabled. Tests at
ordinary aerosol-like indices do not justify a claim of uniform accuracy.

Until maintainers choose and validate an enforced refractive-index/internal-
size regime, use the Mie solver for high-contrast or otherwise unvalidated
materials, even if x passes the fallback check. A new threshold would change
the supported physical input domain and requires explicit review plus boundary
comparisons to Mie; no arbitrary threshold is introduced by this documentation.

PR #105 contains the fallback from #104, so any future numerical validity
guard must be retained when those overlapping changes are combined.

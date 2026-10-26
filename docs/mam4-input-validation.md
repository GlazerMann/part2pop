# MAM4 reconstruction inputs

Modes with exactly zero number concentration are omitted from initialization
and output reconstruction. All per-mode parameters are selected by the same
original mode indices. Unused composition in zero-number modes is ignored.
Every original number concentration is checked for missing, masked, nonfinite
or negative values before filtering. Invalid numbers are not inactive modes.

Active modes require finite nonnegative represented SO4/OC components and a
strictly positive dry total, including when wet sizes imply positive water.
The reader raises ValueError rather than inventing sulfate for an active mode
with no represented dry material. Files whose dry material exists only in
unrepresented native species need an extended reader, not a fallback mixture.

Initialization reads sulfate from `mfso4N` and organic material from `mfpomN`
in the namelist. Output reconstruction reads only `so4_aer` and `soa_aer`,
with optional water inferred from the diameter moments. Native BC, dust, salt
or other organic pools are not automatically represented by those two fields.
Positive number in the primary-carbon mode is therefore not sufficient to
determine compatibility: the represented dry fields must also have a positive
total. A mode containing only unrepresented material raises even if it is a
valid native-model state.
Conversely, a positive represented total does not prove complete species mass
closure: additional unrepresented material is not detected by this guard.
Validate the full native composition before treating a successful read as a
physical round trip.

This intentionally differs from replacing an empty represented mixture with
unit sulfate. Such replacement invents chemical composition and can affect
downstream optics, hygroscopicity and freezing. For an unsupported native
composition, stop reconstruction and explicitly extend the species mapping
with reviewed mass, number and size checks. Do not silently drop a mode with
positive particle number. Existing example files passing validation do not
establish support for every native composition.

Negative values and wet diameters smaller than dry diameters are rejected
without clipping or an invented measurement tolerance. This validation does
not resolve native record timing, unsupported species or the same-GSD wet/dry
moment assumption. The public timestep contract remains 1 for initialization
and n >= 2 for output record n - 2.

For active output modes, pressure (Pa), temperature (K) and derived dry-air
density must be finite positive real scalars. Required dry diameters must be
finite positive and GSD must exceed one. Wet reconstruction additionally
requires finite positive wet diameters at least as large as dry diameters and
finite nonnegative reconstructed water. Dry-only reconstruction does not use
wet-diameter values. Initial fractions and output masses are checked before
normalization; required masked values and complex pressure/temperature inputs
are rejected rather than converted to plausible-looking real values.

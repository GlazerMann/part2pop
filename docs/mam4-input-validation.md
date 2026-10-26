# MAM4 reconstruction inputs

Modes with exactly zero number concentration are omitted from initialization
and output reconstruction. All per-mode parameters are selected by the same
original mode indices. Unused composition in zero-number modes is ignored.

Active modes require finite nonnegative represented SO4/OC components and a
strictly positive dry total, including when wet sizes imply positive water.
The reader raises ValueError rather than inventing sulfate for an active mode
with no represented dry material. Files whose dry material exists only in
unrepresented native species need an extended reader, not a fallback mixture.

Negative values and wet diameters smaller than dry diameters are rejected
without clipping or an invented measurement tolerance. This validation does
not resolve native record timing, unsupported species or the same-GSD wet/dry
moment assumption. The public timestep contract remains 1 for initialization
and n >= 2 for output record n - 2.

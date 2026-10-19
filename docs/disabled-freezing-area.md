# Disabled freezing material and active area

Tabulated material contributes to active surface area only when both effective
`m_log10_Jhet` and `b_log10_Jhet` are non-NaN. Setting either coefficient to NaN
removes that material's volume from the equivalent-sphere area, matching its
exclusion from the heterogeneous-rate mixture. A zero/zero pair remains active.
Table coefficients are converted to floats for the check because the table
reader returns strings.

This standalone P12 patch applies directly to the pinned upstream main and
does not require P05. It retains main's existing treatment of custom materials;
P05 separately adds area support for complete custom coefficient pairs. If P05
is already applied, use the archive's after-p05.patch variant instead of the
standalone change.patch. The full series selects the after-P05 variant.

The equivalent-sphere morphology, freezing coefficients and water-activity
classification are unchanged. P01 separately handles rates when no positive
active volume is present. Both apply variants carry the same nine regressions.

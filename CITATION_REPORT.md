# Inline Citation Report

## Summary
- Citations to fix (in XML directives or `!!{ ... !!}` LaTeX docblocks): 9
- Affected files: 6
- Citations with NO BibTeX match (still need new entries): 0

Scope: only contexts that are rendered to the LaTeX documentation are
included — i.e. XML directives like `<description>`, `<longDescription>`,
`<defaultSource>`, `<defaultValue>`, etc. embedded inside `!![ ... !!]`
directive blocks, and free-text inside `!!{ ... !!}` LaTeX docstring
blocks. Inline `!`/`!!`/`!%` comments and runtime output strings (HDF5
attributes, command-line `description=` strings) are excluded because they
do not get rendered to LaTeX, so a `\cite{}` would not produce a real
citation.

## Citations to fix (grouped by file)

### source/dark_matter_profiles.SIDM.isothermal.F90
- Line 37: "Jiang et al. (2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [latex-docblock !!{...!!}]: `A dark matter halo profile class implementing profiles for self-interacting dark matter following the ``isothermal'' model of Jiang et al. (2022).`

### source/dark_matter_profiles_DMO.SIDM.coreNFW.F90
- Line 38: "Jiang et al. (2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [latex-docblock !!{...!!}]: `on the model of Jiang et al. (2022). The profile is defined by the enclosed mass, with (Jiang et al. 2022):`
- Line 81: "Jiang et al. (2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [xml:defaultSource]: `<defaultSource>Jiang et al. (2022)</defaultSource>`

### source/mass_distributions.spherical.SIDM.coreNFW.F90
- Line 28: "Jiang et al. (2022)" -> `\cite{jiang_semi-analytic_2023}`
  - Context [xml:description]: `on the model of Jiang et al. (2022). The profile is defined by the enclosed mass, with \citep{jiang_semi-analytic_2023}:`

### source/stellar_feedback.outflows.Creasey2013.F90
- Line 216: "Creasey et al. (2012)" -> `\cite{creasey_how_2013}`
  - Context [latex-docblock !!{...!!}]: `Integrand function for the ``Creasey et al. (2012)'' supernovae feedback calculation.`

### source/tests.mass_accretion_history.Hearin2021.F90
- Line 21: "Hearin (2021)" -> `\cite{hearin_differentiable_2021}`
  - Context [latex-docblock !!{...!!}]: `Contains a program which tests the Hearin (2021) halo mass formation history.`
- Line 26: "Hearin (2021)" -> `\cite{hearin_differentiable_2021}`
  - Context [latex-docblock !!{...!!}]: `Tests the Hearin (2021) halo mass formation history algorithm.`

### source/tests.mass_accretion_history.Hearin2021_stochastic.F90
- Line 21: "Hearin (2021)" -> `\cite{hearin_differentiable_2021}`
  - Context [latex-docblock !!{...!!}]: `Contains a program which tests the Hearin (2021) stochastic halo mass formation history.`
- Line 26: "Hearin (2021)" -> `\cite{hearin_differentiable_2021}`
  - Context [latex-docblock !!{...!!}]: `Tests the Hearin (2021) halo mass formation history algorithm.`

## Citations with NO BibTeX match (need new entries)

_None._ (The three apparent unmatched citations — `Bryan & Norman (1998)`,
`van der Wel et al. (2014)`, and `Leauthaud et al. (2012)` inside
`output.analyses.*.F90` `<description>` blocks — are not academic
citations: they appear inside `\mono{...}` as example string values that
users must enter into the `haloMassDefinition` and `reference` attributes
of input data files. These get rendered as literal labels on plots, not as
citations, so they should remain plain text.)

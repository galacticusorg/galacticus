#!/usr/bin/env python3
"""Demangle the Galacticus bibliography to convert accents (and other symbols) to LaTeX syntax."""
# Andrew Benson (27 February 2011); ported to Python

import re

with open('Galacticus.bib', 'r', encoding='utf-8', errors='replace') as ifile, \
     open('GalacticusAccented.bib', 'w', encoding='utf-8') as ofile:
    for line in ifile:

        # Special cases.
        # zwaan_hipass_2005 title
        line = line.replace(r'{\${\textbackslash}Omega\$_{{\textbackslash}rm} {HI}}', r'$\Omega_{\rm HI}$')

        # Obsoleted commands.
        line = re.sub(r'\{\\rm\s', r'\\mathrm{', line)

        # Unknown symbols (Unicode replacement character).
        line = line.replace('\ufffd', ' ')

        # Percentages.
        line = re.sub(r'(?<!\\)%', r'\\%', line)

        # Backslashes.
        line = re.sub(r'\{\\textbackslash\}([a-zA-Z]+)', r'{\\\1}', line)
        line = re.sub(r'\{\\textbackslash([a-zA-Z]+)\}', r'{\\\1}', line)
        line = line.replace(r'{\textbackslash}', '\\')

        # Tildes.
        line = line.replace(r'{\textasciitilde}', r'$\sim$')

        # Greek letters.
        line = line.replace('α', r'$\alpha$')
        line = line.replace('β', r'$\beta$')
        line = line.replace('γ', r'$\gamma$')
        line = line.replace('δ', r'$\delta$')
        line = line.replace('ɛ', r'$\epsilon$')
        line = line.replace('ϵ', r'$\epsilon$')
        line = line.replace('ε', r'$\epsilon$')
        line = line.replace('ζ', r'$\zeta$')
        line = line.replace('η', r'$\eta$')
        line = line.replace('θ', r'$\theta$')
        line = line.replace('ϑ', r'$\theta$')
        line = line.replace('ι', r'$\iota$')
        line = line.replace('κ', r'$\kappa$')
        line = line.replace('ϰ', r'$\kappa$')
        line = line.replace('λ', r'$\lambda$')
        line = line.replace('μ', r'$\mu$')
        line = line.replace('ν', r'$\nu$')
        line = line.replace('ξ', r'$\xi$')
        line = line.replace('ℴ', r'$o$')
        line = line.replace('π', r'$\pi$')
        line = line.replace('ϖ', r'$\pi$')
        line = line.replace('ρ', r'$\rho$')
        line = line.replace('ϱ', r'$\rho$')
        line = line.replace('σ', r'$\sigma$')
        line = line.replace('ς', r'$\sigma$')
        line = line.replace('τ', r'$\tau$')
        line = line.replace('υ', r'$\upsilon$')
        line = line.replace('ϕ', r'$\phi$')
        line = line.replace('φ', r'$\phi$')
        line = line.replace('χ', r'$\chi$')
        line = line.replace('ψ', r'$\psi$')
        line = line.replace('ω', r'$\omega$')
        line = line.replace('Γ', r'$\Gamma$')
        line = line.replace('Δ', r'$\Delta$')
        line = line.replace('Θ', r'$\Theta$')
        line = line.replace('Λ', r'$\Lambda$')
        line = line.replace('Ξ', r'$\Xi$')
        line = line.replace('Π', r'$\Pi$')
        line = line.replace('Σ', r'$\Sigma$')
        line = line.replace('Φ', r'$\Phi$')
        line = line.replace('Ψ', r'$\Psi$')
        line = line.replace('Ω', r'$\Omega$')

        # Astronomical symbols.
        line = line.replace('⊙', r'$\odot$')
        line = line.replace('☉', r'$\odot$')
        line = line.replace('★', r'$\star$')
        line = line.replace('⋆', r'$\star$')

        # Mathematical symbols.
        line = line.replace(r'{\textless}', r'$<$')
        line = line.replace(r'{\textgreater}', r'$>$')
        line = re.sub(r'\{\\textasciicircum\}([+\-\d.]+)', r'$^{\1}$', line)
        line = line.replace(r'{\textasciicircum}', '^')
        line = line.replace('×', r'$\times$')
        line = line.replace('±', r'$\pm$')
        line = line.replace('∓', r'$\mp$')
        line = line.replace('°', r'$^\circ$')
        line = line.replace('˜', r'$\sim$')
        line = line.replace('∼', r'$\sim$')
        line = line.replace('\u0303', r'$\sim$')
        line = line.replace('≈', r'$\approx$')
        line = line.replace('≃', r'$\approx$')
        line = line.replace('≡', r'$\equiv$', 1)
        line = line.replace('≳', r'$\gtrapprox$')
        line = line.replace('≲', r'$\lessapprox$')
        line = line.replace('√', r'$\sqrt{}$')
        line = line.replace('∗', r'$\star$')
        line = line.replace('Ṁ', r'$\dot{M}$')
        line = line.replace('≤', r'$\le$')
        line = line.replace('≥', r'$\ge$')
        line = line.replace('≠', r'$\ne$')
        line = line.replace('∝', r'$\propto$')
        line = line.replace('ℳ', r'$\mathcal{M}$')
        line = line.replace('ˆ', ' ')

        # Quotes.
        line = line.replace('\u2018', '`')
        line = line.replace('\u2019', "'")
        line = line.replace('\u201c', '``')
        line = line.replace('\u201d', "''")

        # Dashes.
        line = line.replace('\u2013', '--')
        line = line.replace('\u2212', '-')
        line = line.replace('\u2014', '-')

        # Accents.
        line = line.replace('Ç', r'\c{C}')
        line = line.replace('ü', r'\"u')
        line = line.replace('é', r"\'e")
        line = line.replace('â', r'\^a')
        line = line.replace('ä', r'\"a')
        line = line.replace('à', r'\`a')
        line = line.replace('ç', r'\c{c}')
        line = line.replace('ê', r'\^e')
        line = line.replace('ë', r'\"e')
        line = line.replace('è', r'\`e')
        line = line.replace('ï', r'\"i')
        line = line.replace('î', r'\^i')
        line = line.replace('ì', r'\`i')
        line = line.replace('Ä', r'\"A')
        line = line.replace('Å', r'\AA')
        line = line.replace('Á', r"\'A")
        line = line.replace('É', r"\'E")
        line = line.replace('ȯ', r'\.o')
        line = line.replace('ô', r'\^o')
        line = line.replace('ö', r'\"o')
        line = line.replace('ò', r'\`o')
        line = line.replace('û', r'\^u')
        line = line.replace('ù', r'\`u')
        line = line.replace('ÿ', r'\"y')
        line = line.replace('Ö', r'\"O')
        line = line.replace('Ü', r'\"U')
        line = line.replace('á', r"\'a")
        line = line.replace('í', r"\'i")
        line = line.replace('ó', r"\'o")
        line = line.replace('ú', r"\'u")
        line = line.replace('ñ', r'\~n')
        line = line.replace('Ñ', r'\~N')
        line = line.replace('ß', r'{\ss}')
        line = line.replace('Ã', r'\~A')
        line = line.replace('ø', r'\o{}')
        line = line.replace('š', r'\v{s}')
        line = line.replace('ć', r"\'c")
        line = line.replace('Ž', r'\v{Z}')
        line = line.replace('ł', r'\L{}')

        # Ligatures.
        line = line.replace('æ', 'ae')
        line = line.replace('Æ', 'Ae')
        line = line.replace('Œ', 'Oe')
        line = line.replace('ﬃ', 'ffi')
        line = line.replace('ﬀ', 'ff')

        ofile.write(line)

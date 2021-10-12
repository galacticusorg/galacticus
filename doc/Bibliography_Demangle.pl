#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";

# Demangle the Galacticus bibliography to convert accents (and other symbols) to LaTeX syntax.
# Andrew Benson (27 February 2011)

open(iHndl,"Galacticus.bib");
open(oHndl,">GalacticusAccented.bib");
while ( my $line = <iHndl> ) {

    # Special cases.
    # zwaan_hipass_2005 title
    $line =~ s/\{\\\$\{\\textbackslash\}Omega\\\$\\_\{\{\\textbackslash\}rm\} \{HI\}\}/\$\\Omega_{\\rm HI}\$/;

    # Unknown symbols.
    $line =~ s/�/ /g;
    $line =~ s/�/ /g;
    $line =~ s/�/ /g;

    # Percentages.
    $line =~ s/[^\\]%/\\%/g;

    # Backslashes.
    $line =~ s/{\\textbackslash}([a-zA-Z]+)/{\\$1}/g;
    $line =~ s/{\\textbackslash([a-zA-Z]+)}/{\\$1}/g;
    $line =~ s/{\\textbackslash}/\\/g;

    # Tildes.
    $line =~ s/{\\textasciitilde}/\$\\sim\$/g;

    # Greek letters.
    $line =~ s/α/\$\\alpha\$/g;
    $line =~ s/β/\$\\beta\$/g;
    $line =~ s/γ/\$\\gamma\$/g;
    $line =~ s/δ/\$\\delta\$/g;
    $line =~ s/ɛ/\$\\epsilon\$/g;
    $line =~ s/ϵ/\$\\epsilon\$/g;
    $line =~ s/ε/\$\\epsilon\$/g;
    $line =~ s/ζ/\$\\zeta\$/g;
    $line =~ s/η/\$\\eta\$/g;
    $line =~ s/θ/\$\\theta\$/g;
    $line =~ s/ϑ/\$\\theta\$/g;
    $line =~ s/ι/\$\\iota\$/g;
    $line =~ s/κ/\$\\kappa\$/g;
    $line =~ s/ϰ/\$\\kappa\$/g;
    $line =~ s/λ/\$\\lambda\$/g;
    $line =~ s/μ/\$\\mu\$/g;
    $line =~ s/ν/\$\\nu\$/g;
    $line =~ s/ξ/\$\\xi\$/g;
    $line =~ s/ℴ/\$\\omicron\$/g;
    $line =~ s/π/\$\\pi\$/g;
    $line =~ s/ϖ/\$\\pi\$/g;
    $line =~ s/ρ/\$\\rho\$/g;
    $line =~ s/ϱ/\$\\rho\$/g;
    $line =~ s/σ/\$\\sigma\$/g;
    $line =~ s/ς/\$\\sigma\$/g;
    $line =~ s/τ/\$\\tau\$/g;
    $line =~ s/υ/\$\\upsilon\$/g;
    $line =~ s/ϕ/\$\\phi\$/g;
    $line =~ s/φ/\$\\phi\$/g;
    $line =~ s/χ/\$\\chi\$/g;
    $line =~ s/ψ/\$\\psi\$/g;
    $line =~ s/ω/\$\\omega\$/g;
    $line =~ s/Γ/\$\\Gamma\$/g;
    $line =~ s/Δ/\$\\Delta\$/g;
    $line =~ s/Θ/\$\\Theta\$/g;
    $line =~ s/Λ/\$\\Lambda\$/g;
    $line =~ s/Ξ/\$\\Xi\$/g;
    $line =~ s/Π/\$\\Pi\$/g;
    $line =~ s/Σ/\$\\Sigma\$/g;
    $line =~ s/Φ/\$\\Phi\$/g;
    $line =~ s/Ψ/\$\\Psi\$/g;
    $line =~ s/Ω/\$\\Omega\$/g;
    
    # Astronomical symbols.
    $line =~ s/⊙/\$\\odot\$/g;
    $line =~ s/☉/\$\\odot\$/g;
    $line =~ s/★/\$\\star\$/g;
    $line =~ s/⋆/\$\\star\$/g;

    # Mathematical symbols.
    $line =~ s/{\\textless}/\$<\$/g;
    $line =~ s/{\\textgreater}/\$>\$/g;
    $line =~ s/{\\textasciicircum}([\+\-\d\.]+)/\$^{$1}\$/g;
    $line =~ s/{\\textasciicircum}/^/g;
    $line =~ s/×/\$\\times\$/g;
    $line =~ s/±/\$\\pm\$/g;
    $line =~ s/∓/\$\\mp\$/g;
    $line =~ s/°/\$^\\circ\$/g;
    $line =~ s/˜/\$\\sim\$/g;
    $line =~ s/∼/\$\\sim\$/g;
    $line =~ s/̃/\$\\sim\$/g;
    $line =~ s/≈/\$\\aproox\$/g;
    $line =~ s/≃/\$\\aproox\$/g;
    $line =~ s/≡/\$\\equiv\$/;
    $line =~ s/≳/\$\\gtrapprox\$/g;
    $line =~ s/≲/\$\\lessapprox\$/g;
    $line =~ s/√/\$\\sqrt{}\$/g;
    $line =~ s/∗/\$\\star\$/g;
    $line =~ s/Ṁ/\$\\dot{M}\$/g;
    $line =~ s/≤/\$\\ge\$/g;
    $line =~ s/≥/\$\\le\$/g;
    $line =~ s/≠/\$\\ne\$/g;
    $line =~ s/∝/\$\\propto\$/g;
    $line =~ s/ℳ/\$\\mathcal{M}\$/g;
    $line =~ s/ˆ/ /g;
     
    # Quotes.
    $line =~ s/‘/`/g;
    $line =~ s/’/'/g;
    $line =~ s/“/``/g;
    $line =~ s/”/''/g;
    
    # Dashes.
    $line =~ s/–/--/g;
    $line =~ s/−/-/g;
    $line =~ s/—/-/g;

    # Accents.
    $line =~ s/Ç/\\c{C}/g;
    $line =~ s/ü/\\"u/g;
    $line =~ s/é/\\'e/g;
    $line =~ s/â/\\^a/g;
    $line =~ s/ä/\\"a/g;
    $line =~ s/à/\\`a/g;
    $line =~ s/ç/\\c{c}/g;
    $line =~ s/ê/\\^e/g;
    $line =~ s/ë/\\"e/g;
    $line =~ s/è/\\`e/g;
    $line =~ s/ï/\\"i/g;
    $line =~ s/î/\\^i/g;
    $line =~ s/ì/\\`i/g;
    $line =~ s/Ä/\\"A/g;
    $line =~ s/Å/\\AA/g;
    $line =~ s/Á/\\'A/g;
    $line =~ s/É/\\'E/g;
    $line =~ s/ȯ/\\.o/g;
    $line =~ s/ô/\\^o/g;
    $line =~ s/ö/\\"o/g;
    $line =~ s/ò/\\`o/g;
    $line =~ s/û/\\^u/g;
    $line =~ s/ù/\\`u/g;
    $line =~ s/ÿ/\\"y/g;
    $line =~ s/Ö/\\"O/g;
    $line =~ s/Ü/\\"U/g;
    $line =~ s/á/\\'a/g;
    $line =~ s/í/\\'i/g;
    $line =~ s/ó/\\'o/g;
    $line =~ s/ú/\\'u/g;
    $line =~ s/ñ/\\~n/g;
    $line =~ s/Ñ/\\~N/g;
    $line =~ s/ß/{\\ss}/g;
    $line =~ s/Ã/\\~A/g;
    $line =~ s/ø/\\o{}/g;
    $line =~ s/š/\\v{s}/g;
    $line =~ s/ć/\\'c/g;
    $line =~ s/Ž/\\v{Z}/g;
    $line =~ s/ł/\\L{}/g;

    # Ligatures.
    $line =~ s/æ/ae/g;
    $line =~ s/Æ/Ae/g;
    $line =~ s/Œ/Oe/g;
    $line =~ s/ﬃ/ffi/g;
    $line =~ s/ﬀ/ff/g;

    # Write the modified line to file.
    print oHndl $line;
}
close(iHndl);
close(oHndl);

exit;

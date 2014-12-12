#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 

# Demangle the Galacticus bibliography to convert accents (and other symbols) to LaTeX syntax.
# Andrew Benson (27 February 2011)

open(iHndl,"Galacticus.bib");
open(oHndl,">GalacticusAccented.bib");
while ( my $line = <iHndl> ) {

    # Special cases.
    # zwaan_hipass_2005 title
    $line =~ s/{\\\${\\textbackslash}Omega\\\$\\_{{\\textbackslash}rm} {HI}}/\$\\Omega_{\\rm HI}\$/;

    # Percentages.
    $line =~ s/[^\\]%/\\%/g;

    # Backslashes.
    $line =~ s/{\\textbackslash}([a-zA-Z]+)/{\\$1}/g;
    $line =~ s/{\\textbackslash([a-zA-Z]+)}/{\\$1}/g;
    $line =~ s/{\\textbackslash}/\\/g;

    # Tildes.
    $line =~ s/{\\textasciitilde}/\$\\sim\$/g;

    # Astronomical symbols.
    $line =~ s/⊙/\$\\odot\$/g;

    # Less than symbol.
    $line =~ s/{\\textless}/\$<\$/g;

    # Greater than symbol.
    $line =~ s/{\\textgreater}/\$>\$/g;

    # Exponentiation.
    $line =~ s/{\\textasciicircum}([\+\-\d\.]+)/\$^{$1}\$/g;
    $line =~ s/{\\textasciicircum}/^/g;

    # Em dashes.
    $line =~ s/–/--/g;

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

    # Ligatures.
    $line =~ s/æ/ae/g;
    $line =~ s/Æ/Ae/g;

    # Write the modified line to file.
    print oHndl $line;
}
close(iHndl);
close(oHndl);

exit;

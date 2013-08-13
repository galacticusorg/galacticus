#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use List::Util qw(first);
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Plot errors from the Li & White (2009) stellar mass function along with a fitting function.
# Andrew Benson (10-December-2011)

# Read the XML data file.
my $xml       = new XML::Simple;
my $data      = $xml->XMLin($galacticusPath."data/observations/massFunctionsStellar/Stellar_Mass_Function_Li_White_2009.xml");
my $columns   = $data->{'massFunction'}->{'columns'};
my $x         = pdl @{$columns->{'mass'        }->{'datum'}};
my $y         = pdl @{$columns->{'massFunction'}->{'datum'}};
my $errorUp   = pdl @{$columns->{'upperError'  }->{'datum'}};
my $errorDown = pdl @{$columns->{'lowerError'  }->{'datum'}};
$errorUp      = 10.0**($y+$errorUp  );
$errorDown    = 10.0**($y-$errorDown);
$x            = 10.0** $x            ;
my $error     = 0.5*($errorUp-$errorDown);
my $errorFit  = 1.0e-3/(($x/4.5e10)**0.3)*exp(-$x/4.5e10)+1.0e-7;

# Declare variables for GnuPlot;
my $plotFile = "plots/stellarMassFunctionErrors_z01.pdf";
my ($gnuPlot, $plotFileEPS, $plot);

# Open a pipe to GnuPlot.
($plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.4,0.2\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e8:1.0e12]\n";
print $gnuPlot "set yrange [1.0e-8:1.0e-2]\n";
print $gnuPlot "set title 'Errors on the stellar mass function at \$z=0.1\$'\n";
print $gnuPlot "set xlabel '\$M_\\star\$ [\$M_\\odot\$]'\n";
print $gnuPlot "set ylabel '\$\\sigma_{{\\rm d}n/{\\rm d}\\log_{10}M_\\star}\$'\n";
&PrettyPlots::Prepare_Dataset(\$plot,
			      $x,$errorFit,
			      style     => "line",
			      weight    => [5,3],
			      color     => $PrettyPlots::colorPairs{'redYellow'},
			      title     => "Fit"
    );
&PrettyPlots::Prepare_Dataset(\$plot,
			      $x,$error,
			      style     => "point",
			      symbol    => [6,7], 
			      weight    => [5,3],
			      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
			      title     => $data->{'stellarMassFunction'}->{'label'}
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;

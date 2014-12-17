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
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
use XML::Simple;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Estimate completeness as a function of mass for the UKIDSS UDS survey of Caputi et al. (2011).
# Andrew Benson (28-April-2014)

# Number of sigma above sky for detection (arbitrary - will scale out of the results).
my $sigma = pdl 3.0;

# Argument in error function required for 80% completeness.
my $y80   = pdl sqrt(2.0)*erfi(2.0*(0.5-0.8));

# Specify mass completeness limits.
my @bins =
    (
     {
	 redshiftMinimum =>  3.00,
	 redshiftMaximum =>  3.50,
	 completeness50  => 10.30,
	 completeness80  => 10.93
     },
     {
	 redshiftMinimum =>  3.50,
	 redshiftMaximum =>  4.25,
	 completeness50  => 10.41,
	 completeness80  => 11.06
     },
     {
	 redshiftMinimum =>  4.25,
	 redshiftMaximum =>  5.00,
	 completeness50  => 10.51,
	 completeness80  => 11.18
     }
    );

# Load the existing data file.
my $xml = new XML::Simple;
my $observations = $xml->XMLin($galacticusPath."data/observations/massFunctionsStellar/Stellar_Mass_Functions_UKIDSS_UDS_2011.xml");

# Begin constructing the plot.
my $plot;
my $gnuPlot;
my $plotFile = $galacticusPath."constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/completeness.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set xlabel 'Galaxy stellar mass; \$M_\\star\\;[M_\\odot]\$\n";
print $gnuPlot "set ylabel 'Completeness; []'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale x\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e10:1.0e12]\n";
print $gnuPlot "set yrange [-0.05:1.05]\n";
print $gnuPlot "set pointsize 2.0\n";

# Iterate over redshift bins.
my $iBin = -1;
foreach my $bin ( @bins ) {
    ++$iBin;
    # Specify the 50% and 80% completeness masses from Caputi et al. (2011; Figure 4).
    my $m50 = pdl 10.0**$bin->{'completeness50'};
    my $m80 = pdl 10.0**$bin->{'completeness80'};
    # Solve for required parameters in the completeness model.
    my $flux = $y80**2*$m80/(($m50-$m80)**2-$y80**2*$m50/$sigma);
    my $sky  = ($flux*$m50/$sigma)**2;
    # Construct masses.
    my $logMass = pdl sequence(1000)*2.0/999.0+10.0;
    my $mass    = 10.0**$logMass;
    # Determine signal and noise as a function of mass.
    my $signal  = $flux*$mass;
    my $noise   = sqrt($sky);
    # Evaluate the completeness model.
    my $x            = ($sigma*$noise-$signal)/sqrt($noise**2+$signal);
    my $completeness = 0.5*(1.0-erf($x/sqrt(2.0)));
    # Construct label.
    my $label = sprintf("\$%4.2f<z<%4.2f\$",$bin->{'redshiftMinimum'},$bin->{'redshiftMaximum'});
    # Plot the completeness.
    &PrettyPlots::Prepare_Dataset(
	\$plot,
	$mass,
	$completeness,
	style      => "line",
	weight     => [5,3],
	title      => $label,
	color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[$iBin]}
	);
    # Compute completeness in each observed mass bin.
    my $observedMass         = pdl @{${$observations->{'massFunction'}->{'columns'}}[$iBin]->{'mass'}->{'datum'}};
    my $observedSignal       = $flux*10.0**$observedMass;
    my $observedX            = ($sigma*$noise-$observedSignal)/sqrt($noise**2+$observedSignal);
    my $observedCompleteness = 0.5*(1.0-erf($observedX/sqrt(2.0)));
    @{${$observations->{'massFunction'}->{'columns'}}[$iBin]->{'completeness'}->{'datum'}} = $observedCompleteness->list();
    ${$observations->{'massFunction'}->{'columns'}}[$iBin]->{'completeness'}->{'scaling'} = "linear";
    ${$observations->{'massFunction'}->{'columns'}}[$iBin]->{'completeness'}->{'description'} = "Completeness in this mass bin.";
}
# Finalize plot.
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

# Write the augmented data back to file.
open(my $outputFile,">".$galacticusPath."data/observations/massFunctionsStellar/Stellar_Mass_Functions_UKIDSS_UDS_2011.xml");
print $outputFile $xml->XMLout($observations,NoAttr=>1, RootName=>"dataSets");
close($outputFile);

exit;

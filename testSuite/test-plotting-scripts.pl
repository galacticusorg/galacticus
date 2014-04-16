#!/usr/bin/env perl
use strict;
use warnings;

# Run a set of short Galacticus models and make plots from them to test the plotting scripts.
# Andrew Benson (10-Oct-2012)

# First run the models.
system("cd ..; scripts/aux/launch.pl testSuite/parameters/test-plotting-scripts.xml");
system("cd ..; bunzip2 testSuite/outputs/test-plotting-scripts/galacticus_*/galacticus.hdf5.bz2");

# Run plotting commands.
my @plottingScripts =
    (
     {
	 script     => "Plot_Black_Hole_vs_Bulge_Mass.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "Plot_HI_Mass_Function.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "Plot_Star_Formation_History.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "stellarMassFunction_SDSS_z0.07.pl",
	 model      => "0:1",
	 constraint => 1
     },
     {
	 script     => "Plot_K_Luminosity_Function.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "Plot_Morphological_Luminosity_Function.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "Plot_SDSS_Color_Distribution.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "Plot_SDSS_Gas_Metallicity.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "Plot_Disk_Scalelengths.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "diskGalaxySizes_SDSS_z0.07.pl",
	 model      => "0:1",
	 constraint => 1
     },
     {
	 script     => "Plot_SDSS_Tully_Fisher.pl",
	 model      => "0:1",
	 constraint => 0
     },
     {
	 script     => "Plot_bJ_Luminosity_Function.pl",
	 model      => "0:1",
	 constraint => 0
     }
    );
foreach ( @plottingScripts ) {
    if ( $_->{'constraint'} == 0 ) {
	system("cd ..; scripts/plotting/".$_->{'script'}." testSuite/outputs/test-plotting-scripts/galacticus_".$_->{'model'}."/galacticus.hdf5 testSuite/outputs/test-plotting-scripts/galacticus_".$_->{'model'}." 1");
	print "FAILED: plotting script ".$_->{'script'}." failed\n"
	    unless ( $? == 0 );
    } else {
	system("cd ..; constraints/scripts/".$_->{'script'}." testSuite/outputs/test-plotting-scripts/galacticus_".$_->{'model'}."/galacticus.hdf5 --plotFile testSuite/outputs/test-plotting-scripts/galacticus_".$_->{'model'}."/tmp.pdf");
	print "FAILED: plotting script ".$_->{'script'}." failed\n"
	    unless ( $? == 0 );
    }
}

exit;

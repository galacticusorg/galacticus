#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use Text::Template qw(fill_in_string);
use PDL;
use PDL::IO::HDF5;
use Galacticus::Options;

# Script to generate content for halo mass function pipeline.
# Andrew Benson 8-April-2022

# Generates XML fragments for inclusion in haloMassFunctionConfig.xml
# Also generates haloMassFunctionBase_*.xml files for all completed halo mass functions (including inserting the measured environment properties)

# Get command line options.
my %options;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Validate required parameters are present.
die('simulationDataPath is required but is not present')
    unless ( exists($options{'simulationDataPath'}) );

my @types =
    (
     {
	 label        => "MilkyWay",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "CDM",
	 name         => "Milky Way",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984 ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_Axion20",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "Axion20",
	 name         => "Milky Way 10^-20 eV axion",
	 redshifts    => [ 0.000  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_Axion21",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "Axion21",
	 name         => "Milky Way 10^-21 eV axion",
	 redshifts    => [ 0.000  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_Axion22",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "Axion22",
	 name         => "Milky Way 10^-22 eV axion",
	 redshifts    => [ 0.000  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_WDM1",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM1",
	 name         => "Milky Way 1 keV WDM",
	 redshifts    => [ 0.000  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_WDM5",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM5",
	 name         => "Milky Way 5 keV WDM",
	 redshifts    => [ 0.000  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_WDM10",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM10",
	 name         => "Milky Way 10 keV WDM",
	 redshifts    => [ 0.000  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "LMC",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "CDM",
	 name         => "LMC",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984 ],
	 massParticle => 5.03538e04
     }
    );
open(my $configFile,">haloMassFunctionConfig.xml");
foreach my $type ( @types ) {
    opendir(my $dir,$options{'simulationDataPath'}."/ZoomIns/".$type->{'label'});
    my @halos = sort map {$_ =~ m/^Halo\d+/ ? $_ : ()} readdir($dir);
closedir($dir);    
    foreach my $halo ( @halos ) {
	$code::label       = $type->{'label'};
	$code::name        = $type->{'name'};
	$code::suite       = $type->{'suite'};
	$code::simulation  = $type->{'simulation'};
	$code::transfer    = $type->{'transfer'};
	$code::halo        = $halo;
	$code::massMinimum = sprintf("%11.5e",3000.0*$type->{'massParticle'});
	foreach my $redshift ( @{$type->{'redshifts'}} ) {
	    $code::redshift      = sprintf("%11.5e",$redshift);
	    $code::redshiftShort =  sprintf("%5.3f",$redshift);
	    # Make the XML for the config file.
	    my $config = fill_in_string(<<'CODE', PACKAGE => 'code');
    <!-- Zoom-in: {$name} {$halo} -->
    <parameterMap value="haloMassFunctionParameters::a             haloMassFunctionParameters::b
			 haloMassFunctionParameters::p             haloMassFunctionParameters::q 
			 haloMassFunctionParameters::normalization haloMassFunctionParameters::c
			 haloMassFunctionParameters::d
			 haloMassFunctionParameters::cW            haloMassFunctionParameters::beta"/>
    <parameterInactiveMap value=""/>
    <posteriorSampleLikelihood value="haloMassFunction">
      <!-- Options matched to those of Benson (2017; https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.3454B) -->
      <baseParametersFileName value="constraints/pipelines/darkMatter/haloMassFunctionBase_{$label}_{$halo}_z{$redshiftShort}.xml"/>
      <fileName               value="%DATASTATICPATH%/darkMatter/haloMassFunction_{$label}_{$halo}_z{$redshiftShort}.hdf5"/>
      <redshift               value="{$redshift}"/>
      <massRangeMinimum       value="{$massMinimum}"/> <!-- 3000 times zoom-in {$name} particle mass -->
      <binCountMinimum        value="0"          />    
    </posteriorSampleLikelihood>
CODE
	    print $configFile $config;
	    # Make the base file.
	    my $hmfFileName = $ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$type->{'label'}."_".$halo."_z".$code::redshiftShort.".hdf5";
	    if ( -e $hmfFileName ) {
		my $hmfFile = new PDL::IO::HDF5($hmfFileName);
		my $hmf = $hmfFile->group('simulation0001');
		(my $massRegion, my $overdensityRegion) = $hmf->attrGet('massRegion', 'overdensityRegion');
		$code::massRegion = sprintf("%+13.6e",$massRegion);
		$code::overdensityRegion = sprintf("%+9.6f    ",$overdensityRegion);
		my $base = fill_in_string(<<'CODE', PACKAGE => 'code');
<?xml version="1.0" encoding="UTF-8"?>
<parameters>
  <formatVersion>2</formatVersion>
  <version>0.9.4</version>

  <!-- Output control -->
  <outputFileName value="haloMassFunction_{$label}_{$halo}_z{$redshiftShort}.hdf5"/>
  <outputTimes value="list">
    <redshifts value="{$redshiftShort}"/>
  </outputTimes>  

  <!-- Halo environments -->
  <haloEnvironment value="fixed">
    <massEnvironment value="{$massRegion}"/>
    <overdensity     value="{$overdensityRegion}"/>
    <redshift        value="{$redshiftShort}" ignoreWarnings="true"/>
  </haloEnvironment>

  <!-- Include Milky Way cosmology and mass function parameters -->
  <xi:include href="haloMassFunctionParameters.xml"            xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="simulation_{$suite}.xml"                   xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="cosmology_{$suite}.xml"                    xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="haloMassFunction_{$suite}.xml"             xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="transferFunction_{$suite}_{$transfer}.xml" xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>

</parameters>
CODE
		open(my $baseFile,">constraints/pipelines/darkMatter/haloMassFunctionBase_".$type->{'label'}."_".$halo."_z".$code::redshiftShort.".xml");
		print $baseFile $base;
		close($baseFile);
	    }
	}
    }
}
close($configFile);

exit;

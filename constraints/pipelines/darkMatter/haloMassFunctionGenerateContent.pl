#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use Text::Template qw(fill_in_string);
use XML::Simple;
use PDL;
use PDL::IO::HDF5;
use Galacticus::Options;
use List::Util;

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
	 label        => "MilkyWay_WDM3",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM3",
	 name         => "Milky Way 3 keV WDM",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984 ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_WDM6.5",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM6.5",
	 name         => "Milky Way 6.5 keV WDM",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984 ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_WDM10",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM10",
	 name         => "Milky Way 10 keV WDM",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_IDM1GeV_envelope",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1GeV_envelope",
	 name         => "Milky Way IDM 1 GeV envelope",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_IDM1GeV_halfmode",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1GeV_halfmode",
	 name         => "Milky Way IDM 1 GeV half-mode",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_IDM1e-2GeV_envelope",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1e-2GeV_envelope",
	 name         => "Milky Way IDM 10^{-2} GeV envelope",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_IDM1e-2GeV_halfmode",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1e-2GeV_halfmode",
	 name         => "Milky Way IDM 10^{-2} GeV half-mode",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_IDM1e-4GeV_envelope",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1e-4GeV_envelope",
	 name         => "Milky Way IDM 10^{-4} GeV envelope",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_IDM1e-4GeV_halfmode",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1e-4GeV_halfmode",
	 name         => "Milky Way IDM 10^{-4} GeV half-mode",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_fdm_25.9e-22eV",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "fdm_25.9e-22eV",
	 name         => "Milky Way FDM 25.9*10^{-22} eV",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_fdm_185e-22eV",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "fdm_185e-22eV",
	 name         => "Milky Way FDM 185*10^{-22} eV",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_fdm_490e-22eV",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "fdm_490e-22eV",
	 name         => "Milky Way FDM 490*10^{-22} eV",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 4.02830e05
     },
     {
	 label        => "MilkyWay_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "CDM",
	 name         => "Milky Way",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984 ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_WDM3_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM3",
	 name         => "Milky Way 3 keV WDM",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984 ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_WDM6.5_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM6.5",
	 name         => "Milky Way 6.5 keV WDM",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984 ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_WDM10_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "WDM10",
	 name         => "Milky Way 10 keV WDM",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_IDM1GeV_envelope_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1GeV_envelope",
	 name         => "Milky Way IDM 1 GeV envelope",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_IDM1GeV_halfmode_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1GeV_halfmode",
	 name         => "Milky Way IDM 1 GeV half-mode",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_IDM1e-2GeV_envelope_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1e-2GeV_envelope",
	 name         => "Milky Way IDM 10^{-2} GeV envelope",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_IDM1e-2GeV_halfmode_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1e-2GeV_halfmode",
	 name         => "Milky Way IDM 10^{-2} GeV half-mode",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_IDM1e-4GeV_envelope_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1e-4GeV_envelope",
	 name         => "Milky Way IDM 10^{-4} GeV envelope",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_IDM1e-4GeV_halfmode_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "IDM1e-4GeV_halfmode",
	 name         => "Milky Way IDM 10^{-4} GeV half-mode",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_fdm_25.9e-22eV_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "fdm_25.9e-22eV",
	 name         => "Milky Way FDM 25.9*10^{-22} eV",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_fdm_185e-22eV_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "fdm_185e-22eV",
	 name         => "Milky Way FDM 185*10^{-22} eV",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "MilkyWay_fdm_490e-22eV_hires",
	 simulation   => "MilkyWay",
	 suite        => "Symphony",
	 transfer     => "fdm_490e-22eV",
	 name         => "Milky Way FDM 490*10^{-22} eV",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984  ],
	 massParticle => 5.03538e04
     },
     {
	 label        => "LMC",
	 simulation   => "LMC",
	 suite        => "Symphony",
	 transfer     => "CDM",
	 name         => "LMC",
	 redshifts    => [ 0.000, 0.504, 0.990, 2.031, 3.984 ],
	 massParticle => 5.03538e04
     }
    );
my @lengths = map {length($_->{'label'})} @types;
my $labelLengthMaximum = &List::Util::max(@lengths);
open(my $configFile  ,">haloMassFunctionConfig.xml");
my $pipelineBase   = "";
my $pipelineSuffix = "";
foreach my $type ( @types ) {
    die("Directory '".$options{'simulationDataPath'}."/ZoomIns/".$type->{'label'}."' does not exist")
	unless ( -e $options{'simulationDataPath'}."/ZoomIns/".$type->{'label'} );
    opendir(my $dir,$options{'simulationDataPath'}."/ZoomIns/".$type->{'label'});
    my @halos = sort map {$_ =~ m/^Halo\d+/ ? $_ : ()} readdir($dir);
    closedir($dir);    
    foreach my $halo ( @halos ) {
	$code::label       = $type->{'label'     }                              ;
	$code::name        = $type->{'name'      }                              ;
	$code::suite       = $type->{'suite'     }                              ;
	$code::simulation  = $type->{'simulation'}                              ;
	$code::hires       = $type->{'label'     } =~ m/_hires$/ ? "_hires" : "";
	$code::transfer    = $type->{'transfer'  }                              ;
	$code::halo        = $halo;
	$code::massMinimum = sprintf("%11.5e",3000.0*$type->{'massParticle'});
	foreach my $redshift ( @{$type->{'redshifts'}} ) {
	    $code::redshift      = sprintf("%11.5e",$redshift);
	    $code::redshiftShort = sprintf( "%5.3f",$redshift);
	    next
		unless ( -e $options{'simulationDataPath'}."/ZoomIns/".$type->{'label'}."/".$halo."/primaryHalo_z".$code::redshiftShort.".xml" );
	    # Extract the target halo properties.
	    my $xml = new XML::Simple();
	    my $haloTarget = $xml->XMLin($options{'simulationDataPath'}."/ZoomIns/".$type->{'label'}."/".$halo."/primaryHalo_z".$code::redshiftShort.".xml");
	    $code::massMaximum = sprintf("%11.5e",$haloTarget->{'mc'}/10.0);
	    # Add entries for the pipeline file.
	    $pipelineBase   .= '	      "haloMassFunctionBase_'.$type->{'label'}.'_'.$halo.'_z'.$code::redshiftShort.'.xml"'.(" " x (8+$labelLengthMaximum-length($type->{'label'})-length($halo))).",\n";
	    $pipelineSuffix .= '	      "'.$type->{'label'}.'_'.$halo.'_z'.$code::redshiftShort.'"'.(" " x (8+$labelLengthMaximum-length($type->{'label'})-length($halo))).",\n";
	    # Make the base file.
	    my $hmfFileName = $ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$type->{'label'}."_".$halo."_z".$code::redshiftShort.".hdf5";
	    if ( -e $hmfFileName ) {
		my $hmfFile = new PDL::IO::HDF5($hmfFileName);
		my $hmf     = $hmfFile->group('simulation0001');
		my $status  = 0;
		foreach my $attributeName ( 'massEnvironment', 'overdensityEnvironment' ) {
		    unless ( grep {$_ eq $attributeName} $hmf->attrs() ) {
			print "Attribute '".$attributeName."' is missing in file '".$hmfFileName."'\n";
			$status = 1;
		    }
		}
		if ( $status == 0 ) {
		    (my $massEnvironment, my $overdensityEnvironment) = $hmf->attrGet('massEnvironment', 'overdensityEnvironment');
		    $code::massEnvironment        = sprintf("%13.6e",       $massEnvironment);
		    $code::overdensityEnvironment = sprintf("%+9.6f",$overdensityEnvironment);
		    # Make the XML for the config file.
		    my $config = fill_in_string(<<'CODE', PACKAGE => 'code');
    <!-- Zoom-in: {$name} {$halo} -->
    <parameterMap value="haloMassFunctionParameters::a             haloMassFunctionParameters::p
                         haloMassFunctionParameters::normalization haloMassFunctionParameters::q
                         haloMassFunctionParameters::b
                         haloMassFunctionParameters::cW            haloMassFunctionParameters::beta
                         haloMassFunctionParameters::alpha         varianceFractionalModelDiscrepancy"/>
    <parameterInactiveMap value="" ignoreWarnings="true"/>
    <posteriorSampleLikelihood value="haloMassFunction">
      <!-- Options matched to those of Benson (2017; https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.3454B) -->
      <baseParametersFileName value="constraints/pipelines/darkMatter/haloMassFunctionBase_{$label}_{$halo}_z{$redshiftShort}.xml"/>
      <fileName               value="%DATASTATICPATH%/darkMatter/haloMassFunction_{$label}_{$halo}_z{$redshiftShort}.hdf5"/>
      <redshift               value="{$redshift}"       />
      <massRangeMinimum       value="{$massMinimum}"    /> <!-- 3000 times zoom-in {$name} particle mass -->
      <massRangeMaximum       value="{$massMaximum}"    /> <!-- 1/10 of the target halo mass             -->
      <binCountMinimum        value="0"                 />    
      <likelihoodPoisson      value="true"              />
      <likelihoodModel        value="simulatonSphere"   />
      <massSphere             value="{$massEnvironment}"/>
      <truncatePower          value="true"              />
    </posteriorSampleLikelihood>
CODE
		    print $configFile $config;
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
    <massEnvironment value="{$massEnvironment}"/>
    <overdensity     value="{$overdensityEnvironment}"/>
    <redshift        value="{$redshiftShort}" ignoreWarnings="true"/>
  </haloEnvironment>

  <!-- Include Milky Way cosmology and mass function parameters -->
  <xi:include href="haloMassFunctionParameters.xml"                xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="simulation_{$suite}_{$simulation}{$hires}.xml" xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="cosmology_{$suite}.xml"                        xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="haloMassFunction_{$suite}.xml"                 xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="transferFunction_{$suite}_{$transfer}.xml"     xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>

</parameters>
CODE
		    open(my $baseFile,">constraints/pipelines/darkMatter/haloMassFunctionBase_".$type->{'label'}."_".$halo."_z".$code::redshiftShort.".xml");
		    print $baseFile $base;
		    close($baseFile);
		}
	    }
	}
    }
}
close($configFile);

# Write base array for pipeline file.
open(my $pipelineBaseFile,">pipeline.base.pl");
print $pipelineBaseFile $pipelineBase;
close($pipelineBaseFile);

# Write suffix array for pipeline file.
open(my $pipelineSuffixFile,">pipeline.suffix.pl");
print $pipelineSuffixFile $pipelineSuffix;
close($pipelineSuffixFile);

exit;

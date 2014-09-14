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
unshift(@INC,$galacticusPath."perl"); 
use Scalar::Util 'reftype';
use XML::LibXML qw();
use Data::Dumper;
require Galacticus::Options;
require List::ExtraUtils;

# Update a Galacticus parameter file from one version to a later version.
# Andrew Benson (13-September-2014)

# Get arguments.
die("Usage: parametersMigrate.pl <inputFile> <outputFile> [options]")
    unless ( scalar(@ARGV) >= 2 );
my $inputFileName  = $ARGV[0];
my $outputFileName = $ARGV[1];
my %options =
    (
     inputVersion  => "0.9.2",
     outputVersion => "0.9.3"
    );

# Define translations.
my @translations =
    (
      {
	 inputVersion  => "0.9.0",
	 outputVersion => "0.9.1",
 	 names         =>
	 {
 	     Omega_0 => "Omega_Matter",
	 },
	 values        =>
	 {
	     adafEnergyOption                      =>
	     {
		 "pure ADAF"                => "pureADAF"               ,
	     },
	     coolingFunctionMethods                =>
	     {
		 "atomic_CIE_Cloudy"        => "atomicCIECloudy"        ,
		 "CIE_from_file"            => "cieFromFile"            ,
	     },
	     coolingSpecificAngularMomentumMethod  => 
	     {
		 "constant rotation"        => "constantRotation"       ,
	     },
	     coolingTimeAvailableMethod            => 
	     {
		 "White-Frenk"              => "White-Frenk1991"        ,
	     },
	     coolingRateMethod                     => 
	     {
		 "White + Frenk"            => "White-Frenk1991"        ,
	     },
	     cosmologyMethod                       =>
	     {
		 "matter + lambda"          => "matter-lambda"          ,
	     },
	     criticalOverdensityMethod             =>
	     {
		 "spherical top hat"        => "sphericalTopHat"        ,
	     },
	     darkMatterConcentrationMethod         =>
	     {
		 "Gao 2008"                 => "Gao2008"                ,
	     },
	     hotHaloCoolingFromNode                =>
	     {
		 "current node"             => "currentNode"            ,
		 "formation node"           => "formationNode"
	     },
	     hotHaloDensityMethod                  =>
	     {
		 "cored isothermal"         => "coredIsothermal"
	     },
	     infallRadiusMethod                    =>
	     {
		 "cooling and freefall"     => "coolingAndFreefall"     ,
	     },
	     nodeMergersMethod                     => 
	     {
		 "single level hierarchy"   => "singleLevelHierarchy"   ,
	     },
	     powerSpectrumMethod                   => 
	     {
		 "power law"                => "powerLaw"               ,
	     },
	     satelliteMergingMethod                => 
	     {
		 "Lacey-Cole + Tormen"      => "Lacey-Cole+Tormen"      ,
	     },
	     starFormationHistoriesMethod          => 
	     {
		 "metallicity split"        => "metallicitySplit"       ,
	     },
	     starFormationHistoriesMethod          => 
	     {
		 "metallicity split"        => "metallicitySplit"       ,
	     },
	     starFormationFeedbackDisksMethod      => 
	     {
		 "power law"                => "powerLaw"               ,
	     },
	     starFormationFeedbackSpheroidsMethod  => 
	     {
		 "power law"                => "powerLaw"               ,
	     },
	     starFormationTimescaleDisksMethod     => 
	     {
		 "dynamical time"           => "dynamicalTime"          ,
	     },
	     starFormationTimescaleSpheroidsMethod => 
	     {
		 "dynamical time"           => "dynamicalTime"          ,
	     },
	     transferFunctionMethod                => 
	     {
		 "Eisenstein + Hu"          => "Eisenstein-Hu1999"      ,
	     },
	     treeBranchingMethod                   => 
	     {
		 "modified Press-Schechter" => "modifiedPress-Schechter",
	     },
	     virialDensityContrastMethod           =>
	     {
		 "spherical top hat"        => "sphericalTopHat"        ,
	     },
	 },
      },
      {
	 inputVersion  => "0.9.1",
	 outputVersion => "0.9.2",
 	 names         =>
	 {
 	     treeNodeMethodSatelliteOrbit                => "treeNodeMethodSatellite",
	 },
	 values        =>
	 {
 	     treeNodeMethodSatellite              =>
	     {
		 "simple"                    => "standard"
	     },
	     treeNodeMethodSpheroid               =>
	     {
		 "Hernquist"                 => {
		                                 value => "standard",
						 new   => [
						           {
							    name  => "spheroidMassDistribution",
							    value => "hernquist"
						           }
                                                          ]
		                                },
		 "sersic"                    => {
		                                 value => "standard",
						 new   => [
						           {
							    name  => "spheroidMassDistribution",
							    value => "sersic"
						           }
                                                          ]
		                                }
	     },
	 }
     },
     {
	 inputVersion  => "0.9.2",
	 outputVersion => "0.9.3",
	 names         =>
	 {
	     accretionHalosMethod                        => "accretionHaloMethod"                 ,
	     cosmologyMethod                             => "cosmologyFunctionsMethod"            ,
	     darkMatterConcentrationMethod               => "darkMatterProfileConcentrationMethod",
	     hotHaloCoredIsothermalCoreRadiiMethod       => "hotHaloColdModeCoredIsothermalCoreRadiiMethod",
	     hotHaloDensityMethod                        => "hotHaloMassDistributionMethod"       ,
	     ionizationStateFile                         => "chemicalStateFile",
	     isothermalCoreRadiusOverScaleRadius         => "hotHaloCoreRadiusOverScaleRadius",
	     isothermalCoreRadiusOverVirialRadius        => "hotHaloCoreRadiusOverVirialRadius"   ,
	     isothermalCoreRadiusOverVirialRadiusMaximum => "hotHaloCoreRadiusOverVirialRadiusMaximum",
	     mergerTreeBuildCole2000MassResolution       => "mergerTreeHaloMassResolution",
	     nfw96ConcentrationC                         => "nfw1996ConcentrationC",
	     satelliteMergingMethod                      => "satelliteMergingTimescalesMethod"    ,
	     treeNodeMethodFormationTimes                => "treeNodeMethodFormationTime",
	     virialOrbitsMethod                          => "virialOrbitMethod"
	 },
	 values        =>
         {
	     cosmologyFunctionsMethod             =>
	     {
		 "matter-lambda"             => "matterLambda"                 ,
		 "matter-darkEnergy"         => "matterDarkEnergy"             ,
	     },
	     darkMatterProfileConcentrationMethod =>
	     {
		 "Gao2008"                   => "gao2008"                      ,
		 "Munoz-Cuartas2011"         => "munozCuartas2011"             ,
		 "Prada2011"                 => "prada2011"                    ,
		 "Zhao2009"                  => "zhao2009"
	     },
	     hotHaloMassDistributionMethod        =>
	     {
		 "coredIsothermal"           => "betaProfile"
	     },
	     satelliteMergingTimescalesMethod     => 
	     {
		 "BoylanKolchin2008"         => "boylanKolchin2008"            ,
		 "Jiang2008"                 => "jiang2008"                    ,
		 "Lacey-Cole"                => "laceyCole1993"                ,
		 "Lacey-Cole+Tormen"         => "laceyCole1993Tormen"          ,
		 "Wetzel-White2010"          => "wetzelWhite2010"
	     },
	     virialDensityContrastMethod          =>
	     {
		 "Bryan-Norman1998"          => "bryanNorman1998"              ,
		 "sphericalTopHatDarkEnergy" => "sphericalCollapseMatterDE"    ,
		 "sphericalTopHat"           => "sphericalCollapseMatterLambda",
		 "Kitayama-Suto1996"         => "kitayamaSuto1996"
	     },
	     virialOrbitMethod                    =>
	     {
		 "Benson2005"                => "benson2005"                   ,
		 "Wetzel2010"                => "wetzel2010"
 
	     },
	 }
     }
    );

# Validate the paramter file.
system($galacticusPath."scripts/aux/validateParameters.pl ".$inputFileName);
die('input file is not a valid Galacticus parameter file')
    unless ( $? == 0 );

# Parse the input file.
my $parser     = XML::LibXML->new();
my $input      = $parser->parse_file($inputFileName);
my $parameters = $input->findnodes('parameters')->[0];

# Check for a version element.
my $version    = $parameters->findnodes('version')->[0];
if ( defined($version) ) {
    $options{'inputVersion'} = $version->textContent();
    $version->firstChild()->setData("0.9.3");
} else {
    my $versionNode = $input->createElement ("version");
    my $newVersion  = $input->createTextNode("0.9.3"  );
    my $newBreak    = $input->createTextNode("\n  "   );
    $versionNode->addChild($newVersion);
    $parameters->insertBefore($versionNode,$parameters->firstChild());
    $parameters->insertBefore($newBreak   ,$parameters->firstChild());
}

# Parse options.
&Options::Parse_Options(\@ARGV,\%options);

# Iterate over translations.
foreach my $translation ( @translations ) {
    # Skip irrelevant translation.
    next
	unless ( $translation->{'inputVersion'} eq $options{'inputVersion'} );
    # Report.
    print "Translating from v".$translation->{'inputVersion'}." to v".$translation->{'outputVersion'}."\n";
    $options{'inputVersion'} = $translation->{'outputVersion'};
    # Apply translations.
    for my $parameter ( $parameters->findnodes('parameter') ) {
	# Get name and value text elements.
	my $name  = $parameter->findnodes('name' )->[0]->firstChild();
	my $value = $parameter->findnodes('value')->[0]->firstChild();
	# Translate names.
	if ( exists($translation->{'names'}->{$name->data()}) ) {
	    print "   translate parameter name: ".$name->data()." --> ".$translation->{'names'}->{$name->data()}."\n";
	    $name->setData($translation->{'names'}->{$name->data()});;
	}
	# Translate values.
	if ( exists($translation->{'values'}->{$name->data()}) ) {
	    # Split values.
	    my $valuesText = $value->data();
	    $valuesText =~ s/^\s*//;
	    $valuesText =~ s/\s*$//;
	    my @values;
	    if ( $translation->{'inputVersion'} eq "0.9.0" ) {
		# For v0.9.0, we cannot split as method names were permitted to contain spaces.
		push(@values,$valuesText);
	    } else {
		@values = split(/\s+/,$valuesText);
	    }
	    foreach my $thisValue ( @values ) {
		if ( exists($translation->{'values'}->{$name->data()}->{$thisValue}) ) {
		    print "   translate parameter value: ".$name->data()."\n";
		    if ( ref($translation->{'values'}->{$name->data()}->{$thisValue}) ) {
			my $newValue = $translation->{'values'}->{$name->data()}->{$thisValue};
			print "                                 ".$thisValue." --> ".$newValue->{'value'}."\n";
			$thisValue = $newValue->{'value'};
			if ( exists($newValue->{'new'}) ) {
			    foreach my $newParameter ( &ExtraUtils::as_array($newValue->{'new'}) ) {
				print "      add parameter: ".$newParameter->{'name'}." = ".$newParameter->{'value'}."\n";
				my $parameterNode = $input->createElement (                "parameter" );
				my $name          = $input->createElement (                "name"      );
				my $value         = $input->createElement (                "value"     );
				my $nameText      = $input->createTextNode($newParameter->{'name'     });
				my $valueText     = $input->createTextNode($newParameter->{'value'    });
				$name ->addChild($nameText );
				$value->addChild($valueText);
				$parameterNode->addChild($input->createTextNode("\n    "));
				$parameterNode->addChild($name  );
				$parameterNode->addChild($input->createTextNode("\n    "));
				$parameterNode->addChild($value );
				$parameterNode->addChild($input->createTextNode("\n  "));
				$parameters->insertAfter($parameterNode,$parameter);
				$parameters->insertAfter($input->createTextNode("\n  "),$parameter);
			    }
			}
		    } else {
			print "                                 ".$thisValue." --> ".$translation->{'values'}->{$name->data()}->{$thisValue}."\n";
			$thisValue = $translation->{'values'}->{$name->data()}->{$thisValue};
		    }
		}
		$value->setData(join(" ",@values));
	    }
	}
    }
    # Finished if output version is reached.
    last
	if ( $options{'outputVersion'} eq $options{'inputVersion'} );
}

# Output the resulting file.
$input->toFile($outputFileName);

exit;

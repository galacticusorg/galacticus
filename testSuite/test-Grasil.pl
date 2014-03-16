#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "../";
 $ENV{"GALACTICUS_ROOT_V093"} = getcwd()."/../";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
require Galacticus::HDF5;
require Galacticus::SubMmFluxesHayward;
require Galacticus::Grasil;
require Galacticus::Survey;

# Run simple Galacticus+Grasil calculations.
# Andrew Benson (12-January-2013)

# First run the models.
system("rm -rf ../aux/Grasil");
system("cd ..; scripts/aux/launch.pl testSuite/parameters/test-Grasil.xml");
system("cd ..; bunzip2 testSuite/outputs/test-Grasil/galacticus_*/galacticus.hdf5.bz2");

# Open the model output.
my $model = new PDL::IO::HDF5("outputs/test-Grasil/galacticus_0:1/galacticus.hdf5");

# Iterate over output times.
my $failuresFound = 0;
my @outputs = $model->group("Outputs")->groups();
foreach my $output ( @outputs ) {
    # Get the output time,
    (my $outputTime) = $model->group("Outputs")->group($output)->attrGet('outputTime');
    # Get the star formation history group for this output.
    my $outputGroup = $model->group("starFormationHistories")->group($output);
    # Get all merger trees in the group.
    my @mergerTrees = $outputGroup->groups();
    # Iterate over merger trees.
    foreach my $mergerTree ( @mergerTrees ) {
	# Get all datasets in the group.
	my @datasets = $outputGroup->group($mergerTree)->datasets();
	# Iterate over datasets.
	foreach my $dataset ( @datasets ) {
	    # Find time datasets.
	    if ( $dataset =~ m /Time\d+$/ ) {
		# Read the list of times.
		my $times = $outputGroup->group($mergerTree)->dataset($dataset)->get();
		# Ensure that each output star formation history extends to the output time.
		if ( $times((-1)) < $outputTime && $failuresFound == 0 ) {
		    print "FAIL: star formation history time series does not extend to output time\n";
		    $failuresFound = 1;
		}
	    }
	}
    }
}
undef($model);
exit
    if ( $failuresFound == 1 );

# Specify minimum flux and maximum redshift.
my $maximumRedshift = 1.5;
my $minimumFlux     = 1.0e-3;

# Get available output times for this model.
my $galacticus;
$galacticus->{'file'} = "outputs/test-Grasil/galacticus_0:1/galacticus.hdf5";
&HDF5::Get_Times ($galacticus);

# Loop over all available redshifts.
my $outputNumber = 0;
foreach my $outputRedshift ( $galacticus->{'outputs'}->{'redshift'}->list() ) {
    ++$outputNumber;

    # Get the star formation history group for this output.
    my $outputGroup = $galacticus->{'hdf5File'}->group("starFormationHistories")->group("Output".$outputNumber);

    # Find all trees at this output.
    my @treeGroups = $outputGroup->groups();

    # Specify model.
    my $galacticus;
    $galacticus->{'file' } = "outputs/test-Grasil/galacticus_0:1/galacticus.hdf5";
    $galacticus->{'store'} = 0;
    $galacticus->{'tree' } = "all";
  
    # Specify Grasil options.
    $galacticus->{'grasilOptions'}->{'includePAHs'            } = 0;
    $galacticus->{'grasilOptions'}->{'fluctuatingTemperatures'} = 0;
    $galacticus->{'grasilOptions'}->{'wavelengthCount'        } = 100;
    $galacticus->{'grasilOptions'}->{'radialGridCount'        } = 10;
    $galacticus->{'grasilOptions'}->{'cpuLimit'               } = 600;
    
    # Set parameters of the Hayward et al. fitting formula for 850um flux to those which best fit Grasil results.
    $galacticus->{'haywardSubMmFit'}->{'dustToMetalsRatio'        } = pdl 0.61;
    $galacticus->{'haywardSubMmFit'}->{'fitNormalization'         } = pdl 1.00e-3;
    $galacticus->{'haywardSubMmFit'}->{'starFormationRateExponent'} = pdl 0.42;
    $galacticus->{'haywardSubMmFit'}->{'dustmassExponent'         } = pdl 0.68;

    # Read results from model.
    &HDF5::Get_Parameters        ($galacticus                );
    &HDF5::Select_Output         ($galacticus,$outputRedshift);
    &HDF5::Get_Datasets_Available($galacticus                );
    if ( exists($galacticus->{'dataSetsAvailable'}->{'nodeIndex'}) ) {
	&HDF5::Get_Dataset($galacticus,
			   [
			    'nodeIndex'           ,
			    'mergerTreeIndex'     ,
			    'redshift'            ,
			    'flux850micronHayward',
			    'diskMassStellar'     ,
			    'spheroidMassStellar'
			   ]
	    );
	my $dataSets = $galacticus->{'dataSets'};
	$galacticus->{'selection'} = 
	    which(
		( $dataSets->{'flux850micronHayward'} > $minimumFlux     )
		& 
		( $dataSets->{'redshift'            } < $maximumRedshift )
	    );


	
	# Extract recycling parameter for the model.
	die("test-Grasil.pl: FAIL - fixed IMF and instantaneous recycling expected for this model")
	    unless ( $galacticus->{'parameters'}->{'imfSelectionMethod'} eq "fixed" && $galacticus->{'parameters'}->{'stellarPopulationPropertiesMethod'} eq "instantaneous" );
	my $imfRecycledAttributeName = "imf".$galacticus->{'parameters'}->{'imfSelectionFixed'}."RecycledInstantaneous";
	my $recycledFraction =$galacticus->{'parameters'}->{$imfRecycledAttributeName};

	# Determine the number of galaxies selected and write a message reporting this.
	my $numberSelected = nelem($galacticus->{'selection'});
	print "At z=".$outputRedshift.", ".$numberSelected." galaxies were selected.\n";

	# Process the galaxies through Grasil.
	&HDF5::Get_Dataset($galacticus,
	  		   [
	  		    'grasilFlux850microns',
	  		    'grasilFlux250microns',
	 		    'grasilFlux350microns',
	  		    'grasilFlux500microns',
	  		    'grasilInfraredLuminosity'
	  		   ]
	    );

	# Write out the results.
	for(my $i=0;$i<nelem($galacticus->{'selection'});++$i) {
	    my $j = $galacticus->{'selection'}->(($i));

	    # Assume zero mass for disk and spheroid by default.
	    my $diskMass     = pdl 0.0;
	    my $spheroidMass = pdl 0.0;
	    # Find the corresponding star formation histories.
	    my $treeIndex = $galacticus->{'dataSets'}->{'mergerTreeIndex'}->(($j));
	    my $nodeIndex = $galacticus->{'dataSets'}->{'nodeIndex'      }->(($j));
	    if ( grep $_ eq "mergerTree".$treeIndex, @treeGroups ) {
		# Find all available datasets.
		my $treeGroup = $outputGroup->group("mergerTree".$treeIndex);
		my @dataSets = $treeGroup->datasets();
		if ( grep $_ eq "diskSFH".$nodeIndex, @dataSets ) {
		    my $sfh    = $treeGroup->dataset("diskSFH".$nodeIndex)->get();
		    $diskMass = sum($sfh)*(1.0-$recycledFraction);

		}
		if ( grep $_ eq "spheroidSFH".$nodeIndex, @dataSets ) {
		    my $sfh    = $treeGroup->dataset("spheroidSFH".$nodeIndex)->get();
		    $spheroidMass = sum($sfh)*(1.0-$recycledFraction);

		}
	    }
	    # Check that the output star formation history integrates to the expected stellar mass.
	    if ( $diskMass == 0.0 ) {
		print "test-Grasil.pl: FAIL - zero mass in star formation history for non-zero mass disk\n"
		    unless ( $galacticus->{'dataSets'}->{'diskMassStellar'}->(($j)) == 0.0 );
	    } else {
		my $fractionalError = abs(($galacticus->{'dataSets'}->{'diskMassStellar'}->(($j))-$diskMass)/$diskMass);
		print "test-Grasil.pl: FAIL - star formation history does not integrate to expected stellar mass for disk [".$fractionalError."]\n"
		    unless ( $fractionalError < 1.0e-2 );
	    }
	    if ( $spheroidMass == 0.0 ) {
		print "test-Grasil.pl: FAIL - zero mass in star formation history for non-zero mass spheroid\n"
		    unless ( $galacticus->{'dataSets'}->{'spheroidMassStellar'}->(($j)) == 0.0 );
	    } else {
		my $fractionalError = abs(($galacticus->{'dataSets'}->{'spheroidMassStellar'}->(($j))-$spheroidMass)/$spheroidMass);
		print "test-Grasil.pl: FAIL - star formation history does not integrate to expected stellar mass for spheroid [".$fractionalError."]\n"
		    unless ( $fractionalError < 1.0e-2 );
	    }

	    foreach my $dataSet ( "redshift", "flux850micronHayward", "grasilFlux850microns", "grasilFlux250microns", "grasilFlux350microns", "grasilFlux500microns", "grasilInfraredLuminosity" ) {
		print $galacticus->{'dataSets'}->{$dataSet}->(($j))."\t";
	    }
	    print "\n";
	}

    }
    
}

exit;

# Contains a Perl module which implements calculations of SEDs using Grasil.

package Galacticus::Grasil;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::GSL::INTERP;
use Text::Table;
use Sys::CPU;
use Data::Dumper;
use List::Util qw(first);
use Galacticus::HDF5;
use Galacticus::Inclination;
use Galacticus::Filters;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::Local;
use Galacticus::Launch::PBS;
use Galacticus::Path;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
			       "^grasilFlux[\\d\\.]+microns\$"        => \&Galacticus::Grasil::Get_Flux                 ,
			       "^grasilFlux:([^:]+)\$"                => \&Galacticus::Grasil::Compute_Flux_Under_Filter,
			       "^grasilInfraredLuminosity\$"          => \&Galacticus::Grasil::Get_Flux                 ,
			       "^luminosity:grasil:([^:]+):([^:]+)\$" => \&Galacticus::Grasil::Get_Flux
    );

# Compute conversion factor for units from luminosity to flux.
my $ergs                = pdl 1.0e-7;
my $Angstrom            = pdl 1.0e-10;
my $Jansky              = pdl 1.0e-26;
my $speedOfLight        = pdl 2.998e8;
my $megaParsec          = pdl 3.086e22;
my $Pi                  = pdl 3.1415927;
my $conversionFactor    = pdl 1.0e30*$ergs*$Angstrom/$speedOfLight/$megaParsec**2/$Jansky;
my $solarLuminosity     = pdl 3.845e3; # Solar luminosity in units of 10^30 ergs/s.
my $luminosityAB        = pdl 4.46592015000e+13; # W/Hz
my $grasilNormalization = pdl 1.00000000000e+23; # W/Angstrom

# Define the range used for total IR luminosity.
my $irWavelengthMinimum =    8.0e4; # Angstroms.
my $irWavelengthMaximum = 1000.0e4; # Angstroms.

sub Get_Flux {
    # Get the flux at a specified wavelength.
    my $dataSet     = shift;
    my $dataSetName = $_[0];

    # Extract a basic property so we can figure out how many galaxies exist.
    &Galacticus::HDF5::Get_Dataset($dataSet,["nodeIndex","mergerTreeIndex","redshift","inclination","luminosityDistance"]);
    
    # Get a reference to the datasets.
    my $dataSets = $dataSet->{'dataSets'};

    # Determine for which selection of galaxies we will compute fluxes.
    my $selection;
    if ( exists($dataSet->{'selection'}) ) {
   	$selection = $dataSet->{'selection'};
    } else {
    	$selection = sequence(nelem($dataSets->{'nodeIndex'}));
    }
    
    # Create an empty version of the flux property.
    $dataSets->{$dataSetName} = pdl zeroes(nelem($dataSets->{'nodeIndex'}));
    
    # Extract the output index.
    my $outputIndex = $dataSet->{'output'};
    
    # Determine Grasil parameter options.
    my $includePAHs             = "1.0";
    my $fluctuatingTemperatures = "1.0";
    my $wavelengthCount         = 400;
    my $radialGridCount         = 40;
    my $dustToMetalsRatio       = 0.61;
    my $recomputeSEDs           = 0;
    my $cpuLimit                = 3600;
    if ( exists($dataSet->{'grasilOptions'}->{'includePAHs'}) ) {
    	if ( $dataSet->{'grasilOptions'}->{'includePAHs'} == 1 ) {
    	    $includePAHs = "1.0";
    	} else {
    	    $includePAHs = "0.0";
    	}
    }
    if ( exists($dataSet->{'grasilOptions'}->{'fluctuatingTemperatures'}) ) {
    	if ( $dataSet->{'grasilOptions'}->{'fluctuatingTemperatures'} == 1 ) {
    	    $fluctuatingTemperatures = "1.0";
    	} else {
    	    $fluctuatingTemperatures = "0.0";
    	}
    }
    $wavelengthCount   = $dataSet->{'grasilOptions'}->{'wavelengthCount'  } if ( exists($dataSet->{'grasilOptions'}->{'wavelengthCount'  }) );
    $radialGridCount   = $dataSet->{'grasilOptions'}->{'radialGridCount'  } if ( exists($dataSet->{'grasilOptions'}->{'radialGridCount'  }) );
    $dustToMetalsRatio = $dataSet->{'grasilOptions'}->{'dustToMetalsRatio'} if ( exists($dataSet->{'grasilOptions'}->{'dustToMetalsRatio'}) );
    $recomputeSEDs     = $dataSet->{'grasilOptions'}->{'recomputeSEDs'    } if ( exists($dataSet->{'grasilOptions'}->{'recomputeSEDs'    }) );
    $cpuLimit          = $dataSet->{'grasilOptions'}->{'cpuLimit'         } if ( exists($dataSet->{'grasilOptions'}->{'cpuLimit'         }) );

    # Determine the number of CPUs available.
    my $launchMethod = exists($dataSet->{'grasilOptions'}->{'launchMethod'}) ? $dataSet->{'grasilOptions'}->{'launchMethod'} : "local";
    my $grasilJobsMaximum = 100;
    $grasilJobsMaximum = Sys::CPU::cpu_count()
	if ( $launchMethod eq "local" );
    $grasilJobsMaximum = $dataSet->{'grasilOptions'}->{'grasilJobsMaximum'} 
        if ( exists($dataSet->{'grasilOptions'}->{'grasilJobsMaximum'}) );
    
    # Open the file for reading.
    &Galacticus::HDF5::Open_File($dataSet);
    my @groupsList;
    my @datasetsList;

    # Create groups as necessary.
    my $grasilSEDsGroup;
    @groupsList = $dataSet->{'hdf5File'}->groups();
    if ( defined(first { $_ eq "grasilSEDs" } @groupsList) ) {
    	$grasilSEDsGroup = $dataSet->{'hdf5File'}->group("grasilSEDs");
    } else {
    	$grasilSEDsGroup = new PDL::IO::HDF5::Group( name => "grasilSEDs", parent => $dataSet->{'hdf5File'},
    						     fileObj => $dataSet->{'hdf5File'} );
    }
    @groupsList = $dataSet->{'hdf5File'}->groups();
    unless ( defined(first { $_ eq "grasilSEDs" } @groupsList) ) {
	print "Failed to open grasilSEDs in ".$dataSet->{'file'}."\n";
	return;
    }

    my $outputGroup;
    @groupsList = $grasilSEDsGroup->groups();
    if ( defined(first { $_ eq "Output".$outputIndex } @groupsList) ) {
    	$outputGroup = $grasilSEDsGroup->group("Output".$outputIndex);
    } else {
    	$outputGroup = new PDL::IO::HDF5::Group( name => "Output".$outputIndex, parent => $grasilSEDsGroup,
    						 fileObj => $dataSet->{'hdf5File'} );
    }
    @groupsList = $grasilSEDsGroup->groups();
    unless ( defined(first { $_ eq "Output".$outputIndex } @groupsList) ) {
	print "Failed to open grasilSEDs/Output".$outputIndex." in ".$dataSet->{'file'}."\n";
	return;
    }

    # Initialize a queue of galaxies to process through Grasil.
    my @grasilQueue;

    # Initialize an array of inclinations.
    my $inclinations;

    # Loop over all nodes to process.
    foreach my $i ( $selection->list() ) {
	
	# Check for the existance of a Grasil SED for this galaxy.
	my $nodeIndex       = $dataSets->{'nodeIndex'      }->index($i);
	my $mergerTreeIndex = $dataSets->{'mergerTreeIndex'}->index($i);
	my $mergerTreeGroup;
	@groupsList = $outputGroup->groups();
	if ( defined(first { $_ eq "mergerTree".$mergerTreeIndex } @groupsList) ) {
	    $mergerTreeGroup = $outputGroup->group("mergerTree".$mergerTreeIndex);
	} else {
	    $mergerTreeGroup = new PDL::IO::HDF5::Group( name => "mergerTree".$mergerTreeIndex, parent => $outputGroup,
							 fileObj => $dataSet->{'hdf5File'} );
	}
	@groupsList = $outputGroup->groups();
	unless ( defined(first { $_ eq "mergerTree".$mergerTreeIndex } @groupsList) ) {
	    print "Failed to open grasilSEDs/Output".$outputIndex."/mergerTree".$mergerTreeIndex." in ".$dataSet->{'file'}."\n";
	    return;
	}
	my $nodeGroup;
	@groupsList = $mergerTreeGroup->groups();
	if ( defined(first { $_ eq "node".$nodeIndex } @groupsList) ) {
	    $nodeGroup = $mergerTreeGroup->group("node".$nodeIndex);
	} else {
	    $nodeGroup = new PDL::IO::HDF5::Group( name => "node".$nodeIndex, parent => $mergerTreeGroup,
						   fileObj => $dataSet->{'hdf5File'} );
	}
	@groupsList = $mergerTreeGroup->groups();
	unless ( defined(first { $_ eq "node".$nodeIndex } @groupsList) ) {
	    print "Failed to open grasilSEDs/Output".$outputIndex."/mergerTree".$mergerTreeIndex."/node".$nodeIndex." in ".$dataSet->{'file'}."\n";
	    return;
	}
	my $wavelengthDatasetName  = "grasilSEDs/Output".$outputIndex."/mergerTree".$mergerTreeIndex."/node".$nodeIndex."/wavelength";
	my $grasilDatasetName      = "grasilSEDs/Output".$outputIndex."/mergerTree".$mergerTreeIndex."/node".$nodeIndex."/SED";
	my $inclinationDatasetName = "grasilSEDs/Output".$outputIndex."/mergerTree".$mergerTreeIndex."/node".$nodeIndex."/inclination";
	my $wavelength;
	my $SED;
	@datasetsList = $nodeGroup->datasets();
	unless ( defined(first { $_ eq "SED" } @datasetsList) && $recomputeSEDs == 0 ) {
	    # Get a queue number for this galaxy.
	    my $queueNumber = scalar(@grasilQueue);

	    # Extract the star formation history for this galaxy.
	    system("mkdir -p ".$dataSet->{'file'}.".grasilTmp".$$.".".$queueNumber);
	    &Galacticus::Grasil::Extract_Star_Formation_History($dataSet,$outputIndex,$mergerTreeIndex,$nodeIndex,$dataSet->{'file'}.".grasilTmp".$$.".".$queueNumber."/grasil".$queueNumber.".dat");

	    # Generate a parameter file for Grasil.
	    open(iHndl,&galacticusPath()."/data/grasil/grasilBaseParameters.txt");
	    open(oHndl,">".$dataSet->{'file'}.".grasilTmp".$$.".".$queueNumber."/grasil".$queueNumber.".par");
	    my $inInclinations = 0;
	    $inclinations = pdl [];
	    while ( my $line = <iHndl> ) {
		$line =~ s/^(\s*pahflag\s+)[\d\.]+/$1$includePAHs/;
		$line =~ s/^(\s*flutflag\s+)[\d\.]+/$1$fluctuatingTemperatures/;
		$line =~ s/^(\s*nlf\s+)\d+/$1$wavelengthCount/;
		$line =~ s/^(\s*ndr\s+)\d+/$1$radialGridCount/;
		$line =~ s/^(\s*dsug\s+)[\d\.]+/$1$dustToMetalsRatio/;
		if ( $inInclinations == 1 ) {
		    if ( $line =~ m/^([\d\.]+)/ ) { 
			$inclinations = $inclinations->append($1);
		    } else {
			$inInclinations = 0;
		    }
		}
		$inInclinations = 1 if ( $line =~ m/^\d+\s+number of directions/ );
		print oHndl $line;
	    }
	    close(iHndl);
	    close(oHndl);
	    
	    # Add this galaxy to the queue for processing by Grasil.
	    $grasilQueue[$queueNumber] = {
		grasilFilesRoot => $dataSet->{'file'}.".grasilTmp".$$.".".$queueNumber."/grasil".$queueNumber,
		galaxyIndex     => $i,
		nodeGroup       => $nodeGroup,
		cpuLimit        => $cpuLimit
	    };

	    # Process through Grasil if the queue is full.
	    if ( scalar(@grasilQueue) >= $grasilJobsMaximum ) {
		&Process_Through_Grasil(\@grasilQueue,$inclinations,$dataSet,$dataSets,$dataSetName,$cpuLimit,$launchMethod);
		undef(@grasilQueue);
	    }

	} else {
	    # Grasil SED already exists - read it.
	    $wavelength   = $dataSet->{'hdf5File'}->dataset($wavelengthDatasetName )->get();
	    $SED          = $dataSet->{'hdf5File'}->dataset($grasilDatasetName     )->get();
	    $inclinations = $dataSet->{'hdf5File'}->dataset($inclinationDatasetName)->get();

	    # Compute the required Grasil property.
	    &Compute_Grasil_Property($dataSet,$wavelength,$SED,$inclinations,$dataSets,$dataSetName,$i);

	}
	
    }
    
    # Process through Grasil if any galaxies remain in the queue.
    if ( scalar(@grasilQueue) > 0 ) {
	&Process_Through_Grasil(\@grasilQueue,$inclinations,$dataSet,$dataSets,$dataSetName,$cpuLimit,$launchMethod);
	undef(@grasilQueue);
    }

}

sub Compute_Grasil_Property {
    # Compute the requested Grasil property from the SED.
    my $dataSet            = shift;
    my $wavelength         = shift;
    my $SED                = shift;
    my $inclinations       = shift;
    my $dataSets           = shift;
    my $dataSetName        = shift;
    my $i                  = shift;    
    # Using the SED, compute the flux and store it.
    if ( $dataSetName eq "grasilInfraredLuminosity" ) {
	# The total infrared luminosity was requested.
	&Compute_Infrared_Luminosity($wavelength,$SED,$inclinations,$dataSets,$dataSetName,$i);
    } elsif ( $dataSetName =~ m/^grasilFlux[\d\.]+microns$/ ) {
	# An individual flux was requested.
	my $observedWavelength;
	if ( $dataSetName =~ m/^grasilFlux([\d\.]+)microns/ ) {
	    $observedWavelength = pdl $1*1.0e4;
	}
	&Compute_Flux($wavelength,$SED,$inclinations,$observedWavelength,$dataSets,$dataSetName,$i);
    } elsif ( $dataSetName =~ m/^luminosity:grasil:([^:]+):([^:]+)$/ ) {
	&Compute_Luminosity_Under_Filter($wavelength,$SED,$inclinations,$dataSets,$dataSetName,$i);
    } else {
	# Unrecognized property name.
	die("Galacticus::Grasil::Compute_Grasil_Property: unknown property name");
    }
}

sub Extract_Star_Formation_History {
    my $dataSet     = shift;
    my $outputIndex = shift;
    my $treeIndex   = shift;
    my $nodeIndex   = shift;
    my $grasilFile  = shift;
    my $plotFile    = shift;
    
    # Define prefixes.
    my $giga        = 1.0e9;

    # Define Solar metallicity.
    my $metallicitySolar = 0.0188;

    # Open groups.
    my $nodeDataGroup            = $dataSet->{'hdf5File'}->group("Outputs")->group("Output".$outputIndex)->group("nodeData");
    my $starFormationGroup       = $dataSet->{'hdf5File'}->group("starFormationHistories");
    my $starFormationOutputGroup = $starFormationGroup->group("Output".$outputIndex); 
    my $starFormationTreeGroup   = $starFormationOutputGroup->group("mergerTree".$treeIndex);

    # Get the output time.
    my @outputTime          = $dataSet->{'hdf5File'}->group("Outputs")->group("Output".$outputIndex)->attrGet("outputTime");
    
    # Get the metallicities at which star formation rates are tabulated.
    my $metallicities       = $starFormationGroup->dataset("metallicities")->get();

    # Read the tree indices, offsets and node counts.
    my $treeIndices         = $dataSet->{'hdf5File'}->group("Outputs")->group("Output".$outputIndex)->dataset("mergerTreeIndex"     )->get();
    my $treeNodeStart       = $dataSet->{'hdf5File'}->group("Outputs")->group("Output".$outputIndex)->dataset("mergerTreeStartIndex")->get();
    my $treeNodeCount       = $dataSet->{'hdf5File'}->group("Outputs")->group("Output".$outputIndex)->dataset("mergerTreeCount"     )->get();

    # Get the offset positions for the tree.
    my $selection = which($treeIndices == $treeIndex);
    my $start = $treeNodeStart->index($selection);
    my $count = $treeNodeCount->index($selection);
    my $end   = $start+$count-1;
  
    # Read in galaxy data.
    my $nodeIndices         = $nodeDataGroup->dataset("nodeIndex"                  )->get($start,$end);
    my $diskScaleLength     = $nodeDataGroup->dataset("diskRadius"                 )->get($start,$end);
    my $spheroidScaleLength = $nodeDataGroup->dataset("spheroidRadius"             )->get($start,$end);
    my $diskStellarMass     = $nodeDataGroup->dataset("diskMassStellar"            )->get($start,$end);
    my $spheroidStellarMass = $nodeDataGroup->dataset("spheroidMassStellar"        )->get($start,$end);
    my $diskGasMass         = $nodeDataGroup->dataset("diskMassGas"                )->get($start,$end);
    my $spheroidGasMass     = $nodeDataGroup->dataset("spheroidMassGas"            )->get($start,$end);
    my $diskGasMetals       = $nodeDataGroup->dataset("diskAbundancesGasMetals"    )->get($start,$end);
    my $spheroidGasMetals   = $nodeDataGroup->dataset("spheroidAbundancesGasMetals")->get($start,$end);

    # Find the node in question.
    my $selected = which($nodeIndices == $nodeIndex);

    # Convert from total metals to metallicity.
    my $diskGasMetallicity;
    if ( $diskGasMass    (($selected)) > 0.0 ) {
 	$diskGasMetallicity      = $diskGasMetals    (($selected))/$diskGasMass    (($selected));
 	$diskGasMetallicity     .= 1.0 if ( $diskGasMetallicity > 1.0 );
 	$diskGasMetallicity     /= $metallicitySolar;
    } else {
 	$diskGasMetallicity      = 0.0;
    }
    my $spheroidGasMetallicity;
    if ( $spheroidGasMass(($selected)) > 0.0 ) {
     	$spheroidGasMetallicity  = $spheroidGasMetals(($selected))/$spheroidGasMass(($selected));
     	$spheroidGasMetallicity .= 1.0 if ( $spheroidGasMetallicity > 1.0 );
     	$spheroidGasMetallicity /= $metallicitySolar;
    } else {
     	$spheroidGasMetallicity  = 0.0;
    }
    
    # Get a list of available datasets and convert to a hash.
    my @dataSets = $starFormationTreeGroup->datasets();
    my %availableDatasets;
    foreach my $dataSet ( @dataSets ) {
 	$availableDatasets{$dataSet} = 1;
    }

    # Read in the star formation data.
    my ($diskTime, $diskSFH);
    if ( exists($availableDatasets{"diskTime".$nodeIndex}) ) {
 	$diskTime            = $dataSet->{'hdf5File'}->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->dataset("diskTime"    .$nodeIndex)->get;
 	$diskSFH             = $dataSet->{'hdf5File'}->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->dataset("diskSFH"     .$nodeIndex)->get;
 	$diskSFH->where($diskSFH < 0.0) .= 0.0;
    } else {
 	$diskTime            = ones  (1);
 	$diskSFH             = zeroes(1,nelem($metallicities));
    }
    my ($spheroidTime, $spheroidSFH);
    if ( exists($availableDatasets{"spheroidTime".$nodeIndex}) ) {
 	$spheroidTime        = $dataSet->{'hdf5File'}->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->dataset("spheroidTime".$nodeIndex)->get;
 	$spheroidSFH         = $dataSet->{'hdf5File'}->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->dataset("spheroidSFH" .$nodeIndex)->get;
 	$spheroidSFH->where($spheroidSFH < 0.0) .= 0.0;
    } else {
 	$spheroidTime        = ones  (1);
 	$spheroidSFH         = zeroes(1,nelem($metallicities));
    }

    # Compute time steps.
    my $diskTimeBegin       = pdl [0.0];
    if ( nelem($diskTime) > 1 ) {
 	$diskTimeBegin       = $diskTimeBegin->append($diskTime->index(sequence(nelem($diskTime)-1)));
    } else {
 	$diskTimeBegin       = $diskTimeBegin;
    }
    my $diskTimeStep        = $diskTime-$diskTimeBegin;
    my $diskTimeCentral     = ($diskTime+$diskTimeBegin)/2.0;
    my $spheroidTimeBegin   = pdl [0.0];
    if ( nelem($spheroidTime) > 1 ) {
 	$spheroidTimeBegin   = $spheroidTimeBegin->append($spheroidTime->index(sequence(nelem($spheroidTime)-1)));
    } else {
 	$spheroidTimeBegin   = $spheroidTimeBegin;
    }
    my $spheroidTimeStep    = $spheroidTime-$spheroidTimeBegin;
    my $spheroidTimeCentral = ($spheroidTime+$spheroidTimeBegin)/2.0;

    # Open the Grasil output file.
    open(gHndl,">".$grasilFile);

    # Output file header.
    print gHndl "# Output index     :\t".$outputIndex."\n";
    print gHndl "# Tree   index     :\t".$treeIndex."\n";
    print gHndl "# Node   index     :\t".$nodeIndex."\n";
    print gHndl "# Output time [Gyr]:\t".$outputTime[0]."\n";
    print gHndl "#\n";
    print gHndl "# Galaxy properties:\n";
    print gHndl "#  Stellar mass    (disk, spheroid) [M_Solar]:\t".$diskStellarMass   (($selected))."\t".$spheroidStellarMass   (($selected))."\n";
    print gHndl "#  Gas     mass    (disk, spheroid) [M_Solar]:\t".$diskGasMass       (($selected))."\t".$spheroidGasMass       (($selected))."\n";
    print gHndl "#  Scale length    (disk, spheroid) [Mpc    ]:\t".$diskScaleLength   (($selected))."\t".$spheroidScaleLength   (($selected))."\n";
    print gHndl "#  Gas metallicity (disk, spheroid) [       ]:\t".$diskGasMetallicity             ."\t".$spheroidGasMetallicity             ."\n";

    # See if we can check the integration of star formation histories.
    print gHndl "#\n";
    my @starFormationParameters = $dataSet->{'hdf5File'}->group("Parameters")->attrGet("imfSelectionMethod","stellarPopulationPropertiesMethod");
    if ( $starFormationParameters[0] eq "fixed" && $starFormationParameters[1] eq "instantaneous" ) {
 	my @imfSelectionFixed = $dataSet->{'hdf5File'}->group("Parameters")->attrGet("imfSelectionFixed");
 	my $imfRecycledAttributeName = "imf".$imfSelectionFixed[0]."RecycledInstantaneous";
 	my @recycledFraction = $dataSet->{'hdf5File'}->group("Parameters")->attrGet($imfRecycledAttributeName);
 	my $diskSFHIntegrated     = (1.0-$recycledFraction[0])*$diskSFH    ->sum;
 	my $spheroidSFHIntegrated = (1.0-$recycledFraction[0])*$spheroidSFH->sum;
 	print gHndl "# Fractional error in SFH integration:\n";
 	if ( $diskStellarMass    (($selected)) > 0.0 ) {
 	    print gHndl "#  Disk    :\t".abs($diskSFHIntegrated    -$diskStellarMass    (($selected)))/$diskStellarMass    (($selected))."\n";
 	} else {
 	    print gHndl "#  Disk    :\tN/A\n";
 	}
 	if ( $spheroidStellarMass(($selected)) > 0.0 ) {
 	    print gHndl "#  Spheroid:\t".abs($spheroidSFHIntegrated-$spheroidStellarMass(($selected)))/$spheroidStellarMass(($selected))."\n";
 	} else {
 	    print gHndl "#  Spheroid:\tN/A\n";
 	}
    } else {
 	print gHndl "# Checks of star formation history integration are disabled for this model.\n";
    }

    # Output the metallicities.
    print gHndl "#\n";
    print gHndl "# Metallicities: ".join("\t",$metallicities->list)."\n";
    
    # Compute mean star formation rates.
    my $diskSFR             = $diskSFH    /$diskTimeStep    /$giga;
    my $spheroidSFR         = $spheroidSFH/$spheroidTimeStep/$giga;

    # Make a plot.
    my $sfrMinimum = 1.0e-2;
    if ( defined($plotFile) && ( any($diskSFR > $sfrMinimum) || any($spheroidSFR > $sfrMinimum) ) ) {
 	# Declare variables for GnuPlot;
 	my ($gnuPlot, $outputFile, $outputFileEPS, $plot);

 	$outputFile = $plotFile;
 	($outputFileEPS = $outputFile) =~ s/\.pdf$/.eps/;
 	open($gnuPlot,"|gnuplot");
 	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
 	print $gnuPlot "set output '".$outputFileEPS."'\n";
 	print $gnuPlot "set lmargin screen 0.15\n";
 	print $gnuPlot "set rmargin screen 0.95\n";
 	print $gnuPlot "set bmargin screen 0.15\n";
 	print $gnuPlot "set tmargin screen 0.95\n";
 	print $gnuPlot "set key spacing 1.2\n";
 	print $gnuPlot "set key at screen 0.45,0.5\n";
 	print $gnuPlot "set key left\n";
 	print $gnuPlot "set key bottom\n";
 	print $gnuPlot "set logscale y\n";
 	print $gnuPlot "set mytics 10\n";
 	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
 	print $gnuPlot "set xlabel '\$t\$ [Gyr]'\n";
 	print $gnuPlot "set ylabel '\$\\dot{M}_\\star\$ \$[M_\\odot \\hbox{yr}^{-1}]\$'\n";
 	print $gnuPlot "set xrange [0.0:15.0]\n";
 	print $gnuPlot "set yrange [1.0e-2:1.0e2]\n";
 	my $iColor = -1;
 	for(my $i=0;$i<$spheroidSFR->dim(1);++$i) {
 	    ++$iColor;
 	    if (any($spheroidSFR(:,($i)) > $sfrMinimum)) {
 		my ($metallicityLow, $metallicityHigh);
 		if ( $i == 0 ) {
 		    $metallicityLow = 0.0;
 		} else {
 		    $metallicityLow = FormatSigFigs($metallicities->index($i-1),2);
 		}
 		if ( $i == $spheroidSFR->dim(1)-1 ) {
 		    $metallicityHigh = "\\\\infty";
 		} else {
 		    $metallicityHigh = FormatSigFigs($metallicities->index($i),2);
 		}
 		my $label = "\\\\small Spheroid: \$".$metallicityLow."<Z<".$metallicityHigh."\$";
 		&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
 					      $spheroidTimeCentral, $spheroidSFR(:,($i)),
 					      style => "point", symbol => [4,5], weight => [5,3],
 					      color => $GnuPlot::PrettyPlots::colorPairs{${$GnuPlot::PrettyPlots::colorPairSequences{'sequence1'}}[$iColor]},
 					      title => $label);
 	    }
 	}
 	$iColor = -1;
 	for(my $i=0;$i<$diskSFR->dim(1);++$i) {
 	    ++$iColor;
 	    if (any($diskSFR(:,($i)) > $sfrMinimum)) {
 		my ($metallicityLow, $metallicityHigh);
 		if ( $i == 0 ) {
 		    $metallicityLow = 0.0;
 		} else {
 		    $metallicityLow = FormatSigFigs($metallicities->index($i-1),2);
 		}
 		if ( $i == $diskSFR->dim(1)-1 ) {
 		    $metallicityHigh = "\\\\infty";
 		} else {
 		    $metallicityHigh = FormatSigFigs($metallicities->index($i),2);
 		}
 		my $label = "\\\\small Disk: \$".$metallicityLow."<Z<".$metallicityHigh."\$";
 		&GnuPlot::PrettyPlots::Prepare_Dataset(\$plot,
 					      $diskTimeCentral, $diskSFR(:,($i)),
 					      style => "point", symbol => [6,7], weight => [5,3],
 					      color => $GnuPlot::PrettyPlots::colorPairs{${$GnuPlot::PrettyPlots::colorPairSequences{'sequence1'}}[$iColor]},
 					      title => $label);
 	    }
 	}
 	&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
 	close($gnuPlot);
 	&GnuPlot::LaTeX::GnuPlot2PDF($outputFileEPS);
	
    }

    # Write disk SFR data.
    print gHndl "#\n";
    print gHndl "# Disk star formation rate:\n";
    print gHndl "# Time [Gyr] | SFR [M_Solar/yr]\n";
    my $table = Text::Table->new();
    for(my $j=0;$j<$diskTimeCentral->nelem;++$j) {
 	my @rowData;
 	$rowData[0] = $diskTimeCentral->index($j);
 	push(@rowData,$diskSFR(($j),:)->list);
 	$table->add(@rowData);
    }
    print gHndl $table;
    
    # Write spheroid SFR data.
    print gHndl "#\n";
    print gHndl "# Spheroid star formation rate:\n";
    print gHndl "# Time [Gyr] | SFR [M_Solar/yr]\n";
    $table = Text::Table->new();
    for(my $j=0;$j<$spheroidTimeCentral->nelem;++$j) {
	my @rowData;
 	$rowData[0] = $spheroidTimeCentral->index($j);
 	push(@rowData,$spheroidSFR(($j),:)->list);
 	$table->add(@rowData);
    }
    print gHndl $table;
    close(gHndl);
}

sub Compute_Flux {
    # Compute the flux at a given wavelength and inclination given the SED of a galaxy.
    my $wavelength         = shift;
    my $SED                = shift;
    my $inclinations       = shift;
    my $observedWavelength = shift;
    my $dataSets           = shift;
    my $dataSetName        = shift;
    my $i                  = shift;

    # Interpolate in wavelength and inclination to get luminosity.
    my $inclination        = $dataSets->{"inclination"       }->index($i);
    my $redshift           = $dataSets->{"redshift"          }->index($i);
    my $luminosityDistance = $dataSets->{"luminosityDistance"}->index($i);
    my $restWavelength     = pdl $observedWavelength/(1.0+$redshift);
    my $inclinationIndices = sequence(nelem($inclinations));
    (my $inclinationIndex, my $inclinationError) = interpolate($inclination,$inclinations,$inclinationIndices);
    my $wavelengthIndices = sequence(nelem($wavelength));
    (my $wavelengthIndex, my $wavelengthError) = interpolate($restWavelength,$wavelength,$wavelengthIndices);
    my $index = pdl([$inclinationIndex,$wavelengthIndex]);
    my $luminosity = $SED->interpND($index);
    
    # Convert luminosity to flux.
    my $flux = $conversionFactor*$luminosity*(1.0+$redshift)*($restWavelength**2)/4.0/$Pi/$luminosityDistance**2;
    
    # Store in the dataset.
    $dataSets->{$dataSetName}->index($i) .= $flux;
}

sub Compute_Infrared_Luminosity {
    # Compute the total infrared (8-1000 microns) luminosity given the SED of a galaxy.
    my $wavelength         = shift;
    my $SED                = shift;
    my $inclinations       = shift;
    my $dataSets           = shift;
    my $dataSetName        = shift;
    my $i                  = shift;

    # Find unique wavelengths.
    my $uniqueWavelengths = $wavelength->uniqind();

    # Set up an interpolator and integrate over the luminosity range.
    my $interpolator = PDL::GSL::INTERP->init('linear',$wavelength->index($uniqueWavelengths),$SED((0),:)->index($uniqueWavelengths),{Sort => 0}); 
    my $luminosity = $interpolator->integ($irWavelengthMinimum,$irWavelengthMaximum,{Extrapolate => 0});

    # Convert to Solar luminosities.
    $luminosity /= $solarLuminosity;

    # Store in the dataset.
    $dataSets->{$dataSetName}->index($i) .= $luminosity;

}

sub Compute_Flux_Under_Filter {
    # Compute the flux in a filter (in Janskys).
    my $dataSet     = shift;
    my $dataSetName = $_[0];
    # Construct the corresponding luminosity name.
    (my $luminosityDataSet = $dataSetName) =~ s/^grasilFlux:([^:]+)$/luminosity:grasil:$1:observed/;
    # Get the luminosity dataset.
    &Galacticus::HDF5::Get_Dataset($dataSet,[$luminosityDataSet,"redshift","luminosityDistance"]);
    # Compute the flux in this filter.
    $dataSet->{'dataSets'}->{$dataSetName} = $luminosityAB*$dataSet->{'dataSets'}->{$luminosityDataSet}*(1.0+$dataSet->{'dataSets'}->{"redshift"})/4.0/$Pi/$dataSet->{'dataSets'}->{"luminosityDistance"}**2/$megaParsec**2/$Jansky;
}

sub Compute_Luminosity_Under_Filter {
    # Compute the luminosity in a filter (in units of "maggies", i.e. the zero-point of the AB magnitude system).
    my $wavelength   = shift;
    my $SED          = shift;
    my $inclinations = shift;
    my $dataSets     = shift;
    my $dataSetName  = shift;
    my $i            = shift;

    # Match the dataset name.
    if ( $dataSetName =~ m/^luminosity:grasil:([^:]+):([^:]+)$/ ) {
	# Extract the name of the filter and the required frame.
	my $filterName   = $1;
	my $frame        = $2;
	# Get the redshift for this galaxy.
	my $redshift     = $dataSets->{'redshift'   }->(($i));
	# Get the inclination for this galaxy.
	my $inclination  = $dataSets->{'inclination'}->(($i));
	# Interpolate the SED in inclination.
	my $inclinationIndices = sequence(nelem($inclinations));
	(my $inclinationIndex, my $inclinationError) = interpolate($inclination,$inclinations,$inclinationIndices);
	my $indices = zeroes(2,nelem($wavelength));
	$indices->((0),:) .= $inclinationIndex;
	$indices->((1),:) .= sequence(nelem($wavelength));
	my $sedInclined = $SED->interpND($indices);
	# Load the filter.
	(my $filterWavelengths, my $filterResponse) = &Galacticus::Filters::Load($filterName);	
	# Make a joint set of filter and SED wavelengths.
	my $filterWavelengthRestFrame = $filterWavelengths->copy();
	$filterWavelengthRestFrame /= (1.0+$redshift)
	    if ( $frame eq "observed" );
	my $jointWavelengths = $filterWavelengthRestFrame->copy();
	$jointWavelengths  = $jointWavelengths->append($wavelength);
	my $nonZero        = which(
	    ($jointWavelengths >= $filterWavelengthRestFrame((0)))
	    &
	    ($jointWavelengths <= $filterWavelengthRestFrame((-1)))
	    );
	$jointWavelengths  = $jointWavelengths->($nonZero)->qsort();
	# Interpolate the filter response onto the joint wavelengths.
	(my $jointResponse, my $responseError) = interpolate($jointWavelengths,$filterWavelengthRestFrame,$filterResponse);
	# Interpolate the SED onto the joint wavelengths.
	(my $jointSED     , my $sedError     ) = interpolate($jointWavelengths,$wavelength               ,$sedInclined   );
	# Generate a set of delta wavelengths for use in integrations.
	my $deltaWavelengths = pdl [];
	for(my $i=0;$i<nelem($jointWavelengths);++$i) {
	    my $deltaWavelength;
	    if ( $i == 0 ) {
		$deltaWavelength = $jointWavelengths((1))-$jointWavelengths((0));
	    } elsif ( $i == nelem($jointWavelengths)-1 ) {
		$deltaWavelength = $jointWavelengths((-1))-$jointWavelengths((-2));
	    } else {
		$deltaWavelength = $jointWavelengths(($i+1))-$jointWavelengths(($i-1));
	    }
	    $deltaWavelengths = $deltaWavelengths->append($deltaWavelength/2.0);
	}	
	# Integrate the SED to get the luminosity.
	my $luminosity = 
	    sum($grasilNormalization*$jointSED    *$jointResponse*$deltaWavelengths*$jointWavelengths*$Angstrom/$speedOfLight)
	    /
	    sum(                     $luminosityAB*$jointResponse*$deltaWavelengths/$jointWavelengths                        );
	# Store in the dataset.
	$dataSets->{$dataSetName}->index($i) .= $luminosity;
    }
}

sub Process_Through_Grasil {
    # Process a list of galaxies through Grasil.
    my @grasilQueue        = @{shift()};
    my $inclinations       =   shift() ;
    my $dataSet            =   shift() ;
    my $dataSets           =   shift() ;
    my $dataSetName        =   shift() ;
    my $cpuLimit           =   shift() ;
    my $launchMethod       =   shift() ;
    # Ensure that we have Grasil.
    &grasilGet();
    # Iterate over galaxies, generating jobs.
    my @submitQueue;
    foreach my $grasilJob ( @grasilQueue ) {
	(my $grasilDirectoryName = $grasilJob->{'grasilFilesRoot'}) =~ s/\/([^\/]+)$//;
	my $grasilRun = "cd `dirname ".$grasilJob->{'grasilFilesRoot'}."`; ".&galacticusPath()."/aux/Grasil/grasil `basename ".$grasilJob->{'grasilFilesRoot'}."`";
	my %job =
	    (
	     launchFile   => $grasilDirectoryName."/grasil.script",
	     logFile      => $grasilDirectoryName."/grasil.log"   ,
	     label        => "grasil",
	     command      => $grasilRun,
	     onCompletion =>
	     {
		 function  => \&grasilPostProcess,
		 arguments => [ $grasilJob, $dataSet, $dataSets, $dataSetName, $inclinations ]
	     }
	    );
	foreach ( 'ppn', 'walltime', 'memory' ) {
	    $job{$_} = $dataSet->{'grasilOptions'}->{$_}
	       if ( exists($dataSet->{'grasilOptions'}->{$_}) );
	}
	push(@submitQueue,\%job);
    }
    # Launch the jobs.
    $dataSet->{'grasilOptions'}->{'allowThreads'} = 0;
    &{$Hooks::moduleHooks{$launchMethod}->{'jobArrayLaunch'}}($dataSet->{'grasilOptions'},@submitQueue);
}

sub grasilPostProcess {
    # Postprocess output of a Grasil job.
    my $grasilJob    = shift();
    my $dataSet      = shift();
    my $dataSets     = shift();
    my $dataSetName  = shift();
    my $inclinations = shift();
    my $jobID        = shift();
    my $exitStatus   = shift();
    # Check that Grasil completed.	
    if ( $exitStatus ==0 && -e $grasilJob->{'grasilFilesRoot'}.".spe" ) {
	# Count number of wavelengths in output file.
	my $wavelengthCount = 0;
	open(my $grasilSED,$grasilJob->{'grasilFilesRoot'}.".spe");
	while ( my $line = <$grasilSED> ) {
	    ++$wavelengthCount unless ( $line =~ m/^\#/ );
	}
	close($grasilSED);
	# Parse and store the output SEDs.
	my $wavelength   = pdl zeroes(                     $wavelengthCount);
	my $SED          = pdl zeroes(nelem($inclinations),$wavelengthCount);
	my $iInclination = -1;
	foreach my $inclination ( $inclinations->list() ) {
	    ++$iInclination;
	    my $iWavelength = -1;
	    open($grasilSED,$grasilJob->{'grasilFilesRoot'}.".".$inclination);
	    while ( my $line = <$grasilSED> ) {
		unless ( $line =~ m/^\#/ ) {
		    ++$iWavelength;
		    $line =~ s/^\s*//;
		    $line =~ s/\s*$//;
		    my @columns = split(/\s+/,$line);
		    $wavelength(                ($iWavelength)) .= $columns[0];
		    $SED       (($iInclination),($iWavelength)) .= $columns[5];
		}
	    }
	    close($grasilSED);
	}
	# Check if the datasets already exist.
	my @existingDatasets = $grasilJob->{'nodeGroup'}->datasets();
	if ( defined(first { $_ eq "wavelength" } @existingDatasets) ) {
	    # Determine size of stored wavelength dataset.
	    my $storedWavelength = $grasilJob->{'nodeGroup'}->dataset('wavelength')->get();
	    if ( nelem($storedWavelength) != nelem($wavelength) ) {
		# New dataset differs in size from the previously stored on. We therefore need to unlink the stored datasets
		# so that we can make new ones.
		$grasilJob->{'nodeGroup'}->unlink($_)
		    foreach ( 'wavelength', 'inclination', 'SED' );
	    }
	}
	# Write datasets.
	my $wavelengthDataset  = $grasilJob->{'nodeGroup'}->dataset('wavelength' );
	my $sedDataset         = $grasilJob->{'nodeGroup'}->dataset('SED'        );
	my $inclinationDataset = $grasilJob->{'nodeGroup'}->dataset('inclination');
	$wavelengthDataset ->set($wavelength  );
	$sedDataset        ->set($SED         );
	$inclinationDataset->set($inclinations);	    
	# Compute the required Grasil property.
	&Compute_Grasil_Property($dataSet,$wavelength,$SED,$inclinations,$dataSets,$dataSetName,$grasilJob->{'galaxyIndex'});
	# Remove the folder.
	system("rm -rf `dirname ".$grasilJob->{'grasilFilesRoot'}."`");
    } else {
	# Grasil did not complete.
	my $workDirectory = $grasilJob->{'grasilFilesRoot'};
	$workDirectory =~ s/\/[^\/]+$//;
	my $baseDirectory = $workDirectory;
	$baseDirectory =~ s/\/[^\/]+$//;
	system("mv ".$workDirectory." ".$baseDirectory."/".$grasilJob->{'galaxyIndex'});
	print "Galacticus::Grasil::grasilPostprocess(): Grasil appears to have FAILED for ".$baseDirectory."/".$grasilJob->{'galaxyIndex'}."\n";
	exit;
    }
}

sub grasilGet {
    # Ensure that we have the Grasil executable and data files.
    my $grasilPath = &galacticusPath()."/aux/Grasil/";
    system("mkdir -p ".$grasilPath);
    system("wget http://adlibitum.oats.inaf.it/silva/grasil/download/gslib.tar.gz -O ".$grasilPath."/gslib.tar.gz")
	unless ( -e $grasilPath."/gslib.tar.gz" );
    die("Galacticus::Grasil::grasilGet(): : unable to download gslib.tar.gz")
	unless ( $? == 0 );
    unless ( -e $grasilPath."/TAU96.OF" ) {
	system("cd ".$grasilPath."; tar xvfz gslib.tar.gz");
	die("Galacticus::Grasil::grasilGet(): : unable to extract gslib.tar.gz")
	    unless ( $? == 0 );
	unlink($grasilPath."/grasil")
	    if ( -e $grasilPath."/grasil" );
    }
    system("wget http://adlibitum.oats.inaf.it/silva/grasil/SSP_zip/grasil.tar.gz -O ".$grasilPath."/grasil.tar.gz")
	unless ( -e $grasilPath."/grasil.tar.gz" );
    die("Galacticus::Grasil::grasilGet(): : unable to download grasil.tar.gz")
	unless ( $? == 0 );
    unless ( -e $grasilPath."/template_grasil.sf" ) {
	system("cd ".$grasilPath."; tar xvfz grasil.tar.gz");
	die("Galacticus::Grasil::grasilGet(): : unable to extract gslib.tar.gz")
	    unless ( $? == 0 );
	unlink($grasilPath."/grasil")
	    if ( -e $grasilPath."/grasil" );
    }
    system("wget http://users.obs.carnegiescience.edu/abenson/galacticus/tools/grasil -O ".$grasilPath."/grasil; chmod u=wrx ".$grasilPath."/grasil")
	unless ( -e $grasilPath."/grasil" );
    die("Galacticus::Grasil::grasilGet(): : unable to download grasil")
	unless ( $? == 0 );
}

1;

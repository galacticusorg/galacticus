#!/usr/bin/env perl
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use PDL::NiceSlice;
use Text::Table;
use Math::SigFigs;

# Plot the star formation history of a galaxy split by metallicity and output the data in a form suitable for input to Grasil.
# Andrew Benson (06-September-2010)

# Read arguments.
die ("Usage: Extract_Star_Formation_History_for_Grasil.pl <inputFile> <outputIndex> <treeIndex> <nodeIndex> <grasilFile> [<plotFile>]") unless ( $#ARGV == 4 || $#ARGV == 5 );
$inputFile   = $ARGV[0];
$outputIndex = $ARGV[1];
$treeIndex   = $ARGV[2];
$nodeIndex   = $ARGV[3];
$grasilFile  = $ARGV[4];
$plotFile    = $ARGV[5] if ( $#ARGV == 5 );

# Define prefixes.
$giga        = 1.0e9;

# Define Solar metallicity.
$metallicitySolar = 0.0188;

# Open the Galacticus file.
$HDFfile = new PDL::IO::HDF5($inputFile);

# Get the output time.
@outputTime          = $HDFfile->group("Outputs/Output".$outputIndex)->attrGet("outputTime");

# Get the metallicities at which star formation rates are tabulated.
$metallicities       = $HDFfile->group("starFormationHistories")->dataset("metallicities")->get;

# Read the tree indices, offsets and node counts.
$treeIndices         = $HDFfile->group("Outputs/Output".$outputIndex)->dataset("mergerTreeIndex"     )->get;
$treeNodeStart       = $HDFfile->group("Outputs/Output".$outputIndex)->dataset("mergerTreeStartIndex")->get;
$treeNodeCount       = $HDFfile->group("Outputs/Output".$outputIndex)->dataset("mergerTreeCount"     )->get;

# Get the offset positions for the tree.
$selection = which($treeIndices == $treeIndex);
$start = $treeNodeStart->index($selection);
$count = $treeNodeCount->index($selection);
$end   = $start+$count-1;

# Read in galaxy data.
$nodeIndices         = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("nodeIndex"          )->get($start,$end);
$diskScaleLength     = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("diskScaleLength"    )->get($start,$end);
$spheroidScaleLength = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("spheroidScaleLength")->get($start,$end);
$diskStellarMass     = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("diskStellarMass"    )->get($start,$end);
$spheroidStellarMass = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("spheroidStellarMass")->get($start,$end);
$diskGasMass         = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("diskGasMass"        )->get($start,$end);
$spheroidGasMass     = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("spheroidGasMass"    )->get($start,$end);
$diskGasMetals       = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("diskGasMetals"      )->get($start,$end);
$spheroidGasMetals   = $HDFfile->group("Outputs/Output".$outputIndex."/nodeData")->dataset("spheroidGasMetals"  )->get($start,$end);

# Find the node in question.
$selected = which($nodeIndices == $nodeIndex);

# Convert from total metals to metallicity.
if ( $diskGasMass    (($selected)) > 0.0 ) {
    $diskGasMetallicity      = $diskGasMetals    (($selected))/$diskGasMass    (($selected));
    $diskGasMetallicity     .= 1.0 if ( $diskGasMetallicity > 1.0 );
    $diskGasMetallicity     /= $metallicitySolar;
} else {
    $diskGasMetallicity      = 0.0;
}
if ( $spheroidGasMass(($selected)) > 0.0 ) {
    $spheroidGasMetallicity  = $spheroidGasMetals(($selected))/$spheroidGasMass(($selected));
    $spheroidGasMetallicity .= 1.0 if ( $spheroidGasMetallicity > 1.0 );
    $spheroidGasMetallicity /= $metallicitySolar;
} else {
    $spheroidGasMetallicity  = 0.0;
}

# Get a list of available datasets and convert to a hash.
@dataSets = $HDFfile->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->datasets;
foreach $dataSet ( @dataSets ) {
    $availableDatasets{$dataSet} = 1;
}

# Read in the star formation data.
if ( exists($availableDatasets{"diskTime".$nodeIndex}) ) {
    $diskTime            = $HDFfile->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->dataset("diskTime"    .$nodeIndex)->get;
    $diskSFH             = $HDFfile->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->dataset("diskSFH"     .$nodeIndex)->get;
} else {
    $diskTime            = ones  (1);
    $diskSFH             = zeroes(1,nelem($metallicities));
}
if ( exists($availableDatasets{"spheroidTime".$nodeIndex}) ) {
    $spheroidTime        = $HDFfile->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->dataset("spheroidTime".$nodeIndex)->get;
    $spheroidSFH         = $HDFfile->group("starFormationHistories/Output".$outputIndex."/mergerTree".$treeIndex)->dataset("spheroidSFH" .$nodeIndex)->get;
} else {
    $spheroidTime        = ones  (1);
    $spheroidSFH         = zeroes(1,nelem($metallicities));
}

# Compute time steps.
$diskTimeBegin       = pdl [0.0];
if ( nelem($diskTime) > 1 ) {
    $diskTimeBegin       = $diskTimeBegin->append($diskTime->index(sequence(nelem($diskTime)-1)));
} else {
    $diskTimeBegin       = $diskTimeBegin;
}
$diskTimeStep        = $diskTime-$diskTimeBegin;
$diskTimeCentral     = ($diskTime+$diskTimeBegin)/2.0;
$spheroidTimeBegin   = pdl [0.0];
if ( nelem($spheroidTime) > 1 ) {
    $spheroidTimeBegin   = $spheroidTimeBegin->append($spheroidTime->index(sequence(nelem($spheroidTime)-1)));
} else {
    $spheroidTimeBegin   = $spheroidTimeBegin;
}
$spheroidTimeStep    = $spheroidTime-$spheroidTimeBegin;
$spheroidTimeCentral = ($spheroidTime+$spheroidTimeBegin)/2.0;

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
@starFormationParameters = $HDFfile->group("Parameters")->attrGet("imfSelectionMethod","stellarPopulationPropertiesMethod");
if ( $starFormationParameters[0] eq "fixed" && $starFormationParameters[1] eq "instantaneous" ) {
    @imfSelectionFixed = $HDFfile->group("Parameters")->attrGet("imfSelectionFixed");
    $imfRecycledAttributeName = "imf".$imfSelectionFixed[0]."RecycledInstantaneous";
    @recycledFraction = $HDFfile->group("Parameters")->attrGet($imfRecycledAttributeName);
    $diskSFHIntegrated     = (1.0-$recycledFraction[0])*$diskSFH    ->sum;
    $spheroidSFHIntegrated = (1.0-$recycledFraction[0])*$spheroidSFH->sum;
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
$diskSFR             = $diskSFH    /$diskTimeStep    /$giga;
$spheroidSFR         = $spheroidSFH/$spheroidTimeStep/$giga;

# Make a plot.
$sfrMinimum = 1.0e-2;
if ( defined($plotFile) && ( any($diskSFR > $sfrMinimum) || any($spheroidSFR > $sfrMinimum) ) ) {
    open(gnuPlot,"|gnuplot");
    print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
    print gnuPlot "set output \"tmp.ps\"\n";
    print gnuPlot "set xlabel \"t [Gyr]\"\n";
    print gnuPlot "set ylabel \"SFR [M yr^{-1}]\"\n";
    print gnuPlot "set title \"Star Formation Rate\"\n";
    print gnuPlot "set logscale y\n";
    print gnuPlot "set format y \ \"10^{\%L}\"\n";
    print gnuPlot "set key left\n";
    print gnuPlot "set pointsize 1.0\n";
    $plotCommand  = "plot";
    $separator = " ";
    for($i=0;$i<$diskSFR->dim(1);++$i) {
	if (any($diskSFR(:,($i)) > $sfrMinimum)) {
	    if ( $i == 0 ) {
		$metallicityLow = 0.0;
	    } else {
		$metallicityLow = FormatSigFigs($metallicities->index($i-1),3);
	    }
	    if ( $i == $diskSFR->dim(1)-1 ) {
		$metallicityHigh = "\infty";
	    } else {
		$metallicityHigh = FormatSigFigs($metallicities->index($i),3);
	    }
	    $label = "Disk: ".$metallicityLow."<Z<".$metallicityHigh;
	    $plotCommand .= $separator."'-' pt 6 title \"".$label."\"";
	    $separator = ", ";
	}
    }
    for($i=0;$i<$spheroidSFR->dim(1);++$i) {
	if (any($spheroidSFR(:,($i)) > $sfrMinimum)) {
	    if ( $i == 0 ) {
		$metallicityLow = 0.0;
	    } else {
		$metallicityLow = FormatSigFigs($metallicities->index($i-1),3);
	    }
	    if ( $i == $diskSFR->dim(1)-1 ) {
		$metallicityHigh = "\infty";
	    } else {
		$metallicityHigh = FormatSigFigs($metallicities->index($i),3);
	    }
	    $label = "Spheroid: ".$metallicityLow."<Z<".$metallicityHigh;
	    $plotCommand .= $separator."'-' pt 4 title \"".$label."\"";
	    $separator = ", ";
	}
    }
    print gnuPlot $plotCommand."\n";
    for($i=0;$i<$diskSFR->dim(1);++$i) {
	if (any($diskSFR(:,($i)) > $sfrMinimum)) {
	    for($j=0;$j<$diskTimeCentral->nelem;++$j) {
		print gnuPlot $diskTimeCentral->index($j)." ".$diskSFR(($j),($i))."\n" unless ( $diskSFR(($j),($i)) <= $sfrMinimum );
	    }
	    print gnuPlot "e\n";
	}
    }
    for($i=0;$i<$spheroidSFR->dim(1);++$i) {
	if (any($spheroidSFR(:,($i)) > $sfrMinimum)) {
	    for($j=0;$j<$spheroidTimeCentral->nelem;++$j) {
		print gnuPlot $spheroidTimeCentral->index($j)." ".$spheroidSFR(($j),($i))."\n" unless ( $spheroidSFR(($j),($i)) <= $sfrMinimum );
	    }
	    print gnuPlot "e\n";
	}
    }
    close(gnuPlot);
    
    # Convert to PDF.
    system("ps2pdf tmp.ps ".$plotFile);
    
    # Clean up files.
    unlink("tmp.ps");
}

# Write disk SFR data.
print gHndl "#\n";
print gHndl "# Disk star formation rate:\n";
print gHndl "# Time [Gyr] | SFR [M_Solar/yr]\n";
$table = Text::Table->new();
for($j=0;$j<$diskTimeCentral->nelem;++$j) {
    undef(@rowData);
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
for($j=0;$j<$spheroidTimeCentral->nelem;++$j) {
    undef(@rowData);
    $rowData[0] = $spheroidTimeCentral->index($j);
    push(@rowData,$spheroidSFR(($j),:)->list);
    $table->add(@rowData);
}
print gHndl $table;

exit;

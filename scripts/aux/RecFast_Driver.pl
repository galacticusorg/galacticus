#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use File::Copy;
use Data::Dumper;
use DateTime;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Download, compile and run RecFast.
# Andrew Benson (18-January-2011)

# Get arguments.
die "Usage: RecFast_Driver.pl <parameterFile> <outputFile>"
    unless ( scalar(@ARGV) == 2 );
my $parameterFile  = $ARGV[0];
my $outputFileName = $ARGV[1];

# Parse the parameter file.
my $xml = new XML::Simple;
my $parameterData = $xml->XMLin($parameterFile);
unlink($parameterFile);

# Extract current file format version.
my $fileFormat        = $parameterData->{'fileFormat'}->{'value'};
my $fileFormatCurrent = 1;
die('RecFast_Driver.pl: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
    unless ( $fileFormat == $fileFormatCurrent );


# Download the code.
unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/recfast.for" ) {
    print "RecFast_Driver.pl: downloading RecFast code.\n";
    system("mkdir -p ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast; wget http://www.astro.ubc.ca/people/scott/recfast.for -O ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/recfast.for");
    die("RecFast_Driver.pl: FATAL - failed to download RecFast code.") unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/recfast.for" );
}

# Patch the code.
unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/patched" ) {
    print "RecFast_Driver.pl: patching RecFast code.\n";
    foreach my $file ( "recfast.for.patch" ) {
	copy($ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast_Galacticus_Modifications/".$file,$ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/".$file);
	if ( $file =~ m/\.patch$/ ) {system("cd ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast; patch < $file")};
	print "$file\n";
    }
    system("touch ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/patched");
}

# Build the code.
unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/recfast.exe" ) {
    print "RecFast_Driver.pl: compiling RecFast code.\n";
    system("cd ".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/; gfortran recfast.for -o recfast.exe -O3 -ffixed-form -ffixed-line-length-none");
    die("RecFast_Driver.pl: FATAL - failed to build RecFast code.") unless ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/recfast.exe" );
}

# Run the RecFast code.
my $buildFile = 0;
system("mkdir -p `dirname ".$outputFileName."`");
if ( -e $outputFileName ) {
    my $existingFile      = new PDL::IO::HDF5($outputFileName);
    (my $previousVersion) = $existingFile->attrGet('fileFormat');
    $buildFile = 1 
	if ( $previousVersion != $fileFormatCurrent );
} else {
    $buildFile = 1;
}
if ( $buildFile == 1 ) {
    # Open output data file.
    my $outputFile = new PDL::IO::HDF5(">".$outputFileName);
    # Check that required parameters exist.
    my @parameters = ( "OmegaBaryon", "OmegaMatter", "OmegaDarkEnergy", "HubbleConstant", "temperatureCMB", "Y_He" );
    my $parametersGroup = $outputFile->group('Parameters');
    foreach my $parameter ( @parameters ) {
	die("CMBFast_Driver.pl: FATAL - parameter ".$parameter." can not be found.") unless ( exists($parameterData->{$parameter}) );
	$parametersGroup->attrSet($parameter => pdl $parameterData->{$parameter}->{'value'});
    }
    # Extract variables.
    my $OmegaB = $parameterData->{'OmegaBaryon'    }->{'value'};
    my $OmegaM = $parameterData->{'OmegaMatter'    }->{'value'};
    my $OmegaL = $parameterData->{'OmegaDarkEnergy'}->{'value'};
    my $H0     = $parameterData->{'HubbleConstant' }->{'value'};
    my $T0     = $parameterData->{'temperatureCMB' }->{'value'};
    my $Yp     = $parameterData->{'Y_He'           }->{'value'};
    # Compute derived quantities.
    my $OmegaDM = $OmegaM-$OmegaB;
    # Drive RecFast.
    my $recfastOutput = "recFastOutput.data";
    open(pHndl,"|".$ENV{'GALACTICUS_EXEC_PATH'}."/aux/RecFast/recfast.exe");
    print pHndl $recfastOutput."\n";
    print pHndl $OmegaB." ".$OmegaDM." ".$OmegaL."\n";
    print pHndl $H0." ".$T0." ".$Yp."\n";
    print pHndl "1\n";
    print pHndl "6\n";
    close(pHndl);
    # Parse the output file.
    my @redshift          = ();
    my @electronFraction  = ();
    my @hIonizedFraction  = ();
    my @heIonizedFraction = ();
    my @matterTemperature = ();
    open(iHndl,$recfastOutput);
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/^\s*\d/ ) {
	    $line =~ s/^\s*//;
	    $line =~ s/\s*$//;
	    my @columns = split(/\s+/,$line);
	    push(@redshift         ,$columns[0]);
	    push(@electronFraction ,$columns[1]);
	    push(@hIonizedFraction ,$columns[2]);
	    push(@heIonizedFraction,$columns[3]);
	    push(@matterTemperature,$columns[4]);
	}
    }
    close(iHndl);
    unlink($recfastOutput);
    # Write arrays to the output file.
    $outputFile->dataset('redshift'         )->set(pdl @redshift         );
    $outputFile->dataset('electronFraction' )->set(pdl @electronFraction );
    $outputFile->dataset('hIonizedFraction' )->set(pdl @hIonizedFraction );
    $outputFile->dataset('heIonizedFraction')->set(pdl @heIonizedFraction);
    $outputFile->dataset('matterTemperature')->set(pdl @matterTemperature);
    # Add units data to output structure.
    $outputFile->dataset('matterTemperature')->attrSet(units     => "Kelvin");
    $outputFile->dataset('matterTemperature')->attrSet(unitsInSI => pdl 1.0 );
    # Add description and provenance to output structure.
    $outputFile->attrSet(description => "IGM ionization/thermal state computed using RecFast");
    my $date = `date`;
    chomp($date);
    my $provenanceGroup = $outputFile->group('provenance');
    $provenanceGroup->attrSet(date   => $date                   );
    $provenanceGroup->attrSet(source => "Galacticus via RecFast");
    my $version = "unknown";
    open(iHndl,"aux/RecFast/recfast.for");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/^CV\s+Version:\s+([\d\.]+)\s*$/ ) {$version = $1};
    }
    close(iHndl);
    my $recFastProvenance = $provenanceGroup->group('recFast');
    $recFastProvenance->attrSet(version => $version);
    $recFastProvenance->attrSet(notes   => "Includes modification of H recombination.\nIncludes all modifications for HeI recombination");
    # Add file format.
    $outputFile->attrSet(fileFormat => pdl long($fileFormatCurrent));
    # Add timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
    $outputFile->attrSet(timeStamp => $now);
}

exit;

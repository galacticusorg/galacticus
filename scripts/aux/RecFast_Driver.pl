#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
use File::Copy;
use Data::Dumper;
use DateTime;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V092'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V092'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");

# Download, compile and run RecFast.
# Andrew Benson (18-January-2011)

# Get arguments.
if ( $#ARGV != 1 ) {die "Usage: RecFast_Driver.pl <parameterFile> <outputFile>"};
my $parameterFile = $ARGV[0];
my $outputFile    = $ARGV[1];

# Parse the parameter file.
my $xml = new XML::Simple;
my $data = $xml->XMLin($parameterFile);
my $parameterHash = $data->{'parameter'};
unlink($parameterFile);

# Declare data structure.
my $output;

# Check that required parameters exist.
my @parameters = ( "Omega_b", "Omega_Matter", "Omega_DE", "H_0", "T_CMB", "Y_He" );
foreach my $parameter ( @parameters ) {
    die("CMBFast_Driver.pl: FATAL - parameter ".$parameter." can not be found.") unless ( exists($data->{'parameter'}->{$parameter}) );
    $output->{'provenance'}->{'recFast'}->{'parameters'}->{$parameter} = $data->{'parameter'}->{$parameter};
}

# Extract variables.
my $OmegaB = $parameterHash->{'Omega_b'      }->{'value'};
my $OmegaM = $parameterHash->{'Omega_Matter' }->{'value'};
my $OmegaL = $parameterHash->{'Omega_DE'     }->{'value'};
my $H0     = $parameterHash->{'H_0'          }->{'value'};
my $T0     = $parameterHash->{'T_CMB'        }->{'value'};
my $Yp     = $parameterHash->{'Y_He'         }->{'value'};

# Extract current file format version.
my $fileFormat        = $parameterHash->{'fileFormat'}->{'value'};
my $fileFormatCurrent = 1;
die('RecFast_Driver.pl: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
    unless ( $fileFormat == $fileFormatCurrent );

# Compute derived quantities.
my $OmegaDM = $OmegaM-$OmegaB;

# Download the code.
unless ( -e $galacticusPath."aux/RecFast/recfast.for" ) {
    print "RecFast_Driver.pl: downloading RecFast code.\n";
    system("mkdir -p ".$galacticusPath."aux/RecFast; wget http://www.astro.ubc.ca/people/scott/recfast.for -O ".$galacticusPath."aux/RecFast/recfast.for");
    die("RecFast_Driver.pl: FATAL - failed to download RecFast code.") unless ( -e $galacticusPath."aux/RecFast/recfast.for" );
}

# Patch the code.
unless ( -e $galacticusPath."aux/RecFast/patched" ) {
    print "RecFast_Driver.pl: patching RecFast code.\n";
    foreach my $file ( "recfast.for.patch" ) {
	copy($galacticusPath."aux/RecFast_Galacticus_Modifications/".$file,$galacticusPath."aux/RecFast/".$file);
	if ( $file =~ m/\.patch$/ ) {system("cd ".$galacticusPath."aux/RecFast; patch < $file")};
	print "$file\n";
    }
    system("touch ".$galacticusPath."aux/RecFast/patched");
}

# Build the code.
unless ( -e $galacticusPath."aux/RecFast/recfast.exe" ) {
    print "RecFast_Driver.pl: compiling RecFast code.\n";
    system("cd ".$galacticusPath."aux/RecFast/; gfortran recfast.for -o recfast.exe -O3 -ffixed-form -ffixed-line-length-none");
    die("RecFast_Driver.pl: FATAL - failed to build RecFast code.") unless ( -e $galacticusPath."aux/RecFast/recfast.exe" );
}

# Run the RecFast code.
my $buildFile = 0;
system("mkdir -p `dirname ".$outputFile."`");
if ( -e $outputFile ) {
    my $xmlFile         = new XML::Simple;
    my $previousFile    = $xmlFile->XMLin($outputFile);
    my $previousVersion = $previousFile->{'fileFormat'};
    $buildFile = 1 if ( $previousVersion != $fileFormatCurrent );
} else {
    $buildFile = 1;
}
if ( $buildFile == 1 ) {
    my $recfastOutput = "recFastOutput.data";
    open(pHndl,"|".$galacticusPath."aux/RecFast/recfast.exe");
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
    
    # Add arrays to output structure.
    @{$output->{'redshift'         }->{'datum'}} = @redshift         ;
    @{$output->{'electronFraction' }->{'datum'}} = @electronFraction ;
    @{$output->{'hIonizedFraction' }->{'datum'}} = @hIonizedFraction ;
    @{$output->{'heIonizedFraction'}->{'datum'}} = @heIonizedFraction;
    @{$output->{'matterTemperature'}->{'datum'}} = @matterTemperature;
    
    # Add units data to output structure.
    $output->{'matterTemperature'}->{'units'    } = "Kelvin";
    $output->{'matterTemperature'}->{'unitsInSI'} = 1.0;

    # Add description and provenance to output structure.
    $output->{'description'} = "IGM ionization/thermal state computed using RecFast";
    my $date = `date`;
    chomp($date);
    $output->{'provenance'}->{'date'} = $date;
    $output->{'provenance'}->{'source'} = "Galacticus via RecFast";
    my $version = "unknown";
    open(iHndl,"aux/RecFast/recfast.for");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/^CV\s+Version:\s+([\d\.]+)\s*$/ ) {$version = $1};
    }
    close(iHndl);
    $output->{'provenance'}->{'recFast'}->{'version'} = $version;
    @{$output->{'provenance'}->{'recFast'}->{'note'}} = ( 
	"Includes modification of H recombination",
	"Includes all modifications for HeI recombination"
	);
    
    # Add file format.
    $output->{'fileFormat'} = $fileFormatCurrent;

    # Add timestamp.
    my $dt = DateTime->now->set_time_zone('local');
    (my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
    my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
    $output->{'timeStamp'} = $now;

    # Output as XML.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"igm");
    open(outHndl,">".$outputFile);
    print outHndl $xmlOutput->XMLout($output);
    close(outHndl);
}

exit;

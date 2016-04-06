#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
use XML::Simple;
require Galacticus::Options;
require Galacticus::Constraints::Parameters;

# Collect systematic shifts arising from model discrepancies and combine them into offsets on systematics model parameter priors.
# Andrew Benson (01-April-2016)

# Get arguments.
die("Usage: collectSystematics.pl <configFile> <configFileNew> [options]")
    unless ( scalar(@ARGV) >= 2 );
my $configFileName    = $ARGV[0];
my $configFileNameNew = $ARGV[1];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Parse the configuration file.
my $config = &Parameters::Parse_Config($configFileName,useStored => 0);

# Validate the config file.
die("monteCarloTrees.pl: workDirectory must be specified in config file" )
    unless ( exists($config->{'likelihood'}->{'workDirectory' }) );

# Determine the work directory.
my $workDirectory = $config->{'likelihood'}->{'workDirectory'};

# Initialize mass shift systematics.
my $massShifts;

# Search for model discrepancies.
my $discrepanciesDirectoryName = $workDirectory."/modelDiscrepancy/";
opendir(my $discrepanciesDirectory,$discrepanciesDirectoryName);
while ( my $discrepancyName = readdir($discrepanciesDirectory) ) {
    next
	if ( $discrepancyName =~ m/^\./ );
    print " -> Processing discrepancies for: ".$discrepancyName."\n";
    # Open the discrepancy directory and search for discrepancy files.
    opendir(my $discrepancy,$discrepanciesDirectoryName.$discrepancyName);
    while ( my $discrepancyFileName = readdir($discrepancy) ) {
	next
	    unless ( $discrepancyFileName =~ m/\.hdf5$/ );
	print "    -> Processing discrepancy for: ".$discrepancyFileName."\n";
	(my $constraintName = $discrepancyFileName) =~ s/^discrepancy(.*)\.hdf5$/$1/;
	$constraintName = lcfirst($constraintName);
	# Open the file and check for any systematics.
	my $discrepancyFile = new PDL::IO::HDF5($discrepanciesDirectoryName.$discrepancyName."/".$discrepancyFileName);
	next
	    unless ( grep {$_ eq "systematicModels"} $discrepancyFile->groups() );
	# Iterate over available systematics.
	my $systematicModels = $discrepancyFile->group('systematicModels');
	foreach my $systematicType ( $systematicModels->groups() ) {
	    print "       -> Processing systematic: ".$systematicType."\n";
	    my $systematic = $systematicModels->group($systematicType);
	    # Handle individual known systematic types.
	    if ( $systematicType eq "MassShift" ) {
		# A simple systematic model for shifts in mass.
		my @attributes = $systematic->attrs();
		my @coefficientList;
		foreach my $attribute ( @attributes ) {
		    if ( $attribute =~ m/^systematicMassShiftCoefficient(\d+)$/ ) {
			my $order = $1;
			($coefficientList[$order]) = $systematic->attrGet($attribute);
		    }
		}
		my $coefficients = pdl @coefficientList;
		if ( exists($massShifts->{$constraintName}) ) {
		    $massShifts->{$constraintName} += $coefficients;
		} else {
		    $massShifts->{$constraintName}  = $coefficients;
		}
	    }
	}
    }
    close($discrepancy);
}
close($discrepanciesDirectory);

# Apply effects of systematics models.
# MassShift systematics.
my %massShiftNameMappings =
    (
     "alfalfaHiMassFunctionZ0.00" => "hiMassFunctionZ0.00"
    );
foreach my $parameter ( @{$config->{'parameters'}->{'parameter'}} ) {
    if ( $parameter->{'name'} =~ m/^(.*)MassFunctionZ([\d\.]+)MassSystematic(\d+)$/ ) {
	# Construct the name of the constraint.
	my $sampleName     = $1;
	my $redshift       = $2;
	my $order          = $3;
	my $constraintName = $sampleName."MassFunctionZ".$redshift;
	$constraintName = $massShiftNameMappings{$constraintName}
	    if ( exists($massShiftNameMappings{$constraintName}) );
	# Check for existance of a mass shift for this systematic.
	if ( exists($massShifts->{$constraintName}) ) {
	    if ( exists($parameter->{'define'}) ) {
		print "Applying shift to parameter definition for : ".$parameter->{'name'}." of ".$massShifts->{$constraintName}->(($order))."\n";	
		$parameter->{'define'} = $massShifts->{$constraintName}->(($order))."+".$parameter->{'define'};
	    } elsif ( exists($parameter->{'prior'}->{'distribution'}) ) {
		print "Applying shift to parameter prior for : ".$parameter->{'name'}." of ".$massShifts->{$constraintName}->(($order))."\n";
		if ( $parameter->{'prior'}->{'distribution'}->{'type'} eq "normal" ) {
		    my $shift = $massShifts->{$constraintName}->(($order))->sclr();
		    foreach ( "mean", "minimum", "maximum" ) {
			if ( exists($parameter->{'prior'}->{'distribution'}->{$_}) ) {
			    my $applyShift = 
				(
				 ( $_ eq "maximum" && $shift < 0.0 )
				 ||
				 ( $_ eq "minimum" && $shift > 0.0 )
				)
				?
				0
				:
				1;
			    $parameter->{'prior'}->{'distribution'}->{$_} += $shift
				if ( $applyShift );
			}
		    }
		} else {
		    die('collectSystematics.pl: do not know how to handle this prior distribution');
		}
	    } else {
		print Dumper($parameter);
		die('collectSystematics.pl: parameter has no prior nor is it "defined"');
	    }
	}	
    }
}

# Output the modified config file.
open(my $configFileNew,">".$configFileNameNew);
my $xml = new XML::Simple();
print $configFileNew $xml->XMLout($config, NoAttr => 1, KeyAttr => [], RootName => "simulationConfig");
close($configFileNew);

exit;

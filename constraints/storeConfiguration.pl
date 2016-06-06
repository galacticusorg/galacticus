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
use Storable;
use XML::Simple;
use Data::Dumper;
require Galacticus::Constraints::Parameters;

# Parse constraints configuration and store in a fast binary format for rapid reuse during MCMC simulation.
# Andrew Benson (20-March-2014)

# Get arguments.
die("Usage: storeConfiguration.pl <configFile>")
    unless ( scalar(@ARGV) == 1 );
my $configFile = $ARGV[0];

# Parse the configuration file.
my $config = &Parameters::Parse_Config($configFile,useStored => 0);

# Store the configuration file.
store($config,$configFile.".store");

# Get an XML worker object.
my $xml = new XML::Simple();

# Determine the compilation file in use.
my $compilationFileName = "constraints/compilations/".$config->{'likelihood'}->{'compilation'};

# Determine the base parameters file in use.
my $parametersFileName = $config->{'likelihood'}->{'baseParameters'};

# Parse and store the compilation file.
my $compilation = $xml->XMLin($compilationFileName);
store($compilation,$compilationFileName.".store");

# Parse and store the base parameters file.
my $parameters = $xml->XMLin($parametersFileName);
store($parameters,$parametersFileName.".store");

# Scan through all constraints.
my @constraints;
if ( ref($compilation->{'constraint'}) eq "ARRAY" ) {
    @constraints = @{$compilation->{'constraint'}};
} else {
    push(@constraints,$compilation->{'constraint'});
}
foreach my $constraint ( @constraints ) {
 my $constraintDefinition = $xml->XMLin($constraint->{'definition'},KeyAttr => "");
 store($constraintDefinition,$constraint->{'definition'}.".store");
}

exit;

#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use XML::SAX::ParserFactory;
use XML::Validator::Schema;
use XML::Simple;
use Scalar::Util 'reftype';
require List::ExtraUtils;

# Validate a Galacticus XML parameter file.
# Andrew Benson (1-December-2013)

# Get command line arguments.
die("Usage: validateParameters.pl <file>")
    unless ( scalar(@ARGV) == 1 );
my $file = $ARGV[0];

# Validate the parameter file using XML schema.
my $validator = XML::Validator::Schema->new(file => 'schema/parameters.xsd');
my $parser    = XML::SAX::ParserFactory->parser(Handler => $validator); 
eval { $parser->parse_file($file) };
die "Parameter file fails XML schema validation\n".$@
    if $@;

# Parse the file.
my $xml = new XML::Simple;
my $parameters = $xml->XMLin($file, KeyAttr => "");

# Check for duplicated entries.
my $valid = 0;
my %names;
++$names{$_->{'name'}}
  foreach ( &ExtraUtils::as_array($parameters->{'parameter'}) );
foreach ( keys(%names) ) {
    if ( $names{$_} > 1 ) {
	$valid = 1;
	print "Parameter '".$_."' appears ".$names{$_}." times - should appear only once\n";
    }
}
exit $valid;

#!/usr/bin/env perl
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use strict;
use warnings;
use XML::SAX::ParserFactory;
use XML::Validator::Schema;
use XML::Simple;
use Scalar::Util 'reftype';
use Data::Dumper;
use List::ExtraUtils;

# Validate a Galacticus XML parameter file.
# Andrew Benson (1-December-2013)

# Get command line arguments.
die("Usage: validateParameters.pl <file>")
    unless ( scalar(@ARGV) == 1 );
my $file = $ARGV[0];
# Parse the file.
my $xml = new XML::Simple;
my $parameters = $xml->XMLin($file, KeyAttr => "");
# Determine the format.
my $format = 2; # Best guess if not other information.
if ( exists($parameters->{'formatVersion'}) ) {
    # Use the format declared in the file itself.
    $format = $parameters->{'formatVersion'};
} elsif ( exists($parameters->{'parameter'}) ) {
    # Does it look like format version 1?
    $format = 1;
}
# Assume valid by default.
my $valid = 0;
# Handle different file formats.
if ( $format == 1 ) {
    # Handle format version 1.
    # Validate the parameter file using XML schema.
    my $validator = XML::Validator::Schema->new(file => $ENV{'GALACTICUS_EXEC_PATH'}.'/schema/parameters.xsd');
    my $parser    = XML::SAX::ParserFactory->parser(Handler => $validator); 
    eval { $parser->parse_file($file) };
    die "Parameter file fails XML schema validation\n".$@
	if $@;
    # Check for duplicated entries.
    my %names;
    ++$names{$_->{'name'}}
        foreach ( &List::ExtraUtils::as_array($parameters->{'parameter'}) );
    foreach ( keys(%names) ) {
	if ( $names{$_} > 1 ) {
	    $valid = 1;
	    print "Parameter '".$_."' appears ".$names{$_}." times - should appear only once\n";
	}
    }
} elsif ( $format == 2 ) {
    # Handle format version 2.
    my $hasFormat;
    my $hasVersion;
    my $hasLastModified;
   # Create a stack of elements to check.
    my @stack = map {{name => $_, node => $parameters->{$_}}} keys(%{$parameters});
    while ( scalar(@stack) > 0 ) {
	my $element = pop(@stack);
	if ( $element->{'name'} eq 'formatVersion' ) {
	    $hasFormat = 1;
	} elsif ( $element->{'name'} eq 'version' ) {
	    $hasVersion = 1;
	} elsif ( $element->{'name'} eq 'lastModified' ) {
	    $hasLastModified = 1;
	} elsif ( reftype($element->{'node'}) eq "ARRAY" ) {
	    if ( $element->{'name'} =~ m/\-\>/ ) { # Duplicates are allowed only in subparameters.
		foreach my $node ( @{$element->{'node'}} ) {
		    if ( ! exists($node->{'value'}) ) {
			unless ( %{$node} ) {
			    $valid = 1;
			    print "Parameter '".$element->{'name'}."' has no value\n";
			}
		    } elsif ( reftype($node->{'value'}) && reftype($node->{'value'}) eq "ARRAY" ) {
			$valid = 1;
			print "Parameter '".$element->{'name'}."' has multiple values\n";
		    }
		    push
			(
			 @stack,
			 map {$_ eq "value" || ! reftype($node->{$_}) ? () : {name => $element->{'name'}."->".$_, node => $node->{$_}}} keys(%{$node})
			);
		}
	    } else {
		$valid = 1;
		print "Parameter '".$element->{'name'}."' appears ".scalar(@{$element->{'node'}})." times - should appear only once\n";
	    }
	} else {
	    if ( ! exists($element->{'node'}->{'value'}) ) {
		unless ( %{$element->{'node'}} ) {
		    $valid = 1;
		    print "Parameter '".$element->{'name'}."' has no value\n";
		}
	    } elsif ( reftype($element->{'node'}->{'value'}) && reftype($element->{'node'}->{'value'}) eq "ARRAY" ) {
		$valid = 1;
		print "Parameter '".$element->{'name'}."' has multiple values\n";
	    }
	    push
		(
		 @stack,
		 map {$_ eq "value" || ! reftype($element->{'node'}->{$_}) ? () : {name => $element->{'name'}."->".$_, node => $element->{'node'}->{$_}}} keys(%{$element->{'node'}})
		);
	}
    }
}

exit $valid;

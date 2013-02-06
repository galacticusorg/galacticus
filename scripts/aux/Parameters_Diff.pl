#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
use Scalar::Util qw(looks_like_number);
use Data::Dumper;

# Read the command line arguments.
if ( $#ARGV != 1 ) {die("Usage: Parameters_Diff.pl <file1> <file2>")};
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

# Get an XML object.
my $xml = new XML::Simple;

# Read the XML for both parameter files.
my $parameters1 = $xml->XMLin($file1);
my $parameters2 = $xml->XMLin($file2);

# Initially assume that files match.
my $filesMatch = 1;

# Loop over all parameters in the first file.
foreach my $parameter1 ( keys(%{$parameters1->{'parameter'}}) ) {
    # Test if the parameter exists in the second file.
    unless ( exists($parameters2->{'parameter'}->{$parameter1}) ) {
	# It does not.
	print "Parameter ".$parameter1." in file 1 is not present in file 2.\n";
	$filesMatch = 0;
    } else {
	# It does, so test for equality. Have to figure out if the values are numeric or strings.
	my $areEqual;
	if ( looks_like_number($parameters1->{'parameter'}->{$parameter1}->{'value'}) == 0 ) {
	    # First value is a string.
	    if ( looks_like_number($parameters2->{'parameter'}->{$parameter1}->{'value'}) == 0 ) {
		# Both values are strings - test equality.
		if ( $parameters1->{'parameter'}->{$parameter1}->{'value'} eq $parameters2->{'parameter'}->{$parameter1}->{'value'} ) {
		    $areEqual = 1;
		} else {
		    $areEqual = 0;
		}
	    } else {
		# Second value is numeric - the values can't be equal.
		$areEqual = 0;
	    }
	} else {
	    # First value is numeric.
	    if ( looks_like_number($parameters2->{'parameter'}->{$parameter1}->{'value'}) != 0 ) {
		# Both values are numeric - test for equality.
		if ( $parameters1->{'parameter'}->{$parameter1}->{'value'} == $parameters2->{'parameter'}->{$parameter1}->{'value'} ) {
		    $areEqual = 1;
		} else {
		    $areEqual = 0;
		}
	    } else {
		# Second value is not numeric - the two values can't be equal.
		$areEqual = 0;
	    }
	}
	# Report unequal values.
	unless ( $areEqual == 1 ) {
	    print "Values for ".$parameter1." differ:\n";
	    print " -> File 1: ".$parameters1->{'parameter'}->{$parameter1}->{'value'}."\n";
	    print " -> File 2: ".$parameters2->{'parameter'}->{$parameter1}->{'value'}."\n";
	    $filesMatch = 0;
	}
    }
}

# Loop over all parameters in the second file.
foreach my $parameter2 ( keys(%{$parameters2->{'parameter'}}) ) {
    # Check that the parameter also exists in the first file.
    unless ( exists($parameters1->{'parameter'}->{$parameter2}) ) {
	print "Parameter ".$parameter2." in file 2 is not present in file 1.\n";
	$filesMatch = 0;
    }
}

# Write final report.
if ( $filesMatch == 1 ) {
    print "Files match.\n";
} else {
    print "Files do not match.\n";
}

exit;

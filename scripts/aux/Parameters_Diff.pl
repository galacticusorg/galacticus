#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
use Scalar::Util qw(looks_like_number reftype);
use Data::Dumper;

# Read the command line arguments.
die("Usage: Parameters_Diff.pl <file1> <file2>")
    unless ( scalar(@ARGV) == 2 );
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

# Get an XML object.
my $xml = new XML::Simple;

# Read the XML for both parameter files.
my $parameters1 = $xml->XMLin($file1);
my $parameters2 = $xml->XMLin($file2);

# Do the comparison.
my $filesMatch = &Compare($parameters1,$parameters2,1,0) && &Compare($parameters2,$parameters1,0,1);

# Write final report.
if ( $filesMatch == 1 ) {
    print "Files match.\n";
} else {
    print "Files do not match.\n";
}

exit;

sub Compare {
    my $parameters1   = shift();
    my $parameters2   = shift();
    my $checkValues   = shift();
    my $reverseLabels = shift();
    # File labels.
    my $label1 = "1";
    my $label2 = "2";
    if ( $reverseLabels ) {
	$label1 = "2";
	$label2 = "1";
    }
    # Assume files match.
    my $filesMatch = 1;
    # Iterate over all parameters in the first file.
    my @parameters1Stack = map {{name => $_, parameter => $parameters1->{$_}}} keys(%{$parameters1});
    while ( scalar(@parameters1Stack) > 0 ) {
	my $parameter1 = pop(@parameters1Stack);
	next
	    unless ( reftype($parameter1->{'parameter'}) );
	# Test if the parameter exists in the second file.
	my $parameter2 = 
	{
	    name      => $parameter1->{'name'},
	    parameter => $parameters2
	};
	my $searchPath = $parameter1->{'name'};
	while ( $searchPath ne "" ) {
	    my $searchRoot;
	    if ( $searchPath =~ m/([^\-\>]*)\-\-\>(.*)/ ) {
		$searchRoot = $1;
		$searchPath = $2;
	    } else {
		$searchRoot = $searchPath;
		$searchPath = "";
	    }
	    if ( exists($parameter2->{'parameter'}->{$searchRoot}) ) {
		$parameter2->{'parameter'} = $parameter2->{'parameter'}->{$searchRoot};
	    } else {
		undef($parameter2);
		last;
	    }
	}
	unless ( $parameter2 ) {
	    # It does not.
	    print "Parameter ".$parameter1->{'name'}." [value: ".$parameter1->{'parameter'}->{'value'}."] in file ".$label1." is not present in file ".$label2.".\n";
	    $filesMatch = 0;
	} elsif ( $checkValues ) {
	    # Check if there is a value for this parameter.
	    unless ( reftype($parameter1->{'parameter'}) eq "ARRAY" ) {
		# It does, so test for equality. Have to figure out if the values are numeric or strings.
		my $areEqual;
		if ( looks_like_number($parameter1->{'parameter'}->{'value'}) == 0 ) {
		    # First value is a string.
		    if ( looks_like_number($parameter2->{'parameter'}->{'value'}) == 0 ) {
			# Both values are strings - test equality.
			if ( $parameter1->{'parameter'}->{'value'} eq $parameter2->{'parameter'}->{'value'} ) {
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
		    if ( looks_like_number($parameter2->{'parameter'}->{'value'}) != 0 ) {
			# Both values are numeric - test for equality.
			if ( $parameter1->{'parameter'}->{'value'} == $parameter2->{'parameter'}->{'value'} ) {
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
		    print "Values for ".$parameter1->{'name'}." differ:\n";
		    if ( $reverseLabels ) {
			print " -> File 1: ".$parameter2->{'parameter'}->{'value'}."\n";
			print " -> File 2: ".$parameter1->{'parameter'}->{'value'}."\n";
		    } else {
			print " -> File 1: ".$parameter1->{'parameter'}->{'value'}."\n";
			print " -> File 2: ".$parameter2->{'parameter'}->{'value'}."\n";
		    }
		    $filesMatch = 0;
		}
	    }
	}
	# Push any subparameters onto the stack.
	if ( reftype($parameter1->{'parameter'}) eq "HASH" ) {
	    foreach ( keys(%{$parameter1->{'parameter'}}) ) {
		push
		    (
		     @parameters1Stack,
		     {
			 name      => $parameter1->{'name'     }."-->".$_ ,
			 parameter => $parameter1->{'parameter'}   -> {$_}
		     }
		    )
		    unless ( $_ eq "value" );
	    }
	} elsif  ( reftype($parameter1->{'parameter'}) eq "ARRAY" ) {
	    # Arrays are currently not handled.
	}
    }
    return $filesMatch;
}

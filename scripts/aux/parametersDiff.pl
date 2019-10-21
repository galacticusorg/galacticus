#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Scalar::Util qw(looks_like_number reftype);
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Options;

# Read the command line arguments.
die("Usage: parametersDiff.pl <file1> <file2>")
    unless ( scalar(@ARGV) >= 2 );
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my %options =
    (
     numericalTolerance => 0.0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Get an XML object.
my $xml = new XML::Simple;

# Read the XML for both parameter files.
my $parameters1 = $xml->XMLin($file1);
my $parameters2 = $xml->XMLin($file2);

# Do the comparison.
my $filesMatch = &Compare($parameters1,$parameters2,\%options);

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
    my %options       = %{shift()};
    # Assume files match.
    my $filesMatch = 1;
    # Initialize stack of subtrees.
    my @subTrees1 = ( $parameters1 );
    my @subTrees2 = ( $parameters2 );
    while ( scalar(@subTrees1) > 0 ) {
	# Pop subtrees off the stacks.
	my $subTree1 = pop(@subTrees1);
	my $subTree2 = pop(@subTrees2);
	my $stage = -1;
	foreach my $subTrees ( {tree1 => {tree => $subTree1, name => "file 1"}, tree2 => {tree => $subTree2, name => "file 2"}},  {tree2 => {tree => $subTree1, name => "file 1"}, tree1 => {tree => $subTree2, name => "file 2"}} ) {
	    ++$stage;
	    my $workTree1 = $subTrees->{'tree1'};
	    my $workTree2 = $subTrees->{'tree2'};
	    foreach my $parameterName ( keys(%{$workTree1->{'tree'}}) ) {
		# Skip non-parameters.
		next
		    unless ( reftype($workTree1->{'tree'}->{$parameterName}) );
		# Test for comparison.
		my $parametersMatch = 1;
		if ( exists($workTree2->{'tree'}->{$parameterName}) ) {
		    if ( $stage == 1 ) {
			# Get array of all parameters.
			my @workParameters1 = &List::ExtraUtils::as_array($workTree1->{'tree'}->{$parameterName});
			my @workParameters2 = &List::ExtraUtils::as_array($workTree2->{'tree'}->{$parameterName});
			if ( scalar(@workParameters1) != scalar(@workParameters2) ) {
			    # Number of parameters do not agree.
			    print "Parameter '".$parameterName."' array length differs from that  in ".$workTree2->{'name'}."\n";
			    $parametersMatch = 0;
			} else {
			    for(my $i=0;$i<scalar(@workParameters1);++$i) {
				# Parameter present in file 2. Check values.
				my $haveValue1 = exists($workParameters1[$i]->{'value'});
				my $haveValue2 = exists($workParameters2[$i]->{'value'});
				if ( $haveValue1 && $haveValue2 ) {
				    # Values are present.
				    my $isNumber1 = looks_like_number($workParameters1[$i]->{'value'});
				    my $isNumber2 = looks_like_number($workParameters2[$i]->{'value'});
				    if ( $isNumber1 && $isNumber2 ) {
					# Both values are numbers - do a numerical comparison.
					if ( ! &approximatelyEqual($workParameters1[$i]->{'value'},$workParameters2[$i]->{'value'},$options{'numericalTolerance'}) ) {
					    print "Parameter '".$parameterName."' value (".$workParameters1[$i]->{'value'}.") in ".$workTree1->{'name'}." differs from that (".$workParameters2[$i]->{'value'}.") in ".$workTree2->{'name'}."\n";
					    $parametersMatch = 0;
					}
				    } else {
					# One or both values are not numbers - do a string comparison.
					if ( $workParameters1[$i]->{'value'} ne $workParameters2[$i]->{'value'} ) {
					    print "Parameter '".$parameterName."' value (".$workParameters1[$i]->{'value'}.") in ".$workTree1->{'name'}." differs from that (".$workParameters2[$i]->{'value'}.") in ".$workTree2->{'name'}."\n";
					    $parametersMatch = 0;
					}
				    }
				} elsif ( $haveValue1 && ! $haveValue2 ) {
				    # Value missing in file 2.
				    print "Parameter '".$parameterName."' has value in ".$workTree1->{'name'}." but not in ".$workTree2->{'name'}."\n";
				    $parametersMatch = 0;
				} elsif ( ! $haveValue1 && $haveValue2 ) {
				    # Value missing in file 2.
				    print "Parameter '".$parameterName."' has no value in ".$workTree1->{'name'}." but has value in ".$workTree2->{'name'}."\n";
				    $parametersMatch = 0;
				}
			    }
			}
			# If parameters match then push any subtrees to the stack.
			if ( $parametersMatch ) {
			    for(my $i=0;$i<scalar(@workParameters1);++$i) {
				delete($workParameters1[$i]->{'value'});
				delete($workParameters2[$i]->{'value'});
			    }
			    push(@subTrees1,@workParameters1);
			    push(@subTrees2,@workParameters2);
			}
		    }
		} else {
		    # Parameter missing in tree 2.
		    print "Parameter '".$parameterName."' in ".$workTree1->{'name'}." not present in ".$workTree2->{'name'}."\n";
		    $parametersMatch = 0;
		}
		# If parameters do not match, files do not match.
		$filesMatch = 0
		    unless ( $parametersMatch );
	    }
	}
    }
    return $filesMatch;
}

sub approximatelyEqual {
    my $number1           = shift();
    my $number2           = shift();
    my $toleranceRelative = shift();
    my $toleranceAbsolute = 0.5*$toleranceRelative*(abs($number1)+abs($number2));
    return abs($number1-$number2) <= $toleranceAbsolute;
}

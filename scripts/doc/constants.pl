#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use File::Find;
use List::ExtraUtils;
use Scalar::Util qw(looks_like_number);

# Format constant definitions for incorporation into the documentation.
# Andrew Benson (29-January-2025)

# Get arguments.
die("Usage: constants.pl <buildPath> <outputFile>")
    unless ( scalar(@ARGV) == 2 );
my $buildPath  = $ARGV[0];
my $outputFile = $ARGV[1];

# Initialize a list to contain constants.
my @constants;

# Process all files.
my @directories = ( $buildPath );
find(\&process,@directories);

# Define known groups.
my %knownGroups =
    (
     "astrophysical" => "Astrophysical constants",
     "atomic"        => "Atomic physics constants",
     "physical"      => "Physical constants",
     "math"          => "Mathematical constants",
     "units"         => "Units",
     "prefixes"      => "SI Prefixes",
     "GSL"           => "Gnu Scientific Library constants",
     "Kernel"        => "Kernel constants",
     "misc"          => "Miscellaneous constants"
    );

# Find all groups.
my $groups;
foreach my $constant ( @constants ) {
    next
	unless ( exists($constant->{'value'}) );
    if ( exists($constant->{'group'}) ) {
	my @groupNames = split(":",$constant->{'group'});
	foreach my $groupName ( @groupNames ) {
	    die("Unknown group name '".$groupName."'")
		unless ( exists($knownGroups{$groupName}) );
	    push(@{$groups->{$groupName}->{'constant'}},$constant);
	}
    } else {
	push(@{$groups->{'misc'}->{'constant'}},$constant);
    }
}

# Iterate over groups and constants, formatting into LaTeX.
open(my $output,">",$outputFile);
foreach my $groupName ( sort(keys(%{$groups})) ) {
    print $output "\\subsection{".$knownGroups{$groupName}."}\n";
    my @constants = sort {$a->{'variable'} cmp $b->{'variable'}} @{$groups->{$groupName}->{'constant'}};
    foreach my $constant ( @constants ) {
	next
	    unless ( exists($constant->{'value'}) );
	my $reference;
	if ( exists($constant->{'referenceURL'}) ) {
	    $reference = "\\href{".$constant->{'referenceURL'}."}{".$constant->{'reference'}."}";
	} else {
	    $reference = $constant->{'reference'};
	}
	$reference =~ s/&/\\&/;
	(my $variable = $constant->{'variable'}) =~ s/_/\\_/g;
	my $external  = exists($constant->{'externalDescription'}) ? " (See \\href{".$constant->{'externalDescription'}."}{here}.)" : "";
	my $symbol    = exists($constant->{'symbol'}) ? "\$".$constant->{'symbol'}."\$: " : "";
	my $units     = exists($constant->{'units' }) && $constant->{'units'} ne "dimensionless" ? " [".$constant->{'units'}."]" : "";
	(my $value    = $constant->{'value'}) =~ s/d([\+\-0-9]+)$/e$1/;
	$value =~ s/_[_a-zA-Z0-9]+$//;
	(my $module   = $constant->{'module'}) =~ s/_/\\_/g;
	(my $fileName = $constant->{'fileName'}) =~ s/\./_/g;
	my $moduleURL = "https://github.com/galacticusorg/galacticus/releases/download/masterRelease/Galacticus_Source.pdf\\#source.".$fileName.":".lc($constant->{'module'});
	print $output "\\noindent {\\normalfont \\ttfamily ".$variable." = ".$value."}".$units."\\\\\n";
	print $output "\\indent ".$symbol.$constant->{'description'}.$external."\\\\\n";
	print $output "\\indent Module: \\href{".$moduleURL."}{\\normalfont \\ttfamily ".$module."}\\\\\n";
	print $output "\\indent Reference: ".$reference."\n\n\\medskip\n";
    }
}
close($output);

exit;

sub process {
    # Process constant definition files, accumulating all constants.
    my $fileName = $_;
    my $fullName = $File::Find::name;
    # Ignore files not containing constant definitions.
    return
	unless ( $fileName =~ m/\.constants.xml$/ );
    # Parse the file.
    my $xml               = new XML::Simple;
    my $constantsFromFile = $xml->XMLin($fileName);
    # Append to our list.
    push(@constants,&List::ExtraUtils::as_array($constantsFromFile->{'constant'}));
}




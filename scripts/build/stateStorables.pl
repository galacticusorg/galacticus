#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Galacticus::Build::Directives;
use Galacticus::Build::SourceTree;
use List::ExtraUtils;
use Data::Dumper;

# Builds a list of classes which support store/restore of their state.
# Andrew Benson (10-April-2018)

# Read command line arguments.
die "Usage: stateStoreables.pl <installDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $installDirectoryName = $ARGV[0];
# Get an XML parser.
my $xml                  = new XML::Simple();
# Parse the directive locations file.
my $directiveLocations   = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
# Initialize a structure for output of classes.
my $stateStorables;
# Find all files which contain functionClass objects - these all support state store/restore.
foreach my $functionClassFileName ( &List::ExtraUtils::as_array($directiveLocations->{'functionClass'}->{'file'}) ) {
    # Extract a functionClass directives from this file.
    foreach my $functionClass ( &Galacticus::Build::Directives::Extract_Directives($functionClassFileName,'functionClass') ) {
	push(@{$stateStorables->{'functionClasses'}},$functionClass->{'name'}."Class");
	# Find all files which contain instances of this functionClass.
	foreach my $instanceFileName ( &List::ExtraUtils::as_array($directiveLocations->{$functionClass->{'name'}}->{'file'}) ) {
	    # Extract the relevant directives from this file.
	    foreach my $instance ( &Galacticus::Build::Directives::Extract_Directives($instanceFileName,$functionClass->{'name'}) ) {
		push(@{$stateStorables->{'functionClassInstances'}},$instance->{'name'});
	    }
	}
    }   
}
# Find all files which contain functionClassType objects - these all support state store/restore.
foreach my $functionClassFileName ( &List::ExtraUtils::as_array($directiveLocations->{'functionClassType'}->{'file'}) ) {
    # Extract a functionClassType directives from this file.
    push(@{$stateStorables->{'functionClassTypes'}},map {{name => $_->{'name'}, file => $functionClassFileName}} &Galacticus::Build::Directives::Extract_Directives($functionClassFileName,'functionClassType'));
}
# Find all files which contain stateStorable objects - these explicitly support state store/restore.
foreach my $stateStorableFileName ( &List::ExtraUtils::as_array($directiveLocations->{'stateStorable'}->{'file'}) ) {
    # Parse the source of this file.
    my $tree = &Galacticus::Build::SourceTree::ParseFile($stateStorableFileName);
    # Walk the tree.
    my $node           = $tree;
    my $depth          = 0    ;
    my @directiveNodes        ;
    my %classes               ;
    while ( $node ) {
	# Capture stateStorable directives.
	push(@directiveNodes,$node)
	    if ( $node->{'type'} eq "stateStorable" );
	# Capture derived type definitions.
	if ( $node->{'type'} eq "type" ) {
	    # Parse class openers to find dependencies.
	    if ( $node->{'opener'} =~ m/^\s*type\s*(,\s*(abstract)\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*([a-z0-9_]+)\s*$/i ) {
		my $type     = $5;
		my $extends  = $3;
		my $abstract = defined($2);
		$classes{$type} =
		{
		     node     => $node   ,
		     extends  => $extends,
		     abstract => $abstract
		};
	    } else {
		die("stateStorables.pl: unable to parse type definition opener");
	    }
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
    # Process any stateSotrable directives.
    foreach my $directiveNode ( @directiveNodes ) {
	my $directive = $directiveNode->{'directive'};
	# Scan all known classes, finding all which derive from the base class.
	foreach my $className ( sort(keys(%classes)) ) {
	    my $matches         = 0;
	    my $parentClassName = $className;
	    while ( defined($parentClassName) ) {
		if ( $parentClassName eq $directive->{'class'} ) {
		    $matches = 1;
		    last;
		}
		$parentClassName = $classes{$parentClassName}->{'extends'};
	    }
	    push(@{$stateStorables->{'stateStorables'}},{type => $className, class => $directive->{'class'}})
		if ( $matches );
	}
    }
}
# Sort results.
@{$stateStorables->{'functionClasses'       }} = sort                                 @{$stateStorables->{'functionClasses'       }};
@{$stateStorables->{'functionClassInstances'}} = sort                                 @{$stateStorables->{'functionClassInstances'}};
@{$stateStorables->{'stateStorables'        }} = sort {$a->{'type'} cmp $b->{'type'}} @{$stateStorables->{'stateStorables'        }};
# Output the results.
open(my $outputFile,">".$ENV{'BUILDPATH'}."/stateStorables.xml");
print $outputFile $xml->XMLout($stateStorables, RootName => "storables");
close($outputFile);

exit;

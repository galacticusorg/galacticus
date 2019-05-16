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

# Builds a list of classes which support actions during deep copy.
# Andrew Benson (15-May-2019)

# Read command line arguments.
die "Usage: deepCopyActions.pl <installDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $installDirectoryName = $ARGV[0];
# Get an XML parser.
my $xml                  = new XML::Simple();
# Parse the directive locations file.
my $directiveLocations   = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
# Initialize a structure for output of classes.
my $deepCopyActions;
# Find all files which contain deepCopyActions objects - these explicitly support state store/restore.
foreach my $deepCopyActionsFileName ( &List::ExtraUtils::as_array($directiveLocations->{'deepCopyActions'}->{'file'}) ) {
    # Parse the source of this file.
    my $tree = &Galacticus::Build::SourceTree::ParseFile($deepCopyActionsFileName);
    # Walk the tree.
    my $node           = $tree;
    my $depth          = 0    ;
    my @directiveNodes        ;
    my %classes               ;
    while ( $node ) {
	# Capture deepCopyActions directives.
	push(@directiveNodes,$node)
	    if ( $node->{'type'} eq "deepCopyActions" );
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
		die("deepCopyActions.pl: unable to parse type definition opener");
	    }
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
    # Process any deepCopyActions directives.
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
	    push(@{$deepCopyActions->{'deepCopyActions'}},{type => $className, class => $directive->{'class'}})
		if ( $matches );
	}
    }
}
# Sort results.
@{$deepCopyActions->{'deepCopyActions'}} = sort {$a->{'type'} cmp $b->{'type'}} @{$deepCopyActions->{'deepCopyActions'}};
# Output the results.
open(my $outputFile,">".$ENV{'BUILDPATH'}."/deepCopyActions.xml");
print $outputFile $xml->XMLout($deepCopyActions, RootName => "deepCopyActions");
close($outputFile);

exit;

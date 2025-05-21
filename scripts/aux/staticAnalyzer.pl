#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::SourceTree;

# Perform static analysis of Fortran files.
# Andrew Benson (28-February-2023)

# Get the file to process.
die('Usage: staticAnalyzer.pl <fileName>')
    unless ( scalar(@ARGV) == 1 );
my $fileName = $ARGV[0];

# Parse the file.
my $tree = &Galacticus::Build::SourceTree::ParseFile($fileName);
# Walk the tree.
my $node   = $tree;
my $depth  = 0;
my $status = 0;
while ( $node ) {
    # Class/type pointers in derived types should be null initialized.
    if ( $node->{'type'} eq "declaration" ) {
	foreach my $declaration ( @{$node->{'declarations'}} ) {
	    if ( $node->{'parent'}->{'type'} eq "type" && ( $declaration->{'intrinsic'} eq "type" || $declaration->{'intrinsic'} eq "class" ) && grep {$_ eq "pointer"} @{$declaration->{'attributes'}} ) {
		for(my $i=0;$i<scalar(@{$declaration->{'variables'}});++$i) {
		    next
			if ( $declaration->{'variables'}->[$i] =~ m/=>null\(\)$/ );
		    (my $typeName = $node->{'parent'}->{'opener'}) =~ s/.*::\s*([a-zA-Z0-9_]+).*/$1/;
		    chomp($typeName);
		    print "Pointer variable '".$declaration->{'variableNames'}->[$i]."' in type '".$typeName."' in file '".$fileName."' is not null initialized\n";
		    $status = 1;
		}
	    }
	}
    }
    # Look for duplicated assignments in `constructorAssign` directives.
    if ( $node->{'type'} eq "constructorAssign" ) {
	my @variables = split(/\s*,\s*/,$node->{'directive'}->{'variables'});
	my %countAssignments;
	foreach my $variable ( @variables ) {
	    ++$countAssignments{$variable};
	}
	foreach my $variable ( keys(%countAssignments) ) {
	    if ( $countAssignments{$variable} > 1 ) {
		print "Duplicated assignment of '".$variable."' in `constructorAssign` directive in file '".$fileName."'\n";
		$status = 1;
	    }
	}
    }
    # Walk to the next node in the tree.
    $node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
}

exit $status;

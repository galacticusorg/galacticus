# Contains a Perl module which analyzes pre-processed source code for duplication of "use" statements.

package Galacticus::Build::SourceTree::Analyze::UseDuplication;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::analyzeHooks{'useDuplication'} = \&Analyze_UseDuplication;

sub Analyze_UseDuplication {
    # Get the tree.
    my $tree  = shift();
    # Walk the tree.
    my $node  = $tree;
    my $depth = 0;
    my $fileName;
    while ( $node ) {
	if ( $node->{'type'} eq "moduleUse" ) {
	    # Process only module uses in functions or subroutines.
	    if ( $node->{'parent'}->{'type'} eq "function" || $node->{'parent'}->{'type'} eq "subroutine" ) {
		# Walk up the hierarchy finding already used modules.
		my @usedModules = ();
		my $nodeParent = $node->{'parent'}->{'parent'};
		while ( $nodeParent ) {
		    my $nodeChild = $nodeParent->{'firstChild'};
		    while ( $nodeChild ) {
			if ( $nodeChild->{'type'} eq "moduleUse" ) {
			    push(@usedModules,keys(%{$nodeChild->{'moduleUse'}}));
			}
			$nodeChild = $nodeChild->{'sibling'};
		    }
		    $nodeParent = $nodeParent->{'parent'};
		}
		foreach my $moduleName ( keys(%{$node->{'moduleUse'}}) ) {
		    if ( grep {$moduleName eq $_} @usedModules ) {
			print "Warning: module '".$moduleName."' used in ".$node->{'parent'}->{'name'}."() was already used in container [file: ".$fileName."]\n";
		    }
		}
	    }
	} elsif ( $node->{'type'} eq "file" ) {
	    $fileName = $node->{'name'};
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

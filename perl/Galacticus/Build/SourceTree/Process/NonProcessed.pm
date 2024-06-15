# Contains a Perl module which implements allocation of variables with some extra, useful functionality.

package Galacticus::Build::SourceTree::Process::NonProcessed;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks       {'nonProcessed'} = \&Process_NonProcessed;
$Galacticus::Build::SourceTree::Hooks::processDependencies{'nonProcessed'} = [ "generics" ];

sub Process_NonProcessed {
    # Get the tree.
    my $tree = shift();
    # Non-processed directives that we simply mark as processed to avoid warnings.
    my @nonProcessedDirectives = ( "methods", "workaround", "include", "functionGlobal", "component", "radiusSolverPlausibility", "interTreePositionInsert", "expiry", "scoping" );
    # Walk the tree, looking for our directive.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if (
	    $node->{'type'} =~ m/Task$/
	    ||
	    grep {$node->{'type'} eq $_} @nonProcessedDirectives
	    ) {
	    # Record that node is processed.
	    $node->{'directive'}->{'processed'} = 1;
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

# Contains a Perl module which implements processing of object builder directives.

package ObjectBuilder;
use strict;
use warnings;
use utf8;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use Data::Dumper;
require List::ExtraUtils;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Hooks::processHooks{'objectBuilder'} = \&Process_ObjectBuilder;

sub Process_ObjectBuilder {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "objectBuilder" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate source code for the object builder. The logic here is that, if the passed
	    # parameter set is *not* the global set, and *does* contain a definition of the
	    # relevant class, then use that definition to build an object of the
	    # class. Otherwise, we want to use the default (global) implementation of the class,
	    # so we simply get a pointer to the default object. This avoids building the default
	    # implementation more than once.
	    my $builderCode;
	    $builderCode .= "if (".$node->{'directive'}->{'source'}."%isPresent('".$node->{'directive'}->{'class'}."Method').and..not.".$node->{'directive'}->{'source'}."%isGlobal()) then\n";
	    $builderCode .= "   ! ...construct the object from the provided parameters.\n";
	    $builderCode .= "   ".$node->{'directive'}->{'name'}." => ".$node->{'directive'}->{'class'}."(".$node->{'directive'}->{'source'}.")\n";
	    $builderCode .= "else\n";
	    $builderCode .= "   ! ...otherwise, use the default object.\n";
	    $builderCode .= "   ".$node->{'directive'}->{'name'}." => ".$node->{'directive'}->{'class'}."()\n";
	    $builderCode .= "end if\n";
	    # Build a code node.
	    my $newNode =
	    {
		type       => "code"      ,
		content    => $builderCode,
		firstChild => undef()
	    };
	    # Insert the node.
	    &SourceTree::InsertAfterNode($node,[$newNode]);
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	}
	# Walk to the next node in the tree.
	$node = &SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

# Contains a Perl module which implements reseting to begin a new sequence of deepCopy actions.

package Galacticus::Build::SourceTree::Process::DeepCopyReset;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use XML::Simple;
use Encode;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'deepCopyReset'} = \&Process_DeepCopyReset;

sub Process_DeepCopyReset {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "deepCopyReset" && ! $node->{'directive'}->{'processed'} ) {
	    # Trigger deep-copy reset actions on the named objects.
	    my @objects = split(" ",$node->{'directive'}->{'variables'});
	    my $code;
	    foreach my $object ( @objects ) {
		$code .= "call ".$object."\%deepCopyReset()\n";
	    }
	    my $newNode =
	    {
		type       => "code",
		content    => $code,
		firstChild => undef()
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

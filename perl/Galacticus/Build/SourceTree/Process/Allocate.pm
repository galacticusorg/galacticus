# Contains a Perl module which implements allocation of variables with some extra, useful functionality.

package Galacticus::Build::SourceTree::Process::Allocate;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks       {'allocate'} = \&Process_Allocate;
$Galacticus::Build::SourceTree::Hooks::processDependencies{'allocate'} = [ "generics" ];

sub Process_Allocate {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for our directive.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "allocate" && ! $node->{'processed'} ) {
	    # Record that node is processed.
	    $node->{'processed'} = 1;
	    # Get the declaration and determine rank of the shape variable.
	    my $rank;
	    if ( exists($node->{'directive'}->{'rank'}) ) {
		$rank = $node->{'directive'}->{'rank'};
	    } else {
		my $declaration = &Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration($node->{'parent'},$node->{'directive'}->{'variable'});
		foreach my $attribute ( @{$declaration->{'attributes'}} ) {
		    if ( $attribute =~ m/^dimension/ ) {
			$rank = ($attribute =~ tr/,//)+1;
		    }
		}
	    }
	    # Generate code.
	    my $allocator  = "! Auto-generated allocation\n";
	    if (exists($node->{'directive'}->{'shape'})) {
		$allocator    .= "allocate(".$node->{'directive'}->{'variable'}.($rank == 0 ? "" : "(".join(",",map {"".$node->{'directive'}->{'shape'}."(".$_.")"} 1..$rank).")").")\n";
	    } elsif (exists($node->{'directive'}->{'size'}) ) {
		$allocator    .= "allocate(".$node->{'directive'}->{'variable'}.($rank == 0 ? "" : "(".join(",",map {"size(".$node->{'directive'}->{'size' }.",dim=".$_.")"} 1..$rank).")").")\n";
	    } else {
		die("No source given for allocation");
	    }
	    $allocator    .= "! End auto-generated allocation\n";
	    # Create a new node.
	    my $allocatorNode =
	    {
		type       => "code" ,
		content    => $allocator,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$allocatorNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

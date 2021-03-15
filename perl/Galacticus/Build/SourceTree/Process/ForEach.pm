# Contains a Perl module which implements iteration over all elements of arrays of arbitrary dimension.

package Galacticus::Build::SourceTree::Process::ForEach;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks       {'forEach'} = \&Process_ForEach;
$Galacticus::Build::SourceTree::Hooks::processDependencies{'forEach'} = [  "generics" ];

sub Process_ForEach {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for our directive.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "forEach" && ! $node->{'processed'} ) {
	    # Record that node is processed.
	    $node->{'processed'} = 1;
	    # Get the declaration and determine rank.
	    my $declaration = &Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration($node->{'parent'},$node->{'directive'}->{'variable'});
	    my $rank = 0;
	    foreach my $attribute ( @{$declaration->{'attributes'}} ) {
		if ( $attribute =~ m/^dimension/ ) {
		    $rank = ($attribute =~ tr/,//)+1;
		}
	    }
	    # Insert declarations for the indexing variables.
	    if ( $rank > 0 ) {
		my $indices =
		{
		    intrinsic => "integer",
		    type => undef(),
		    openMP => 0,
		    variables => [map {"foreach__".$_} 1..$rank],
		    attributes => []
		};
		&Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},[$indices]);
	    }
	    # Generate code.
	    my $iterator = "! Auto-generated iteration over elements\n";
	    for(my $i=1;$i<=$rank;++$i) {
		$iterator .= "do foreach__".$i."=1,size(".$node->{'directive'}->{'variable'}.",dim=".$i.")\n";
	    }
	    my $indexesFormat = $rank == 0 ? "a1"  : "a1,".join(",",map {"i1"} 1..$rank).",a1";
	    my $indexes       = $rank == 0 ? "'.'" : "'[',".join(",",map {"foreach__".$_} 1..$rank).",']'";
	    my $indexer       = $rank == 0 ? ""    : "(".join(",",map {"foreach__".$_} 1..$rank).")";
	    open(my $code,"<",\$node->{'directive'}->{'content'});
	    while ( my $line = <$code> ) {
		$line      =~ s/%index%/$indexesFormat/g;
		$line      =~ s/\{\{index\}\}/$indexes/g;
		$line      =~ s/\{index\}/$indexer/g;
		$iterator .=  $line;
	    }
	    for(my $i=1;$i<=$rank;++$i) {
		$iterator .= "end do\n";
	    }
	    $iterator   .= "! End auto-generated iteration over elements\n";
	    # Create a new node.
	    my $iteratorNode =
	    {
		type       => "code" ,
		content    => $iterator,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$iteratorNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

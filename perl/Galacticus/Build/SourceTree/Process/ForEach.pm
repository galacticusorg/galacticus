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
	if ( $node->{'type'} eq "forEach" && ! $node->{'directive'}->{'processed'} ) {
	    # Record that node is processed.
	    $node->{'directive'}->{'processed'} = 1;
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
	    my $indexesFormat = $rank == 0 ? "'(a1)'"  : "indexesFormat__";
	    my $indexes       = $rank == 0 ? "'.'" :     join(",",map {"foreach__".$_} 1..$rank)    ;
	    my $indexer       = $rank == 0 ? ""    : "(".join(",",map {"foreach__".$_} 1..$rank).")";
	    my $needFormat    = 0;
	    open(my $code,"<",\$node->{'directive'}->{'content'});
	    while ( my $line = <$code> ) {
		if ( $line =~ m/%index%/ ) {
		    $needFormat = 1;
		}
		$line      =~ s/%index%/$indexesFormat/g;
		$line      =~ s/\{\{index\}\}/$indexes/g;
		$line      =~ s/\{index\}/$indexer/g;
		$iterator .=  $line;
	    }
	    for(my $i=1;$i<=$rank;++$i) {
		$iterator .= "end do\n";
	    }
	    # Add formatter if needed.
	    if ( $needFormat && $rank > 0 ) {
		my @indexesFormatVariables =
		    (
		     {
			 intrinsic => "character",
			 type => "len=5+".$rank."*(7+int(log10(dble(huge(0_c_size_t)))))",
			 openMP => 0,
			 variables => ["indexesFormat__"],
			 attributes => []
		     },
		     {
			 intrinsic => "character",
			 type => "len=12",
			 openMP => 0,
			 variables => ["indexesFormatMeta_"],
			 attributes => []
		     },
		     {
			 intrinsic => "character",
			 type => "len=5+".$rank."*(6+int(log10(dble(huge(0_c_size_t)))))",
			 openMP => 0,
			 variables => ["indexesFormatMeta__"],
			 attributes => []
		     }
		     );
		&Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@indexesFormatVariables);
		my $usesNode =
		{
		    type      => "moduleUse",
		    moduleUse =>
		    {
			ISO_C_Binding =>
			{
			    intrinsic => 1,
			    only      => {c_size_t => 1}
			}
		    },
		    source     => "Galacticus::Build::SourceTree::Process::ForEach::Process_ForEach()",
		    line       => 1
		};
		&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);

		my $format;
		$format .= "write (indexesFormatMeta_,'(\"(\"\"i\"\",i\",i2.2,\".\",i2.2,\")\")') 1+int(log10(dble(huge(0_c_size_t)))),1+int(log10(dble(huge(0_c_size_t))))\n";
		$format .= "indexesFormat__='(\"[\",'\n";
		for(my $i=1;$i<=$rank;++$i) {
		    $format .= "write (indexesFormatMeta__,indexesFormatMeta_) 1+int(log10(dble(size(".$node->{'directive'}->{'variable'}.",dim=".$i."))))\n";
		    $format .= "indexesFormat__=trim(indexesFormat__)//trim(indexesFormatMeta__)".($i < $rank ? "//',\",\",'" : "")."\n";
		}
		$format .= "indexesFormat__=trim(indexesFormat__)//',\"]\")'\n";
		
		$iterator = $format.$iterator;
	    }
	    # Add start and end comments.
	    $iterator   = "! Auto-generated iteration over elements\n"    .
		          $iterator                                       .
	    	          "! End auto-generated iteration over elements\n";
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

# Contains a Perl module which implements some of the missing !GCC$ attribute functionality.

package Galacticus::Build::SourceTree::Process::GCCAttributes;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::SourceTree::Hooks;
use IO::Scalar;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'gccAttributes'} = \&Process_GCCAttributes;

sub Process_GCCAttributes {
    # Get the tree.
    my $tree = shift(); 
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "code" ) {
	    my $newContent;
	    my $code = new IO::Scalar \$node->{'content'};
	    while ( defined(my $line = $code->getline()) ) {
		if ( $line =~ m /^\s*!GCC\$\s+attributes\s+unused\s*::\s*(.*)/ ) {
		    # <workaround type="gfortran" PR="41209" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=41209"/>		    
		    # This is an unused variable attribute. To avoid compiler warnings we test the location (or sizeof for pure
		    # functions where loc() can't be used) of the unused variables and do nothing as a result.
		    # Check that we are in a function or subroutine.
		    (my $variableList = $1) =~ s/\s//g;
		    my $parent = $node->{'parent'};
		    die("Process_GCCAttributes(): 'unused' attribute is relevant only in a function or subroutine")
			unless ( $parent->{'type'} eq "function" || $parent->{'type'} eq "subroutine" );
		    my $functionToUse = $parent->{'opener'} =~ m/^\s*(elemental|pure)\s/ ? "sizeof" : "loc";
		    # Find the final node in the function, or a contains node.
		    my $finalNode = $parent->{'firstChild'};
		    while ( exists($finalNode->{'sibling'}) && ref($finalNode->{'sibling'}) ) {
			$finalNode = $finalNode->{'sibling'};
		    }
		    # Generate a new node.
		    my $newCode = join("",map {"if (".$functionToUse."(".$_.")<0.and.".$functionToUse."(".$_.")>0) then\nend if\n"} split(/,/,$variableList));
		    my $newNode =
		    {
			type       => "code"  ,
			content    => $newCode,
			firstChild => undef()
		    };
		    # Insert at the end of the function, or before any "contains".
		    if ( $finalNode->{'type'} eq "contains" ) {
			&Galacticus::Build::SourceTree::InsertBeforeNode($finalNode,[$newNode]);
		    } else {
			&Galacticus::Build::SourceTree::InsertAfterNode ($finalNode,[$newNode]);
		    }
		    # Neutralize the directive.
		    $line =~ s/^\s*!GCC\$/!GCC/;
		}
		$newContent .= $line;
	    }
	    $code->close();
	    $node->{'content'} = $newContent;
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

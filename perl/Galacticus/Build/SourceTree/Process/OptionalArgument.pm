# Contains a Perl module which implements processing of optional argument directives.

package Galacticus::Build::SourceTree::Process::OptionalArgument;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'optionalArguments'} = \&Process_OptionalArguments;

sub Process_OptionalArguments {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "optionalArgument" && ! $node->{'directive'}->{'processed'} ) {
	    # Record that node is processed.
	    $node->{'directive'}->{'processed'} = 1;
	    # Generate declaration and setting code for the optional argument.
	    my $declaration = &Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration($node->{'parent'},$node->{'directive'}->{'name'});
	    # Change the name of the declared variable, strip argument-based attributes.
	    $declaration->{'variables'} = [ $node->{'directive'}->{'name'}."_" ];
	    @{$declaration->{'attributes'}} =
		map
	    {
		(
		 $_ eq "optional"
		 ||
		 $_ =~ m/^intent\((in)?(out)?\)/
		)
		    ?
		    ()
		    :
		    $_	 
	    }
	    @{$declaration->{'attributes'}};
	    $declaration->{'preprocessor'} = $node->{'directive'}->{'if'}
	        if ( exists($node->{'directive'}->{'if'}) );
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},[$declaration]);
	    # Generate setting code for the optional argument.
	    my $setter       = "   ! Auto-generated optional argument setter\n";
	    $setter         .= "   ".$node->{'directive'}->{'name'}."_=".$node->{'directive'}->{'defaultsTo'}."\n";
	    $setter         .= "   if (present(".$node->{'directive'}->{'name'}.")) ".$node->{'directive'}->{'name'}."_=".$node->{'directive'}->{'name'}."\n";
	    $setter         .= "   ! End auto-generated optional argument setter\n";
	    # Create a new node.
	    my $setterNode =
	    {
		type       => "code" ,
		content    => $setter,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$setterNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

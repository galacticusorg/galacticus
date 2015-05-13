# Contains a Perl module which implements processing of optional argument directives.

package OptionalArguments;
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
require Galacticus::Build::SourceTree::Parse::Declarations;

# Insert hooks for our functions.
$Hooks::processHooks{'optionalArguments'} = \&Process_OptionalArguments;

sub Process_OptionalArguments {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "optionalArgument" && ! $node->{'processed'} ) {
	    # Record that node is processed.
	    $node->{'processed'} = 1;
	    # Generate declaration and setting code for the optional argument.
	    my $declaration = &Declarations::GetDeclaration($node->{'parent'},$node->{'directive'}->{'name'});
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
	    &Declarations::AddDeclarations($node->{'parent'},[$declaration]);
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
	    &SourceTree::InsertAfterNode($node,[$setterNode]);
	}
	$node = &SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

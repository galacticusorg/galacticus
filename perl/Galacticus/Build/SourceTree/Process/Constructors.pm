# Contains a Perl module which implements processing of constructor directives.

package Galacticus::Build::SourceTree::Process::Constructors;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use XML::Simple;
use List::ExtraUtils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'constructors'} = \&Process_Constructors;

sub Process_Constructors {
    # Get the tree.
    my $tree  = shift();
    # Get an XML parser.
    my $xml   = new XML::Simple();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "constructorAssign" && ! $node->{'directive'}->{'processed'} ) {
	    # Assert that our parent is a function.
	    die("Process_Constructors: parent node must be a function")
		unless ( $node->{'parent'}->{'type'} eq "function" );
	    # Determine function return value name.
	    my $returnValueLabel;
	    if ( $node->{'parent'}->{'opener'} =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$/ ) {
		$returnValueLabel = $1;
	    } else {
		$returnValueLabel = $node->{'parent'}->{'name'}; 
	    }
	    # Generate source code for the enumeration.
	    $node->{'directive'}->{'processed'} = 1;
	    my $assignmentSource = "  ! Auto-generated constructor assignment\n";
	    foreach ( grep {$_ ne ""} split(/\s*,\s*/,$node->{'directive'}->{'variables'}) ) {
		my $assigner     = "=";
		my $argumentName = $_;
		if ( $_ =~ m/^\*(.*)/ ) {
		    $assigner     = " => ";
		    $argumentName = $1;
		} 
		$assignmentSource   .= "   ".$returnValueLabel."%".$argumentName.$assigner.$argumentName."\n";
	    }
	    $assignmentSource   .= "  ! End auto-generated constructor assignment.\n\n";
	    # Create a new node.
	    my $newNode =
	    {
		type       => "code"           ,
		content    => $assignmentSource,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

# Contains a Perl module which implements parsing of visibilities in the Galacticus preprocessor system.

package Visibilities;
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
use XML::Simple;
require Fortran::Utils;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Hooks::parseHooks{'visibilities'} = \&Parse_Visibilities;

sub Parse_Visibilities {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Find code blocks.
	if ( $node->{'type'} eq "code" ) {
	    # Initialize a set of new nodes.
	    my @newNodes;
	    # Read the code block, accumulating visibilities as we go.
	    my $rawCode;
	    my $rawVisibility;
	    my $visibilities;
	    open(my $code,"<",\$node->{'content'});
	    do {
		# Get a line.
		&Fortran_Utils::Get_Fortran_Line($code,my $rawLine, my $processedLine, my $bufferedComments);
		# Determine if line is a visibility line.
		my $isVisibility = 0;
		$isVisibility    = 1
		    if ( $processedLine =~ m/^\s*(public|private)\s*(::)?\s*(.*?)\s*$/ );
		# Accumulate raw text.
		if ( $isVisibility == 1 ) {
		    $rawVisibility .= $rawLine;
		    my $visibility = $1;
		    $visibilities->{$visibility} = undef();
		    map {$visibilities->{$visibility}->{$_} = 1} split(/\s*,\s*/,$3)
			if ( $3 );
		} else {
		    $rawCode      .= $rawLine;
		}
		# Process code and visibility blocks as necessary.		
		if ( ( $isVisibility == 1 || eof($code) ) && $rawCode      ) {
		    # Create a new node.
		    my $newNode =
		    {
			type       => "code"  ,
			content    => $rawCode,
			firstChild => undef()
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    push(
			@newNodes,
			$newNode
			);
		    # Reset the raw code text.
		    undef($rawCode);
		}
		if ( ( $isVisibility == 0 || eof($code) ) && $rawVisibility ) {
		    # Create a new node.
		    my $newNode =
		    {
			type       => "visibility"  ,
			visibility => $visibilities
		    };		    
		    $newNode->{'firstChild'} =
		    {
			type       => "code"        ,
			content    => $rawVisibility,
			parent     => $newNode      ,
			sibling    => undef()       ,
			firstChild => undef()
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    # If the end of the code has been reached and we're in a code block, pop that code block from the children
		    # array before we push our visbility node.
		    my $codeNode = pop(@newNodes)
			if ( eof($code) && $isVisibility == 0 );
		    push(
			@newNodes,
			$newNode
			);
		    push(@newNodes,$codeNode)
			if ( eof($code) && $isVisibility == 0 );		    
		    # Reset the raw visibility text.
		    undef($rawVisibility);
		}
	    } until ( eof($code) );
	    close($code);
	    # If we have a single code block, nothing needs to change.
	    unless ( scalar(@newNodes) == 1 && $newNodes[0]->{'type'} eq "code" ) {
		# New nodes created, insert them, replacing the old node.
		&SourceTree::ReplaceNode($node,\@newNodes);
	    }
	}
	$node = &SourceTree::Walk_Tree($node,\$depth);
    }    
}

1;

# Contains a Perl module which implements parsing of module procedures in the Galacticus preprocessor system.

package Galacticus::Build::SourceTree::Parse::ModuleProcedures;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use Fortran::Utils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::parseHooks{'moduleProcedures'} = \&Parse_ModuleProcedures;

sub Parse_ModuleProcedures {
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
	    # Read the code block, accumulating module procedures as we go.
	    my $rawCode;
	    my $rawModuleProcedure;
	    my @moduleProcedures;
	    my $lineNumber             = exists($node->{'line'  }) ? $node->{'line'  } : 0        ;
	    my $source                 = exists($node->{'source'}) ? $node->{'source'} : "unknown";
	    my $rawCodeLine            = $lineNumber;
	    my $rawModuleProcedureLine = $lineNumber;
	    open(my $code,"<",\$node->{'content'});
	    do {
		# Get a line.
		&Fortran::Utils::Get_Fortran_Line($code,my $rawLine, my $processedLine, my $bufferedComments);
		# Determine if line is a module procedure line.
		my $isModuleProcedure = 0;
		$isModuleProcedure    = 1
		    if ( $processedLine =~ m/^\s*module\s+procedure\s*(::)??\s*(\S*?)\s*$/ );
		# Accumulate raw text.
		if ( $isModuleProcedure == 1 ) {
		    $rawModuleProcedure .= $rawLine;
		    push(@moduleProcedures,$2);
		} else {
		    $rawCode            .= $rawLine;
		}
		# Process code and module procedure blocks as necessary.		
		if ( ( $isModuleProcedure == 1 || eof($code) ) && $rawCode            ) {
		    # Create a new node.
		    my $newNode =
		    {
			type       => "code"  ,
			content    => $rawCode,
			firstChild => undef(),
			source     => $source,
			line       => $rawCodeLine
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    push(
			@newNodes,
			$newNode
			);
		    # Reset the raw code text.
		    undef($rawCode);
		    $rawCodeLine            = $lineNumber;
		    $rawModuleProcedureLine = $lineNumber;
		}
		if ( ( $isModuleProcedure == 0 || eof($code) ) && $rawModuleProcedure ) {
		    # Create a new node.
		    my $newNode =
		    {
			type  => "moduleProcedure",
			names => \@moduleProcedures
		    };		    
		    $newNode->{'firstChild'} =
		    {
			type       => "code"        ,
			content    => $rawModuleProcedure,
			parent     => $newNode      ,
			sibling    => undef()       ,
			firstChild => undef(),
			source     => $source,
			line       => $rawModuleProcedureLine
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    # If the end of the code has been reached and we're in a code block, pop that code block from the children
		    # array before we push our module procedure node.
		    my $codeNode = pop(@newNodes)
			if ( eof($code) && $isModuleProcedure == 0 );
		    push(
			@newNodes,
			$newNode
			);
		    push(@newNodes,$codeNode)
			if ( eof($code) && $isModuleProcedure == 0 );		    
		    # Reset the raw module procedure text.
		    undef($rawModuleProcedure);
		    $rawCodeLine            = $lineNumber;
		    $rawModuleProcedureLine = $lineNumber;
		}
		$lineNumber += $rawLine =~ tr/\n//;
	    } until ( eof($code) );
	    close($code);
	    # If we have a single code block, nothing needs to change.
	    unless ( scalar(@newNodes) == 1 && $newNodes[0]->{'type'} eq "code" ) {
		# New nodes created, insert them, replacing the old node.
		&Galacticus::Build::SourceTree::ReplaceNode($node,\@newNodes);
	    }
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }    
}

1;

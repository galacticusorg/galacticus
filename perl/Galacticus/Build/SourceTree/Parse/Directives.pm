# Contains a Perl module which implements parsing of directives in the Galacticus preprocessor system.

package Galacticus::Build::SourceTree::Parse::Directives;
use strict;
use warnings;
use utf8;
use Data::Dumper;
use XML::Simple;
## AJB HACK use Galacticus::Build::SourceTree::Hooks;
## AJB HACK use Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::parseHooks{'directives'} = \&Parse_Directives;

sub Parse_Directives {
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
	    # Read the code block, accumulating directives as we go.
	    my $rawCode;
	    my $rawDirective;
	    my $strippedDirective;
	    my $inDirective = 0;
	    my $directiveRoot;
	    open(my $code,"<",\$node->{'content'});
	    while ( my $line = <$code> ) {
		# Determine if line is a directive line.
		my $isDirective = 0;
		$isDirective    = 1
		    if ( $line =~ m/^\s*\!\#\s+\<([^\s\>]+)/ || $inDirective == 1 );
		$directiveRoot = $1
		    if ( $isDirective == 1 && $inDirective == 0 );		
		# Catch the end of directives.
		my $endDirective = 0;
		$endDirective = 1
		    if ( $isDirective == 1 && $line =~ m/\s*\!\#\s+\<\/$directiveRoot\>/ );
		$endDirective = 1
		    if ( $isDirective == 1 && $inDirective == 0 && $line =~ m/\s*\!\#\s+\<$directiveRoot\s.*\/\>/ );
		# Record whether we are currently in or out of a directive.
		$inDirective = 1
		    if ( $isDirective == 1 );
		# Accumulate raw text.
		if ( $inDirective == 1 ) {
		    (my $strippedLine = $line) =~ s/^\s*\!\#\s*//;
		    $rawDirective      .= $line;
		    $strippedDirective .= $strippedLine;
		} else {
		    $rawCode      .= $line
		};
		# Process code and directive blocks as necessary.
		if ( ( $inDirective == 1 || eof($code) ) && $rawCode      ) {
		    # Create a new node.
		    my $newNode =
		    {
			type       => "code"            ,
			content    => $rawCode          ,
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
		if ( ( $inDirective == 0 || eof($code) || $endDirective ) && $rawDirective ) {
		    # Attempt to parse the directive XML.
		    my $directive = eval{$xml->XMLin($strippedDirective, keepRoot => 1)};
		    if ( $@ ) {
			print "Parse_Directives: failed parsing with message:\n".$@."\n";
			print $strippedDirective;
			die();
		    }
		    my $directiveName = (keys %{$directive})[0];
		    # Create a new node.
		    my $newNode =
		    {
			type       => $directiveName              ,
			directive  => $directive->{$directiveName},
			processed  => 0
		    };
		    $newNode->{'firstChild'} =
		    {
			type       => "code"       ,
			content    => $rawDirective,
			parent     => $newNode     ,
			sibling    => undef()      ,
			firstChild => undef()
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    # If the end of the code has been reached and we're in a code block, pop that code block from the children
		    # array before we push our directive node.
		    my $codeNode = pop(@newNodes)
			if ( eof($code) && $isDirective == 0 );
		    push(
			@newNodes,
			$newNode
			);
		    push(@newNodes,$codeNode)
			if ( eof($code) && $isDirective == 0 );
		    # Reset the raw directive text.
		    $inDirective = 0;
		    undef($rawDirective     );
		    undef($strippedDirective);
		}
	    }
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

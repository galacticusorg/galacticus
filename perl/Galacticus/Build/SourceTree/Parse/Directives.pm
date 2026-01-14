# Contains a Perl module which implements parsing of directives in the Galacticus preprocessor system.

package Galacticus::Build::SourceTree::Parse::Directives;
use strict;
use warnings;
use utf8;
use Data::Dumper;
use XML::Simple;
use XML::SAX::ParserFactory;
use XML::Validator::Schema;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::parseHooks      {'directives'} = \&Parse_Directives;
$Galacticus::Build::SourceTree::Hooks::postprocessHooks{'directives'} = \&PostProcess_Directives;

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
	    my $rawOpener;
	    my $rawDirective;
	    my $strippedDirective;
	    my $inDirective      = 0;
	    my $inXML            = 0;
	    my $directiveRoot;
	    my $lineNumber       = exists($node->{'line'  }) ? $node->{'line'  } : 0        ;
	    my $source           = exists($node->{'source'}) ? $node->{'source'} : "unknown";
	    my $rawCodeLine      = $lineNumber;
	    my $rawDirectiveLine = $lineNumber;
	    open(my $code,"<",\$node->{'content'});
	    while ( my $line = <$code> ) {
		# Detect the end of an XML section and change state.
		$inXML = 0
		    if ( $line =~ m/^\s*!!\]/ );
		# Process XML blocks.
		my $isDirective  = 0;
		my $endDirective = 0;
		# Strip the initial "!<" and any following whitespace (but not tabs - this allows us to use tabs for formatting
		# purposes).
		(my $strippedLine = $line) =~ s/^\s*\!<\s*//;
		if ( $inXML ) {
		    # Determine if line is a directive line.
		    $isDirective    = 1
			if ( $strippedLine =~ m/^\s*\<([^\s\>\/]+)/ || $inDirective == 1 );
		    $directiveRoot = $1
			if ( $isDirective == 1 && $inDirective == 0 );		
		    # Catch the end of directives.
		    $endDirective = 1
			if ( $isDirective == 1 && $strippedLine =~ m/\s*\<\/$directiveRoot\>/ );
		    $endDirective = 1
			if ( $isDirective == 1 && $inDirective == 0 && ( $strippedLine =~ m/\s*\<$directiveRoot\s.*\/\>/ || $strippedLine =~ m/\s*\<$directiveRoot\/\>/ ) );
		    # Record whether we are currently in or out of a directive.
		    $inDirective = 1
			if ( $isDirective == 1 );
		}
		# Accumulate raw text.
		if ( $inDirective ) {
		    # Process non-breaking spaces as a special case.
		    $strippedLine =~ s/&nbsp;/ /g;
		    # Accumulate the line.
		    $rawDirective      .= $line;
		    $strippedDirective .= $strippedLine;
		} elsif ( $line !~ m/^\s*!!(\[|\])/ ) {
		    $rawCode           .= $line;
		} elsif ( $line =~ m/^\s*!!\[/ ) {
		    $rawOpener          = $line;
		} else {
		    $rawCodeLine      = $lineNumber+1;
		    $rawDirectiveLine = $lineNumber+1;
		}
		# Process code and directive blocks as necessary.
		if ( ( $inDirective == 1 || eof($code) ) && $rawCode      ) {
		    # Create a new node.
		    my $newNode =
		    {
			type       => "code"            ,
			content    => $rawCode          ,
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
		    $rawCodeLine      = $lineNumber;
		    $rawDirectiveLine = $lineNumber;
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
		    # Validate the directive if possible.
		    if ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/schema/".$directiveName.".xsd" ) {
			my $validator = XML::Validator::Schema->new(file => $ENV{'GALACTICUS_EXEC_PATH'}."/schema/".$directiveName.".xsd");
			my $parser    = XML::SAX::ParserFactory->parser(Handler => $validator); 
			eval { $parser->parse_string($strippedDirective) };
			if ( $@ ) {
			    my $nodeParent = $node;
			    while ( $nodeParent->{'type'} ne "file" ){
				$nodeParent = $nodeParent->{'parent'};
			    }
			    die "Galacticus::Build::SourceTree::Parse::Directives(): validation failed in file ".$nodeParent->{'name'}." at line ".$node->{'line'}.":\n".$@;
			}
		    }
		    # Create a new node.
		    my $newNode =
		    {
			type       => $directiveName              ,
			directive  => $directive->{$directiveName},
			line       => $node     ->{'line'        },
			processed  => 0                           ,
			source     => $source                     ,
			line       => $rawDirectiveLine
		    };
		    (my $rawCloser = $rawOpener) =~ s/\[/\]/;
		    $rawDirective = $rawOpener.$rawDirective.$rawCloser;
		    $newNode->{'firstChild'} =
		    {
			type       => "code"           ,
			content    => $rawDirective    ,
			parent     => $newNode         ,
			sibling    => undef()          ,
			firstChild => undef()          ,
			source     => $source          ,
			line       => $rawDirectiveLine
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
		    $rawCodeLine      = $lineNumber;
		    $rawDirectiveLine = $lineNumber;
		}
		# Detect the start of an XML section and change state.
		$inXML = 1
		    if ( $line =~ m/^\s*!!\[/ );
		# Increment line number count.
		++$lineNumber;
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

sub PostProcess_Directives {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Find directives.
	if ( exists($node->{'directive'}) ) {
	    die("directive '".$node->{'type'}."' was not processed at line ".$node->{'line'}." in ".$tree->{'name'})
		unless ( exists($node->{'directive'}->{'processed'}) && $node->{'directive'}->{'processed'} );
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }    
}

1;

# Contains a Perl module which implements parsing of module uses in the Galacticus preprocessor system.

package Galacticus::Build::SourceTree::Parse::ModuleUses;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Fortran::Utils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::parseHooks{'moduleUses'} = \&Parse_ModuleUses;

sub Parse_ModuleUses {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Find code blocks.
	if ( $node->{'type'} eq "code" ) {
	    # Initialize a set of new nodes.
	    my @newNodes;
	    # Read the code block, accumulating "uses" as we go.
	    my $rawCode;
	    my $rawModuleUse;
	    my $moduleUses;
	    my @preprocessorStack;
	    my $lineNumber       = exists($node->{'line'  }) ? $node->{'line'  } : 0        ;
	    my $source           = exists($node->{'source'}) ? $node->{'source'} : "unknown";
	    my $rawCodeLine      = $lineNumber;
	    my $rawModuleUseLine = $lineNumber;
	    open(my $code,"<",\$node->{'content'});
	    do {
		# Get a line.
		&Fortran::Utils::Get_Fortran_Line($code,my $rawLine, my $processedLine, my $bufferedComments);
		# Determine if line is a module use line.
		my $isModuleUse = 0;
		$isModuleUse    = 1
		    if ( $processedLine =~ m/^\s*(!\$)?\s*use\s*(\s+|,\s*(intrinsic))\s*(::)?\s*([a-zA-Z0-9_]+)\s*(,\s*only\s*:)?\s*([a-zA-Z0-9_\(\)=,\s]+)?\s*$/ );
		# Accumulate raw text.
		if ( $isModuleUse == 1 ) {
		    $rawModuleUse .= $rawLine;
		    my $isOpenMP    = $1;
		    my $isIntrinsic = $3;
		    my $moduleName  = $5;
		    my $only        = $7;
		    $moduleUses->{$moduleName}->{'openMP'} = $isOpenMP ? 1 : 0;
		    if ( $isIntrinsic ) {
			$moduleUses->{$moduleName}->{'intrinsic'} = 1;
		    } else {
			$moduleUses->{$moduleName}->{'intrinsic'} = 0;
		    }
		    if ( $only && ! exists($moduleUses->{$moduleName}->{'all'}) ) {
			chomp($only);
			map {$moduleUses->{$moduleName}->{'only'}->{$_} = 1} split(/\s*,\s*/,$only);
		    } else {
			$moduleUses->{$moduleName}->{'all'} = 1;
		    }
		    @{$moduleUses->{$moduleName}->{'conditions'}} = @preprocessorStack
			if ( scalar(@preprocessorStack) > 0 );
		} elsif ( $processedLine =~ m/^#/ && $rawModuleUse ) {
		    $isModuleUse = 1;
		    $rawModuleUse .= $rawLine
		} else {
		    $rawCode      .= $rawLine;
		}
		# Handle preprocessor lines. Note that for #else and #endif directives we check that we have elements in the
		# stack. Since we process this on a node-by-node basis it's possible that the matching #if directive was in a
		# previous node, so isn't caught. Ideally we need to do a better job of parsing the preprocessor directives by
		# walking through sibling nodes when necessary.
		push(@preprocessorStack,{name => $1, invert => 0})
		    if ( $processedLine =~ m/^#ifdef\s+([A-Z0-9]+)/ );
		push(@preprocessorStack,{name => $1, invert => 1})
		    if ( $processedLine =~ m/^#ifndef\s+([A-Z0-9]+)/ );
		pop(@preprocessorStack)
		    if ( $processedLine =~ m/^#endif/ && scalar(@preprocessorStack) > 0 );
		if ( $processedLine =~ m/^#else/ && scalar(@preprocessorStack) > 0 ) {
		    my $descriptor = pop(@preprocessorStack);
		    die("Galacticus::Build::SourceTree::Parse::ModuleUses::Parse_ModuleUses(): invert is not defined in file '".$node->{'source'}."' in node at line ".$node->{'line'})
			unless ( exists($descriptor->{'invert'}) && defined($descriptor->{'invert'}) );
		    $descriptor->{'invert'} = 1-$descriptor->{'invert'};
		    push(@preprocessorStack,$descriptor);
		}
		# Process code and module use blocks as necessary.		
		if ( ( $isModuleUse == 1 || eof($code) ) && $rawCode      ) {
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
		    $rawCodeLine      = $lineNumber;
		    $rawModuleUseLine = $lineNumber;
		}
		if ( ( $isModuleUse == 0 || eof($code) ) && $rawModuleUse ) {
		    # Create a new node.
		    my $newNode =
		    {
			type      => "moduleUse"  ,
			moduleUse => $moduleUses
		    };
		    $newNode->{'firstChild'} =
		    {
			type       => "code"        ,
			content    => $rawModuleUse,
			parent     => $newNode      ,
			sibling    => undef()       ,
			firstChild => undef()       ,
			source       => $source            ,
			line         => $rawModuleUseLine
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    # If the end of the code has been reached and we're in a code block, pop that code block from the children
		    # array before we push our visbility node.
		    my $codeNode = pop(@newNodes)
			if ( eof($code) && $isModuleUse == 0 );
		    push(
			@newNodes,
			$newNode
			);
		    push(@newNodes,$codeNode)
			if ( eof($code) && $isModuleUse == 0 );		    
		    # Reset the raw module use text.
		    undef($rawModuleUse);
		    $rawCodeLine      = $lineNumber;
		    $rawModuleUseLine = $lineNumber;
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

sub AddUses {
    # Grab the node to add uses to, and the new uses to add.
    my $node       = shift();
    my $moduleUses = shift();
    # Locate a moduleUses node.
    my $usesNode;
    my $childNode = $node->{'firstChild'};
    while ( $childNode ) {
	$usesNode = $childNode
	    if ( $childNode->{'type'} eq "moduleUse" );
	$childNode = $childNode->{'sibling'};
    }
    unless ( $usesNode ) {
	# No module use node exists - insert one.
	$usesNode =
	{
	    type       => "moduleUse",
	    sibling    => undef()    ,
	    parent     => undef()    ,
	    moduleUse  => undef()
	};
	$usesNode->{'firstChild'} =
	{
	    type       => "code"   ,
	    content    => ""       ,
	    sibling    => undef()  ,
	    parent     => $usesNode,
	    firstChild => undef()
	};
	&Galacticus::Build::SourceTree::InsertBeforeNode($node->{'firstChild'},[$usesNode]);
    }
    # Merge the module usages.
    foreach my $moduleName ( keys(%{$moduleUses->{'moduleUse'}}) ) {
	$usesNode->{'moduleUse'}->{$moduleName}->{'openMP'   } = $moduleUses->{'moduleUse'}->{$moduleName}->{'openMP'   };
	$usesNode->{'moduleUse'}->{$moduleName}->{'intrinsic'} = $moduleUses->{'moduleUse'}->{$moduleName}->{'intrinsic'};
	unless ( exists($usesNode->{'moduleUse'}->{$moduleName}->{'all'}) ) {
	    if ( exists($moduleUses->{'moduleUse'}->{$moduleName}->{'all'}) ) {
		$usesNode->{'moduleUse'}->{$moduleName}->{'all'} = 1;
		delete($usesNode->{'moduleUse'}->{$moduleName}->{'only'});
	    } else {
		$usesNode->{'moduleUse'}->{$moduleName}->{'only'}->{$_} = 1x2
		    foreach ( keys(%{$moduleUses->{'moduleUse'}->{$moduleName}->{'only'}}) );
	    }
	}
    }
    # Update the contained code.
    $usesNode->{'firstChild'}->{'content'} = undef();
    foreach my $moduleName ( keys(%{$usesNode->{'moduleUse'}}) ) {
	if ( exists($usesNode->{'moduleUse'}->{$moduleName}->{'conditions'}) ) {
	    foreach ( @{$usesNode->{'moduleUse'}->{$moduleName}->{'conditions'}} ) {
		$usesNode->{'firstChild'}->{'content'} .= ($_->{'invert'} ? "#ifndef" : "#ifdef")." ".$_->{'name'}."\n";
	    }
	}
	$usesNode->{'firstChild'}->{'content'} .= "   ";
	$usesNode->{'firstChild'}->{'content'} .= "!\$ "
	    if ( $usesNode->{'moduleUse'}->{$moduleName}->{'openMP'} );
	$usesNode->{'firstChild'}->{'content'} .= "use";
	$usesNode->{'firstChild'}->{'content'} .= ", intrinsic"
	    if ( $usesNode->{'moduleUse'}->{$moduleName}->{'intrinsic'} );
	$usesNode->{'firstChild'}->{'content'} .= " :: ".$moduleName;
	$usesNode->{'firstChild'}->{'content'} .= ", only : ".join(", ",keys(%{$usesNode->{'moduleUse'}->{$moduleName}->{'only'}}))
	    if ( $usesNode->{'moduleUse'}->{$moduleName}->{'only'} );
  	$usesNode->{'firstChild'}->{'content'} .= "\n";
	if ( exists($usesNode->{'moduleUse'}->{$moduleName}->{'conditions'}) ) {
	    foreach ( @{$usesNode->{'moduleUse'}->{$moduleName}->{'conditions'}} ) {
		$usesNode->{'firstChild'}->{'content'} .= "#endif\n";
	    }
	}
    }
}

1;

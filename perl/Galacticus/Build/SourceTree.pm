# Contains a Perl module which implements a tree-based preprocessor for Galacticus Fortran source files.

package SourceTree;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");
use File::Slurp qw(slurp);
use Data::Dumper;
use Scalar::Util qw(reftype);
require Fortran::Utils;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree::Parse::Directives;
require Galacticus::Build::SourceTree::Parse::Visibilities;
require Galacticus::Build::SourceTree::Parse::ModuleUses;
require Galacticus::Build::SourceTree::Parse::Declarations;
require Galacticus::Build::SourceTree::Process::Enumeration;
require Galacticus::Build::SourceTree::Process::InputParameter;
require Galacticus::Build::SourceTree::Process::InputParameterList;
require Galacticus::Build::SourceTree::Process::FunctionClass;
require Galacticus::Build::SourceTree::Process::OptionalArgument;
require Galacticus::Build::SourceTree::Process::Generics;
require Galacticus::Build::SourceTree::Process::SourceDigest;                                                                                
require Galacticus::Build::SourceTree::Process::ObjectBuilder;                                                                               
require Galacticus::Build::SourceTree::Process::DebugHDF5;

sub ParseFile {
    # Grab the file name.
    my $fileName = shift();
    # Read the code.
    my $code = slurp($fileName);
    # Initialize root object.
    (my $fileRootName = $fileName) =~ s/^.*\/([^\/]+)$/$1/;
    my $tree =
    {
	type       => "file"       ,
	name       => $fileRootName,
	content    => $code        ,
	parent     => undef()      ,
	firstChild => undef()      ,
	sibling    => undef()
    };
    &BuildTree($tree);
    return $tree;
}

sub BuildTree {
    # Get the root of the tree;
    my $tree = shift();
    # Build the tree.
    my @stack = ( $tree );
    my $i = 0;
    while ( scalar(@stack) > 0 ) {
	++$i;
	my $current = pop(@stack);
	&Parse_Unit($current);
	my $child = $current->{'firstChild'};
	while ( $child ) {
	    push(@stack,$child);
	    $child = $child->{'sibling'};
	}
    }
    # Run all defined parsers on the tree.
    my @unitParsers = ( "directives" );
    for(my $stage=0;$stage<2;++$stage) {
	foreach my $parser ( keys(%Hooks::parseHooks) ) {
	    &{$Hooks::parseHooks{$parser}}($tree)
		if ( 
		    ( $stage == 0 &&   grep {$_ eq $parser} @unitParsers )
		    ||
		    ( $stage == 1 && ! grep {$_ eq $parser} @unitParsers )
		);
	}
    }
}

sub ProcessTree {
    # Get the tree.
    my $tree      = shift();
    my (%options) = @_;
    # Run all defined processors on the tree.
    &{$Hooks::processHooks{$_}}($tree,\%options)
	foreach ( keys(%Hooks::processHooks) );
}

sub Parse_Unit {
    my $unit = shift();
    # Process any content unless this is a code block.
    if ( exists($unit->{'content'}) && $unit->{'type'} ne "code" ) {
	my @children = &Build_Children($unit->{'content'});
	$unit->{'firstChild'} = $children[0];
	for(my $i=0;$i<scalar(@children);++$i) {
	    $children[$i]->{'parent'  } = $unit;
	    $children[$i]->{'sibling' } = $children[$i+1]
		unless ( $i == scalar(@children)-1 );
	}
	delete($unit->{'content'});
    }
    return;
}

sub Build_Children {
    # Grab the code passed to us.
    my $codeText = shift();
    # Initialize an empty array of children.
    my @children = ();
    # Initialize current raw code buffer.
    my $rawCode;
    # Initialize custom openers;
    my %unitOpeners = %Fortran_Utils::unitOpeners;
    $unitOpeners{'contains'} = { unitName => -1, regEx => qr/^\s*contains\s*$/ };
    # Connect a file handle to the code text.
    open(my $code,"<",\$codeText);
    # Read lines.
    do {
	# Get a line.
	my ($rawLine, $processedLine, $bufferedComments);
	&Fortran_Utils::Get_Fortran_Line($code,$rawLine,$processedLine,$bufferedComments);
	# Check for unit opening.
	my $unitFound;	
	foreach my $unitType ( keys(%unitOpeners) ) {
	    if ( my @matches = ( $processedLine =~ $unitOpeners{$unitType}->{'regEx'} ) ) {
		# A match is found, extract the name of the unit.
	     	my $unitName    = $matches[$unitOpeners{$unitType}->{'unitName'}]
		    if ( $unitOpeners{$unitType}->{'unitName'} >= 0 );
		$unitFound      = 1;
		# Store the unit opener.
		my $opener      = $rawLine;
		my $closer                ;
		my $unitRawCode           ;
		# Push any existing raw code into a "code" child.
		push(
		    @children,
		    {
			type       => "code"  ,
			content    => $rawCode,
			parent     => undef() ,
			firstChild => undef() ,
			sibling    => undef()
		    }
		    )
		    if ( $rawCode );
		undef($rawCode);
		# Read ahead until the matching unit closer is found, storing the unit code.
		my $closerRegExp;
		if ( $unitType eq "contains" ) {
		    $closerRegExp = qr/(?!)/; # Never match.
		} else {
		    $closerRegExp = $Fortran_Utils::unitClosers{$unitType}->{'regEx'};
		}
		do {
		    # Get the next line.
		    if ( eof($code) ) {
			$rawLine       = "";
			$processedLine = "";
		    } else {
			&Fortran_Utils::Get_Fortran_Line($code,$rawLine, $processedLine, $bufferedComments);
		    }
		    if ( ( $unitType eq "contains" && eof($code) ) || ( my @matches = ( $processedLine =~ $closerRegExp ) ) ) {
			my $closeUnitName = $matches[$Fortran_Utils::unitClosers{$unitType}->{'unitName'}]
			    unless ( $unitType eq "contains" );
			if ( $unitType eq "contains" || ( defined($unitName) && $closeUnitName eq $unitName ) ) {
			    # Matching unit closer has been found.
			    if ( $unitType eq "contains" ) {
				$unitRawCode .= $rawLine;
			    } else {
				$closer = $rawLine;
			    }
			    # Create a new node.
			    my $newNode = 
			    {
				type       => $unitType   ,
				opener     => $opener     ,
				parent     => undef()     ,
				firstChild => undef()     ,
				sibling    => undef()
			    };
			    $newNode->{'content'} = $unitRawCode
				if ( $unitRawCode );
			    $newNode->{'name'  } = $unitName
				if ( $unitName );
			    $newNode->{'closer'} = $closer
				if ( $closer   );	
			    # Append a child.
			    push(
				@children,
				$newNode
				);
			    last;			    
			}
		    }
		    $unitRawCode .= $rawLine
			unless ( $closer );
		} until ( eof($code) );
		die("Build_Children: no matching unit closer was found for opener '".$opener."'")
		    unless ( $closer );
 	    }
	    last
		if ( $unitFound );
	}
	# If no unit was found, treat this as raw code.
	unless ( $unitFound ) {
	    # Append this line to the current raw code buffer.
	    $rawCode .= $rawLine;
	}	
    } until ( eof($code) );
    close($code);
    # If code exists in the raw code buffer, append it as a child node.
    push(
	@children,
	{
     	    type       => "code"  ,
     	    content    => $rawCode,
     	    parent     => undef() ,
     	    firstChild => undef() ,
	    sibling    => undef()
    	}
     	)
	if ( $rawCode );
    # Return the child nodes.
    return @children;
}

sub Show_Tree {
    # Show the structure of the source code tree.
    my $tree = shift();
    # Build a breadth first stack of the tree.
    my @stack = ();
    @stack = &StackIt($tree,\@stack,-1);
    my $depthPrevious = 0;
    foreach ( @stack ) {
	my $indent = "  " x $_->{'depth'};
	print $indent."|\n"
	    if ( $_->{'depth'} > $depthPrevious );
	print $indent."--> ".$_->{'node'}->{'type'};
	print " {".$_->{'node'}->{'name'}."}"
	    if ( exists($_->{'node'}->{'name'}) );
	print"\n";
	$depthPrevious = $_->{'depth'};
    }
}

sub StackIt {
    my $node  =   shift() ;
    my @stack = @{shift()};
    my $depth =   shift() ;
    ++$depth;
    push(@stack,{node => $node, depth => $depth});
    if ( $node->{'firstChild'} ) {
	my $child = $node->{'firstChild'};
	while ( $child ) {			      
	    @stack = &StackIt($child,\@stack,$depth);
	    $child = $child->{'sibling'};
	}
    }
    return @stack;
}

sub Walk_Tree {
    # Do a depth-first walk of the tree.
    my $node     = shift();
    my $depthRef = shift();
    # Walk the tree.
    unless ( $node->{'parent'} ) {
    	while ( $node->{'firstChild'} ) {
	    ++${$depthRef};
	    $node = $node->{'firstChild'};
	}
    } else {
	if ( $node->{'sibling'} ) {
	    $node = $node->{'sibling'};
	    while ( $node->{'firstChild'} ) {		
		++${$depthRef};
		$node = $node->{'firstChild'};
	    }
	} else {
	    --${$depthRef};
	    $node = $node->{'parent'};
	    undef($node)
		unless ( $node->{'parent'} );
	}
    }
    return $node;
}

sub ReplaceNode {
    # Grab the node to replace and the replacements.
    my $node     =   shift() ;
    my @newNodes = @{shift()};
    # Find the prior sibling to the node.
    if ( $node->{'parent'} ) {
	my $child = $node->{'parent'}->{'firstChild'};
	if ( $child == $node ) {
	    $newNodes[$#newNodes]->{'sibling'} = $node->{'parent'}->{'firstChild'}->{'sibling'};
	    $node->{'parent'}->{'firstChild'} = $newNodes[0];
	} else {
	    while ( $child ) {
		my $nextChild = $child->{'sibling'};
		if ( defined($nextChild) && $nextChild == $node ) {
		    $newNodes[$#newNodes]->{'sibling'} = $nextChild->{'sibling'};
		    $child               ->{'sibling'} = $newNodes[0];
		    last;
		}
		$child = $nextChild;
	    }
	}
	for(my $i=0;$i<scalar(@newNodes);++$i) {
	    $newNodes[$i]->{'parent'  } = $node->{'parent'};
	    $newNodes[$i]->{'sibling' } = $newNodes[$i+1]
		unless ( $i == scalar(@newNodes)-1 );
	}
	return $newNodes[0];
    } else {
	die("ReplaceNode: cannot replace root node");
    }
}

sub Serialize {
    my $node = shift();
    my $serialization;
    my $currentNode = $node;
    while ( $currentNode ) {
	if ( $currentNode->{'type'} eq "code" ) {
	    $serialization .= $currentNode->{'content'}
	} else {
	    $serialization .= $currentNode->{'opener'}
	        if ( exists($currentNode->{'opener'}) );
	    $serialization .= &Serialize($currentNode->{'firstChild'})
		if ( $currentNode->{'firstChild'} );
	    $serialization .= $currentNode->{'closer'}
	        if ( exists($currentNode->{'closer'}) );
	}
	if ( $currentNode->{'sibling'   } ) {
	    $currentNode = $currentNode->{'sibling'};
	} else {
	    undef($currentNode);
	}
    }
    return $serialization;
}

sub InsertAfterNode {
    # Grab the node to insert after and the insertion.
    my $node     =   shift() ;
    my @newNodes = @{shift()};
    # Check that the node has a parent.
    if ( $node->{'parent'} ) {
	$newNodes[$#newNodes]->{'sibling'} = $node->{'sibling'};
	for(my $i=0;$i<scalar(@newNodes);++$i) {
	    $newNodes[$i]->{'parent' } = $node->{'parent'};	    
	    $newNodes[$i]->{'sibling'} = $newNodes[$i+1]
		unless ( $i == scalar(@newNodes)-1 );
	}
	$node   ->{'sibling'} = $newNodes[0];
    } else {
	die("InsertAfterNode: cannot insert after root node");
    }
}

sub InsertBeforeNode {
    # Grab the node to insert before and the insertion.
    my $node     =   shift() ;
    my @newNodes = @{shift()};
    # Check that the node has a parent.
    if ( $node->{'parent'} ) {
	# Find where to insert.
	my $child = $node->{'parent'}->{'firstChild'};
	if ( $child == $node ) {
	    $newNodes[$#newNodes]->{'sibling'} = $child;
	    $node->{'parent'}->{'firstChild'} = $newNodes[0];
	    for(my $i=0;$i<scalar(@newNodes);++$i) {
	    	$newNodes[$i]->{'parent' } = $node->{'parent'};	    
	    	$newNodes[$i]->{'sibling'} = $newNodes[$i+1]
	    	    unless ( $i == scalar(@newNodes)-1 );
	    }
	} else {
	    while ( $child->{'sibling'} != $node ) {
		$child = $child->{'sibling'};
	    }
	    &InsertAfterNode($child,\@newNodes);
	}
    } else {
	die("InsertBeforeNode: cannot insert before root node");
    }
}

sub PrependChildToNode {
    # Insert nodes as children of the given node, at the start of the list of children.
    my $node     =   shift() ;
    my @newNodes = @{shift()};
    # Check if the node already has children.
    if ( $node->{'firstChild'} ) {
	# It does, insert before the first child.
	&InsertBeforeNode($node->{'firstChild'},\@newNodes);
    } else {
	# It does not, create children.
	$node->{'firstChild'} = $newNodes[0];
	for(my $i=0;$i<scalar(@newNodes);++$i) {
	    $newNodes[$i]->{'parent' } = $node;
	    $newNodes[$i]->{'sibling'} = $newNodes[$i+1]
		unless ( $i == scalar(@newNodes)-1 );
	}
   }
}

sub InsertPreContains {
    # Grab the node and nodes to insert.
    my $node     = shift();
    my $newNodes = shift();
    # Find the container if one exists.
    my $child    = $node->{'firstChild'};
    my $containsNode;
    my $lastChild;
    while ( $child ) {
	last
	    if ( $child->{'type'} eq "contains" );
	$lastChild = $child;
	$child     = $child->{'sibling'};
    }
    &InsertAfterNode($lastChild,$newNodes);
}

sub InsertPostContains {
    # Grab the node and nodes to insert.
    my $node     = shift();
    my $newNodes = shift();
    # Find the container if one exists.
    my $child    = $node->{'firstChild'};
    my $containsNode;
    my $lastChild;
    while ( $child ) {
	$containsNode = $child	
	    if ( $child->{'type'} eq "contains" );
	$lastChild = $child;
	$child     = $child->{'sibling'};
    }
    # If no container exists, add one now.
    unless ( $containsNode ) {
	$containsNode =
	{
	    type       => "contains"  ,
	    opener     => "contains\n",
	    closer     => "\n"        ,
	    firstChild => undef()
	};
	&InsertAfterNode($lastChild,[$containsNode]);
    }
    # Insert the new nodes.
    &PrependChildToNode($containsNode,$newNodes);
}

sub SetVisibility {
    # Grab the node, name, and visibility.
    my $node       = shift();
    my $unitName   = shift();
    my $visibility = shift();
    # Verify visibility.
    die("SetVisibility: visibility must be 'public' or 'private'")
	unless ( $visibility eq "public" || $visibility eq "private" );
    # Find the visibility node if one exists.
    my $child      = $node->{'firstChild'};
    my $visibilityNode;
    my $lastChild;
    while ( $child ) {
	$visibilityNode = $child	
	    if ( $child->{'type'} eq "visibility" );
	$lastChild = $child;
	$child     = $child->{'sibling'};
    }
    # If no visibility exists, add one now.
    unless ( $visibilityNode ) {
	$visibilityNode =
	{
	    type       => "visibility",
	    firstChild =>
	    {
		type       => "code"         ,
		firstChild => undef()        ,
		sibling    => undef()
	    }
	};
	$visibilityNode->{'firstChild'}->{'parent'} = $visibilityNode;

	# If the node has a module use block, visibilities must appear after it.
	my $child         = $node->{'firstChild'};
	my $moduleUseNode;
	while ( $child ) {
	    $moduleUseNode = $child
		if ( $child->{'type'} eq "moduleUse" );
	    $child = $child->{'sibling'};
	}
	if ( $moduleUseNode ) {
	    # A module use block exists, place visibilities after it.
	    &InsertAfterNode   ($moduleUseNode,[$visibilityNode]);
	} else {	
	    # No module use block exists, place visibilities as first child of the node.
	    &PrependChildToNode($node         ,[$visibilityNode]);
	}
    }
    # Set visibility.
    $visibilityNode->{'visibility'}->{$visibility}->{$unitName} = 1;
    # Rebuild content.
    my $newContent;
    foreach ( 'private', 'public' ) {
	if ( exists($visibilityNode->{'visibility'}->{$_}) ) {
	    $newContent .= "  ".$_;
	    $newContent .= " :: ".join(", ",keys(%{$visibilityNode->{'visibility'}->{$_}}))
		if ( $visibilityNode->{'visibility'}->{$_} );
	    $newContent .= "\n";
	}
    }
    $visibilityNode->{'firstChild'}->{'content'} = $newContent."\n";
}

1;

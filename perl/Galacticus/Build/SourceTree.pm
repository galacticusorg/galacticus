# Contains a Perl module which implements a tree-based preprocessor for Galacticus Fortran source files.

package Galacticus::Build::SourceTree;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Scalar::Util qw(reftype);
use Fortran::Utils;
use File::Slurp;
use Sort::Topo;
use Galacticus::Build::SourceTree::Hooks;
use Galacticus::Build::SourceTree::Parse::Directives;
use Galacticus::Build::SourceTree::Parse::Visibilities;
use Galacticus::Build::SourceTree::Parse::ModuleUses;
use Galacticus::Build::SourceTree::Parse::Declarations;
use Galacticus::Build::SourceTree::Parse::ModuleProcedures;
use Galacticus::Build::SourceTree::Process::Enumeration;
use Galacticus::Build::SourceTree::Process::InputParameter;
use Galacticus::Build::SourceTree::Process::InputParametersValidate;
use Galacticus::Build::SourceTree::Process::FunctionClass;
use Galacticus::Build::SourceTree::Process::StateStorable;
use Galacticus::Build::SourceTree::Process::DeepCopyActions;
use Galacticus::Build::SourceTree::Process::OptionalArgument;
use Galacticus::Build::SourceTree::Process::ForEach;
use Galacticus::Build::SourceTree::Process::Allocate;
use Galacticus::Build::SourceTree::Process::Generics;
use Galacticus::Build::SourceTree::Process::SourceDigest;
use Galacticus::Build::SourceTree::Process::SourceIntrospection;
use Galacticus::Build::SourceTree::Process::ObjectBuilder;
use Galacticus::Build::SourceTree::Process::DeepCopyReset;
use Galacticus::Build::SourceTree::Process::DeepCopyFinalize;
use Galacticus::Build::SourceTree::Process::DebugHDF5;
use Galacticus::Build::SourceTree::Process::DebugMPI;
use Galacticus::Build::SourceTree::Process::ProfileOpenMP;
use Galacticus::Build::SourceTree::Process::GSLConstants;
use Galacticus::Build::SourceTree::Process::HDF5FCInterop;
use Galacticus::Build::SourceTree::Process::Constructors;
use Galacticus::Build::SourceTree::Process::ConditionalCall;
use Galacticus::Build::SourceTree::Process::EventHooks;
use Galacticus::Build::SourceTree::Process::ClassDocumentation;
use Galacticus::Build::SourceTree::Analyze::UseDuplication;
use Encode;

sub ParseFile {    
    # Grab the file name.
    my $fileName = shift();
    my (%options) = @_;
    # Read the code.
    my $code = &Galacticus::Build::SourceTree::Process::SourceIntrospection::ReadFile($fileName,%options);
    # Initialize root object.
    (my $fileRootName = $fileName) =~ s/^.*\/([^\/]+)$/$1/;
    my $tree =
    {
	type       => "file"       ,
	name       => $fileRootName,
	content    => $code        ,
	parent     => undef()      ,
	firstChild => undef()      ,
	sibling    => undef()      ,
	source     => $fileName    ,
	line       => 0
    };
    &BuildTree($tree);
    return $tree;
}

sub ParseCode {
    # Grab the source to parse.
    my $code     = shift();
    my $fileName = shift();
    my (%options) = @_;
    # Instrument code.
    $code = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Instrument($code)
	unless ( exists($options{'instrument'}) && ! $options{'instrument'} );
    # Initialize root object.
    (my $fileRootName = $fileName) =~ s/^.*\/([^\/]+)$/$1/;
    my $tree =
    {
	type       => "file"       ,
	name       => $fileRootName,
	content    => $code        ,
	parent     => undef()      ,
	firstChild => undef()      ,
	sibling    => undef()      ,
	source     => $fileName    ,
	line       => 0
    };
    &BuildTree($tree);
    return $tree;
}

sub Children {
    # Return an array of children.
    my $node = shift();
    my @children;
    my $child = $node->{'firstChild'};
    while ( $child ) {
	push(@children,$child);
	$child = $child->{'sibling'};
    }
    return @children;
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
	foreach my $parser ( keys(%Galacticus::Build::SourceTree::Hooks::parseHooks) ) {
	    &{$Galacticus::Build::SourceTree::Hooks::parseHooks{$parser}}($tree)
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
    # Get a list of all defined processors.
    my @processors = sort(keys(%Galacticus::Build::SourceTree::Hooks::processHooks));
    # Perform a topological sort to enforce dependencies.
    my %dependencies;
    foreach my $processor ( keys(%Galacticus::Build::SourceTree::Hooks::processDependencies) ) {
	foreach my $dependent ( @{$Galacticus::Build::SourceTree::Hooks::processDependencies{$processor}} ) {
	    push(@{$dependencies{$dependent}},$processor);
	}
    }
    my @processorsOrdered = &Sort::Topo::sort(\@processors,\%dependencies);
    # Run all defined processors on the tree.
    &{$Galacticus::Build::SourceTree::Hooks::processHooks{$_}}($tree,\%options)
	foreach ( @processorsOrdered );
    return $tree;
}

sub AnalyzeTree {
    # Get the tree.
    my $tree      = shift();
    my (%options) = @_;
    # Run all defined processors on the tree.
    &{$Galacticus::Build::SourceTree::Hooks::analyzeHooks{$_}}($tree,\%options)
	foreach ( sort(keys(%Galacticus::Build::SourceTree::Hooks::analyzeHooks)) );
    return $tree;
}

sub Parse_Unit {
    my $unit = shift();
    # Process any content unless this is a code block.
    if ( exists($unit->{'content'}) && $unit->{'type'} ne "code" ) {
	my @children = &Build_Children($unit->{'type'},$unit->{'content'},$unit->{'source'},$unit->{'line'});
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
    my $type       = shift();
    my $codeText   = shift();
    my $source     = shift();
    my $lineNumber = shift();
    # Initialize an empty array of children.
    my @children = ();
    # Initialize current raw code buffer.
    my $rawCode;
    my $rawLineNumber = $lineNumber+1;
    # Initialize custom openers;
    my %unitOpeners = %Fortran::Utils::unitOpeners;
    $unitOpeners{'contains'} = { unitName => -1, regEx => qr/^\s*contains\s*$/ };
    # Only parse "module procedure" openers within a "contains" node.
    ## Within a contains section of a submodule "module procedure" is the opener of a subroutine/function whose interface was declared in the corresponding module.
    ## Outside of a contains section "module procedure" can appear within an interface. 
    delete($unitOpeners{'moduleProcedure'})
	unless ( $type eq "contains" );
    # Connect a file handle to the code text.
    ## NOTE: There seems to be a change in how Unicode charcaters are handled here between different Perl versions. The following
    ## is a hack which switches between the methods that work depending on Perl version. Note that the transition version is a
    ## guess. Known working versions are as follows:
    ##
    ### Version  Method
    ### 5.10.1   Old
    ### 5.26.3   New
    my $code;
    if ( $^V > v5.18 ) {
	# New method.
 	open($code,"<",\$codeText);
    } else {
	# Old method.
	my $codeTextDecoded = decode(q{utf8},$codeText);
	open($code,"<",\$codeTextDecoded);
    }
    # Read lines.
    do {
	# Get a line.
	my ($rawLine, $processedLine, $bufferedComments);
	&Fortran::Utils::Get_Fortran_Line($code,$rawLine,$processedLine,$bufferedComments);
	$lineNumber += $rawLine =~ tr/\n//;
	# Check for unit opening.
	my $unitFound;	
	foreach my $unitType ( keys(%unitOpeners) ) {
	    if ( my @matches = ( $processedLine =~ $unitOpeners{$unitType}->{'regEx'} ) ) {		
		# A match is found, extract the name of the unit, store current line number.
		my $lineNumberOpener = $lineNumber;
	     	my $unitName         = $matches[$unitOpeners{$unitType}->{'unitName'}]
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
			type       => "code"        ,
			content    => $rawCode      ,
			parent     => undef()       ,
			firstChild => undef()       ,
			sibling    => undef()       ,
			source     => $source       ,
			line       => $rawLineNumber
		    }
		    )
		    if ( $rawCode );
		undef($rawCode);
		$rawLineNumber = $lineNumber+1;
		# Read ahead until the matching unit closer is found, storing the unit code.
		my $closerRegExp;
		if ( $unitType eq "contains" ) {
		    $closerRegExp = qr/(?!)/; # Never match.
		} else {
		    $closerRegExp = $Fortran::Utils::unitClosers{$unitType}->{'regEx'};
		}
		do {
		    # Get the next line.
		    if ( eof($code) ) {
			$rawLine       = "";
			$processedLine = "";
		    } else {
			&Fortran::Utils::Get_Fortran_Line($code,$rawLine, $processedLine, $bufferedComments);
			$lineNumber += $rawLine =~ tr/\n//;
		    }
		    if ( ( $unitType eq "contains" && eof($code) ) || ( my @matches = ( $processedLine =~ $closerRegExp ) ) ) {
			my $closeUnitName = $matches[$Fortran::Utils::unitClosers{$unitType}->{'unitName'}]
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
				type       => $unitType        ,
				opener     => $opener          ,
				parent     => undef()          ,
				firstChild => undef()          ,
				sibling    => undef()          ,
				source     => $source          ,
				line       => $lineNumberOpener
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
			    $rawLineNumber = $lineNumber+1;
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
     	    type       => "code"        ,
     	    content    => $rawCode      ,
     	    parent     => undef()       ,
     	    firstChild => undef()       ,
	    sibling    => undef()       ,
	    source     => $source       ,
	    line       => $rawLineNumber
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
	# print " [".$_->{'node'}->{'source'}.":".$_->{'node'}->{'line'}."]"
	#     if ( exists($_->{'node'}->{'source'}) && exists($_->{'node'}->{'line'}) );
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

sub Walk_Branch {
    # Do a depth-first walk of the tree.
    my $branchHead = shift();
    my $node       = shift();
    my $depthRef   = shift();
    # Walk the tree.
    if ( $node == $branchHead ) {
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
		if ( $node == $branchHead );
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
    my (%options) = @_;
    $options{'annotate'} = 1
	unless ( exists($options{'annotate'}) );
    my $serialization;
    my $currentNode = $node;
    while ( $currentNode ) {
	$serialization .= "!--> ".$currentNode->{'line'}." \"".$currentNode->{'source'}."\"\n"
	    if ( exists($currentNode->{'source'}) && exists($currentNode->{'line'}) && $options{'annotate'} );
	if ( $currentNode->{'type'} eq "code" ) {
	    $serialization .= $currentNode->{'content'}
	} else {
	    $serialization .= $currentNode->{'opener'}
	        if ( exists($currentNode->{'opener'}) );
	    $serialization .= &Serialize($currentNode->{'firstChild'},%options)
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

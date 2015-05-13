# Contains a Perl module which implements parsing of module uses in the Galacticus preprocessor system.

package ModuleUses;
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
require Fortran::Utils;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Hooks::parseHooks{'moduleUses'} = \&Parse_ModuleUses;

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
	    open(my $code,"<",\$node->{'content'});
	    do {
		# Get a line.
		&Fortran_Utils::Get_Fortran_Line($code,my $rawLine, my $processedLine, my $bufferedComments);
		# Determine if line is a module use line.
		my $isModuleUse = 0;
		$isModuleUse    = 1
		    if ( $processedLine =~ m/^\s*use\s*(,\s*intrinsic)?\s*(::)?\s*([a-zA-Z0-9_]+)\s*(,\s*only\s*:)?\s*([a-zA-Z0-9_,\s]+)?\s*$/ );
		# Accumulate raw text.
		if ( $isModuleUse == 1 ) {
		    $rawModuleUse .= $rawLine;
		    my $isIntrinsic = $1;
		    my $moduleName  = $3;
		    my $only        = $5;
		    if ( $isIntrinsic ) {
			$moduleUses->{$moduleName}->{'intrinsic'} = 1;
		    } else {
			$moduleUses->{$moduleName}->{'intrinsic'} = 0;
		    }
		    if ( $only && ! exists($moduleUses->{$moduleName}->{'all'}) ) {
			map {$moduleUses->{$moduleName}->{'only'}->{$_} = 1} split(/\s*,\s*/,$only);
		    } else {
			$moduleUses->{$moduleName}->{'all'} = 1;
		    }   
		} else {
		    $rawCode      .= $rawLine;
		}
		# Process code and module use blocks as necessary.		
		if ( ( $isModuleUse == 1 || eof($code) ) && $rawCode      ) {
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
			firstChild => undef()
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
	&SourceTree::InsertBeforeNode($node->{'firstChild'},[$usesNode]);
    }
    # Merge the module usages.
    foreach my $moduleName ( keys(%{$moduleUses->{'moduleUse'}}) ) {
	$usesNode->{'moduleUse'}->{$moduleName}->{'intrinsic'} = $moduleUses->{'moduleUse'}->{$moduleName}->{'intrinsic'};
	unless ( exists($usesNode->{'moduleUse'}->{$moduleName}->{'all'}) ) {
	    if ( exists($moduleUses->{'moduleUse'}->{$moduleName}->{'all'}) ) {
		$usesNode->{'moduleUse'}->{$moduleName}->{'all'} = 1;
		delete($usesNode->{'moduleUse'}->{$moduleName}->{'only'});
	    } else {
		$usesNode->{'moduleUse'}->{$moduleName}->{'only'}->{$_} = 1
		    foreach ( keys(%{$moduleUses->{'moduleUse'}->{$moduleName}->{'only'}}) );
	    }
	}
    }
    # Update the contained code.
    $usesNode->{'firstChild'}->{'content'} = undef();
    foreach my $moduleName ( keys(%{$usesNode->{'moduleUse'}}) ) {
	$usesNode->{'firstChild'}->{'content'} .= "   use";
	$usesNode->{'firstChild'}->{'content'} .= ", intrinsic"
	    if ( $usesNode->{'moduleUse'}->{$moduleName}->{'intrinsic'} );
	$usesNode->{'firstChild'}->{'content'} .= " :: ".$moduleName;
	$usesNode->{'firstChild'}->{'content'} .= ", only : ".join(", ",keys(%{$usesNode->{'moduleUse'}->{$moduleName}->{'only'}}))
	    if ( $usesNode->{'moduleUse'}->{$moduleName}->{'only'} );
  	$usesNode->{'firstChild'}->{'content'} .= "\n";
    }
}

1;

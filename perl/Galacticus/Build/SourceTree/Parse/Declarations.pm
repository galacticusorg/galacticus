# Contains a Perl module which implements parsing of variable declarations in the Galacticus preprocessor system.

package Galacticus::Build::SourceTree::Parse::Declarations;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Storable qw(dclone);
use Fortran::Utils;
use List::MoreUtils qw(first_index);

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::parseHooks{'declarations'} = \&Parse_Declarations;

sub Parse_Declarations {
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
	    # Read the code block, accumulating declarations as we go.
	    my $rawCode;
	    my $rawDeclaration;
	    my @declarations;
	    my $lineNumber         = exists($node->{'line'  }) ? $node->{'line'  } : 0        ;
	    my $source             = exists($node->{'source'}) ? $node->{'source'} : "unknown";
	    my $rawCodeLine        = $lineNumber;
	    my $rawDeclarationLine = $lineNumber;
	    my $isImplicitNone     = 0;
	    open(my $code,"<",\$node->{'content'});
	    do {
		# Get a line.
		&Fortran::Utils::Get_Fortran_Line($code,my $rawLine, my $processedLine, my $bufferedComments);
		# Determine if line is a declaration line.
		my $isDeclaration = 0;
		if ( $processedLine =~ m/^\s*implicit\s+none\s*$/ ) {
		    $isImplicitNone   = 1;
		    $isDeclaration    = 1;
		}
		my $declaration = &parseDeclaration($processedLine);
		$isDeclaration = 1
		    if ( defined($declaration) );
		# Accumulate raw text.
		if ( $isDeclaration == 1 ) {
		    $rawDeclaration .= $rawLine;
		    if ( $declaration ) {
			$declaration->{'line'} = $lineNumber;
			push(@declarations,$declaration);
		    }
		} else {
		    $rawCode        .= $rawLine;
		}
		# Process code and declaration blocks as necessary.		
		if ( ( $isDeclaration == 1 || eof($code) ) && $rawCode ) {
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
		    $rawCodeLine        = $lineNumber;
		    $rawDeclarationLine = $lineNumber;
		}
		if ( ( $isDeclaration == 0 || eof($code) ) && $rawDeclaration ) {
		    # Create a new node.
		    my $newNode =
		    {
			type         => "declaration",
			implicitNone => $isImplicitNone
		    };
		    @{$newNode->{'declarations'}} = @declarations;
		    $newNode->{'firstChild'} =
		    {
			type         => "code"             ,
			content      => $rawDeclaration    ,
			parent       => $newNode           ,
			sibling      => undef()            ,
			firstChild   => undef()            ,
			source       => $source            ,
			line         => $rawDeclarationLine
		    };
		    $newNodes[$#newNodes]->{'sibling'} = $newNode
			if ( scalar(@newNodes) > 0 );
		    # If the end of the code has been reached and we're in a code block, pop
		    # that code block from the children array before we push our declaration node.
		    my $codeNode = pop(@newNodes)
			if ( eof($code) && $isDeclaration == 0 );
		    push(
			@newNodes,
			$newNode
			);
		    push(@newNodes,$codeNode)
			if ( eof($code) && $isDeclaration == 0 );
		    # Reset implicit none status.
		    $isImplicitNone = 0;
		    # Reset the raw declaration text.
		    undef($rawDeclaration);
		    undef(@declarations  );
		    $rawCodeLine        = $lineNumber;
		    $rawDeclarationLine = $lineNumber;
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

sub BuildDeclarations {
    # Build Fortran code from declarations list.
    my $node =   shift() ;
    $node->{'firstChild'}->{'content'} = $node->{'implicitNone'} ? "implicit none\n" : "";
    foreach my $declaration ( @{$node->{'declarations'}} ) {
	my $declarationCode  = "";
	$declarationCode    .= "#ifdef ".$declaration->{'preprocessor'}."\n"
	    if ( exists($declaration->{'preprocessor'}) );
	$declarationCode    .= "  ";
	$declarationCode    .= "!\$ "
	    if ( exists($declaration->{'openMP'}) && $declaration->{'openMP'} );
	$declarationCode    .= $declaration->{'intrinsic'};
	if ( exists($declaration->{'type'}) && defined($declaration->{'type'}) ) {
	    my $hasParentheses = $declaration->{'type'} =~ m/^\(/ && $declaration->{'type'} =~ m/\)$/;
	    $declarationCode .= "("
		unless ( $hasParentheses );
	    $declarationCode .= $declaration->{'type'};
	    $declarationCode .= ")"
		unless ( $hasParentheses );
	}
	$declarationCode    .= ", ".join(", ",@{$declaration->{'attributes'}})
	    if ( exists($declaration->{'attributes'}) && $declaration->{'attributes'} && scalar(@{$declaration->{'attributes'}}) > 0 );
	$declarationCode    .= " :: ".join(", ",@{$declaration->{'variables'}})."\n";
	$declarationCode    .= " !\$omp threadprivate(".join(",",map {$_ =~ s/([a-zA-Z0-9_]+).*/$1/; $_} @{$declaration->{'variables'}}).")\n"
	    if ( exists($declaration->{'threadprivate'}) && $declaration->{'threadprivate'} );
	$declarationCode    .= "#endif\n"
	    if ( exists($declaration->{'preprocessor'}) );
	$node->{'firstChild'}->{'content'} .= $declarationCode;
    }
}

sub AddDeclarations {
    # Grab the node to add declarations to, and the new declarations to add.
    my $node         =   shift() ;
    my @declarations = @{shift()};
    # Locate a declarations node.
    my $declarationsNode;
    my $childNode = $node->{'firstChild'};
    my $usesNode;
    while ( $childNode ) {
	$declarationsNode = $childNode
	    if ( $childNode->{'type'} eq "declaration" && ! defined($declarationsNode) );
	$usesNode = $childNode
	    if ( $childNode->{'type'} eq "moduleUse"   && ! defined($usesNode        ) );
	$childNode = $childNode->{'sibling'};
    }
    unless ( $declarationsNode ) {
	# No declarations node exists - create one.
	$declarationsNode =
	{
	    type       => "declaration",
	    sibling    => undef()      ,
	    parent     => undef()
	};
	$declarationsNode->{'firstChild'} =
	{
	    type       => "code"   ,
	    content    => ""       ,
	    sibling    => undef()  ,
	    parent     => $declarationsNode,
	    firstChild => undef()
	};
	# Insert the node, after any module use node if one exists.
	if ( $usesNode ) {
	    &Galacticus::Build::SourceTree::InsertAfterNode ($usesNode            ,[$declarationsNode]);
	} else {	
	    &Galacticus::Build::SourceTree::InsertBeforeNode($node->{'firstChild'},[$declarationsNode]);
	}
    }
    # Add the declarations.
    push(@{$declarationsNode->{'declarations'}},@declarations);
    &BuildDeclarations($declarationsNode);
}

sub AddAttributes {
    # Grab the node to add declarations to, and the new declarations to add.
    my $node         =   shift() ;
    my $variableName =   shift() ;
    my @attributes   = @{shift()};
    # Locate the declarations node.
    my $childNode  = $node->{'firstChild'};
    my $declarationsFound;
    my $declarationFound;
    while ( $childNode ) {
	if ( $childNode->{'type'} eq "declaration" ) {
	    $declarationsFound = $childNode;
	    # Locate the variable in the list of declarations.
	    foreach my $declaration ( @{$childNode->{'declarations'}} ) {
		if ( grep {$_ eq lc($variableName)} @{$declaration->{'variables'}} ) {
		    $declarationFound = $declaration;
		    last;
		}
	    }
	    last
		if ( $declarationFound );
	}
 	$childNode = $childNode->{'sibling'};
    }
    die('Galacticus::Build::SourceTree::Parse::Declarations::AddAttributes: no declarations present in '.$node->{'type'}.' "'.$node->{'name'}.'"'       )
	unless ( $declarationsFound );
    die('Galacticus::Build::SourceTree::Parse::Declarations::AddAttributes: variable declaration ['.$variableName.'] not found in '.$node->{'type'}.' "'.$node->{'name'}.'"')
	unless ( $declarationFound  );
    # Modify the attributes.
    if ( scalar(@{$declarationFound->{'variables'}}) > 1 ) {
	# Extract out other variables and push into their own declaration.
	my $declarationCopy = dclone($declarationFound);
	my $index           = first_index {$_ eq lc($variableName)} @{$declarationCopy->{'variables'}};
	$declarationFound->{'variables'    } = [ $declarationCopy->{'variables'    }->[$index] ];
	$declarationFound->{'variableNames'} = [ $declarationCopy->{'variableNames'}->[$index] ];
	splice(@{$declarationCopy->{'variables'    }},$index,1);
	splice(@{$declarationCopy->{'variableNames'}},$index,1);
	push(@{$declarationsFound->{'declarations'}},$declarationCopy);
    }
    foreach my $attribute ( @attributes ) {
	unless ( grep {$_ eq $attribute} @{$declarationFound->{'attributes'}} ) {
	    push(@{$declarationFound->{'attributes'}},$attribute);
	}
    }
    # Rebuild the declarations.
    &BuildDeclarations($declarationsFound);
}

sub GetDeclaration {
    # Return a descriptor of a variable declaration.
    my $node         = shift();
    my $variableName = shift();
    # Locate the declarations node.
    my $childNode         = $node->{'firstChild'};
    my $declarationsFound;
    my $declarationFound;
    while ( $childNode ) {
	if ( $childNode->{'type'} eq "declaration" ) {
    	    $declarationsFound = 1;
	    # Locate the variable in the list of declarations.
	    foreach my $declaration ( @{$childNode->{'declarations'}} ) {
		my @variableNames = map {$_ =~ m/^([^=]+)\s*=/ ? $1 : $_} @{$declaration->{'variables'}};
		if ( grep {$_ eq lc($variableName)} @variableNames ) {
		    $declarationFound = dclone($declaration);
		    $declarationFound->{'variables'} = [ $variableName ];
		    last;
		}
	    }
	}
 	$childNode = $childNode->{'sibling'};
    }
    die('Galacticus::Build::SourceTree::Process::Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration: no declarations present'       )
	unless ( $declarationsFound );
    die('Galacticus::Build::SourceTree::Process::Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration: variable declaration for "'.$variableName.'" not found in node "'.$node->{'opener'}.'"')
	unless ( $declarationFound  );
    return $declarationFound;
}

sub DeclarationExists {
    # Test whether a named variable is declared
    my $node         = shift();
    my $variableName = shift();
    # Locate the declarations node.
    my $exists       = 0;
    my $childNode    = $node->{'firstChild'};
    while ( $childNode ) {
	if ( $childNode->{'type'} eq "declaration" ) {
	    # Locate the variable in the list of declarations.
	    foreach my $declaration ( @{$childNode->{'declarations'}} ) {
		if ( grep {lc($_) eq lc($variableName)} @{$declaration->{'variables'}} ) {
		    $exists = 1;
		    last;
		}
	    }
	}
 	$childNode = $childNode->{'sibling'};
    }
    return $exists;
}

sub parseDeclaration {
    # Parse an individual declaration line.
    my $line = shift();
    my $declaration;
    foreach ( sort(keys(%Fortran::Utils::intrinsicDeclarations)) ) {
	if ( my @matches = ( $line =~ $Fortran::Utils::intrinsicDeclarations{$_}->{'regEx'} ) ) {
	    my $intrinsic  = $Fortran::Utils::intrinsicDeclarations{$_}->{'intrinsic'};
	    my $type;
	    ($type         = $matches[$Fortran::Utils::intrinsicDeclarations{$_}->{'type'}]) =~ s/\(\s*([^\s]*)\s*\)/$1/
               if ( $matches[$Fortran::Utils::intrinsicDeclarations{$_}->{'type'}] );
	    my $openMP = $matches[$Fortran::Utils::intrinsicDeclarations{$_}->{'openmp'}] ? 1 : 0;
	    my $attributesText = $matches[$Fortran::Utils::intrinsicDeclarations{$_}->{'attributes'}];
	    $attributesText =~ s/^\s*,\s*//
		if  ( $attributesText );
	    my @attributes    = &Fortran::Utils::Extract_Variables($attributesText                                                    ,keepQualifiers => 1                );
	    my @variables     = &Fortran::Utils::Extract_Variables($matches[$Fortran::Utils::intrinsicDeclarations{$_}->{'variables'}],keepQualifiers => 1                );
	    my @variableNames = &Fortran::Utils::Extract_Variables($matches[$Fortran::Utils::intrinsicDeclarations{$_}->{'variables'}],keepQualifiers => 0, lowerCase => 0);
	    $declaration = {
		intrinsic     => $intrinsic     ,
		type          => $type          ,
		openMP        => $openMP        ,
		attributes    => \@attributes   ,
		variables     => \@variables    ,
		variableNames => \@variableNames
	    };
	    last;
	}
    }
    return $declaration;
}

1;

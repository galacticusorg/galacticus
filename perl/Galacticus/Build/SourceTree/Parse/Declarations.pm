# Contains a Perl module which implements parsing of variable declarations in the Galacticus preprocessor system.

package Declarations;
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
$Hooks::parseHooks{'declarations'} = \&Parse_Declarations;

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
	    open(my $code,"<",\$node->{'content'});
	    do {
		# Get a line.
		&Fortran_Utils::Get_Fortran_Line($code,my $rawLine, my $processedLine, my $bufferedComments);
		# Determine if line is a declaration line.
		my $isDeclaration = 0;
		$isDeclaration    = 1
		    if ( $processedLine =~ m/^\s*implicit\s+none\s*$/ );
		my $declaration;
		foreach ( keys(%Fortran_Utils::intrinsicDeclarations) ) {
		    if ( my @matches = ( $processedLine =~ $Fortran_Utils::intrinsicDeclarations{$_}->{'regEx'} ) ) {
			my $intrinsic  = $Fortran_Utils::intrinsicDeclarations{$_}->{'intrinsic'};
			my $type;
			($type         = $matches[$Fortran_Utils::intrinsicDeclarations{$_}->{'type'}]) =~ s/\((.*)\)/$1/
			    if ( $matches[$Fortran_Utils::intrinsicDeclarations{$_}->{'type'}] );
			my $attributesText = $matches[$Fortran_Utils::intrinsicDeclarations{$_}->{'attributes'}];
			$attributesText =~ s/^\s*,\s*//
			    if  ( $attributesText );
			my @attributes = &Fortran_Utils::Extract_Variables($attributesText,keepQualifiers => 1);
			my @variables  = &Fortran_Utils::Extract_Variables($matches[$Fortran_Utils::intrinsicDeclarations{$_}->{'variables'}]);
			$declaration = {
			    intrinsic  => $intrinsic  ,
			    type       => $type       ,
			    attributes => \@attributes,
			    variables  => \@variables 
			};
			$isDeclaration = 1;
		    }
		}
		# Accumulate raw text.
		if ( $isDeclaration == 1 ) {
		    $rawDeclaration .= $rawLine;
		    push(@declarations,$declaration)
			if ( $declaration );
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
		if ( ( $isDeclaration == 0 || eof($code) ) && $rawDeclaration ) {
		    # Create a new node.
		    my $newNode =
		    {
			type         => "declaration"
		    };
		    @{$newNode->{'declarations'}} = @declarations;
		    $newNode->{'firstChild'} =
		    {
			type         => "code"         ,
			content      => $rawDeclaration,
			parent       => $newNode       ,
			sibling      => undef()        ,
			firstChild   => undef()
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
		    # Reset the raw module use text.
		    undef($rawDeclaration);
		    undef(@declarations  );
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
	    if ( $childNode->{'type'} eq "declaration" );
	$usesNode = $childNode
	    if ( $childNode->{'type'} eq "moduleUse"   );
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
	# Inert the node, after any module use node if one exists.
	if ( $usesNode ) {
	    &SourceTree::InsertAfterNode ($usesNode            ,[$declarationsNode]);
	} else {	
	    &SourceTree::InsertBeforeNode($node->{'firstChild'},[$declarationsNode]);
	}
    }
    # Add the declarations.
    push(@{$declarationsNode->{'declarations'}},@declarations);
    foreach ( @declarations ) {
	my $declarationCode  = "  ".$_->{'intrinsic'};
	$declarationCode    .= "(".$_->{'type'}.")"
	    if ( exists($_->{'type'}) && $_->{'type'} );
	$declarationCode    .= ", ".join(", ",@{$_->{'attributes'}})
	    if ( exists($_->{'attributes'}) && $_->{'attributes'} && scalar(@{$_->{'attributes'}}) > 0 );
	$declarationCode    .= " :: ".join(", ",@{$_->{'variables'}})."\n";
	$declarationsNode->{'firstChild'}->{'content'} .= $declarationCode;
    }
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
		if ( grep {$_ eq lc($variableName)} @{$declaration->{'variables'}} ) {
		    $declarationFound = $declaration;
		    $declarationFound->{'variables'} = [ $variableName ];
		    last;
		}
	    }
	}
 	$childNode = $childNode->{'sibling'};
    }
    die('Galacticus::Build::SourceTree::Process::Declarations::GetDeclaration: no declarations present'      )
	unless ( $declarationsFound );
    die('Galacticus::Build::SourceTree::Process::Declarations::GetDeclaration: variable declaration not found')
	unless ( $declarationFound  );
    return $declarationFound;
}

1;

# Contains a Perl module which implements processing of generic directives.

package Galacticus::Build::SourceTree::Process::Generics;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use Storable qw(dclone);
use Scalar::Util qw(reftype);
use List::ExtraUtils;
use Galacticus::Build::SourceTree::Parse::Declarations;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'generics'} = \&Process_Generics;

sub Process_Generics {
    # Get the tree.
    my $tree = shift();
    # Walk the tree.
    my $node  = $tree;
    my $depth = 0;    
    while ( $node ) {
	if ( $node->{'type'} eq "generic" ) {
	    $node->{'directive'}->{'processed'} = 1;
	    # Build regExs.
	    my $genericRegEx   = qr/\{$node->{'directive'}->{'identifier'}¦.*\}/;
	    # Iterate over all sibling nodes.
	    my $sibling = $node->{'sibling'};
	    while ( $sibling ) {
		# Get a breadth-first stack of the tree.
		my @stack = ();
		@stack = &Galacticus::Build::SourceTree::StackIt($sibling,\@stack,-1);		
		# Walk the tree, breadth first.
		my $depthPrevious = 0;
		while ( scalar(@stack) > 0 ) {
		    my $subTreeNode = pop(@stack);
		    if ( exists($subTreeNode->{'node'}->{'opener'}) && $subTreeNode->{'node'}->{'opener'} =~ $genericRegEx ) {
			# Initialize copies of the subtree.
			my @copies;
			# Iterate over instances.
			foreach my $instance ( &List::ExtraUtils::as_array($node->{'directive'}->{'instance'}) ) {
			    # Make a copy.
			    my $copy = dclone($subTreeNode->{'node'});
			    $copy->{'parent' } = undef();
			    $copy->{'sibling'} = undef();
			    # Walk the copied subtree.
			    my $copyNode  = $copy;
			    my $copyDepth = -1;			    
			    while ( $copyNode ) {
				# Replace generic tags in opener, closer, and name.
				foreach my $element ( 'opener', 'closer', 'name' ) {
				    if ( exists($copyNode->{$element}) ) {
					$copyNode->{$element} = &ReplaceGeneric           ($copyNode->{$element},$node->{'directive'}->{'identifier'},$instance,$_)
					    foreach ( sort(keys(%{$instance})) );
					$copyNode->{$element} = &ReplaceGenericConditional($copyNode->{$element},$node->{'directive'}->{'identifier'},$instance   );
				    }
				}
				# Replace generic tags in code.
				if ( exists($copyNode->{'content'}) ) {
				    my $newCode;
				    open(my $code,"<",\$copyNode->{'content'});
				    while ( my $line = <$code> ) {
					$line = &ReplaceGeneric           ($line,$node->{'directive'}->{'identifier'},$instance,$_)
					    foreach ( sort(keys(%{$instance})) );
					$line = &ReplaceGenericConditional($line,$node->{'directive'}->{'identifier'},$instance   );
					$newCode .= $line;
				    }
				    close($code);
				    $copyNode->{'content'} = $newCode;
				    # Check if the code is part of a variable declaration.
				    if ( $copyNode->{'parent'}->{'type'} eq "declaration" ) {
					my $nodeCopy = dclone($copyNode);
					my $treeTmp  = {type => "null", firstChild => $nodeCopy, sibling => undef(), parent => undef()};
					$nodeCopy->{'parent'} = $treeTmp;
					&Galacticus::Build::SourceTree::Parse::Declarations::Parse_Declarations($treeTmp);
					&Galacticus::Build::SourceTree::ReplaceNode($copyNode->{'parent'},[$treeTmp->{'firstChild'}]);
					$copyNode = $treeTmp->{'firstChild'};
				    }
				}
				# Move to the next node in the copied tree.
				$copyNode = &Galacticus::Build::SourceTree::Walk_Tree($copyNode,\$copyDepth);
			    }
			    # Reparse the new content.
			    my $copyReparsed = &Galacticus::Build::SourceTree::ParseCode(&Galacticus::Build::SourceTree::Serialize($copy),$tree->{'name'}, instrument => 0 ,reinstateBlocks => 1);
			    # Push copy to list of copies.
			    push(@copies,$copyReparsed);
			}
			# Replace the subtree with the array of copies.
			&Galacticus::Build::SourceTree::ReplaceNode($subTreeNode->{'node'},\@copies);
			# Skip over entries in the stack that belonged to the subtree.
			my $depthTarget = $subTreeNode->{'depth'};
			if ( scalar(@stack) > 0 ) {
			    do {
				$subTreeNode = pop(@stack);
			    } until ( $subTreeNode->{'depth'} == $depthTarget || scalar(@stack) == 0 );
			    push(@stack,$subTreeNode)
				if ( $subTreeNode->{'depth'} == $depthTarget );
			}
		    } elsif ( exists($subTreeNode->{'node'}->{'content'}) ) {
			# Check that the code is not contained in a node with generic opener.
			my $processContent = 1;
			my $parent         = $subTreeNode->{'node'}->{'parent'};			
			while ( $parent ) {
			    $processContent = 0
				if ( exists($parent->{'opener'}) && $parent->{'opener'} =~ $genericRegEx );
			    $parent = $parent->{'parent'};
			}
			# Replace generic tags in the code.
			if ( $processContent ) {
			    my $newCode;
			    open(my $code,"<",\$subTreeNode->{'node'}->{'content'});
			    while ( my $line = <$code> ) {
				if ( $line =~ $genericRegEx ) {
				    # Iterate over instances.
				    foreach my $instance ( &List::ExtraUtils::as_array($node->{'directive'}->{'instance'}) ) {
					my $copiedLine = $line;
					$copiedLine = &ReplaceGeneric           ($copiedLine,$node->{'directive'}->{'identifier'},$instance,$_)
					    foreach ( sort(keys(%{$instance})) );
					$copiedLine = &ReplaceGenericConditional($copiedLine,$node->{'directive'}->{'identifier'},$instance   );
					$newCode .= $copiedLine;
				    }
				} else {
				    $newCode .= $line;
				}
			    }
			    close($code);
			    $subTreeNode->{'node'}->{'content'} = $newCode;
			    # Check if the code is part of a variable declaration.
			    if ( $subTreeNode->{'node'}->{'parent'}->{'type'} eq "declaration" ) {
				my $nodeCopy = dclone($subTreeNode->{'node'});
				my $treeTmp  = {type => "null", firstChild => $nodeCopy, sibling => undef(), parent => undef()};
				$nodeCopy->{'parent'} = $treeTmp;
				&Galacticus::Build::SourceTree::Parse::Declarations::Parse_Declarations($treeTmp);
				&Galacticus::Build::SourceTree::ReplaceNode($subTreeNode->{'node'}->{'parent'},[$treeTmp->{'firstChild'}]);
			    }
			}
		    }
		}
		# Next sibling.
		$sibling = $sibling->{'sibling'};
	    }
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

sub ReplaceGeneric {
    # Replace a generic directive.
    my $line       = shift();
    my $identifier = shift();
    my $instance   = shift();
    my $modifier   = shift();
    # Replace defined modifiers.
    if ( $instance->{$modifier} =~ m/^regEx\x{0A6}(.*)\x{0A6}(.*)\x{0A6}/ ) {
	my $from  = $1;
	my $to    = $2;
	my $regEx = qr/\{$identifier¦$modifier¦$from\}/s;
	$line      =~ s/$regEx/eval qq{"$to"}/eg;
    } else {
	$line      =~ s/\{$identifier(¦$modifier)?\}/$instance->{$modifier}/g;
    }
    # Return the modified line.   
    return $line;
}

sub ReplaceGenericConditional {
    # Replace a generic conditional directive.
    my $line       = shift();
    my $identifier = shift();
    my $instance   = shift();
    # Replace any inline conditionals.
    while ( $line =~ m/\{$identifier¦match¦([^¦]*)¦([^¦]*)¦([^¦]*)\}/ ) {
	my $instanceRegEx = qr/$1/;
	my $match         = $2;
	my $noMatch       = $3;
	if ( $instance->{'label'} =~ $instanceRegEx ) {
	    $line =~ s/\{$identifier¦match¦([^¦]*)¦([^¦]*)¦([^¦]*)\}/$match/;
	} else {
	    $line =~ s/\{$identifier¦match¦([^¦]*)¦([^¦]*)¦([^¦]*)\}/$noMatch/;
	}
    }
    # Return the modified line.   
    return $line;
}

1;

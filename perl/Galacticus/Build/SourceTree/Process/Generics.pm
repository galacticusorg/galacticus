# Contains a Perl module which implements processing of generic directives.

package Generics;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use Data::Dumper;
use XML::Simple;
use LaTeX::Encode;
use Storable qw(dclone);
use Scalar::Util qw(reftype);
require List::ExtraUtils;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Hooks::processHooks{'generics'} = \&Process_Generics;

sub Process_Generics {
    # Get the tree.
    my $tree = shift();
    # Walk the tree.
    my $node  = $tree;
    my $depth = 0;    
    while ( $node ) {
	if ( $node->{'type'} eq "generic" ) {
	    # Build regExs.
	    my $genericRegEx   = qr/\{$node->{'directive'}->{'identifier'}¦.*\}/;
	    # Iterate over all sibling nodes.
	    my $sibling = $node->{'sibling'};
	    while ( $sibling ) {
		# Get a breadth-first stack of the tree.
		my @stack = ();
		@stack = &SourceTree::StackIt($sibling,\@stack,-1);		
		# Walk the tree, breadth first.
		my $depthPrevious = 0;
		while ( scalar(@stack) > 0 ) {
		    my $subTreeNode = pop(@stack);
		    if ( exists($subTreeNode->{'node'}->{'opener'}) && $subTreeNode->{'node'}->{'opener'} =~ $genericRegEx ) {
			# Initialize copies of the subtree.
			my @copies;
			# Iterate over instances.
			foreach my $instance ( &ExtraUtils::as_array($node->{'directive'}->{'instance'}) ) {
			    # Make a copy.
			    my $copy = dclone($subTreeNode->{'node'});
			    $copy->{'parent' } = undef();
			    $copy->{'sibling'} = undef();
			    # Walk the copied subtree.
			    my $copyNode  = $copy;
			    my $copyDepth = -1;			    
			    while ( $copyNode ) {
				# Replace generic tags in opener and closer.
				foreach my $element ( 'opener', 'closer' ) {
				    if ( exists($copyNode->{$element}) ) {
					$copyNode->{$element} = &ReplaceGeneric($copyNode->{$element},$node->{'directive'}->{'identifier'},$instance,$_)
					    foreach ( keys(%{$instance}) );
				    }
				}
				# Replace generic tags in code.
				if ( exists($copyNode->{'content'}) ) {
				    my $newCode;
				    open(my $code,"<",\$copyNode->{'content'});
				    while ( my $line = <$code> ) {
					$line = &ReplaceGeneric($line,$node->{'directive'}->{'identifier'},$instance,$_)
					    foreach ( keys(%{$instance}) );
					$newCode .= $line;
				    }
				    close($code);				    
				    $copyNode->{'content'} = $newCode;
				}
				# Move to the next node in the copied tree.
				$copyNode = &SourceTree::Walk_Tree($copyNode,\$copyDepth);
			    }
			    # Push copy to list of copies.
			    push(@copies,$copy);
			}
			# Replace the subtree with the array of copies.
			&SourceTree::ReplaceNode($subTreeNode->{'node'},\@copies);
			# Skip over entries in the stack that belonged to the subtree.
			my $depthTarget = $subTreeNode->{'depth'};
			do {
			    $subTreeNode = pop(@stack);
			} until ( $subTreeNode->{'depth'} == $depthTarget || scalar(@stack) == 0 );
			push(@stack,$subTreeNode)
			    if ( $subTreeNode->{'depth'} == $depthTarget );
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
				    foreach my $instance ( &ExtraUtils::as_array($node->{'directive'}->{'instance'}) ) {
					my $copiedLine = $line;
					$copiedLine = &ReplaceGeneric($copiedLine,$node->{'directive'}->{'identifier'},$instance,$_)
					    foreach ( keys(%{$instance}) );
					$newCode .= $copiedLine;
				    }
				} else {
				    $newCode .= $line;
				}
			    }
			    close($code);
			    $subTreeNode->{'node'}->{'content'} = $newCode;
			}
		    }
		}
		# Next sibling.
		$sibling = $sibling->{'sibling'};
	    }
	}
	$node = &SourceTree::Walk_Tree($node,\$depth);
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

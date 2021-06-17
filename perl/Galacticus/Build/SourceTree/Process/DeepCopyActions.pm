# Contains a Perl module which implements processing of deepCopyActions directives.

package Galacticus::Build::SourceTree::Process::DeepCopyActions;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Text::Template 'fill_in_string';

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'deepCopyActions'} = \&Process_DeepCopyActions;

sub Process_DeepCopyActions {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Initialize the directives and classes.
    my @directiveNodes;
    my %classes;
    # Initialize deep copy actions database.
    my $deepCopyActions;
    # Walk the tree.
    my $node  = $tree;
    my $depth = 0    ;
    while ( $node ) {
	# Capture deepCopyActions directives.
	if ( $node->{'type'} eq "deepCopyActions" ) {
	    # Extract the directive.
	    push(@directiveNodes,$node);	    
	    # Get state storables database if we do not have it.
	    $deepCopyActions = $xml->XMLin($ENV{'BUILDPATH'}."/deepCopyActions.xml")
		unless ( $deepCopyActions );	    
	}
	# Capture derived type definitions.
	if ( $node->{'type'} eq "type" ) {
	    # Parse class openers to find dependencies.
	    if ( $node->{'opener'} =~ m/^\s*type\s*(,\s*(abstract)\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*([a-z0-9_]+)\s*$/i ) {
		my $type     = $5;
		my $extends  = $3;
		my $abstract = defined($2);
		$classes{$type} =
		{
		     node     => $node   ,
		     extends  => $extends,
		     abstract => $abstract
		};
	    } elsif ( $node->{'opener'} !~ /\{/ ) { # Ignore types that use generics as this module does not currently support them.
		die("Galacticus::Build::SourceTree::Process::DeepCopyActions::Process_DeepCopyActions: unable to parse type definition opener");
	    }
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
    # Process any deepCopyActions directives.
    foreach my $directiveNode ( @directiveNodes ) {
	my $directive = $directiveNode->{'directive'};
	# Assert that class must exist.
	die("Galacticus::Build::SourceTree::Process::DeepCopyActions::Process_DeepCopyActions: class '".$directive->{'class'}."' not found")
	    unless ( exists($classes{$directive->{'class'}}) );
	# Begin building deep copy actions function code.
	$code::className   = $directive->{'class'};
	my $deepCopyAction = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine {$className}DeepCopyActions(self)
 !!\{
 Perform actions needed for deep copy of this object.
 !!\}
 implicit none
 class({$className}), intent(inout) :: self
 select type (self)
CODE
	# Scan all known classes, finding all which derive from the base class.
	foreach my $className ( sort(keys(%classes)) ) {
	    my $matches = 0;
	    my $parentClassName = $className;
	    while ( defined($parentClassName) ) {
		if ( $parentClassName eq $directive->{'class'} ) {
		    $matches = 1;
		    last;
		}
		$parentClassName = $classes{$parentClassName}->{'extends'};
	    }
	    if ( $matches && ! $classes{$className}->{'abstract'} ) {
		# This class is derived from the class of interest, so perform actions.
		$deepCopyAction .= " type is (".$className.")\n";
		# Search the class node for declarations.
		my @methodCalls;
		my $parentClassName = $className;
		while ( defined($parentClassName) ) {
		    my $classNode = $classes{$parentClassName}->{'node'}->{'firstChild'};
		    while ( $classNode ) {
			if ( $classNode->{'type'} eq "declaration" ) {
			    # "setTo" actions - simply set a variable to a value.
			    if ( exists($directive->{$parentClassName}) && exists($directive->{$parentClassName}->{'setTo'}) ) {
				foreach my $setTo ( &List::ExtraUtils::as_array($directive->{$parentClassName}->{'setTo'}) ) {
				    my @variables = split(/\s*,\s*/,$setTo->{'variables'});
				    $deepCopyAction .= " self%".$_."=".$setTo->{'state'}."\n"
					foreach ( @variables );
				}
			    }
			    # "methodCall" actions - simply call a method.
			    if ( exists($directive->{$parentClassName}) && exists($directive->{$parentClassName}->{'methodCall'}) ) {
				foreach my $methodCall ( &List::ExtraUtils::as_array($directive->{$parentClassName}->{'methodCall'}) ) {				    
				    push(@methodCalls," call self%".$methodCall->{'method'}."(".(exists($methodCall->{'arguments'}) ? $methodCall->{'arguments'} : "").")");
				}
			    }
			}
			$classNode = $classNode->{'sibling'};
		    }
		    $parentClassName = $classes{$parentClassName}->{'extends'};
		}
		$deepCopyAction .= join("\n",@methodCalls)."\n"
		    if ( @methodCalls );
	    }
	}
	# Close function.
	$deepCopyAction .= fill_in_string(<<'CODE', PACKAGE => 'code');
  end select
  return
end subroutine {$className}DeepCopyActions
CODE
	# Insert type-bindings.
	my $content = fill_in_string(<<'CODE', PACKAGE => 'code');
    !![
    <methods>
     <method method="deepCopyActions" description="Perform actions needed for deep copy of this object."/>
    </methods>
    !!]
    procedure :: deepCopyActions => {$className}DeepCopyActions
CODE
	my $treeContent = &Galacticus::Build::SourceTree::ParseCode($content,'null');
	my @childrenContent = &Galacticus::Build::SourceTree::Children($treeContent);
	&Galacticus::Build::SourceTree::InsertPostContains($classes{$directive->{'class'}}->{'node'},\@childrenContent);

	# Insert code.
	my $treeDeepCopyActions = &Galacticus::Build::SourceTree::ParseCode($deepCopyAction,'null');
	&Galacticus::Build::SourceTree::ProcessTree($treeDeepCopyActions);
	my @childrenDeepCopyActions = &Galacticus::Build::SourceTree::Children($treeDeepCopyActions);
	&Galacticus::Build::SourceTree::InsertPostContains($directiveNode->{'parent'},\@childrenDeepCopyActions);

    }
}

1;

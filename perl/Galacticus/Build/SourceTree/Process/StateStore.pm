# Contains a Perl module which implements processing of stateStore directives.

package Galacticus::Build::SourceTree::Process::StateStore;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Build::SourceTree::Process::SourceIntrospection;
use Galacticus::Build::SourceTree::Parse::Declarations;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'stateStore'} = \&Process_StateStore;

sub Process_StateStore {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Walk the tree.
    my $node       = $tree;
    my $moduleNode        ;
    my $depth      = 0    ;
    while ( $node ) {
	# Capture stateStore directives.
	if ( $node->{'type'} eq "stateStore" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} = 1;
	    my @variables = split(" ",$node->{'directive'}->{'variables'});
	    my $code;
	    foreach my $variable ( @variables ) {
		$code .= "write (stateFile) associated(".$variable.")\n";
		$code .= "if (associated(".$variable.")) call ".$variable."%stateStore(stateFile,gslStateFile,stateOperationID)\n";
	    }
	    # Create a new node.
	    my $newNode =
	    {
		type       => "code",
		content    => $code ,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	}
	# Capture stateRestore directives.
	if ( $node->{'type'} eq "stateRestore" && ! $node->{'directive'}->{'processed'} ) {
	    $node->{'directive'}->{'processed'} = 1;
	    my @variables = split(" ",$node->{'directive'}->{'variables'});
	    my $code;
	    foreach my $variable ( @variables ) {
		$code .= "read (stateFile) wasAllocated_\n";
		$code .= "if (wasAllocated_) then\n";
		$code .= "   if (.not.associated(".$variable.")) call Error_Report(\"'".$variable."' was stored, but is now not allocated\"//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		$code .= "   call ".$variable."%stateRestore(stateFile,gslStateFile,stateOperationID)\n";
		$code .= "end if\n";
	    }
	    # Create a new node.
	    my $newNode =
	    {
		type       => "code",
		content    => $code ,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Add state variable.
	    my @stateVariable =
		(
		 {
		     intrinsic => "logical",
		     variables => [ "wasAllocated_" ]
		 }
		);
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@stateVariable);

	    # Add required modules.
	    my $moduleUses =
	    {
		type      => "moduleUse",
		moduleUse =>
		{
		    Error => {intrinsic => 0, only => { "Error_Report" => 1}}
		}
	    };
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$moduleUses);
	}
	# Walk to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

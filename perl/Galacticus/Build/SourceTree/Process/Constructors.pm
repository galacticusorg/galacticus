# Contains a Perl module which implements processing of constructor directives.

package Galacticus::Build::SourceTree::Process::Constructors;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use List::ExtraUtils;
use Galacticus::Build::SourceTree::Parse::Declarations;
$XML::Simple::PREFERRED_PARSER = "XML::Parser";

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'constructors'} = \&Process_Constructors;

sub Process_Constructors {
    # Get the tree.
    my $tree  = shift();
    # Get an XML parser.
    my $xml   = new XML::Simple();
    # Initialize state storables database.
    my $stateStorables;
    # Determine if debugging output is required.
    my $debugging = exists($ENV{'GALACTICUS_OBJECTS_DEBUG'}) && $ENV{'GALACTICUS_OBJECTS_DEBUG'} eq "yes";
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "constructorAssign" && ! $node->{'directive'}->{'processed'} ) {
	    # Assert that our parent is a function.
	    die("Process_Constructors: parent node must be a function")
		unless ( $node->{'parent'}->{'type'} eq "function" || $node->{'parent'}->{'type'} eq "moduleProcedure" );
	    # Get state storables database if we do not have it.
	    $stateStorables = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml")
		unless ( $stateStorables );
	    # Determine if automatic allocation of variables should be performed.
	    my $allocate = exists($node->{'directive'}->{'allocate'}) ? $node->{'directive'}->{'allocate'} : "yes";
	    # Determine function return value name.
	    my $returnValueLabel;
	    if ( $node->{'parent'}->{'opener'} =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$/ ) {
		$returnValueLabel = $1;
	    } else {
		$returnValueLabel = $node->{'parent'}->{'name'}; 
	    }
	    # Generate source code for the assignment.
	    $node->{'directive'}->{'processed'} = 1;
	    my $assignmentSource = "  ! Auto-generated constructor assignment\n";
	    (my $variables = $node->{'directive'}->{'variables'}) =~ s/^\s*(.*?)\s*$/$1/;
	    foreach ( grep {$_ ne ""} split(/\s*,\s*/,$variables) ) {
		my $matches      = $_ =~ m/^(\*??)([a-zA-Z0-9_]+)/;
		my $isPointer    = $1 eq "*";
		my $argumentName = $2;
		my $assigner     = $isPointer ? " => " : "=";
		my $hasDefault   = $_ =~ m/=\s*(.+)/;
		my $default      = $hasDefault ? $1 : undef();
		die("Galacticus::Build::SourceTree::Process::Constructor::Process_Constructors(): syntax error")
		    unless ( $matches );
		# Get the variable declaration.
		my $declaration;
		if ( $node->{'parent'}->{'type'} eq "moduleProcedure" ) {
		    
		} else {
		    $declaration = &Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration($node->{'parent'},$argumentName);
		}
		# Detect optional arguments.
		my $optional    = (grep {$_ eq "optional"} @{$declaration->{'attributes'}}) ? "if (present(".$argumentName.")) " : "";
		# Detect allocatable objects.
		if ( $allocate eq "yes" && grep {$_ =~ m/dimension\s*\([:,]+\)/} @{$declaration->{'attributes'}} ) {
		    # Determine the rank of the variable.
		    my $rank = join("",map {$_ =~ m/dimension\s*\(([:,]+)\)/ ? $1 : ""} @{$declaration->{'attributes'}}) =~ tr/://;
		    $assignmentSource .= "   ".$optional."allocate(".$returnValueLabel."%".$argumentName."(".join(",",map {"size(".$argumentName.",dim=".$_.")"} 1..$rank)."))\n";
		}
		# Build the assignment.
		if ( $optional eq "" ) {
		    $assignmentSource   .= "   "                    .$returnValueLabel."%".$argumentName.$assigner.$argumentName."\n";
		} elsif ( $hasDefault ) {
		    $assignmentSource   .= "   ".$optional." then\n".$returnValueLabel."%".$argumentName.$assigner.$argumentName."\nelse\n".$returnValueLabel."%".$argumentName.$assigner.$default."\nend if\n";
		} else {
		    $assignmentSource   .= "   ".$optional          .$returnValueLabel."%".$argumentName.$assigner.$argumentName."\n";
		}
		# Detect functionClass objects and increment their reference count.
		if ( 
		    ( $declaration->{'intrinsic'} eq "type" || $declaration->{'intrinsic'} eq "class" )
		    &&
		    $isPointer
		    ) {
		    my $type = lc($declaration->{'type'});
		    $type =~ s/\s//g;
		    if ( grep {lc($_) eq $type} (keys(%{$stateStorables->{'functionClasses'}}),@{$stateStorables->{'functionClassInstances'}})) {
			$assignmentSource .= "   ".$optional." then\n"
			    unless ( $optional eq "" );
			$assignmentSource .= "   if (associated(".$returnValueLabel."%".$argumentName.")) "
			    if ( $isPointer );
			$assignmentSource .= " call ".$returnValueLabel."%".$argumentName."%referenceCountIncrement()\n";
			if ( $debugging ) {
			    $assignmentSource .= "   if (debugReporting.and.mpiSelf\%isMaster()) then\n";
			    $assignmentSource .= "   ".$optional." call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): [".$argumentName."] : ".$returnValueLabel." : ')//debugStackGet()//' : '//loc(".$returnValueLabel."%".$argumentName.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'},compact => 1).",verbosityLevelSilent)\n";
			    $assignmentSource .= "   end if\n";
			}
			$assignmentSource .= "   end if\n"
			    unless ( $optional eq "" );
		    } elsif ( $declaration->{'type'} eq "*" ) {
			$assignmentSource .= "   ".$optional." then\n"
			    unless ( $optional eq "" );
			$assignmentSource .= "   if (associated(".$returnValueLabel."%".$argumentName.")) then\n"
			    if ( $isPointer );
			$assignmentSource .= "select type(s__ => ".$returnValueLabel."%".$argumentName.")\n";
			$assignmentSource .= "class is (functionClass)\n";
			$assignmentSource .= " call s__%referenceCountIncrement()\n";
			if ( $debugging ) {
			    $assignmentSource .= "   if (debugReporting.and.mpiSelf\%isMaster()) then\n";
			    $assignmentSource .= "   ".$optional." call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): [".$argumentName."] : ".$returnValueLabel." : ')//debugStackGet()//' : '//loc(".$returnValueLabel."%".$argumentName.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'},compact => 1).",verbosityLevelSilent)\n";
			    $assignmentSource .= "   end if\n";
			}
			$assignmentSource .= "end select\n";
			$assignmentSource .= "   end if\n"
			    if ( $isPointer );
			$assignmentSource .= "   end if\n"
			    unless ( $optional eq "" );
	
		    }
		}
	    }
	    $assignmentSource   .= "  ! End auto-generated constructor assignment.\n\n";
	    # Create a new node.
	    my $newNode =
	    {
		type       => "code"           ,
		content    => $assignmentSource,
		firstChild => undef()          ,
		source     => "Galacticus::Build::SourceTree::Process::Constructor::Process_Constructors()",
		line       => 1
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Add modules for debugging.
	    if ( $debugging ) {
		my $usesNode =
		{
		    type      => "moduleUse",
		    moduleUse =>
		    {
			MPI_Utilities =>
			{
			    intrinsic => 0,
			    all       => 1
			},
			Display =>
			{
			    intrinsic => 0,
			    all       => 1
			},
			String_Handling    =>
			{
			    intrinsic => 0,
			    all       => 1
			},
			ISO_Varying_String =>
			{
			    intrinsic => 0,
			    all       => 1
			},
			Function_Classes   =>
			{
			    intrinsic => 0,
			    all       => 1
			},
		    }
		};
		&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    }
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

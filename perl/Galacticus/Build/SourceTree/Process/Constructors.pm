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

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'constructors'} = \&Process_Constructors;

sub Process_Constructors {
    # Get the tree.
    my $tree  = shift();
    # Get an XML parser.
    my $xml   = new XML::Simple();
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "constructorAssign" && ! $node->{'directive'}->{'processed'} ) {
	    # Assert that our parent is a function.
	    die("Process_Constructors: parent node must be a function")
		unless ( $node->{'parent'}->{'type'} eq "function" );
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
		my $assigner     = "=";
		my $argumentName = $_;
		if ( $_ =~ m/^\*(.*)/ ) {
		    $assigner     = " => ";
		    $argumentName = $1;
		} 
		# Get the variable declaration.
		my $declaration = &Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration($node->{'parent'},$argumentName);
		# Detect optional arguments.
		my $optional    = (grep {$_ eq "optional"} @{$declaration->{'attributes'}}) ? "if (present(".$argumentName.")) " : "";
		# Detect allocatable objects.
		if ( $allocate eq "yes" && grep {$_ =~ m/dimension\s*\([:,]+\)/} @{$declaration->{'attributes'}} ) {
		    # Determine the rank of the variable.
		    my $rank = join("",map {$_ =~ m/dimension\s*\(([:,]+)\)/ ? $1 : ""} @{$declaration->{'attributes'}}) =~ tr/://;
		    $assignmentSource .= "   allocate(".$returnValueLabel."%".$argumentName."(".join(",",map {"size(".$argumentName.",dim=".$_.")"} 1..$rank)."))\n";
		}
		# Build the assignment.
		$assignmentSource   .= "   ".$optional.$returnValueLabel."%".$argumentName.$assigner.$argumentName."\n";
	    }
	    $assignmentSource   .= "  ! End auto-generated constructor assignment.\n\n";
	    # Create a new node.
	    my $newNode =
	    {
		type       => "code"           ,
		content    => $assignmentSource,
		firstChild => undef()
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

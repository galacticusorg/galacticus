# Contains a Perl module which implements validation of parameter names in parameter-based constructors.

package Galacticus::Build::SourceTree::Process::InputParametersValidate;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use XML::Simple;
use LaTeX::Encode;
use List::ExtraUtils;
use Fortran::Utils;
use Galacticus::Build::Directives;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'inputParametersValidate'} = \&Process_InputParametersValidate;

sub Process_InputParametersValidate {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Get the file name.
    my $fileName;
    $fileName = $tree->{'name'}
        if ( $tree->{'type'} eq "file" );
    # Get code directive locations.
    my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
     # Walk the tree, looking for input parameter validation directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	# Look for inputParametersValidate directives and process them.
	if ( $node->{'type'} eq "inputParametersValidate" && ! $node->{'directive'}->{'processed'} ) {	    
	    # Record that this directive has been processed.
	    $node->{'directive'}->{'processed'} =  1;
	    # Determine the parameter source.
	    my $source = exists($node->{'directive'}->{'source'}) ? $node->{'directive'}->{'source'} : "globalParameters";
	    # Step through sibling nodes looking for input parameter directives.
	    my @objectBuilderNames;
	    my $sibling = $node->{'parent'}->{'firstChild'};
	    while ( $sibling ) {
		if ( $sibling->{'type'} eq "objectBuilder" ) {
		    my $objectBuilderSource = exists($sibling->{'directive'}->{'source'}) ? $sibling->{'directive'}->{'source'} : "globalParameters";
		    push(@objectBuilderNames,$sibling->{'directive'}->{'name'})
			if ( $objectBuilderSource eq $source );
		}
		$sibling = $sibling->{'sibling'};
	    }
	    # Generate the variable declaration.
	    my $variableName = exists($node->{'directive'}->{'label'}) ? $node->{'directive'}->{'label'} : "allowedParameterNames_";
	    unless ( &Galacticus::Build::SourceTree::Parse::Declarations::DeclarationExists($node->{'parent'},$variableName) ) {
		my $declaration =
		{
		    intrinsic  => "type",
		    type       => "varying_string",
		    attributes => [ "dimension(:)", "allocatable" ],
		    variables  => [ $variableName ]
		};
		&Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},[$declaration]);
	    }
	    # Add module usage.
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},{moduleUse => {ISO_Varying_String => {all => 1}}});
	    # Generate the validation code.
	    my $code;
	    my $result;
	    if      ( $node->{'parent'}->{'type'} eq "subroutine" ) {
		if ( exists($node->{'directive'}->{'target'}) ) {
		    $result = $node->{'directive'}->{'target'};
		} else {
		    die('Galacticus::Build::SourceTree::Process::InputParametersValidate::Process_InputParametersValidate: target must be specified');
		}
	    } elsif ( $node->{'parent'}->{'type'} eq "function"   ) {
		if ( $node->{'parent'}->{'opener'} =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$/ ) {
		    $result = $1;
		} else {
		    $result = $node->{'parent'}->{'name'};	       
		}
	    } else {
		die('Galacticus::Build::SourceTree::Process::InputParametersValidate::Process_InputParametersValidate: parent is neither function nor subroutine');
	    }
	    $code .= "   call ".$result."%allowedParameters(".$variableName.",'".$source."')\n";
	    $code .= "   call ".$_     ."%allowedParameters(".$variableName.",'parameters')\n"
		foreach ( @objectBuilderNames);
	    $code .= "   call ".$source."%checkParameters(".$variableName.")\n";
	    $code .= "   if (allocated(".$variableName.")) deallocate(".$variableName.")\n";
	    # Insert new code.
	    my $codeNode =
	    {
		type       => "code"           ,
		content    => $code            ,
		sibling    => undef()          ,
		parent     => $node->{'parent'},
		firstChild => undef()
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$codeNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

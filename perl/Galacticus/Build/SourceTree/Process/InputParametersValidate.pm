# Contains a Perl module which implements validation of parameter names in parameter-based constructors.

package Galacticus::Build::SourceTree::Process::InputParametersValidate;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
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
    # Initialize state storables database.
    my $stateStorables = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml");
    # Get the file name.
    my $fileName;
    $fileName = $tree->{'name'}
        if ( $tree->{'type'} eq "file" );
    # Walk the tree, looking for input parameter validation directives.
    my $node  = $tree;
    my $depth = 0;
    my $functionClassName;
    while ( $node ) {
	# Look for function class definitions.
	$functionClassName = $node->{'type'}
	    if ( grep {$node->{'type'}."Class" eq $_} keys(%{$stateStorables->{'functionClasses'}}) );
	# Look for inputParametersValidate directives and process them.
	if ( $node->{'type'} eq "inputParametersValidate" && ! $node->{'directive'}->{'processed'} ) {	    
	    # Record that this directive has been processed.
	    $node->{'directive'}->{'processed'} =  1;
	    # Determine the parameter source.
	    die("Galacticus::Build::SourceTree::Process::InputParametersValidate::Process_InputParametersValidate(): no source given")
		unless ( exists($node->{'directive'}->{'source'}) );
	    my $source = $node->{'directive'}->{'source'};
	    # Step through sibling nodes looking for input parameter directives.
	    my @objectBuilders;
	    my $sibling = $node->{'parent'}->{'firstChild'};
	    while ( $sibling ) {
		if ( $sibling->{'type'} eq "objectBuilder" ) {
		    die("Galacticus::Build::SourceTree::Process::InputParametersValidate::Process_InputParametersValidate(): no source given")
			unless ( exists($sibling->{'directive'}->{'source'}) );
		    push(@objectBuilders,$sibling->{'directive'})
			if ( $sibling->{'directive'}->{'source'} eq $source );
		}
		$sibling = $sibling->{'sibling'};
	    }
	    # Generate the variable declarations.
	    my $variableName = "allowedParameterNames_";
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
	    # Initialize new code.
	    my $code;
	    # Generate allowed multi-parameter names.
	    my $multiNames = "allowedMultiParameterNames_";
	    if ( exists($node->{'directive'}->{'multiParameters'}) ) {
		my @multiParameterNames = split(/\s*,\s*/,$node->{'directive'}->{'multiParameters'});
		unless ( &Galacticus::Build::SourceTree::Parse::Declarations::DeclarationExists($node->{'parent'},$multiNames) ) {
		    my $declaration =
		    {
			intrinsic  => "type",
			type       => "varying_string",
			attributes => [ "dimension(".scalar(@multiParameterNames).")" ],
			variables  => [ $multiNames ]
		    };
		    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},[$declaration]);
		}
		my $i = 0;
		foreach my $multiParameterName ( @multiParameterNames ) {
		    ++$i;
		    $code .= $multiNames."(".$i.")='".$multiParameterName."'\n";
		}
	    }
	    # Add module usage.
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},{moduleUse => {ISO_Varying_String => {all => 1}}});
	    # Add any extra allowed names.
	    if ( exists($node->{'directive'}->{'extraAllowedNames'}) ) {
		my @extraAllowedNames = split(" ",$node->{'directive'}->{'extraAllowedNames'});
		$code .= "allocate(".$variableName."(".scalar(@extraAllowedNames)."))\n";
		for(my $i=0;$i<scalar(@extraAllowedNames);++$i) {
		    $code .= $variableName."(".($i+1).")='".$extraAllowedNames[$i]."'\n";
		}
	    }
	    # Generate the validation code.
	    my $result;
	    if ( $node->{'parent'}->{'type'} eq "function" ) {
		if ( $node->{'parent'}->{'opener'} =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$/ ) {
		    $result = $1;
		} else {
		    $result = $node->{'parent'}->{'name'};	       
		}
	    } else {
		die('Galacticus::Build::SourceTree::Process::InputParametersValidate::Process_InputParametersValidate: parent is not a function');
	    }
	    $code .= "   call ".$result."%allowedParameters(".$variableName.",'".$source."',.false.)\n";
	    # Perform the check.
	    $code .= "   if (".$functionClassName."DsblVldtn == 0) call ".$source."%checkParameters(allowedParameterNames=".$variableName.(exists($node->{'directive'}->{'multiParameters'}) ? ",allowedMultiParameterNames=".$multiNames : "").")\n";
	    $code .= "   if (allocated(".$variableName.")) deallocate(".$variableName.")\n";
	    # Insert new code.
	    my $codeNode =
	    {
		type       => "code"           ,
		content    => $code            ,
		sibling    => undef()          ,
		parent     => $node->{'parent'},
		firstChild => undef()          ,
		source     => "Galacticus::Build::SourceTree::Process::InputParametersValidate::Process_InputParametersValidate()",
		line       => 1
	    };
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$codeNode]);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

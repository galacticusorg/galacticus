# Contains a Perl module which implements processing of input parameter directives.

package Galacticus::Build::SourceTree::Process::InputParameter;
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
use Digest::MD5 qw(md5 md5_hex md5_base64);

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'inputParameters'} = \&Process_InputParameters;

sub Process_InputParameters {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();

    # Get a list of executables (excluding test suite codes).
    my @executables;
    open(my $makeFile,$ENV{'BUILDPATH'}."/Makefile_All_Execs");
    while ( my $line = <$makeFile> ) {
	if ( $line =~ m/^all_exes\s+=\s+(.*)/ ) {
	    my $executableText = $1;
	    @executables = map {$_ =~ m/^tests\./ ? () : $_} split(" ",$executableText);
	}
    }
    close($makeFile);
    # For each executable find a list of dependent files.
    my $dependencies;
    foreach my $executableName ( @executables ) {
	(my $dependencyFileName = $executableName) =~ s/\.exe/.d/;
	if ( -e $ENV{'BUILDPATH'}."/".$dependencyFileName ) {
	    open(my $dependencyFile,$ENV{'BUILDPATH'}."/".$dependencyFileName);
	    while ( my $line = <$dependencyFile> ) {
		if ( $line =~ m/([^\/]+)\.o$/ ) {
		    $dependencies->{$executableName}->{$1} = 1;
		}
	    }
	    close($dependencyFile);
	}
    }
    # Get code directive locations.
    my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "inputParameter" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate source code for the input parameter.
	    $node->{'directive'}->{'processed'} =  1;
	    my $inputParameterSource;
	    $inputParameterSource .= "  ! Auto-generated input parameter\n";
	    if ( exists($node->{'directive'}->{'name'}) ) {
		# Simple parameter defined by a name.
		$inputParameterSource .= "  call ";
		if ( exists($node->{'directive'}->{'source'}) ) {
		    $inputParameterSource .= $node->{'directive'}->{'source'};
		} else {
		    $inputParameterSource .= "globalParameters";
		}
		my $parameterName = $node->{'directive'}->{'name'};
		$parameterName = "'".$parameterName."'" # Add delimiters to name unless the name is actually a function.
		    unless ( $parameterName =~ m/\(/ );
		$inputParameterSource .= "%value(".$parameterName.",";
		if ( exists($node->{'directive'}->{'variable'}) ) {
		    $inputParameterSource .= $node->{'directive'}->{'variable'};
		} else {
		    $inputParameterSource .= $node->{'directive'}->{'name'    };
		}
		$inputParameterSource .= ",defaultValue=".$node->{'directive'}->{'defaultValue'}
		if ( exists($node->{'directive'}->{'defaultValue'}) );
		$inputParameterSource .= ",writeOutput=".($node->{'directive'}->{'writeOutput'} eq "no" ? ".false." : ".true.")
		    if ( exists($node->{'directive'}->{'writeOutput'}) );
		$inputParameterSource .= ")\n";
	    }
	    $inputParameterSource .= "  ! End auto-generated input parameter\n\n";
	    # Create a new node.
	    my $newNode =
	    {
		type       => "code"            ,
		content    => $inputParameterSource,
		firstChild => undef(),
		source     => "Galacticus::Build::SourceTree::Process::InputParameter::Process_InputParameters()",
		line       => 1
	    };
	    # Insert the node.
	    &Galacticus::Build::SourceTree::InsertAfterNode($node,[$newNode]);
	    # Ensure input parameters module is used.
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse =>
		{
		    Input_Parameters =>
		    {
			intrinsic => 0,
			all       => 1
		    }
		}
	    };
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    # Construct list of executables which this parameter influences.
	    (my $fileName = $tree->{'name'}) =~ s/\.F90$//
		if ( $tree->{'type'} eq "file" );
	    my @influencedExecutableNames = map {exists($dependencies->{$_}->{$fileName}) ? $_ : ()} @executables
		if ( $fileName );
	    # Test for required properties.
	    foreach my $property ( "cardinality" ) {
		unless ( exists($node->{'directive'}->{$property}) ) {
		    print Dumper($node->{'directive'});
		    die("Process_InputParameters(): missing property '".$property."'");
		}
	    }
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

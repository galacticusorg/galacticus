# Contains a Perl module which implements processing of input parameter directives.

package InputParameters;
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
use XML::Simple;
use LaTeX::Encode;
require List::ExtraUtils;
require Fortran::Utils;
require Galacticus::Build::Directives;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Hooks::processHooks{'inputParameters'} = \&Process_InputParameters;

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
    my $directiveLocations = $xml->XMLin($galacticusPath."work/build/Code_Directive_Locations.xml");
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "inputParameter" && ! $node->{'directive'}->{'processed'} ) {
	    # Generate source code for the input parameter.
	    $node->{'directive'}->{'processed'} =  1;
	    my $nameForFile;
	    my $nameForDocumentation;
	    my $inputParameterSource                ;
	    $inputParameterSource .= "  ! Auto-generated input parameter\n";
	    if ( exists($node->{'directive'}->{'name'}) ) {
		# Simple parameter defined by a name.
		$inputParameterSource .= "  call ";
		if ( exists($node->{'directive'}->{'source'}) ) {
		    $inputParameterSource .= $node->{'directive'}->{'source'};
		} else {
		    $inputParameterSource .= "globalParameters";
		}
		$inputParameterSource .= "%value('".$node->{'directive'}->{'name'}."',";
		if ( exists($node->{'directive'}->{'variable'}) ) {
		    $inputParameterSource .= $node->{'directive'}->{'variable'};
		} else {
		    $inputParameterSource .= $node->{'directive'}->{'name'    };
		}
		$inputParameterSource .= ",defaultValue=".$node->{'directive'}->{'defaultValue'}
	            if ( exists($node->{'directive'}->{'defaultValue'}) );
		$inputParameterSource .= ")\n";
		# Use raw name for file and documentation.
		$nameForFile          = $node->{'directive'}->{'name'};
		$nameForDocumentation = $node->{'directive'}->{'name'};
	    } elsif ( exists($node->{'directive'}->{'iterator'})) {
		# A parameter whose name iterates over a set of possible names.
		if ( $node->{'directive'}->{'iterator'} =~ m/\(\#([a-zA-Z0-9]+)\-\>([a-zA-Z0-9]+)\)/ ) {
		    my $directiveName = $1;
		    my $attributeName = $2;
		    die('Process_InputParameter(): locations not found for directives')
			unless ( exists($directiveLocations->{$directiveName}) );
		    foreach my $fileName ( &ExtraUtils::as_array($directiveLocations->{$directiveName}->{'file'}) ) {
			foreach ( &Directives::Extract_Directives($fileName,$directiveName) ) {
			    (my $parameterName = $node->{'directive'}->{'iterator'}) =~ s/\(\#$directiveName\-\>$attributeName\)/$_->{$attributeName}/;
			    # Generate code.
			    $inputParameterSource .= "  call ";
			    if ( exists($node->{'directive'}->{'source'}) ) {
				$inputParameterSource .= $node->{'directive'}->{'source'};
			    } else {
				$inputParameterSource .= "globalParameters";
			    }
			    $inputParameterSource .= "%value('".$parameterName."',";
			    if ( exists($node->{'directive'}->{'variable'}) ) {
				(my $variableName = $node->{'directive'}->{'variable'}) =~ s/\$1/$_->{$attributeName}/;
				$inputParameterSource .= $variableName;
			    } else {
				$inputParameterSource .= $parameterName;
			    }
			    if ( exists($node->{'directive'}->{'defaultValue'}) ) {
				(my $defaultValue = $node->{'directive'}->{'defaultValue'}) =~ s/\$1/$_->{$attributeName}/;
				$inputParameterSource .= ",defaultValue=".$defaultValue;
			    }
			    $inputParameterSource .= ")\n";
			}
		    }
		} else {
		    die('Process_InputParameter(): nothing to iterate over');
		}
		# Construct names for file and documentation.
		$nameForFile          = $node->{'directive'}->{'iterator'};
		$nameForDocumentation = $node->{'directive'}->{'iterator'};
	    }	    
	    $inputParameterSource .= "  ! End auto-generated input parameter\n\n";
	    # Create a new node.
	    my $newNode =
	    {
		type       => "code"            ,
		content    => $inputParameterSource,
		firstChild => undef()
	    };
	    # Insert the node.
	    &SourceTree::InsertAfterNode($node,[$newNode]);
	    # Ensure input parameters module is used.
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse =>
		{
		    Input_Parameters2 =>
		    {
			intrinsic => 0,
			all       => 1
		    }
		}
	    };
	    &ModuleUses::AddUses($node->{'parent'},$usesNode);
	    # Construct list of executables which this parameter influences.
	    (my $fileName = $tree->{'name'}) =~ s/\.F90$//
		if ( $tree->{'type'} eq "file" );
	    my @influencedExecutableNames = map {exists($dependencies->{$_}->{$fileName}) ? $_ : ()} @executables
		if ( $fileName );
	    # Create documentation.
	    system("mkdir -p doc/inputParameters");
	    open(my $defHndl,">doc/inputParameters/".$nameForFile.".tex");
	    print $defHndl "\\noindent {\\normalfont \\bfseries Name:} {\\normalfont \\ttfamily ".latex_encode($nameForDocumentation)."}\\\\\n";
	    my $definedIn;
	    my @hyperTarget;
	    my $fileIn;
	    my $parent = $node->{'parent'};
	    while ( $parent ) {
		if ( exists($parent->{'name'}) ) {
		    $definedIn = $parent->{'type'}.":".latex_encode($parent->{'name'})
			if (
			    ! $definedIn
			    &&
			    (
			     $parent->{'type'} eq "module" 
			     ||
			     $parent->{'type'} eq "subroutine"
			     ||
			     $parent->{'type'} eq "function"
			     ||
			     $parent->{'type'} eq "program"
			    )
			);
		    $fileIn = $parent->{'name'}
		        if ( $parent->{'type'} eq "file" );
		    if ( $definedIn ) {
			my $hyperTargetName = $parent->{'name'};
			$hyperTargetName = lc($hyperTargetName)
			    unless ( $parent->{'type'} eq "file" );
			unshift(@hyperTarget,$hyperTargetName);
		    }
		}
		$parent = $parent->{'parent'};
	    }
	    print $defHndl "{\\normalfont \\bfseries Type:} ".$node->{'directive'}->{'type'}." ".$node->{'directive'}->{'cardinality'}."\\\\\n";	    
	    if ( exists($node->{'directive'}->{'defaultValue'}) ) {
		my $defaultValue = $node->{'directive'}->{'defaultValue'};
		if ( $defaultValue =~ m/^\s*([\d\.]+)d(\+|\-)?(\d+)\s*$/ ) {
		    my $number     = $1;
		    my $sign       = $2;
		    my $exponent   = $3;
		    $defaultValue  = "\$".$number;
		    if ( $exponent != 0 ) {
			$defaultValue .= " \\times 10^{";
			$defaultValue .= "-"
			    if ( $sign && $sign eq "-" );
			$defaultValue .= $exponent."}";
		    }
		    $defaultValue   .= "\$";
		}
		print $defHndl "{\\normalfont \\bfseries Default value:} ".$defaultValue;
		print $defHndl " ".$node->{'directive'}->{'defaultSource'}
		    if ( exists($node->{'directive'}->{'defaultSource'}) );
		print $defHndl "\\\\\n";
	    }
	    print $defHndl "{\\normalfont \\bfseries Description:} ".$node->{'directive'}->{'description'}." \\\\\n";	    
	    print $defHndl "{\\normalfont \\bfseries Defined in:} \\hyperlink{".join(":",@hyperTarget)."}{{\\normalfont \\ttfamily ".$definedIn."}}\\\\\n";
	    print $defHndl "{\\normalfont \\bfseries File:} \\hyperlink{".$fileIn."}{{\\normalfont \\ttfamily ".latex_encode($fileIn)."}}\\\\\n";
	    print $defHndl "{\\normalfont \\bfseries Used by:} ".join(", ",map {"\\hyperlink{".&replace($_,qr/\.exe$/s,".F90")."}{\\normalfont \\ttfamily ".latex_encode($_)."}"} @influencedExecutableNames)."\\\\\n\n";
	    close($defHndl);
	}
	$node = &SourceTree::Walk_Tree($node,\$depth);
    }
}

sub replace {
    my $text    = shift();
    my $regEx   = shift();
    my $replace = shift();
    $text =~ s/$regEx/$replace/;
    return $text;
}

1;

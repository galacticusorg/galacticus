#!/usr/bin/env perl
use lib './perl';
use XML::Simple;
use Data::Dumper;
use Switch;
use Sort::Topological qw(toposort);
use Fortran::Utils;
use UNIVERSAL 'isa';

# Scans source code for "!#" directives and generates an include file.

# Define the source directory.
if ( $#ARGV != 1 ) {die "Usage: Build_Include_File.pl sourcedir xmlFile"};
$sourcedir = $ARGV[0];
$sourcedirs[0] = $sourcedir."/source";
$xmlFile = $ARGV[1];

# Specify verbosity.
$verbosity = 0;

# Specify evolvability of different object types.
my %isEvolvable = (
    scalar      => 1,
    array       => 1,
    history     => 1,
    keplerOrbit => 0
    );

# Load the file of directive locations.
if ( -e "./work/build/Code_Directive_Locations.xml" ) {
    $xml = new XML::Simple;
    $locations = $xml->XMLin("./work/build/Code_Directive_Locations.xml");
    $gotLocations = 1;
} else {
    $gotLocations = 0;
}

# Create XML object and process the XML file.
$xml = new XML::Simple;
$instructions = $xml->XMLin($xmlFile);

# Initially not inside any module.
$moduleName = "";

# Find files to scan.
if ( $gotLocations == 1 ) {
    if ( isa($locations->{$instructions->{'directive'}}->{'file'},'ARRAY') ) {
	@filesToScan    = @{$locations->{$instructions->{'directive'}}->{'file'}};
    } else {

	if (exists($locations->{$instructions->{'directive'}}->{'file'})) {
	    $filesToScan[0] =   $locations->{$instructions->{'directive'}}->{'file'} ;
	}
    }
} else {
    foreach $srcdir ( @sourcedirs ) {
	opendir(indir,$srcdir) or die "Can't open the source directory: #!";
	while ($fname = readdir indir) {	
	    if ( $fname =~ m/\.[fF](90)??t??$/ && $fname !~ m/^\.\#/ ) {
		$fullname = $srcdir."/".$fname;
		push(@filesToScan,$fullname);
	    }
	}
    }
}

# Scan all necessary files.
foreach $fullname ( @filesToScan ) {
    # Determine source code type.
    $codeType = "fortran";
    $codeType = "c" if ( $fullname =~ m/\.c(pp)??$/ );

    open($infile,$fullname) or die "Can't open input file: #!";
    until ( eof($infile) ) {

	switch ( $codeType) {
	    case ( "fortran" ) {
		# Get next line from the Fortran source.
		&Fortran_Utils::Get_Fortran_Line($infile,$rawLine,$processedLine,$bufferedComments);
	    }
	    case ( "c" ) {
		# Get the next line from a C(++) source.
		$processedLine = <$infile>;
		$rawLine = $processedLine;	
	    }
	}
		
	# Check if we've entered or left a module.
	if ( $processedLine =~ /^\s*module\s*([a-z0-9_]+)\s*$/i ) {$moduleName = $1};
	if ( $processedLine =~ /^\s*end\s+module\s*([a-z0-9_]+)\s*$/i ) {$moduleName = ""};
		
	if ( $rawLine =~ m/^\s*(!|\/\/)\#\s+(<\s*([a-zA-Z]+)+.*>)\s*$/ ) {
	    $xmlCode = $2."\n";
	    $xmlTag  = $3;
	    # Read ahead until a matching close tag is found.
	    unless ( $xmlCode =~  m/\/>/ ) {
		$nextLine = "";
		until ( $nextLine =~ m/<\/$xmlTag>/ || eof($infile) ) {
		    switch ( $codeType) {
			case ( "fortran" ) {
			    &Fortran_Utils::Get_Fortran_Line($infile,$nextLine,$processedLine,$bufferedComments);
			}
			case ( "c" ) {
			    $nextLine = <$infile>;
			}
		    }
		    $nextLine =~ s/^\s*(!|\/\/)\#\s+//;
		    $xmlCode .= $nextLine;
		}
	    }
	    # Check if this dircetive matches that which we are currently processing.
	    if ( $xmlTag eq $instructions->{'directive'} ) {
		$data = $xml->XMLin($xmlCode);
		if ( $verbosity == 1 ) {
		    print Dumper($data);
		}
		# Process the directive based on the type of instruction we have.
		switch ( $instructions->{'type'} ) {
		    # Instruction is to insert some type of code.				
		    case ( "code" ) {
			switch ( $instructions->{'action'} ) {
			    # Action is to assign a procedure pointer.
			    case ( "procPointer" ) {
				$inserts{$data->{'unitName'}} = $instructions->{'pointerName'}." => ".$data->{'unitName'}."\n";
				if ( exists($instructions->{'pointerAction'}) ) {
				    $inserts{$data->{'unitName'}} .= $instructions->{'pointerAction'}."\n";
				}
			    }
			    # Action is to call a subroutine.
			    case ( "subroutine" ) {
				$inserts{$data->{'unitName'}} = "call ".$data->{'unitName'}."(";
				# Check for existance of an option name. This will be used as a parameter to select between
				# different implementations of a component.
				if ( exists($data->{'optionName'}) ) {
				    if ( exists($data->{'optionName'}->{'content'}) ) {
					$content = $data->{'optionName'}->{'content'};
				    } else {
					$content = $data->{'optionName'};
				    }
				    $inserts{$data->{'unitName'}} .= $content;
				}
				# If we have subroutine arguments, append them to the call.
				if ( exists($instructions->{'subroutineArgs'}) ) {
				    if ( exists($instructions->{'subroutineArgs'}->{$codeType}) ) {
					$arguments = $instructions->{'subroutineArgs'}->{$codeType};
				    } else {
					die ("Build_Include_File.pl: subroutine arguments not derived for this language") unless ( $codeType eq "fortran" );
					$arguments = $instructions->{'subroutineArgs'};
				    }
				} else {
				    $arguments = "";
				}
				unless ( $arguments eq "" ) {
				    if ( $inserts{$data->{'unitName'}} !~ m/\($/ ) {$inserts{$data->{'unitName'}} .= ","};
				    $arguments =~ s/\#label/$data->{'label'}/g;
				    $inserts{$data->{'unitName'}} .= $arguments;
				}
				$inserts{$data->{'unitName'}} .= ")\n";
				# If some action is specified, perform this action after the subroutine call.
				if ( exists($instructions->{'subroutineAction'}) ) {
				    $inserts{$data->{'unitName'}} .= $instructions->{'subroutineAction'}."\n";
				}
			    }
			}
		    }
		    # Instruction is to insert a "use" statement for the module which the directive appears in.
		    case ( "moduleUse" ) {
			switch ( $codeType ) {
			    case ( "fortran" ) {
				$inserts{$moduleName} = "use ".$moduleName."\n";	
			    }
			    case ( "c" ) {
				($leafName = $fullname) =~ s/.*?([^\/]+)\.c(pp)??$/\1.o/;
				$inserts{$moduleName} = "!: ./work/build/".$leafName."\n";
			    }
			}
		    }
		    # Instruction is to declare C bindings as appropriate.
		    case ( "cBinding" ) {
			switch ( $codeType ) {
			    case ( "c" ) {
				$inserts{$moduleName} = "bind(c,name='".$data->{'unitName'}."') :: ".$data->{'unitName'}."\n";
			    }
			}
		    }
		    # Instruction is to process a new component property method.
		    case ( "methods" ) {
			# Get the type of this method.
			$inserts{$data->{'methodName'}} = &Get_Type($data);
		    }
		    # Instruction is to generate a C-bound wrapper for a component property method.
		    case ( "cWrapper" ) {
			# Get the type of this method.
			$inserts{$data->{'methodName'}} = &Get_Type($data);
		    }
		    # Instruction is to generate templates for C wrappers for a component property method.
		    case ( "cTemplate" ) {
			# Get the type of this method.
			$inserts{$data->{'methodName'}} = &Get_Type($data);
		    }
		    # Instruction is to process a pipe to/from a component property.
		    case ( "pipe" ) {
			# Get the type of this pipe.
			$inserts{$data->{'pipeName'}} = &Get_Type($data);
		    }
		    # Instruction is to build a set of calls to evaluate property derivatives.
		    case ( "derivatives" ) {
			$inserts{$data->{'methodName'}} = &Get_Type($data);
		    }
		    # Instruction is to build a set of calls that initialize methods.
		    case ( "initializeMethods" ) {
			# Get the type of this method.
			$inserts{$data->{'methodName'}} = &Get_Type($data);
		    }
		    # Process labels.
		    case ( "label" ) {
			# Get the prefix for this method.
			++$labelCount{$instructions->{'prefix'}};
			$inserts{$data->{'unitName'}} .= "integer, public, parameter :: ".$instructions->{'prefix'}.$data->{'label'}."=".$labelCount{$instructions->{'prefix'}}."\n";
		    }
		    # Process option names associated with component implementations.
		    case ( [ "optionNames", "optionDefinitions" ] ) {
			if ( exists($data->{'optionName'}) ) {
			    if ( exists($data->{'optionName'}->{'content'}) ) {
				$content = $data->{'optionName'}->{'content'};
			    } else {
				$content = $data->{'optionName'};
			    }
			    # Check if a default option is specified
			    if ( exists($data->{'optionName'}->{'default'}) ) {
				# A default has already been specified, report and error.
				if ( exists($inserts{$content}) && $inserts{$content} ne "none" ) {die("Build_Include_File.pl: FATAL - multiple defaults found")};
				# Store this default.
				$inserts{$content} = $data->{'optionName'}->{'default'};
			    } else {
				# If no default has already been specified, then set to "none".
				unless ( exists($inserts{$content}) ) {
				    $inserts{$content} = "none";
				}
			    }
			}
		    }
		}
		# Check for "after" tags, which influence the sort order of commands.
		if ( exists($data->{'after'}) ) {
		    # Store the names of procedures named in after tags.
		    if ( UNIVERSAL::isa($data->{'after'},"ARRAY" ) ) {
			foreach $after ( @{$data->{'after'}} ) {
			    ${$Dependencies{$after}}[++$#{$Dependencies{$after}}] = $data->{'unitName'};
			}
		    } else {
			${$Dependencies{$data->{'after'}}}[++$#{$Dependencies{$data->{'after'}}}] = $data->{'unitName'};
		    }
		}
		# Check for "before" tags, which influence the sort order of commands.
		if ( exists($data->{'before'}) ) {
		    # Store the names of procedures named in before tags.
		    if ( UNIVERSAL::isa($data->{'before'},"ARRAY" ) ) {
			push(@{$Dependencies{$data->{'unitName'}}},@{$data->{'before'}});
		    } else {
			${$Dependencies{$data->{'unitName'}}}[++$#{$Dependencies{$data->{'unitName'}}}] = $data->{'before'};
		    }
		}
		# Check if a specific sorting name has been given for this directive.
		if ( exists($data->{'sortName'}) ) {
		    # One has, so store it.
		    $sortNames{$data->{'unitName'}} = $data->{'sortName'};
		}
	    }
	}
    }
    close($infile);
}

# Sort the list of commands for dependencies (from <before> and <after> tags).
sub Dependencies { @{$Dependencies{$_[0]} || []}; } # Comparison function for dependency sort.
# Create a list of names for sorting.
foreach $BlockName ( keys(%inserts) ) {
    # Sort on the given sorting name or on the block name otherwise.
    if ( exists($sortName{$BlockName}) ) {
	$UnsortedBlocks[++$#UnsortedBlocks] = $sortName{$BlockName};
    } else {
	$UnsortedBlocks[++$#UnsortedBlocks] = $BlockName;
    }
}
# Perform initial alphanumerical sort.
@PresortedBlocks = sort(@UnsortedBlocks);
# Perform dependency sort.
@SortedBlocks    = toposort(\&Dependencies, \@PresortedBlocks);

# Open the output file.
open(includeFile,">".$instructions->{'fileName'});
switch ( $instructions->{'type'} ) {
    # For option defintions, output code to define the variables.
    case ( "optionDefinitions" ) {
	# Create variables.
	foreach $optionName ( @SortedBlocks ) {
	    print includeFile "type(varying_string) :: $optionName\n";
	}
    }
    # For option names, output code to read in the values from the parameter file.
    case ( "optionNames" ) {
	# Read in parameter value.
	foreach $optionName ( @SortedBlocks ) {
	    print includeFile "call Get_Input_Parameter('".$optionName."',".$optionName;
	    unless ( $inserts{$insert} eq "none" ) {
		print includeFile ",defaultValue='".$inserts{$optionName}."'";
	    }
	    print includeFile ")\n";
	}
    }
    # For methods, create code which creates procedure pointers for the given method.
    case ( "methods" ) {
	foreach $method ( @SortedBlocks ) {
	    $suffix = &Get_Suffix($inserts{$method});
	    print includeFile "procedure(Get_Template".$suffix."), pointer, public :: $method => null()\n";
	    print includeFile "procedure(Set_Template".$suffix."), pointer, public :: $method\_Set => null()\n";
	    print includeFile "procedure(Rate_Adjust_Template".$suffix."), pointer, public :: $method\_Rate_Adjust => null()\n"
		unless ( $isEvolvable{$inserts{$method}} == 0 );
	    print includeFile "procedure(Rate_Compute_Template), pointer, public :: $method\_Rate_Compute => null()\n"
		unless ( $isEvolvable{$inserts{$method}} == 0 );
	}
    }
    # For cWrapper, create C-bound wrappers for the given method.
    case ( "cWrapper" ) {
	foreach $method ( @SortedBlocks ) {
	    $suffix = &Get_Suffix($inserts{$method});
	    switch ( $suffix ) {
		case ( "" ) {
		    # Simple scalar method.
		    print includeFile "function c".$method."(cNode) bind(c,name='".$method."')\n";
		    print includeFile "  use, intrinsic :: ISO_C_Binding\n";
		    print includeFile "  use Tree_Nodes\n";
		    print includeFile "  implicit none\n";
		    print includeFile "  real(c_double)          :: c".$method."\n";
		    print includeFile "  type(c_ptr),    value   :: cNode\n";
		    print includeFile "  type(treeNode), pointer :: thisNode\n";
		    print includeFile "  call c_f_pointer(cNode,thisNode)\n";
		    print includeFile "  c".$method."=".$method."(thisNode)\n";
		    print includeFile "  return\n";
		    print includeFile "end function c".$method."\n\n";
		}
	    }
	}
    }
    # For cTemplate, create C templates for the given method.
    case ( "cTemplate" ) {
	print includeFile "extern \"C\"\n"; 
	print includeFile "{\n";
	foreach $method ( @SortedBlocks ) {
	    $suffix = &Get_Suffix($inserts{$method});
	    switch ( $suffix ) {
		case ( "" ) {
		    # Simple scalar method.
		    print includeFile "  double ".$method."(void *thisNode);\n";
		}
	    }
	}
	print includeFile "}\n";
    }
    # Create code to compute derivatives for properties.
    case ( "derivatives" ) {
	foreach $method ( @SortedBlocks ) {
	    print includeFile "if (.not.interrupt) call $method\_Rate_Compute(thisNode,interrupt,interruptProcedure)\n"
		unless ( $isEvolvable{$inserts{$method}} == 0 );
	}
    }
    # Create code to intialize all methods.
    case ( "initializeMethods" ) {
	foreach $method ( @SortedBlocks ) {
	    $suffix = &Get_Suffix($inserts{$method});
	    print includeFile "$method\_Rate_Adjust  => Tree_Node_Rate_Adjust".$suffix."_Dummy\n"
		unless ( $isEvolvable{$inserts{$method}} == 0 );
	    print includeFile "$method\_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy\n"
		unless ( $isEvolvable{$inserts{$method}} == 0 );
	}
    }
    # Create code to define procedure pointers for pipes.
    case ( "pipe" ) {
	foreach $pipe ( @SortedBlocks ) {
	    $suffix = &Get_Suffix($inserts{$pipe});
	    print includeFile "procedure(Rate_Adjust_Template".$suffix."),  pointer, public :: $pipe  => null()\n";
	}
    }
    # Generic process: simply output the contents of the insert hash.
    else {
	foreach $insert ( @SortedBlocks ) {
	    print includeFile $inserts{$insert};
	}
    }
}
close(includeFile);
exit;

sub Get_Type {
    # Returns the type of a method of pipe.
    my $data = shift;
    # Assume scalar type by default
    my $type = "scalar";
    # If a type is specified, then return it instead.
    if ( exists($data->{'type'}) ) {$type = $data->{'type'}};
    return $type;
}

sub Get_Suffix {
    # Returns the suffix for a method.
    my $methodType = shift;
    # Determine the suffix to use.
    my $suffix = "";
    switch ( $methodType ) {
	case ( "scalar" ) {
	    $suffix = "";
	}
	case ( "array" ) {
	    $suffix = "_Array";
	}
	case ( "history" ) {
	    $suffix = "_History";
	}
	case ( "keplerOrbit" ) {
	    $suffix = "_Kepler_Orbit";
	}
	else {
	    die("Build_Include_File.pl: unrecognized method type");
	}
    }
    return $suffix;
}

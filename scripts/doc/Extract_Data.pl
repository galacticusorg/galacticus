#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use lib './perl';
use strict;
use warnings;
use Switch;
use XML::Simple;
use Data::Dumper;
use LaTeX::Encode;
use Fcntl qw(SEEK_SET);
use UNIVERSAL;
require Fortran::Utils;
require Galacticus::Doc::Parameters;

# Scan Fortran90 source code and extract various useful data from "!@" lines.
# Andrew Benson 12-Mar-2010

# Get source directory
die 'Usage: Extract_Data.pl <sourceDir> <outputRoot>'
    unless ( scalar(@ARGV) == 2 );
my $sourceDir  = $ARGV[0];
my $outputRoot = $ARGV[1];

# Specify regex for derived-type declarations.
my $derivedTypeOpenRegex  = "^\\s*type\\s*(,\\s*abstract\\s*|,\\s*public\\s*|,\\s*private\\s*|,\\s*extends\\s*\\([a-zA-Z0-9_]+\\)\\s*)*(::)??\\s*([a-z0-9_]+)\\s*\$";
my $derivedTypeCloseRegex = "^\\s*end\\s+type\\s+([a-z0-9_]+)";

# Specify regex for type-bound procedures.
my $typeBoundRegex = "^\\s*(procedure|generic)\\s*(,\\s*nopass\\s*)*::\\s*([a-z0-9_]+)\\s*=>\\s*([a-z0-9_,\\s]+)\$";

# Create XML object.
my $xml = new XML::Simple;

# Open the source directory.
opendir(dirHndl,$sourceDir);
# Read files from the source directory.
my @fileNames = ( "work/build/objects.nodes.components.Inc" );
while ( my $fileName = readdir(dirHndl) ) {
    # Find Fortran 90 and C++ source files.
    push(@fileNames,$sourceDir."/".$fileName)
	if ( ( $fileName =~ m/\.F90$/ || $fileName =~ m/\.cpp$/ ) && $fileName !~ m/^\.\#/ );
}
closedir(dirHndl);

# Initialize data hashes.
my %parametersRead;
my %parametersData;
my %objects;
my %enumerations;
my @emptyDefaults;

foreach my $fileName ( @fileNames ) {
    # Get a printable file name.
    (my $fileNamePrint = $fileName) =~ s/_/\\_/g;
    (my $leafName = $fileName) =~ s/^$sourceDir\///;
    (my $leafNamePrint = $leafName) =~ s/_/\\_/g;

    # Add the file to the list of filenames to process.
    my @fileNames     = ( $fileName );
    my @filePositions = (        -1 );

    # Initialize the data directive buffer.
    my $dataDirective  = "";
    my $openingElement = "";

    # Initialize the program unit stack.
    my @programUnits;
    undef(@programUnits);

    # Initialize derived-type boolean.
    my $inDerivedType = "";

    # Process files until none remain.
    while ( scalar(@fileNames) > 0 ) {
	# Open the file.
	open(my $fileHandle,$fileNames[0]);
	seek($fileHandle,$filePositions[0],SEEK_SET) unless ( $filePositions[0] == -1 );
	
	# Scan the file.
	until ( eof($fileHandle) ) {
	    
	    my $rawLine;
	    my $processedLine;
	    my $bufferedComments;
	    if ( $fileName =~ m/\.(F90|Inc)$/ ) {
		&Fortran_Utils::Get_Fortran_Line($fileHandle,$rawLine,$processedLine,$bufferedComments);
	    } else {
		$rawLine = <$fileHandle>;
	    }
	    
	    # Detect include files, and recurse into them.
	    if ( defined($processedLine) && $processedLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
		my $includeFile = $1;
		if ( -e $sourceDir."/../work/build/".$includeFile ) {
		    $filePositions[0] = tell($fileHandle);
		    unshift(@fileNames,$sourceDir."/../work/build/".$includeFile);
		    unshift(@filePositions,-1);
		    last;
		} elsif ( -e $sourceDir."/".$includeFile ) {
		    $filePositions[0] = tell($fileHandle);
		    unshift(@fileNames,$sourceDir."/".$includeFile);
		    unshift(@filePositions,-1);
		    last;
		}
	    }
	    
	    # Check for entering, leaving program units.
	    if ( $rawLine =~ m/^\s*module\s+([a-zA-Z0-9_]+)/ ) {
		my $moduleName = $1;
		unless ( $moduleName eq "procedure" ) {push(@programUnits,"module:".$moduleName)};
	    }
	    if ( $rawLine =~ m/^\s*program\s+([a-zA-Z0-9_]+)/ ) {push(@programUnits,"program:".$1)};
	    if ( $rawLine =~ m/^\s*(pure\s+|elemental\s+|recursive\s+)*\s*subroutine\s+([a-zA-Z0-9_]+)/ ) {push(@programUnits,"subroutine:".$2)};
	    if ( $rawLine =~ m/^\s*(pure\s+|elemental\s+|recursive\s+)*\s*(real|integer|double precision|character|logical)*\s*(\((kind|len)=[\w\d]*\))*\s*function\s+([a-zA-Z0-9_]+)/ ) {push(@programUnits,"function:".$5)};
	    if ( $rawLine =~ m/^\s*end\s+(program|module|subroutine|function)\s/ ) {pop(@programUnits)};

	    # Check for derived-type definitions.
	    if ( defined($processedLine) ) {
		if ( $processedLine =~ m/$derivedTypeOpenRegex/i ) {
		    $inDerivedType = $3;
		    $objects{lc($inDerivedType)}->{'name'} = $inDerivedType;
		    if ( defined($1) ) {
			my $attributes = $1;
			if ( $attributes =~ m/extends\s*\(\s*([a-zA-Z0-9_]+)\s*\)/ ) {
			    $objects{lc($inDerivedType)}->{'extends'} = lc($1);
			}
		    }
		}
		if ( $processedLine =~ m/$derivedTypeCloseRegex/i ) {
		    $inDerivedType = "";
		}
	    }
	    
	    # Check for type-bound procedures.
	    if ( $inDerivedType ne "" ) {
		if ( $processedLine =~ m/$typeBoundRegex/i ) {
		    my $method = lc($3);
		    $objects{lc($inDerivedType)}->{'methods'}->{$method}->{'description'} = "UNDEFINED"
			unless ( exists($objects{lc($inDerivedType)}->{'methods'}->{$method}) );
		}
	    }
	    
	    # Check for parameter reading lines - used to check that all parameters are documented.
	    if ( $rawLine =~ m/^\s*call\s+get_input_parameter\s*\(\s*[\"\']\s*(\w+)\s*[\"\']/i ) {
		$parametersRead{$1} .= $fileName." ";
	    }
	    
	    # Search for "!@".
	    my $process = 0;
	    if ( $rawLine =~ m/^\s*(\!|\/\/)\@\s/ ) {
		# Found directive - add it to the buffer.
		if ( $dataDirective eq "" && $rawLine =~ m/<([^>]+)>/ ) {
		    $openingElement = $1;
		} elsif ( $rawLine =~ m/<\/$openingElement>/ ) {
		    $process = 1;
		}
		$rawLine =~ s/^\s*(\!|\/\/)\@\s*//;
		$rawLine =~ s/\s*$//;
		$dataDirective .= $rawLine." ";
	    } elsif ( $dataDirective ne "" ) {
		$process = 1;
	    }
	    if ( $process == 1 ) {
		# Process a directive.
		$openingElement = "";
		my $data = eval{$xml->XMLin($dataDirective, KeepRoot => 1)};
		my $lineNumber = $.;
		die("Extract_Data.pl failed in file ".$fileName." at line ".$lineNumber." with message:\n".$@."and data \n".$dataDirective) if ($@);
		
		# Loop over all data types found.
		foreach my $dataType ( keys(%{$data}) ) {
		    my $contents = $data->{$dataType};
		    switch ( $dataType ) {
			case ( "enumeration" ) {
			    foreach ( @{$contents->{'entry'}}) {
				push(@{$enumerations{$contents->{'name'}}->{'entry'}},$_->{'label'});
			    }
			    my $programUnitIndex = scalar(@programUnits)-1;
			    $programUnitIndex = -1 
				unless ( defined($programUnits[$programUnitIndex]) );
			    my $regEx = "^module";
			    until ( $programUnitIndex == -1 || $programUnits[$programUnitIndex] =~ m/$regEx/ ) {
				--$programUnitIndex
			    }
			    ($enumerations{$contents->{'name'}}->{'module'     } = $programUnits[$programUnitIndex]) =~ s/^module://;
			    $enumerations {$contents->{'name'}}->{'file'       } = $leafName;
			    $enumerations {$contents->{'name'}}->{'description'} = $contents->{'description'};
			}
			case ( "inputParameter" ) {
			    # Handle parameters defined as regular expressions.
			    if ( exists($contents->{'regEx'}) ) {
				$contents->{'name'} = "[regEx] ".latex_encode(&Parameters::ExpandRegEx($contents->{'regEx'},$sourceDir));
			    }
			    # Construct output data for this parameter.
			    (my $printName = $contents->{'name'}) =~ s/([^\\])_/$1\\_/g;
			    my $buffer  = "\\noindent {\\bf Name:} {\\tt ".$printName."}\\\\\n";
			    $buffer .= "{\\bf Attached to:} ";
			    my $attachedTo;
			    my $attachedAt;
			    # Detect empty "defaultValue" elements.
			    if ( exists($contents->{'defaultValue'}) && UNIVERSAL::isa($contents->{'defaultValue'},"HASH") ) {
				push(
				    @emptyDefaults,
				    {
					name => $contents->{'name'},
					file => $fileNames[0]
				    }
				    );
				delete($contents->{'defaultValue'});
			    }
			    # Determine to what the parameter is attached.
			    if ( exists($contents->{'attachedTo'}) ) {
				my $programUnitIndex = scalar(@programUnits)-1;
				my $regEx = $contents->{'attachedTo'}.":";
				$programUnitIndex = -1 
				    unless ( defined($programUnits[$programUnitIndex]) );
				until ( $programUnitIndex == -1 || $programUnits[$programUnitIndex] =~ m/$regEx/ ) {
				    --$programUnitIndex
				}
				if ( $programUnitIndex >= 0 ) {
				    $attachedTo = "{\\tt ".$programUnits[$programUnitIndex]."}";
				    $attachedAt = $programUnitIndex;
				} else {
				    $attachedTo = "unknown";
				    $attachedAt = -1;
				}
			    } else {
				$attachedTo = "{\\tt ".$programUnits[-1]."}";
				$attachedAt = scalar(@programUnits)-1;
			    }
			    my $attachedLink = $leafName.":";
			    for(my $iAttach=0;$iAttach<=$attachedAt;++$iAttach) {
				if ( $programUnits[$iAttach] =~ m/^[^:]+:(.*)/ ) {
				    my $unitName = lc($1);
				    $unitName =~ s/\\_/_/g;
				    $attachedLink .= $unitName.":";
				}
			    }
			    $attachedLink =~ s/:$//;
			    my $targetPrefix;
			    my $targetSuffix;
			    if ( $attachedTo =~ m/:/ ) {
				$targetPrefix = "\\hyperlink{".$attachedLink."}{";
				$targetPrefix =~ s/\{\\tt\s+([^\}]+)\}/$1/;
				$targetPrefix =~ s/program:/prog:/;
				$targetPrefix =~ s/module:/mod:/;
				$targetPrefix =~ s/subroutine:/sub:/;
				$targetPrefix =~ s/function:/func:/;
				$targetSuffix = "}";
			    } else {
				$targetPrefix = "";
				$targetSuffix = "";
			    }
			    $attachedTo =~ s/_/\\_/g;
			    $buffer .= $targetPrefix.$attachedTo.$targetSuffix."\\\\\n";
			    $buffer .= "{\\bf File:} \\hyperlink{".$leafName."}{{\\tt ".$leafNamePrint."}}\\\\\n";
			    $buffer .= "{\\bf Default value:} ";
			    if ( exists($contents->{'defaultValue'}) ) {
				$buffer .= $contents->{'defaultValue'};
			    } else {
				$buffer .= "none";
			    }
			    $buffer .= "\\\\\n";
			    $buffer .= "{\\bf Description:} ".$contents->{'description'};
			    $buffer .= "\\\\" unless ( $buffer =~ m/\}\s*$/ );
			    $buffer .= "\n\n";
			    my $descriptor = 
			    {
				description      => $buffer      ,
				file             => $fileNames[0],
				attachmentStatus => $attachedTo 
			    };
			    # Use this instance of this parameter description if it is the first
			    # instance seen, or if previous instances had unknown attachment status.
			    my $useThisInstance = 1;
			    $useThisInstance = 0
				if ( exists($parametersData{$contents->{'name'}}) && $parametersData{$contents->{'name'}}->{'attachmentStatus'} ne "unknown" );
			    $parametersData{$contents->{'name'}} = $descriptor
				if ( $useThisInstance == 1 );
			}
			case ( "objectMethod" ) {
			    my $object = lc($contents->{'object'});
			    $objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'method'     } = $contents->{'method'     };
			    $objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'description'} = $contents->{'description'};
			    $objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'arguments'  } = $contents->{'arguments'  }
			    if ( exists($contents->{'arguments'}) );
			    $objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'type'       } = $contents->{'type'       }
			    if ( exists($contents->{'type'}) );
			}
			case ( "objectMethods" ) {
			    my $object = lc($contents->{'object'});
			    my @methods;
			    if ( UNIVERSAL::isa($contents->{'objectMethod'},"ARRAY") ) {
				@methods = @{$contents->{'objectMethod'}};
			    } else {
				push(@methods,$contents->{'objectMethod'});
			    }
			    foreach my $method ( @methods ) {
				$objects{$object}->{'methods'}->{lc($method->{'method'})}->{'method'     } = $method->{'method'     };
				$objects{$object}->{'methods'}->{lc($method->{'method'})}->{'description'} = $method->{'description'};
				$objects{$object}->{'methods'}->{lc($method->{'method'})}->{'arguments'  } = $method->{'arguments'  }
				if ( exists($method->{'arguments'}) );
				$objects{$object}->{'methods'}->{lc($method->{'method'})}->{'type'       } = $method->{'type'       }
				if ( exists($method->{'type'}) );
			    }
			}
		    }
		}
		
		# Reset the directive buffer.
		$dataDirective = "";
	    }
	    
	}
	
	# Close the file and shift the list of filenames.
	if ( eof($fileHandle) ) {
	    shift(@fileNames);
	    shift(@filePositions);
	}
	close($fileHandle);
    }
}

# Copy any methods inherited from parent classes.
my $foundExtensions = 1;
while ( $foundExtensions == 1 ) {
    $foundExtensions = 0;
    foreach my $object ( sort(keys(%objects)) ) {
	if ( defined($objects{$object}->{'extends'}) ) {
	    $foundExtensions = 1;
	    my $parent = $objects{$object}->{'extends'};
	    unless ( defined($objects{$parent}->{'extends'}) ) {
		foreach my $method ( keys(%{$objects{$parent}->{'methods'}}) ) {
		    if ( ( ! exists($objects{$object}->{'methods'}->{$method}) ) || $objects{$object}->{'methods'}->{$method}->{'description'} eq "UNDEFINED" ) {
			$objects{$object}->{'methods'}->{$method} = $objects{$parent}->{'methods'}->{$method};
		    }
		}
		undef($objects{$object}->{'extends'});
	    }
	}
    }
}

# Open output files.
open(parametersHndl,">".$outputRoot."Parameters.tex");
open(methodHndl    ,">".$outputRoot."Methods.tex");

# Write parameter descriptions.
my @sortedParameters = sort(keys(%parametersData));
foreach my $parameterName ( @sortedParameters ) {
    print parametersHndl $parametersData{$parameterName}->{'description'};
}

# Write method descriptions.
foreach my $object ( sort(keys(%objects)) ) {
    my $hasEntries = 0;
    foreach my $method ( sort(keys(%{$objects{$object}->{'methods'}})) ) {
	if ( $objects{$object}->{'methods'}->{$method}->{'description'} eq "UNDEFINED" ) {
	    print "Warning: missing data for method ".$method." of ".$object." object\n";
	} else {
	    $hasEntries = 1;
	}
    }
    if ( $hasEntries == 1 ) {
	print methodHndl "\\subsubsection{\\large {\\tt ".$objects{$object}->{'name'}."}}\\label{sec:AutoMethods".ucfirst($objects{$object}->{'name'})."}\n\n";
	print methodHndl "\\begin{description}\n";
	foreach my $method ( sort(keys(%{$objects{$object}->{'methods'}})) ) {
	    if ( $objects{$object}->{'methods'}->{$method}->{'description'} ne "UNDEFINED" ) {
		print methodHndl "\\item[]{\\tt ";
		if ( exists($objects{$object}->{'methods'}->{$method}->{'type'}) ) {
		    print methodHndl latex_encode($objects{$object}->{'methods'}->{$method}->{'type'})."\\ ";
		} else {
		    print "Warning: missing type for method ".$method." of ".$object." object\n";
		}
		print methodHndl $objects{$object}->{'methods'}->{$method}->{'method'}."(";
		if ( exists($objects{$object}->{'methods'}->{$method}->{'arguments'}) ) {
		    print methodHndl $objects{$object}->{'methods'}->{$method}->{'arguments'}
		       unless ( UNIVERSAL::isa($objects{$object}->{'methods'}->{$method}->{'arguments'},"HASH") );
		} else {
		    print "Warning: missing arguments for method ".$method." of ".$object." object\n";
		}
		print methodHndl ")} ".$objects{$object}->{'methods'}->{$method}->{'description'}."\n";
	    }
	}
	print methodHndl "\\end{description}\n";
    }
}

# Close the output files.
close(parametersHndl);
close(methodHndl);

# Write warning messages about missing data.
foreach my $parameterRead ( keys(%parametersRead) ) {
    unless ( exists($parametersData{$parameterRead}) ) {
	print "Warning: missing data for input parameter [".$parameterRead."] in file ".$parametersRead{$parameterRead}."\n";
    }
}

# Write warning messages about unknown attachments.
foreach ( keys(%parametersRead) ) {
    print "Warning: unknown attachment for input parameter [".$_."] in file ".$parametersData{$_}->{'file'}."\n"
	if ( $parametersData{$_}->{'attachmentStatus'} eq "unknown" );
}

# Write warning messages about empty default values.
foreach ( @emptyDefaults ) {
    print "Warning: empty defaultValue for input parameter [".$_->{'name'}."] in file ".$_->{'file'}."\n";
}

# Construct enumeration documentation.
open(oHndl,">".$outputRoot."Enumerations.tex");
open(sHndl,">".$outputRoot."EnumerationSpecifiers.tex");
foreach ( sort(keys(%enumerations)) ) {
    print oHndl "\\subsubsection{\\large {\\tt ".$_."}}\\hypertarget{ht:AutoEnumerations".ucfirst($_)."}{}\\label{sec:AutoEnumerations".ucfirst($_)."}\\index{enumerations!".$_."\@{\\tt ".$_."}}\n\n";
    print oHndl "\\begin{tabular}{rp{130mm}}\n";
    print oHndl "Description: & ".$enumerations{$_}->{'description'}." \\\\\n";
    print oHndl "Provided by: & {\\tt module} \\hyperlink{".$enumerations{$_}->{'file'}.":".lc($enumerations{$_}->{'module'})."}{\\tt ".latex_encode($enumerations{$_}->{'module'})."} \\\\\n";
    my $first = 1;
    foreach my $entry ( @{$enumerations{$_}->{'entry'}} ) {
	print oHndl "Members:"
	    if ( $first == 1 );
	print oHndl " & {\\tt ".$entry."}\\\\\n";
	$first = 0;
    }
    print oHndl "\\end{tabular}\n";
    print sHndl "\\def\\enum".ucfirst($_)."{\\textcolor{red}{\\hyperlink{ht:AutoEnumerations".ucfirst($_)."}{\\textless ".$_."\\textgreater}}}\n";
}
close(oHndl);
close(sHndl);

exit;

#!/usr/bin/env perl
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use lib './perl';
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use LaTeX::Encode;
use Fcntl qw(SEEK_SET);
use UNIVERSAL;
use Fortran::Utils;
use Storable qw(dclone);

# Scan Fortran90 source code and extract various useful data from "!@" lines.
# Andrew Benson 12-Mar-2010

# Get source directory
die 'Usage: Extract_Data.pl <sourceDir> <outputRoot>'
    unless ( scalar(@ARGV) == 2 );
my $sourceDir  = $ARGV[0];
my $outputRoot = $ARGV[1];

# Specify regex for derived-type declarations.
my $derivedTypeOpenRegex  = "^\\s*type\\s*((,\\s*abstract\\s*|,\\s*public\\s*|,\\s*private\\s*|,\\s*extends\\s*\\([a-zA-Z0-9_]+\\)\\s*)*)(::)??\\s*([a-zA-Z0-9_]+)\\s*\$";
my $derivedTypeCloseRegex = "^\\s*end\\s+type\\s+([a-z0-9_]+)";

# Specify regex for type-bound procedures.
my $typeBoundRegex = "^\\s*(procedure|generic)\\s*(,\\s*nopass\\s*)*::\\s*([a-z0-9_]+)\\s*=>\\s*([a-z0-9_,\\s]+)\$";
my $genericRegex   = "^\\s*procedure\\s*(,\\s*nopass\\s*)*::\\s*([a-z0-9_]+)\\s*\$";

# Create XML object.
my $xml = new XML::Simple;

# Open the source directory.
opendir(dirHndl,$sourceDir);
# Read files from the source directory.
my @fileNames = ( "work/build/objects.nodes.components.Inc" );
while ( my $fileName = readdir(dirHndl) ) {
    # Find Fortran 90 and C++ source files.
    if ( ( $fileName =~ m/\.F90$/ || $fileName =~ m/\.cpp$/ ) && $fileName !~ m/^\.\#/ ) {
	(my $preProcessedName = $fileName) =~ s/\.F90$/.p.F90/;
	if ( -e "work/build/".$preProcessedName ) {
	    push(@fileNames,"work/build/".$preProcessedName);
	} else {
	    push(@fileNames,$sourceDir."/".$fileName);
	}
    }
}
closedir(dirHndl);

# Initialize data hashes.
my %objects;
my %enumerations;
my @emptyDefaults;

foreach my $fileName ( @fileNames ) {
    # Get a printable file name.
    (my $fileNamePrint = $fileName) =~ s/_/\\_/g;
    (my $leafName = $fileName) =~ s/^$sourceDir\///;
    (my $leafNamePrint = $leafName) =~ s/_/\\_/g;
    $leafNamePrint =~ s/^work\/build\/(.*)\.p\.F90/$1.F90/;
    
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
		&Fortran::Utils::Get_Fortran_Line($fileHandle,$rawLine,$processedLine,$bufferedComments);
	    } else {
		$rawLine = <$fileHandle>;
	    }
	    
	    # Detect include files, and recurse into them.
	    if ( defined($processedLine) && $processedLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
		my $includeFile = $1;
		$includeFile =~ s/\.inc/.Inc/;
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
		    $inDerivedType = $4;
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
                    if ( $1 eq "generic" ) {
			my $methodList = lc($4);
			$methodList =~ s/^\s*//;
			$methodList =~ s/\s*$//;
			my @methods = split(/\s*,\s*/,$methodList);
			foreach ( @methods ) {
			    if ( exists($objects{lc($inDerivedType)}->{'methods'}->{$_}) ) {
				$objects{lc($inDerivedType)}->{'methods'}->{$_}->{'description'} = "GENERIC";
			    }
			}
                    }


		}
		if ( $processedLine =~ m/$genericRegex/i ) {
		    my $procedure = lc($2);
		    $objects{lc($inDerivedType)}->{'procedures'}->{$procedure} = 1;
		}
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
		    if ( $dataType eq "enumeration" ) {
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
		    } elsif ( $dataType eq "objectMethod" ) {
			my $object = lc($contents->{'object'});
			$objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'method'     } = $contents->{'method'     };
			$objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'description'} = $contents->{'description'};
			$objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'arguments'  } = $contents->{'arguments'  }
			if ( exists($contents->{'arguments'}) );
			$objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'type'       } = $contents->{'type'       }
			if ( exists($contents->{'type'}) );
			$objects{$object}->{'methods'}->{lc($contents->{'method'})}->{'inherited'  } = 0;
		    } elsif ( $dataType eq "objectMethods" ) {
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
			    $objects{$object}->{'methods'}->{lc($method->{'method'})}->{'inherited'  } = 0;
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
foreach my $object ( sort(keys(%objects)) ) {
    $objects{$object}->{'isFunctionClass'} = 0;
}
my $foundExtensions = 1;
while ( $foundExtensions == 1 ) {
    $foundExtensions = 0;
    foreach my $object ( sort(keys(%objects)) ) {
	if ( defined($objects{$object}->{'extends'}) ) {
	    $foundExtensions = 1;
	    my $parent = $objects{$object}->{'extends'};
	    $objects{$object}->{'isFunctionClass'} = 1
		if ( $parent eq "functionclass" || $objects{$parent}->{'isFunctionClass'} );
	    $objects{$object}->{'parent'} = $objects{$parent}->{'name'};
	    push(@{$objects{$parent}->{'children'}},$objects{$object}->{'name'});
	    unless ( defined($objects{$parent}->{'extends'}) ) {
		foreach my $method ( keys(%{$objects{$parent}->{'methods'}}) ) {
		    if ( ( ! exists($objects{$object}->{'methods'}->{$method}) ) || $objects{$object}->{'methods'}->{$method}->{'description'} eq "UNDEFINED" ) {
			$objects{$object}->{'methods'}->{$method} = dclone($objects{$parent}->{'methods'}->{$method});
		    }
		    $objects{$object}->{'methods'}->{$method}->{'inherited'} = 1;
		}
		foreach my $procedure ( keys(%{$objects{$parent}->{'procedures'}}) ) {
		    delete($objects{$object}->{'methods'}->{$procedure})
			if (
			    exists($objects{$object}->{'methods'}->{$procedure})
			    &&
			           $objects{$object}->{'methods'}->{$procedure}->{'description'} eq "UNDEFINED"
			);
		}
		undef($objects{$object}->{'extends'});
	    }
	}
    }
}

# Open output files.
open(methodHndl    ,">".$outputRoot."Methods.tex");

# Write method descriptions.
my @functionClassExcludes = ( "allowedParameters", "autoHook", "deepCopy", "deepCopyReset", "descriptor", "hashedDescriptor", "objectType", "stateRestore", "stateStore" );
foreach my $object ( sort(keys(%objects)) ) {
    my $hasEntries = 0;
    foreach my $method ( sort(keys(%{$objects{$object}->{'methods'}})) ) {
	if ( $objects{$object}->{'methods'}->{$method}->{'description'} eq "UNDEFINED" ) {
	    print "Warning: missing data for method ".$method." of ".$object." object\n";
	} else {
	    $hasEntries = 1;
	}
    }
    print methodHndl "\\subsection{\\large {\\normalfont \\ttfamily ".$objects{$object}->{'name'}."}}\\label{class:".$objects{$object}->{'name'}."}\\hyperdef{class}{".$objects{$object}->{'name'}."}{}\n\n";
    if ( $hasEntries == 1 ) {
	my @methodNames = 
	    map 
	{
	      $objects{$object}->{'methods'}->{$_}->{'description'} ne "UNDEFINED" 
	    &&
	      $objects{$object}->{'methods'}->{$_}->{'description'} ne "GENERIC" 
	    &&
	    (
	     ! exists($objects{$object}->{'methods'}->{$_}->{'inherited'  })
	     ||
	     !        $objects{$object}->{'methods'}->{$_}->{'inherited'  }
	    )
	    ?
	    $_
	    :
	    ()
	}
	sort(keys(%{$objects{$object}->{'methods'}}));
	if ( $objects{$object}->{'isFunctionClass'} ) {
	    my @methodNamesUnique;
	    foreach my $methodName ( @methodNames ) {
		push(@methodNamesUnique,$methodName)
		    unless ( grep {lc($_) eq lc($methodName)} @functionClassExcludes );
	    }
	    @methodNames = @methodNamesUnique;
	}
	print methodHndl "\\emph{Physics model:} \\refPhysics{".$objects{$object}->{'name'}."}\n\n"
	    if ( $objects{$object}->{'isFunctionClass'} );
	print methodHndl "\\noindent\\emph{Parent class:} \\refClass{".$objects{$object}->{'parent'}."}\n\n"
	    if ( exists($objects{$object}->{'parent'}) );
	if ( exists($objects{$object}->{'children'}) ) {
	    my @sortedChildren = sort(@{$objects{$object}->{'children'}});
		    print methodHndl "\\noindent\\emph{Child classes:}\n\n\\begin{tabular}{ll}\n";
		    for(my $i=0;$i<scalar(@sortedChildren);$i+=2) {
			print methodHndl    "\\refClass{".$sortedChildren[$i  ]."}";
			print methodHndl " & \\refClass{".$sortedChildren[$i+1]."}"
			    if ( $i+1 < scalar(@sortedChildren) );
			print methodHndl "\\\\\n";
		    }
		    print methodHndl "\\end{tabular}\n\n";
	}
	if ( scalar(@methodNames) > 0 || $objects{$object}->{'isFunctionClass'} ) {
	    print methodHndl "\\begin{description}\n";
	    print methodHndl "\\item[] All standard \\refClass{functionClass} methods (see \\S\\ref{sec:functionClassAll})\n"
		if ( $objects{$object}->{'isFunctionClass'} );
	    foreach my $method ( @methodNames ) {
		print methodHndl "\\item[]{\\normalfont \\ttfamily ";
	        if ( exists($objects{$object}->{'methods'}->{$method}->{'type'}) ) {
		    (my $methodLabel = $objects{$object}->{'methods'}->{$method}->{'type'}) =~ s/([^\\])_/$1\\_/g;
		    print methodHndl $methodLabel."\\ ";
		} else {
		    print "Warning: missing type for method ".$method." of ".$object." object\n";
		}
	        print methodHndl latex_encode($objects{$object}->{'methods'}->{$method}->{'method'})."(";
	        if ( exists($objects{$object}->{'methods'}->{$method}->{'arguments'}) ) {
		    unless ( UNIVERSAL::isa($objects{$object}->{'methods'}->{$method}->{'arguments'},"HASH") ) {
			(my $arguments = $objects{$object}->{'methods'}->{$method}->{'arguments'}) =~ s/([^\\])_/$1\\_/g;
			print methodHndl $arguments;
		    }
		} else {
		    print "Warning: missing arguments for method ".$method." of ".$object." object\n";
		}
	        print methodHndl ")} ".$objects{$object}->{'methods'}->{$method}->{'description'}."\n";
            }
	    print methodHndl "\\end{description}\n";
        }
    }
}

# Close the output files.
close(methodHndl);

# Write warning messages about empty default values.
foreach ( @emptyDefaults ) {
    print "Warning: empty defaultValue for input parameter [".$_->{'name'}."] in file ".$_->{'file'}."\n";
}

# Construct enumeration documentation.
open(oHndl,">".$outputRoot."Enumerations.tex");
open(sHndl,">".$outputRoot."EnumerationSpecifiers.tex");
foreach ( sort(keys(%enumerations)) ) {
    print oHndl "\\subsection{\\large {\\normalfont \\ttfamily ".$_."}}\\hypertarget{ht:AutoEnumerations".ucfirst($_)."}{}\\label{sec:AutoEnumerations".ucfirst($_)."}\\index{enumerations!".$_."\@{\\normalfont \\ttfamily ".$_."}}\n\n";
    print oHndl "\\begin{tabular}{rp{130mm}}\n";
    print oHndl "Description: & ".$enumerations{$_}->{'description'}." \\\\\n";
    (my $enumerationFile = $enumerations{$_}->{'file'}) =~ s/^work\/build\/(.*)\.p\.F90/$1.F90/;
    $enumerationFile =~ s/\./_/g;
    print oHndl "Provided by: & {\\normalfont \\ttfamily module} \\href{https://github.com/galacticusorg/galacticus/releases/download/masterRelease/Galacticus_Source.pdf\\#source.".$enumerationFile.":".lc($enumerations{$_}->{'module'})."}{\\normalfont \\ttfamily ".latex_encode($enumerations{$_}->{'module'})."} \\\\\n";
    my $first = 1;
    foreach my $entry ( @{$enumerations{$_}->{'entry'}} ) {
	print oHndl "Members:"
	    if ( $first == 1 );
	print oHndl " & {\\normalfont \\ttfamily ".$entry."}\\\\\n";
	$first = 0;
    }
    print oHndl "\\end{tabular}\n";
    print sHndl "\\def\\enum".ucfirst($_)."{\\textcolor{red}{\\hyperlink{ht:AutoEnumerations".ucfirst($_)."}{\\textless ".$_."\\textgreater}}}\n";
}
close(oHndl);
close(sHndl);

exit;

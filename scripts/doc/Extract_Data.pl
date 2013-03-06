#!/usr/bin/env perl
use lib './perl';
use Switch;
use XML::Simple;
use Data::Dumper;

# Scan Fortran90 source code and extract various useful data from "!@" lines.
# Andrew Benson 12-Mar-2010

# Get source directory
unless ($#ARGV == 1) {die 'Usage: Extract_Data.pl <sourceDir> <outputRoot>'};
$sourceDir  = $ARGV[0];
$outputRoot = $ARGV[1];

# Create XML object.
$xml = new XML::Simple;

# Open the source directory.
opendir(dirHndl,$sourceDir);

# Read files from the source directory.
while ( $fileName = readdir(dirHndl) ) {

    # Find Fortran 90 and C++ source files.
    if ( ( $fileName =~ m/\.F90$/ || $fileName =~ m/\.cpp$/ ) && $fileName !~ m/^\.\#/ ) {

	# Get a printable file name.
	($fileNamePrint = $fileName) =~ s/_/\\_/g;

	# Open the file.
	open(fileHndl,$sourceDir."/".$fileName);

	# Initialize the data directive buffer.
	$dataDirective  = "";
	$openingElement = "";

	# Initialize the program unit stack.
	undef(@programUnits);

	# Scan the file.
	while ( $line = <fileHndl> ) {

	    # Check for entering, leaving program units.
	    if ( $line =~ m/^\s*module\s+([a-zA-Z0-9_]+)/ ) {
		$moduleName = $1;
		unless ( $moduleName eq "procedure" ) {$programUnits[++$#programUnits] = "module:".$moduleName};
	    }
	    if ( $line =~ m/^\s*program\s+([a-zA-Z0-9_]+)/ ) {$programUnits[++$#programUnits] = "program:".$1};
	    if ( $line =~ m/^\s*(pure\s+|elemental\s+|recursive\s+)*\s*subroutine\s+([a-zA-Z0-9_]+)/ ) {$programUnits[++$#programUnits] = "subroutine:".$2};
	    if ( $line =~ m/^\s*(pure\s+|elemental\s+|recursive\s+)*\s*(real|integer|double precision|character|logical)*\s*(\((kind|len)=[\w\d]*\))*\s*function\s+([a-zA-Z0-9_]+)/ ) {$programUnits[++$#programUnits] = "function:".$5};
	    if ( $line =~ m/^\s*end\s+(program|module|subroutine|function)\s/ ) {--$#programUnits};

	    # Check for parameter reading lines - used to check that all parameters are documented.
	    if ( $line =~ m/^\s*call\s+get_input_parameter\s*\(\s*[\"\']\s*(\w+)\s*[\"\']/i ) {
		$parametersRead{$1} .= $fileName." ";
	    }

	    # Search for "!@".
	    $process = 0;
	    if ( $line =~ m/^\s*(\!|\/\/)\@\s/ ) {
		# Found directive - add it to the buffer.
		if ( $dataDirective eq "" && $line =~ m/<([^>]+)>/ ) {
		    $openingElement = $1;
		} elsif ( $line =~ m/<\/$openingElement>/ ) {
		    $process = 1;
		}
		$line =~ s/^\s*(\!|\/\/)\@\s*//;
		$line =~ s/\s*$//;
		$dataDirective .= $line." ";
	    } elsif ( $dataDirective ne "" ) {
		$process = 1;
	    }
	    if ( $process == 1 ) {
		# Process a directive.
		$openingElement = "";
		$data = eval{$xml->XMLin($dataDirective, KeepRoot => 1)};
		$lineNumber = $.;
		die("Extract_Data.pl failed in file ".$sourceDir."/".$fileName." at line ".$lineNumber." with message:\n".$@."and data \n".$dataDirective) if ($@);

		# Loop over all data types found.
		foreach $dataType ( keys(%{$data}) ) {
		    $contents = $data->{$dataType};
		    switch ( $dataType ) {
			case ( "inputParameter" ) {
			    ($printName = $contents->{'name'}) =~ s/_/\\_/g;
			    $buffer  = "\\noindent {\\bf Name:} {\\tt ".$printName."}\\\\\n";
			    $buffer .= "{\\bf Attached to:} ";
			    if ( exists($contents->{'attachedTo'}) ) {
				$programUnitIndex = $#programUnits;
				$regEx = $contents->{'attachedTo'}.":";
				until ( $programUnits[$programUnitIndex] =~ m/$regEx/ || $programUnitIndex == -1 ) {
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
				$attachedTo = "{\\tt ".$programUnits[$#programUnits]."}";
				$attachedAt = $#programUnits;
			    }

			    $attachedLink = $fileName.":";
			    for($iAttach=0;$iAttach<=$attachedAt;++$iAttach) {
				if ( $programUnits[$iAttach] =~ m/^[^:]+:(.*)/ ) {
				    $unitName = lc($1);
				    $unitName =~ s/\\_/_/g;
				    $attachedLink .= $unitName.":";
				}
			    }
			    $attachedLink =~ s/:$//;

			    if ( $attachedTo =~ m/:/ ) {
				$targetPrefix = "\\hyperlink{".$attachedLink."}{";
				$targetPrefix =~ s/\{\\tt\s+([^\}]+)\}/\1/;
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
			    $buffer .= "{\\bf File:} \\hyperlink{".$fileName."}{{\\tt ".$fileNamePrint."}}\\\\\n";
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
			    $parametersData{$contents->{'name'}} = $buffer;
			}
			case ( "objectMethod" ) {
			    $object = $contents->{'object'};
			    $objects{$object}->{$contents->{'method'}} = $contents->{'description'};
			}
			case ( "objectMethods" ) {
			    $object = $contents->{'object'};
			    foreach $method ( @{$contents->{'objectMethod'}} ) {
				$objects{$object}->{$method->{'method'}} = $method->{'description'};
			    }
			}
		    }
		}

		# Reset the directive buffer.
		$dataDirective = "";
	    }

	}

	# Close the file.
	close(fileHndl);

    }

}

# Close the source directory.
closedir(dirHndl);

# Open output files.
open(parametersHndl,">".$outputRoot."Parameters.tex");
open(methodHndl    ,">".$outputRoot."Methods.tex");

# Write parameter descriptions.
@sortedParameters = sort { $parametersData{$a} cmp $parametersData{$b} } keys %parametersData;
foreach $parameterName ( @sortedParameters ) {
    print parametersHndl $parametersData{$parameterName};
}

# Write method descriptions.
foreach $object ( sort(keys(%objects)) ) {
    print methodHndl "\\subsubsection{{\\tt ".$object."}}\\label{sec:AutoMethods".ucfirst($object)."}\n\n";
    print methodHndl "\\begin{description}\n";
    foreach $method ( sort(keys(%{$objects{$object}})) ) {
	print methodHndl "\\item[{\\tt ".$method."}] ".$objects{$object}->{$method}."\n";
    }
    print methodHndl "\\end{description}\n";
}

# Close the output files.
close(parametersHndl);
close(methodHndl);

# Write warning messages about missing data.
foreach $parameterRead ( keys(%parametersRead) ) {
    unless ( exists($parametersData{$parameterRead}) ) {
	print "Warning: missing data for input parameter [".$parameterRead."] in file ".$parametersRead{$parameterRead}."\n";
    }
}

exit;

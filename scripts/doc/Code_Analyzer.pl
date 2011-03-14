#!/usr/bin/env perl
use lib './perl';
use File::Find;
use Fortran::Utils;
use Data::Dumper;
use Cwd;
use Text::Balanced qw (extract_bracketed);
use LaTeX::Encode;

# Scan Fortran 90 source code and output a Latex file describing it.
#
# Currently outputs lists of modules used, and used by with hyperlinks, along with description and code line count. Collects much
# more information on subroutine calls, function calls, module procedures etc. Could be made to output information on these if
# required.
#
# Andrew Benson (01-May-2010)

# Get arguments.
die("Usage: Code_Analyzer.pl <sourceDir> <outputFile>") unless ( $#ARGV == 1 );
push(@sourceDir,$ARGV[0]);
$outputFile = $ARGV[1];

# Get current working directory.
$currentDirectory = cwd;

# Specify location of include files.
$includeFileDirectory = $currentDirectory."/".$sourceDir[0]."/../work/build";

# Specify regex for source code comments.
$commentRegex = "^\\s*\\%";

# Specify regex for use statements.
$useRegex = "^\\s*use\\s*(,\\s*intrinsic\\s*)*(::)??\\s*([a-z0-9_]+)";

# Specify regex for subroutine calls.
$callRegex = "(^|\\W)call\\s+([a-z0-9_]+)\\s*(\\(|\\s*\$)";

# Specify regex for calls to type-bound subroutines.
$callTypeBoundRegex = "(^|\\W)call\\s+([a-z0-9_]+)\\s*%\\s*([a-z0-9_]+)\\s*(\\(|\\s*\$)";

# Specify regex for module procedures.
$moduleProcedureRegex = "^\\s*module\\s+procedure\\s+([a-z0-9_]+)\\s*\$";

# Specify regex for derived-type declarations.
$derivedTypeRegex = "^\\s*type\\s*\\(\\s*([a-z0-9_]+)\\s*\\)([\\sa-z0-9_,:\\+\\-\\*\\/\\(\\)]*::)*\\s*([a-z0-9_,:\\+\\-\\*\\/\\(\\)]+)\\s*\$";

# Specify regex for type-bound procedures.
$typeBoundRegex = "^\\s*(procedure|generic)\\s*::\\s*([a-z0-9_]+)\\s*=>\\s*([a-z0-9_,\\s]+)\$";

# Specify regexs for intrinsic variable declarations.
%intrinsicDeclarations = (
    integer   => { variables => 3, regEx => "^\\s*integer\\s*(\\(\\s*kind\\s*=\\s*[a-z0-9_]+\\s*\\))*([\\sa-z0-9_,:\\+\\-\\*\\/\\(\\)]*::)*\\s*([\\sa-z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    real      => { variables => 3, regEx => "^\\s*real\\s*(\\(\\s*kind\\s*=\\s*[a-z0-9_]+\\s*\\))*([\\sa-z0-9_,:\\+\\-\\*\\/\\(\\)]*::)*\\s*([\\sa-z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    double    => { variables => 3, regEx => "^\\s*double\\s+precision\\s*(\\(\\s*kind\\s*=\\s*[a-z0-9_]+\\s*\\))*([\\sa-z0-9_,:\\+\\-\\*\\/\\(\\)]*::)*\\s*([\\sa-z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    logical   => { variables => 3, regEx => "^\\s*logical\\s*(\\(\\s*kind\\s*=\\s*[a-z0-9_]+\\s*\\))*([\\sa-z0-9_,:\\+\\-\\*\\/\\(\\)]*::)*\\s*([\\sa-z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    character => { variables => 4, regEx => "^\\s*character\\s*(\\(\\s*(len|kind)\\s*=\\s*[a-z0-9_\\*\\(\\)]+\\s*\\))*([\\sa-z0-9_,:\\+\\-\\*\\/\\(\\)]*::)*\\s*([\\sa-z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    procedure => { variables => 3, regEx => "^\\s*procedure\\s*(\\(\\s*[a-z0-9_]*\\s*\\))*([\\sa-z0-9_,:\\+\\-\\*\\/\\(\\)]*::)*\\s*([\\sa-z0-9_,:=>\\+\\-\\*\\/\\(\\)]+)\\s*\$" }
    );

# Specify unit opening regexs.
%unitOpeners = (
    # Find module openings, avoiding module procedures.
    module             => { unitName => 1, regEx => "^\\s*module\\s+(?!procedure\\s)([a-z0-9_]+)" },
    # Find program openings.
    program            => { unitName => 1, regEx => "^\\s*program\\s+([a-z0-9_]+)" },
    # Find subroutine openings, allowing for pure, elemental and recursive subroutines.
    subroutine         => { unitName => 2, regEx => "^\\s*(pure\\s+|elemental\\s+|recursive\\s+)*\\s*subroutine\\s+([a-z0-9_]+)"},
    # Find function openings, allowing for pure, elemental and recursive functions, and different function types.
    function           => { unitName => 5, regEx => "^\\s*(pure\\s+|elemental\\s+|recursive\\s+)*\\s*(real|integer|double\\s+precision|character|logical)*\\s*(\\((kind|len)=[\\w\\d]*\\))*\\s*function\\s+([a-z0-9_]+)"},
    # Find interfaces.
    interface          => { unitName => 2, regEx => "^\\s*(abstract\\s+)??interface\\s+([a-z0-9_\\(\\)\\/\\+\\-\\*\\.=]*)"},
    # Find types.
    type               => { unitName => 3, regEx => "^\\s*type\\s*(,\\s*public\\s*|,\\s*private\\s*)*(::)??\\s*([a-z0-9_]+)\\s*\$"}
    );

# Specify unit closing regexs.
%unitClosers = (
    module             => { unitName => 1, regEx => "^\\s*end\\s+module\\s+([a-z0-9_]+)" },
    program            => { unitName => 1, regEx => "^\\s*end\\s+program\\s+([a-z0-9_]+)" },
    subroutine         => { unitName => 1, regEx => "^\\s*end\\s+subroutine\\s+([a-z0-9_]+)"},
    function           => { unitName => 1, regEx => "^\\s*end\\s+function\\s+([a-z0-9_]+)"},
    interface          => { unitName => 1, regEx => "^\\s*end\\s+interface"},
    type               => { unitName => 1, regEx => "^\\s*end\\s+type\\s+([a-z0-9_]+)"}
    );

# Find and process all files in the source directory tree.
find(\&processFile, @sourceDir);

#print Dumper(%units); 

# Create list of module IDs.
&Build_Modules_Hash;

# Output the data.
&Output_Data;

exit;

sub processFile {

    # Process a file in the source directory tree. Scans the file and extracts information on units, variables, types, calls etc.
    $fileName = $_;
    chomp($fileName);

    # Check if this is a Fortran source file.
    if ( $fileName =~ m/\.[fF]90$/ ) {
	# Initialize the unitIdList array and the units hash.
	undef(@unitIdList);
	$unitIdList[0] = $fileName;
	$units{$fileName} = {
		unitType => "file",
		unitName => $fileName
	};

	# Add the file to the list of filenames to process.
	undef(@fileNames);
	undef(@filePositions);
	unshift(@fileNames,$fileName);
	unshift(@filePositions,-1);

	# Process files until none remain.
	while ( $#fileNames >= 0 ) {

	    # Open the file.
	    open($fileHandle,$fileNames[0]);
	    seek($fileHandle,$filePositions[0],SEEK_SET) unless ( $filePositions[0] == -1 );

	    # Process until end of file is reached.
	    LINE: until ( eof($fileHandle) ) {

		# Specify that line has yet to be processed.
		$lineProcessed = 0;
		
		# Grab the next Fortran line.
		&Fortran_Utils::Get_Fortran_Line($fileHandle,$rawLine,$processedLine,$bufferedComments);
		foreach $unitID ( @unitIdList ) {
		    ++$units{$unitID}->{"codeLines"};
		}

                # Detect include files, and recurse into them.
		if ( $processedLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
		    $includeFile = $includeFileDirectory."/".$1;
		    if ( -e $includeFile ) {
			$filePositions[0] = tell($fileHandle);
			unshift(@fileNames,$includeFile);
			unshift(@filePositions,-1);
			last;
		    }
		}

		# Detect unit opening.
		foreach $unitType ( keys(%unitOpeners) ) {
		    
		    # Check for a match to a unit opening regex.
		    if ( $processedLine =~ m/$unitOpeners{$unitType}->{"regEx"}/i ) {
			$matchIndex = $unitOpeners{$unitType}->{"unitName"};
			$unitName   = lc($$matchIndex);
			# Get the ID of the parent unit.
			$parentID = $unitIdList[$#unitIdList];
			# Generate identifier for this unit.
			$unitID = $parentID.":".$unitName;
			# Add this unit to the "contains" list for its parent.
			$units{$parentID}->{"contains"}->{$unitID} = 1;
			# Create an entry for this new unit.
			push(@unitIdList,$unitID);
			$units{$unitID} = {
			    unitType => $unitType,
			    unitName => $unitName,
			    belongsTo => $parentID
			};
			# Mark line as processed.
			$lineProcessed = 1;
		    }
		}
		if ( $lineProcessed == 1 ) {next LINE};
		
		# Detect unit closing.
		foreach $unitType ( keys(%unitClosers) ) {
		    
		    # Check for a match to a unit closing regex.
		    if ( $processedLine =~ m/$unitClosers{$unitType}->{"regEx"}/i ) {
			$matchIndex = $unitClosers{$unitType}->{"unitName"};
			$unitName   = lc($$matchIndex);
			# Get ID of last unit opening.
			$openerID = $unitIdList[$#unitIdList];
			# Check we have a matching opening.
			unless ( $unitType eq $units{$openerID}->{"unitType"} && ( $unitName eq $units{$openerID}->{"unitName"} || $unitType =~ m/interface/i ) ) {
			    print "Unit close does not match unit open:\n";
			    print " Closing with: ".$unitType." ".$unitName."\n";
			    print "  Opened with: ".$units{$openerID}->{"unitType"}." ".$units{$openerID}->{"unitName"}."\n";
			    print "      In file: ".$fileName."\n";
			    die;
			}
			# Remove entry from ID list.
			--$#unitIdList;
			# Mark line as processed.
			$lineProcessed = 1;
		    }
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Check for source code comments.
		if ( $bufferedComments =~ m/$commentRegex/i ) {
		    $unitId = $unitIdList[$#unitIdList];
		    ($trimmedComments = $bufferedComments) =~ s/^\s*\%//;
		    $units{$unitId}->{"comments"} .= $trimmedComments;
		    # Mark line as processed.
		    $lineProcessed = 1;
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Check for use statements.
		if ( $processedLine =~ m/$useRegex/i ) {
		    $moduleName = lc($3);
		    $unitId = $unitIdList[$#unitIdList];
		    $units{$unitId}->{"modulesUsed"}->{$moduleName} = 1;
		    # Mark line as processed.
		    $lineProcessed = 1;
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Check for direct subroutine calls.
		if ( $processedLine =~ m/$callRegex/i ) {
		    $subroutineName = lc($2);
		    $unitId = $unitIdList[$#unitIdList];
 		    $units{$unitId}->{"subroutinesCalled"}->{$subroutineName} = -1;
		    # Mark line as processed.
		    $lineProcessed = 1;
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Check for subroutine calls via type-bound procedures.
		if ( $processedLine =~ m/$callTypeBoundRegex/i ) {
		    $subroutineName = lc($3);
		    $variableName   = lc($2);
		    $unitId = $unitIdList[$#unitIdList];
		    # Store the name of the type-bound subroutine call, and the name of the variable it was applied to.
 		    $units{$unitId}->{"subroutinesCalled"}->{$subroutineName} = $variableName;
		    # Mark line as processed.
		    $lineProcessed = 1;
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Check for type-bound procedures.
		if ( $units{$unitIdList[$#unitIdList]}->{"unitType"} eq "type" ) {
		    if ( $processedLine =~ m/$typeBoundRegex/i ) {
			$methodName = lc($2);
			($procedureList = lc($3)) =~ s/\s//g;
			@procedures = split(/,/,$procedureList);
			$unitId = $unitIdList[$#unitIdList];
			@{$units{$unitId}->{"methods"}->{$methodName}} = @procedures;
			# Mark line as processed.
			$lineProcessed = 1;
		    }
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Detect intrinsic variable declarations.
		foreach $intrinsicType ( keys(%intrinsicDeclarations) ) {
		    
		    # Check for a match to an intrinsic declaration regex.
		    if ( $processedLine =~ m/$intrinsicDeclarations{$intrinsicType}->{"regEx"}/i ) {
			$matchIndex = $intrinsicDeclarations{$intrinsicType}->{"variables"};
			$variablesList = lc($$matchIndex);
			# Get ID of unit.
			$unitId = $unitIdList[$#unitIdList];
			@variables = &Extract_Variables($variablesList);
			# Store the variable list.
			push(@{$units{$unitId}->{$intrinsicType}},@variables);
			# Mark line as processed.
			$lineProcessed = 1;
		    }
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Check for derived-type declarations.
		if ( $processedLine =~ m/$derivedTypeRegex/i ) {
		    $derivedType = $1;
		    $variablesList = $3;
		    @variables = &Extract_Variables($variablesList);		    
		    $unitId = $unitIdList[$#unitIdList];
		    push(@{$units{$unitId}->{"derivedTypesUsed"}->{$derivedType}},@variables);
		    # Mark line as processed.
		    $lineProcessed = 1;
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Check for module procedures.
		if ( $processedLine =~ m/$moduleProcedureRegex/i ) {
		    $procedureName = lc($1);
		    $unitId = $unitIdList[$#unitIdList];
		    $units{$unitId}->{"moduleProcedures"}->{$procedureName} = 1;
		    # Mark line as processed.
		    $lineProcessed = 1;
		}
		if ( $lineProcessed == 1 ) {next LINE};

		# Check for function calls.
		$functionSeek = lc($processedLine);
		$functionSeek =~ s/'[^']+'//g;
		while ( $functionSeek =~ m/\(/ ) {
		    ($extracted,$remainder,$prefix) = extract_bracketed($functionSeek,"()","[^\\(]+");
		    if ( $prefix =~ m/(([a-z0-9_]+)\s*%\s*)*([a-z0-9_]+)$/ ) {
			$functionName = lc($3);
			($derivedType = lc($1)) =~ s/%\s*$//;
			if ( $derivedType eq "" ) {$derivedType = -1};
			$unitId = $unitIdList[$#unitIdList];
			$units{$unitId}->{"functionCalls"}->{$functionName} = $derivedType;
		    }
		    $functionSeek = $remainder;
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
}

sub Extract_Variables {
    # Given the post-"::" section of a variable declaration line, return an array of all variable names.
    $variableList = shift;
    # Convert to lower case.
    $variableList = lc($variableList);
    # Remove whitespace.
    $variableList =~ s/\s//g;
    # Remove *'s (can appear for character variables.
    $variableList =~ s/\*//g;
    # Remove text within matching () pairs.
    while ( $variableList =~ m/\(/ ) {
	($extracted,$remainder,$prefix) = extract_bracketed($variableList,"()","[a-z0-9_,=>\\+\\-\\*\\/]+");
	$variableList = $prefix.$remainder;
    }
    # Remove any definitions or associations.
    $variableList =~ s/=[^,]*(,|$)//g;
    # Split variables into an array and store.
    @variables = split(/,/,$variableList);
    return @variables;
}

sub Build_Modules_Hash {
    # Create a hash giving unit IDs as a function of module names.
    foreach $unitID ( keys(%units) ) {
	if ( $units{$unitID}->{"unitType"} eq "module" ) {
	    $modules{$units{$unitID}->{"unitName"}} = $unitID;
	}
    }
    # Create records of which units use each module.
    foreach $unitID ( keys(%units) ) {
	if ( exists($units{$unitID}->{"modulesUsed"}) ) {
	    foreach $moduleName ( keys(%{$units{$unitID}->{"modulesUsed"}}) ) {
		if ( exists($modules{$moduleName}) ) {
		    $units{$modules{$moduleName}}->{"usedBy"}->{$unitID} = 1;
		}
	    }
	}
    }
}

sub Output_Data {

    # Output the output file and begin a section.
    open($outputHandle,">".$outputFile);
    print $outputHandle "\\section{Program units}\n";

    # Get a list of all units.
    @unsortedIDs = keys(%units);
    # Get a sorted list of all units.
    @sortedIDs = sort(@unsortedIDs);
    # Loop over all units.
    foreach $unitID ( @sortedIDs ) {
	# Get unit type.
	$unitType = $units{$unitID}->{"unitType"};
	# Get unit name in LaTeX encoding.
	$unitName = &latex_encode($units{$unitID}->{"unitName"});
	$unitName =~ s/\\_/\\\-\\_/g;
	# Get ID of parent.
	if ( exists($units{$unitID}->{"belongsTo"}) ) {
	    $parentID = $units{$unitID}->{"belongsTo"};
	} else {
	    $parentID = "";
	}
	# Check if this is an abstract interface.
	if ( $unitType eq "interface" && $unitName eq "" ) {
	    $unitIsAbstract = 1;
	} else {
	    $unitIsAbstract = 0;
	}
	# Check if parent is an abstract interface.
	if ( $units{$parentID}->{"unitType"} eq "interface" && $units{$parentID}->{"unitName"} eq "" ) {
	    $parentUnitIsAbstract = 1;
	} else {
	    $parentUnitIsAbstract = 0;
	}
	
	# Output header line for unit, skipping nameless interfaces (i.e. abstract interfaces) and any children of such interfaces.
	unless ( $unitIsAbstract == 1 || $parentUnitIsAbstract == 1 ) {
	    print $outputHandle "\\noindent{\\bf ".$unitType.":} \\hypertarget{".$unitID."}{{\\tt ".$unitName."}}\\index[code]{".$unitName."\@\{\\tt ".$unitName."} (".$unitType.")}\n\n";

	    # Begin tabulated output.
	    $tableOpen = "\\begin{supertabular}{lp{70mm}p{70mm}}\n";
	    $tableIsOpen = 0;
	    
	    # Output any description.
	    if ( exists($units{$unitID}->{"comments"}) ) {
		unless ( $units{$unitID}->{"comments"} =~ m/{verbatim}/ ) {
		    if ( $tableIsOpen == 0 ) {
			print $outputHandle $tableOpen;
			$tableIsOpen = 1;
		    }
		    print $outputHandle "\\emph{Description:} & \\multicolumn{2}{l}{\n\\begin{minipage}[t]{140mm}\n".$units{$unitID}->{"comments"}."\\end{minipage}\n}\\\\\n";
		} else {
		    print $outputHandle "\\emph{Description:} ".$units{$unitID}->{"comments"};
		    print $outputHandle "\\\\" unless ( $units{$unitID}->{"comments"} =~ m/\}\s*$/ );
		    print $outputHandle "\n";
		}
	    }

	    # Output count of code lines.
	    if ( $tableIsOpen == 0 ) {
		print $outputHandle $tableOpen;
		$tableIsOpen = 1;
	    }
	    print $outputHandle "\\emph{Code lines:} & \\multicolumn{2}{l}{".$units{$unitID}->{"codeLines"}."} \\\\\n";

	    # Output any parent.
	    unless ( $parentID eq "" ) {
		if ( $tableIsOpen == 0 ) {
		    print $outputHandle $tableOpen;
		    $tableIsOpen = 1;
		}
		print $outputHandle "\\emph{Contained by:} & \\multicolumn{2}{l}{".$units{$parentID}->{"unitType"}." \\hyperlink{".$parentID."}{{\\tt ".&latex_encode($units{$parentID}->{"unitName"})."}}} \\\\ \n" 
	    }
	    
	    # Output lists of modules used.
	    if ( exists($units{$unitID}->{"modulesUsed"}) ) {
		if ( $tableIsOpen == 0 ) {
		    print $outputHandle $tableOpen;
		    $tableIsOpen = 1;
		}
		print $outputHandle "\\emph{Modules used:} ";
		@unsortedModules = keys(%{$units{$unitID}->{"modulesUsed"}});
		@sortedModules = sort(@unsortedModules);
		foreach ( @sortedModules ) {
		    $_ = "\\hyperlink{".$modules{$_}."}{{\\tt ".&latex_encode($_)."}}";
		}
		&printTwoColumn($outputHandle,\@sortedModules);
	    }

	    # Output lists of units which use us.
	    if ( exists($units{$unitID}->{"usedBy"}) ) {
		if ( $tableIsOpen == 0 ) {
		    print $outputHandle $tableOpen;
		    $tableIsOpen = 1;
		}
		print $outputHandle "\\emph{Used by:} ";
		@unsortedUnits = keys(%{$units{$unitID}->{"usedBy"}});
		@sortedUnits = sort(@unsortedUnits);
		foreach ( @sortedUnits ) {
		    $_ = $units{$_}->{"unitType"}." \\hyperlink{".$_."}{{\\tt ".&latex_encode($units{$_}->{"unitName"})."}}";
		}
		&printTwoColumn($outputHandle,\@sortedUnits);
	    }

	    # Close tabulated output.
	    if ( $tableIsOpen == 1 ) {
		print $outputHandle " & & \\\\\n\\end{supertabular}\n\\\\\n";
	    } else {
		print $outputHandle "\n";
	    }

	}

    }

    close($outputHandle);

}

sub printTwoColumn {
    # Print contents of an array in a two column table.
    $outputHandle = $_[0];
    @toPrint = @{$_[1]};
    $iColumn = 0;
    foreach $item ( @toPrint ) {
	print $outputHandle " & \\RaggedRight ".$item;
	print $outputHandle "\\\\\n" if ( $iColumn == 1);
	$iColumn = 1-$iColumn;
    }
    print $outputHandle " & \\\\\n" if ( $iColumn == 1 );
}

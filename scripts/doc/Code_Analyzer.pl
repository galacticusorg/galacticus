#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use File::Find;
use Data::Dumper;
use Cwd;
use Text::Balanced qw (extract_bracketed);
use LaTeX::Encode;
use Fcntl qw(SEEK_SET);
use Fortran::Utils;

# Scan Fortran 90 source code and output a LaTeX file describing it.
#
# Currently outputs lists of modules used, and used by with hyperlinks, along with description and code line count. Collects much
# more information on subroutine calls, function calls, module procedures etc. Could be made to output information on these if
# required.
#
# Andrew Benson (01-May-2010)

# Get arguments.
die("Usage: Code_Analyzer.pl <sourceDir> <outputFile>") 
    unless ( scalar(@ARGV) == 2 );
my @sourceDir;
push(@sourceDir,$ARGV[0]);
my $outputFile = $ARGV[1];

# Global data structures.
my %units;
my %modules;

# Get current working directory.
my $currentDirectory = cwd;

# Specify location of include files.
my $includeFileDirectory = $currentDirectory."/".$sourceDir[0]."/../".$ENV{'BUILDPATH'};

# Specify regex for use statements.
my $useRegex = "^\\s*use\\s*(,\\s*intrinsic\\s*)*(::)??\\s*([a-z0-9_]+)";

# Specify regex for subroutine calls.
my $callRegex = "(^|\\W)call\\s+([a-z0-9_]+)\\s*(\\(|\\s*\$)";

# Specify regex for calls to type-bound subroutines.
my $callTypeBoundRegex = "(^|\\W)call\\s+([a-z0-9_]+)\\s*%\\s*([a-z0-9_]+)\\s*(\\(|\\s*\$)";

# Specify regex for module procedures.
my $moduleProcedureRegex = "^\\s*module\\s+procedure\\s+([a-z0-9_]+)\\s*\$";

# Specify regex for derived-type declarations.
my $derivedTypeRegex = "^\\s*type\\s*\\(\\s*([a-z0-9_]+)\\s*\\)([\\sa-z0-9_,:\\+\\-\\*\\/\\(\\)]*::)*\\s*([a-z0-9_,:\\+\\-\\*\\/\\(\\)]+)\\s*\$";

# Specify regex for type-bound procedures.
my $typeBoundRegex = "^\\s*(procedure|generic)\\s*::\\s*([a-z0-9_]+)\\s*=>\\s*([a-z0-9_,\\s]+)\$";

# Specify regexs for intrinsic variable declarations.
my %intrinsicDeclarations = (
    integer   => { intrinsic => "integer", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)integer(?-i)\\s*(\\(\\s*kind\\s*=\\s*[a-zA-Z0-9_]+\\s*\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    real      => { intrinsic => "real", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)real(?-i)\\s*(\\(\\s*kind\\s*=\\s*[a-zA-Z0-9_]+\\s*\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9\\._,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    double    => { intrinsic => "double precision", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)double\\s+precision(?-i)\\s*(\\(\\s*kind\\s*=\\s*[a-zA-Z0-9_]+\\s*\\))*([\\sa-zA-Z0-9_,:=\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9\\._,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    logical   => { intrinsic => "logical", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)logical(?-i)\\s*(\\(\\s*kind\\s*=\\s*[a-zA-Z0-9_]+\\s*\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    character => { intrinsic => "character", type => 1, attributes => 4, variables => 5, regEx => "^\\s*(?i)character(?-i)\\s*(\\((\\s*(len|kind)\\s*=\\s*[a-zA-Z0-9_,\\+\\-\\*\\(\\)]+\\s*)+\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9_,:=>\\+\\-\\*\\/\\(\\)\\[\\]]+)\\s*\$" },
    procedure => { intrinsic => "procedure", type => 1, attributes => 2, variables => 3, regEx => "^\\s*(?i)procedure(?-i)\\s*(\\(\\s*[a-zA-Z0-9_]*\\s*\\))*([\\sa-zA-Z0-9_,:\\+\\-\\*\\/\\(\\)]*)??::\\s*([\\sa-zA-Z0-9_,:=>\\+\\-\\*\\/\\(\\)]+)\\s*\$" },
    );

# Specify unit opening regexs.
my %unitOpeners = (
    # Find module openings, avoiding module procedures.
    module             => { unitName => 1, regEx => "^\\s*module\\s+(?!procedure\\s)([a-z0-9_]+)\$" },
    # Find program openings.
    program            => { unitName => 1, regEx => "^\\s*program\\s+([a-z0-9_]+)" },
    # Find subroutine openings, allowing for impure, pure, elemental and recursive subroutines.
    subroutine         => { unitName => 2, regEx => "^\\s*(impure\\s+|pure\\s+|elemental\\s+|recursive\\s+|module\\s+)*\\s*subroutine\\s+([a-z0-9_]+)"},
    # Find function openings, allowing for impure, pure, elemental and recursive functions, and different function types.
    function           => { unitName => 6, regEx => "^\\s*(impure\\s+|pure\\s+|elemental\\s+|recursive\\s+)*\\s*(real\\s*|integer\\s*|double\\s+precision\\s*|double\\s+complex\\s*|character\\s*|logical\\s*|module\\s*)*\\s*(\\(((kind|len)=)??[\\w\\d]*\\))*\\s*function\\s+([a-z0-9_]+)"},
    # Find interfaces.
    interface          => { unitName => 2, regEx => "^\\s*(abstract\\s+)??interface\\s+([a-z0-9_\\(\\)\\/\\+\\-\\*\\.=]*)"},
    # Find types.
    type               => { unitName => 3, regEx => "^\\s*type\\s*(,\\s*abstract\\s*|,\\s*public\\s*|,\\s*private\\s*|,\\s*extends\\s*\\([a-zA-Z0-9_]+\\)\\s*|,\\s*bind\\(c\\)\\s*)*(::)??\\s*([a-z0-9_\{\}¦]+)\\s*\$"}
    );

# Specify unit closing regexs.
my %unitClosers = (
    module             => { unitName => 1, regEx => "^\\s*end\\s+module\\s+([a-z0-9_]+)" },
    program            => { unitName => 1, regEx => "^\\s*end\\s+program\\s+([a-z0-9_]+)" },
    subroutine         => { unitName => 1, regEx => "^\\s*end\\s+subroutine\\s+([a-z0-9_]+)"},
    function           => { unitName => 1, regEx => "^\\s*end\\s+function\\s+([a-z0-9_]+)"},
    interface          => { unitName => 1, regEx => "^\\s*end\\s+interface"},
    type               => { unitName => 1, regEx => "^\\s*end\\s+type\\s+([a-z0-9_\{\}¦]+)"}
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
    my $fileName = $_;
    chomp($fileName);

    # Check if this is a Fortran or C++ source file.
    if ( ( $fileName =~ m/\.[fF]90$/ || $fileName =~ m/\.cpp$/ ) && $fileName !~ m/^\.\#/ ) {

	# Initialize the unitIdList array and the units hash.
	my @unitIdList = ( $fileName );
	$units{$fileName} = {
	    unitType => "file",
	    unitName => $fileName
	};

	# Process further only for Fortran.
	if ( $fileName =~ m/\.[fF]90$/ ) {
	    
	    # Add the file to the list of filenames to process.
	    (my $preProcessedFile = $fileName) =~ s/\.F90$/.p.F90/;
	    my @fileStack =
		(
		 {
		     name     => $fileName,
		     position => -1,
		     inXML    => 0,
		     inLaTeX  => 0
		 }
		);
	    
	    # Process files until none remain.
	    while ( scalar(@fileStack) > 0 ) {

		# Open the file.
		open(my $fileHandle,$fileStack[0]->{'name'});
		seek($fileHandle,$fileStack[0]->{'position'},SEEK_SET) unless ( $fileStack[0]->{'position'} == -1 );
		
		# Process until end of file is reached.
	      LINE: until ( eof($fileHandle) ) {
		  
		  # Specify that line has yet to be processed.
		  my $lineProcessed = 0;
		  
		  # Grab the next Fortran line.
		  &Fortran::Utils::Get_Fortran_Line($fileHandle,my $rawLine,my $processedLine,my $bufferedComments);
		  # Detect leaving LaTeX and XML blocks.
		  $fileStack[0]->{'inXML'  } = 0
		      if ( $rawLine =~ m/^\s*!!\]/ );
		  $fileStack[0]->{'inLaTeX'} = 0
		      if ( $rawLine =~ m/^\s*!!\}/ );
		  foreach my $unitID ( @unitIdList ) {
		      ++$units{$unitID}->{"codeLines"};
		  }

		  # Detect include files, and recurse into them.
		  if ( $processedLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
		      my $includeFile = $includeFileDirectory."/".$1;
		      if ( -e $includeFile ) {
			  $fileStack[0]->{'position'} = tell($fileHandle);
			  unshift(
			      @fileStack,
			      {
				  name     => $includeFile,
				  position => -1,
				  inXML    => 0,
				  inLaTeX  => 0
			      }
			      );
			  last;
		      }
		  }

		  # Detect unit opening.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      foreach my $unitType ( keys(%unitOpeners) ) {			  
			  # Check for a match to a unit opening regex.
			  if ( my @matches = $processedLine =~ m/$unitOpeners{$unitType}->{"regEx"}/i ) {
			      my $matchIndex = $unitOpeners{$unitType}->{"unitName"};
			      my $unitName   = lc($matches[$matchIndex-1]);
			      # Get the ID of the parent unit.
			      my $parentID = $unitIdList[-1];
			      # Generate identifier for this unit.
			      my $unitID = $parentID.":".$unitName;
			      $unitID =~ s/[¦\{\}]//g;
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
		  }

		  # Detect unit closing.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      foreach my $unitType ( keys(%unitClosers) ) {
			  # Check for a match to a unit closing regex.
			  if ( my @matches = $processedLine =~ m/$unitClosers{$unitType}->{"regEx"}/i ) {
			      my $matchIndex = $unitClosers{$unitType}->{"unitName"};
			      my $unitName   = lc($matches[$matchIndex-1]);
			      # Get ID of last unit opening.
			      my $openerID = $unitIdList[-1];
			      # Check we have a matching opening.
			      unless ( $unitType eq $units{$openerID}->{"unitType"} && ( $unitName eq $units{$openerID}->{"unitName"} || $unitType =~ m/interface/i ) ) {
				  print "Unit close does not match unit open:\n";
				  print " Closing with: ".$unitType." ".$unitName."\n";
				  print "  Opened with: ".$units{$openerID}->{"unitType"}." ".$units{$openerID}->{"unitName"}."\n";
				  print "      In file: ".$fileName."\n";
				  die;
			      }
			      # Remove entry from ID list.
			      pop(@unitIdList);
			      # Mark line as processed.
			      $lineProcessed = 1;
			  }
		      }
		  }

		  # Check for source code comments.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} ) {
		      if ( $bufferedComments =~ m/^!\{/ ) {
			  my $unitId = $unitIdList[-1];
			  until ( eof($fileHandle) ) {
			      my $commentLine = <$fileHandle>;
			      last
				  if ( $commentLine =~ m/^\s*!!\}/ );
			      $commentLine =~ s/^\s*!//;
			      $units{$unitId}->{"comments"} .= $commentLine;
			  }
			  $fileStack[0]->{'inLaTeX'} = 0;
			  # Mark line as processed.
			  $lineProcessed = 1;
		      }
		  }

		  # Check for use statements.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      if ( $processedLine =~ m/$useRegex/i ) {
			  my $moduleName = lc($3);
			  my $unitId = $unitIdList[-1];
			  $units{$unitId}->{"modulesUsed"}->{$moduleName} = 1;
			  # Mark line as processed.
			  $lineProcessed = 1;
		      }
		  }

		  # Check for direct subroutine calls.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      if ( $processedLine =~ m/$callRegex/i ) {
			  my $subroutineName = lc($2);
			  my $unitId = $unitIdList[-1];
			  $units{$unitId}->{"subroutinesCalled"}->{$subroutineName} = -1;
			  # Mark line as processed.
			  $lineProcessed = 1;
		      }
		  }

		  # Check for subroutine calls via type-bound procedures.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      if ( $processedLine =~ m/$callTypeBoundRegex/i ) {
			  my $subroutineName = lc($3);
			  my $variableName   = lc($2);
			  my $unitId = $unitIdList[-1];
			  # Store the name of the type-bound subroutine call, and the name of the variable it was applied to.
			  $units{$unitId}->{"subroutinesCalled"}->{$subroutineName} = $variableName;
			  # Mark line as processed.
			  $lineProcessed = 1;
		      }
		  }

		  # Check for type-bound procedures.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      if ( $units{$unitIdList[-1]}->{"unitType"} eq "type" ) {
			  if ( $processedLine =~ m/$typeBoundRegex/i ) {
			      my $methodName = lc($2);
			      (my $procedureList = lc($3)) =~ s/\s//g;
			      my @procedures = split(/,/,$procedureList);
			      my $unitId = $unitIdList[-1];
			      @{$units{$unitId}->{"methods"}->{$methodName}} = @procedures;
			      # Mark line as processed.
			      $lineProcessed = 1;
			  }
		      }
		  }

		  # Detect intrinsic variable declarations.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      foreach my $intrinsicType ( keys(%intrinsicDeclarations) ) {
			  # Check for a match to an intrinsic declaration regex.
			  if ( my @matches = $processedLine =~ m/$intrinsicDeclarations{$intrinsicType}->{"regEx"}/i ) {
			      my $matchIndex = $intrinsicDeclarations{$intrinsicType}->{"variables"};
			      my $variablesList = lc($matches[$matchIndex-1]);
			      # Get ID of unit.
			      my $unitId = $unitIdList[-1];
			      my @variables = &Fortran::Utils::Extract_Variables($variablesList);
			      # Store the variable list.
			      push(@{$units{$unitId}->{$intrinsicType}},@variables);
			      # Mark line as processed.
			      $lineProcessed = 1;
			  }
		      }
		  }

		  # Check for derived-type declarations.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      if ( $processedLine =~ m/$derivedTypeRegex/i ) {
			  my $derivedType = $1;
			  my $variablesList = $3;
			  my @variables = &Fortran::Utils::Extract_Variables($variablesList);		    
			  my $unitId = $unitIdList[-1];
			  push(@{$units{$unitId}->{"derivedTypesUsed"}->{$derivedType}},@variables);
			  # Mark line as processed.
			  $lineProcessed = 1;
		      }
		  }

		  # Check for module procedures.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      if ( $processedLine =~ m/$moduleProcedureRegex/i ) {
			  my $procedureName = lc($1);
			  my $unitId = $unitIdList[-1];
			  $units{$unitId}->{"moduleProcedures"}->{$procedureName} = 1;
			  # Mark line as processed.
			  $lineProcessed = 1;
		      }
		  }

		  # Check for function calls.
		  unless ( $lineProcessed || $fileStack[0]->{'inXML'} || $fileStack[0]->{'inLaTeX'} ) {
		      my $functionSeek = lc($processedLine);
		      $functionSeek =~ s/''//g;
		      $functionSeek =~ s/""//g;
		      $functionSeek =~ s/'[^']+'//g;
		      $functionSeek =~ s/"[^"]+"//g;
		      my $countIterations = 0;
		      while ( $functionSeek =~ m/\(/ ) {
			  ++$countIterations;
			  die("Code_Analyzer.pl: exceeded 1000 iterations attempting to extract function call - probably indicates a failure")
			      if ( $countIterations > 1000 );
			  (my $extracted, my $remainder, my $prefix) = extract_bracketed($functionSeek,"()","[^\\(]+");
			  unless ( defined($prefix) ) {
			      print "failed to find function name:\n";
			      print "\traw line: ".$rawLine;
			      print "\tprocessed line: ".$processedLine;
			      print "\tbuffered comments: ".$bufferedComments;
			      print "\tremaining line: ".$functionSeek;
			      exit;
			  }
			  if ( $prefix =~ m/(([a-z0-9_]+)\s*\%\s*)([a-z0-9_]+)$/ ) {
			      my $functionName = lc($3);
			      (my $derivedType = lc($1)) =~ s/%\s*$//;
			      if ( $derivedType eq "" ) {$derivedType = -1};
			      my $unitId = $unitIdList[-1];
			      $units{$unitId}->{"functionCalls"}->{$functionName} = $derivedType;
			  }
			  $functionSeek = $remainder;
		      }
		  }
		  
		  # Detect entering LaTeX and XML blocks.
		  $fileStack[0]->{'inXML'  } = 1
		      if ( $rawLine =~ m/^\s*!!\[/ );
		  $fileStack[0]->{'inLaTeX'} = 1
		      if ( ! $lineProcessed && $rawLine =~ m/^\s*!!\{/ );
	      }
		
	      # Close the file and shift the file stack.
	      shift(@fileStack)
		  if ( eof($fileHandle) );
	      close($fileHandle);

	    }
	}
    }
}

sub Build_Modules_Hash {
    # Create a hash giving unit IDs as a function of module names.
    foreach my $unitID ( keys(%units) ) {
	if ( $units{$unitID}->{"unitType"} eq "module" ) {
	    $modules{$units{$unitID}->{"unitName"}} = $unitID;
	}
    }
    # Create records of which units use each module.
    foreach my $unitID ( keys(%units) ) {
	if ( exists($units{$unitID}->{"modulesUsed"}) ) {
	    foreach my $moduleName ( keys(%{$units{$unitID}->{"modulesUsed"}}) ) {
		if ( exists($modules{$moduleName}) ) {
		    $units{$modules{$moduleName}}->{"usedBy"}->{$unitID} = 1;
		}
	    }
	}
    }
}

sub Output_Data {

    # Output the output file and begin a section.
    open(my $outputHandle,">".$outputFile);
    print $outputHandle "\\section{Program units}\n";

    # Get a list of all units.
    my @unsortedIDs = keys(%units);
    # Get a sorted list of all units.
    my @sortedIDs = sort(@unsortedIDs);
    # Loop over all units.
    foreach my $unitID ( @sortedIDs ) {
	# Get unit type.
	my $unitType = $units{$unitID}->{"unitType"};
	# Get unit name in LaTeX encoding.
	my $unitName = &latex_encode($units{$unitID}->{"unitName"});
	$unitName =~ s/\\_/\\\-\\_/g;
	$unitName =~ s/¦/\\textbrokenbar/g;
	# Get ID of parent.
	my $parentID;
	if ( exists($units{$unitID}->{"belongsTo"}) ) {
	    $parentID = $units{$unitID}->{"belongsTo"};
	}
	# Check if this is an abstract interface.
	my $unitIsAbstract;
	if ( $unitType eq "interface" && $unitName eq "" ) {
	    $unitIsAbstract = 1;
	} else {
	    $unitIsAbstract = 0;
	}
	# Check if parent is an abstract interface.
	my $parentUnitIsAbstract = 0;
	if ( defined($parentID) ) {
	    if ( $units{$parentID}->{"unitType"} eq "interface" && $units{$parentID}->{"unitName"} eq "" ) {
		$parentUnitIsAbstract = 1;
	    }
	}
	
	# Output header line for unit, skipping nameless interfaces (i.e. abstract interfaces) and any children of such interfaces.
	unless ( $unitIsAbstract == 1 || $parentUnitIsAbstract == 1 ) {
	    (my $hyperdefTarget = $unitID) =~ s/\./_/g;
	    print $outputHandle "\\noindent{\\normalfont \\bfseries ".$unitType.":} \\hypertarget{".$unitID."}{{\\normalfont \\ttfamily ".$unitName."}}\\hyperdef{source}{".$hyperdefTarget."}{}\\index[code]{".$unitName."\@\{\\normalfont \\ttfamily ".$unitName."} (".$unitType.")}\n\n";

	    # Begin tabulated output.
	    my $tableOpen = "\\begin{supertabular}{lp{70mm}p{70mm}}\n";
	    my $tableIsOpen = 0;
	    
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
	    if ( exists($units{$unitID}->{"codeLines"}) ) {
		print $outputHandle "\\emph{Code lines:} & \\multicolumn{2}{l}{".$units{$unitID}->{"codeLines"}."} \\\\\n";
	    }

	    # Output any parent.
	    if  ( defined($parentID) ) {
		if ( $tableIsOpen == 0 ) {
		    print $outputHandle $tableOpen;
		    $tableIsOpen = 1;
		}
		print $outputHandle "\\emph{Contained by:} & \\multicolumn{2}{l}{".$units{$parentID}->{"unitType"}." \\hyperlink{".$parentID."}{{\\normalfont \\ttfamily ".&latex_encode($units{$parentID}->{"unitName"})."}}} \\\\ \n" 
	    }
	    
	    # Output lists of modules used.
	    if ( exists($units{$unitID}->{"modulesUsed"}) ) {
		if ( $tableIsOpen == 0 ) {
		    print $outputHandle $tableOpen;
		    $tableIsOpen = 1;
		}
		print $outputHandle "\\emph{Modules used:} ";
		my @unsortedModules = keys(%{$units{$unitID}->{"modulesUsed"}});
		my @sortedModules = sort(@unsortedModules);
		foreach ( @sortedModules ) {
		    if ( exists($modules{$_}) ) {
			if ( $modules{$_} eq "" ) {
			    $_ = "{\\normalfont \\ttfamily ".&latex_encode($_)."}";
			} else {
			    $_ = "\\hyperlink{".$modules{$_}."}{{\\normalfont \\ttfamily ".&latex_encode($_)."}}";
			}
		    } else {
			$_ = "{\\normalfont \\ttfamily ".&latex_encode($_)."}";
		    }
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
		my @unsortedUnits = keys(%{$units{$unitID}->{"usedBy"}});
		my @sortedUnits = sort(@unsortedUnits);
		foreach ( @sortedUnits ) {
		    if ( $_ eq "" ) {
			$_ = $units{$_}->{"unitType"}." {\\normalfont \\ttfamily ".&latex_encode($units{$_}->{"unitName"})."}";
		    } else {
			$_ = $units{$_}->{"unitType"}." \\hyperlink{".$_."}{{\\normalfont \\ttfamily ".&latex_encode($units{$_}->{"unitName"})."}}";
		    }
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
    my $outputHandle = $_[0];
    my @toPrint = @{$_[1]};
    my $iColumn = 0;
    foreach my $item ( @toPrint ) {
	print $outputHandle " & \\RaggedRight ".$item;
	print $outputHandle "\\\\\n" if ( $iColumn == 1);
	$iColumn = 1-$iColumn;
    }
    print $outputHandle " & \\\\\n" if ( $iColumn == 1 );
}

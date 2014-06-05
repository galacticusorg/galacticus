# Contains a Perl module which implements various useful functionality for handling Fortran source code.

package Fortran_Utils;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use File::Copy;
use Text::Balanced qw (extract_bracketed);
use Text::Table;
use Switch;
use Data::Dumper;
use Fcntl qw(SEEK_SET);

# RegEx's useful for matching Fortran code.
our $classDeclarationRegEx = qr/^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*([a-z0-9_]+)\s*$/i;
our $variableDeclarationRegEx = qr/^\s*(?i)(integer|real|double precision|logical|character|type|class|complex)(?-i)\s*(\(\s*[a-zA-Z0-9_=]+\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$/;

sub Truncate_Fortran_Lines {
    # Scans a Fortran file and truncates source lines to be less than 132 characters in length as (still) required by some compilers.
    # Includes intelligent handling of OpenMP directives.

    # Get file name.
    my $inFile = $_[0];

    # Make a backup copy.
    copy($inFile,$inFile."~");

    # Specify input and output file names.
    my $outFile = $inFile;
    $inFile     = $inFile."~";

    # Specify maximum line length (short enough to allow inclusion of continuation character).
    my $lineLengthMaximum = 128;
    
    # Open input and output files.
    open(inHandle ,    $inFile );
    open(outHandle,">".$outFile);
    # Loop through the input file.
    while ( ! eof(inHandle) ) {
	my $truncateDo = 0;    # Set to not truncate by default.
	my $comments;          # Clear comments buffer.
	my $indent;            # Clear indentation string.
	my $line = <inHandle>; # Get the next line.
	my $buffer = $line;    # Put into the buffer.
	chomp($line);          # Remove newline.
	my $lineLength = length($line);
	if ( $line =~ s/^(\s*)//       ) {$indent   = $1};
	if ( $line =~ s/(\![^\$].*)$// ) {$comments = $1};         # Remove comments and store in comments buffer.
	(my $linesJoined = $line) =~ s/&\s*$//;                    # Put line into joined lines buffer.
	if ( $lineLength > $lineLengthMaximum ) {$truncateDo = 1}; # Check if line is overlong and flag if necessary.

	# Test for OpenMP directives.
	my $isOpenMpConditional = 0;
	my $isOpenMpDirective   = 0;
	if ( $line =~ m/^\s*\!\$\s/  ) {$isOpenMpConditional = 1}; # Line is an OpenMP conditional.
	if ( $line =~ m/^\s*\!\$omp/ ) {$isOpenMpDirective   = 1}; # Line is an OpenMP directive.
	my $endedOnPreprocessorDirective = 0;
	while ( $line =~ m/&\s*$/ ) { # While continuation lines are present, read next line and process in same way.
	    my $inFileStoredPosition = tell(inHandle);
	    $line = <inHandle>;
	    if ( $line =~ m/^\#/ ) {
		seek(inHandle,$inFileStoredPosition,0);
		$endedOnPreprocessorDirective = 1;
		last;
	    }
	    $buffer .= $line;
	    chomp($line);
	    $lineLength = length($line); 
	    if ( $line =~ s/^(\s*)//  ) {$indent    = $1};
	    if ( $line =~ s/(\!.*)$// ) {$comments .= $1};
	    (my $lineTemporary  = $line)          =~ s/^\s*&//;
	    ($linesJoined      .= $lineTemporary) =~ s/&\s*$//;
	    if ( $lineLength > $lineLengthMaximum ) {$truncateDo = 1};
	}
	# Check if we need to truncate.
	unless ( $truncateDo == 0 ) {
	    # Yes we do - split the line intelligently.
	    undef($buffer); # Clear the buffer.
	    my $linePrefix = "";
	    my $continuationSuffix = "&";
	    my $continuationPrefix = "&";
	    my $targetLineLength = $lineLengthMaximum-length($indent)-length($continuationSuffix);
	    # Check if we need to add directive prefixes.
	    if ( $isOpenMpConditional == 1 ) {
		$linesJoined =~ s/\s*\!\$\s*//g;
		$targetLineLength -= 3;
		$linePrefix = "!\$ ";
	    }
	    if ( $isOpenMpDirective == 1 ) {
		$linesJoined =~ s/\s*\!\$omp\s*//g;
		$targetLineLength -= 6;
		$linePrefix = "!\$omp ";
		$continuationSuffix = "&";
		$continuationPrefix  = "";
	    }
	    my $lineJoiner = "";
	    while ( length($linesJoined) > $targetLineLength ) {
		# Find a point to break the line.
		my $breakPoint = $targetLineLength;
		while (
		       ( substr($linesJoined,$breakPoint,1) !~ m/[\s\+\-\*\/\(\)\,\%]/
			 || substr($linesJoined,0,$breakPoint+1) =~ m/e\s*[\+\-]$/
			 || substr($linesJoined,0,$breakPoint+1) =~ m/\*\*$/
			 || substr($linesJoined,0,$breakPoint+1) =~ m/\/\/$/
			 )
		       &&
		       $breakPoint >= 0
		       ) {
		    # Ensure that we jump over "e+" or "e-" exponents.
		    if ( substr($linesJoined,0,$breakPoint+1) =~ m/e\s*[\+\-]$/ ) {--$breakPoint};
		    --$breakPoint;
		}
		if ( $breakPoint <= 0 ) { # Didn't find a suitable breakpoint.
		    copy($inFile,$outFile); # Restore the original file.
		    die "Truncate_Fortran_Lines.pl: failed to find a break point in line (original file - ".$inFile." - restored)";
		}
		# Add the pre-break line to the buffer.
		$buffer .= $indent.$linePrefix.$lineJoiner.substr($linesJoined,0,$breakPoint);
		if ( $linesJoined !~ m/^\s*$/ ) {$buffer .= $continuationSuffix};
		$buffer .= "\n";
		# Remove pre-break line from the joined line buffer.
		$linesJoined = substr($linesJoined,$breakPoint,length($linesJoined)-$breakPoint);
		# After first pass, set the joining string to the continutation string.
		$lineJoiner = $continuationPrefix;
	    }
	    # Append the remaining text of the joined buffer and any comments.
	    $buffer .= $indent.$linePrefix.$lineJoiner.$linesJoined;
	    if ( $endedOnPreprocessorDirective == 1 ) {$buffer .= $continuationSuffix};
	    $buffer .= $comments
		if ( defined($comments) );
	    $buffer .= "\n";
	}
	print outHandle $buffer;
    }
    close(inHandle );
    close(outHandle);

}

sub Get_Matching_Lines {
    # Return a list of all lines in a file matching a supplied regular expression.
    my $fileName = shift;
    my $regEx    = shift;
    # Open the file, and read each line.
    my @matches;
    open(my $fileHandle,$fileName);
    until ( eof($fileHandle) ) {
	&Get_Fortran_Line($fileHandle,my $rawLine,my $processedLine,my $bufferedComments);
	if ( my @submatches = $processedLine =~ $regEx ) {
	    push(
		@matches,
		{
		    line => $processedLine,
		    submatches => \@submatches
		}
		);
	}
    }
    close($fileName);
    return @matches;
}

sub read_file {
    # Return a complete listing of a source file, optionally returning only the processed lines.
    my $fileName = shift;
    (my %options) = @_
	if ( scalar(@_) > 1 );
    # Determine what to return.
    my $returnType = "raw";
    $returnType = $options{'state'}
       if ( exists($options{'state'}) );
    # Determine whitespace stripping options.
    my $stripEmptyLines    = 0;
    my $stripLeadingSpace  = 0;
    my $stripTrailingSpace = 0;
    my $stripRegEx;
    $stripEmptyLines    = $options{'stripEmpty'   }
       if ( exists($options{'stripEmpty'   }) );
    $stripLeadingSpace  = $options{'stripLeading' }
       if ( exists($options{'stripLeading' }) );
    $stripTrailingSpace = $options{'stripTrailing'}
       if ( exists($options{'stripTrailing'}) );
    $stripRegEx         = $options{'stripRegEx'   }
       if ( exists($options{'stripRegEx'   }) );
    # Determine whether or not to follow included files.
    my $followIncludes = 0;
    $followIncludes = $options{'followIncludes'}
       if ( exists($options{'followIncludes'}) );
    my @includeLocations = ( "" );
    push(@includeLocations,@{$options{'includeLocations'}})
	if ( exists($options{'includeLocations'}) );
    # Initialize the file name stack.
    my @fileNames     = ( $fileName );
    my @filePositions = (        -1 );  
    # Initialize the code buffer.
    my $codeBuffer;
    # Iterate until all files are processed.
    while ( scalar(@fileNames) > 0 ) {
	# Open the file.
	open(my $fileHandle,$fileNames[0]);
	seek($fileHandle,$filePositions[0],SEEK_SET) 
	    unless ( $filePositions[0] == -1 );
	until ( eof($fileHandle) ) {
	    &Get_Fortran_Line($fileHandle,my $rawLine,my $processedLine,my $bufferedComments);
	    # Detect include files, and recurse into them.
	    if ( $followIncludes == 1 && $processedLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
		my $includeFileLeaf  = $1;
		my $includeFileFound = 0;
		foreach my $suffix ( ".inc", ".Inc" ) {
		    foreach my $includeLocation ( @includeLocations ) {
			(my $includePath = $fileNames[0]) =~ s/\/[^\/]+$/\//;
			(my $includeFile = $includePath.$includeLocation."/".$includeFileLeaf) =~ s/\.inc$/$suffix/;
			if ( -e $includeFile ) {
			    $filePositions[0] = tell($fileHandle);
			    unshift(@fileNames,$includeFile);
			    unshift(@filePositions,-1);
			    $includeFileFound = 1;
			    last;
			}
		    }
		    last
			if ( $includeFileFound == 1 );
		}
		last
		    if ( $includeFileFound == 1 );
	    }
	    # Process the line.
	    my $line;
	    switch ( $returnType ) {
		case ( "raw"       ) {$line = $rawLine         }
		case ( "processed" ) {$line = $processedLine   }
		case ( "comments"  ) {$line = $bufferedComments}
	    }
	    $line =~ s/$stripRegEx//
		if ( defined($stripRegEx) );
	    $line =~ s/^\s*//
		if ( $stripLeadingSpace  == 1 );
	    $line =~ s/\s*$//g 
		if ( $stripTrailingSpace == 1 );
	    chomp($line);	
	    $codeBuffer .= $line."\n"
		unless ( $line =~ m/^\s*$/ && $stripEmptyLines == 1 );
	}
	# Close the file and shift the list of filenames.
	if ( eof($fileHandle) ) {
	    shift(@fileNames    );
	    shift(@filePositions);
	}
	close($fileHandle);
    }
    return $codeBuffer;
}

sub Get_Fortran_Line {
    # Reads a line (including all continuation lines) from a Fortran source. Returns the line, a version with comments stripped
    # and all continuations combined into a single line and a buffer of removed comments.

    # Get file handle.
    my $inHndl = $_[0];

    my $processedFullLine = 0;
    my $rawLine           = "";
    my $processedLine     = "";
    my $bufferedComments  = "";
    my $firstLine         = 1;
    while ( $processedFullLine == 0 ) {
	# Get a line;
	my $line = <$inHndl>;
	# Strip comments and grab any continuation lines.
	my $tmpLine = $line;
	my $inDoubleQuotes  =  0;
	my $inSingleQuotes  =  0;
	my $commentPosition = -1;
	for(my $iChar=0;$iChar<length($tmpLine);++$iChar) {
	    my $char = substr($tmpLine,$iChar,1);
	    if ( $char eq "'" ) {
		if ( $inDoubleQuotes == 0 ) {
		    $inDoubleQuotes = 1;
		} else {
		    $inDoubleQuotes = 0;
		}
	    }
	    if ( $char eq '"' ) {
		if ( $inSingleQuotes == 0 ) {
		    $inSingleQuotes = 1;
		} else {
		    $inSingleQuotes = 0;
		}
	    }
	    if ( $commentPosition == -1 && $char eq "!" && $inDoubleQuotes == 0 && $inSingleQuotes == 0 ) {$commentPosition = $iChar};
	}
	$rawLine .= $line;
	chomp($processedLine);
	if ( $commentPosition == -1 ) {
	    $tmpLine = $line;
	} else {
	    $tmpLine = substr($line,0,$commentPosition)."\n";	       
	    $bufferedComments .= substr($line,$commentPosition+1,length($line)-$commentPosition)."\n";
	    chomp($bufferedComments);
	}
	$tmpLine =~ s/^\s*&\s*//;
	$processedLine .= " " unless ( $processedLine eq "" );
	$processedLine .= $tmpLine;
	if ( $processedLine =~ m/&\s*$/ ) {
	    $processedLine =~ s/\s*&\s*$//;
	} elsif ( $firstLine == 0 && $line =~ m/^\#/ ) {
	    # This is a preprocessor directive in the middle of continuation lines. Just concatenate it.
	} else {
	    $processedFullLine = 1;
	}
	$firstLine = 0;
    }

    # Return date.
    $_[1] = $rawLine;
    $_[2] = $processedLine;
    $_[3] = $bufferedComments;
}

sub Format_Variable_Defintions {
    # Generate formatted variable definitions.
    my $variables = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};
    $options{'variableWidth'} = -1
	unless ( exists($options{'variableWidth'}) );
    $options{'alignVariables'} = 0
	unless ( exists($options{'alignVariables'}) );
    # Record of inline comments.
    my @inlineComments;
    # Scan data content searching for repeated attributes.
    my %attributes;
    foreach my $datum ( @{$variables} ) {
	if ( exists($datum->{'attributes'}) ) {
	    foreach ( @{$datum->{'attributes'}} ) {
		$_ =~ s/intent\(\s*in\s*\)/intent(in   )/;
		$_ =~ s/intent\(\s*out\s*\)/intent(  out)/;
		$_ =~ s/intent\(\s*inout\s*\)/intent(inout)/;
		(my $attributeName = $_) =~ s/^([^\(]+).*/$1/;
		++$attributes{$attributeName}->{'count'};
		$attributes{$attributeName}->{'column'} = -1;
	    }
	}
    }
    # Find column for aligned attributes.
    my $columnCountMaximum = -1;
    foreach my $datum ( @{$variables} ) {
	if ( exists($datum->{'attributes'}) ) {
	    my $columnCount  = -1;
	    foreach ( sort(@{$datum->{'attributes'}}) ) {
		(my $attributeName = $_) =~ s/^([^\(]+).*/$1/;
		++$columnCount;
		if ( $attributes{$attributeName}->{'count'} > 1 ) {
		    if ( $columnCount > $attributes{$attributeName}->{'column'} ) {
			foreach my $otherAttribute ( sort(keys(%attributes)) ) {
			    ++$attributes{$otherAttribute}->{'column'}
			    if (
				$attributes{$otherAttribute}->{'column'} >= $columnCount &&
				$attributes{$otherAttribute}->{'count' } > 1             &&
				$otherAttribute ne $attributeName
				);
			}
			$attributes{$attributeName}->{'column'} = $columnCount;
		    }
		    $columnCount = $attributes{$attributeName}->{'column'};
		}
	    }
	    $columnCountMaximum = $columnCount+1
		if ( $columnCount+1 > $columnCountMaximum );
	}
    }
    foreach ( keys(%attributes) ) {
	$columnCountMaximum = $attributes{$_}->{'column'}
	if ( $attributes{$_}->{'column'} > $columnCountMaximum);
    }
    ++$columnCountMaximum;
    my @attributeColumns;
    push @attributeColumns, {is_sep => 1, body => ""},{align => "left"} foreach (1..$columnCountMaximum);

    # Find the number of columns to use for variables.
    my @variablesColumns;
    my $variableColumnCount;
    if ( $options{'alignVariables'} == 1 ) {
	my $variableLengthMaximum = 0;
	my $variableColumnCountMaximum    = 0;
	foreach my $definition ( @{$variables} ) {
	    $variableColumnCountMaximum = scalar(@{$definition->{'variables'}})
		if ( scalar(@{$definition->{'variables'}}) > $variableColumnCountMaximum );
	    foreach ( @{$definition->{'variables'}} ) {
		$variableLengthMaximum = length($_)
		    if ( length($_) > $variableLengthMaximum );
	    }
	}
	if ( $options{'variableWidth'} > 0 ) {
	    $variableColumnCount = int($options{'variableWidth'}/($variableLengthMaximum+2));
	    $variableColumnCount = 2
		if ( $variableColumnCount < 2 );
	} else {
	    $variableColumnCount = $variableColumnCountMaximum;
	}
	push @variablesColumns, {is_sep => 1, body => ""}, {align => "left"} foreach (1..5*$variableColumnCount);
    } else {
	push @variablesColumns, {is_sep => 1, body => " :: "}, {align => "left"};
    }

    # Construct indentation.
    my $indent = "    ";
    $indent = " " x $options{'indent'}
    if ( exists($options{'indent'}) );
    # Create a table for data content.
    my @columnsDef = 
	(
	 {
	     is_sep => 1,
	     body   => $indent
	 },
	 {
	     align  => "left"
	 },
	 {
	     is_sep => 1,
	     body   => ""
	 },
	 {
	     align  => "left"
	 },
	 {
	     is_sep => 1,
	     body   => ""
	 },
	 {
	     align  => "left"
	 },
	 {
	     is_sep => 1,
	     body   => ""
	 },
	 {
	     align  => "left"
	 },
	 @attributeColumns
	);
    if ( $options{'alignVariables'} == 1 ) {
	push(
	    @columnsDef,   
	    {
		align  => "left"
	    },
	    @variablesColumns,
	    {
		is_sep => 1,
		body   => ""
	    },
	    {
		align  => "left"
	    },
	    {
		align  => "left"
	    }
	    );
    } else {
	push(
	    @columnsDef,   
	    @variablesColumns,
	    {
		align  => "left"
	    }
	    );
    }
    my $dataTable = Text::Table->new(@columnsDef);

    # Iterate over all data content.
    foreach ( @{$variables} ) {
	# Construct the type definition.
	my @typeDefinition = ( "", "", "" );
	@typeDefinition = ( "(", $_->{'type'}, ")" )
	    if ( exists($_->{'type'}) );
	# Add attributes.
	my @attributeList;
	if ( exists($_->{'attributes'}) ) {
	    foreach ( sort(@{$_->{'attributes'}}) ) {
		(my $attributeName = $_) =~ s/^([^\(]+).*/$1/;
		if ( $attributes{$attributeName}->{'column'} >= 0 ) {
		    push(@attributeList,"")
			while ( scalar(@attributeList) < $attributes{$attributeName}->{'column'} );
		}
		push(@attributeList,", ".$_);
	    }
	    push(@attributeList,"")
		while ( scalar(@attributeList) < $columnCountMaximum);
	} else {
	    @attributeList = ("") x $columnCountMaximum;
	}
	# Construct any comment.
	my $comment = "";
	$comment = " ! ".$_->{'comment'}
	    if ( exists($_->{'comment'}) );
	# Add a row to the table.
	if ( $options{'alignVariables'} == 1 ) {
	    for(my $i=0;$i<scalar(@{$_->{'variables'}});$i+=$variableColumnCount) {
		my $i0 = $i;
		my $i1 = $i+$variableColumnCount-1;
		$i1 = $#{$_->{'variables'}}
		   if ( $i1 > $#{$_->{'variables'}} );
		my @thisVariables = @{$_->{'variables'}}[$i0..$i1];
		my @splitVariables;
		for(my $j=0;$j<scalar(@thisVariables);++$j) {
		    if ( $thisVariables[$j] =~ m/([a-zA-Z0-9_]+)\s*(\([a-zA-Z0-9:,\(\)]+\))??(\s*=\s*(.*?)\s*)??$/ ) {
			my $variableName = $1;
			$variableName = ", ".$variableName
			    if ( $j > 0 );
			my $dimensions   = $2;
			my $assignment   = $4;
			$assignment = "=".$assignment
			    if ( defined($assignment) );
			my @dimensioner;
			if ( defined($dimensions) ) {
			    $dimensions =~ s/^\(//;
			    $dimensions =~ s/\)$//;
			    @dimensioner = ( "(", $dimensions, ")" );
			} else {
			    @dimensioner = ( "", "", "" );
			}
			push(@splitVariables,$variableName,@dimensioner,$assignment);
		    } else {
			print "Variable ".$thisVariables[$j]." is unmatched\n";
			die;
		    }
		}
		my $continuation = "";
		$continuation = ", &"
		    if ( $i+$variableColumnCount < scalar(@{$_->{'variables'}}) );
		if ( $i == 0 ) {
		    $dataTable->add(
			$_->{'intrinsic'},
			@typeDefinition,
			@attributeList,
			":: ",
			@splitVariables,
			$continuation,
			$comment
			);
		} else {
		    $dataTable->add(
			"     &",
			map("",@typeDefinition),
			map("",@attributeList ),
			"",
			@splitVariables,
			$continuation,
			""
			);
		}
	    }
	} else {
	    $dataTable->add(
		$_->{'intrinsic'},
		@typeDefinition,
		@attributeList,
		join(", ",@{$_->{'variables'}}),
		$comment
		);
	}
	# Add a comment after this row.
	if ( exists($_->{'commentAfter'}) ) {
	    push(@inlineComments,$_->{'commentAfter'});
	    $dataTable->add("%C".scalar(@inlineComments));
   	}
    }
    # Get the formatted table.
    my $formattedVariables = $dataTable->table();
    # Reinsert inline comments.
    for(my $i=scalar(@inlineComments);$i>=1;--$i) {	
	chomp($inlineComments[$i-1]);
	$inlineComments[$i-1] =~ s/^\s*//;
	$inlineComments[$i-1] =~ s/\s*$//;
	$formattedVariables =~ s/\%C$i\s*\n/$inlineComments[$i-1]\n/;
    }
    $formattedVariables =~ s/\%BLANKLINE//g;
    # Return the table.
    return $formattedVariables;
}

sub Extract_Variables {
    # Given the post-"::" section of a variable declaration line, return an array of all variable names.
    my $variableList = shift;
    return
	unless ( defined($variableList) );
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};
    $options{'lowerCase'} = 1
	unless ( exists($options{'lowerCase'}) );
    $options{'keepQualifiers'} = 0
	unless ( exists($options{'keepQualifiers'}) );
    $options{'removeSpaces'} = 1
	unless ( exists($options{'removeSpaces'}) );
    die("Fortran::Utils::Extract_Variables: variable list contains '::' - most likely regex matching failed") if ( $variableList =~ m/::/ );
    # Convert to lower case.
    $variableList = lc($variableList)
	if ( $options{'lowerCase'} == 1 );
    # Remove whitespace.
    if ( $options{'removeSpaces'} == 1 ) {
	$variableList =~ s/\s//g;
    } else {
	$variableList =~ s/\s*$//;
    }
    # Remove *'s (can appear for character variables).
    unless ( $options{'removeSpaces'} == 0 ) {
	if ( $options{'keepQualifiers'} == 0 ) {
	    $variableList =~ s/\*//g;
	} else {
            $variableList =~ s/\*/\%\%ASTERISK\%\%/g;
	}
    }
    # Remove text within matching () pairs.
    while ( $variableList =~ m/[\(\[]/ ) {
	(my $extracted, my $remainder, my $prefix) = extract_bracketed($variableList,"()[]","[\\sa-zA-Z0-9_,:=>\%\\+\\-\\*\\/\.]+");
	if ( $options{'keepQualifiers'} == 0 ) {
	    $variableList = $prefix.$remainder;
	} else {
	    $extracted =~ s/\(/\%\%OPEN\%\%/g;
	    $extracted =~ s/\)/\%\%CLOSE\%\%/g;
	    $extracted =~ s/\[/\%\%OPENSQ\%\%/g;
	    $extracted =~ s/\]/\%\%CLOSESQ\%\%/g;
	    $extracted =~ s/,/\%\%COMMA\%\%/g;
	    $variableList = $prefix.$extracted.$remainder;
	}
    }
    # Remove any definitions or associations.
    $variableList =~ s/=[^,]*(,|$)//g
	if ( $options{'keepQualifiers'} == 0 );
    # Split variables into an array and store.
    my @variables = split(/\s*,\s*/,$variableList);
    if ( $options{'keepQualifiers'} == 1 ) {
	foreach ( @variables ) {
	    $_ =~ s/\%\%OPEN\%\%/\(/g;
	    $_ =~ s/\%\%CLOSE\%\%/\)/g;
	    $_ =~ s/\%\%OPENSQ\%\%/\[/g;
	    $_ =~ s/\%\%CLOSESQ\%\%/\]/g;
	    $_ =~ s/\%\%COMMA\%\%/,/g;
	    $_ =~ s/\%\%ASTERISK\%\%/\*/g;
	}
    }
    # Return the list.
    return @variables;
}

1;

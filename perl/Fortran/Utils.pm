# Contains a Perl module which implements various useful functionality for handling Fortran source code.

package Fortran::Utils;
use strict;
use warnings;
use File::Copy;
use Text::Balanced qw (extract_bracketed);
use Text::Table;
use Data::Dumper;
use List::Uniq qw(uniq);
use Fcntl qw(SEEK_SET);

# RegEx's useful for matching Fortran code.
our $label = qr/[a-zA-Z0-9_\{\}¦]+/;
our $argumentList = qr/[a-zA-Z0-9_\{\}¦,\s]*/;
our $classDeclarationRegEx = qr/^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\((${label})\)\s*)*(::)??\s*([a-z0-9_]+)\s*$/i;
our $variableDeclarationRegEx = qr/^\s*(!\$\s*)??(?i)(integer|real|double precision|logical|character|type|class|complex|procedure|generic)(?-i)\s*(\(\s*[a-zA-Z0-9_=\*]+\s*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$/;

# Specify unit opening regexs.
our %unitOpeners = (
    # Find module openings, avoiding module procedures, functions, and subroutines.
    module             => { unitName => 0                                           , regEx => qr/^\s*module\s+(${label})\s*$/ },
    # Find submodule openings.
    submodule          => { unitName => 0                                           , regEx => qr/^\s*submodule\s+\(\s*[a-zA-Z0-9_:\{\}¦]+\s*\)\s+(${label})\s*$/ },
    # Find program openings.
    program            => { unitName => 0                                           , regEx => qr/^\s*program\s+(${label})/ },
    # Find subroutine openings, allowing for pure, elemental and recursive subroutines.
    subroutine         => { unitName => 1                           , arguments => 3, regEx => qr/^\s*(impure\s+|pure\s+|elemental\s+|recursive\s+|module\s+)*\s*subroutine\s+(${label})\s*(\(\s*(${argumentList})\))*/},
    # Find function openings, allowing for pure, elemental, and recursive functions, and different function types.
    function           => { unitName => 5, intrinsic => 1, kind => 2, arguments => 7, regEx => qr/^\s*(impure\s+|pure\s+|elemental\s+|recursive\s+|module\s+)*\s*(real|integer|double\s+precision|double\s+complex|character|logical)*\s*(\(((kind|len)=)??[\w\d]*\))*\s*function\s+(${label})\s*(\(\s*(${argumentList})\))*/},
     # Find submodule module procedure openings.
    moduleProcedure    => { unitName => 0                                           , regEx => qr/^\s*module\s+procedure\s+(${label})/},
    # Find interfaces.
    interface          => { unitName => 1                                           , regEx => qr/^\s*(abstract\s+)??interface\s+([a-zA-Z0-9_\(\)\/\+\-\*\.=]*)/},
    # Find types.
    type               => { unitName => 2                                           , regEx => qr/^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(${label}\)\s*)*(::)??\s*(${label})\s*$/}
    );

# Specify unit closing regexs.
our %unitClosers = (
    module             => { unitName => 0, regEx => qr/^\s*end\s+module\s+(${label})/ },
    submodule          => { unitName => 0, regEx => qr/^\s*end\s+submodule\s+(${label})/ },
    program            => { unitName => 0, regEx => qr/^\s*end\s+program\s+(${label})/ },
    subroutine         => { unitName => 0, regEx => qr/^\s*end\s+subroutine\s+(${label})/},
    function           => { unitName => 0, regEx => qr/^\s*end\s+function\s+(${label})/},
    moduleProcedure    => { unitName => 0, regEx => qr/^\s*end\s+procedure\s+(${label})/},
    interface          => { unitName => 0, regEx => qr/^\s*end\s+interface\s*([a-zA-Z0-9_\(\)\/\+\-\*\.=]*)/},
    type               => { unitName => 0, regEx => qr/^\s*end\s+type\s+(${label})/}
    );

# Specify regexs for intrinsic variable declarations.
our %intrinsicDeclarations = (
    integer       => { intrinsic => "integer"         , openmp => 0, type => 1, attributes => 3, variables => 4, regEx => qr/^\s*(!\$)??\s*(?i)integer(?-i)\s*(\(\s*([a-zA-Z0-9_=]+|kind\s*\(\s*[a-zA-Z0-9_]+\s*\))\s*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    real          => { intrinsic => "real"            , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)real(?-i)\s*(\(\s*[a-zA-Z0-9_=]+\s*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    double        => { intrinsic => "double precision", openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)double\s+precision(?-i)\s*(\(\s*[a-zA-Z0-9_=]+\s*\))*([\sa-zA-Z0-9_,%:=\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    complex       => { intrinsic => "complex"         , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)complex(?-i)\s*(\(\s*[a-zA-Z0-9_=]+\s*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    doubleComplex => { intrinsic => "double complex"  , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)double\s+complex(?-i)\s*(\(\s*[a-zA-Z0-9_=]+\s*\))*([\sa-zA-Z0-9_,%:=\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    logical       => { intrinsic => "logical"         , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)logical(?-i)\s*(\(\s*[a-zA-Z0-9_=]+\s*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_\.,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    character     => { intrinsic => "character"       , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)character(?-i)\s*(\(\s*[a-zA-Z0-9_=,\+\-\*\(\)]+\s*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    type          => { intrinsic => "type"            , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)type(?-i)\s*(\(\s*${label}\s*\))?([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    class         => { intrinsic => "class"           , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)class(?-i)\s*(\(\s*[a-zA-Z0-9_\*]+\s*\))?([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$/ },
    procedure     => { intrinsic => "procedure"       , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)procedure(?-i)\s*(\([a-zA-Z0-9_\s]*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)<>\.\{\}¦]+)\s*$/ },
    generic       => { intrinsic => "generic"         , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)generic(?-i)\s*(\([a-zA-Z0-9_\s]*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=\+\-\*\/\(\)<>\.\{\}¦]+)\s*$/ },
    final         => { intrinsic => "final"           , openmp => 0, type => 1, attributes => 2, variables => 3, regEx => qr/^\s*(!\$)??\s*(?i)final(?-i)\s*(\([a-zA-Z0-9_\s]*\))*([\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,]+)\s*$/ },
    );

# Hash of files which have been read and processed.
my %processedFiles;

# Database of symbols exported by modules.
my %symbolsExported;

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
    my $fileName = shift();
    my $regEx    = shift();
    # Determine if we need to read the file.
    unless ( exists($processedFiles{$fileName}) ) {
	open(my $fileHandle,$fileName);
	until ( eof($fileHandle) ) {
	    &Get_Fortran_Line($fileHandle,my $rawLine,my $processedLine,my $bufferedComments);
	    push(@{$processedFiles{$fileName}},$processedLine);
	}
	close($fileName);
    }
    # Open the file, and read each line.
    my @matches;
    if ( defined($processedFiles{$fileName}) ) {
	foreach my $processedLine ( @{$processedFiles{$fileName}} ) {
	    if ( my @submatches = $processedLine =~ $regEx ) {
		push(
		    @matches,
		    {
			line       => $processedLine,
			submatches => \@submatches
		    }
		    );
	    }
	}
    }
    return @matches;
}

sub read_file {
    # Return a complete listing of a source file, optionally returning only the processed lines.
    my $fileName = shift;
    (my %options) = @_
	if ( scalar(@_) > 1 );
    # Set default set of include files to exclude.
    @{$options{'includeFilesExcluded'}} = ()
	unless ( exists($options{'includeFilesExcluded'}) );
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
		unless ( grep {$_ eq $includeFileLeaf} @{$options{'includeFilesExcluded'}} ) {
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
		}
		last
		    if ( $includeFileFound == 1 );
	    }
	    # Process the line.
	    my $line;
	    if      ( $returnType eq "raw"       ) {
		$line = $rawLine;
	    } elsif ( $returnType eq "processed" ) {
		$line = $processedLine;
	    } elsif ( $returnType eq "comments"  ) {
		$line = $bufferedComments;		
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
	my $tmpLine         = $line;
	my $inDoubleQuotes  =  0;
	my $inSingleQuotes  =  0;
	my $inBraces        =  0;
	my $commentPosition;
	# Search for start of comment in the line, if any. Regex is fast here.
	if ( $tmpLine !~ m/!/ ) {
	    ## This regex tests for the possible existance of any comment. If no "!" is present there can be no comment, and most lines don't contain comments, so we can quickly ignored them.
	    $commentPosition = -1;
	} elsif ( $tmpLine =~ m/^([^'"\{]+?)![^\$]/ ) {
	    ## This regex checks for the first "!" prior to any of "'{, and not followed by a "$" (which indicates an OpenMP directive). If found, this must be the start of a comment.
	    $commentPosition = length($1);
	} else {
	    # General case. Comments begin with a "!" symbol, but not if it is found inside a string literal, or inside braces, "{}".
	    my $removedCount = 0;
	    while ( $tmpLine =~ m/([!'"\{])/ ) {
	    	if ( $1 eq "!" ) {
	    	    # Comment character is prior to any other, so we're done.
		    if ( $tmpLine =~ s/^([^!]*!\$)// ) {
			$removedCount += length($1);
		    } else {
			last;
		    }
	    	} elsif ( $1 eq "{" ) {
	    	    # Braces can be nested so must be counted.
	    	    my $iChar = index($tmpLine,"{");
	    	    my $depth = 0;
	    	    while ( $iChar < length($tmpLine) ) {
	    		++$depth
	    		    if ( substr($tmpLine,$iChar,1) eq "{" );
	    		--$depth
	    		    if ( substr($tmpLine,$iChar,1) eq "}" );
	    		last
	    		    if ( $depth == 0 );
	    		++$iChar;
	    	    }
	    	    $removedCount += $iChar;
	    	    $tmpLine = substr($tmpLine,$iChar);
	    	} else {
	    	    # Quotes are not nested, so can be removed by regex.
	    	    my $quote = $1;
	    	    $tmpLine =~ s/^(.*?$quote[^$quote]*$quote)//;
	    	    $removedCount += length($1);
	    	}
	    }
	    $commentPosition = index($tmpLine,"!");
	    $commentPosition += $removedCount
	    	unless ( $commentPosition == -1 );
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
	$tmpLine =~ s/^\s*!\$\s+&\s*//;
	$processedLine .= " " unless ( $processedLine eq "" );
	$processedLine .= $tmpLine;
	if ( $processedLine =~ m/&\s*$/ ) {
	    $processedLine =~ s/\s*&\s*$//;
	} elsif ( $firstLine == 0 && ( $line =~ m/^\#/ || $line =~ m/^\s*![^\$]/ ) ) {
	    # This is a preprocessor directive or comment in the middle of continuation lines. Just concatenate it.
	} else {
	    $processedFullLine = 1;
	}
	$firstLine = 0;
    }

    # Return data.
    $_[1] = $rawLine;
    $_[2] = $processedLine;
    $_[3] = $bufferedComments;
}

sub Format_Variable_Definitions {
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
		if ( $columnCount > $attributes{$attributeName}->{'column'} ) {
		    foreach my $otherAttribute ( sort(keys(%attributes)) ) {
			++$attributes{$otherAttribute}->{'column'}
			if (
			    $attributes{$otherAttribute}->{'column'} >= $columnCount &&
			    $otherAttribute ne $attributeName
			    );
		    }
		    $attributes{$attributeName}->{'column'} = $columnCount;
		}
		$columnCount = $attributes{$attributeName}->{'column'};
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
    my $dataTable       = Text::Table->new(@columnsDef);
    my $ompPrivateTable = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "  !\$omp threadprivate("
	},
	{
	    align  => "left"
	},
	{
	    is_sep  => 1,
	    body    => ")"
	}
	);
    # Iterate over all data content.
    foreach ( @{$variables} ) {
	# Construct the type definition.
	my @typeDefinition = ( "", "", "" );
	@typeDefinition = ( "(", $_->{'type'}, ")" )
	    if ( exists($_->{'type'}) && defined($_->{'type'}) );
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
		    if ( $thisVariables[$j] =~ m/(${label})\s*(\([a-zA-Z0-9:,\(\)]+\))??(\s*=\s*(.*?)\s*)??$/ ) {
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
			(exists($_->{'openMP'}) && $_->{'openMP'} ? "!\$ " : "").$_->{'intrinsic'},
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
		(exists($_->{'openMP'}) && $_->{'openMP'} ? "!\$ " : "").$_->{'intrinsic'},
		@typeDefinition,
		@attributeList,
		join(", ",@{$_->{'variables'}}),
		$comment
		);
	}
	# Add any OpenMP threadprivate statement.
	$ompPrivateTable->add(join(", ",@{$_->{'variables'}}))
	    if ( exists($_->{'ompPrivate'}) && $_->{'ompPrivate'} );
	# Add a comment after this row.
	if ( exists($_->{'commentAfter'}) ) {
	    push(@inlineComments,$_->{'commentAfter'});
	    $dataTable->add("%C".scalar(@inlineComments));
   	}
    }
    # Get the formatted table.
    my $formattedVariables = 
	$dataTable      ->table().
	$ompPrivateTable->table();    
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

sub Unformat_Variables {
    # Given a Fortran-formatted variable string, decode it and return a standard variable structure.
    my $variableString = shift();
    # Iterate over intrinsic declaration regexes.
    foreach my $intrinsicType ( keys(%intrinsicDeclarations) ) {
	# Check for a match to an intrinsic declaration regex.
	if ( my @matches = $variableString =~ m/$intrinsicDeclarations{$intrinsicType}->{"regEx"}/i ) {
	    my $type               = $matches[$intrinsicDeclarations{$intrinsicType}->{"type"      }];
	    my $variablesString    = $matches[$intrinsicDeclarations{$intrinsicType}->{"variables" }];
	    my $attributesString   = $matches[$intrinsicDeclarations{$intrinsicType}->{"attributes"}];
	    $type                  =~ s/^\((.*)\)$/$1/
		if ( defined($type            ) );
	    $type                  =~ s/\s//g
		if ( defined($type            ) );
	    $attributesString      =~ s/^\s*,\s*//
		if ( defined($attributesString) );	    
	    my @variables          =  &Extract_Variables($variablesString ,keepQualifiers => 1,removeSpaces => 1);
	    my @variableNames      =  &Extract_Variables($variablesString ,keepQualifiers => 0,removeSpaces => 1);
	    my @attributes         =  &Extract_Variables($attributesString,keepQualifiers => 1,removeSpaces => 1);
	    my $variableDefinition =
	    {
		intrinsic     => $intrinsicDeclarations{$intrinsicType}->{'intrinsic'},
		variables     => \@variables                                          ,
		variableNames => \@variableNames
	    };
	    $variableDefinition->{'type'      } = $type
		if ( defined($type      )     );
	    $variableDefinition->{'attributes'} = \@attributes
		if ( scalar (@attributes) > 0 );
	    return $variableDefinition;
	}
    }
    return undef();
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
    die("Fortran::Utils::Extract_Variables: variable list '".$variableList."' contains '::' - most likely regex matching failed")
	if ( $variableList =~ m/::/ );
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
    my $iteration = 0;
    while ( $variableList =~ m/[\(\[]/ ) {
	++$iteration;
	if ( $iteration > 10000 ) {
	    print "Fortran::Utils::Extract_Variables(): maximum iterations exceeded for input:\n";
	    print " --> '".$variableList."'\n";
	    exit 1;
	}
	(my $extracted, my $remainder, my $prefix) = extract_bracketed($variableList,"()[]","[\\sa-zA-Z0-9_,:=>\%\\+\\-\\*\\/\.]+");
	if ( $options{'keepQualifiers'} == 0 ) {
	    die('Extract_Variables: failed to find prefix in "'.$variableList.'"')
		unless ( defined($prefix) );
	    $variableList = $prefix.$remainder;
	} else {
	    die("failed to remove bracketed text in '".$variableList."'")
		unless ( defined($extracted) );
	    $extracted =~ s/\(/\%\%OPEN\%\%/g;
	    $extracted =~ s/\)/\%\%CLOSE\%\%/g;
	    $extracted =~ s/\[/\%\%OPENSQ\%\%/g;
	    $extracted =~ s/\]/\%\%CLOSESQ\%\%/g;
	    $extracted =~ s/,/\%\%COMMA\%\%/g;
	    $variableList = $prefix.$extracted.$remainder;
	}
    }
    # Remove any definitions or associations.
    $variableList =~ s/=[^,]*(,|$)/$1/g
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

sub moduleSymbols {
    # Return a list of symbols exported by a named module.
    my $moduleName = shift();
    # Build list of paths that we can search for modules.
    my @moduleSearchPaths = ( $ENV{'BUILDPATH'} );
    if ( exists($ENV{'GALACTICUS_FCFLAGS'}) ) {
	my @options = split(" ",$ENV{'GALACTICUS_FCFLAGS'});
	while ( scalar(@options) > 0 ) {
	    my $option = shift(@options);
	    push(@moduleSearchPaths,shift(@options))
		if ( $option eq "-fintrinsic-modules-path" );
	}
    }
    # Parse the file if we have not yet done so.
    unless ( exists($symbolsExported{$moduleName}) ) {
	if ( $moduleName eq "iso_c_binding" ) {
	    @{$symbolsExported{$moduleName}} = map {lc($_)}
		(
		 "C_ASSOCIATED",
		 "C_F_POINTER",
		 "C_F_PROCPOINTER",
		 "C_FUNLOC",
		 "C_LOC",
		 "C_INT",
		 "C_SHORT",
		 "C_LONG",
		 "C_LONG_LONG",
		 "C_SIGNED_CHAR",
		 "C_SIZE_T",
		 "C_INT8_T",
		 "C_INT16_T",
		 "C_INT32_T",
		 "C_INT64_T",
		 "C_INT128_T",
		 "C_INT_LEAST8_T",
		 "C_INT_LEAST16_T",
		 "C_INT_LEAST32_T",
		 "C_INT_LEAST64_T",
		 "C_INT_LEAST128_T",
		 "C_INT_FAST8_T",
		 "C_INT_FAST16_T",
		 "C_INT_FAST32_T",
		 "C_INT_FAST64_T",
		 "C_INT_FAST128_T",
		 "C_INTMAX_T",
		 "C_INTPTR_T",
		 "C_FLOAT",
		 "C_DOUBLE",
		 "C_LONG_DOUBLE",
		 "C_FLOAT_COMPLEX",
		 "C_DOUBLE_COMPLEX",
		 "C_LONG_DOUBLE_COMPLEX",
		 "C_BOOL",
		 "C_CHAR",
		 "C_NULL_CHAR",
		 "C_ALERT",
		 "C_BACKSPACE",
		 "C_FORM_FEED",
		 "C_NEW_LINE",
		 "C_CARRIAGE_RETURN",
		 "C_HORIZONTAL_TAB",
		 "C_VERTICAL_TAB"
		);
	} elsif ( $moduleName eq "omp_lib" ) {
	    @{$symbolsExported{$moduleName}} = map {lc($_)}
	    (
	     "omp_get_active_level",
	     "omp_get_ancestor_thread_num",
	     "omp_get_dynamic",
	     "omp_get_level",
	     "omp_get_max_active_levels",
	     "omp_get_max_threads",
	     "omp_get_nested",
	     "omp_get_num_procs",
	     "omp_get_num_threads",
	     "omp_get_schedule",
	     "omp_get_team_size",
	     "omp_get_thread_limit",
	     "omp_get_thread_num",
	     "omp_in_parallel",
	     "omp_in_final",
	     "omp_set_dynamic",
	     "omp_set_max_active_levels",
	     "omp_set_nested",
	     "omp_set_num_threads",
	     "omp_set_schedule",
	     "omp_init_lock",
	     "omp_set_lock",
	     "omp_test_lock",
	     "omp_unset_lock",
	     "omp_destroy_lock",
	     "omp_init_nest_lock",
	     "omp_set_nest_lock",
	     "omp_test_nest_lock",
	     "omp_unset_nest_lock",
	     "omp_destroy_nest_lock",
	     "omp_get_wtick",
	     "omp_get_wtime"
	    );
	} else {
	    @{$symbolsExported{$moduleName}} = ();
	}
	# Locate the module.
	my $modulePath;
	foreach my $moduleSearchPath ( @moduleSearchPaths ) {
	    if ( -e $moduleSearchPath."/".$moduleName.".mod" ) {
		$modulePath = $moduleSearchPath."/".$moduleName.".mod";
		last;
	    }
	}
	if ( defined($modulePath) ) {
	    my $symbolBuffer = "";
	    open(my $module,"gunzip -c ".$modulePath."|");
	    do {
		my $char = getc($module);
		if ( $char eq "(" || $char eq ")" ) {
		    # Parentheses - and on buffer and clear.
		    if ( $symbolBuffer =~ m/^\s*(\d*\s+)??'([a-z0-9_]+)'\s+'$moduleName'/ ) {
			push(@{$symbolsExported{$moduleName}},$2);
		    }
		    $symbolBuffer = "";
		} else {
		    # Other - accumulate to buffer.
		    $char = " "
			if ( $char eq "\n" );
		    $symbolBuffer .= $char;
		}
	    } until ( eof($module) );
	    close($module);
	    # Uniqueify symbols.
	    @{$symbolsExported{$moduleName}} = uniq(sort(@{$symbolsExported{$moduleName}}));
	}
    }
    return @{$symbolsExported{$moduleName}};
}

1;

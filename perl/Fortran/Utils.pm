# Contains a Perl module which implements various useful functionality for handling Fortran source code.

package Fortran_Utils;
use File::Copy;
my $status = 1;
$status;

sub Truncate_Fortran_Lines {
    # Scans a Fortran file and truncates source lines to be less than 132 characters in length as (still) required by some compilers.
    # Includes intelligent handling of OpenMP directives.

    # Get file name.
    my $inFile = $_[0];

    # Make a backup copy.
    copy($inFile,$inFile."~");

    # Specify input and output file names.
    $outFile = $inFile;
    $inFile  = $inFile."~";

    # Specify maximum line length (short enough to allow inclusion of continuation character).
    $lineLengthMaximum = 128;

    # Open input and output files.
    open(inHandle ,    $inFile );
    open(outHandle,">".$outFile);
    # Loop through the input file.
    while ( ! eof(inHandle) ) {
	$truncateDo = 0;   # Set to not truncate by default.
	undef($comments);   # Clear comments buffer.
	undef($indent);     # Clear indentation string.
	$line = <inHandle>; # Get the next line.
	$buffer = $line;    # Put into the buffer.
	chomp($line);       # Remove newline.
	$lineLength = length($line);
	if ( $line =~ s/^(\s*)//       ) {$indent   = $1};
	if ( $line =~ s/(\![^\$].*)$// ) {$comments = $1};         # Remove comments and store in comments buffer.
	($linesJoined = $line) =~ s/&\s*$//;                       # Put line into joined lines buffer.
	if ( $lineLength > $lineLengthMaximum ) {$truncateDo = 1}; # Check if line is overlong and flag if necessary.

	# Test for OpenMP directives.
	$isOpenMpConditional = 0;
	$isOpenMpDirective   = 0;
	if ( $line =~ m/^\s*\!\$\s/  ) {$isOpenMpConditional = 1}; # Line is an OpenMP conditional.
	if ( $line =~ m/^\s*\!\$omp/ ) {$isOpenMpDirective   = 1}; # Line is an OpenMP directive.
	$endedOnPreprocessorDirective = 0;
	while ( $line =~ m/&\s*$/ ) { # While continuation lines are present, read next line and process in same way.
	    $inFileStoredPosition = tell(inHandle);
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
	    ($lineTemporary  = $line)          =~ s/^\s*&//;
	    ($linesJoined   .= $lineTemporary) =~ s/&\s*$//;
	    if ( $lineLength > $lineLengthMaximum ) {$truncateDo = 1};
	}
	# Check if we need to truncate.
	unless ( $truncateDo == 0 ) {
	    # Yes we do - split the line intelligently.
	    undef($buffer); # Clear the buffer.
	    $linePrefix = "";
	    $continuationSuffix = "&";
	    $continuationPrefix = "&";
	    $targetLineLength = $lineLengthMaximum-length($indent)-length($continuationSuffix);
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
	    $lineJoiner = "";
	    while ( length($linesJoined) > $targetLineLength ) {
		# Find a point to break the line.
		$breakPoint = $targetLineLength;
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
	    $buffer .= $comments."\n";
	}
	print outHandle $buffer;
    }
    close(inHandle );
    close(outHandle);

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
	} else {
	    $processedFullLine = 1;
	}
    }

    # Return date.
    $_[1] = $rawLine;
    $_[2] = $processedLine;
    $_[3] = $bufferedComments;
}

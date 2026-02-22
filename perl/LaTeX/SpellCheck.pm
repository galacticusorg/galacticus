package LaTeX::SpellCheck;
use XML::Simple;
use File::Temp;
use Regexp::Common;
use List::Uniq qw(uniq);

# Spell-checking functions for LaTeX.
# Andrew Benson (28-February-2023)

sub spellCheck {
    # Perform spell-checking on a text fragment.
    my $text             = shift();
    my $type             = shift();
    my $fileNameOriginal = shift();
    # Determine if this is LaTeX source.
    my $isLaTeX = $type eq "latex";
    my $suffix  = $isLaTeX ? ".tex" : ".txt";
    # Build a temporary file.
    my $tmpFile = File::Temp->new( UNLINK => 0, SUFFIX => $suffix);
    my $fileName = $tmpFile->filename();
    print $tmpFile $text;
    close($tmpFile);
    # Do the spell check.
    my $warnings = &spellCheckFile($fileName,$fileNameOriginal);
    # Clean up.
    unlink($fileName);
    return $warnings;
}

sub spellCheckFile {
    # Perform spell-checking on a file.
    # Get file name and determine if this is a LaTeX file.
    my $fileName         = shift();
    my $fileNameOriginal = shift();
    my $isLaTeX          = $fileName =~ m/\.tex$/;
    
    # Load functionClass names.
    unless ( @spellWords ) {
	if ( -e "work/build/stateStorables.xml" ) {
	    my $xml            = new XML::Simple();
	    my $stateStorables = $xml->XMLin("work/build/stateStorables.xml");
	    my @instanceNames;
	    push(@spellWords,                                                         @{$stateStorables->{'functionClassInstances'}} );
	    push(@spellWords,map {(my $prefix = $_) =~ s/Class$//; ($_,$prefix)} keys(%{$stateStorables->{'functionClasses'       }}));
	    foreach my $instanceName ( @{$stateStorables->{'functionClassInstances'}} ) {
		my $matchLongest   = 0;
		my $instanceSuffix    ;
		foreach my $className ( keys(%{$stateStorables->{'functionClasses'}}) ) {
		    (my $classPrefix = $className) =~ s/Class$//;
		    if ( $instanceName =~ m/^$classPrefix/ ) {
			if ( length($classPrefix) > $matchLongest ) {
			    $matchLongest = length($classPrefix);
			    ($instanceSuffix = $instanceName) =~ s/^$classPrefix//;
			    $instanceSuffix = lcfirst($instanceSuffix);
			}
		    }
		}
		push(@spellWords,$instanceSuffix)
		    if ( defined($instanceSuffix) );
	    }
	}
	open(my $words,"<","aux/words.dict");
	while ( my $line = <$words> ) {
	    chomp($line);
	    push(@spellWords,$line);
	}
	close($words);
    }
    open(my $dictionary,">",$fileName.".dic");
    print $dictionary scalar(@spellWords)."\n";
    print $dictionary join("\n",uniq(sort(@spellWords)))."\n";
    close($dictionary);
    system("cp aux/words.aff ".$fileName.".aff");

    # Pre-process file.
    open(my $fileIn ,"<",$fileName         );
    open(my $fileOut,">",$fileName.".spell");
    while ( my $line = <$fileIn> ) {
	# Split camelCase words into their components.
	while ( my @parts = $line =~ m/^(.*)\b(?<!\\)([a-zA-Z0-9][a-zA-Z]*(?:[a-z0-9][a-zA-Z]*[A-Z]|[A-Z][a-zA-Z]*[a-z])[a-zA-Z0-9]*)\b(.*)$/g ) {
	    my $linePrevious = $line;
	    if ( grep {$_ eq $parts[1]} @spellWords ) {
		$parts[1] = $parts[1];
	    } elsif ( $parts[1] eq "FoX" ) {
		$parts[1] = "fox";
	    } else {
		$parts[1] =~ s/([a-z0-9])([A-Z0-9])/$1 $2/g;
	    }
	    $line = $parts[0].$parts[1].$parts[2]."\n";
	    last
		if ( $line eq $linePrevious );
	}
	# Remove non-roman text in LaTeX subscripts.
	## Cases where the subscript starts with a "{".
	 if ( $isLaTeX ) {
	    my $lineNew;
	    my $bp = $RE{balanced}{-parens=>'{}'};
	    while ( $line =~ m/^(.*?_)($bp)/ ) {
		$line     = substr($line,length($1)+length($2));
		$lineNew .= $1."{";
		my $subscript = substr($2,1,length($2)-2);
		my $subscriptNew = "";
		while ( $subscript =~ m/^(.*?)(\\mathrm($bp))/ ) {
		    $subscriptNew .= "\\mathrm{".(length($3) > 4 ? $3 : "{}")."}";
		    $subscript     = substr($subscript,length($1)+length($2));
		}
		$subscript = $subscriptNew;
		$lineNew  .= $subscript."}";
	    }
	    $lineNew .= $line;
	    $line     = $lineNew;
	}
	## Cases where the subscript starts with a "\mathrm{".
	if ( $isLaTeX ) {
	    my $lineNew;
	    my $bp = $RE{balanced}{-parens=>'{}'};
	    while ( $line =~ m/^(.*?_\\mathrm)($bp)/ ) {
		$line     = substr($line,length($1)+length($2));
		$lineNew .= $1."{";
		my $subscript = substr($2,1,length($2)-2);
		$subscript = length($subscript) > 2 ? $subscript : "";
		$lineNew      .= $subscript."}";
	    }
	    $lineNew .= $line;
	    $line     = $lineNew;
	}
	# Remove glossary labels.
	 if ( $isLaTeX ) {
	    my $bp = $RE{balanced}{-parens=>'{}'};
	    $line =~ s/\\gls$bp/\\gls\{\}/g;
	    $line =~ s/\\glslink$bp/\\glslink\{\}/g;
	    $line =~ s/\\newacronym$bp$bp/\\newacronym\{\}\{\}/g;
	    $line =~ s/\\newglossaryentry$bp\{name=$bp/\\newglossaryentry\{\}\{/g;
	    $line =~ s/firstplural=/first plural=/g;
	}
	# Translate accents. Note that we translate to unaccented characters, as hunspell seems to not handle accented characters in
	# personal dictionaries correctly.
	if ( $isLaTeX ) {
	    $line =~ s/\\'e/e/g;
	    $line =~ s/\\"o/o/g;
	}
	# Remove URLs.
	if ( $isLaTeX ) {
	    $line =~ s/\\href\{[^\}]+\}//g;
	} else {
	    $line =~ s/http:\/\/.+(\s|\))/$1/g;
	}
	# Write the processed line to output.
	print $fileOut $line;
    }
    close($fileIn );
    close($fileOut);

    # Spell check.
    my %words;
    my $lineNumber = 0;
    open(my $spell,"hunspell -l ".($isLaTeX ? "-t " : "")."-i utf-8 -d en_US,".$fileName." ".$fileName.".spell |");
    while ( my $word = <$spell> ) {
	++$lineNumber; 
	chomp($word);
	$words{lc($word)}++;
    }
    close($spell);   
    unlink($fileName.".spell",$fileName.".dic",$fileName.".aff");
    my $warnings = "";
    foreach my $word ( sort(keys(%words)) ) {
	$warnings .= ":warning: Possible misspelled word '".$word."' ".($words{$word} > 1 ? "(".$words{$word}." instances)" : "")." in file '".$fileNameOriginal."'\n";
    }
    return $warnings;
}

1;

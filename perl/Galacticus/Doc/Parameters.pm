# Contains a Perl module which implements various useful functionality for extracting data on
# parameters from Galacticus source code.

package Parameters;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use XML::Simple;
use Fcntl qw(SEEK_SET);
require Fortran::Utils;

sub ExpandRegEx {
    my $regEx           = shift;
    my $sourceDirectory = shift;
    # Load the file of directive locations.
    my $locations;
    if ( -e "./work/build/Code_Directive_Locations.xml" ) {
	my $xml    = new XML::Simple;
	$locations = $xml->XMLin("./work/build/Code_Directive_Locations.xml");
    }
    # Process the regEx until it no longer contains anything to expand.
    while ( $regEx =~ m/\(\#([a-zA-Z0-9]+)\-\>([a-zA-Z0-9]+)\)/ ) {
	my $directive = $1;
	my $element   = $2;
	# Initialize list of acceptable values.
	my @acceptableValues;
	# Find files to scan.
	my @filesToScan;
	if ( defined($locations) ) {
	    if ( UNIVERSAL::isa($locations->{$directive}->{'file'},'ARRAY') ) {
		@filesToScan = @{$locations->{$directive}->{'file'}};
	    } else {	
		@filesToScan = ( $locations->{$directive}->{'file'} )
		    if (exists($locations->{$directive}->{'file'}));
	    }
	} else {
	    opendir(my $indir,$sourceDirectory) or die "Can't open the source directory: #!";
	    while ( my $fname = readdir $indir) {	
		if ( $fname =~ m/\.[fF](90)??t??$/ && $fname !~ m/^\.\#/ ) {
		    my $fullname = $sourceDirectory."/".$fname;
		    push(@filesToScan,$fullname);
		}
	    }
	}
	# Scan all necessary files.
	foreach my $currentFileName ( @filesToScan ) {
	    # Determine source code type.
	    my $codeType = "fortran";
	    $codeType = "c"
		if ( $currentFileName =~ m/\.c(pp)??$/ );
	    # Add the file to the list of filenames to process.
	    my @fileNames;
	    my @filePositions;
	    unshift(@fileNames,$currentFileName);
	    unshift(@filePositions,-1);
	    # Process files until none remain.
	    while ( scalar(@fileNames) > 0 ) {
		# Open the file.
		open(my $infile,$fileNames[0]) or die "Can't open input file: #!";
		seek($infile,$filePositions[0],SEEK_SET) unless ( $filePositions[0] == -1 );
		until ( eof($infile) ) {
		    my $rawLine;
		    my $processedLine;
		    my $bufferedComments;
		    if ( $codeType eq "fortran" ) {
			# Get next line from the Fortran source.
			&Fortran_Utils::Get_Fortran_Line($infile,$rawLine,$processedLine,$bufferedComments);
		    } elsif ( $codeType eq "c" ) {
			# Get the next line from a C(++) source.
			$processedLine = <$infile>;
			$rawLine       = $processedLine;	
		    }
		    my $lineNumber = $.;	    
		    # Detect include files, and recurse into them.
		    if ( $codeType eq "fortran" && $processedLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
			my $includeFile = $sourceDirectory."/source/".$1;
			$includeFile =~ s/\.inc$/.Inc/;
			if ( -e $includeFile ) {
			    $filePositions[0] = tell($infile);
			    unshift(@fileNames,$includeFile);
			    unshift(@filePositions,-1);
			    last;
			}
		    }	    
		    if ( $rawLine =~ m/^\s*(!|\/\/)\#\s+(<\s*([a-zA-Z]+)+.*>)\s*$/ ) {
			my $xmlCode = $2."\n";
			my $xmlTag  = $3;
			# Read ahead until a matching close tag is found.
			unless ( $xmlCode =~  m/\/>/ ) {
			    my $nextLine = "";
			    until ( $nextLine =~ m/<\/$xmlTag>/ || eof($infile) ) {
				if ( $codeType eq "fortran" ) {
				    &Fortran_Utils::Get_Fortran_Line($infile,$nextLine,$processedLine,$bufferedComments);
				} elsif ( $codeType eq "c" ) {
				    $nextLine = <$infile>;
				}
				$nextLine =~ s/^\s*(!|\/\/)\#\s+//;
				$xmlCode .= $nextLine;
			    }
			}
			# Check if this directive matches that which we are currently processing.
			if ( $xmlTag eq $directive ) {
			    my $xml = new XML::Simple;
			    my $xmlData = eval{$xml->XMLin($xmlCode, ForceArray => ["data","property","binding"])};
			    die("Find_Parameter_Dependencies.pl failed in ".$currentFileName." at line ".$lineNumber." with message:\n".$@)
				if ( $@ );
			    push(@acceptableValues,$xmlData->{$element})
				if ( exists($xmlData->{$element}) );
			}
		    }
		}
		if ( eof($infile) ) {
		    shift(@fileNames);
		    shift(@filePositions);
		}
		close($infile);
	    }
	}
	my $acceptableRegex = "(".join("|",sort(@acceptableValues)).")";
	$regEx =~ s/\(\#$directive\-\>$element\)/$acceptableRegex/g;
    }
    return $regEx;
}

1;

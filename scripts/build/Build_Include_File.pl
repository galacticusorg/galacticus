#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use XML::Simple;
use Data::Dumper;
use Switch;
use Scalar::Util 'reftype';
require Fortran::Utils;
require Galacticus::Build::Hooks;
require Galacticus::Build::ModuleUse;
require Galacticus::Build::Labels;
require Galacticus::Build::FunctionCall;
require Galacticus::Build::Components;
require Galacticus::Build::BindingsC;

# Scans source code for "!#" directives and generates an include file.
# Andrew Benson (18-November-2011)

# Read command line arguments.
die "Usage: Build_Include_File.pl <sourceDirectory> <xmlFile>"
    unless ( scalar(@ARGV) == 2 );
my $sourceDirectory = $ARGV[0];
my $xmlFile         = $ARGV[1];

# Specify verbosity.
$verbosity = 0;

# Load the file of directive locations.
my $locations;
if ( -e "./work/build/Code_Directive_Locations.xml" ) {
    my $xml    = new XML::Simple;
    $locations = $xml->XMLin("./work/build/Code_Directive_Locations.xml");
}

# Create XML object and process the XML file.
my $xml       = new XML::Simple;
my $buildData = $xml->XMLin($xmlFile);

# Initially not inside any module.
$buildData->{'moduleName'} = "";

# Find files to scan.
my @filesToScan;
if ( defined($locations) ) {
    if ( reftype($locations->{$buildData->{'directive'}}->{'file'}) eq 'ARRAY' ) {
	@filesToScan = @{$locations->{$buildData->{'directive'}}->{'file'}};
    } else {	
	@filesToScan = ( $locations->{$buildData->{'directive'}}->{'file'} )
	    if (exists($locations->{$buildData->{'directive'}}->{'file'}));
    }
} else {
    opendir(indir,$sourceDirectory) or die "Can't open the source directory: #!";
    while ( my $fname = readdir indir) {	
	if ( $fname =~ m/\.[fF](90)??t??$/ && $fname !~ m/^\.\#/ ) {
	    my $fullname = $sourceDirectory."/".$fname;
	    push(@filesToScan,$fullname);
	}
    }
}

# Scan all necessary files.
foreach my $currentFileName ( @filesToScan ) {
    # Store the current file name.
    $buildData->{'currentFileName'} = $currentFileName;

    # Determine source code type.
    $buildData->{'codeType'} = "fortran";
    $buildData->{'codeType'} = "c"
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
	    switch ( $buildData->{'codeType'} ) {
		case ( "fortran" ) {
		    # Get next line from the Fortran source.
		    &Fortran_Utils::Get_Fortran_Line($infile,$rawLine,$processedLine,$bufferedComments);
		}
		case ( "c" ) {
		    # Get the next line from a C(++) source.
		    $processedLine = <$infile>;
		    $rawLine       = $processedLine;	
		}
	    }
	    my $lineNumber = $.;
	    
	    # Detect include files, and recurse into them.
	    if ( $buildData->{'codeType'} eq "fortran" && $processedLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
		my $includeFile = $sourceDirectory."/source/".$1;
		$includeFile =~ s/\.inc$/.Inc/;
		if ( -e $includeFile ) {
		    $filePositions[0] = tell($infile);
		    unshift(@fileNames,$includeFile);
		    unshift(@filePositions,-1);
		    last;
		}
	    }
	    
	    # Check if we've entered or left a module.
	    if ( $processedLine =~ /^\s*module\s*([a-z0-9_]+)\s*$/i       ) {$buildData->{'moduleName'} = $1};
	    if ( $processedLine =~ /^\s*end\s+module\s*([a-z0-9_]+)\s*$/i ) {$buildData->{'moduleName'} = ""};
	    
	    if ( $rawLine =~ m/^\s*(!|\/\/)\#\s+(<\s*([a-zA-Z]+)+.*>)\s*$/ ) {
		my $xmlCode = $2."\n";
		my $xmlTag  = $3;
		# Read ahead until a matching close tag is found.
		unless ( $xmlCode =~  m/\/>/ ) {
		    my $nextLine = "";
		    until ( $nextLine =~ m/<\/$xmlTag>/ || eof($infile) ) {
			switch ( $buildData->{'codeType'}) {
			    case ( "fortran" ) {
				&Fortran_Utils::Get_Fortran_Line($infile,$nextLine,$processedLine,$bufferedComments);
			    }
			    case ( "c" ) {
				$nextLine = <$infile>;
			    }
			}
			$nextLine =~ s/^\s*(!|\/\/)\#\s+//;
			$xmlCode .= $nextLine;
		    }
		}
		# Check if this directive matches that which we are currently processing.
		if ( $xmlTag eq $buildData->{'directive'} ) {
		    $buildData->{'currentDocument'} = eval{$xml->XMLin($xmlCode, ForceArray => ["data","method","binding"])};
		    die("Build_Include_File.pl failed in ".$currentFileName." at line ".$lineNumber." with message:\n".$@)
			if ( $@              );
		    print Dumper($buildData->{'currentDocument'})
			if ( $verbosity == 1 );

		    # Look for a match for this action type and call the relevant function to parse it.
		    my $foundMatch = 0;
		    foreach my $hook ( keys(%Hooks::moduleHooks) ) {
			if ( $buildData->{'type'} eq $hook ) {
			    $foundMatch = 1;
			    my $parseFunction = $Hooks::moduleHooks{$hook}->{'parse'};
			    &{$parseFunction}($buildData);
			}
		    }
		    die("Build_Include_File.pl: failed to find a function to parse ".$buildData->{'type'}." action")
			if ( $foundMatch == 0 );
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

# Look for a match for this action type and call the relevant function to generate content.
    my $foundMatch = 0;
    foreach my $hook ( keys(%Hooks::moduleHooks) ) {
	if ( $buildData->{'type'} eq $hook ) {
	    $foundMatch = 1;
	    my $generateFunction = $Hooks::moduleHooks{$hook}->{'generate'};
	    &{$generateFunction}($buildData);
	}
    }
    die("Build_Include_File.pl: failed to find a function to generate ".$buildData->{'type'}." action")
	if ( $foundMatch == 0 );
    
# Output the generated content.
    open(includeFile,">".$buildData->{'fileName'});
    print includeFile $buildData->{'content'};
    close(includeFile);

    exit;

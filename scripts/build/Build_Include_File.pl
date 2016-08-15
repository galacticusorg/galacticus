#!/usr/bin/env perl
use strict;
use warnings;
our $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use XML::Simple;
use Data::Dumper;
use Scalar::Util 'reftype';
use Fcntl qw(SEEK_SET);
use UNIVERSAL;
require Fortran::Utils;
require Galacticus::Build::Hooks;
require Galacticus::Build::ModuleUse;
require Galacticus::Build::MethodNames;
require Galacticus::Build::Labels;
require Galacticus::Build::Function;
require Galacticus::Build::FunctionCall;
require Galacticus::Build::BindingsC;
require Galacticus::Build::FunctionGlobal;
require Galacticus::Build::SourceTree;

# Scans source code for "!#" directives and generates an include file.
# Andrew Benson (18-November-2011)

# Read command line arguments.
die "Usage: Build_Include_File.pl <sourceDirectory> <xmlFile>"
    unless ( scalar(@ARGV) == 2 );
my $sourceDirectory = $ARGV[0];
my $xmlFile         = $ARGV[1];

# Specify verbosity.
my $verbosity = 0;

# Set load status of large modules.
my $componentsLoaded = 0;

# Load the file of directive locations.
my $locations;
if ( -e $ENV{'BUILDPATH'}."/Code_Directive_Locations.xml" ) {
    my $xml    = new XML::Simple;
    $locations = $xml->XMLin($ENV{'BUILDPATH'}."/Code_Directive_Locations.xml");
}

# Create XML object and process the XML file.
my $xml       = new XML::Simple;
my $buildData = $xml->XMLin($xmlFile, KeyAttr => []);

# Initially not inside any module.
$buildData->{'moduleName'} = "";

# Find files to scan.
my @filesToScan;
if ( defined($locations) ) {
    if ( UNIVERSAL::isa($locations->{$buildData->{'directive'}}->{'file'},'ARRAY') ) {
	@filesToScan = @{$locations->{$buildData->{'directive'}}->{'file'}};
    } else {	
	@filesToScan = ( $locations->{$buildData->{'directive'}}->{'file'} )
	    if (exists($locations->{$buildData->{'directive'}}->{'file'}));
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
	    if ( $buildData->{'codeType'} eq "fortran" ) {
		# Get next line from the Fortran source.
		&Fortran_Utils::Get_Fortran_Line($infile,$rawLine,$processedLine,$bufferedComments);
	    } elsif ( $buildData->{'codeType'} eq "c" ) {
		# Get the next line from a C(++) source.
		$processedLine = <$infile>;
		$rawLine       = $processedLine;	
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
	    
	    my $includeFile;
	    if ( $rawLine =~ m/^\s*(!|\/\/)\#\s+(<\s*([a-zA-Z]+)+.*>)\s*$/ ) {
		my $xmlCode = $2."\n";
		my $xmlTag  = $3;
		# Read ahead until a matching close tag is found.
		unless ( $xmlCode =~  m/\/>/ ) {
		    my $nextLine = "";
		    until ( $nextLine =~ m/<\/$xmlTag>/ || eof($infile) ) {
			if ( $buildData->{'codeType'} eq "fortran" ) {
				&Fortran_Utils::Get_Fortran_Line($infile,$nextLine,$processedLine,$bufferedComments);
			} elsif ( $buildData->{'codeType'} eq "c" ) {
			    $nextLine = <$infile>;
			}
			# Check for included files.
			if ( $buildData->{'codeType'} eq "fortran" && $nextLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ ) {
			    $includeFile = $sourceDirectory."/".$ENV{'BUILDPATH'}."/".$1;
			    $includeFile =~ s/\.inc$/.Inc/;
			}
			# Add the line to our XML.
			$nextLine =~ s/^\s*(!|\/\/)\#\s+//;
			$xmlCode .= $nextLine;
		    }
		}
		# Check if this directive matches that which we are currently processing.
		if ( $xmlTag eq $buildData->{'directive'} ) {
		    $buildData->{'currentDocument'} = eval{$xml->XMLin($xmlCode, ForceArray => ["data","property","binding"])};
		    die("Build_Include_File.pl failed in ".$currentFileName." at line ".$lineNumber." with message:\n".$@)
			if ( $@              );
		    print Dumper($buildData->{'currentDocument'})
			if ( $verbosity == 1 );
		    # Load large modules needed for this action type.
		    if ( $componentsLoaded == 0 && $buildData->{'type'} eq "component" ) {
			require Galacticus::Build::Components;
			$componentsLoaded = 1;
		    }
		    # Look for a match for this action type and call the relevant function to parse it.
		    my $foundMatch = 0;
		    foreach my $hook ( keys(%Hooks::moduleHooks) ) {
			if ( $buildData->{'type'} eq $hook ) {
			    $foundMatch = 1;
			    if ( exists($Hooks::moduleHooks{$hook}->{'validate'}) ) {
				my $validateFunction = $Hooks::moduleHooks{$hook}->{'validate'};
				&{$validateFunction}($xmlCode,$currentFileName);
			    }
			    my $parseFunction = $Hooks::moduleHooks{$hook}->{'parse'};
			    &{$parseFunction}($buildData);
			}
		    }
		    die("Build_Include_File.pl: failed to find a function to parse ".$buildData->{'type'}." action")
			if ( $foundMatch == 0 );
		}
	    }
	    # Process any include file that was found.
	    if ( defined($includeFile) && -e $includeFile ) {
		$filePositions[0] = tell($infile);
		unshift(@fileNames,$includeFile);
		unshift(@filePositions,-1);
		last;
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

# Parse Fortran files, simply output other files.
if ( $buildData->{'fileName'} =~ m/\.Inc$/ ) {
    # Parse the file to build a tree.
    my $tree = &SourceTree::ParseCode($buildData->{'content'},$buildData->{'fileName'});
    # Process the tree.
    &SourceTree::ProcessTree($tree);
    # Serialize back to source code.
    open(my $outputFile,">",$buildData->{'fileName'});
    print $outputFile &SourceTree::Serialize($tree);
    close($outputFile);
} else {
    open(my $outputFile,">",$buildData->{'fileName'});
    print $outputFile $buildData->{'content'};
    close($outputFile);
}

exit;

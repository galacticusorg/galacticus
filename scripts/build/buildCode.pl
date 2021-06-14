#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Data::Dumper;
use Scalar::Util 'reftype';
use Fcntl qw(SEEK_SET);
use UNIVERSAL;
use Fortran::Utils;
use List::ExtraUtils;
use File::Changes;
use Storable;
use Galacticus::Build::Hooks;
use Galacticus::Build::ModuleUse;
use Galacticus::Build::FunctionCall;
use Galacticus::Build::FunctionGlobal;
use Galacticus::Build::SourceTree;

# Scans source code for "!#" directives and generates code from these.
# Andrew Benson (18-November-2011)

# Read command line arguments.
die "Usage: buildCode.pl <sourceDirectory> <xmlFile>"
    unless ( scalar(@ARGV) == 2 );
my $sourceDirectoryName = $ARGV[0];
my $xmlFile             = $ARGV[1];
# Specify verbosity.
my $verbosity           = 0;
# Set load status of large modules.
my $componentsLoaded    = 0;
# Get an XML parser.
my $xml                 = new XML::Simple;
# Load the file of directive locations if available.
my $locations           = -e $ENV{'BUILDPATH'}."/directiveLocations.xml" ? $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml") : undef();
# Process the XML file.
my $build               = $xml->XMLin($xmlFile, KeyAttr => []);
# Initialize structure to hold record of directives from each source file.
my $directivesPerFile;
my $havePerFile = -e $build->{'fileName'}.".blob";
my $updateTime;
if ( $havePerFile ) {
    $directivesPerFile = retrieve($build->{'fileName'}.".blob");
    $updateTime        = -M       $build->{'fileName'}.".blob" ;
}
# Initialize to be not inside any module.
$build->{'moduleName'}   = "";
# Find files to scan.
my @fileNamesToScan;
if ( defined($locations) ) {
    # A list of files containing the directive of interest is available - use it.
    @fileNamesToScan = &List::ExtraUtils::as_array($locations->{$build->{'directive'}}->{'file'});
} else {
    # Alternative, just scan all source files.
    opendir(my $sourceDirectory,$sourceDirectoryName) 
	or die "buildCode.pl: can not open the source directory: #!";
    @fileNamesToScan = map {$_ =~ m/\.[fF](90)??t??$/ && $_ !~ m/^\.\#/ ? $sourceDirectoryName."/".$_ : ()} readdir($sourceDirectory);
    closedir($sourceDirectory);
}
# Scan all files of interest.
foreach my $fileName ( @fileNamesToScan ) {
    # Store the current file name.
    $build->{'currentFileName'} = $fileName;
    # Determine source code type.
    $build->{'codeType'       } = $fileName =~ m/\.c(pp)??$/ ? "c" : "fortran";
    # Add the file name to the stack of file names to process.
    my @fileStack = ( { name => $fileName, position => -1 } );
    (my $fileIdentifier = $fileName) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;
    # Check if file is updated. If it is not, skip processing it. If it is, remove previous record of directives and rescan.
    if ( $havePerFile && exists($directivesPerFile->{$fileIdentifier}) ) {
	next
	    unless ( grep {-M $_ < $updateTime} &List::ExtraUtils::as_array($directivesPerFile->{$fileIdentifier}->{'files'}) );
    }
    delete($directivesPerFile->{$fileIdentifier})
    	    if ( $havePerFile && exists($directivesPerFile->{$fileIdentifier}) );
    push(@{$directivesPerFile->{$fileIdentifier}->{'files'}},$fileName);
    # Process files until none remain.
    while ( scalar(@fileStack) > 0 ) {	
	# Pop a file from the stack and move back to the previous position.
	my $fileProcess = pop(@fileStack);
	open(my $file,$fileProcess->{'name'}) 
	    or die "buildCode.pl: can not open input file: ".$fileProcess->{'name'};
	seek($file,$fileProcess->{'position'},SEEK_SET) 
	    unless ( $fileProcess->{'position'} == -1 );
	# Read the file until end of file is reached.
	until ( eof($file) ) {
	    # Read a line from the file.
	    my $rawLine;
	    my $processedLine;
	    my $bufferedComments;
	    my $includeFile;
	    if ( $build->{'codeType'} eq "fortran" ) {
		# Get next line from the Fortran source.
		&Fortran::Utils::Get_Fortran_Line($file,$rawLine,$processedLine,$bufferedComments);
	    } elsif ( $build->{'codeType'} eq "c" ) {
		# Get the next line from a C(++) source.
		$processedLine = <$file>;
		$rawLine       = $processedLine;	
	    }
	    # Record the current line number for any subsequent error reporting.
	    my $lineNumber = $.;
	    # Fortran-specific processing.
	    if ( $build->{'codeType'} eq "fortran" ) {
		# Detect include files.
		($includeFile = $sourceDirectoryName."/source/".$1) =~ s/\.inc$/.Inc/
		    if ( $processedLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ );
		# Check if we've entered or left a module.
		if ( my @matches = $processedLine =~ $Fortran::Utils::unitOpeners{'module'}->{'regEx'} ) 
	            {$build->{'moduleName'} = $matches[$Fortran::Utils::unitOpeners{'module'}->{'unitName'}]};
		if (               $processedLine =~ $Fortran::Utils::unitClosers{'module'}->{'regEx'} )
	            {$build->{'moduleName'} = ""                                                            };
	    }
	    if ( $rawLine =~ m/^\s*(!|\/\/)\#\s+(<\s*([a-zA-Z]+)+.*>)\s*$/ ) {
		my $xmlCode = $2."\n";
		my $xmlTag  = $3;
		# Read ahead until a matching close tag is found.
		unless ( $xmlCode =~  m/\/>/ ) {
		    my $nextLine = "";
		    until ( $nextLine =~ m/<\/$xmlTag>/ || eof($file) ) {
			if ( $build->{'codeType'} eq "fortran" ) {
			    &Fortran::Utils::Get_Fortran_Line($file,$nextLine,$processedLine,$bufferedComments);
			} elsif ( $build->{'codeType'} eq "c" ) {
			    $nextLine = <$file>;
			}
			# Check for included files.
			($includeFile = $sourceDirectoryName."/".$ENV{'BUILDPATH'}."/".$1) =~ s/\.inc$/.Inc/
			    if ( $build->{'codeType'} eq "fortran" && $nextLine =~ m/^\s*include\s*['"]([^'"]+)['"]\s*$/ );
			# Add the line to our XML.
			$nextLine =~ s/^\s*(!|\/\/)\#\s+//;
			$xmlCode .= $nextLine;
		    }
		}
		# Check if this directive matches that which we are currently processing.
		if ( $xmlTag eq $build->{'directive'} ) {
		    my $directive;
		    $directive->{'fileName'  } = $fileName;
		    $directive->{'xmlCode'   } = $xmlCode;
		    $directive->{$_          } = $build->{$_}
		        foreach ( "moduleName", "currentFileName", "codeType" );
		    $directive->{'directive' } = eval{$xml->XMLin($xmlCode, ForceArray => ["data","property","binding"])};
		    die("buildCode.pl: failed in ".$fileProcess->{'name'}." at line ".$lineNumber." with message:\n".$@)
			if ( $@ );
		    push(@{$directivesPerFile->{$fileIdentifier}->{'directives'}},$directive);
		    print Dumper($build->{'currentDocument'})
			if ( $verbosity == 1 );
		}
	    }
	    # If an include file was found, push the current file back onto the stack, followed by the include file, and finish
	    # processing the current file. We process files in this way (rather than simply pushing all include files to the stack
	    # and processing them after the main file) because we want to know the module context within which the file is
	    # included.
	    if ( defined($includeFile) && -e $includeFile ) {
		$fileProcess->{'position'} = tell($file);
		push(@fileStack,$fileProcess,{name => $includeFile, position => -1});
		push(@{$directivesPerFile->{$fileIdentifier}->{'files'}},$includeFile);
		last;
	    }
	}
	close($file);
    }
}
# Output the per file directives.
store($directivesPerFile,$build->{'fileName'}.".blob");
# Load large modules needed for this action type.
if ( ! $componentsLoaded && $build->{'type'} eq "component" ) {
    require Galacticus::Build::Components;
    $componentsLoaded = 1;
}
# Validate and parse all directives.
foreach my $fileIdentifier ( keys(%{$directivesPerFile}) ) {
    foreach my $directive ( &List::ExtraUtils::as_array($directivesPerFile->{$fileIdentifier}->{'directives'}) ) {
	$build->{'currentDocument'} = $directive->{'directive' };
	$build->{$_               } = $directive->{$_          }
	    foreach ( "moduleName", "currentFileName", "codeType" );
	# Validate and parse this directive using the appropriate handler.
	die("buildCode.pl: failed to find a function to parse '".$build->{'type'}."' action")
	    unless ( exists($Galacticus::Build::Hooks::moduleHooks{$build->{'type'}}) );
	&{$Galacticus::Build::Hooks::moduleHooks{$build->{'type'}}->{'validate'}}($directive->{'xmlCode'},$directive->{'fileName'})
	    if ( exists($Galacticus::Build::Hooks::moduleHooks{$build->{'type'}}->{'validate'}) );
	&{$Galacticus::Build::Hooks::moduleHooks{$build->{'type'}}->{'parse'   }}($build                                          );	
    }
}
# Call the relevant function to generate content for this action.
die("buildCode.pl: failed to find a function to generate '".$build->{'type'}."' action")
    unless ( exists($Galacticus::Build::Hooks::moduleHooks{$build->{'type'}}) );
&{$Galacticus::Build::Hooks::moduleHooks{$build->{'type'}}->{'generate'}}($build);
# Generate output. For Fortran source we run the code through the processor first. Otherwise it is simply output.
open(my $outputFile,">",$build->{'fileName'}.".tmp");
# Parse Fortran files, simply output other files.
print $outputFile 
    $build->{'fileName'} =~ m/\.Inc$/
    ?
    &Galacticus::Build::SourceTree::Serialize(
	&Galacticus::Build::SourceTree::ProcessTree(
	     &Galacticus::Build::SourceTree::ParseCode($build->{'content'},$build->{'fileName'})
	)
    )
    :
    $build->{'content'};
close($outputFile);
&File::Changes::Update($build->{'fileName'},$build->{'fileName'}.".tmp", proveUpdate => "yes");
exit;

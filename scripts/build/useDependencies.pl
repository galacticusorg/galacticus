#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use List::Uniq ':all';
use Data::Dumper;
use XML::Simple;
use Galacticus::Build::Directives;
use Fortran::Utils;
use List::ExtraUtils;

# Locate source files which have dependencies on modules.
# Andrew Benson (06-September-2016)

# Define the source directory
die "Usage: findUseDependencies.pl sourcedir"
    unless ( scalar(@ARGV) == 1 );
my $rootSourceDirectoryName = $ARGV[0];
# Specify work directory.
my $workDirectoryName       = $ENV{'BUILDPATH'}."/";
# Get an XML parser.  
my $xml                     = new XML::Simple();
# Load the file of directive locations.
my $locations               = -e $workDirectoryName."directiveLocations.xml" ? $xml->XMLin($workDirectoryName."directiveLocations.xml") : undef();
# List of external modules (which will be ignored for dependency analysis of the source code).
my @externalModules = ( "omp_lib", "hdf5", "h5tb", "h5lt", "h5global", "h5fortran_types", "fox_common", "fox_dom", "fox_wxml", "fox_utils",
			"fgsl", "mpi", "yeplibrary", "yepcore", "yepmath");
# Modules that require a library to be linked. These are key-value pairs with the key being the module name, and the value the
# name of the required library.
my %moduleLibararies = (
    nearest_neighbors => "ANN",
    fftw3             => "fftw3",
    fgsl              => "fgsl_gfortran",
    fox_common        => "FoX_common",
    fox_dom           => "FoX_dom",
    fox_wxml          => "FoX_wxml",
    fox_utils         => "FoX_utils",
    hdf5              => "hdf5_fortran",
    h5tb              => "hdf5hl_fortran",
    vectors           => "blas",
    yeplibrary        => "yeppp",
    yepcore           => "yeppp",
    yepmath           => "yeppp"
    );
# C includes that require a library to be linked. These are key-value pairs with the key being the include name, and the value the
# name of the required library.
my %includeLibararies = (
    crypt             => "crypt"
    );
# Parse the Makefile to find preprocessor macros that are set.
my @preprocessorDirectives;
my @conditionsStack;
open(my $makefile,"Makefile");
while ( my $line = <$makefile> ) {
    # Build a stack of true/false values indicating whether the current line is active or not based on pre-processor arguments.
    push(@conditionsStack,exists($ENV{$1}))
	if ( $line =~ m/^ifdef ([A-Z]+)/ );
    push(@conditionsStack,1)
	if ( $line =~ m/^ifeq /          );
    pop (@conditionsStack  )
	if ( $line =~ m/^endif/          );
    # If this line sets compiler flags, and no false entries exist in the preprocessor conditions stack, extract all preprocessor
    # directives from the line.
    push(@preprocessorDirectives,map {$_ =~ m/\-D([0-9A-Z]+)/ ? $1 : ()} split(" ",$line))
	if ( ! grep {! $_} @conditionsStack && $line =~ m/^\s*FCFLAGS\s*\+??=/ );
}
close($makefile);
# Extract any preprocessor directives specified via the GALACTICUS_FCFLAGS environment variable.
push(@preprocessorDirectives,map {$_ =~ m/\-D([0-9A-Z]+)/ ? $1 : ()} split(" ",$ENV{"GALACTICUS_FCFLAGS"}));
# Open the output dependencies file
open(my $dependenciesFile,">".$workDirectoryName."Makefile_Use_Dependencies");
# Build a list of source directories (including all sub-directories of the main directory).
my @sourceDirectoryNames  = ( $rootSourceDirectoryName."/source" );
if ( -e $sourceDirectoryNames[0] ) {
    opendir(my $sourceDirectory,$sourceDirectoryNames[0]);
    push(@sourceDirectoryNames,map {chomp($_);my $path = $sourceDirectoryNames[0]."/".$_; (-d $path && $_ !~ m/^\.+$/) ? $path : ()} readdir($sourceDirectory) );	 
    closedir($sourceDirectory);
}
# Iterate over source directories.
my @sourceFilesToProcess;
foreach my $sourceDirectoryName ( @sourceDirectoryNames ) {
    # Extract the subdirectory name.
    (my $subDirectoryName = $sourceDirectoryName) =~ s/^$rootSourceDirectoryName\/source\/?//;
    # Find all source files to process.
    opendir(my $sourceDirectory,$sourceDirectoryName) 
	or die "findUseDependencies.pl: can not open the source directory: #!";
    push
	(
	 @sourceFilesToProcess,
	 map 
	 {
	     ($_ =~ m/\.(f|f90|c|cpp|inc)$/i && $_ !~ m/^\.\#/) 
	     ?
	     {
		 fileName         =>                          $_, 
		 fullPathFileName => $sourceDirectoryName."/".$_,
		 subDirectoryName => $subDirectoryName
	     }
	     :
	     ()
	 }
	 readdir($sourceDirectory));
    closedir($sourceDirectory);
}
# Iterate over files to process.
foreach my $sourceFile ( @sourceFilesToProcess ) {
    # Initialize records of modules used, provided, explicit dependencies, and any library dependencies.
    my @modulesUsed;
    my %modulesProvided;
    my @dependenciesExplicit;
    my %libraryDependencies;
    # Push the main file onto the scan stack.
    my @fileNamesToProcess = ( $sourceFile->{'fullPathFileName'} );
    # Extract lists of directives from this file which require special handling.
    my $directives;
    @{$directives->{$_}} = &Galacticus::Build::Directives::Extract_Directives($sourceFile->{'fullPathFileName'},$_)
	foreach ( "functionClass", "inputParameter", "enumeration" );
    # Special handling for functionClass directives - add implementation files to the list of files to scan.
    &List::ExtraUtils::smart_push(\@fileNamesToProcess,$locations->{${$directives->{'functionClass'}}[0]->{'name'}}->{'file'})
	if ( scalar(@{$directives->{'functionClass'}}) > 0 );
    # Add dependence on input parameters module if necessary.
    push(@modulesUsed,$workDirectoryName."input_parameters.mod")
	if ( scalar(@{$directives->{'functionClass'}}) > 0 ||  scalar(@{$directives->{'inputParameter'}}) > 0 );
    # Add dependence on error reporting module if necessary.
    push(@modulesUsed,$workDirectoryName."galacticus_error.mod")
	if ( grep {exists($_->{'encodeFunction'}) && $_->{'encodeFunction'} eq "yes"} @{$directives->{'enumeration'}} );
    # Find modules used in functionClass directives.
    foreach my $functionClass ( @{$directives->{'functionClass'}} ) {
	next 
	    unless ( exists($functionClass->{'method'}) );
	foreach my $method ( exists($functionClass->{'method'}->{'name'}) ? $functionClass->{'method'} : map {$functionClass->{'method'}->{$_}} keys(%{$functionClass->{'method'}}) ) {
	    next
		unless ( exists($method->{'modules'}) );
	    push
		(
		 @modulesUsed,
		 map {$_ eq "hdf5" ? () : $workDirectoryName.$_.".mod"} split(" ",lc($method->{'modules'}))
		);
	}
    }
    # Scan files on the stack until stack is empty.
    while ( scalar(@fileNamesToProcess) > 0 ) {
	my $fullPathFileName = pop(@fileNamesToProcess);
	# Make the file if necessary.
	unless ( -e $fullPathFileName ) {
	    my $leaf = $fullPathFileName =~ m/\/([\w\.]+?)$/ ? $1 : $fullPathFileName;
	    system("make ".$leaf);
	}
	# Add a dependency on libstdc++ for any C++ source file.
	$libraryDependencies{"stdc++"} = 1
	    if ( $fullPathFileName =~ m/\.cpp$/i );
	# Initialize preprocessor conditional compilation state and state stack.
	my @preprocessorConditionalsStack;
	my $conditionallyCompile = 1;
	open(my $file,$fullPathFileName) or die "Can't open input file: $fullPathFileName";
	while (my $line = <$file>) {
	    if ( $line =~ m/^\s*\!;\s*([a-zA-Z0-9_]+)\s*$/ ) {
		$libraryDependencies{$1} = 1;	
	    }
	    if ( $line =~ m/^\s*\#include\s+<([a-zA-Z0-9_]+)\.h>/ ) {
		my $includeFile = $1;
		$libraryDependencies{$includeLibararies{lc($includeFile)}} = 1
		    if ( exists($includeLibararies{lc($includeFile)}) );
	    }
	    # Detect preprocessor lines.
	    if ( $line =~ m/^\#/ ) {
		# Build a stack of preprocessor conditional directives. The stack stores the directive name (where possible) and
		# the sign of the logic (1 for "compile-if-true", 0 for "compile-if-false"). Preprocessor "else" directives simply
		# invert the state of the last entry on the stack.
		push(@preprocessorConditionalsStack,{name => $1           , state => 1})
		    if ( $line =~ m/^\#ifdef\s+([0-9A-Za-z_]+)\s*$/ );
		push(@preprocessorConditionalsStack,{name => "conditional", state => 1})
		    if ( $line =~ m/^\#if\s/                        );
		push(@preprocessorConditionalsStack,{name => $1           , state => 0})
		    if ( $line =~ m/^\#ifndef\s+([0-9A-Z_]+)\s*$/   );
		pop (@preprocessorConditionalsStack                                    )
		    if ( $line =~ m/^\#endif\s*$/                   );
		$preprocessorConditionalsStack[-1]->{'state'} = 1-$preprocessorConditionalsStack[-1]->{'state'}
		if ( $line =~ m/^\#else\s*$/                    );
		# Determine whether or not the current code will be conditionally compiled.
		$conditionallyCompile = 1;
		foreach my $preprocessorConditional ( @preprocessorConditionalsStack ) {
		    my $conditionalActive = 
			(grep {$_ eq $preprocessorConditional->{'name'}} @preprocessorDirectives)
			?
			  $preprocessorConditional->{'state'} 
		        : 
			1-$preprocessorConditional->{'state'};
		    $conditionallyCompile = 0	
			if ( $conditionalActive == 0 );
		}
	    }
	    # Process line only if conditional compilation state is active.
	    if ( $conditionallyCompile == 1 ) {
		# Locate any lines which use the "use" statement and extract the name of the file they use. Any externally
		# provided modules are excluded.
		if ( $line =~ m/^\s*use\s+([a-zA-Z0-9_]+)/i ) {
		    my $usedModule = $1;
		    # Add any library dependency for this module.
		    $libraryDependencies{$moduleLibararies{lc($usedModule)}} = 1
			if ( exists($moduleLibararies{lc($usedModule)}) );
		    push(@modulesUsed,$workDirectoryName.lc($usedModule).".mod")
			unless ( grep {$_ eq lc($usedModule)} @externalModules );
		}
		# Locate explicit dependencies.
		push(@dependenciesExplicit,split(" ",$2))
		    if ( $line =~ m/^\s*(\!|\/\/):\s*(.*)$/ );
		# Locate any modules provided by this file (we do not need to include an explicit dependence on any modules which
		# are self-provided).
		$modulesProvided{$1.".mod"} = 1
		    if ( $line =~ m/^\s*module\s+([a-zA-Z0-9_]+) / );
		# Locate included files and push them to the stack of files to process.
		if ( $line =~ m/^\s*include\s+(\'|\")([\w\.\-]+)(\'|\")/i ) {
		    my $preprocessedIncludedFile = $2;
		    (my $rawIncludedFile = $preprocessedIncludedFile) =~ s/\.inc$/.Inc/;
		    if ( -e $rootSourceDirectoryName."/source/".$rawIncludedFile ) {
			# A raw (unpreprocessed) matching file exists in the source directory. Use this as the dependency.
			push(@fileNamesToProcess,$rootSourceDirectoryName."/source/".$rawIncludedFile);
		    } elsif ( -e $workDirectoryName.$preprocessedIncludedFile ) {
			# A preprocessed matching file exists in the work directory. Use it as the dependency.
			push(@fileNamesToProcess,$workDirectoryName.$preprocessedIncludedFile);
		    } elsif ( $preprocessedIncludedFile =~ m/(.*)\.type\.inc/ ) {
			# For old-style method include files, add a dependency on all files which contain the associated
			# directive.
			&List::ExtraUtils::smart_push(\@fileNamesToProcess,$locations->{$1}->{'file'});
		    }
		}
	    }
	}
	close($file);
    }
    # Construct name of assoicated object and dependency file.
    (my $objectFileName      = $sourceFile->{'fileName'}) =~ s/\.(f|f90|c|cpp)$/.o/i;
    (my $dependencyFileName  = $sourceFile->{'fileName'}) =~ s/\.(f|f90|c|cpp|inc)$/.d/i;
    # Construct name of work subdirectory.
    my $workSubDirectoryName = $workDirectoryName.$sourceFile->{'subDirectoryName'}.($sourceFile->{'subDirectoryName'} eq "" ? "" : "/");
    # Output library file rule.
    unless ( $objectFileName =~ m/\.Inc$/ ) {
	(my $libraryFileName = $objectFileName) =~ s/.o$/.fl/;
	print $dependenciesFile $workSubDirectoryName.$libraryFileName,":\n";
	print $dependenciesFile "\t\@echo -n      > ".$workSubDirectoryName.$libraryFileName."\n";
	print $dependenciesFile "\t\@echo ".$_." >> ".$workSubDirectoryName.$libraryFileName."\n"
	    foreach ( keys(%libraryDependencies) );
    }    
    # For files which used any modules, generate dependency rules.
    if ( scalar(@modulesUsed) > 0 || scalar(@dependenciesExplicit) > 0 ) {
	# Sort the list of modules used, and remove any duplicate entries.
	@modulesUsed = grep {! exists($modulesProvided{$_})} uniq(sort(@modulesUsed));
	# Output the dependencies.
	print $dependenciesFile $workSubDirectoryName.$objectFileName,": ".join(" ",@modulesUsed,@dependenciesExplicit)." Makefile\n";
	# Generate rules for dependency files - we first append a ".d" to any module file names used.
	my @dependenciesUsed = map {$_ =~ m/\.mod$/ ? $_.".d" : $_} @modulesUsed;	
	print $dependenciesFile $workSubDirectoryName.$dependencyFileName,": ".join(" ",@dependenciesUsed,map {(my $modifiedName = $_) =~ s/\.o$/.d/; $modifiedName} @dependenciesExplicit)."\n";
	print $dependenciesFile "\t\@echo ".$workSubDirectoryName.$objectFileName." > ".$workSubDirectoryName.$dependencyFileName."\n";
	foreach my $dependencyExplicit ( @dependenciesExplicit ) {
	    (my $dependencyExplicitFileName = $dependencyExplicit) =~ s/\.o$/.d/;
	    print $dependenciesFile "\t\@cat ".($dependencyExplicit =~ m/\// ? "" : $workDirectoryName).$dependencyExplicitFileName." >> ".$workSubDirectoryName.$dependencyFileName."\n";
	}
	print $dependenciesFile "\t\@cat ".$_." >> ".$workSubDirectoryName.$dependencyFileName."\n"
	    foreach ( @dependenciesUsed );
	print $dependenciesFile "\t\@sort -u ".$workSubDirectoryName.$dependencyFileName." -o ".$workSubDirectoryName.$dependencyFileName."\n\n";
	# Create rules for making dependency trees with GraphViz.
	my @graphVizesUsed   = map {$_ =~ m/\.d$/ ? $_.".gv" : $_} @modulesUsed;		 
	my $graphVizFileName = $sourceFile->{'fileName'}.".gv";
	print $dependenciesFile $workSubDirectoryName.$graphVizFileName,": ".$workSubDirectoryName.$dependencyFileName." ".join(" ",@graphVizesUsed)."\n";
	print $dependenciesFile "\t\@echo \\\"".$sourceFile->{'subDirectoryName'}.$sourceFile->{'fileName'}."\\\" > ".$workSubDirectoryName.$graphVizFileName."\n";
	foreach my $dependencyExplicit ( @dependenciesExplicit ) {
	    (my $dependencyExplicitFileName = $dependencyExplicit) =~ s/\.o$/.d/;
	    print $dependenciesFile "\t\@awk '{print \"\\\"".$sourceFile->{'subDirectoryName'}.$sourceFile->{'fileName'}."\\\" -> \\\"\"\$\$1\"\\\"\"}' ".($dependencyExplicit =~ m/\// ? "" : $workDirectoryName).$dependencyExplicitFileName." >> ".$workSubDirectoryName.$graphVizFileName."\n";
	}
	foreach my $graphVizUsed ( @graphVizesUsed ) {
	    print $dependenciesFile "\t\@awk '{print \"\\\"".$sourceFile->{'subDirectoryName'}.$sourceFile->{'fileName'}."\\\" -> \\\"\"\$\$1\"\\\"\"}' ".$graphVizUsed." >> ".$workSubDirectoryName.$graphVizFileName."\n";
	    print $dependenciesFile "\t\@cat `awk '{print \"".$workDirectoryName."\"\$\$1\".gv\"}' ".$graphVizUsed."` >> ".$workSubDirectoryName.$graphVizFileName."\n";
	}
	print $dependenciesFile "\t\@sort -u ".$workSubDirectoryName.$graphVizFileName." -o ".$workSubDirectoryName.$graphVizFileName."\n\n";
    }
}

exit;

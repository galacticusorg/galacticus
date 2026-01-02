#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use List::Uniq ':all';
use Data::Dumper;
use XML::Simple;
use Galacticus::Build::Directives;
use Fortran::Utils;
use List::ExtraUtils;
use Storable;
use Scalar::Util qw(reftype);

# Locate source files which have dependencies on modules.
# Andrew Benson (06-September-2016)

# Define the source directory
die "Usage: useDependencies.pl sourcedir"
    unless ( scalar(@ARGV) == 1 );
my $rootSourceDirectoryName = $ARGV[0];
# Specify work directory.
my $workDirectoryName       = $ENV{'BUILDPATH'}."/";
# Get an XML parser.  
my $xml                     = new XML::Simple();
# Load the file of directive locations.
my $locations               = -e $workDirectoryName."directiveLocations.xml" ? $xml->XMLin($workDirectoryName."directiveLocations.xml") : undef();
# List of external modules (which will be ignored for dependency analysis of the source code).
my @externalModules = ( "omp_lib", "hdf5", "h5tb", "h5lt", "h5global", "h5fortran_types", "fox_common", "fox_dom", "fox_wxml", "fox_utils", "mpi", "mpi_f08" );
# Modules that require a library to be linked. These are key-value pairs with the key being the module name, and the value the
# name of the required library.
my %moduleLibraries = (
    nearest_neighbors   => "ANN"           ,
    points_convex_hull  => "qhullcpp"      ,
    fftw3               => "fftw3"         ,
    fox_common          => "FoX_common"    ,
    fox_dom             => "FoX_dom"       ,
    fox_wxml            => "FoX_wxml"      ,
    fox_utils           => "FoX_utils"     ,
    hdf5                => "hdf5_fortran"  ,
    h5tb                => "hdf5hl_fortran",
    vectors             => "blas"          ,
    models_likelihoods  => "matheval"      ,
    input_parameters    => "matheval"      ,
    interface_gsl       => "gsl"           ,
    output_versioning   => "git2"
    );
# C includes that require a library to be linked. These are key-value pairs with the key being the include name, and the value the
# name of the required library.
my %includeLibraries = (
    crypt             => "crypt"
    );
# Parse the Makefile to find preprocessor macros that are set.
my @preprocessorDirectives;
my @makefiles = ( "Makefile", glob($ENV{'BUILDPATH'}."/Makefile*") );
foreach my $makefileName ( @makefiles ) { 
    my @conditionsStack = ( 1 );;
    open(my $makefile,$makefileName);
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
		if ( $line =~ m/^\s*FCFLAGS\s*\+??=/ && ! grep {! $_} @conditionsStack );
    }
    close($makefile);
}

# Extract any preprocessor directives specified via the GALACTICUS_FCFLAGS environment variable.
push(@preprocessorDirectives,map {$_ =~ m/\-D([0-9A-Z]+)/ ? $1 : ()} split(" ",$ENV{"GALACTICUS_FCFLAGS"}))
    if ( exists($ENV{"GALACTICUS_FCFLAGS"}) );
my $compiler = exists($ENV{'CCOMPILER'}) ? $ENV{'CCOMPILER'} : "gcc";
open(my $compilerDefs,$compiler." -dM -E - < /dev/null |");
while ( my $line = <$compilerDefs> ) {
    my @columns = split(" ",$line);
    push(@preprocessorDirectives,$columns[1]);
}

# Initialize structure to hold record of directives from each source file.
my $usesPerFile;
my $havePerFile = -e $workDirectoryName."Makefile_Use_Dependencies.blob";
my $updateTime;
if ( $havePerFile ) {
    $usesPerFile = retrieve($workDirectoryName."Makefile_Use_Dependencies.blob");
    $updateTime  = -M       $workDirectoryName."Makefile_Use_Dependencies.blob" ;
}
# Open the output dependencies file
open(my $dependenciesFile,">".$workDirectoryName."Makefile_Use_Dependencies");
# Build a list of source directories (including all sub-directories of the main directory).
my @sourceDirectoryNames  = ( $rootSourceDirectoryName."/source" );
if ( -e $sourceDirectoryNames[0] ) {
    opendir(my $sourceDirectory,$sourceDirectoryNames[0]);
    push(@sourceDirectoryNames,map {chomp($_);my $path = $sourceDirectoryNames[0]."/".$_; (-d $path && $_ !~ m/^\.+$/) ? $path : ()} readdir($sourceDirectory) );	 
    closedir($sourceDirectory);
}
push(@sourceDirectoryNames,$ENV{'BUILDPATH'}."/libgalacticus");

# Iterate over source directories.
my @sourceFilesToProcess;
foreach my $sourceDirectoryName ( @sourceDirectoryNames ) {
    # Extract the subdirectory name.
    (my $subDirectoryName = $sourceDirectoryName) =~ s/^($rootSourceDirectoryName\/source|$ENV{'BUILDPATH'})\/?//;
    # Find all source files to process.
    opendir(my $sourceDirectory,$sourceDirectoryName) 
	or die "useDependencies.pl: can not open the source directory: #!";
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

# Look for changes in source files. Force a rescan of all files if anything has changed.
my @fileIdentifiers;
my $forceRescan = 0;
# Iterate over source directories.
foreach my $sourceFile ( @sourceFilesToProcess ) {
    (my $fileIdentifier = $sourceFile->{'fullPathFileName'}) =~ s/\//_/g;
    push(@fileIdentifiers,$fileIdentifier);
}
# Check for new files.
foreach my $fileIdentifier ( @fileIdentifiers ) {
    unless ( exists($usesPerFile->{$fileIdentifier}) ) {
	$forceRescan = 1;
    }
}
# Check for removed files.
foreach my $fileIdentifier ( keys(%{$usesPerFile}) ) {
    unless ( grep {$_ eq $fileIdentifier} @fileIdentifiers ) {
	$forceRescan = 1;
    }
}

# Initialize list of modules needed for event hooks.
my @eventHookModules;
# Iterate over files to process.
foreach my $sourceFile ( @sourceFilesToProcess ) {
    # Push the main file onto the scan stack.
    my @fileNamesToProcess = ( $sourceFile->{'fullPathFileName'} );
    (my $fileIdentifier = $sourceFile->{'fullPathFileName'}) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;
    # Check if file is updated. If it is not, skip processing it. If it is, remove previous record of uses and rescan.
    my $rescan = 1;
    if ( $havePerFile && exists($usesPerFile->{$fileIdentifier}) ) {
	$rescan = 0
	    unless ( grep {-M $_ < $updateTime} &List::ExtraUtils::as_array($usesPerFile->{$fileIdentifier}->{'files'}) );
    }
    next
	unless ( $rescan || $forceRescan );
    delete($usesPerFile->{$fileIdentifier})
	if ( $havePerFile && exists($usesPerFile->{$fileIdentifier}) );
    push(@{$usesPerFile->{$fileIdentifier}->{'files'}},$sourceFile->{'fullPathFileName'});
    # Initialize records of modules used, provided, explicit dependencies, and any library dependencies.
    @{$usesPerFile->{$fileIdentifier}->{'modulesUsed'         }} = ();
    @{$usesPerFile->{$fileIdentifier}->{'dependenciesExplicit'}} = ();
    %{$usesPerFile->{$fileIdentifier}->{'modulesProvided'     }} = ();
    @{$usesPerFile->{$fileIdentifier}->{'submodules'          }} = ();
    %{$usesPerFile->{$fileIdentifier}->{'libraryDependencies' }} = ();
    # Extract lists of directives from this file which require special handling.
    my $directives;
    @{$directives->{$_}} = &Galacticus::Build::Directives::Extract_Directives($sourceFile->{'fullPathFileName'},$_)
	foreach ( "functionClass", "inputParameter", "enumeration", "eventHook", "eventHookStatic", "eventHookManager" );
    # Special handling for functionClass directives - add implementation files to the list of files to scan.
    if ( scalar(@{$directives->{'functionClass'}}) > 0 ) {
	foreach my $functionClass ( @{$directives->{'functionClass'}} ) {
	    &List::ExtraUtils::smart_push(\@fileNamesToProcess,$locations->{$functionClass->{'name'}}->{'file'});
	}
    }
    # Add dependence on functionClass module if necessary.
    push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName."function_classes.mod")
	if ( scalar(@{$directives->{'functionClass'}}) > 0 );
    # Accumulate dependencies for event hook modules.
    if ( scalar(@{$directives->{'eventHook'}}) > 0  ) {
	push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName."events_hooks.mod");
	foreach my $eventHook ( @{$directives->{'eventHook'}} ) {
	    next
		unless ( exists($eventHook->{'import'}) );
	    my @modules;
	    if ( exists($eventHook->{'import'}->{'module'}->{'name'}) ) {
		@modules = &List::ExtraUtils::as_array($eventHook->{'import'}->{'module'}                 );
	    } else {
		@modules = &List::ExtraUtils::hashList($eventHook->{'import'}->{'module'}, keyAs => "name");
	    }
	    foreach my $module ( @modules ) {
		push(@eventHookModules,$workDirectoryName.lc($module->{'name'}).".mod");
	    }
	}
    }
    if ( scalar(@{$directives->{'eventHookManager'}}) > 0  ) {
	my $workSubDirectoryName = $workDirectoryName.$sourceFile->{'subDirectoryName'}.($sourceFile->{'subDirectoryName'} eq "" ? "" : "/");
	($usesPerFile->{'eventHooksManager'}->{'objectFileName'} = $workSubDirectoryName.$sourceFile    ->{'fileName'}) =~ s/\.(f|f90|c|cpp)$/.o/i;
	$usesPerFile ->{'eventHooksManager'}->{'fileIdentifier'} =                       $fileIdentifier                                          ;
    }
    # Accumulate dependencies for static event hook modules.
    if ( scalar(@{$directives->{'eventHookStatic'}}) > 0  ) {
	foreach my $eventHookStatic ( @{$directives->{'eventHookStatic'}} ) {
	    foreach my $eventHookedStaticFile ( &List::ExtraUtils::as_array($locations->{$eventHookStatic->{'name'}}->{'file'}) ) {
		my $moduleName;
		open(my $file,$eventHookedStaticFile) or die "Can't open input file: $eventHookedStaticFile";
		while (my $line = <$file>) {
		    if ( $line =~ m/^\s*module\s+([a-zA-Z0-9_]+)/ ) {
			$moduleName = $1;
			last;
		    }
		}
		close($file);
		die("useDependencies.pl: unable to locate containing module for static event '".$eventHookStatic->{'name'}."' in file '".$eventHookedStaticFile."'")
		    unless ( defined($moduleName) );
		push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName.lc($moduleName).".mod");
	    }
	}
    }
    # Add dependence on input parameters module if necessary.
    push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName."input_parameters.mod")
	if ( scalar(@{$directives->{'functionClass'}}) > 0 ||  scalar(@{$directives->{'inputParameter'}}) > 0 );
    # Add dependence on error reporting module if necessary.
    push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName."error.mod")
	if (
	    (grep {exists($_->{'encodeFunction'}) && $_->{'encodeFunction'} eq "yes" && ! exists($_->{'errorValue'})} @{$directives->{'enumeration'}})
	    ||
	    (grep {exists($_->{'decodeFunction'}) && $_->{'decodeFunction'} eq "yes" && ! exists($_->{'errorValue'})} @{$directives->{'enumeration'}})
	);
    # Add dependence on enumerations module if necessary.
    push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName."enumerations.mod")
	if ( scalar(@{$directives->{'enumeration'}}) > 0 );
    # Find modules used in functionClass directives.
    foreach my $functionClass ( @{$directives->{'functionClass'}} ) {
	next 
	    unless ( exists($functionClass->{'method'}) );
	foreach my $method ( exists($functionClass->{'method'}->{'name'}) ? $functionClass->{'method'} : map {$functionClass->{'method'}->{$_}} keys(%{$functionClass->{'method'}}) ) {
	    next
		unless ( exists($method->{'modules'}) );
	    if ( reftype($method->{'modules'}) ) {
		# Array of modules.
		foreach my $moduleName ( map {lc($_)} sort(keys(%{$method->{'modules'}}))  ) {
		    push
			(
			 @{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},
			 $workDirectoryName.$moduleName.".mod"
			)
			unless ( grep {$_ eq $moduleName} @externalModules );
		}
	    } else {
		# Simple space-separated list of modules.
		foreach my $moduleName ( split(" ",lc($method->{'modules'})) ) {
		    push
			(
			 @{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},
			 $workDirectoryName.$moduleName.".mod"
			)
			unless ( grep {$_ eq $moduleName} @externalModules );
		}
	    }
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
	$usesPerFile->{$fileIdentifier}->{'libraryDependencies'}->{"stdc++"} = 1
	    if ( $fullPathFileName =~ m/\.cpp$/i );
	# Initialize preprocessor conditional compilation state and state stack.
	my @preprocessorConditionalsStack;
	my $conditionallyCompile = 1;
	my $inXML                = 0;
	my $inLaTeX              = 0;
	open(my $file,$fullPathFileName) or die "Can't open input file: $fullPathFileName";
	while (my $line = <$file>) {
	    # Detect leaving LaTeX and XML blocks.
	    $inXML   = 0
		if ( $line =~ m/^\s*!!\]/ );
	    $inLaTeX = 0
		if ( $line =~ m/^\s*!!\}/ );
	    next
		if ( $inXML || $inLaTeX );
	    if ( $line =~ m/^\s*\!;\s*([a-zA-Z0-9_]+)\s*$/ ) {
		$usesPerFile->{$fileIdentifier}->{'libraryDependencies'}->{$1} = 1;	
	    }
	    if ( $line =~ m/^\s*\#include\s+<([a-zA-Z0-9_]+)\.h>/ ) {
		my $includeFile = $1;
		$usesPerFile->{$fileIdentifier}->{'libraryDependencies'}->{$includeLibraries{lc($includeFile)}} = 1
		    if ( exists($includeLibraries{lc($includeFile)}) && $conditionallyCompile );
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
		if ( $line =~ m/^\s*(!\$\s)??\s*use\s*(::|\s)\s*([a-zA-Z0-9_]+)/i ) {
		    my $usedModule = $3;
		    # Add any library dependency for this module.
		    $usesPerFile->{$fileIdentifier}->{'libraryDependencies'}->{$moduleLibraries{lc($usedModule)}} = 1
			if ( exists($moduleLibraries{lc($usedModule)}) );
		    push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName.lc($usedModule).".mod")
			unless ( grep {$_ eq lc($usedModule)} @externalModules );
		}
		# Find any OpenMP parallel directives - these require a dependence on the event hooks module.
		push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName."events_filters.mod")
		    if ( $line =~ m/^\s*\!\$omp\s+parallel/i );
		# Find any OpenMP critical directives - these require a dependence on the OpenMP utilities data module.
		push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName."openmp_utilities_data.mod")
		    if ( $line =~ m/^\s*\!\$omp\s+critical\s*\([a-z0-9_]+\)/i );
		# Locate explicit dependencies.
		push(@{$usesPerFile->{$fileIdentifier}->{'dependenciesExplicit'}},split(" ",$2))
		    if ( $line =~ m/^\s*(\!|\/\/):\s*(.*)$/ );
		# Locate programs - this require a dependency on the "ISO_Varying_String" module for auto-generated allowed parameter names.
		push(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},$workDirectoryName."iso_varying_string.mod")
		    if ( $line =~ m/^\s*program\s/i );
		# Locate explicit dependencies.
		# Locate any modules provided by this file (we do not need to include an explicit dependence on any modules which
		# are self-provided).
		if ( $line =~ m/^\s*module\s+([a-zA-Z0-9_]+)/ ) {
		    my $moduleName = $1;
		    unless ( $moduleName eq "procedure" ) {
			$usesPerFile->{$fileIdentifier}->{'modulesProvided'}->{$moduleName.".mod"} = 1;
			if ( scalar(@{$directives->{'functionClass'}}) > 0 ) {
			    foreach my $functionClass ( @{$directives->{'functionClass'}} ) {
				foreach my $functionClassFileName ( &List::ExtraUtils::as_array($locations->{$functionClass->{'name'}}->{'file'}) ) {
				    my @functionClasses = &Galacticus::Build::Directives::Extract_Directives($functionClassFileName,$functionClass->{'name'});
				    my @submoduleNames = map {$_->{'name'}."_"} @functionClasses;
				    push(@{$usesPerFile->{$fileIdentifier}->{'submodules'}},@submoduleNames);
				}
			    }
			}
		    }
		}
		# Locate any submodules provided by this file.
		if ( $line =~ m/^\s*submodule\s*\(\s*([a-z0-9_]+)(:([a-z0-9_]+))??\s*\)\s*([a-zA-Z0-9_]+)/i ) {
		    # Construct the names of the module file and object file generated from this source file.
		    my $moduleFileName = lc($1).".mod";
		    my $submoduleName  = lc($4);
		    push(@{$usesPerFile->{$fileIdentifier}->{'submodulesProvided'}},{submoduleName => $submoduleName, moduleFileName => $moduleFileName});
		}
		# Locate included files and push them to the stack of files to process.
		if ( $line =~ m/^\s*include\s+(\'|\")([\w\.\-]+)(\'|\")/i ) {
		    my $preprocessedIncludedFile = $2;
		    (my $rawIncludedFile = $preprocessedIncludedFile) =~ s/\.inc$/.Inc/;
		    if ( -e $rootSourceDirectoryName."/source/".$rawIncludedFile ) {
			# A raw (unpreprocessed) matching file exists in the source directory. Use this as the dependency.
			push(@fileNamesToProcess,$rootSourceDirectoryName."/source/".$rawIncludedFile);
			push(@{$usesPerFile->{$fileIdentifier}->{'files'}},$rootSourceDirectoryName."/source/".$rawIncludedFile);
		    } elsif ( -e $workDirectoryName.$preprocessedIncludedFile ) {
			# A preprocessed matching file exists in the work directory. Use it as the dependency.
			push(@fileNamesToProcess,$workDirectoryName.$preprocessedIncludedFile);
			push(@{$usesPerFile->{$fileIdentifier}->{'files'}},$workDirectoryName.$preprocessedIncludedFile);
		    } elsif ( $preprocessedIncludedFile =~ m/(.*)\.type\.inc/ ) {
			# For old-style method include files, add a dependency on all files which contain the associated
			# directive.
			&List::ExtraUtils::smart_push(\@fileNamesToProcess,$locations->{$1}->{'file'});
			&List::ExtraUtils::smart_push($usesPerFile->{$fileIdentifier}->{'files'},$locations->{$1}->{'file'});
		    }
		}
	    }
	    # Detect entering LaTeX and XML blocks.
	    $inXML   = 1
		if ( $line =~ m/^\s*!!\[/ );
	    $inLaTeX = 1
		if ( $line =~ m/^\s*!!\{/ );
	}
	close($file);
    }
}
# Add any extra dependencies for the eventHooksManager.
push(@{$usesPerFile->{'eventHooksManager'}->{'modules'}},@eventHookModules);
@{$usesPerFile->{'eventHooksManager'}->{'modules'}} = uniq(sort(@{$usesPerFile->{'eventHooksManager'}->{'modules'}}));
my $eventHooksFileIdentifier = $usesPerFile->{'eventHooksManager'}->{'fileIdentifier'};
$usesPerFile->{$eventHooksFileIdentifier}->{'modulesUsed'} = []
    unless ( exists($usesPerFile->{$eventHooksFileIdentifier}->{'modulesUsed'}) );
push(@{$usesPerFile->{$eventHooksFileIdentifier}->{'modulesUsed'}},@{$usesPerFile->{'eventHooksManager'}->{'modules'}});
@{$usesPerFile->{$eventHooksFileIdentifier}->{'modulesUsed'}} = uniq(sort(@{$usesPerFile->{$eventHooksFileIdentifier}->{'modulesUsed'}}));
# Build a map of submodules associated with each module.
my %submodules;
## First include submodules generated from functionClasses.
foreach my $fileIdentifier ( keys(%{$usesPerFile}) ) {
    next
	unless ( exists($usesPerFile->{$fileIdentifier}->{'submodules'}) && scalar(@{$usesPerFile->{$fileIdentifier}->{'submodules'}}) > 0 );
    my @modulesProvided = keys(%{$usesPerFile->{$fileIdentifier}->{'modulesProvided'}});
    die("useDependencies.pl: submodules [".join(", ",@{$usesPerFile->{$fileIdentifier}->{'submodules'}})."] associated with multiple modules [".join(", ",@modulesProvided)."]")
	unless ( scalar(@modulesProvided) == 1 );
    $submodules{$workDirectoryName.lc($modulesProvided[0])} = $usesPerFile->{$fileIdentifier}->{'submodules'};
}
## Next include explictly defined submodules.
foreach my $fileIdentifier ( keys(%{$usesPerFile}) ) {
    next
	unless ( exists($usesPerFile->{$fileIdentifier}->{'submodulesProvided'}) );
    foreach my $submodule ( @{$usesPerFile->{$fileIdentifier}->{'submodulesProvided'}} ) {
	push(@{$submodules{$workDirectoryName.$submodule->{'moduleFileName'}}},$submodule->{'submoduleName'});
    }
}

# Iterate over files to generate make rules.
foreach my $sourceFile ( @sourceFilesToProcess ) {
    # Push the main file onto the scan stack.
    my @fileNamesToProcess = ( $sourceFile->{'fullPathFileName'} );
    (my $fileIdentifier = $sourceFile->{'fullPathFileName'}) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;
    # Construct name of associated object and dependency file.
    (my $objectFileName      = $sourceFile->{'fileName'}) =~ s/\.(f|f90|c|cpp)$/.o/i;
    (my $dependencyFileName  = $sourceFile->{'fileName'}) =~ s/\.(f|f90|c|cpp|inc)$/.d/i;
    # Construct name of work subdirectory.
    my $workSubDirectoryName = $workDirectoryName.$sourceFile->{'subDirectoryName'}.($sourceFile->{'subDirectoryName'} eq "" ? "" : "/");
    # Output library file rule.
    unless ( $objectFileName =~ m/\.Inc$/ ) {
	(my $libraryFileName = $objectFileName) =~ s/.o$/.fl/;
	print $dependenciesFile $workSubDirectoryName.$libraryFileName,":\n";
	my $redirect = ">";
	foreach ( sort(keys(%{$usesPerFile->{$fileIdentifier}->{'libraryDependencies'}})) ) {
	    print $dependenciesFile "\t\@echo ".$_." ".$redirect." ".$workSubDirectoryName.$libraryFileName."\n";
	    $redirect = ">>";
	}
    }    
    # For files which used any modules, generate dependency rules.
    if ( scalar(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}}) > 0 || scalar(@{$usesPerFile->{$fileIdentifier}->{'dependenciesExplicit'}}) > 0 ) {
	# Sort the list of modules used, and remove any duplicate entries.
	@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}} = grep {! exists($usesPerFile->{$fileIdentifier}->{'modulesProvided'}->{$_})} uniq(sort(@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}}));
	# Output the dependencies.
	my @submodulesUsed;
	foreach my $moduleUsed ( @{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}} ) {
	    next
		unless ( exists($submodules{$moduleUsed}) );
	    (my $moduleName = $moduleUsed) =~ s/\.mod$//;
	    push(@submodulesUsed,map {$moduleName."\@".lc($_).".smod"} @{$submodules{$moduleUsed}});
	}
	print $dependenciesFile $workSubDirectoryName.$objectFileName,": ".$workDirectoryName."utility.OpenMP.workaround.o ".join(" ",@{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}},@submodulesUsed,@{$usesPerFile->{$fileIdentifier}->{'dependenciesExplicit'}})." Makefile\n";
	# Generate rules for dependency files - we first append a ".d" to any module file names used.
	my @dependenciesUsed = map {$_ =~ m/\.mod$/ ? $_.".d" : $_} @{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}};	
	print $dependenciesFile $workSubDirectoryName.$dependencyFileName,": ".join(" ",@dependenciesUsed,map {(my $modifiedName = $_) =~ s/\.o$/.d/; $modifiedName} @{$usesPerFile->{$fileIdentifier}->{'dependenciesExplicit'}})."\n";
	print $dependenciesFile "\t\@echo ".$workSubDirectoryName.$objectFileName." > " .$workSubDirectoryName.$dependencyFileName."~\n";
	print $dependenciesFile "\t\@echo ".$workDirectoryName."utility.OpenMP.workaround.o >> ".$workSubDirectoryName.$dependencyFileName."~\n";
	foreach my $dependencyExplicit ( @{$usesPerFile->{$fileIdentifier}->{'dependenciesExplicit'}} ) {
	    (my $dependencyExplicitFileName = $dependencyExplicit) =~ s/\.o$/.d/;
	    print $dependenciesFile "\t\@cat ".($dependencyExplicit =~ m/\// ? "" : $workDirectoryName).$dependencyExplicitFileName." >> ".$workSubDirectoryName.$dependencyFileName."~\n";
	}
	print $dependenciesFile "\t\@cat ".$_." >> ".$workSubDirectoryName.$dependencyFileName."~\n"
	    foreach ( @dependenciesUsed );
	print $dependenciesFile "\t\@sort -u ".$workSubDirectoryName.$dependencyFileName."~ -o ".$workSubDirectoryName.$dependencyFileName."~\n";
	print $dependenciesFile "\t\@if cmp -s ".$workSubDirectoryName.$dependencyFileName." ".$workSubDirectoryName.$dependencyFileName."~ ; then \\\n";
	print $dependenciesFile "\t rm ".$workSubDirectoryName.$dependencyFileName."~ ; \\\n";
	print $dependenciesFile "\telse \\\n";
	print $dependenciesFile "\t mv ".$workSubDirectoryName.$dependencyFileName."~ ".$workSubDirectoryName.$dependencyFileName." ; \\\n";
	print $dependenciesFile "\tfi\n\n";
	# Create rules for making dependency trees with GraphViz.
	my @graphVizesUsed   = map {$_ =~ m/\.d$/ ? $_.".gv" : $_} @{$usesPerFile->{$fileIdentifier}->{'modulesUsed'}};		 
	my $graphVizFileName = $sourceFile->{'fileName'}.".gv";
	print $dependenciesFile $workSubDirectoryName.$graphVizFileName,": ".$workSubDirectoryName.$dependencyFileName." ".join(" ",@graphVizesUsed)."\n";
	print $dependenciesFile "\t\@echo \\\"".$sourceFile->{'subDirectoryName'}.$sourceFile->{'fileName'}."\\\" > ".$workSubDirectoryName.$graphVizFileName."\n";
	foreach my $dependencyExplicit ( @{$usesPerFile->{$fileIdentifier}->{'dependenciesExplicit'}} ) {
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
# Output the per file module use data.
store($usesPerFile,$workDirectoryName."Makefile_Use_Dependencies.blob");

exit;

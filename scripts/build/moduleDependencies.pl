#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Storable;
use List::ExtraUtils;
use List::Uniq ':all';
use Galacticus::Build::Directives;
use Galacticus::Build::SourceTree;
use Fortran::Utils;
use XML::Simple;
use Data::Dumper;

# Build Makefile rules which describe how to construct module files from source files.
# Andrew Benson (12-September-2016)

# Get arguments.
die "Usage: moduleDependencies.pl <sourceDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $sourceDirectoryName = $ARGV[0];
# Specify a work directory.
my $workDirectoryName = $ENV{'BUILDPATH'}."/";
# Get an XML parser.
my $xml                     = new XML::Simple();
# Load the file of directive locations.
my $locations               = -e $workDirectoryName."directiveLocations.xml" ? $xml->XMLin($workDirectoryName."directiveLocations.xml") : undef();
# Build a list of source directories to search.
my @sourceDirectoryDescriptors = ( {path => $sourceDirectoryName."/source", leaf => ""} );
if ( -e $sourceDirectoryDescriptors[0]->{'path'} ) {
    opendir(my $sourceDirectory,$sourceDirectoryDescriptors[0]->{'path'});
    push(@sourceDirectoryDescriptors,map {chomp($_); ( -d $sourceDirectoryDescriptors[0]->{'path'}."/".$_ && $_ !~ m/^\.+$/) ? {path => $sourceDirectoryDescriptors[0]->{'path'}."/".$_, leaf => $_."/"} : ()} readdir($sourceDirectory));
    system("mkdir -p ".$workDirectoryName."/".$_->{'leaf'})
	foreach ( @sourceDirectoryDescriptors );    
}
# Initialize data structure to hold per-file information.
my $modulesPerFile;
my $havePerFile = -e $workDirectoryName."Makefile_Module_Dependencies.blob";
my $forceRescan = 0;
my $updateTime;
if ( $havePerFile ) {
    $modulesPerFile = retrieve($workDirectoryName."Makefile_Module_Dependencies.blob");
    $updateTime     = -M       $workDirectoryName."Makefile_Module_Dependencies.blob" ;
    # Look for changes in source files. Force a rescan of all files if anything has changed.
    my @fileIdentifiers;
    # Iterate over source directories.
    foreach my $sourceDirectoryDescriptor ( @sourceDirectoryDescriptors ) {
	# Find all source files to process in this directory.
	opendir(my $sourceDirectory,$sourceDirectoryDescriptor->{'path'})
	    or die "moduleDepedencies.pl: can not open the source directory: #!";
	my @sourceFileNames = grep {$_ =~ m/\.(f|f90)$/i && $_ !~ m/^\.\#/} readdir($sourceDirectory);
	closedir($sourceDirectory);
	foreach my $sourceFileName ( @sourceFileNames ) {
	    (my $fileIdentifier = $sourceDirectoryDescriptor->{'path'}."/".$sourceFileName) =~ s/\//_/g;
	    push(@fileIdentifiers,$fileIdentifier);
	}
    }
    # Check for new files.
    foreach my $fileIdentifier ( @fileIdentifiers ) {
	unless ( exists($modulesPerFile->{$fileIdentifier}) ) {
	    $forceRescan = 1;
	}
    }
    # Check for removed files.
    foreach my $fileIdentifier ( keys(%{$modulesPerFile}) ) {
	unless ( grep {$_ eq $fileIdentifier} @fileIdentifiers ) {
	    $forceRescan = 1;
	}
    }
}

# Iterate over source directories.
foreach my $sourceDirectoryDescriptor ( @sourceDirectoryDescriptors ) {
    # Find all source files to process in this directory.
    opendir(my $sourceDirectory,$sourceDirectoryDescriptor->{'path'})
	or die "moduleDepedencies.pl: can not open the source directory: #!";
    my @sourceFileNames = grep {$_ =~ m/\.(f|f90)$/i && $_ !~ m/^\.\#/} readdir($sourceDirectory);
    closedir($sourceDirectory);
    # Iterate over source files.
    foreach my $sourceFileName ( @sourceFileNames ) {
	# Initialize list of modules provided, and a stack of file names to process.
	my @sourceFileNameStack = ( $sourceDirectoryDescriptor->{'path'}."/".$sourceFileName );
	(my $fileIdentifier = $sourceDirectoryDescriptor->{'path'}."/".$sourceFileName) =~ s/\//_/g;
	$fileIdentifier =~ s/^\._??//;
	# Check if file is updated. If it is not, skip processing it. If it is, remove previous record of uses and rescan.
	my $rescan = 1;
	if ( $havePerFile && exists($modulesPerFile->{$fileIdentifier}) ) {
	    $rescan = 0
		unless ( grep {-M $_ < $updateTime} &List::ExtraUtils::as_array($modulesPerFile->{$fileIdentifier}->{'files'}) );
	}
	next
	    unless ( $rescan || $forceRescan );
	delete($modulesPerFile->{$fileIdentifier})
    	    if ( $havePerFile && exists($modulesPerFile->{$fileIdentifier}) );
	# Extract list of functionClass directives from this file which require special handling.
	my @functionClasses = &Galacticus::Build::Directives::Extract_Directives($sourceDirectoryDescriptor->{'path'}."/".$sourceFileName,'functionClass');
	if ( scalar(@functionClasses) > 0 ) {
	    foreach my $functionClass ( @functionClasses ) {
		foreach my $functionClassFileName ( &List::ExtraUtils::as_array($locations->{$functionClass->{'name'}}->{'file'}) ) {
		    my $classTree  = &Galacticus::Build::SourceTree::ParseFile($functionClassFileName);
		    my $classNode  = $classTree;
		    my $classDepth = 0;
		    my %classes;
		    while ( $classNode ) {
			if ( $classNode->{'type'} eq $functionClass->{'name'} ) {
			    $classes{$classNode->{'directive'}->{'name'}}->{'name'} = $classNode->{'directive'}->{'name'};
			} elsif ( $classNode->{'type'} eq "type" ) {
			    if (
				$classNode->{'opener'} =~ m/^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*$functionClass->{'name'}([a-z0-9_]+)\s*$/i 
				&&
				defined($2)
				) {
				$classes{$functionClass->{'name'}.$4}->{'extends'} = $2;
			    }
			}
			$classNode = &Galacticus::Build::SourceTree::Walk_Tree($classNode,\$classDepth);
		    }
		    my @submodules;
		    foreach my $className ( keys(%classes) ) {
			(my $fileName     = $functionClassFileName) =~ s/^.*\/(.*)\.F90$/$1/;
			my $submoduleName =                                                                                   $classes{$className}->{'name'   }."_";
			my $extends       = $classes{$className}->{'extends'} eq $functionClass->{'name'}."Class" ? undef() : $classes{$className}->{'extends'}."_";
			my $submodule =
			{
			    name          => $submoduleName        ,
			    fileName      => $fileName             ,
			    extends       => $extends              ,
			    source        => $functionClassFileName,
			    functionClass => 1
			};
			push(@{$modulesPerFile->{$fileIdentifier}->{'submodules'}},$submodule);
		    }
		}
	    }
	}
	push(@{$modulesPerFile->{$fileIdentifier}->{'files'}},$sourceDirectoryDescriptor->{'path'}."/".$sourceFileName);
        $modulesPerFile->{$fileIdentifier}->{'sourceFileName'           } = $sourceFileName;
        $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'} = $sourceDirectoryDescriptor;
        @{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}} = ();
	while ( scalar(@sourceFileNameStack) > 0 ) {
	    # Open and read the source file.
	    my $sourceFilePathName = pop(@sourceFileNameStack);
	    my $inXML              = 0;
	    my $inLaTeX            = 0;
	    open(my $sourceFile,$sourceFilePathName) 
		or die "moduleDepedencies.pl: can not open input file: ".$sourceFilePathName;	    
	    while ( my $line = <$sourceFile> ) {
		# Detect leaving LaTeX and XML blocks.
		$inXML   = 0
		    if ( $line =~ m/^\s*!!\]/ );
		$inLaTeX = 0
		    if ( $line =~ m/^\s*!!\}/ );
		next
		    if ( $inXML || $inLaTeX );
		# Strip comments and trailing space.
		$line =~ s/\s*(!.*)?$//;
		# Locate any lines which use the "module" statement and extract the name of that module.
		if ( $line =~ m/^\s*module\s+([a-zA-Z0-9_]+)$/i ) {
		    # Construct the names of the module file and object file generated from this source file.
		    my $moduleFileName  = lc($1).".mod";
		    push(@{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}},$moduleFileName);
		}
		# Locate any lines which use the "submodule" statement and extract the name of that submodule.
		if ( $line =~ m/^\s*submodule\s*\(\s*([a-z0-9_]+)(:([a-z0-9_]+))??\s*\)\s*([a-zA-Z0-9_]+)$/i ) {
		    # Construct the names of the module file and object file generated from this source file.
		    my $moduleFileName = lc($1).".mod";
		    my $submoduleName  = lc($4);
		    (my $fileName      = $sourceFilePathName) =~ s/^.*\/([^\/]+)\.F90$/$1/;
		    push(@{$modulesPerFile->{$fileIdentifier}->{'submodulesProvided'}},{moduleName => $moduleFileName, submodule => {name => $submoduleName, fileName => $fileName, source => $sourceFilePathName, extends => undef(), functionClass => 0}});
		}
		# Detect include files and add to the stack of files to process.
		if ( $line =~ m/include\s+\'(\w+)\'/i ) {
		    my $includeFileName = $sourceDirectoryName."/".$1;
		    push(@sourceFileNameStack,$includeFileName);
		    push(@{$modulesPerFile->{$fileIdentifier}->{'files'}},$includeFileName);
		}
		# Detect entering LaTeX and XML blocks.
		$inXML   = 1
		    if ( $line =~ m/^\s*!!\[/ );
		$inLaTeX = 1
		    if ( $line =~ m/^\s*!!\{/ );
	    }
	    close($sourceFile);
	}
    }
}
# Build a map of submodules associated with each module.
my %submodules;
foreach my $fileIdentifier ( keys(%{$modulesPerFile}) ) {
    next
	unless ( exists($modulesPerFile->{$fileIdentifier}->{'submodules'}) && scalar(@{$modulesPerFile->{$fileIdentifier}->{'submodules'}}) > 0 );
    die("moduleDependencies.pl: submodules associated with multiple modules")
	unless ( scalar(@{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}}) == 1 );
    (my $moduleName = ${$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}}[0]) =~ s/\.mod$//;
    $submodules{lc($moduleName).".mod"} = $modulesPerFile->{$fileIdentifier}->{'submodules'};
}
foreach my $fileIdentifier ( keys(%{$modulesPerFile}) ) {
    next
	unless ( exists($modulesPerFile->{$fileIdentifier}->{'submodulesProvided'}) );
    foreach my $submodule ( @{$modulesPerFile->{$fileIdentifier}->{'submodulesProvided'}} ) {
	@{$submodules{$submodule->{'moduleName'}}} = ()
	    unless ( exists($submodules{$submodule->{'moduleName'}}) );
	push(@{$submodules{$submodule->{'moduleName'}}},$submodule->{'submodule'});	
    }
}

# Create the output Makefile.
open(my $makefile,">".$workDirectoryName."Makefile_Module_Dependencies");
foreach my $fileIdentifier ( sort(keys(%{$modulesPerFile})) ) {
    if ( exists($modulesPerFile->{$fileIdentifier}->{'modulesProvided'}) && scalar(@{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}}) > 0 ) {
	foreach my $moduleFileName ( @{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}} ) {
	    # Construct a Makefile rule for this module file.
	    (my $objectFileName = $modulesPerFile->{$fileIdentifier}->{'sourceFileName'}) =~ s/\.(f|f90)$/.o/i;
	    print $makefile $workDirectoryName.$moduleFileName.": ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName."\n";
	    print $makefile "\t\@if [ ! -f ".$workDirectoryName.$moduleFileName." ]; then \\\n";
	    print $makefile "\t  rm "      .$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName." ; \\\n";
	    print $makefile "\t  \$(MAKE) ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName." ; \\\n";
	    print $makefile "\tfi\n\n";
	    (my $moduleName = $moduleFileName) =~ s/\.mod$//;
	    foreach my $submodule ( @{$submodules{lc($moduleFileName)}} ) {
		 print $makefile $workDirectoryName.$moduleName."\@".lc($submodule->{'name'}).".smod: ".$workDirectoryName.$submodule->{'fileName'}.".o\n";
		 print $makefile "\t\@if [ ! -f ".$workDirectoryName.$moduleName."\@".lc($submodule->{'name'}).".smod ]; then \\\n";
		 print $makefile "\t  rm "      .$workDirectoryName.$submodule->{'fileName'}.".o ; \\\n";
		 print $makefile "\t  \$(MAKE) ".$workDirectoryName.$submodule->{'fileName'}.".o ; \\\n";
		 print $makefile "\tfi\n\n";
		 my $dependsOn;
		 if ( defined($submodule->{'extends'}) ) {
		     my @match = grep {$_->{'name'} eq $submodule->{'extends'}} @{$submodules{lc($moduleFileName)}};
		     if ( scalar(@match) == 0 ) {
			 print "no matching submodule found:\n";
			 print "\t'".$submodule->{'name'}."' extends '".$submodule->{'extends'}."'\n";
			 print "\tavailable submodules are:\n";
			 foreach my $availableSubmodule ( @{$submodules{lc($moduleFileName)}} ) {
			     print "\t\t'".$availableSubmodule->{'name'}."'\n";
			 }
			 die('ERROR: moduleDependencies.pl: no matching submodule found');
		     }
		     if ( scalar(@match) >  1 ) {
			 print "ERROR: class '".$submodule->{'name'}."' is an extension of class '".$submodule->{'extends'}."' which is multiply defined in files:\n";
			 foreach ( @match ) {
			     print "\t".$_->{'source'}."\n";
			 }
			 die('moduleDependencies.pl: multiple matching submodules found');
		     }
		     $dependsOn = $match[0]->{'fileName'}.".o";
		 } else {
		     $dependsOn = $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName;
		 }
		 print $makefile $workDirectoryName.$submodule->{'fileName'}.".o: ".$workDirectoryName.$dependsOn."\n\n";
		 # For submodules from functionClass objects the preprocessed file is made by preprocessing the parent functionClass type.
		 if ( $submodule->{'functionClass'} ) {		 
		     (my $preprocessedFileName = $modulesPerFile->{$fileIdentifier}->{'sourceFileName'}) =~ s/\.F90$/.p.F90/;
		     print $makefile $workDirectoryName.$submodule->{'fileName'}.".p.F90.up: ".$workDirectoryName.$modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$preprocessedFileName." ".$submodule->{'source'}."\n";
		     print $makefile "\t\@if [ ! -f ".$workDirectoryName.$submodule->{'fileName'}.".p.F90 ]; then \\\n";
		     print $makefile "\t  rm "      .$workDirectoryName.$modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$preprocessedFileName." ; \\\n";
		     print $makefile "\t  \$(MAKE) ".$workDirectoryName.$modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$preprocessedFileName." ; \\\n";
		     print $makefile "\tfi\n";
		     print $makefile $workDirectoryName.$submodule->{'fileName'}.".p.F90: ".$workDirectoryName.$submodule->{'fileName'}.".p.F90.up\n\n";
		 }
	    }
	    (my $objectDependencyFileName = $objectFileName) =~ s/.o$/.d/;
	    print $makefile $workDirectoryName.$moduleFileName.".d: ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectDependencyFileName." ".join(" ",map {$workDirectoryName.$_->{'fileName'}.".d"} @{$submodules{lc($moduleFileName)}})."\n";
	    print $makefile "\t\@echo ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName          ." > " .$workDirectoryName.$moduleFileName.".d~\n";
	    print $makefile "\t\@cat " .$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectDependencyFileName." >> ".$workDirectoryName.$moduleFileName.".d~\n";
	    foreach my $submodule ( @{$submodules{lc($moduleFileName)}} ) {
		print $makefile "\t\@echo ".$workDirectoryName.$submodule->{'fileName'}.".o >> ".$workDirectoryName.$moduleFileName.".d~\n";
		print $makefile "\t\@cat " .$workDirectoryName.$submodule->{'fileName'}.".d >> ".$workDirectoryName.$moduleFileName.".d~\n";
	    }
	    print $makefile "\t\@if cmp -s ".$workDirectoryName.$moduleFileName.".d ".$workDirectoryName.$moduleFileName.".d~ ; then \\\n";
	    print $makefile "\t rm ".$workDirectoryName.$moduleFileName.".d~ ; \\\n";
	    print $makefile "\telse \\\n";
	    print $makefile "\t mv ".$workDirectoryName.$moduleFileName.".d~ ".$workDirectoryName.$moduleFileName.".d ; \\\n";
	    print $makefile "\tfi\n\n";
	    # Create rule for making a *.mod.gv file which is used in building GraphViz descriptions of source
	    # file dependencies.
	    print $makefile $workDirectoryName.$moduleFileName.".gv: ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectDependencyFileName." ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$modulesPerFile->{$fileIdentifier}->{'sourceFileName'}.".gv\n";
	    print $makefile "\t\@echo ". $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$modulesPerFile->{$fileIdentifier}->{'sourceFileName'}." > ".$workDirectoryName.lc($moduleFileName).".gv\n";
	}
	# Create a rule for the module list file if needed.
	(my $modulesProvidedFileName = $modulesPerFile->{$fileIdentifier}->{'sourceFileName'})  =~ s/\.(f|f90)$/\.m/i;
	print $makefile $workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$modulesProvidedFileName.":\n";
	my $director = "> ";
	foreach ( @{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}} ) {
	    print $makefile "\t\@echo ".$workDirectoryName.$_." ".$director." ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$modulesProvidedFileName."\n";
	    $director = ">>";
	    if ( scalar(@{$submodules{lc($_)}}) > 0 ) {
		(my $smodName = $_) =~ s/\.mod$/.smod/;
		print $makefile "\t\@echo ".$workDirectoryName.$smodName." ".$director." ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$modulesProvidedFileName."\n";
	    }
	}
	print $makefile "\n";
	foreach ( @{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}} ) {
	    (my $moduleName = $_) =~ s/\.mod$//;
	    foreach my $submodule ( @{$submodules{lc($_)}} ) {
		print $makefile $workDirectoryName.$submodule->{'fileName'}.".m:\n";
		print $makefile "\t\@echo ".$workDirectoryName.$moduleName."\@".lc($submodule->{'name'}).".smod > ".$workDirectoryName.$submodule->{'fileName'}.".m\n\n";
	    }
	}
    }   
}
close($makefile);
# Output the per file module use data.
store($modulesPerFile,$workDirectoryName."Makefile_Module_Dependencies.blob");

exit 0;

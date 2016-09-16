#!/usr/bin/env perl
use strict;
use warnings;

# Build Makefile rules which describe how to construct module files from source files.
# Andrew Benson (12-September-2016)

# Get arguments.
die "Usage: moduleDependencies.pl <sourceDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $sourceDirectoryName = $ARGV[0];
# Specify a work directory.
my $workDirectoryName = $ENV{'BUILDPATH'}."/";
# Build a list of source directories to search.
my @sourceDirectoryDescriptors = ( {path => $sourceDirectoryName."/source", leaf => ""} );
if ( -e $sourceDirectoryDescriptors[0]->{'path'} ) {
    opendir(my $sourceDirectory,$sourceDirectoryDescriptors[0]->{'path'});
    push(@sourceDirectoryDescriptors,map {chomp($_); ( -d $sourceDirectoryDescriptors[0]->{'path'}."/".$_ && $_ !~ m/^\.+$/) ? {path => $sourceDirectoryDescriptors[0]->{'path'}."/".$_, leaf => $_."/"} : ()} readdir($sourceDirectory));
    system("mkdir -p ".$workDirectoryName."/".$_->{'leaf'})
	foreach ( @sourceDirectoryDescriptors );    
}
# Open the output Makefile.
open(my $makefile,">".$sourceDirectoryName."/".$workDirectoryName."Makefile_Module_Dependencies");
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
	my @modulesProvided;
	my @sourceFileNameStack = ( $sourceDirectoryDescriptor->{'path'}."/".$sourceFileName );
	while ( scalar(@sourceFileNameStack) > 0 ) {
	    # Open and read the source file.
	    my $sourceFilePathName = pop(@sourceFileNameStack);
	    open(my $sourceFile,$sourceFilePathName) 
		or die "moduleDepedencies.pl: can not open input file: ".$sourceFilePathName;	    
	    while ( my $line = <$sourceFile> ) {
		# Strip comments and trailing space.
		$line =~ s/\s*(!.*)?$//;
		# Locate any lines which use the "module" statement and extract the name of that module.
		if ( $line !~ m/^\s*module\s+procedure/i && $line =~ m/^\s*module\s+([a-zA-Z0-9_]+)/i ) {
		    # Construct the names of the module file and object file generated from this source file.
		    my $moduleFileName  = lc($1).".mod";
		    (my $objectFileName = $sourceFileName) =~ s/\.(f|f90)$/.o/i;
		    push(@modulesProvided,$moduleFileName);
		    # Construct a Makefile rule for this module file.
		    print $makefile $workDirectoryName.$moduleFileName.": ".$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$objectFileName."\n";
		    print $makefile "\t\@if [ ! -f ".$workDirectoryName.$moduleFileName." ]; then \\\n";
		    print $makefile "\t  rm "      .$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$objectFileName." ; \\\n";
		    print $makefile "\t  \$(MAKE) ".$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$objectFileName." ; \\\n";
		    print $makefile "\tfi\n\n";
		    (my $objectDependencyFileName = $objectFileName) =~ s/.o$/.d/;
		    print $makefile $workDirectoryName.$moduleFileName.".d: ".$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$objectDependencyFileName."\n";
		    print $makefile "\t\@echo ".$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$objectFileName          ." > " .$workDirectoryName.$moduleFileName.".d\n";
		    print $makefile "\t\@cat " .$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$objectDependencyFileName." >> ".$workDirectoryName.$moduleFileName.".d\n";		    
		    # Create rule for making a *.mod.gv file which is used in building GraphViz descriptions of source
		    # file dependencies.
		    print $makefile $workDirectoryName.$moduleFileName.".gv: ".$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$objectDependencyFileName." ".$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$sourceFileName.".gv\n";
		    print $makefile "\t\@echo ".$sourceDirectoryDescriptor->{'leaf'}.$sourceFileName." > ".$workDirectoryName.lc($moduleFileName).".gv\n";

		}
		# Detect include files and add to the stack of files to process.
		push(@sourceFileNameStack,$sourceDirectoryName."/".$1)
		    if ( $line =~ m/include\s+\'(\w+)\'/i );
	    }
	    close($sourceFile);
	    # Create a rule for the module list file if needed.
	    if ( scalar(@modulesProvided) > 0 ) {
		(my $modulesProvidedFileName = $sourceFileName)  =~ s/\.(f|f90)$/\.m/i;
		print $makefile $workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$modulesProvidedFileName.":\n";
		my $director = "> ";
		foreach ( @modulesProvided ) {
		    print $makefile "\t\@echo ".$workDirectoryName.$_." ".$director." ".$workDirectoryName.$sourceDirectoryDescriptor->{'leaf'}.$modulesProvidedFileName."\n";
		    $director = ">>";
		}
		print $makefile "\n";
	    }
	}
    }
}

exit 0;

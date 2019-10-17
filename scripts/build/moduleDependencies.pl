#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Storable;
use List::ExtraUtils;
use XML::Simple;

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
# Initialize data structure to hold per-file information.
my $modulesPerFile;
my $havePerFile = -e $workDirectoryName."Makefile_Module_Dependencies.blob";
my $updateTime;
if ( $havePerFile ) {
    $modulesPerFile = retrieve($workDirectoryName."Makefile_Module_Dependencies.blob");
    $updateTime     = -M       $workDirectoryName."Makefile_Module_Dependencies.blob" ;
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
	    unless ( $rescan );
	delete($modulesPerFile->{$fileIdentifier})
    	    if ( $havePerFile && exists($modulesPerFile->{$fileIdentifier}) );
	push(@{$modulesPerFile->{$fileIdentifier}->{'files'}},$sourceDirectoryDescriptor->{'path'}."/".$sourceFileName);
        $modulesPerFile->{$fileIdentifier}->{'sourceFileName'           } = $sourceFileName;
        $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'} = $sourceDirectoryDescriptor;
        @{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}} = ();
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
		    push(@{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}},$moduleFileName);
		}
		# Detect include files and add to the stack of files to process.
		if ( $line =~ m/include\s+\'(\w+)\'/i ) {
		    my $includeFileName = $sourceDirectoryName."/".$1;
		    push(@sourceFileNameStack,$includeFileName);
		    push(@{$modulesPerFile->{$fileIdentifier}->{'files'}},$includeFileName);
		}
	    }
	    close($sourceFile);
	}
    }
}
# Create the output Makefile.
open(my $makefile,">".$sourceDirectoryName."/".$workDirectoryName."Makefile_Module_Dependencies");
foreach my $fileIdentifier ( sort(keys(%{$modulesPerFile})) ) {
    if ( scalar(@{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}}) > 0 ) {
	foreach my $moduleFileName ( @{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}} ) {
	    # Construct a Makefile rule for this module file.
	    (my $objectFileName = $modulesPerFile->{$fileIdentifier}->{'sourceFileName'}) =~ s/\.(f|f90)$/.o/i;
	    print $makefile $workDirectoryName.$moduleFileName.": ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName."\n";
	    print $makefile "\t\@if [ ! -f ".$workDirectoryName.$moduleFileName." ]; then \\\n";
	    print $makefile "\t  rm "      .$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName." ; \\\n";
	    print $makefile "\t  \$(MAKE) ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName." ; \\\n";
	    print $makefile "\tfi\n\n";
	    (my $objectDependencyFileName = $objectFileName) =~ s/.o$/.d/;
	    print $makefile $workDirectoryName.$moduleFileName.".d: ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectDependencyFileName."\n";
	    print $makefile "\t\@echo ".$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectFileName          ." > " .$workDirectoryName.$moduleFileName.".d~\n";
	    print $makefile "\t\@cat " .$workDirectoryName. $modulesPerFile->{$fileIdentifier}->{'sourceDirectoryDescriptor'}->{'leaf'}.$objectDependencyFileName." >> ".$workDirectoryName.$moduleFileName.".d~\n";
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
	}
	print $makefile "\n";
    }
}
close($makefile);
# Output the per file module use data.
store($modulesPerFile,$workDirectoryName."Makefile_Module_Dependencies.blob");

# Create an output XML file.
my $modulesByFile;
foreach my $fileIdentifier ( sort(keys(%{$modulesPerFile})) ) {
    $modulesByFile->{$modulesPerFile->{$fileIdentifier}->{'sourceFileName'}} = $modulesPerFile->{$fileIdentifier}->{'modulesProvided'}
        if ( scalar(@{$modulesPerFile->{$fileIdentifier}->{'modulesProvided'}}) > 0 );
}
my $xml = new XML::Simple();
open(my $xmlOutput,">".$sourceDirectoryName."/".$workDirectoryName."moduleLocations.xml");
print $xmlOutput $xml->XMLout($modulesByFile, rootName => "moduleLocations");
close($xmlOutput);

exit 0;

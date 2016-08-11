#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use Text::Template 'fill_in_string';

# Locate files which contain programs and append to a list of executables.
# Andrew Benson (08-August-2016)

# Define the source directory
die "Usage: findExecutables.pl <galacticusDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $galacticusDirectoryName = $ARGV[0];

# Determine source directory.
my $sourceDirectoryName = $galacticusDirectoryName."/source";

# Determine work directory.
$make::workDirectoryName = $ENV{'BUILDPATH'}."/";

# Open an output file
open(my $outputFile,">".$galacticusDirectoryName."/".$ENV{'BUILDPATH'}."/Makefile_All_Execs");

# Open the source directory and iterate over all files.
my @executableNames;
opendir(my $sourceDirectory,$sourceDirectoryName) 
    or die "findExecutables.pl: can not open the source directory: #!";
while ( my $fileName = readdir($sourceDirectory) ) {
    # Skip temporary files.
    next
	if ( $fileName =~ m/^\.\#/ );
    # Skip non-source files.
    next
	unless ( $fileName =~ m/\.[fF](90)?$/ );
    # Parse the file.
    my $fileFullName       = $sourceDirectoryName."/".$fileName;
    my $excludeFromMakeAll = 0;
    open(my $file,$fileFullName) 
	or die "findExecuitables.pl: can not open input file: #!";	    
    while (my $line = <$file>) {
	# Record if file should be excluded from "make all".
	$excludeFromMakeAll = 1
	    if ( $line =~ m/^\s*!\/\s+exclude/ );
	# Test for line indicating the start of a "program".
	if ( $line =~ m/^\s*program\s/i ) {
	    ($make::fileNameRoot = $fileName) =~ s/\.[fF](90)?t?$//;
	    push(@executableNames,$make::fileNameRoot.".exe")
		unless ( $excludeFromMakeAll == 1 );
	    print $outputFile fill_in_string(<<'MAKE', PACKAGE => 'make');
{$fileNameRoot}.exe: {$workDirectoryName.$fileNameRoot}.o {$workDirectoryName.$fileNameRoot}.d $(MAKE_DEPS)
	$(FCCOMPILER) `cat {$workDirectoryName.$fileNameRoot}.d` -o {$fileNameRoot}.exe $(FCFLAGS) `./scripts/build/Library_Dependencies.pl {$fileNameRoot}.exe $(FCFLAGS)`
	./scripts/build/executableSize.pl {$fileNameRoot}.exe {$workDirectoryName.$fileNameRoot}.size
	./scripts/build/parameterDependencies.pl `pwd` {$fileNameRoot}.exe

MAKE
	}
    }
    close($file);
}
closedir($sourceDirectory);

# Generate a Makefile rule to build all executables.
print $outputFile "all_exes = ".join(" ",@executableNames)."\n"
    if ( scalar(@executableNames) > 0 );
close($outputFile);

exit;


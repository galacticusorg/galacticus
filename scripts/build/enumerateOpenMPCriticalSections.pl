#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use XML::Simple;
use Fortran::Utils;
use File::Changes;
use Data::Dumper;

# Locate all OpenMP critical sections, and build an enumeration of them for use in source code instrumentation.
# Andrew Benson (29-November-2016)

# Get command line arguments.
die "Usage: enumerateOpenMPCriticalSections.pl <sourceDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $sourceDirectoryName = $ARGV[0];
# Dictionary of all named critical sections.
my %criticalSectionNames;
# Open the source directorys and scan
opendir(my $sourceDirectory,$sourceDirectoryName."/source") 
    or die "Can't open the source directory: #!";
my @sourceFileNames = grep {$_ =~ m/\.f(90)??$/i} readdir($sourceDirectory);
closedir($sourceDirectory);
foreach my $sourceFileName ( @sourceFileNames ) {
    # Open and read the file.
    my $sourceFilePathName = $sourceDirectoryName."/source/".$sourceFileName;
    open(my $sourceFile,$sourceFilePathName)
	or die "Can't open input file: ".$sourceFilePathName;
    until ( eof($sourceFile) ) {
	# Get the next line of this file.
	&Fortran::Utils::Get_Fortran_Line($sourceFile,my $rawLine,my $processedLine,my $bufferedComments);

	if ( $processedLine =~ m/^\s*\!\$omp\s+critical\s*\(([a-z0-9_]+)\)/i ) {
	    my $criticalSectionName = lc($1);
	    ++$criticalSectionNames{$criticalSectionName};
	}
    }
}
# Build a data structure describing the critical sections.
my $descriptor    ;
my $id         = 0;
foreach ( sort(keys(%criticalSectionNames)) ) {
    push(@{$descriptor->{'critical'}},
	 {
	     name      => $_,
	     frequency => $criticalSectionNames{$_},
	     id        => ++$id
	 }
	);
}
# Output the critical section enumeration to XML, replacing a preexisting file if and only if the new file differs.
my $xmlOutput = new XML::Simple (RootName=>"criticalSections");
open(my $enumerationFile,">".$ENV{'BUILDPATH'}."/openMPCriticalSections.xml.tmp");
print $enumerationFile $xmlOutput->XMLout($descriptor);
close($enumerationFile);
&File::Changes::Update($ENV{'BUILDPATH'}."/openMPCriticalSections.xml",$ENV{'BUILDPATH'}."/openMPCriticalSections.xml.tmp");
# Build an include file that specifies the number of named sections.
open(my $countFile,">".$ENV{'BUILDPATH'}."/openMPCriticalSections.count.inc.tmp");
print $countFile "! Number of named OpenMP critical sections in the source.\n";
print $countFile "integer, public, parameter :: criticalSectionCount=".scalar(keys(%criticalSectionNames))."\n";
close($countFile);
&File::Changes::Update($ENV{'BUILDPATH'}."/openMPCriticalSections.count.inc",$ENV{'BUILDPATH'}."/openMPCriticalSections.count.inc.tmp");
# Build an include file that enumerates critical section names.
open(my $enumerateFile,">".$ENV{'BUILDPATH'}."/openMPCriticalSections.enumerate.inc.tmp");
print $enumerateFile "type(varying_string), dimension(criticalSectionCount) :: criticalSectionNames\n";
print $enumerateFile "criticalSectionNames=[ &\n";
print $enumerateFile " & var_str('".join("'), &\n & var_str('",sort(keys(%criticalSectionNames)))."') &\n";
print $enumerateFile " & ]\n";
close($enumerateFile);
&File::Changes::Update($ENV{'BUILDPATH'}."/openMPCriticalSections.enumerate.inc",$ENV{'BUILDPATH'}."/openMPCriticalSections.enumerate.inc.tmp");

exit;

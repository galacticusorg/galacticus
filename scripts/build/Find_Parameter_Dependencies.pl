#!/usr/bin/env perl
use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use List::MoreUtils qw{ any };

# Scan source files for input parameter definitions for a given executable.
# Andrew Benson (18-October-2011)

# Get the name of the executable to find parameters for.
die("Usage: Find_Parameter_Dependencies.pl <executable>") unless ( scalar(@ARGV) == 1 );
my $executable = $ARGV[0];

# Build a hash of dependencies.
my $dependencies;
(my $dependencyFile = "work/build/".$executable) =~ s/\.exe/\.d/;
open(iHndl,$dependencyFile);
while ( my $line = <iHndl> ) {
    if ( $line =~ m/\.\/work\/build\/(.+\.o)$/ ) {
	my $sourceFile = $1;
	push(@{$dependencies->{'objectFiles'}},$sourceFile);
    }
}
close(iHndl);

# Open the source diretory, finding F90 and cpp files.
my %parametersListed;
opendir(sDir,"source");
while ( my $fileName = readdir(sDir) ) {
    if ( $fileName =~ m/\.F90$/ || $fileName =~ m/\.cpp$/ ) {
	(my $objectFile = $fileName) =~ s/\.[a-zA-Z0-9]+$/\.o/;
	# Open the file and scan for parameters.
	my $xmlBuffer;
	my @fileStack = ( "source/".$fileName );
	while ( scalar(@fileStack) > 0 ) {
	    my $fileToProcess = shift(@fileStack);
	    open(sFile,$fileToProcess);
	    while ( my $line = <sFile> ) {
		# Detect included files and push onto the stack.
		if ( exists {map { $_ => 1 } @{$dependencies->{'objectFiles'}}}->{$objectFile} ) {
		    if ( $line =~ m/^\s*include\s*[\'\"](.*)[\'\"]/ ) {
			push(@fileStack,"work/build/".$1);
		    }
		}
		# Search for XML.
		if ( $line =~ m/^\s*(\!|\/\/)@(.*)/ ) {
		    my $xmlLine = $2;
		    $xmlBuffer = "" if ( $xmlLine =~ m/^\s*<inputParameter>\s*$/ );
		    $xmlBuffer .= $xmlLine;
		    if ( $xmlLine =~ m/^\s*<\/inputParameter>\s*$/ ) {
			# Parse the XML.
			my $xml = new XML::Simple;
			my $inputParameter = $xml->XMLin($xmlBuffer);
			if ( exists {map { $_ => 1 } @{$dependencies->{'objectFiles'}}}->{$objectFile} ) {
			    unless ( exists {map { $_ => 1 } @{$dependencies->{'parameters' }}}->{$inputParameter->{'name'}} ) {
				push(@{$dependencies->{'parameters'}},$inputParameter->{'name'});
			    }
			}
		    }
		}
	    }
	    close(sFile);
	}
    }
}
close(sDir);

# Check for any extra parameters that are not discoverable from the source files.
(my $extraFile = "./work/build/".$executable) =~ s/\.exe/.parameters.extra.xml/;
if ( -e $extraFile ) {
    my $xml = new XML::Simple;
    my $extraParameters = $xml->XMLin($extraFile);
    push(@{$dependencies->{'parameters'}},@{$extraParameters->{'parameter'}});
}

# Store in a structure for output.
(my $outputFile = $executable) =~ s/\.exe/.parameters.xml/;
my $outputList;
@{$outputList->{'parameter'}} = @{$dependencies->{'parameters'}} if ( exists($dependencies->{'parameters'}) );
my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
open(oHndl,">./work/build/".$outputFile);
print oHndl $xmlOutput->XMLout($outputList);
close(oHndl);

exit;

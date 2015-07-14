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
use XML::Simple;
use Data::Dumper;
use List::MoreUtils qw{ any };
require Galacticus::Doc::Parameters;

# Scan source files for input parameter definitions for a given executable.
# Andrew Benson (18-October-2011)

# Get the name of the executable to find parameters for.
die("Usage: Find_Parameter_Dependencies.pl <sourceDirectory> <executable>")
    unless ( scalar(@ARGV) == 2 );
my $sourceDirectory = $ARGV[0];
my $executable      = $ARGV[1];

# Build a hash of dependencies.
my $dependencies;
(my $dependencyFile = $ENV{'BUILDPATH'}."/".$executable) =~ s/\.exe/\.d/;
open(iHndl,$dependencyFile);
while ( my $line = <iHndl> ) {
    if ( $line =~ m/$ENV{'BUILDPATH'}\/(.+\.o)$/ ) {
	my $sourceFile = $1;
	push(@{$dependencies->{'objectFiles'}},$sourceFile);
    }
}
close(iHndl);

# Open the source diretory, finding F90 and cpp files.
my %parametersListed;
opendir(sDir,"source");
while ( my $fileName = readdir(sDir) ) {
    # Skip junk files.
    next
	if ( $fileName =~ m/^\.\#/ );
    # For Fortran files, check for a ".p" parameter file in the build directory.
    if ( $fileName =~ m/\.F90$/) {
	(my $parameterFile = $fileName) =~ s/\.F90$/.p/;
	if ( -e $ENV{'BUILDPATH'}."/".$parameterFile ) {
	    open(my $parameterFile,$ENV{'BUILDPATH'}."/".$parameterFile);
	    while ( my $parameterName = <$parameterFile> ) {
		chomp($parameterName);
		unless ( grep {$_ eq $parameterName} @{$dependencies->{'parameters' }} ) {
		    push(@{$dependencies->{'parameters'}},$parameterName);
		}
	    }
	    close($parameterFile);
	}
    }
    # Process Fortran and C++ files.
    if ( $fileName =~ m/\.(F90|cpp)$/ ) {
	(my $objectFile = $fileName) =~ s/\.(F90|cpp)$/\.o/;
	# Open the file and scan for parameters.
	my $xmlBuffer;
	my @fileStack = ( "source/".$fileName );
	while ( scalar(@fileStack) > 0 ) {
	    my $fileToProcess = shift(@fileStack);
	    open(sFile,$fileToProcess);
	    while ( my $line = <sFile> ) {
		# Detect included files and push onto the stack.
		if ( $line =~ m/^\s*include\s*[\'\"](.*)[\'\"]/ ) {
		    if ( exists {map { $_ => 1 } @{$dependencies->{'objectFiles'}}}->{$objectFile} ) {
			my $fileName = $1;
			push(@fileStack,$ENV{'BUILDPATH'}."/".$1)
			    unless ( $fileName eq "fftw3.f03" );
		    }
		}
		# Search for XML.		
		if ( $line =~ m/^\s*(!|\/\/)\@(.*)/ ) {
		    my $xmlLine = $2;
		    $xmlBuffer = "" if ( $xmlLine =~ m/^\s*<inputParameter>\s*$/ );
		    $xmlBuffer .= $xmlLine;
		    if ( $xmlLine =~ m/^\s*<\/inputParameter>\s*$/ ) {
			# Parse the XML.
			my $xml = new XML::Simple;
			my $inputParameter = $xml->XMLin($xmlBuffer);
			$inputParameter->{'name'} = "regEx:".&Parameters::ExpandRegEx($inputParameter->{'regEx'},$sourceDirectory)
			    if ( exists($inputParameter->{'regEx'}) );
			if ( grep {$_ eq $objectFile}  @{$dependencies->{'objectFiles'}} ) {
			    unless ( grep {$_ eq $inputParameter->{'name'}} @{$dependencies->{'parameters' }} ) {
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
(my $extraFile = $ENV{'BUILDPATH'}."/".$executable) =~ s/\.exe/.parameters.extra.xml/;
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
my $output = $xmlOutput->XMLout($outputList);
# Escape square brackets in the output to that they get correctly parsed by FoX.
$output =~ s/\[/&#x005B;/g;
$output =~ s/\]/&#x005D;/g;
open(oHndl,">".$ENV{'BUILDPATH'}."/".$outputFile);
print oHndl $output;
close(oHndl);

exit;

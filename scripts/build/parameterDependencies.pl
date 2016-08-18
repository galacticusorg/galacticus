#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use XML::Simple;
use Data::Dumper;
use File::Slurp;
use List::Uniq ':all';
use List::MoreUtils qw{ any };
use Galacticus::Doc::Parameters;
use Galacticus::Build::Directives;

# Scan source files for input parameter definitions for a given executable.
# Andrew Benson (18-October-2011)

# Get the name of the executable to find parameters for.
die("Usage: parameterDependencies.pl <sourceDirectory> <executable>")
    unless ( scalar(@ARGV) == 2 );
my $sourceDirectoryName = $ARGV[0];
my $executableName      = $ARGV[1];

# Include files to exclude from parameter search.
my @includeFilesExcluded = ( "fftw3.f03" );

# Build a list of object file dependencies.
(my $dependencyFileName = $ENV{'BUILDPATH'}."/".$executableName) =~ s/\.exe/\.d/;
my @objectFiles = map { $_ =~ /^$ENV{'BUILDPATH'}\/(.+\.o)$/ ? $1 : () } read_file($dependencyFileName, chomp => 1);

# Open the source diretory, finding F90 and cpp files.
my $output;
opendir(my $sourceDirectory,$sourceDirectoryName);
while ( my $fileName = readdir($sourceDirectory) ) {
    # Skip junk files.
    next
	if ( $fileName =~ m/^\.\#/ );
    # Skip non-F90, non-cpp files
    next
	unless ( $fileName =~ m/\.(F90|cpp)$/ );
    # Find corresponding object file name.    
    (my $objectFileName = $fileName) =~ s/\.(F90|cpp)$/\.o/;
    # Skip non-dependency files.
    next
	unless ( grep {$_ eq $objectFileName} @objectFiles );
    # For Fortran files, check for a ".p" parameter file in the build directory. These files are created by
    # Galacticus::Build::SourceTree::Process::FunctionClass.
    if ( $fileName =~ m/\.F90$/) {
	(my $rootFileName = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/./;
	push(@{$output->{'parameters'}},read_file($rootFileName."p", chomp => 1, err_mode => 'quiet'))
	    if ( grep {$_ eq $rootFileName."o"} @objectFiles );
    }
    # Process Fortran and C++ files.
    # Create a stack of files to process.
    my @fileStack = ( $sourceDirectoryName."/".$fileName );
    while ( scalar(@fileStack) > 0 ) {
	my $fileToProcess = shift(@fileStack);
	# Find "include" lines in the file, extract the name of the included file, filter out any include files which are to be
	# excluded from parameter search, and push those remaining onto the file stack.
	push
	    (
	     @fileStack,
	     map 
	      {($_ =~ m/^\s*include\s*[\'\"](.*)[\'\"]/ && ! grep {$1 eq $_} @includeFilesExcluded) ? $ENV{'BUILDPATH'}."/".$1 : ()} 
	      read_file($fileToProcess, chomp => 1)
	    );
	# Find all "inputParameter" directives, extract names from them, and push to the list of parameters.
	push
	    (
	     @{$output->{'parameters'}},
	     map 
	      {exists($_->{'regEx'}) ? "regEx:".&Galacticus::Doc::Parameters::ExpandRegEx($_->{'regEx'},$sourceDirectoryName) : $_->{'name'}}
	      &Galacticus::Build::Directives::Extract_Directives($fileToProcess,"inputParameter",comment => qr/^\s*(!|\/\/)\@/) 
	    );
    }
}
close($sourceDirectory);

# Remove duplicated parameters.
@{$output->{'parameters'}} = uniq({sort => 1}, @{$output->{'parameters'}});

# Serialize to XML.
my $xmlOutput        = new XML::Simple (NoAttr=>1, RootName=>"parameters");
my $outputSerialized = $xmlOutput->XMLout($output);

# Escape square brackets in the output to that they get correctly parsed by FoX.
$outputSerialized =~ s/\[/&#x005B;/g;
$outputSerialized =~ s/\]/&#x005D;/g;

# Output the results.
(my $outputFileName = $executableName) =~ s/\.exe/.parameters.xml/;
open(oHndl,">".$ENV{'BUILDPATH'}."/".$outputFileName);
print oHndl $outputSerialized;
close(oHndl);

exit;

#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Data::Dumper;
use File::Slurp;
use List::Uniq ':all';
use List::MoreUtils qw{ any };
use List::ExtraUtils;
use Galacticus::Build::Directives;
use Storable;

# Scan source files for input parameter definitions for a given executable.
# Andrew Benson (18-October-2011)

# Get the name of the executable to find parameters for.
die("Usage: parameterDependencies.pl <sourceDirectory> <target>")
    unless ( scalar(@ARGV) == 2 );
my $sourceDirectoryName = $ARGV[0];
my $targetName          = $ARGV[1];

# Include files to exclude from parameter search.
my @includeFilesExcluded = ( "fftw3.f03" );

# Build a list of object file dependencies.
(my $dependencyFileName = $ENV{'BUILDPATH'}."/".$targetName) =~ s/\.(exe|o)$/\.d/;
my @objectFiles = map { $_ =~ /^$ENV{'BUILDPATH'}\/(.+\.o)$/ ? $1 : () } read_file($dependencyFileName, chomp => 1);

# Initialize structure to hold record of parameters from each source file.
(my $blobFileName = $targetName) =~ s/\.(exe|o)$/.blob/;
my $parametersPerFile;
my $havePerFile = -e $ENV{'BUILDPATH'}."/".$blobFileName;
my $updateTime;
if ( $havePerFile ) {
    $parametersPerFile = retrieve($ENV{'BUILDPATH'}."/".$blobFileName);
    $updateTime        = -M       $ENV{'BUILDPATH'}."/".$blobFileName ;
}

# Open the source diretory, finding F90 and cpp files.
opendir(my $sourceDirectory,$sourceDirectoryName."/source");
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
    # Create a stack of files to process.
    my @fileStack;
    (my $rootFileName = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/./;
    if ( $fileName =~ m/\.F90$/ && -e $rootFileName."p.F90" ) {
	@fileStack = ( $rootFileName."p.F90" );
    } else {
	@fileStack = ( $sourceDirectoryName."/source/".$fileName ); 
    }
    # Check if file is updated. If it is not, skip processing it. If it is, remove previous record of uses and rescan.
    (my $fileIdentifier = $fileStack[0]) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;
    my $rescan = 1;
    if ( $havePerFile && exists($parametersPerFile->{$fileIdentifier}) ) {
	$rescan = 0
	    unless ( grep {-M $_ < $updateTime} &List::ExtraUtils::as_array($parametersPerFile->{$fileIdentifier}->{'files'}) );
    }
    if ( $rescan ) {
	delete($parametersPerFile->{$fileIdentifier})
    	    if ( $havePerFile && exists($parametersPerFile->{$fileIdentifier}) );
	push(@{$parametersPerFile->{$fileIdentifier}->{'files'}},$fileStack[0]);
	# For Fortran files, check for a ".p" parameter file in the build directory. These files are created by
	# Galacticus::Build::SourceTree::Process::FunctionClass.
	if ( $fileName =~ m/\.F90$/ ) {
	    push(@{$parametersPerFile->{$fileIdentifier}->{'parameter'}},read_file($rootFileName."p", chomp => 1, err_mode => 'quiet'))
		if ( grep {$_ eq $objectFileName} @objectFiles );
	}    
	# Process the file stack.
	while ( scalar(@fileStack) > 0 ) {
	    my $fileToProcess = shift(@fileStack);
	    # Find "include" lines in the file, extract the name of the included file, filter out any include files which are to be
	    # excluded from parameter search, and push those remaining onto the file stack.
	    my @includedFiles = 
		map 
	    {($_ =~ m/^\s*include\s*[\'\"](.*)[\'\"]/ && ! grep {$1 eq $_} @includeFilesExcluded) ? $ENV{'BUILDPATH'}."/".$1 : ()} 
	    read_file($fileToProcess, chomp => 1);
	    push(@fileStack                                         ,@includedFiles);
	    push(@{$parametersPerFile->{$fileIdentifier}->{'files'}},@includedFiles);
	    # Find all "inputParameter" directives, extract names from them, and push to the list of parameters.
	    push
		(
		 @{$parametersPerFile->{$fileIdentifier}->{'parameter'}},
		 map
		 {$_->{'name'}}
		 &Galacticus::Build::Directives::Extract_Directives($fileToProcess,"inputParameter") 
		);
	    # Find all "objectBuilder" directives with non-standard parameter names, extract names from them, and push to the list of parameters.
	    push
		(
		 @{$parametersPerFile->{$fileIdentifier}->{'parameter'}},
		 map 
		 {
		     exists($_->{'parameterName'}) && $_->{'parameterName'} ne $_->{'class'}
		     ?
		     (
                            $_->{'parameterName'}
		     )
		     :
		     ()
		 }
		 &Galacticus::Build::Directives::Extract_Directives($fileToProcess,"objectBuilder")
		);
	}
    }
}
close($sourceDirectory);

# Reduce over files.
my $output;
@{$output->{'parameter'}} = map {$parametersPerFile->{$_}->{'parameter'}} keys(%{$parametersPerFile});

# Remove duplicated parameters.
@{$output->{'parameter'}} = uniq({sort => 1}, @{$output->{'parameter'}});

# Output the results.
(my $outputFileName = $targetName) =~ s/\.(exe|o)$/.parameters.F90/;
open(my $outputFile,">".$ENV{'BUILDPATH'}."/".$outputFileName);
print $outputFile "subroutine knownParameterNames(names)\n";
print $outputFile "  use ISO_Varying_String\n";
print $outputFile "  implicit none\n";
print $outputFile "  type(varying_string), dimension(:), allocatable, intent(inout) :: names \n";
print $outputFile "  allocate(names(".scalar(@{$output->{'parameter'}})."))\n";
for(my $i=0;$i<scalar(@{$output->{'parameter'}});++$i) {
    print $outputFile "  names(".($i+1).")='".$output->{'parameter'}->[$i]."'\n";
}
print $outputFile "end subroutine knownParameterNames\n";
close($outputFile);
# Output the per file module use data.
store($parametersPerFile,$ENV{'BUILDPATH'}."/".$blobFileName);

exit;

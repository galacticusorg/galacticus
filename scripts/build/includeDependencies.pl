#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use List::Uniq ':all';
use File::Slurp;
use Data::Dumper;

# Construct Makefile rules for source files which have dependencies on include files.
# Andrew Benson (08-August-2016)

# Get arguments.
die "Usage: includeDependencies.pl <galacticusDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $galacticusDirectoryName = $ARGV[0];

# Specify the source directory.
my $sourceDirectoryName = $galacticusDirectoryName."/source";

# Build a list of pre-existing C include file paths. This includes all paths searched by GCC, any files specified in environ, and
# the source directory.
my @cIncludeDirectoryNames =
    (
     $sourceDirectoryName, 
     map {$_ =~ m/^\s*(\/\S+)/ ? $1 : ()} split(/\n/,`echo | cpp -Wp,-v 2>&1`),     
    );
foreach my $environmentVariable ( "GALACTICUS_CFLAGS", "GALACTICUS_CPPFLAGS" ) {
    push(@cIncludeDirectoryNames,map {$_ =~ m/^\-I(.*)/ ? $1 : ()} split(" ",$ENV{$environmentVariable}))
	if ( exists($ENV{$environmentVariable}) );
}

# Open an output file
open(my $makeFile,">".$ENV{'BUILDPATH'}."/Makefile_Include_Deps");

# Initialize a list of files on which Makefile_Use_Deps will depend.
my @dependencyFileNames;

# Scan the source directory.
opendir(my $sourceDirectory,$sourceDirectoryName)
    or die "includeDependencies.pl: can not open the source directory: #!";
while ( my $fileName = readdir($sourceDirectory) ) {
    # Skip temporary files.
    next
	if ( $fileName =~ m/^\.\#/ );
    # Skip non-source files.
    next
	unless
	(
	 $fileName =~ m/\.f(90)??$/i # Fortran and Fortran-90
	 ||
	 $fileName =~ m/\.inc$/i     # Fortran include files
	 ||
	 $fileName =~ m/\.c(pp)??$/  # C and C++ source files
	 ||
	 $fileName =~ m/\.h$/        # C/C++ include files
	);
    # Parse the file. We match include lines, extracting the name of the included file into a list of hashes. Each hash also holds
    # an "automatic" element which is true unless the "NO_USES" qualifier appears in the include statement line (which instructs
    # us to ignore that include file for the purposes of finding module "use" dependencies).
    my $fileFullName  = $sourceDirectoryName."/".$fileName;
    my @includedFiles = 
	map
         {
	     (
	      $_ =~ /^\s*#??include\s*[<'"]((.+)\.(inc|h)\d*)['">](\s*\!\s*NO_USES)?/
	      &&
	      (
	       $3 eq "inc"
	       ||
	       ($3 eq "h" && ! -e $sourceDirectoryName."/".$1 && ! grep {-e $_."/".$1} @cIncludeDirectoryNames) # C-include files must not be pre-existing in source or system directories,
	       )
	     )
          ?
             {fileName => $1, automatic => !defined($4)} 
	  :
             ()
         }
        read_file($fileFullName, chomp => 1);
    # Output the dependencies.
    (my $objectFileName = $fileName)  =~ s/\.[^\.]+$/.o/;
    print $makeFile $ENV{'BUILDPATH'}."/".$objectFileName.": ".join(" ",sort(map {$ENV{'BUILDPATH'}."/".$_->{'fileName'}} @includedFiles))."\n\n"
	if ( scalar(@includedFiles) > 0 );
    # Generate a list of include files on which Makefile_Use_Deps will depend. This is the unpreprocessed include file in the
    # source directory if such exists, or the raw include file in the build directory otherwise.
    foreach my $includedFile ( @includedFiles ) {
	(my $fileNameUnpreprocessed = $includedFile->{'fileName'}) =~ s/\.inc$/.Inc/;
	push(@dependencyFileNames,-e $sourceDirectoryName."/".$fileNameUnpreprocessed ? $fileNameUnpreprocessed : $includedFile->{'fileName'})
	   if ( $includedFile->{'fileName'} =~ m/\.inc$/ && $includedFile->{'automatic'} );
    }
}
closedir($sourceDirectory);
# Add dependencies for Makefile_Use_Deps, such that it will be regenerated after any automatically-built include files are
# created. This will ensure that any module "uses" in the automatically-generated include files will be discovered.
print $makeFile "\n".$ENV{'BUILDPATH'}."/Makefile_Use_Deps: ".join(" ",map {$ENV{'BUILDPATH'}."/".$_} sort(@dependencyFileNames))."\n";
close($makeFile);

exit;

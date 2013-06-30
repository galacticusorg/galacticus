#!/usr/bin/env perl
use strict;
use warnings;
use Sort::Topological qw(toposort);

# Output linker options to link required libraries for building an executable.
# Andrew Benson (24-May-2012)

# Get the name of the executable being built.
die("Usage: Library_Dependencies.pl <executable>")
    unless ( scalar(@ARGV) == 1 );
my $executable = $ARGV[0];

# Specify library dependencies.
my %dependencies = 
    (
     hdf5_fortran  => [ "hdf5"                             ],
     hdf5          => [ "z"                                ],
     fgsl_gfortran => [ "gsl"                              ],
     gsl           => [ "gslcblas"                         ],
     FoX_dom       => [ "FoX_fsys", "FoX_utils", "FoX_sax" ],
     FoX_sax       => [ "FoX_common"                       ]
    );

# Library order dependencies for static linking.
my %staticLinkDependencies =
    (
     hdf5          => [ "z"                                ],
     hdf5_fortran  => [ "hdf5"                             ],
     fgsl_gfortran => [ "gsl"                              ],
     gsl           => [ "gslcblas"                         ],
     FoX_dom       => [ "FoX_fsys", "FoX_utils", "FoX_sax" ],
     FoX_sax       => [ "FoX_common"                       ],
     FoX_common    => [ "FoX_fsys"                         ]
    );

# Initialize a hash of required libraries.
my %libraries;

# Open the file of dependencies for the executable.
(my $dependencyFile = "./work/build/".$executable) =~ s/\.exe$/\.d/;
die("Library_Dependencies.pl: dependency file is missing")
    unless ( -e $dependencyFile );
open(iHndl,$dependencyFile);
while ( my $objectFile = <iHndl> ) {
    # Construct the name of the library dependency file.
    (my $libraryDependenciesFile = $objectFile) =~ s/\.o$/\.fl/;
    chomp($libraryDependenciesFile);
    # Test whether the library dependency file exists.
    if ( -e $libraryDependenciesFile ) {
	open(lHndl,$libraryDependenciesFile);
	while ( my $library = <lHndl> ) {
	    chomp($library);
	    $libraries{$library} += 1;
	}
	close(lHndl);
    }
}
close(iHndl);

# Augment the list of libraries with any dependencies.
my $libraryCount = 0;
while ( scalar(keys(%libraries)) != $libraryCount) {
    $libraryCount = scalar(keys(%libraries));
    foreach my $library ( keys(%libraries) ) {
	if ( exists($dependencies{$library}) ) {
	    foreach ( @{$dependencies{$library}} ) {$libraries{$_} +=1};
	}
    }
}

# Perform a topological sort on libraries to ensure they are in the correct order for static linking.
my @unsortedLibraries = keys(%libraries);
sub staticLinkDependency { @{$staticLinkDependencies{$_[0]} || []}; }
my @sortedLibraries = toposort(\&staticLinkDependency, \@unsortedLibraries);

# Write the linker options to standard output.
print join(" ",map {"-l".$_} @sortedLibraries)."\n";

exit;

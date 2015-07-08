#!/usr/bin/env perl
use strict;
use warnings;
use Sort::Topological qw(toposort);

# Output linker options to link required libraries for building an executable.
# Andrew Benson (24-May-2012)

# Get the name of the executable being built.
die("Usage: Library_Dependencies.pl <executable>")
    unless ( scalar(@ARGV) >= 1 );
my $executable = shift(@ARGV);
my @compilerOptions = @ARGV;

# Specify library dependencies.
my %dependencies = 
    (
     hdf5hl_fortran => [ "hdf5_hl"                          ],
     hdf5_hl        => [ "hdf5"                             ],
     hdf5_fortran   => [ "hdf5"                             ],
     hdf5           => [ "z"                                ],
     fgsl_gfortran  => [ "gsl"                              ],
     gsl            => [ "gslcblas"                         ],
     FoX_dom        => [ "FoX_fsys", "FoX_utils", "FoX_sax" ],
     FoX_sax        => [ "FoX_common"                       ],
     yepLibrary     => [ "yeppp"                            ],
     yepCore        => [ "yeppp"                            ],
     yepMath        => [ "yeppp"                            ]
    );

# Library order dependencies for static linking.
my %staticLinkDependencies =
    (
     hdf5           => [ "z", "dl"                          ],
     hdf5_hl        => [ "hdf5"                             ],
     hdf5_fortran   => [ "hdf5"                             ],
     hdf5hl_fortran => [ "hdf5_hl"                          ],
     fgsl_gfortran  => [ "gsl"                              ],
     gsl            => [ "gslcblas"                         ],
     FoX_dom        => [ "FoX_fsys", "FoX_utils", "FoX_sax" ],
     FoX_sax        => [ "FoX_common"                       ],
     FoX_common     => [ "FoX_fsys"                         ],
     YEPLibrary     => [ "yeppp"                            ],
     YEPCore        => [ "yeppp"                            ],
     YEPMath        => [ "yeppp"                            ]
    );

# Detect static linking.
my $isStatic = 0;
$isStatic = 1
    if ( grep {$_ eq "-static"} @compilerOptions);
push(@{$dependencies{'hdf5'}},"dl")
    if ( $isStatic == 1 );

# Detect if pthread is already included.
my $pthreadIncluded = 0;
$pthreadIncluded = 1
    if ( grep {$_ eq "-lpthread"} @compilerOptions);

# Initialize a hash of required libraries.
my %libraries;

# Open the file of dependencies for the executable.
(my $dependencyFile = $ENV{'BUILDPATH'}."/".$executable) =~ s/\.exe$/\.d/;
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

# Add static link options.
my $staticOptions = "";
$staticOptions = " -Wl,--whole-archive -lpthread -Wl,--no-whole-archive"
    if ( $isStatic == 1 && $pthreadIncluded == 0 );

# Write the linker options to standard output.
print join(" ",map {"-l".$_} @sortedLibraries).$staticOptions;

# If we are linking BLAS, cause GFortran to use the external BLAS library.
print " -fexternal-blas"
    if ( grep{$_ eq "blas"} @sortedLibraries );

# Check glibc version.
my $glibcVersionMajor;
my $glibcVersionMinor;
open(my $pipe,"ldd --version|");
while ( my $line = <$pipe> ) {
    if ( $line =~ m/^ldd \(GNU libc\) (\d+)\.(\d+)/ ) {
	$glibcVersionMajor = $1;
	$glibcVersionMinor = $2;
    }
}
close($pipe);

# Determine if we need to link the realtime library.
print " -lrt"
    if ( $glibcVersionMajor < 2 || $glibcVersionMajor == 2 && $glibcVersionMinor <= 16 );

# Write newline.
print "\n";

exit;

#!/usr/bin/env perl
use strict;
use warnings;
use Sort::Topological qw(toposort);
use Data::Dumper;

# Output linker options to link required libraries for building an executable.
# Andrew Benson (24-May-2012)

# Get the name of the executable being built.
die("Usage: libraryDependencies.pl <executable>")
    unless ( scalar(@ARGV) >= 1 );
my $executable      = shift(@ARGV);
my @compilerOptions =       @ARGV ;
# Specify library dependencies. Each key corresponds to a library name. The associated value is a list of other library names on
# which the key library depends.
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
     FoX_utils      => [ "FoX_wxml"                         ],
     yepLibrary     => [ "yeppp"                            ],
     yepCore        => [ "yeppp"                            ],
     yepMath        => [ "yeppp"                            ]
    );
# Library order dependencies for static linking. Each key specifies a library name. The associated value is a list of libraries
# before which the key library must appear in the link command when static linking is used.
my %staticLinkDependencies =
    (
     hdf5           => [ "z", "dl"                                      ],
     hdf5_hl        => [ "hdf5"                                         ],
     hdf5_fortran   => [ "hdf5"                                         ],
     hdf5hl_fortran => [ "hdf5_hl"                                      ],
     fgsl_gfortran  => [ "gsl"                                          ],
     gsl            => [ "gslcblas"                                     ],
     FoX_dom        => [ "FoX_fsys", "FoX_utils", "FoX_sax", "FoX_wxml" ],
     FoX_sax        => [ "FoX_common"                                   ],
     FoX_wxml       => [ "FoX_utils"                                    ],
     FoX_common     => [ "FoX_fsys"                                     ],
     YEPLibrary     => [ "yeppp"                                        ],
     YEPCore        => [ "yeppp"                                        ],
     YEPMath        => [ "yeppp"                                        ]
    );
# Detect static linking.
my $isStatic = grep {$_ eq "-static"} @compilerOptions;
# Add explicit dependency on libdl (for dynamic linking support) for the HDF5 library when linking statically. This appears to be
# necessary as HDF5 relies on some functions in libdl, but it is not linked by default when static linking is used.
push(@{$dependencies{'hdf5'}},"dl")
    if ( $isStatic );
# Detect if libpthread is already included.
my $pthreadIncluded = grep {$_ eq "-lpthread"} @compilerOptions;
# Initialize a hash of required libraries.
my %libraries;
# Open the file of dependencies for the executable.
(my $dependencyFileName = $ENV{'BUILDPATH'}."/".$executable) =~ s/\.exe$/\.d/;
die("libraryDependencies.pl: dependency file is missing")
    unless ( -e $dependencyFileName );
open(my $dependencyFile,$dependencyFileName);
chomp(my @objectFileNames = <$dependencyFile>);
close($dependencyFile);
# Read library dependency files for each object file.
foreach my $objectFileName ( map {$_ =~ s/\.o$/\.fl/; -e $_ ? $_ : ()} @objectFileNames ) {
    open(my $libraryDependenciesFile,$objectFileName);
    chomp(my @libraries = <$libraryDependenciesFile>);
    close($libraryDependenciesFile);
    $libraries{$_} += 1
	foreach ( @libraries );
}
# Augment the list of libraries with any dependencies. Iterate this operation until no new libraries are added (which we detect by
# testing whether the number of keys in the %libraries hash has changed after the operation).
my $libraryCount = 0;
while ( scalar(keys(%libraries)) != $libraryCount) {
    $libraryCount = scalar(keys(%libraries));
    foreach my $library ( keys(%libraries) ) {
	if ( exists($dependencies{$library}) ) {
	    foreach ( @{$dependencies{$library}} ) {$libraries{$_} +=1};
	}
    }
}
# Remove YEPPP library if not used.
delete($libraries{'yeppp'})
    if ( exists($libraries{'yeppp'}) && ! grep {$_ eq "-DYEPPP"} @compilerOptions );
# Perform a topological sort on libraries to ensure they are in the correct order for static linking.
my @sortedLibraries = toposort(sub { @{$staticLinkDependencies{$_[0]} || []}; }, [keys(%libraries)]);
# Add static link options.
my $staticOptions = ( $isStatic && ! $pthreadIncluded ) ? " -Wl,--whole-archive -lpthread -Wl,--no-whole-archive" : "";
# Determine glibc version.
my $linkLibRT = 1;
open(my $lddPipe,"ldd --version|");
while ( my $line = <$lddPipe> ) {
    if ( $line =~ m/^ldd \(GNU libc\) (\d+)\.(\d+)/ ) {
	$linkLibRT = $1 < 2 || ( $1 == 2 && $2 <= 16 );
    }
}
close($lddPipe);
# Write the linker options to standard output. If we are linking BLAS, cause GFortran to use the external BLAS library. Also link
# the realtime library for sufficiently early glibc's.
print 
    join(" ",map {"-l".$_} @sortedLibraries).
    $staticOptions.
    ((grep {$_ eq "blas"} @sortedLibraries) ? " -fexternal-blas" : "").
    ($linkLibRT                             ? " -lrt"            : "").
    "\n";
exit;

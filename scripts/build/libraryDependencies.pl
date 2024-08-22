#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Config;
use Sort::Topo;
use Data::Dumper;
use File::Which qw(which);

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
     gsl            => [ "gslcblas"                         ],
     FoX_dom        => [ "FoX_fsys", "FoX_utils", "FoX_sax" ],
     FoX_sax        => [ "FoX_common"                       ],
     FoX_utils      => [ "FoX_wxml"                         ],
     qhullcpp       => [ "qhull_r", "stdc++"                ]
    );
# Library order dependencies for static linking. Each key specifies a library name. The associated value is a list of libraries
# before which the key library must appear in the link command when static linking is used.
my %staticLinkDependencies =
    (
     hdf5           => [ "z", "dl"                                      ],
     hdf5_hl        => [ "hdf5"                                         ],
     hdf5_fortran   => [ "hdf5"                                         ],
     hdf5hl_fortran => [ "hdf5_hl"                                      ],
     gsl            => [ "gslcblas"                                     ],
     FoX_dom        => [ "FoX_fsys", "FoX_utils", "FoX_sax", "FoX_wxml" ],
     FoX_sax        => [ "FoX_common"                                   ],
     FoX_wxml       => [ "FoX_utils"                                    ],
     FoX_common     => [ "FoX_fsys"                                     ],
     qhullcpp       => [ "stdc++"                                       ]
    );
# Find default preprocessor directives.
my @preprocessorDirectives;
my $compiler = exists($ENV{'CCOMPILER'}) ? $ENV{'CCOMPILER'} : "gcc";
open(my $compilerDefs,$compiler." -dM -E - < /dev/null |");
while ( my $line = <$compilerDefs> ) {
    my @columns = split(" ",$line);
    push(@preprocessorDirectives,$columns[1]);
}
# Detect MacOS.
my $isMacOS = grep {$_ eq "__APPLE__"} @preprocessorDirectives;
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
# Verify that BUILDPATH is set.
die('libraryDependencies.pl: "BUILDPATH" environment variable is not set')
    unless ( exists($ENV{'BUILDPATH'}) );
# Open the file of dependencies for the executable.
(my $mainDependencyFileName = $ENV{'BUILDPATH'}."/".$executable) =~ s/\.(exe|o)$/\.d/;
die("libraryDependencies.pl: dependency file is missing")
    unless ( -e $mainDependencyFileName );
my @dependencyFileNames = ( $mainDependencyFileName );
push(@dependencyFileNames,$ENV{'BUILDPATH'}."/libgalacticus_classes.d")
    if ( $executable eq "libgalacticus.o" );
my @objectFileNames;
foreach my $dependencyFileName ( @dependencyFileNames ) {
    open(my $dependencyFile,$dependencyFileName);
    chomp(my @thisObjectFileNames = <$dependencyFile>);
    close($dependencyFile);
    push(@objectFileNames,@thisObjectFileNames);
}
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
    foreach my $library ( sort(keys(%libraries)) ) {
	if ( exists($dependencies{$library}) ) {
	    foreach ( @{$dependencies{$library}} ) {$libraries{$_} +=1};
	}
    }
}
# Remove FFTW3 library if not used.
delete($libraries{'fftw3'})
    if ( exists($libraries{'fftw3'}) && ! grep {$_ eq "-DFFTW3AVAIL"} @compilerOptions );
# Remove ANN library if not used.
delete($libraries{'ANN'})
    if ( exists($libraries{'ANN'}) && ! grep {$_ eq "-DANNAVAIL"} @compilerOptions );
# Remove qhull library if not used.
delete($libraries{'qhullcpp'})
    if ( exists($libraries{'qhullcpp'}) && ! grep {$_ eq "-DQHULLAVAIL"} @compilerOptions );
delete($libraries{'qhull_r' })
    if ( exists($libraries{'qhull_r' }) && ! grep {$_ eq "-DQHULLAVAIL"} @compilerOptions );
# Remove libmatheval if not used.
delete($libraries{'matheval'})
    if ( exists($libraries{'matheval'}) && ! grep {$_ eq "-DMATHEVALAVAIL"} @compilerOptions );
# Remove libgit2 if not used.
delete($libraries{'git2'})
    if ( exists($libraries{'git2'}) && ! grep {$_ eq "-DGIT2AVAIL"} @compilerOptions );
# Perform a topological sort on libraries to ensure they are in the correct order for static linking.
my @libraryNames = sort(keys(%libraries));
my @sortedLibraries = &Sort::Topo::sort(\@libraryNames,\%staticLinkDependencies);
# For static linking qhull libraries must be renamed.
@sortedLibraries = map {$_ eq "qhull_r" ? "qhullstatic_r" : $_} @sortedLibraries
    if ( $isStatic );
# Add static link options.
my $staticOptions = ( $isStatic && ! $pthreadIncluded ) ? " -Wl,--whole-archive -lpthread -Wl,--no-whole-archive" : "";
# Add OS specific options.
my $osOptions = ($Config{osname} eq "darwin") ? " -Wl,-commons,use_dylibs" : "";
# Determine glibc version.
my $ldd = which("ldd");
my $linkLibRT = 1;
if ( defined($ldd) ) {
    open(my $lddPipe,"ldd --version|");
    while ( my $line = <$lddPipe> ) {
	if ( $line =~ m/^ldd \(GNU libc\) (\d+)\.(\d+)/ ) {
	    $linkLibRT = $1 < 2 || ( $1 == 2 && $2 <= 16 );
	}
    }
    close($lddPipe);
} else {
    $linkLibRT = 0;
}

# Write the linker options to standard output. If we are linking BLAS, cause GFortran to use the external BLAS library. Also link
# the realtime library for sufficiently early glibc's.
print 
    join(" ",map {"-l".$_} @sortedLibraries).
    $staticOptions.
    $osOptions.
    ((grep {$_ eq "-DDEBUGGING"} @compilerOptions) ? " -rdynamic"       : "").
    ((grep {$_ eq "blas"       } @sortedLibraries) ? " -fexternal-blas" : "").
    ($linkLibRT                                    ? " -lrt"            : "").
    "\n";
exit;

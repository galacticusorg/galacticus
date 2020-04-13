# Contains a Perl module which implements utility functions for source code generation.

package Galacticus::Build::SourceTree::Process::Utils;
use strict;
use warnings;
use utf8;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Exporter qw(import);

our @EXPORT_OK = qw(performIO);

# Determine if thread safe IO is required.
my $threadSafeIO = exists($ENV{'GALACTICUS_FCFLAGS'}) && grep {$_ eq "-DTHREADSAFEIO"} split(" ",$ENV{'GALACTICUS_FCFLAGS'});

sub performIO {
    # Perform IO. This wrapper function is used to allow optional thread control to be included to ensure thread safety.
    my $ioCode = shift();
    my $newIOCode = "";
    # <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
    #  <description>Internal file I/O in gfortran can be non-thread safe.</description>
    # </workaround>
    $newIOCode .= "!\$omp critical(gfortranInternalIO)\n"
	if ( $threadSafeIO );
    $newIOCode .= $ioCode;
    $newIOCode .= "!\$omp end critical(gfortranInternalIO)\n"
	if ( $threadSafeIO );
    return $newIOCode;
}

1;

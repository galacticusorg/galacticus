# Contains a Perl module which implements locking Fortran IO to make it thread-safe.

package Galacticus::Build::SourceTree::Process::ThreadSafeIO;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Fortran::Utils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks       {'threadSafeIO'} = \&Lock_IO;
$Galacticus::Build::SourceTree::Hooks::processDependencies{'threadSafeIO'} = [ "metaPropertyDatabase", "stateStore", "functionClass", "stateStorable" ];

sub Lock_IO {
    ## <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
    ##  <description>Internal file I/O in gfortran can be non-thread safe.</description>
    ## </workaround>
    
    # Get the tree.
    my $tree = shift();
    # Check if debugging is required.
    my $lock = 0;
    if ( exists($ENV{'GALACTICUS_FCFLAGS'}) ) {
	$lock = 1
	    if ( grep {$_ eq "-DTHREADSAFEIO"} split(" ",$ENV{'GALACTICUS_FCFLAGS'}) );
    }
    return
	unless ( $lock );
    # Detect if we should add reporting on lock acquire/release.
    my $report = exists($ENV{'GALACTICUS_REPORT_THREADSAFEIO'}) && $ENV{'GALACTICUS_REPORT_THREADSAFEIO'} eq "yes";
    # Do not apply thread locking inside error reporting functions. These can be triggered from anywhere in the code - including
    # within another I/O critical section which would then lead to deadlocking. Since at this point we are in an error condition
    # anyway, we can just take our chances with thread race conditions.
    my $depth  = 0;
    my $module = $tree;
    while ( $module ) {
	return
	    if ( $module->{'type'} eq "module" && $module->{'name'} eq "Error" );
	$module = &Galacticus::Build::SourceTree::Walk_Tree($module,\$depth);
    }
    # Walk the tree, looking for IO statements.
    my @ioStatements = ( 'open', 'close', 'read', 'write' );
    my $node         = $tree;
    $depth           = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "code" ) {
	    my $newContent    ;
	    my $inIO       = 0;
	    my $inCritical = 0;
	    my $ignoreLine = 0;
	    open(my $content,"<",\$node->{'content'});
	    while ( ! eof($content) ) {
		&Fortran::Utils::Get_Fortran_Line($content,my $rawLine, my $processedLine, my $bufferedComments);
		# Detect entry and exit to existing critical sections - we don't want to add locks inside such a section as that
		# would result in deadlocks.
		if ( $processedLine =~ m/^\s*!\$omp\s+critical\s*\(gfortranInternalIO(_)??\)\s*$/ ) {
		    my $isExplicit = ! defined($1);
		    $processedLine =~ s/\(gfortranInternalIO\)/(gfortranInternalIO_)/;
		    if ( $report && $isExplicit ) {
			$processedLine = $processedLine."    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),' acquired the ''gfortranInternalIO'' lock'\n"               ;
			$ignoreLine    = 1;
			&addUse($node);
		    }
		    $inCritical = 1;
		}
		if ( $processedLine =~ m/^\s*!\$omp\s+end\s+critical\s*\(gfortranInternalIO(_)??\)\s*$/ ) {
		    my $isExplicit = ! defined($1);
		    $processedLine =~ s/\(gfortranInternalIO\)/(gfortranInternalIO_)/;
		    if ( $report && $isExplicit ) {
			$processedLine =                "    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),' acquired the ''gfortranInternalIO'' lock'\n".$processedLine;
			$ignoreLine    = 1;
			&addUse($node);
		    }
		    $inCritical = 0;
		}
		# Convert "FoX_DOM_Access" critical sections to "gfortranInternalIO_". This still prevents multi-threaded access to
		# the FoX library, but also restricts access to gfortran internal IO within the FoX library.
		if ( $rawLine =~ m/^\s*!\$omp\s+(end\s+)??critical\s*\(FoX_DOM_Access\)\s*$/ ) {
		    $rawLine =~ s/FoX_DOM_Access/gfortranInternalIO_/;
		    if ( $report ) {
			&addUse($node);
			if ( defined($1) ) {
			    $rawLine =          "    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),' released the ''gfortranInternalIO'' lock'\n".$rawLine;
			} else {
			    $rawLine = $rawLine."    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),' acquired the ''gfortranInternalIO'' lock'\n"         ;
			}
		    }
		}
		# Detect IO statements and lock if not already in a critical section.
		my $isIO = 0;
		foreach my $ioStatement ( @ioStatements ) {
		    if ( $processedLine =~ m/^\s*(!\$)??\s*$ioStatement\s*\(/ ) {
			$isIO = 1;
			last;
		    }
		}
		$isIO = 1
		    if ( $processedLine =~ m/^\s*(!\$)??\s*call\s+flush\s*\(/i );
		if ( $isIO && ! $inIO && ! $inCritical && ! $ignoreLine ) {
		    $inIO        = 1;
		    $newContent .= "    !\$omp critical(gfortranInternalIO_)\n";
		    if ( $report ) {
			$newContent  .= "    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),' acquired the ''gfortranInternalIO'' lock'\n";
			&addUse($node);
		    }
		}
		if ( ! $isIO && $inIO && ! $inCritical && ! $ignoreLine ) {
		    $inIO        = 0;
		    if ( $report ) {
			$newContent  .= "    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),' released the ''gfortranInternalIO'' lock'\n";
			&addUse($node);
		    }
		    $newContent .= "    !\$omp end critical(gfortranInternalIO_)\n";
		}
		$newContent .= $rawLine;
	    }
	    close($content);
	    if ( $inIO && ! $inCritical && ! $ignoreLine ) {
		if ( $report ) {
		    $newContent  .= "    write (output_unit,*) '*** thread ',OMP_Get_Thread_Num(),' released the ''gfortranInternalIO'' lock'\n";
		    &addUse($node);
		}
		$newContent .= "    !\$omp end critical(gfortranInternalIO_)\n";		
	    }
	    $node->{'content'} = $newContent;
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

sub addUse {
    my $node      = shift();
    my $container = $node;
    while ( defined($container) ) {
	last
	    if (
		$container->{'type'} eq "function"
		||
		$container->{'type'} eq "subroutine"
		||
		$container->{'type'} eq "moduleProcedure"
		||
		$container->{'type'} eq "module"
		||
		$container->{'type'} eq "program"
	    );
	$container = $container->{'parent'};
    }
    # Add use of the OpenMP function.
    my $usesNode =
    {
	type      => "moduleUse",
	moduleUse =>
	{
	    "OMP_Lib"   =>
	    {
		openMP    => 1,
		intrinsic => 0,
		only      => {OMP_Get_Thread_Num => 1}
	    }
	}
    };
    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($container,$usesNode);
}

1;

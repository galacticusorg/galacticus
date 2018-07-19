# Contains a Perl module which implements debugging tools for MPI.

package Galacticus::Build::SourceTree::Process::DebugMPI;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use Fortran::Utils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'debugMPI'} = \&Process_DebugMPI;

sub Process_DebugMPI {
    # Get the tree.
    my $tree = shift();
    # Check if debugging is required.
    my $debug = 0;
    if ( exists($ENV{'GALACTICUS_FCFLAGS'}) ) {
	$debug = 1
	    if ( grep {$_ eq "-DDEBUGMPI"} split(" ",$ENV{'GALACTICUS_FCFLAGS'}) );
    }
    return
	unless ( $debug );    
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "code" ) {
	    my $newContent;
	    my $lineNumber = $node->{'line'};
	    open(my $content,"<",\$node->{'content'});
	    do {
		&Fortran::Utils::Get_Fortran_Line($content,my $rawLine, my $processedLine, my $bufferedComments);
		++$lineNumber;
		if ( $processedLine =~ m/mpiSelf%([a-zA-Z0-9_]+)/ ) {
		    my $method = $1;
		    $rawLine = 
			"write (0,*) 'MPI DEBUG ['//char(mpiSelf%rankLabel())//': mpiSelf call to method \"".$method."\"'//".
			&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber)."\n".
			$rawLine;
		}
		if ( $processedLine =~ m/^\s*call\s+mpiBarrier\s*\(\s*\)/ ) {
		    $rawLine = 
			"write (0,*) 'MPI DEBUG ['//char(mpiSelf%rankLabel())//']: mpiBarrier'//".
			&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber)."\n".
			$rawLine;
		}
		$newContent .= $rawLine;
	    } until ( eof($content) );
	    close($content);
	    $node->{'content'} = $newContent;
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

# Contains a Perl module which implements debugging tools for the HDF5 interface.

package Galacticus::Build::SourceTree::Process::DebugHDF5;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use List::ExtraUtils;
## AJB HACK use Galacticus::Build::SourceTree::Hooks;
## AJB HACK use Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'debugHDF5'} = \&Process_DebugHDF5;

sub Process_DebugHDF5 {
    # Get the tree.
    my $tree = shift();
    # Check if debugging is required.
    my $debug = 0;
    if ( exists($ENV{'GALACTICUS_FCFLAGS'}) ) {
	$debug = 1
	    if ( grep {$_ eq "-DDEBUGHDF5"} split(" ",$ENV{'GALACTICUS_FCFLAGS'}) );
    }
    return
	unless ( $debug );    
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "code" ) {
	    my $newContent;
	    open(my $content,"<",\$node->{'content'});
	    while ( my $line = <$content> ) {
		$line .= "call IO_HDF5_Start_Critical()\n"
		    if ( $line =~ m /^\s*\!\$omp\s+critical\s*\(HDF5_Access\)\s*$/ );
		$line  = "call IO_HDF5_End_Critical()\n".$line
		    if ( $line =~ m /^\s*\!\$omp\s+end\s+critical\s*\(HDF5_Access\)\s*$/ );
		$newContent .= $line;
	    }
	    close($content);
	    $node->{'content'} = $newContent;
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

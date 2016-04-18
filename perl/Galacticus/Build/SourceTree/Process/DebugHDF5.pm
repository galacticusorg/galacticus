# Contains a Perl module which implements debugging tools for the HDF5 interface.

package DebugHDF5;
use strict;
use warnings;
use utf8;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use Data::Dumper;
require List::ExtraUtils;
require Galacticus::Build::SourceTree::Hooks;
require Galacticus::Build::SourceTree;

# Insert hooks for our functions.
$Hooks::processHooks{'debugHDF5'} = \&Process_DebugHDF5;

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
	$node = &SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

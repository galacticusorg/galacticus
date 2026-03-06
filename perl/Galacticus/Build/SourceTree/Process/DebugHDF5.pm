# Contains a Perl module which implements debugging tools for the HDF5 interface.

package Galacticus::Build::SourceTree::Process::DebugHDF5;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;

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
	    my $useAdded = 0;
	    my $newContent;
	    open(my $content,"<",\$node->{'content'});
	    while ( my $line = <$content> ) {
		if ( $line =~ m /call\s+hdf5Access\s*\%\s*set\s*\(\s*\)/ ) {
		    $line =~ s/call\s+hdf5Access\s*\%\s*set\s*\(\s*\)/call IO_HDF5_Start_Locked\(\)/
		        if ( &addUse($node,$useAdded) );
		}
		if ( $line =~ m /call\s+hdf5Access\s*\%\s*unset\s*\(\s*\)/ ) {
		    $line  =~ s/call\s+hdf5Access\s*\%\s*unset\s*\(\s*\)/call IO_HDF5_End_Locked\(\)/
			if ( &addUse($node,$useAdded) );
		}
		$newContent .= $line;
	    }
	    close($content);
	    $node->{'content'} = $newContent;
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

sub addUse {
    # Add required module uses.
    my $node = shift();
    # Find the parent function/subroutine/program.
    while ( $node->{'type'} ne "function" && $node->{'type'} ne "subroutine" && $node->{'type'} ne "program" ) {
	$node = $node->{'parent'};
    }
    # Check if we are in a module.
    my $moduleNode = $node;
    my $skipUse    = 0;
    my $moduleName = "";
    while ( $moduleNode ) {
	if ( $moduleNode->{'type'} eq "module" ) {
	    $moduleName = $moduleNode->{'name'};
	    $skipUse    = $moduleNode->{'name'} eq "IO_HDF5" || $moduleNode->{'name'} eq "Error";
	    last;
	}
	$moduleNode = $moduleNode->{'parent'};
    }
    my $addCall = $moduleName ne "Error" && $moduleName ne "IO_HDF5";
    # No need to add a module use statement if this is the IO_HDF5 or Error module.
    unless ( $skipUse ) {
	&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses(
	     $node,
	     {
		 moduleUse =>
		 {
		     IO_HDF5 =>
		     {
			 openMP => 0,
			 only =>
			 {
			     IO_HDF5_Start_Locked => 1,
			     IO_HDF5_End_Locked   => 1
			 }
		     }
		 }
	     }
	    );
    }
    $_[0] = 1;
    return $addCall;
}

1;

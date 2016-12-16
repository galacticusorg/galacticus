# Contains a Perl module which converts HDF5 Fortran types into C-interoperable types to silence compiler warnings.

package Galacticus::Build::SourceTree::Process::HDF5FCInterop;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use List::ExtraUtils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'hdf5FCInterop'} = \&Process_HDF5FCInterop;

sub Process_HDF5FCInterop {
    # Get the tree.
    my $tree = shift();
    # Build the type map.
    die("Galacticus::Build::SourceTree::Process::HDF5FCInterop::Process_HDF5FCInterop(): type map file does not exist")
	unless ( -e $ENV{'BUILDPATH'}."/hdf5FCInterop.dat" );
    my %typeMap;
    open(my $mapFile,$ENV{'BUILDPATH'}."/hdf5FCInterop.dat");
    while ( my $line = <$mapFile> ) {
	if ( $line =~ m/(\S+)\s*=\s*(\S+)/ ) {
	    $typeMap{$1} = $2;
	}
    }
    close($mapFile);
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "declaration" ) {
	    my $typeChanged = 0;
	    my %typesAdded;
	    foreach ( @{$node->{'declarations'}} ) {
		if ( $_->{'intrinsic'} eq "integer" && defined($_->{'type'}) && $_->{'type'} =~ m/(kind=)?\s*(\S+)/i ) {
		    my $prefix = defined($1) ? $1 : "";
		    my $kind   = lc($2);
		    if ( exists($typeMap{$kind}) ) {
			$_->{'type'} = $prefix.$typeMap{$kind};
			$typeChanged = 1;
			$typesAdded{$typeMap{$kind}} = 1;
		    }
		}
	    }
	    if ( $typeChanged ) {
		# Rebuild declarations.
		&Galacticus::Build::SourceTree::Parse::Declarations::BuildDeclarations($node);
		# Check that the ISO_C_Binding module is used.
		my $sibling = $node->{'parent'}->{'firstChild'};
		while ( $sibling ) {
		    if ( $sibling->{'type'} eq "moduleUse" ) {
			unless ( grep {lc($_) eq "iso_c_binding"} keys(%{$sibling->{'moduleUse'}}) ) {
			    my $moduleUses;
			    $moduleUses->{'moduleUse'}->{'ISO_C_Binding'} =
			    {
				intrinsic => 1,
				only      => \%typesAdded
			    };
			    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$moduleUses);
			}
			last;
		    }
		    $sibling = $sibling->{'sibling'};
		}
	    }
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;

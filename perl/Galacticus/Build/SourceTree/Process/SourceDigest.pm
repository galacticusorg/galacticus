# Contains a Perl module which implements processing of source digest directives.

package Galacticus::Build::SourceTree::Process::SourceDigest;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use Digest::MD5 qw(md5_base64);
use Fortran::Utils;
use List::ExtraUtils;
use Galacticus::Path;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'sourceDigests'} = \&Process_SourceDigests;

# Database of previously computed digests.
our %digests;

sub Process_SourceDigests {
    # Get the tree.
    my $tree = shift();
    # Walk the tree, looking for source digest directives.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "sourceDigest" && ! $node->{'directive'}->{'processed'} ) {
	    # Mark the directive as processed.
	    $node->{'directive'}->{'processed'} =  1;
	    # Generate declaration for the digest variable.
	    my @digestDeclaration = 
		(
		 { 
		     intrinsic  => "character",
		     type       => "len=22",
		     attributes => [ "parameter" ],
		     variables  => [ $node->{'directive'}->{'name'}."=\"".&Find_Hash($tree->{'name'})."\"" ] 
		 }
		);
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@digestDeclaration);
	}
	# Step to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

sub Find_Hash {
    # Get names of files to process.
    my @fileNames = @_;
    # Initialize an MD5 hash.
    my $hasher = Digest::MD5->new();
    # Iterate over files.
    foreach my $fileName ( @fileNames ) {
	# Process all source files upon which this file depends.
	(my $dependencyFileName = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/.d/;
	if ( -e $dependencyFileName ) {
	    open(my $dependencyFile,$dependencyFileName);
	    while ( my $objectFileName = <$dependencyFile> ) {
		chomp($objectFileName);
		(my $sourceFileNamePrefix = $objectFileName) =~ s/$ENV{'BUILDPATH'}\/(.*)\.o$/source\/$1/;
		foreach my $suffix ( "F90", "c", "h", "Inc", "cpp" ) {
		    my $sourceFileName = $sourceFileNamePrefix.".".$suffix;
		    if ( -e $sourceFileName ) {
			unless ( exists($digests{$sourceFileName}) ) {
			    my $fileHasher = Digest::MD5->new();
			    if ( $suffix eq "F90" || $suffix eq "Inc" ) {
				# Parse the file ignoring whitespace and comments.
				$fileHasher->add(&Fortran::Utils::read_file($sourceFileName,state => "raw", followIncludes => 1, includeLocations => [ "../source", "../".$ENV{'BUILDPATH'} ], stripRegEx => qr/^\s*![^\#\@].*$/, stripLeading => 1, stripTrailing => 1));
				# Search for use on any files from the data directory by this source file.
				&Hash_Data_Files(
				    $fileHasher,
				    map 
				    {$_->{'submatches'}->[0]} 
				    &Fortran::Utils::Get_Matching_Lines($sourceFileName,qr/[\"\'](data\/[a-zA-Z0-9_\.\-\/]+\.(xml|hdf5))[\"\']/)
				    );
			    } else {
				# Parse the raw file.
				open(my $sourceFile,$sourceFileName);
				$fileHasher->addfile($sourceFile);
				close($sourceFile);
			    }
			    $digests{$sourceFileName} = $fileHasher->b64digest();
			}
			$hasher->add($digests{$sourceFileName});
		    }
		}
	    }
	    close($dependencyFile);
	}
    }
    return $hasher->b64digest();
}

sub Hash_Data_Files {
    # Run a supplied list of files through a supplied MD5 hash object.
    my $hasher   = shift;
    my @files = @_;
    foreach ( @files ) {
    	# Run each data file through the MD5 hash.
    	my $dataFileName = &galacticusPath()."/".$_;
    	if ( -e $dataFileName ) {
    	    open(my $dataHandle,$dataFileName);
    	    $hasher->addfile($dataHandle);
    	    close($dataHandle);
    	}
    }
}

1;

# Contains a Perl module which implements processing of source digest directives.

package Galacticus::Build::SourceTree::Process::SourceDigest;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Digest::MD5 qw(md5_base64);
use Fortran::Utils;
use List::ExtraUtils;
use List::Uniq ':all';

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'sourceDigests'} = \&Process_SourceDigests;

# Database of previously computed digests.
our %digests;
our %compositeDigests;
our %modificationTimes;

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
	# Check for a pre-existing composite digest.
	if ( exists($compositeDigests{$fileName}) ) {
	    # Use the composite digest.
	    $hasher->add($compositeDigests{$fileName});
	} else {
	    # Process all source files upon which this file depends.
	    my $compositeHasher = Digest::MD5->new();
	    (my $dependencyFileName = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/.d/;
	    if ( -e $dependencyFileName ) {
		(my $hashFileName = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/.md5c/;
		my $useStoredCompositeHash = -e $hashFileName;
		if ( $useStoredCompositeHash ) {
		    open(my $dependencyFile,$dependencyFileName);
		    while ( my $objectFileName = <$dependencyFile> ) {
			chomp($objectFileName);
			(my $sourceFileNamePrefix = $objectFileName) =~ s/$ENV{'BUILDPATH'}\/(.*)\.o$/source\/$1/;
			foreach my $suffix ( "F90", "c", "h", "Inc", "cpp" ) {
			    my $sourceFileName = $sourceFileNamePrefix.".".$suffix;
			    if ( -e $sourceFileName ) {
				unless ( exists($digests{$sourceFileName}) ) {
				    # Determine if we can used a stored value of the hash.
				    (my $md5FileName = $sourceFileName) =~ s/^source//;
				    $md5FileName = $ENV{'BUILDPATH'}.$md5FileName.".md5";
				    $useStoredCompositeHash = 0
					unless ( -e $md5FileName && &modificationTime($hashFileName) > &modificationTime($md5FileName) );
				}
			    }
			}
		    }
		}
		if ( $useStoredCompositeHash ) {
		    # Use the stored composite hash.
		    open(my $md5File,$hashFileName);
		    $compositeDigests{$fileName} = <$md5File>;
		    close($md5File);
		    $hasher->add($compositeDigests{$fileName});
		} else {
		    open(my $dependencyFile,$dependencyFileName);
		    while ( my $objectFileName = <$dependencyFile> ) {
			chomp($objectFileName);
			(my $sourceFileNamePrefix = $objectFileName) =~ s/$ENV{'BUILDPATH'}\/(.*)\.o$/source\/$1/;
			foreach my $suffix ( "F90", "c", "h", "Inc", "cpp" ) {
			    my $sourceFileName = $sourceFileNamePrefix.".".$suffix;
			    if ( -e $sourceFileName ) {
				unless ( exists($digests{$sourceFileName}) ) {
				    # Determine if we can used a stored value of the hash.
				    (my $md5FileName = $sourceFileName) =~ s/^source//;
				    $md5FileName = $ENV{'BUILDPATH'}.$md5FileName.".md5";
				    my $useStoredHash = 0;
				    if ( -e $md5FileName && &modificationTime($md5FileName) > &modificationTime($sourceFileName) ) {
					$useStoredHash = 1;
					if ( $suffix eq "F90" || $suffix eq "Inc" ) {
					    foreach my $dataFileName (
						map 
						{$_->{'submatches'}->[0]} 
						&Fortran::Utils::Get_Matching_Lines($sourceFileName,qr/[\"\'](data\/[a-zA-Z0-9_\.\-\/]+\.(xml|hdf5))[\"\']/)
						) {
						$useStoredHash = 0
						    unless ( -e $md5FileName && &modificationTime($md5FileName) > &modificationTime($dataFileName) );
					    }
					}
				    }
				    if ( $useStoredHash ) {
					# Use the stored hash.
					open(my $md5File,$md5FileName);
					$digests{$sourceFileName} = <$md5File>;
					close($md5File);
				    } else {
					# Stored hash is out of date or does not exist. Compute the hash now and store it.
					my $fileHasher = Digest::MD5->new();
					if ( $suffix eq "F90" || $suffix eq "Inc" ) {
					    # Parse the file ignoring whitespace and comments.
					    $fileHasher->add(&Fortran::Utils::read_file($sourceFileName,state => "raw", followIncludes => 1, includeLocations => [ "../source", "../".$ENV{'BUILDPATH'} ], stripRegEx => qr/^\s*![^\#\@].*$/, stripLeading => 1, stripTrailing => 1));
					    # Search for use on any files from the data directory by this source file.
					    &Hash_Data_Files(
						$fileHasher,
						map 
						{$_->{'submatches'}->[1]} 
						&Fortran::Utils::Get_Matching_Lines($sourceFileName,qr/(char\s*\()??\s*galacticusPath\s*\(\s*pathTypeDataStatic\s*\)\s*\)??\/\/\]\s*[\"\']([a-zA-Z0-9_\.\-\/]+\.(xml|hdf5))[\"\']/)
						);
					} else {
					    # Parse the raw file.
					    open(my $sourceFile,$sourceFileName);
					    $fileHasher->addfile($sourceFile);
					    close($sourceFile);
					}
					$digests{$sourceFileName} = $fileHasher->b64digest();
					open(my $md5File,">".$md5FileName);
					print $md5File $digests{$sourceFileName};
					close($md5File);
					&updateModificationTime($md5FileName);
				    }
				}
				if ( exists($digests{$sourceFileName}) ) {
				    $compositeHasher->add($digests{$sourceFileName});
				} else {
				    die("Galacticus::Build::SourceTree::Process::SourceDigest::Process_SourceDigests: failed to build digest for '".$sourceFileName."'");
				}
			    }
			}
		    }
		    close($dependencyFile);
		    $compositeDigests{$fileName} = $compositeHasher->b64digest();
		    $hasher->add($compositeDigests{$fileName});
		    open(my $md5File,">".$hashFileName);
		    print $md5File $compositeDigests{$fileName};
		    close($md5File);
		    &updateModificationTime($fileName);
		}
	    }
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
    	my $dataFileName = $ENV{'GALACTICUS_DATA_PATH'}."/".$_;
    	if ( -e $dataFileName ) {
    	    open(my $dataHandle,$dataFileName);
    	    $hasher->addfile($dataHandle);
    	    close($dataHandle);
    	}
    }
}

sub modificationTime {
    # Return the modification time of a file.
    my $fileName = shift();
    $modificationTimes{$fileName} = (stat($fileName))[9]
	unless ( exists($modificationTimes{$fileName}) );
    return $modificationTimes{$fileName};
}

sub updateModificationTime {
    # Update the modification time of a file.
    my $fileName = shift();
    $modificationTimes{$fileName} = (stat($fileName))[9];
}

1;

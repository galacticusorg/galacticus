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
use Fcntl ':flock';

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
		     type       => "C_Char",
		     attributes => [ "dimension(23)", "bind(C, name=\"".$node->{'directive'}->{'name'}."MD5\")" ],
		     variables  => [ $node->{'directive'}->{'name'} ] 
		 }
		);
	    &Galacticus::Build::SourceTree::Parse::Declarations::AddDeclarations($node->{'parent'},\@digestDeclaration);
	    my $bindingNode =
	    {
		type       => "moduleUse",
	        sibling    => undef()    ,
	        parent     => undef()    ,
		firstChild => undef()    ,
                moduleUse  => {ISO_C_Binding => {intrinsic => 1, only => {C_Char => 1}}}
	    };
            &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$bindingNode);
	}
	# Step to the next node in the tree.
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

sub Binding {
    # Return code for a binding to the C global variable containing the named hash.
    my $name = shift();
    # Generate the code.
    my $code = "   character(C_Char), dimension(23), bind(C, name=\"".$name."MD5\") :: ".$name."5\n";
    return $code;
}

sub Find_Hash {
    # Get names of files to process.
    my @fileNames = @{shift(@_)};
    my (%options) =         @_
	if ( scalar(@_) > 0 );
    # Determine if files should be locked.
    my $doLock = 1;
    if ( exists($ENV{'LOCKMD5'}) ) {
	$doLock = $ENV{'LOCKMD5'} eq "yes";
    }
    # Set default set of include files to exclude.
    @{$options{'includeFilesExcluded'}} = ()
	unless ( exists($options{'includeFilesExcluded'}) );
    # Initialize reporting.
    $options{'report'} = 0
	unless ( exists($options{'report'}) );
    print "=> Begin computing MD5 hash\n"
	if ( $options{'report'} );
    # Initialize an MD5 hash.
    my $hasher = Digest::MD5->new();
    # Iterate over files.
    foreach my $fileName ( @fileNames ) {
	print " => Process file: ".$fileName."\n"
	    if ( $options{'report'} );
	# Check for a pre-existing composite digest.
	if ( exists($compositeDigests{$fileName}) ) {
	    # Use the composite digest.
	    $hasher->add($compositeDigests{$fileName});
	    print "  => Use pre-existing composite hash: ".$compositeDigests{$fileName}."\n"
		if ( $options{'report'} );
	} else {
	    # Process all source files upon which this file depends, plus the file itself.
	    my $compositeHasher = Digest::MD5->new();
	    (my $hashFileName       = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/.md5c/;
	    (my $dependencyFileName = $ENV{'BUILDPATH'}."/".$fileName) =~ s/\.F90$/.d/;
	    if ( -e $dependencyFileName ) {
		print "  => Processing dependencies\n"
		    if ( $options{'report'} );
		open(my $md5Lock,">".$hashFileName.".lock");
		flock($md5Lock,LOCK_EX) or die "Could not lock '".$hashFileName.".lock' - $!"
		    if ( $doLock );
		my $useStoredCompositeHash = -e $hashFileName;
		if ( $useStoredCompositeHash ) {
		    open(my $dependencyFile,$dependencyFileName);
		    while ( my $objectFileName = <$dependencyFile> ) {
			chomp($objectFileName);
			(my $sourceFileNamePrefix = $objectFileName) =~ s/$ENV{'BUILDPATH'}\/(.*)\.o$/source\/$1/;
			print "   => Dependency file: ".$sourceFileNamePrefix."\n"
			    if ( $options{'report'} );
			foreach my $suffix ( "F90", "c", "h", "Inc", "cpp" ) {
			    my $sourceFileName = $sourceFileNamePrefix.".".$suffix;
			    if ( -e $sourceFileName ) {
				print "    => Dependency file: ".$sourceFileName."\n"
				    if ( $options{'report'} );
				(my $md5FileName = $sourceFileName) =~ s/^source//;
				$md5FileName = $ENV{'BUILDPATH'}.$md5FileName.".md5";
				if ( ! exists($digests{$sourceFileName}) ) {
				    # Determine if we can used a stored value of the hash.
				    $useStoredCompositeHash = 0
					unless ( -e $md5FileName && &modificationTime($hashFileName) > &modificationTime($md5FileName) && &modificationTime($md5FileName) > &modificationTime($sourceFileName) );
				    if ( $options{'report'} ) {
					if ( $useStoredCompositeHash ) {
					    print "     => Can use stored hash: "    .$sourceFileName."\n";
					} else {
					    print "     => Can not use stored hash: ".$sourceFileName."\n";
					}
				    }
				} else {
				    $useStoredCompositeHash = 0
					unless ( &modificationTime($hashFileName) > &modificationTime($md5FileName) );
				    if ( $options{'report'} ) {
					if ( $useStoredCompositeHash ) {
					    print "     => Digest pre-exists\n";
					} else {
					    print "     => Digest pre-exists but is outdated\n";
					}
				    }
				}
			    }
			}
		    }
		    close($dependencyFile);
		}
		if ( $useStoredCompositeHash ) {
		    # Use the stored composite hash.
		    open(my $md5File,$hashFileName);
		    my $compositeDigest = <$md5File>;
		    close($md5File);
		    $useStoredCompositeHash = defined($compositeDigest);
		    if ( $useStoredCompositeHash ) {
			$compositeDigests{$fileName} = $compositeDigest;
			$hasher->add($compositeDigests{$fileName});
			print "   => Reading stored composite hash:\t".$fileName."\t".$hashFileName."\t".$compositeDigests{$fileName}."\n"
			    if ( $options{'report'} );
		    }
		}
		if ( ! $useStoredCompositeHash ) {
		    print "   => Computing composite hash\n"
			if ( $options{'report'} );
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
				    open(my $md5Lock,">".$md5FileName.".lock");
				    flock($md5Lock,LOCK_EX) or die "Could not lock '".$md5FileName.".lock' - $!"
					if ( $doLock );
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
					my $digest = <$md5File>;
					close($md5File);
					$useStoredHash = defined($digest);
					if ( $useStoredHash ) {
					    $digests{$sourceFileName} = $digest;
					    print "   => Reading stored hash: ".$sourceFileName."\t".$digests{$sourceFileName}."\n"
						if ( $options{'report'} );
					}
				    }
				    if ( ! $useStoredHash ) {
					# Stored hash is out of date or does not exist. Compute the hash now and store it.
					print "   => Computing hash\n"
					    if ( $options{'report'} );
					my $fileHasher = Digest::MD5->new();
					if ( $suffix eq "F90" || $suffix eq "Inc" ) {
					    # Parse the file ignoring whitespace and comments.
					    $fileHasher->add(&Fortran::Utils::read_file($sourceFileName,state => "raw", followIncludes => 1, includeLocations => [ "../source", "../".$ENV{'BUILDPATH'} ], includeFilesExcluded => $options{'includeFilesExcluded'}, stripRegEx => qr/^\s*!(?!(!\[|\$)).*$/, stripLeading => 1, stripTrailing => 1, stripEmpty => 1));
					    # Search for use on any files from the data directory by this source file.
					    my @extraFiles = 
						map 
						{$_->{'submatches'}->[1]} 
					        &Fortran::Utils::Get_Matching_Lines($sourceFileName,qr/(char\s*\()??\s*galacticusPath\s*\(\s*pathTypeDataStatic\s*\)\s*\)??\/\/\]\s*[\"\']([a-zA-Z0-9_\.\-\/]+\.(xml|hdf5))[\"\']/);
					    &Hash_Data_Files(
						$fileHasher,
						@extraFiles
						);
					    print "    => Computed hash from: ".$sourceFileName." {".join(", ",@extraFiles)."}\n"
						if ( $options{'report'} );
					} else {
					    # Parse the raw file.
					    open(my $sourceFile,$sourceFileName);
					    $fileHasher->addfile($sourceFile);
					    close($sourceFile);
					    print "    => Computed hash from: ".$sourceFileName." {RAW}\n"
						if ( $options{'report'} );
					}
					$digests{$sourceFileName} = $fileHasher->b64digest();
					open(my $md5File,">".$md5FileName) or die $!;
					print $md5File $digests{$sourceFileName};
					close($md5File);
					&updateModificationTime($md5FileName);
					print "   => Stored hash: ".$sourceFileName."\t".$digests{$sourceFileName}."\n"
					    if ( $options{'report'} );
				    }
				    close($md5Lock);
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
		    print "   => Composite hash stored\t".$fileName."\t".$hashFileName."\t".$compositeDigests{$fileName}."\n"
			if ( $options{'report'} );
		}
		close($md5Lock);
	    }
	}
    }
    my $hash = $hasher->b64digest();
    print "=> MD5 hash: ".$hash."\n"
	if ( $options{'report'} );
    return $hash;
}

sub Hash_Data_Files {
    # Run a supplied list of files through a supplied MD5 hash object.
    my $hasher = shift();
    my @files  = @_     ;
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

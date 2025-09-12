#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Data::Dumper;
use Digest::MD5 qw(md5_base64);
use List::Uniq ':all';
use File::Slurp;
use Galacticus::Build::Directives;
use Galacticus::Build::SourceTree;
use Galacticus::Build::SourceTree::Process::FunctionClass::Utils;
use Storable;

# Build source digests for functionClass objects.
# Andrew Benson (06-February-2020)

# Get the name of the executable to compute source digests for.
die("Usage: sourceDigests.pl <sourceDirectory> <target> <useLocks>")
    unless ( scalar(@ARGV) == 3 );
my $sourceDirectoryName = $ARGV[0];
my $targetName          = $ARGV[1];
my $useLocks            = $ARGV[2];

# Include files to exclude from MD5 construction:
#  fftw3.f03 is part of the FFTW3 library, so not relevant to us;
#  galacticus.output.build.environment.inc - contains compiler information so can change depending on, e.g. MPI vs. non-MPI build
#  galacticus.output.version.revision.inc - contains revision number and build time
my @includeFilesExcluded = ( "fftw3.f03", "galacticus.output.build.environment.inc", "galacticus.output.version.revision.inc" );
# Specify a work directory.
my $workDirectoryName = $ENV{'BUILDPATH'}."/";
# Get an XML parser.
my $xml                     = new XML::Simple();
# Load the state storables file which lists all known functionClasses.
my $stateStorables          = -e $workDirectoryName."stateStorables.xml" ? $xml->XMLin($workDirectoryName."stateStorables.xml") : undef();
# Load the file of directive locations.
die("Error: directiveLocations.xml not found")
    unless ( -e $workDirectoryName."directiveLocations.xml" );
my $locations               = $xml->XMLin($workDirectoryName."directiveLocations.xml");
# Extract names of functionClasses.
my @functionClasses = map {$_ =~ s/Class$//; $_} sort(keys(%{$stateStorables->{'functionClasses'}}));
# Build a list of files containing functionClass implementations.
my @functionClassFileList = map {&List::ExtraUtils::as_array($locations->{$_}->{'file'})} @functionClasses;
my %functionClassFiles = map { $_ => 1 } @functionClassFileList;
# Build list of allowed names.
my @allowedNames = ( "functionClass", @functionClasses, @{$stateStorables->{'functionClassInstances'}} );
if ( exists($stateStorables->{'functionClassTypes'}->{'name'}) ) {
    push(@allowedNames,            $stateStorables->{'functionClassTypes'}->{'name'}   );
} else {
    push(@allowedNames,sort(keys(%{$stateStorables->{'functionClassTypes'}          })));
}
# Build a list of object file dependencies.
(my $dependencyFileName = $ENV{'BUILDPATH'}."/".$targetName) =~ s/\.(exe|o)$/\.d/;
my @objectFiles = map { $_ =~ /^$ENV{'BUILDPATH'}\/(.+\.o)$/ ? $1 : () } read_file($dependencyFileName, chomp => 1);
# Initialize structure to hold record of parameters from each source file.
(my $blobFileName = $targetName) =~ s/\.(exe|o)$/.md5.blob/;
my $digestsPerFile;
my $havePerFile = -e $ENV{'BUILDPATH'}."/".$blobFileName;
my $updateTime;
if ( $havePerFile ) {
    $digestsPerFile = retrieve($ENV{'BUILDPATH'}."/".$blobFileName);
    $updateTime     = -M       $ENV{'BUILDPATH'}."/".$blobFileName ;
}
# Hash of types which were updated.
my %updatedTypes;
# Open the source diretory, finding F90 and cpp files.
opendir(my $sourceDirectory,$sourceDirectoryName."/source");
my @fileNames = sort(readdir($sourceDirectory));
close($sourceDirectory);
foreach my $fileName ( @fileNames ) {
    # Skip junk files.
    next
	if ( $fileName =~ m/^\.\#/ );
    # Skip non-F90, non-cpp files
    next
	unless ( $fileName =~ m/\.(F90|cpp)$/ );
    # Find corresponding object file name.    
    (my $objectFileName = $fileName) =~ s/\.(F90|cpp)$/\.o/;
    # Skip non-dependency files.
    next
	unless ( grep {$_ eq $objectFileName} @objectFiles );
    # Construct full file name.
    my $fileToProcess = $sourceDirectoryName."/source/".$fileName;
    # Check if file is updated. If it is not, skip processing it. If it is, remove previous record of uses and rescan.
    (my $fileIdentifier = $fileToProcess) =~ s/\//_/g;
    $fileIdentifier =~ s/^\._??//;
    my $rescan = 1;
    if ( $havePerFile && exists($digestsPerFile->{$fileIdentifier}) ) {
	$rescan = 0
	    unless ( grep {-M $_ < $updateTime} &List::ExtraUtils::as_array($digestsPerFile->{$fileIdentifier}->{'files'}) );
    }
    if ( $rescan ) {
	delete($digestsPerFile->{$fileIdentifier})
    	    if ( $havePerFile && exists($digestsPerFile->{$fileIdentifier}) );
	push(@{$digestsPerFile->{$fileIdentifier}->{'files'}},$fileToProcess);
	# Find the dependency file.
	(my $dependencyFileName = $fileToProcess) =~ s/^.*\/([^\/]+)\.F90$/$ENV{'BUILDPATH'}\/$1\.d/;
	open(my $dependencyFile,$dependencyFileName);
	while ( my $dependentFileName = <$dependencyFile> ) {
	    chomp($dependentFileName);
	    (my $dependencyRoot = $dependentFileName) =~ s/^.*\/([^\/]+)\.o$/source\/$1\./;
	    foreach my $suffix ( "F90", "Inc", "cpp", "c", "h" ) {
		if ( -e $dependencyRoot.$suffix ) {
		    push(@{$digestsPerFile->{$fileIdentifier}->{'files'}},$dependencyRoot.$suffix);
		    last;
		}
	    }
	}
	close($dependencyFile);
	# Handle any directives in this file.
	## sourceDigest directive.
	if ( grep {$_ eq $fileToProcess} &List::ExtraUtils::as_array($locations->{'sourceDigest'}->{'file'}) ) {
	    my @sourceDigests = &Galacticus::Build::Directives::Extract_Directives($fileToProcess,'sourceDigest');
	    foreach my $sourceDigest ( @sourceDigests ) {
		my $hashName = $sourceDigest->{'name'};
		$digestsPerFile->{'types'}->{$hashName}->{'sourceMD5'} = &Galacticus::Build::SourceTree::Process::SourceDigest::Find_Hash([$fileName],includeFilesExcluded => \@includeFilesExcluded, useLocks => $useLocks);
		@{$digestsPerFile->{'types'}->{$hashName}->{'dependencies'}} = ();
		$updatedTypes{$hashName} = 1;
	    }
	}
	## functionClassType directive.
	if ( grep {$_ eq $fileToProcess} &List::ExtraUtils::as_array($locations->{'functionClassType'}->{'file'}) ) {
	    my @functionClassTypes = &Galacticus::Build::Directives::Extract_Directives($fileToProcess,'functionClassType');
	    foreach my $functionClassType ( @functionClassTypes ) {
		my $hashName = $functionClassType->{'name'};
		$digestsPerFile->{'types'}->{$hashName}->{'sourceMD5'} = &Galacticus::Build::SourceTree::Process::SourceDigest::Find_Hash([$fileName],includeFilesExcluded => \@includeFilesExcluded, useLocks => $useLocks);
		@{$digestsPerFile->{'types'}->{$hashName}->{'dependencies'}} = ( "functionClass" );
		$updatedTypes{$hashName} = 1;
	    }
	}
	## functionClass directive.
	if ( grep {$_ eq $fileToProcess} &List::ExtraUtils::as_array($locations->{'functionClass'}->{'file'}) ) {
	    my @functionClasses = &Galacticus::Build::Directives::Extract_Directives($fileToProcess,'functionClass');
	    foreach my $functionClass ( @functionClasses ) {
		my $hashName = $functionClass->{'name'}."Class";
		$digestsPerFile->{'types'}->{$hashName}->{'sourceMD5'} = &Galacticus::Build::SourceTree::Process::SourceDigest::Find_Hash([$fileName],includeFilesExcluded => \@includeFilesExcluded, useLocks => $useLocks);
		if ( exists($functionClass->{'extends'}) ) {
		    @{$digestsPerFile->{'types'}->{$hashName}->{'dependencies'}} = ( $functionClass->{'extends'} );
		} else {
		    @{$digestsPerFile->{'types'}->{$hashName}->{'dependencies'}} = ( "functionClass"             );
		}
		$updatedTypes{$hashName} = 1;
	    }
	}
	# Walk the source code tree.
	if ( exists($functionClassFiles{$fileToProcess}) ) {
	    my $tree  = &Galacticus::Build::SourceTree::ParseFile($fileToProcess);
	    my $depth = 0;
	    my $node  = $tree;
	    while ( $node ) {
		if ( grep {$node->{'type'} eq $_} @functionClasses ) {
		    my $hashName = $node->{'directive'}->{'name'};
		    # Find which class it extends.
		    my $classNode = $node;
		    (my $class, my @classDependencies) = &Galacticus::Build::SourceTree::Process::FunctionClass::Utils::Class_Dependencies($classNode,$node->{'type'});
		    # For self-referencing types, remove the self-reference here.
		    @classDependencies = map {$_ eq $class->{'type'} ? () : $_} @classDependencies;
		    # Store hash and dependencies.
		    $digestsPerFile->{'types'}->{$hashName}->{'sourceMD5'} = &Galacticus::Build::SourceTree::Process::SourceDigest::Find_Hash([$fileName],includeFilesExcluded => \@includeFilesExcluded, useLocks => $useLocks);
		    @{$digestsPerFile->{'types'}->{$hashName}->{'dependencies'}} = @classDependencies;
		    $updatedTypes{$hashName} = 1;
		}
		$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
	    }
	}
    }
}
# Manually add the base "functionClass" type.
my $functionClassFileName = "objects.function_class.F90";
if ( ! $updateTime || ( -M "source/".$functionClassFileName < $updateTime ) ) {
    $digestsPerFile->{'types'}->{'functionClass'}->{'sourceMD5'} = &Galacticus::Build::SourceTree::Process::SourceDigest::Find_Hash([$functionClassFileName],includeFilesExcluded => \@includeFilesExcluded, useLocks => $useLocks);
    @{$digestsPerFile->{'types'}->{'functionClass'}->{'dependencies'}} = ();
    delete($digestsPerFile->{'types'}->{'functionClass'}->{'compositeMD5'});
    $updatedTypes{'functionClass'} = 1;
}
# Recursively remove the composite hash for this type, and any dependent type.
my %dependenciesInverted;
foreach my $hashName ( keys(%{$digestsPerFile->{'types'}}) ) {
    foreach my $hashNameDependent ( @{$digestsPerFile->{'types'}->{$hashName}->{'dependencies'}} ) {
	push(@{$dependenciesInverted{$hashNameDependent}},$hashName);
    }
}
while ( scalar(keys(%updatedTypes)) > 0 ) {
    my ($hashNameReset) = %updatedTypes;
    delete($updatedTypes{$hashNameReset});
    foreach ( @{$dependenciesInverted{$hashNameReset}} ) {
 	$updatedTypes{$_} = 1;
    }
    delete($digestsPerFile->{'types'}->{$hashNameReset}->{'compositeMD5'})
	if ( exists($digestsPerFile->{'types'}->{$hashNameReset}->{'compositeMD5'}) );
}
# Construct composite hashes which include the hashes of all dependencies.
my $resolved = 0;
while ( ! $resolved ) {
    $resolved = 1;
    my $updated  = 0;
    foreach my $hashName ( sort(keys(%{$digestsPerFile->{'types'}})) ) {
	# If a composite MD5 is already computed for this type, skip it.
	next
	    if ( exists($digestsPerFile->{'types'}->{$hashName}) && exists($digestsPerFile->{'types'}->{$hashName}->{'compositeMD5'}) );
	# Check if all dependencies are resolved.
	if ( grep {exists($digestsPerFile->{'types'}->{$_}) && ! exists($digestsPerFile->{'types'}->{$_}->{'compositeMD5'})} @{$digestsPerFile->{'types'}->{$hashName}->{'dependencies'}} ) {
	    # Not all dependencies yet have a resolved compositeMD5 - record that we're not yet resolved and skip this type.
	    $resolved = 0;
	} else {
	    # Dependencies are resolved - compute our composite MD5;
	    $updated   = 1;
	    my $hasher = Digest::MD5->new();
	    $hasher->add($digestsPerFile->{'types'}->{$hashName}->{'sourceMD5'});
	    foreach my $dependency ( @{$digestsPerFile->{'types'}->{$hashName}->{'dependencies'}} ) {
		next
		    unless ( grep {$_ eq $dependency} @allowedNames );
		$hasher->add($digestsPerFile->{'types'}->{$dependency}->{'compositeMD5'});
	    }
	    $digestsPerFile->{'types'}->{$hashName}->{'compositeMD5'} = $hasher->b64digest();
	}
    }
    die("sourceDigest.pl: failed to resolve dependencies")
	if ( ! $resolved && ! $updated );
}
# Output the per file digest data.
store($digestsPerFile,$ENV{'BUILDPATH'}."/".$blobFileName);
# Output the results.
(my $outputFileName = $targetName) =~ s/\.(exe|o)$/.md5s.c/;
open(my $outputFile,">".$ENV{'BUILDPATH'}."/".$outputFileName);
foreach my $hashName ( sort(keys(%{$digestsPerFile->{'types'}})) ) {
    (my $digest = $digestsPerFile->{'types'}->{$hashName}->{'compositeMD5'}) =~ s/\//\@/g;
    print $outputFile "char ".$hashName."MD5[]=\"".$digest."\";\n"
}
close($outputFile);

exit;

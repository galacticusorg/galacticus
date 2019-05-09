#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use Fortran::Utils;
use File::Changes;
use Storable;

# Locate all allocatable arrays in the code base and determine their type and dimensionality.
# Andrew Benson (07-September-2016)

# Get command line arguments.
die "Usage: allocatableArrays.pl <sourceDirectory>"
    unless ( scalar(@ARGV) == 1 );
my $sourceDirectoryName = $ARGV[0];
# Intrinsics to be excluded.
my @excludedIntrinsics = ( "type", "class" );
# Types to be excluded.
my @excludedTypes = ( "hid", "hsize_t", "c_size_t", "c_char", "omp_lock_kind" );
# Initialize record of allocatable classes.
my %allocatableClasses;
# Initialize per-file allocatable classes.
my $allocatableClassesPerFile;
my $havePerFile = -e $ENV{'BUILDPATH'}."/allocatableArraysPerFile.blob";
my $updateTime;
if ( $havePerFile ) {
    $allocatableClassesPerFile = retrieve($ENV{'BUILDPATH'}."/allocatableArraysPerFile.blob");
    $updateTime                = -M       $ENV{'BUILDPATH'}."/allocatableArraysPerFile.blob" ;
}
# Open the source directorys and scan
opendir(my $sourceDirectory,$sourceDirectoryName."/source") 
    or die "Can't open the source directory: #!";
my @sourceFileNames = grep {$_ =~ m/\.f(90)??$/i} readdir($sourceDirectory);
closedir($sourceDirectory);
foreach my $sourceFileName ( @sourceFileNames ) {
    # Construct source file path.
    my $sourceFilePathName = $sourceDirectoryName."/source/".$sourceFileName;
    # Check if file is updated.
    if ( ! $havePerFile || -M $sourceFilePathName < $updateTime ) {
	delete($allocatableClassesPerFile->{$sourceFileName})
	    if ( $havePerFile && exists($allocatableClassesPerFile->{$sourceFileName}) );
	# Open and read the file.
	open(my $sourceFile,$sourceFilePathName)
	    or die "Can't open input file: ".$sourceFilePathName;
	until ( eof($sourceFile) ) {
	    # Get the next line of this file.
	    &Fortran::Utils::Get_Fortran_Line($sourceFile,my $rawLine,my $processedLine,my $bufferedComments);
	    # Test for variable definition.
	    next
		unless ( grep {$processedLine =~ $Fortran::Utils::intrinsicDeclarations{$_}->{'regEx'}} keys(%Fortran::Utils::intrinsicDeclarations) );
	    # Test for allocatable variable.
	    my $descriptor = &Fortran::Utils::Unformat_Variables($processedLine);
	    next
		unless (
		    exists($descriptor->{'attributes'})
		    &&
		    (grep {$_ eq "allocatable"             } @{$descriptor->{'attributes'}})
		    && 
		    (grep {$_ =~ m/^dimension\s*\([:,]+\)/ } @{$descriptor->{'attributes'}})
		);
	    # Skip excluded intrinsics.
	    next
		if ( exists($descriptor->{'intrinsic'}) && grep {lc($descriptor->{'intrinsic'}) eq $_               } @excludedIntrinsics );
	    # Skip excluded types.
	    next
		if ( exists($descriptor->{'type'     }) && grep {   $descriptor->{'type'     }   =~ m/^(kind=)??$_/i} @excludedTypes      );
	    # Determine the rank of the allocatable.
	    my $rank;
	    foreach ( @{$descriptor->{'attributes'}} ) {
		$rank = ($1 =~ tr/,//)+1
		    if ( $_ =~ m/^dimension\s*\(([:,]+)\)/ );
	    }
	    # Modify type.
	    if ( exists($descriptor->{'type'}) ) {
		# Remove any "kind" specifier, and make the type lower case.
		$descriptor->{'type'} =~ s/^kind=//;
		$descriptor->{'type'} = lc($descriptor->{'type'});
		# For character variables, remove any length.
		delete($descriptor->{'type'})
		    if ( $descriptor->{'type'} =~ m/^len=/ );
	    }
	    # Build an identifier for this allocatable class, and store the descriptor.
	    $descriptor->{'intrinsic'} = lc($descriptor->{'intrinsic'});
	    #  Store the rank in the descriptor, and remove any variable and attribute data.
	    $descriptor->{'rank'} = $rank;
	    delete($descriptor->{$_})
		foreach ( "variables", "attributes" );
	    my $identifier = join(".",$descriptor->{'intrinsic'},exists($descriptor->{'type'}) ? $descriptor->{'type'} : (),$rank);
	    $identifier =~ s/\s//g;
	    $allocatableClassesPerFile->{$sourceFileName}->{$identifier} = $descriptor
		unless ( $allocatableClassesPerFile->{$sourceFileName}->{$identifier} );
	}
    }
}
# Combine allocatable types across all souce file.
foreach my $sourceFileName ( keys(%{$allocatableClassesPerFile}) ) {
    foreach my $identifier ( keys(%{$allocatableClassesPerFile->{$sourceFileName}}) ) {
	$allocatableClasses{$identifier} = $allocatableClassesPerFile->{$sourceFileName}->{$identifier}
	    unless ( exists($allocatableClasses{$identifier}) );
    }
}
# Output the allocatables to XML, replacing a preexisting file if and only if the new file differs.
my $allocatables;
@{$allocatables->{'allocatable'}} = map {$allocatableClasses{$_}} sort(keys(%allocatableClasses));
store($allocatableClassesPerFile,$ENV{'BUILDPATH'}."/allocatableArraysPerFile.blob");
my $xmlOutput = new XML::Simple (NoAttr => 1, RootName=>"allocatables");
open(outHndl,">".$ENV{'BUILDPATH'}."/allocatableArrays.xml.tmp");
print outHndl $xmlOutput->XMLout($allocatables);
close(outHndl);
&File::Changes::Update($ENV{'BUILDPATH'}."/allocatableArrays.xml",$ENV{'BUILDPATH'}."/allocatableArrays.xml.tmp", proveUpdate => "yes");

exit;

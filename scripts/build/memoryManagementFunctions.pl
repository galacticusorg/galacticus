#!/usr/bin/env perl
use strict;
use warnings;
no warnings 'once';
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use List::Util qw(max);
use List::ExtraUtils;
use File::Changes;
use LaTeX::Encode;
use Text::Template 'fill_in_string';

# Build include files which provide memory allocation and deallocation functions for arrays of a variety of types and sizes. These
# are used by utility.memory_management.F90.
# Andrew Benson (24-April-2007)

# Check arguments.
die("Usage: memoryManagementFunctions.pl")
    unless ( scalar(@ARGV) == 0 );
# Get an XML parser.
my $xml = new XML::Simple();
# Load the data structure describing what types of allocatable array are needed.
my $allocatables = $xml->XMLin($ENV{'BUILDPATH'}."/allocatableArrays.xml");
# Find the maximum length of intrinsic type names.
my $intrinsicNameLengthMaximum = max(map {length($_->{'intrinsic'})} @{$allocatables->{'allocatable'}});
# Define data structure to hold code fragments during construction.
my $codeFragments;
# Define list of code sections.
my @codeSections = ( "preContain", "postContain" );
# Open files for pre- and post-"contains" code that will be included into the Memory_Management module.
my $includeFiles;
open($includeFiles->{$_} ,">".$ENV{'BUILDPATH'}."/utility.memory_management.".$_.".inc.tmp" )
    foreach ( @codeSections );
# Write some header information to these files.
print {$includeFiles->{'preContain' }} "!% Contains interface and type definitions for memory management functions along with storage space for pointers and sizes.\n";
print {$includeFiles->{'preContain' }} "!% This file was created automatically by {\\normalfont \\ttfamily memoryUseageFunctions.pl}\n\n";
print {$includeFiles->{'postContain'}} "!% Contains memory management functions.\n";
print {$includeFiles->{'postContain'}} "!% This file was created automatically by {\\normalfont \\ttfamily memoryUseageFunctions.pl}\n\n";
# Begin construction of interfaces to allocate and deallocate functions.
$codeFragments->{'preContain'}->{"allocInterfaceCode"  } = "interface allocateArray\n  !% Generic interface to routines which allocate arrays.\n";
$codeFragments->{'preContain'}->{"deallocInterfaceCode"} = "interface deallocateArray\n  !% Generic interface to routines which deallocate arrays.\n";
# Specify list of kinds acceptable as array shape descriptors.
my @arrayShapeDescriptorKinds = ( "", "kind_int8" );
# Iterate over all classes of allocatable variable and generate code for them.
foreach my $allocatable ( @{$allocatables->{'allocatable'}}  ) {
    # Construct intrinsic name, rank, and a type name.
    $code::intrinsicName  = $allocatable->{'intrinsic'}.(exists($allocatable->{'type'}) ? " (kind=".$allocatable->{'type'}.")" : "").($allocatable->{'intrinsic'} eq "character" ? "(len=*)" : "");
    $code::rank           = $allocatable->{'rank'     }                                                                             ;
    $code::typeName       = $allocatable->{'intrinsic'}.(exists($allocatable->{'type'}) ? "_"      .$allocatable->{'type'}     : "");
    # Create function label based on type and dimensionality.
    ($code::functionLabel  = $code::typeName."_".$code::rank."D") =~ s/\s/_/g;
    # Append code to the allocate and deallocate interfaces for this variable type.
    $codeFragments->{'preContain'}->{"deallocInterfaceCode"} .= "  module procedure deallocateArray_".$code::functionLabel                         ."\n";
    $codeFragments->{'preContain'}->{"allocInterfaceCode"  } .= "  module procedure allocateArray_"  .$code::functionLabel.($_ eq "" ? "" : "_".$_)."\n"
	foreach ( @arrayShapeDescriptorKinds );
    # Add a deallocate function for this allocatable type.
    $codeFragments->{'postContain'}->{"dealloc".$code::functionLabel} = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine deallocateArray_{$functionLabel}(thisArray,memoryType,file,line)
  !% Deallocate a {$rank}D \{\normalfont \ttfamily {&LaTeX::Encode::latex_encode($typeName)}\} array.
  use Display
  use ISO_Varying_String
  use String_Handling
  implicit none
  {$intrinsicName}         , intent(inout), allocatable  , dimension({join(",",split(//,":" x $rank))}) :: thisArray
  integer                  , intent(in   ), optional                                                    :: memoryType
  character(len=*         ), intent(in   ), optional                                                    :: file
  integer                  , intent(in   ), optional                                                    :: line
  type     (varying_string)                                                                             :: message
  integer                                                                                               :: memoryTypeActual
  integer  (kind_int8     )                                                                             :: accumulation
  if (present(memoryType)) then
     memoryTypeActual=memoryType
  else
     memoryTypeActual=memoryTypeMisc
  end if
  accumulation=sizeof(thisArray)+allocationOverhead
  !$omp atomic
  usedMemory%memoryType(memoryTypeActual)%usage=usedMemory%memoryType(memoryTypeActual)%usage-accumulation
  deallocate(thisArray)
  if (displayVerbosity() >= verbosityLevelDebug) then
     if (present(file).and.present(line)) then
      message='memory deallocate: '
      message=message//sizeof(thisArray)+allocationOverhead
      message=message//' ['//file//':'//line//']'
      call displayMessage(message)
    end if
  end if
  return
end subroutine deallocateArray_{$functionLabel}
CODE
    # Add block of code containing allocateArray unit.
    foreach my $arrayShapeDescriptorKind ( @arrayShapeDescriptorKinds ) {
	$code::suffix         = $arrayShapeDescriptorKind eq "" ? "" : "_"     .$arrayShapeDescriptorKind    ;
	$code::typeDefinition = $arrayShapeDescriptorKind eq "" ? "" : "(kind=".$arrayShapeDescriptorKind.")";
	$codeFragments->{'postContain'}->{"alloc".$code::functionLabel.$code::suffix} = fill_in_string(<<'CODE', PACKAGE => 'code');
subroutine allocateArray_{$functionLabel.$suffix}(thisArray,dimensions,lowerBounds,memoryType,file,line)
  !% Allocate a {$rank}D \{\normalfont \ttfamily {&LaTeX::Encode::latex_encode($typeName)}\} array.
  use Display
  use ISO_Varying_String
  use String_Handling
  implicit none
  {$intrinsicName}                       , intent(inout), allocatable  , dimension({join(",",split(//,":" x $rank))}) :: thisArray
  integer                                , intent(in   ), optional                                                    :: memoryType
  integer{$typeDefinition}               , intent(in   )                                                              :: dimensions ({$rank})
  integer{$typeDefinition}               , intent(in   ), optional                                                    :: lowerBounds({$rank})
  character               (len=*        ), intent(in   ), optional                                                    :: file
  integer                                , intent(in   ), optional                                                    :: line
  type                    (varying_string)                                                                            :: message
  integer                                                                                                             :: memoryTypeActual
  integer                 (kind_int8     )                                                                            :: accumulation
  if (allocated(thisArray)) call deallocateArray_{$functionLabel}(thisArray,memoryType)
  if (present(lowerBounds)) then
     allocate(thisArray({join(",",map {"lowerBounds(".$_."):lowerBounds(".$_.")+dimensions(".$_.")-1"} 1..$rank)}))
  else
     allocate(thisArray({join(",",map {"dimensions(".$_.")"} 1..$rank)}))
  end if
  if (present(memoryType)) then
     memoryTypeActual=memoryType
  else
     memoryTypeActual=memoryTypeMisc
  end if
  accumulation=sizeof(thisArray)+allocationOverhead
  !$omp atomic
  usedMemory%memoryType(memoryTypeActual)%usage=usedMemory%memoryType(memoryTypeActual)%usage+accumulation
  if (displayVerbosity() >= verbosityLevelDebug) then
     if (present(file).and.present(line)) then
      message='memory allocate: '
      message=message//sizeof(thisArray)+allocationOverhead
      message=message//' ['//file//':'//line//']'
      call displayMessage(message)
    end if
  end if
  call Memory_Usage_Report
  return
end subroutine allocateArray_{$functionLabel.$suffix}
CODE
    }
}
# Add closing statements to the allocate and deallocate interface blocks.
$codeFragments->{'preContain'}->{'allocInterfaceCode'  } .= "end interface allocateArray\n\n";
$codeFragments->{'preContain'}->{'deallocInterfaceCode'} .= "end interface deallocateArray\n\n";
# Output code blocks to files.
foreach my $section ( @codeSections ) {
    print {$includeFiles->{$section}} $_."\n"
	foreach ( &List::ExtraUtils::hashList($codeFragments->{$section}) )
}
# Close the output files.
close($includeFiles->{$_} )
    foreach ( @codeSections );
# Truncate lines and update old files if necessary.
&File::Changes::Update($ENV{'BUILDPATH'}."/utility.memory_management.".$_.".inc",$ENV{'BUILDPATH'}."/utility.memory_management.".$_.".inc.tmp" )
    foreach ( @codeSections );
exit 0;

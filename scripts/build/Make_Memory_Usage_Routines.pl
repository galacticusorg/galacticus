#!/usr/bin/env perl
use lib './perl';
use Fortran::Utils;
use File::Changes;
use Switch;
use XML::Simple;

# Script which builds the include files used by utility.memory_management.F90
# Andrew Benson (24-Apr-2007)
#
# Essentially we need several near identical routines for each type (dimension and variable type) array that must be handled.
# This script automates the process of writing these routines and ensures that they all conform to the same methodology.

# Create an XML object.
$xml = new XML::Simple;

# Load the data structure describing what types of allocatable array are needed.
$allocatables = $xml->XMLin("./work/build/Allocatable_Arrays.xml");

# Find the longest name length.
$long_name=0;
foreach $allocatable ( @{$allocatables->{'allocatable'}}  ) {
    $name = $allocatable->{'type'};
    if (length($name) > $long_name) {$long_name = length($name)};
}

# Open files for pre and post "contains" code that will be included into the Memory_Management module.
open(preContainHandle ,">./work/build/utility.memory_management.precontain.inc.tmp" );
open(postContainHandle,">./work/build/utility.memory_management.postcontain.inc.tmp");

# Write some header information to these files.
print preContainHandle  "!% Contains interface and type definitions for memory management routines along with storage space for pointers and sizes.\n";
print preContainHandle  "!% This file was created automatically by {\\tt Make\\_Memory\\_Usage\\_Routines.pl}\n\n";
print postContainHandle "!% Contains memory management subroutines.\n";
print postContainHandle "!% This file was created automatically by {\\tt Make\\_Memory\\_Usage\\_Routines.pl}\n\n";

# Open interface to Alloc_Array routines.
$preBlocks{"allocInterfaceCode"} = "interface Alloc_Array\n  !% Generic interface to routines which allocate arrays.\n";

# Open interface to Dealloc_Array routines.
$preBlocks{"deallocInterfaceCode"} = "interface Dealloc_Array\n  !% Generic interface to routines which deallocate arrays.\n";

# Loop over all classes of allocatable variable and generate code for them.
foreach $allocatable ( @{$allocatables->{'allocatable'}}  ) {
    $typeName  = $allocatable->{'type'};      # Type of variable.
    $dimension = $allocatable->{'dimension'}; # Dimensionality.
    unless ( UNIVERSAL::isa($allocatable->{'kind'}, "HASH") ) {
	$kind    = $allocatable->{'kind'};      # Get kind number if present.
    } else {
	$kind = "";
    }
    $type = $typeName;
    $type =~ s/^(\w)/uc($1)/e;           # Capitalize first letter of word.
    $type =~ s/([\s_])(\w)/$1.uc($2)/ge; # Capitalize first letter of each subsequent word.
    $type =~ s/\s/_/g;                   # Convert spaces to underscores.
    $typeLowerCase = lc($type);          # All lower case version of type name.
    switch ( $typeName ) {
	case ( "character"        ) {
	    $typeSize = "len(thisArray)";
	    $typeName .= "(len=*)";
	}
    }

    # Create "kind" label if necessary.
    if ( $kind ne "" ) {
	$typeName .= " (kind=".$kind.")";
	$type     .= "_".$kind;
    }

    # Create version of type name with underscores escaped for LaTeX.
    $variableTypeLatex = $typeLowerCase;
    $variableTypeLatex =~ s/_/\\_/g;

    # Create label based on type and dimensionality.
    $typeLabel = $type."_".$dimension."D";
		
    # Append code to the various interfaces for this variable type.
    $preBlocks{"allocInterfaceCode"}    .= "  module procedure Alloc_Array_"  .$typeLabel.       "\n";
    $preBlocks{"deallocInterfaceCode"}  .= "  module procedure Dealloc_Array_".$typeLabel.       "\n";

    # Add block of code containing Dealloc_Array unit.
    $postBlocks{"deallocCode"} .= "subroutine Dealloc_Array_".$typeLabel."(thisArray,memoryType)\n";
    $postBlocks{"deallocCode"} .= "  !% Deallocate a ".$dimension."D ".$variableTypeLatex." array.\n";
    $postBlocks{"deallocCode"} .= "  use Galacticus_Error\n";
    $postBlocks{"deallocCode"} .= "  implicit none\n";
    $c1_width = 8;
    if ( length($typeName) > $c1_width ) {$c1_width = length($typeName)};
    $pad = " " x ($c1_width-length($typeName));
    $postBlocks{"deallocCode"} .= "  ".$typeName.",".$pad." allocatable    :: thisArray(".join(",",split(//,":" x $dimension)).")\n";
    $pad = " " x ($c1_width-length("integer"));
    $postBlocks{"deallocCode"} .= "  integer,".$pad." intent(in), optional :: memoryType\n";
    $postBlocks{"deallocCode"} .= "  integer ".$pad."                      :: memoryTypeActual\n\n";
    if ( $typeSize =~ /^sizeof\((.*)\)/ ) {
	$pad = " " x ($c1_width-length("$typeName"));
	$postBlocks{"deallocCode"} .= "  ".$typeName.$pad."                :: ".$1."\n";
    }
    $postBlocks{"deallocCode"} .= "  if (present(memoryType)) then\n";
    $postBlocks{"deallocCode"} .= "     memoryTypeActual=memoryType\n";
    $postBlocks{"deallocCode"} .= "  else\n";
    $postBlocks{"deallocCode"} .= "     memoryTypeActual=memoryTypeMisc\n";
    $postBlocks{"deallocCode"} .= "  end if\n";
    $postBlocks{"deallocCode"} .= "  !\$omp critical(Memory_Management_Usage)\n";
    $postBlocks{"deallocCode"} .= "  usedMemory%memoryType(memoryTypeActual)%usage=usedMemory%memoryType(memoryTypeActual)%usage-sizeof(thisArray)-allocationOverhead\n";
    $postBlocks{"deallocCode"} .= "  !\$omp end critical(Memory_Management_Usage)\n";
    $postBlocks{"deallocCode"} .= "  deallocate(thisArray)\n";
    $postBlocks{"deallocCode"} .= "  return\n";
    $postBlocks{"deallocCode"} .= "end subroutine Dealloc_Array_".$typeLabel."\n\n";

    # Add block of code containing Alloc_Array unit.
    $postBlocks{"allocCode"} .= "subroutine Alloc_Array_".$typeLabel."(thisArray,dimensions,lowerBounds,memoryType)\n";
    $postBlocks{"allocCode"} .= "  !% Allocate a ".$dimension."D ".$variableTypeLatex." array.\n";
    $postBlocks{"allocCode"} .= "  implicit none\n";
    $c1_width = 16;
    if ( length($typeName) > $c1_width ) {$c1_width = length($typeName)};
    $pad = " " x ($c1_width-length($typeName));
    $postBlocks{"allocCode"} .= "  ".$typeName.",".$pad." allocatable             :: thisArray(".join(",",split(//,":" x $dimension)).")\n";
    $pad = " " x ($c1_width-length("integer"));
    $postBlocks{"allocCode"} .= "  integer,".$pad." intent(in),    optional :: memoryType\n";
    $postBlocks{"allocCode"} .= "  integer,".$pad." intent(in)              :: dimensions (".$dimension.")\n";
    $postBlocks{"allocCode"} .= "  integer,".$pad." intent(in),    optional :: lowerBounds(".$dimension.")\n";
    $pad = " " x ($c1_width-length("character(len=*)"));
    $pad = " " x ($c1_width-length("integer")+1);
    $postBlocks{"allocCode"} .= "  integer ".$pad."                        :: memoryTypeActual\n\n";
    $postBlocks{"allocCode"} .= "  if (allocated(thisArray)) call Dealloc_Array_".$typeLabel."(thisArray,memoryType)\n";
    $postBlocks{"allocCode"} .= "  if (present(lowerBounds)) then\n";
    $postBlocks{"allocCode"} .= "     allocate(thisArray(".join(",",map {"lowerBounds(".$_."):lowerBounds(".$_.")+dimensions(".$_.")-1"} 1..$dimension)."))\n";
    $postBlocks{"allocCode"} .= "  else\n";
    $postBlocks{"allocCode"} .= "     allocate(thisArray(".join(",",map {"dimensions(".$_.")"} 1..$dimension)."))\n";
    $postBlocks{"allocCode"} .= "  end if\n";
    $postBlocks{"allocCode"} .= "  if (present(memoryType)) then\n";
    $postBlocks{"allocCode"} .= "     memoryTypeActual=memoryType\n";
    $postBlocks{"allocCode"} .= "  else\n";
    $postBlocks{"allocCode"} .= "     memoryTypeActual=memoryTypeMisc\n";
    $postBlocks{"allocCode"} .= "  end if\n";
    $postBlocks{"allocCode"} .= "\n";
    $postBlocks{"allocCode"} .= "  !\$omp critical(Memory_Management_Usage)\n";
    $postBlocks{"allocCode"} .= "  usedMemory%memoryType(memoryTypeActual)%usage=usedMemory%memoryType(memoryTypeActual)%usage+sizeof(thisArray)+allocationOverhead\n";
    $postBlocks{"allocCode"} .= "  !\$omp end critical(Memory_Management_Usage)\n";
    $postBlocks{"allocCode"} .= "  call Memory_Usage_Report\n";
    $postBlocks{"allocCode"} .= "  return\n";
    $postBlocks{"allocCode"} .= "end subroutine Alloc_Array_".$typeLabel."\n\n";
}

# Add closing statements to the various interface blocks.
$preBlocks{"allocInterfaceCode"}    .= "end interface\n\n";
$preBlocks{"deallocInterfaceCode"}  .= "end interface\n\n";

# Output code blocks to files.
foreach $preBlockKey ( sort keys %preBlocks ) {
    print preContainHandle $preBlocks{$preBlockKey};
}
foreach $postBlockKey ( sort keys %postBlocks ) {
    print postContainHandle $postBlocks{$postBlockKey};
}

# Close the output files.
close(preContainHandle);
close(postContainHandle);

# Truncate lines and update old files..
foreach $file ( "utility.memory_management.precontain.inc", "utility.memory_management.postcontain.inc", "utility.memory_management.use.inc" ) {
    &Fortran_Utils::Truncate_Fortran_Lines("./work/build/".$file.".tmp");
    &File_Changes::Update("./work/build/".$file ,"./work/build/".$file.".tmp" );
}

exit;

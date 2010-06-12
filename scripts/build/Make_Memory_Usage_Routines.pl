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

# Load the data structure describing where various derived type variables are defined.
$dependencies = $xml->XMLin("./work/build/Type_Definitions.xml");

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

# Create block of code containing Free_Memory routine.
$postBlocks{"freeMemoryCode"} = "subroutine Free_Memory\n  !% Deallocates all memory that was previously registered. Should be called just prior to program termination.\n  implicit none\n  return\nend subroutine Free_Memory\n";

# Open interface to the Register_Memory routines.
$preBlocks{"registerInterfaceCode"} = "interface Register_Memory\n  !% Generic interface to routines which register memory.\n";

# Open interface to Alloc_Array routines.
$preBlocks{"allocInterfaceCode"} = "interface Alloc_Array\n  !% Generic interface to routines which allocate arrays and register it.\n";

# Open interface to Dealloc_Array routines.
$preBlocks{"deallocInterfaceCode"} = "interface Dealloc_Array\n  !% Generic interface to routines which deallocate arrays and unregister them.\n";

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
    # Determine size of variable.
    switch ( $typeName ) {
	case ( "real"             ) {$typeSize =  4}
	case ( "integer"          ) {$typeSize =  4}
	case ( "double precision" ) {$typeSize =  8}
	case ( "logical"          ) {$typeSize =  4}
	case ( "complex"          ) {$typeSize =  8}
	case ( "double complex"   ) {$typeSize = 16}
	case ( "character"        ) {
	    $typeSize = "len(arr)";
	    $typeName .= "(len=*)";
	}
	case ( "varying_string"  ) {
	    $typeSize = 0;
	    $typeName = "type(varying_string)";
	    $use_modules{'iso_varying_string'} += 1;
	}
	else {
	    $typeSize="sizeof(Dummy_".$type.")";
	    if ( ! exists $use_modules{$typeName} ) {$use_modules{$dependencies->{$typeName}} += 1};
	    $typeName = "type(".$typeName.")";
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
    $preBlocks{"registerInterfaceCode"} .= "  module procedure Register_"     .$typeLabel."_Memory\n";
    $preBlocks{"allocInterfaceCode"}    .= "  module procedure Alloc_Array_"  .$typeLabel.       "\n";
    $preBlocks{"allocInterfaceCode"}    .= "  module procedure Alloc_Array_"  .$typeLabel."_LBound\n";
    $preBlocks{"deallocInterfaceCode"}  .= "  module procedure Dealloc_Array_".$typeLabel.       "\n";

    # Add a block of code containing register memory routine.
    $postBlocks{"regsiterUnitCode"} .= "subroutine Register_".$typeLabel."_Memory(arr,mem_type,isize)\n";
    $postBlocks{"regsiterUnitCode"} .= "  !% Tracks memory usage for ".$variableTypeLatex." ".$dimension."D array memory.\n";
    $postBlocks{"regsiterUnitCode"} .= "  use Galacticus_Error\n";
    $postBlocks{"regsiterUnitCode"} .= "  implicit none\n";
    $c1_width = 8;
    if ( length($typeName) > $c1_width ) {$c1_width = length($typeName)};
    if ( length($type)+23 > $c1_width ) {$c1_width = length($type)+23};
    $pad = " " x ($c1_width-length($typeName));
    $postBlocks{"regsiterUnitCode"} .= "  ".$typeName.",".$pad." allocatable   :: arr(".join(",",split(//,":" x $dimension)).")\n";
    $pad = " " x ($c1_width-length("integer"));
    $postBlocks{"regsiterUnitCode"} .= "  integer,".$pad." intent(in)    :: isize\n";
    $postBlocks{"regsiterUnitCode"} .= "  integer,".$pad." intent(in)    :: mem_type\n";
    if ( $typeSize =~ /^sizeof\((.*)\)/ ) {
	$pad = " " x ($c1_width-length("$typeName"));
	$postBlocks{"regsiterUnitCode"} .= "  ".$typeName.$pad."                :: ".$1."\n";
    }
    $postBlocks{"regsiterUnitCode"} .= "\n";
    $postBlocks{"regsiterUnitCode"} .= "  !\$ if (omp_in_parallel()) then\n";
    $postBlocks{"regsiterUnitCode"} .= "  !\$omp critical (MemAdd)\n";
    $postBlocks{"regsiterUnitCode"} .= "  !\$ Memory_Used%mem_type(mem_type)%usage=Memory_Used%mem_type(mem_type)%usage+dble(".$typeSize."*isize)\n";
    $postBlocks{"regsiterUnitCode"} .= "  !\$omp end critical (MemAdd)\n";
    $postBlocks{"regsiterUnitCode"} .= "  !\$ else\n";
    $postBlocks{"regsiterUnitCode"} .= "  Memory_Used%mem_type(mem_type)%usage=Memory_Used%mem_type(mem_type)%usage+dble(".$typeSize."*isize)\n";
    $postBlocks{"regsiterUnitCode"} .= "  !\$ end if\n";
    $postBlocks{"regsiterUnitCode"} .= "  call Report_Memory_Usage\n";
    $postBlocks{"regsiterUnitCode"} .= "  return\n";
    $postBlocks{"regsiterUnitCode"} .= "end subroutine Register_".$typeLabel."_Memory\n\n";

    # Add block of code containing Dealloc_Array unit.
    $postBlocks{"deallocCode"} .= "subroutine Dealloc_Array_".$typeLabel."(arr,iimem_type)\n";
    $postBlocks{"deallocCode"} .= "  !% Deallocate a ".$dimension."D ".$variableTypeLatex." array.\n";
    $postBlocks{"deallocCode"} .= "  use Galacticus_Error\n";
    $postBlocks{"deallocCode"} .= "  implicit none\n";
    $c1_width = 8;
    if ( length($typeName) > $c1_width ) {$c1_width = length($typeName)};
    $pad = " " x ($c1_width-length($typeName));
    $postBlocks{"deallocCode"} .= "  ".$typeName.",".$pad." allocatable    :: arr(".join(",",split(//,":" x $dimension)).")\n";
    $pad = " " x ($c1_width-length("integer"));
    $postBlocks{"deallocCode"} .= "  integer,".$pad." intent(in), optional :: iimem_type\n";
    $postBlocks{"deallocCode"} .= "  integer ".$pad."                      :: imem_type\n\n";
    if ( $typeSize =~ /^sizeof\((.*)\)/ ) {
	$pad = " " x ($c1_width-length("$typeName"));
	$postBlocks{"deallocCode"} .= "  ".$typeName.$pad."                :: ".$1."\n";
    }
    $postBlocks{"deallocCode"} .= "  if (present(iimem_type)) then\n";
    $postBlocks{"deallocCode"} .= "     imem_type=iimem_type\n";
    $postBlocks{"deallocCode"} .= "  else\n";
    $postBlocks{"deallocCode"} .= "     imem_type=MEM_TYPE_MISC\n";
    $postBlocks{"deallocCode"} .= "  end if\n";
    $postBlocks{"deallocCode"} .= "  !\$ if (omp_in_parallel()) then\n";
    $postBlocks{"deallocCode"} .= "  !\$omp critical (MemAdd)\n";
    $postBlocks{"deallocCode"} .= "  !\$ Memory_Used%mem_type(imem_type)%usage=Memory_Used%mem_type(imem_type)%usage-dble(".$typeSize;
    for ($i=1;$i<=$dimension;++$i) {
	$postBlocks{"deallocCode"} .= "*size(arr,dim=".$i.")";
    }
    $postBlocks{"deallocCode"} .= ")\n";
    $postBlocks{"deallocCode"} .= "  !\$omp end critical (MemAdd)\n";
    $postBlocks{"deallocCode"} .= "  !\$ else\n";
    $postBlocks{"deallocCode"} .= "  Memory_Used%mem_type(imem_type)%usage=Memory_Used%mem_type(imem_type)%usage-dble(".$typeSize;
    for ($i=1;$i<=$dimension;++$i) {
	$postBlocks{"deallocCode"} .= "*size(arr,dim=".$i.")";
    }
    $postBlocks{"deallocCode"} .= ")\n";
    $postBlocks{"deallocCode"} .= "  !\$ end if\n";
    $postBlocks{"deallocCode"} .= "  deallocate(arr)\n";
    $postBlocks{"deallocCode"} .= "  return\n";
    $postBlocks{"deallocCode"} .= "end subroutine Dealloc_Array_".$typeLabel."\n\n";

    # Add block of code containing Alloc_Array unit.
    $postBlocks{"allocCode"} .= "subroutine Alloc_Array_".$typeLabel."(arr,".join(",",map {"size".$_} 1..$dimension).",name,mem_type)\n";
    $postBlocks{"allocCode"} .= "  !% Allocate a ".$dimension."D ".$variableTypeLatex." array and register it.\n";
    $postBlocks{"allocCode"} .= "  use Galacticus_Error\n";
    $postBlocks{"allocCode"} .= "  implicit none\n";
    $c1_width = 16;
    if ( length($typeName) > $c1_width ) {$c1_width = length($typeName)};
    $pad = " " x ($c1_width-length($typeName));
    $postBlocks{"allocCode"} .= "  ".$typeName.",".$pad." allocatable             :: arr(".join(",",split(//,":" x $dimension)).")\n";
    $pad = " " x ($c1_width-length("integer"));
    $postBlocks{"allocCode"} .= "  integer,".$pad." intent(in),    optional :: mem_type\n";
    $postBlocks{"allocCode"} .= "  integer,".$pad." intent(in)              :: ".join(",",map {"size".$_} 1..$dimension)."\n";
    $pad = " " x ($c1_width-length("character(len=*)"));
    $postBlocks{"allocCode"} .= "  character(len=*),".$pad." intent(in)              :: name\n";
    $pad = " " x ($c1_width-length("integer")+1);
    $postBlocks{"allocCode"} .= "  integer,".$pad." save                   :: alloc_err,imem_type,tsize,".join(",",map {"llim".$_} 1..$dimension)."\n";
    $postBlocks{"allocCode"} .= "  !\$omp threadprivate(alloc_err,imem_type,tsize,".join(",",map {"llim".$_} 1..$dimension).")\n\n";
    $postBlocks{"allocCode"} .= join("\n",map {"llim".$_."=1"} 1..$dimension)."\n";
    $postBlocks{"allocCode"} .= "  if (allocated(arr)) call Dealloc_Array_".$typeLabel."(arr,imem_type)\n";
    $postBlocks{"allocCode"} .= "  allocate(arr(".join(",",map {"llim".$_.":size".$_} 1..$dimension)."),stat=alloc_err)\n";
    $postBlocks{"allocCode"} .= "  if (alloc_err.ne.0) then\n";
    $postBlocks{"allocCode"} .= "     write (0,*) 'Alloc_Array_".$typeLabel."(): FATAL - unable to allocate memory [',trim(name),']'\n";
    $postBlocks{"allocCode"} .= "     call Galacticus_Error_Report\n";
    $postBlocks{"allocCode"} .= "  end if\n";
    $postBlocks{"allocCode"} .= "  if (present(mem_type)) then\n";
    $postBlocks{"allocCode"} .= "     imem_type=mem_type\n";
    $postBlocks{"allocCode"} .= "  else\n";
    $postBlocks{"allocCode"} .= "     imem_type=MEM_TYPE_MISC\n";
    $postBlocks{"allocCode"} .= "  end if\n";
    $postBlocks{"allocCode"} .= "  tsize=".join("*",map {"(size".$_."-llim".$_."+1)"} 1..$dimension)."\n";
    $postBlocks{"allocCode"} .= "  call Register_Memory(arr,imem_type,tsize)\n";
    $postBlocks{"allocCode"} .= "  return\n";
    $postBlocks{"allocCode"} .= "end subroutine Alloc_Array_".$typeLabel."\n\n";
    
    $postBlocks{"allocCode"} .= "subroutine Alloc_Array_".$typeLabel."_LBound(arr,".join(",",map {"size".$_} 1..$dimension).",name,mem_type)\n";
    $postBlocks{"allocCode"} .= "  !% Allocate a ".$dimension."D ".$variableTypeLatex." array and register it.\n";
    $postBlocks{"allocCode"} .= "  use Galacticus_Error\n";
    $postBlocks{"allocCode"} .= "  implicit none\n";
    $c1_width = 16;
    if ( length($typeName) > $c1_width ) {$c1_width = length($typeName)};
    $pad = " " x ($c1_width-length($typeName));
    $postBlocks{"allocCode"} .= "  ".$typeName.",".$pad." allocatable                 :: arr(".join(",",split(//,":" x $dimension)).")\n";
    $pad = " " x ($c1_width-length("integer"));
    $postBlocks{"allocCode"} .= "  integer,".$pad." intent(in),    optional     :: mem_type\n";
    $postBlocks{"allocCode"} .= "  integer,".$pad." intent(in),    dimension(2) :: ".join(",",map {"size".$_} 1..$dimension)."\n";
    $pad = " " x ($c1_width-length("character(len=*)"));
    $postBlocks{"allocCode"} .= "  character(len=*),".$pad." intent(in)                  :: name\n";
    $pad = " " x ($c1_width-length("integer")+1);
    $postBlocks{"allocCode"} .= "  integer,".$pad." save                       :: alloc_err,imem_type,tsize\n";
    $postBlocks{"allocCode"} .= "  !\$omp threadprivate(alloc_err,imem_type,tsize)\n\n";
    $postBlocks{"allocCode"} .= "  if (allocated(arr)) call Dealloc_Array_".$typeLabel."(arr,imem_type)\n";
    $postBlocks{"allocCode"} .= "  allocate(arr(".join(",",map {"size".$_."(1):size".$_."(2)"} 1..$dimension)."),stat=alloc_err)\n";
    $postBlocks{"allocCode"} .= "  if (alloc_err.ne.0) then\n";
    $postBlocks{"allocCode"} .= "     write (0,*) 'Alloc_Array_".$typeLabel."_LBound(): FATAL - unable to allocate memory [',trim(name),']'\n";
    $postBlocks{"allocCode"} .= "     call Galacticus_Error_Report\n";
    $postBlocks{"allocCode"} .= "  end if\n";
    $postBlocks{"allocCode"} .= "  if (present(mem_type)) then\n";
    $postBlocks{"allocCode"} .= "     imem_type=mem_type\n";
    $postBlocks{"allocCode"} .= "  else\n";
    $postBlocks{"allocCode"} .= "     imem_type=MEM_TYPE_MISC\n";
    $postBlocks{"allocCode"} .= "  end if\n";
    $postBlocks{"allocCode"} .= "  tsize=".join("*",map {"(size".$_."(2)-size".$_."(1)+1)"} 1..$dimension)."\n";
    $postBlocks{"allocCode"} .= "  call Register_Memory(arr,imem_type,tsize)\n";
    $postBlocks{"allocCode"} .= "  return\n";
    $postBlocks{"allocCode"} .= "end subroutine Alloc_Array_".$typeLabel."_LBound\n\n";
    
}

# Add closing statements to the various interface blocks.
$preBlocks{"registerInterfaceCode"} .= "end interface\n\n";
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

# Create a file containg all necessary "use" statements.
open(useHandle,">./work/build/utility.memory_management.use.inc.tmp");
print useHandle join("\n",map {"   use ".$_} keys(%use_modules))."\n";
close(useHandle);

# Truncate lines and update old files..
foreach $file ( "utility.memory_management.precontain.inc", "utility.memory_management.postcontain.inc", "utility.memory_management.use.inc" ) {
    &Fortran_Utils::Truncate_Fortran_Lines("./work/build/".$file.".tmp");
    &File_Changes::Update("./work/build/".$file ,"./work/build/".$file.".tmp" );
}

exit;

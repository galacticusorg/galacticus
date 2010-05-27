!% Contains a module which implements sorting sequences.

module Sort
  !% Implements sorting.
  use FGSL
  use, intrinsic :: iso_c_binding
  private
  public :: Sort_Do

  interface Sort_Do
     module procedure Sort_Do_Double
     module procedure Sort_Do_Integer
  end interface

contains

  subroutine Sort_Do_Double(array)
    !% Given an unsorted double precision {\tt array}, sorts it in place.
    implicit none
    double precision, intent(inout), dimension(:) :: array
    
    call Sort_Do_Double_C(size(array),array)
    return
  end subroutine Sort_Do_Double

  subroutine Sort_Do_Integer(array)
    !% Given an unsorted integer {\tt array}, sorts it in place.
    implicit none
    integer, intent(inout), dimension(:) :: array
    
    call Sort_Do_Integer_C(size(array),array)
    return
  end subroutine Sort_Do_Integer

  subroutine Sort_Do_Double_C(arraySize,array)
    !% Do a double precision sort.
    implicit none
    integer,          intent(in)            :: arraySize
    real(c_double),   intent(inout), target :: array(arraySize)
    integer(c_size_t)                       :: arraySizeC
    type(c_ptr)                             :: arrayPointer
    
    arrayPointer=c_loc(array)
    arraySizeC=arraySize
    call FGSL_HeapSort(arrayPointer,arraySizeC,FGSL_SizeOf(1.0d0),Compare_Double)
    return
  end subroutine Sort_Do_Double_C

  subroutine Sort_Do_Integer_C(arraySize,array)
    !% Do a integer sort.
    implicit none
    integer,          intent(in)            :: arraySize
    integer(c_int),   intent(inout), target :: array(arraySize)
    integer(c_size_t)                       :: arraySizeC
    type(c_ptr)                             :: arrayPointer
    
    arrayPointer=c_loc(array)
    arraySizeC=arraySize
    call FGSL_HeapSort(arrayPointer,arraySizeC,FGSL_SizeOf(1),Compare_Integer)
    return
  end subroutine Sort_Do_Integer_C

  function Compare_Double(x,y) bind(c)
    !% Comparison function for double precision data.
    type(c_ptr),    value   :: x,y
    integer(c_int)          :: Compare_Double
    real(c_double), pointer :: rx,ry

    call c_f_pointer(x,rx)
    call c_f_pointer(y,ry)
    Compare_Double=0
    if (rx < ry) Compare_Double=-1
  end function Compare_Double

  function Compare_Integer(x,y) bind(c)
    !% Comparison function for integer data.
    type(c_ptr),    value   :: x,y
    integer(c_int)          :: Compare_Integer
    integer(c_int), pointer :: rx,ry

    call c_f_pointer(x,rx)
    call c_f_pointer(y,ry)
    Compare_Integer=0
    if (rx < ry) Compare_Integer=-1
  end function Compare_Integer

end module Sort

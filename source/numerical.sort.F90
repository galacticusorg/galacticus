!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module which implements sorting sequences.

module Sort
  !% Implements sorting.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Sort_Do, Sort_Index_Do

  interface Sort_Index_Do
     !% Generic interface to index sort routines.
     module procedure Sort_Index_Do_Integer8
     module procedure Sort_Index_Do_Integer
     module procedure Sort_Index_Do_Double
  end interface

  interface Sort_Do
     !% Generic interface to in-place sort routines.
     module procedure Sort_Do_Double
     module procedure Sort_Do_Integer
  end interface

contains

  subroutine Sort_Do_Double(array)
    !% Given an unsorted double precision {\tt array}, sorts it in place.
    implicit none
    double precision, dimension(:), intent(inout) :: array

    call Sort_Do_Double_C(size(array),array)
    return
  end subroutine Sort_Do_Double

  subroutine Sort_Do_Integer(array)
    !% Given an unsorted integer {\tt array}, sorts it in place.
    implicit none
    integer, dimension(:), intent(inout) :: array

    call Sort_Do_Integer_C(size(array),array)
    return
  end subroutine Sort_Do_Integer

  function Sort_Index_Do_Integer8(array)
    !% Given an unsorted integer {\tt array}, sorts it in place.
    use Kind_Numbers
    implicit none
    integer(kind=kind_int8), dimension(:)          , intent(in   ) :: array
    integer(kind=c_size_t ), dimension(size(array))                :: Sort_Index_Do_Integer8

    call Sort_Index_Do_Integer8_C(size(array),array,Sort_Index_Do_Integer8)
    Sort_Index_Do_Integer8=Sort_Index_Do_Integer8+1
    return
  end function Sort_Index_Do_Integer8

  function Sort_Index_Do_Integer(array)
    !% Given an unsorted integer {\tt array}, sorts it in place.
    use Kind_Numbers
    implicit none
    integer(kind=kind_int4), dimension(:)          , intent(in   ) :: array
    integer(kind=c_size_t ), dimension(size(array))                :: Sort_Index_Do_Integer

    call Sort_Index_Do_Integer_C(size(array),array,Sort_Index_Do_Integer)
    Sort_Index_Do_Integer=Sort_Index_Do_Integer+1
    return
  end function Sort_Index_Do_Integer

  function Sort_Index_Do_Double(array)
    !% Given an unsorted double {\tt array}, sorts it in place.
    use Kind_Numbers
    implicit none
    real   (kind=c_double), dimension(:)          , intent(in   ) :: array
    integer(kind=c_size_t), dimension(size(array))                :: Sort_Index_Do_Double

    call Sort_Index_Do_Double_C(size(array),array,Sort_Index_Do_Double)
    Sort_Index_Do_Double=Sort_Index_Do_Double+1
    return
  end function Sort_Index_Do_Double

  subroutine Sort_Do_Double_C(arraySize,array)
    !% Do a double precision sort.
    implicit none
    integer               , intent(in   )         :: arraySize
    real   (kind=c_double), intent(inout), target :: array       (arraySize)
    integer(kind=c_size_t)                        :: arraySizeC
    type   (c_ptr        )                        :: arrayPointer

    arrayPointer=c_loc(array)
    arraySizeC=arraySize
    call FGSL_HeapSort(arrayPointer,arraySizeC,FGSL_SizeOf(1.0d0),Compare_Double)
    return
  end subroutine Sort_Do_Double_C

  subroutine Sort_Do_Integer_C(arraySize,array)
    !% Do a integer sort.
    implicit none
    integer               , intent(in   )         :: arraySize
    integer(kind=c_int   ), intent(inout), target :: array       (arraySize)
    integer(kind=c_size_t)                        :: arraySizeC
    type   (c_ptr        )                        :: arrayPointer

    arrayPointer=c_loc(array)
    arraySizeC=arraySize
    call FGSL_HeapSort(arrayPointer,arraySizeC,FGSL_SizeOf(1),Compare_Integer)
    return
  end subroutine Sort_Do_Integer_C

  subroutine Sort_Index_Do_Integer8_C(arraySize,array,idx)
    !% Do a integer sort.
    use Kind_Numbers
    implicit none
    integer                , intent(in   )         :: arraySize
    integer(kind=kind_int8), intent(in   ), target :: array       (arraySize)
    integer(kind=c_size_t ), intent(inout)         :: idx         (arraySize)
    integer(kind=c_size_t )                        :: arraySizeC
    integer                                        :: status
    type   (c_ptr         )                        :: arrayPointer

    arrayPointer=c_loc(array)
    arraySizeC=arraySize
    status=FGSL_HeapSort_Index(idx,arrayPointer,arraySizeC,sizeof(1_kind_int8),Compare_Integer8)
    return
  end subroutine Sort_Index_Do_Integer8_C

  subroutine Sort_Index_Do_Integer_C(arraySize,array,idx)
    !% Do an integer sort.
    use Kind_Numbers
    implicit none
    integer                , intent(in   )         :: arraySize
    integer(kind=kind_int4), intent(in   ), target :: array       (arraySize)
    integer(kind=c_size_t ), intent(inout)         :: idx         (arraySize)
    integer(kind=c_size_t )                        :: arraySizeC
    integer                                        :: status
    type   (c_ptr         )                        :: arrayPointer

    arrayPointer=c_loc(array)
    arraySizeC=arraySize
    status=FGSL_HeapSort_Index(idx,arrayPointer,arraySizeC,sizeof(1_kind_int4),Compare_Integer)
    return
  end subroutine Sort_Index_Do_Integer_C

  subroutine Sort_Index_Do_Double_C(arraySize,array,idx)
    !% Do an double sort.
    use Kind_Numbers
    implicit none
    integer                , intent(in   )         :: arraySize
    real   (c_double)      , intent(in   ), target :: array       (arraySize)
    integer(kind=c_size_t ), intent(inout)         :: idx         (arraySize)
    integer(kind=c_size_t )                        :: arraySizeC
    integer                                        :: status
    type   (c_ptr         )                        :: arrayPointer

    arrayPointer=c_loc(array)
    arraySizeC=arraySize
    status=FGSL_HeapSort_Index(idx,arrayPointer,arraySizeC,sizeof(1_c_double),Compare_Double)
    return
  end subroutine Sort_Index_Do_Double_C

  function Compare_Double(x,y) bind(c)
    !% Comparison function for double precision data.
    type   (c_ptr        ), value   :: x             , y
    integer(kind=c_int   )          :: Compare_Double
    real   (kind=c_double), pointer :: rx            , ry

    call c_f_pointer(x,rx)
    call c_f_pointer(y,ry)
    Compare_Double=0
    if (rx < ry) Compare_Double=-1
  end function Compare_Double

  function Compare_Integer(x,y) bind(c)
    !% Comparison function for integer data.
    type   (c_ptr     ), value   :: x              , y
    integer(kind=c_int)          :: Compare_Integer
    integer(kind=c_int), pointer :: rx             , ry

    call c_f_pointer(x,rx)
    call c_f_pointer(y,ry)
    Compare_Integer=0
    if (rx < ry) Compare_Integer=-1
  end function Compare_Integer

  function Compare_Integer8(x,y) bind(c)
    !% Comparison function for integer data.
    type   (c_ptr      ), value   :: x               , y
    integer(kind=c_int )          :: Compare_Integer8
    integer(kind=c_long), pointer :: rx              , ry

    call c_f_pointer(x,rx)
    call c_f_pointer(y,ry)
    Compare_Integer8=0
    if (rx < ry) Compare_Integer8=-1
  end function Compare_Integer8

end module Sort

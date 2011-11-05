!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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

  function Sort_Index_Do_Integer8(array)
    !% Given an unsorted integer {\tt array}, sorts it in place.
    use Kind_Numbers
    implicit none
    integer(kind=kind_int8), intent(in), dimension(:)           :: array
    integer,                             dimension(size(array)) :: Sort_Index_Do_Integer8
 
    call Sort_Index_Do_Integer8_C(size(array),array,Sort_Index_Do_Integer8)
    Sort_Index_Do_Integer8=Sort_Index_Do_Integer8+1
    return
  end function Sort_Index_Do_Integer8

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

  subroutine Sort_Index_Do_Integer8_C(arraySize,array,idx)
    !% Do a integer sort.
    use Kind_Numbers
    implicit none
    integer,              intent(in)            :: arraySize
    integer(c_long_long), intent(in),    target :: array(arraySize)
    integer(c_int),       intent(inout)         :: idx(arraySize)
    integer(c_size_t)                           :: arraySizeC
    integer                                     :: status
    type(c_ptr)                                 :: arrayPointer
    
    arrayPointer=c_loc(array)
    arraySizeC=arraySize
    status=FGSL_HeapSort_Index(idx,arrayPointer,arraySizeC,sizeof(1_kind_int8),Compare_Integer8)
    return
  end subroutine Sort_Index_Do_Integer8_C

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

  function Compare_Integer8(x,y) bind(c)
    !% Comparison function for integer data.
    type(c_ptr),     value   :: x,y
    integer(c_int)           :: Compare_Integer8
    integer(c_long), pointer :: rx,ry

    call c_f_pointer(x,rx)
    call c_f_pointer(y,ry)
    Compare_Integer8=0
    if (rx < ry) Compare_Integer8=-1
  end function Compare_Integer8

end module Sort

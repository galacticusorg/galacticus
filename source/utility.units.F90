!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!!{
Contains a module that defines the \mono{unitType} derived type for encoding unit metadata
written to HDF5 output datasets.
!!}

module Units_MetaData
  !!{
  Defines the \mono{unitType} derived type and associated constructor for encoding unit metadata as a compound HDF5
  attribute. Each instance stores:
  \begin{description}
    \item[\mono{unitsInSI}] Multiplicative conversion factor to SI units.
    \item[\mono{description}] Human-readable units description (e.g.\ ``Solar masses'').
    \item[\mono{quantity}] \href{https://docs.astropy.org/en/stable/units/}{astropy}-parseable
      units string (e.g.\ \mono{Msun}).
    \item[\mono{isComoving}] \mono{0} for physical units, \mono{1} for comoving units.
  \end{description}
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_char, c_double, c_int
  implicit none
  private
  public :: unitType

  !! Maximum length of string fields inside unitType. This is the length used for the string members of the HDF5 compound
  !! datatype built in IO_HDF5_Write_Attribute_Units_Scalar (utility.IO.HDF5.F90).
  integer, parameter, public :: unitStringLength=512

  type, bind(C) :: unitType
     !!{
     A derived type that holds unit metadata for a single dataset property. The \mono{bind(C)} attribute ensures a predictable
     memory layout so that the HDF5 compound-type member offsets obtained from it (via \mono{H5OFFSETOF}) when writing units
     attributes are exact.
     !!}
     real     (c_double)                              :: unitsInSI  =0.0_c_double
     character(c_char  ), dimension(unitStringLength) :: description=c_char_""
     character(c_char  ), dimension(unitStringLength) :: quantity   =c_char_""
     integer  (c_int   )                              :: isComoving =0_c_int
  end type unitType

  interface unitType
     module procedure unitsConstructor
  end interface unitType

contains

  function unitsConstructor(unitsInSI,description,quantity,isComoving) result(units)
    !!{
    Construct a \mono{unitType} value from its components. All arguments except \mono{unitsInSI} are optional and default to
    empty strings / zero. Uses plain Fortran intrinsic types (not C-interoperable types) in the interface (double precision and
    default logical) for convenience at call sites.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_char, c_int, c_null_char
    implicit none
    type            (unitType)                       :: units
    double precision          , intent(in)           :: unitsInSI
    character       (len=*   ), intent(in), optional :: description
    character       (len=*   ), intent(in), optional :: quantity
    logical                   , intent(in), optional :: isComoving
    integer                                          :: i          , lengthString

    units%unitsInSI =unitsInSI
    units%isComoving=0
    if (present(isComoving)) then
       if (isComoving) then
          units%isComoving=1_c_int
       else
          units%isComoving=0_c_int
       end if
    end if
    ! Initialise character arrays to null bytes then copy content.
    units%description=c_null_char
    units%quantity   =c_null_char
    if (present(description)) then
       lengthString=min(len_trim(description),unitStringLength-1)
       do i=1,lengthString
          units%description(i)=description(i:i)
       end do
       units%description(lengthString+1)=c_null_char
    end if
    if (present(quantity)) then
       lengthString=min(len_trim(quantity),unitStringLength-1)
       do i=1,lengthString
          units%quantity(i)=quantity(i:i)
       end do
       units%quantity(lengthString+1)=c_null_char
    end if
    return
  end function unitsConstructor

end module Units_MetaData

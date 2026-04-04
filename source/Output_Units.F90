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
Contains a module that defines the \reftype{unitType} derived type for encoding unit metadata
written to HDF5 output datasets.
!!}

module Output_Units
  !!{
  Defines the \reftype{unitType} derived type and associated constructor for encoding unit
  metadata as a compound HDF5 attribute.  Each instance stores:
  \begin{description}
    \item[\texttt{unitsInSI}] Multiplicative conversion factor to SI units.
    \item[\texttt{description}] Human-readable units description (e.g.\ ``Solar masses'').
    \item[\texttt{quantity}] \href{https://docs.astropy.org/en/stable/units/}{astropy}-parseable
      units string (e.g.\ \texttt{Msun}).
    \item[\texttt{isComoving}] 0 for physical units, 1 for comoving units.
  \end{description}
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_char, c_double, c_int
  implicit none
  private

  !! Maximum length of string fields inside unitType.  This must match the length
  !! used in hdf5_cTypes.c when the compound HDF5 type is built.
  integer, parameter, public :: unitStringLength = 512

  type, public, bind(C) :: unitType
     !!{
     A \href{https://en.wikipedia.org/wiki/Derived_type}{derived type} that holds unit metadata
     for a single dataset property.  The \texttt{bind(C)} attribute ensures a predictable
     memory layout so that the HDF5 C compound-type offsets computed with \texttt{offsetof}
     in \texttt{hdf5\_cTypes.c} are exact.
     !!}
     real   (c_double)                          :: unitsInSI   = 0.0d0
     character(c_char), dimension(unitStringLength) :: description = c_char_""
     character(c_char), dimension(unitStringLength) :: quantity    = c_char_""
     integer(c_int)                             :: isComoving  = 0_c_int
  end type unitType

  public :: unitsMake

contains

  function unitsMake(unitsInSI, description, quantity, isComoving) result(u)
    !!{
    Construct a \reftype{unitType} value from its components.  All arguments except
    \texttt{unitsInSI} are optional and default to empty strings / zero.
    Uses plain Fortran intrinsic types in the interface (double precision and default integer) for
    convenience at call sites.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_char, c_int, c_null_char
    implicit none
    type            (unitType        )                :: u
    double precision                 , intent(in)     :: unitsInSI
    character       (len=*           ), intent(in), optional :: description
    character       (len=*           ), intent(in), optional :: quantity
    integer                          , intent(in), optional  :: isComoving
    integer                                          :: i, slen

    u%unitsInSI  = unitsInSI
    u%isComoving = 0_c_int
    if (present(isComoving)) u%isComoving = int(isComoving, c_int)
    ! Initialise character arrays to null bytes then copy content.
    u%description = c_null_char
    u%quantity    = c_null_char
    if (present(description)) then
       slen = min(len_trim(description), unitStringLength-1)
       do i = 1, slen
          u%description(i) = description(i:i)
       end do
       u%description(slen+1) = c_null_char
    end if
    if (present(quantity)) then
       slen = min(len_trim(quantity), unitStringLength-1)
       do i = 1, slen
          u%quantity(i) = quantity(i:i)
       end do
       u%quantity(slen+1) = c_null_char
    end if
    return
  end function unitsMake

end module Output_Units

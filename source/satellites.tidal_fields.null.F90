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
  Implements a satellite tidal field class which assumes zero tidal field.
  !!}

  !![
  <satelliteTidalField name="satelliteTidalFieldNull">
   <description>
    A satellite tidal field class which assumes a zero tidal field always.
   </description>
  </satelliteTidalField>
  !!]
  type, extends(satelliteTidalFieldClass) :: satelliteTidalFieldNull
     !!{
     Implementation of a satellite tidal friction class which assumes no tidal field.
     !!}
     private
   contains
     procedure :: tidalTensorRadial => nullTidalTensorRadial
  end type satelliteTidalFieldNull

  interface satelliteTidalFieldNull
     !!{
     Constructors for the null satellite tidal field class.
     !!}
     module procedure nullConstructorParameters
  end interface satelliteTidalFieldNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteTidalFieldNull} satellite tidal field class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteTidalFieldNull)                :: self
    type(inputParameters        ), intent(inout) :: parameters

    self=satelliteTidalFieldNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  double precision function nullTidalTensorRadial(self,node)
    !!{
    Return the radial part of the tidal tensor for satellite halos assumed to be zero.
    !!}
    implicit none
    class(satelliteTidalFieldNull), intent(inout) :: self
    type (treeNode               ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    nullTidalTensorRadial=0.0d0
    return
  end function nullTidalTensorRadial

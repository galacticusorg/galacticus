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
  Implementation of a cooling rate class which modifies another cooling rate by cutting off cooling in satellites.
  !!}


  !![
  <coolingRate name="coolingRateNoCoolingSatellites">
   <description>A cooling rate class which modifies another cooling rate by cutting off cooling in satellites</description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateNoCoolingSatellites
     !!{
     Implementation of cooling rate class which modifies another cooling rate by cutting off cooling in satellites.
     !!}
     private
     class(coolingRateClass), pointer :: coolingRate_ => null()
   contains
     final     ::         noCoolingSatellitesDestructor
     procedure :: rate => noCoolingSatellitesRate
  end type coolingRateNoCoolingSatellites

  interface coolingRateNoCoolingSatellites
     !!{
     Constructors for the cut off cooling rate class.
     !!}
     module procedure noCoolingSatellitesConstructorParameters
     module procedure noCoolingSatellitesConstructorInternal
  end interface coolingRateNoCoolingSatellites

contains

  function noCoolingSatellitesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the cut off cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingRateNoCoolingSatellites)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(coolingRateClass              ), pointer       :: coolingRate_

    !![
    <objectBuilder class="coolingRate" name="coolingRate_" source="parameters"/>
    !!]
    self=coolingRateNoCoolingSatellites(coolingRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRate_"/>
    !!]
    return
  end function noCoolingSatellitesConstructorParameters

  function noCoolingSatellitesConstructorInternal(coolingRate_) result(self)
    !!{
    Internal constructor for the cut off cooling rate class.
    !!}
    type (coolingRateNoCoolingSatellites)                        :: self
    class(coolingRateClass              ), intent(in   ), target :: coolingRate_

    !![
    <constructorAssign variables="*coolingRate_"/>
    !!]
    return
  end function noCoolingSatellitesConstructorInternal

  subroutine noCoolingSatellitesDestructor(self)
    !!{
    Destructor for the cut off cooling rate class.
    !!}
    implicit none
    type(coolingRateNoCoolingSatellites), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRate_"/>
    !!]
    return
  end subroutine noCoolingSatellitesDestructor

  double precision function noCoolingSatellitesRate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for a model in which this rate is cut off
    in satellites.
    !!}
    implicit none
    class(coolingRateNoCoolingSatellites), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node

    if (node%isSatellite()) then
       noCoolingSatellitesRate=0.0d0
    else
       noCoolingSatellitesRate=self%coolingRate_%rate(node)
    end if
    return
  end function noCoolingSatellitesRate


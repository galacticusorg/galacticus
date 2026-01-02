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
  Implementation of a star formation rate in galactic spheroids which sets rates in satellites to zero.
  !!}

  !![
  <starFormationRateSpheroids name="starFormationRateSpheroidsCentralsOnly">
   <description>A star formation rate in galactic spheroids which sets rates in satellites to zero.</description>
  </starFormationRateSpheroids>
  !!]
  type, extends(starFormationRateSpheroidsClass) :: starFormationRateSpheroidsCentralsOnly
     !!{
     Implementation of a rate for star formation in galactic spheroids which computes the rate by integrating a star formation rate
     over the spheroid.
     !!}
     private
     class(starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_ => null()
   contains
     final     ::         centralsOnlyDestructor
     procedure :: rate => centralsOnlyRate
  end type starFormationRateSpheroidsCentralsOnly

  interface starFormationRateSpheroidsCentralsOnly
     !!{
     Constructors for the \refClass{starFormationRateSpheroidsCentralsOnly} star formation rate in spheroids class.
     !!}
     module procedure centralsOnlyConstructorParameters
     module procedure centralsOnlyConstructorInternal
  end interface starFormationRateSpheroidsCentralsOnly

contains

  function centralsOnlyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationRateSpheroidsCentralsOnly} star formation rate in spheroids class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (starFormationRateSpheroidsCentralsOnly)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(starFormationRateSpheroidsClass       ), pointer       :: starFormationRateSpheroids_

    !![
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    !!]
    self=starFormationRateSpheroidsCentralsOnly(starFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateSpheroids_"/>
    !!]
    return
  end function centralsOnlyConstructorParameters

  function centralsOnlyConstructorInternal(starFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the \refClass{starFormationRateSpheroidsCentralsOnly} star formation rate in spheroids class.
    !!}
    implicit none
    type (starFormationRateSpheroidsCentralsOnly)                        :: self
    class(starFormationRateSpheroidsClass       ), intent(in   ), target :: starFormationRateSpheroids_
    !![
    <constructorAssign variables="*starFormationRateSpheroids_"/>
    !!]

    return
  end function centralsOnlyConstructorInternal

  subroutine centralsOnlyDestructor(self)
    !!{
    Destructor for the \refClass{starFormationRateSpheroidsCentralsOnly} star formation rate in spheroids class.
    !!}
    implicit none
    type(starFormationRateSpheroidsCentralsOnly), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateSpheroids_" />
    !!]
    return
  end subroutine centralsOnlyDestructor

  double precision function centralsOnlyRate(self,node)
    !!{
    Returns the star formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic spheroid of {\normalfont \ttfamily
    node}. Assumes zero rate for satellites, falling through to another class for centrals.
    !!}
    implicit none
    class(starFormationRateSpheroidsCentralsOnly), intent(inout), target :: self
    type (treeNode                              ), intent(inout), target :: node
    
    if (node%isSatellite()) then
       centralsOnlyRate=0.0d0
    else
       centralsOnlyRate=self%starFormationRateSpheroids_%rate(node)
    end if
    return
  end function centralsOnlyRate

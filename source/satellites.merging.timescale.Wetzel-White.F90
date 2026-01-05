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

  !+    Contributions to this file made by:  Martin White.

  !!{
  Implements a satellite merging timescale class which uses the \cite{wetzel_what_2010} method.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesWetzelWhite2010">
   <description>
    A satellite merging timescale class which computes merging timescales using the dynamical friction calibration of
    \cite{wetzel_what_2010}.
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesWetzelWhite2010
     !!{
     A class implementing the \cite{wetzel_what_2010} method for satellite merging timescales.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: timescaleMultiplier
   contains
     final     ::                     wetzelWhite2010Destructor
     procedure :: timeUntilMerging => wetzelWhite2010TimeUntilMerging
  end type satelliteMergingTimescalesWetzelWhite2010

  interface satelliteMergingTimescalesWetzelWhite2010
     !!{
     Constructors for the \cite{lacey_merger_1993} merging timescale class.
     !!}
     module procedure wetzelWhite2010ConstructorParameters
     module procedure wetzelWhite2010ConstructorInternal
  end interface satelliteMergingTimescalesWetzelWhite2010

contains

  function wetzelWhite2010ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{wetzel_what_2010} merging timescale class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteMergingTimescalesWetzelWhite2010)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    double precision                                                           :: timescaleMultiplier

    !![
    <inputParameter>
      <name>timescaleMultiplier</name>
      <defaultValue>0.75d0</defaultValue>
      <description>A multiplier for the merging timescale in dynamical friction timescale calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=satelliteMergingTimescalesWetzelWhite2010(timescaleMultiplier,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function wetzelWhite2010ConstructorParameters

  function wetzelWhite2010ConstructorInternal(timescaleMultiplier,cosmologyFunctions_) result(self)
    !!{
    Constructor for the \cite{wetzel_what_2010} merging timescale class.
    !!}
    implicit none
    type            (satelliteMergingTimescalesWetzelWhite2010)                        :: self
    double precision                                           , intent(in   )         :: timescaleMultiplier
    class           (cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="timescaleMultiplier, *cosmologyFunctions_"/>
    !!]

    return
  end function wetzelWhite2010ConstructorInternal

  subroutine wetzelWhite2010Destructor(self)
    !!{
    Destructor for the \refClass{satelliteMergingTimescalesWetzelWhite2010} satellite merging timescale class.
    !!}
    implicit none
    type(satelliteMergingTimescalesWetzelWhite2010), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine wetzelWhite2010Destructor

  double precision function wetzelWhite2010TimeUntilMerging(self,node,orbit)
    !!{
    Return the timescale for merging satellites using the \cite{wetzel_what_2010} method.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Kepler_Orbits   , only : keplerOrbit
    implicit none
    class           (satelliteMergingTimescalesWetzelWhite2010), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    type            (keplerOrbit                              ), intent(inout) :: orbit
    type            (treeNode                                 ), pointer       :: nodeHost
    class           (nodeComponentBasic                       ), pointer       :: basicHost                   , basic
    double precision                                           , parameter     :: timeScaleNormalization=0.2d0        !   C_dyn from Wetzel & White (2010).
    double precision                                                           :: massRatio
    !$GLC attributes unused :: self, orbit

    ! Find the host node.
    if (node%isSatellite()) then
       nodeHost => node%parent
    else
       nodeHost => node%parent%firstChild
    end if
    ! Compute mass ratio.
    basic     => node    %basic ()
    basicHost => nodeHost%basic ()
    if (node%isSatellite()) then
       ! Node is already a satellite in its host - compute the mass ratio directly.
       massRatio=+basicHost%mass () &
            &    /basic    %mass ()
    else
       ! Node is not yet a satellite in its host - correct the host mass to what it will be after the node becomes a satellite in the
       ! host.
       massRatio=+basicHost%mass () &
            &    /basic    %mass () &
            &    +1.0d0
    end if
    ! Compute dynamical friction timescale using eqn. (2) from Wetzel & White (2010).
    wetzelWhite2010TimeUntilMerging=+self%timescaleMultiplier                               &
         &                          *timeScaleNormalization                                 &
         &                          /self%cosmologyFunctions_%expansionRate  (              &
         &                           self%cosmologyFunctions_%expansionFactor (             &
         &                                                                     basic%time() &
         &                                                                    )             &
         &                                                                   )              &
         &                          *           massRatio                                   &
         &                          /log(1.0d0+massRatio)
    return
  end function wetzelWhite2010TimeUntilMerging

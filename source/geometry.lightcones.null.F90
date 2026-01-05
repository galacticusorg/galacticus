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
  A null implementation of the lightcone geometry class.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Output_Times        , only : outputTimesClass

  !![
  <geometryLightcone name="geometryLightconeNull">
   <description>
    A null implementation of the lightcone geometry class. The lightcone has zero solid angle/volume, so no galaxy ever lies within it.
   </description>
  </geometryLightcone>
  !!]
  type, extends(geometryLightconeClass) :: geometryLightconeNull
     !!{
     A null lightcone geometry class.
     !!}
     private
   contains
     procedure :: timeMinimum               => nullTimeMinimum
     procedure :: timeMaximum               => nullTimeMaximum
     procedure :: isInLightcone             => nullIsInLightcone
     procedure :: replicationCount          => nullReplicationCount
     procedure :: solidAngle                => nullSolidAngle
     procedure :: position                  => nullPosition
     procedure :: velocity                  => nullVelocity
     procedure :: timeLightconeCrossing     => nullTimeLightconeCrossing
     procedure :: positionLightconeCrossing => nullPositionLightconeCrossing
     procedure :: velocityLightconeCrossing => nullVelocityLightconeCrossing
  end type geometryLightconeNull

  interface geometryLightconeNull
     !!{
     Constructors for the \refClass{geometryLightconeNull} dark matter halo spin distribution class.
     !!}
     module procedure nullConstructorParameters
  end interface geometryLightconeNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{geometryLightconeNull} lightcone geometry distribution class which takes a parameter list as
    input.
    !!}
    implicit none
    type(geometryLightconeNull)                 :: self
    type(inputParameters      ), intent(inout)  :: parameters
    !$GLC attributes unused :: parameters
    
    self=geometryLightconeNull()
    return
  end function nullConstructorParameters

  function nullReplicationCount(self,node)
    !!{
    Determine the number of times {\normalfont \ttfamily node} appears in the lightcone.
    !!}
    implicit none
    integer(c_size_t             )                :: nullReplicationCount
    class  (geometryLightconeNull), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    nullReplicationCount=0_c_size_t
    return
  end function nullReplicationCount

  double precision function nullTimeMinimum(self)
    !!{
    Return the minimum time in the lightcone.
    !!}
    implicit none
    class(geometryLightconeNull), intent(inout) :: self
    !$GLC attributes unused :: self

    nullTimeMinimum=huge(0.0d0)
    return
  end function nullTimeMinimum

  double precision function nullTimeMaximum(self)
    !!{
    Return the maximum time in the lightcone.
    !!}
    implicit none
    class(geometryLightconeNull), intent(inout) :: self
    !$GLC attributes unused :: self

    nullTimeMaximum=-huge(0.0d0)
    return
  end function nullTimeMaximum

  logical function nullIsInLightcone(self,node,atPresentEpoch,radiusBuffer)
    !!{
    Determine if the given {\normalfont \ttfamily node} lies within the lightcone
    !!}
    implicit none
    class           (geometryLightconeNull), intent(inout)            :: self
    type            (treeNode             ), intent(inout)            :: node
    logical                                , intent(in   ) , optional :: atPresentEpoch
    double precision                       , intent(in   ) , optional :: radiusBuffer
    !$GLC attributes unused :: self, node, atPresentEpoch, radiusBuffer

    nullIsInLightcone=.false.
    return
  end function nullIsInLightcone

  double precision function nullSolidAngle(self)
    !!{
    Return the solid angle (in steradians) of a null lightcone.
    !!}
    implicit none
    class(geometryLightconeNull), intent(inout) :: self
    !$GLC attributes unused :: self

    nullSolidAngle=0.0d0
    return
  end function nullSolidAngle

  function nullPosition(self,node,instance)
    !!{
    Return the position of the node in lightcone coordinates.
    !!}
    implicit none
    double precision                       , dimension(3)          :: nullPosition
    class           (geometryLightconeNull), intent(inout)         :: self
    type            (treeNode             ), intent(inout), target :: node
    integer         (c_size_t             ), intent(in   )         :: instance
    !$GLC attributes unused :: self, node, instance

    nullPosition=0.0d0
    return
  end function nullPosition

  function nullVelocity(self,node,instance)
    !!{
    Return the velocity of the node in lightcone coordinates.
    !!}
    implicit none
    double precision                       , dimension(3)  :: nullvelocity
    class           (geometryLightconeNull), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    integer         (c_size_t             ), intent(in   ) :: instance
    !$GLC attributes unused :: self, node, instance

    nullVelocity=0.0d0
    return
  end function nullVelocity

  double precision function nullTimeLightconeCrossing(self,node,timeStart,timeEnd,timesCrossing)
    !!{
    Return the time of the next lightcone crossing for this node.
    !!}
    implicit none
    class           (geometryLightconeNull), intent(inout), target                              :: self
    type            (treeNode             ), intent(inout), target                              :: node
    double precision                       , intent(in   )                                      :: timeStart    , timeEnd
    double precision                       , intent(inout), dimension(:), allocatable, optional :: timesCrossing
    !$GLC attributes unused :: self, node, timeStart, timeEnd, timesCrossing

    nullTimeLightconeCrossing=huge(0.0d0)
    return
  end function nullTimeLightconeCrossing

  function nullPositionLightconeCrossing(self,node)
    !!{
    Return the position at the next lightcone crossing for this node.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                       , dimension(3)  :: nullPositionLightconeCrossing
    class           (geometryLightconeNull), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    nullPositionLightconeCrossing=0.0d0
    return
  end function nullPositionLightconeCrossing
  
  function nullVelocityLightconeCrossing(self,node)
    !!{
    Return the velocity at the next lightcone crossing for this node.
    !!}
    implicit none
    double precision                       , dimension(3)  :: nullVelocityLightconeCrossing
    class           (geometryLightconeNull), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    nullVelocityLightconeCrossing=0.0d0
    return
  end function nullVelocityLightconeCrossing

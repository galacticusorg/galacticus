!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements an N-body data operator which computes the distance of each particle from a point.
!!}

  !![
  <nbodyOperator name="nbodyOperatorDistanceFromPoint">
   <description>An N-body data operator which computes the distance of each particle from a point.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorDistanceFromPoint
     !!{
     An N-body data operator which computes the distance of each particle from a point.
     !!}
     private
     double precision, dimension(3) :: point
   contains
     procedure :: operate => distanceFromPointOperate
  end type nbodyOperatorDistanceFromPoint

  interface nbodyOperatorDistanceFromPoint
     !!{
     Constructors for the \refClass{nbodyOperatorDistanceFromPoint} N-body operator class.
     !!}
     module procedure distanceFromPointConstructorParameters
     module procedure distanceFromPointConstructorInternal
  end interface nbodyOperatorDistanceFromPoint

contains

  function distanceFromPointConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorDistanceFromPoint} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorDistanceFromPoint)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                , dimension(3)  :: point

    !![
    <inputParameter>
      <name>point</name>
      <source>parameters</source>
      <description>The Cartesian coordinates of the point from which to compute the distance.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorDistanceFromPoint(point)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function distanceFromPointConstructorParameters

  function distanceFromPointConstructorInternal(point) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorDistanceFromPoint} N-body operator class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorDistanceFromPoint)                              :: self
    double precision                                , intent(in   ), dimension(3) :: point
    !![
    <constructorAssign variables="point"/>
    !!]

    return
  end function distanceFromPointConstructorInternal

  subroutine distanceFromPointOperate(self,simulations)
    !!{
    Compute the distance of each particle from a point.
    !!}
    use :: Display, only : displayIndent, displayUnindent, verbosityLevelStandard
    implicit none
    class           (nbodyOperatorDistanceFromPoint), intent(inout)                 :: self
    type            (nBodyData                     ), intent(inout), dimension(:  ) :: simulations
    double precision                                , pointer      , dimension(:,:) :: position
    double precision                                , pointer      , dimension(:  ) :: distanceFromPoint
    integer         (c_size_t                      )                                :: iSimulation

    call displayIndent('compute distance from point',verbosityLevelStandard)
    do iSimulation=1,size(simulations)
       ! Retrieve required properties.
       position => simulations(iSimulation)%propertiesRealRank1%value('position')
       ! Allocate workspace.
       allocate(distanceFromPoint(size(position,dim=2)))
       ! Compute distances.
       !$omp workshare
       distanceFromPoint=sqrt(                                  &
            &                 +(position(1,:)-self%point(1))**2 &
            &                 +(position(2,:)-self%point(2))**2 &
            &                 +(position(3,:)-self%point(3))**2 &
            &                )
       !$omp end workshare
       ! Store results.
       call simulations(iSimulation)%propertiesReal%set('distanceFromPoint',distanceFromPoint)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine distanceFromPointOperate

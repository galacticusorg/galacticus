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
Implements an N-body data operator which constructs the convex hull of the particles.
!!}
  
  !![
  <nbodyOperator name="nbodyOperatorConvexHull">
   <description>An N-body data operator which constructs the convex hull of the particles.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorConvexHull
     !!{
     An N-body data operator which constructs the convex hull of the particles.
     !!}
     private
   contains
     procedure :: operate => convexHullOperate
  end type nbodyOperatorConvexHull

  interface nbodyOperatorConvexHull
     !!{
     Constructors for the {\normalfont \ttfamily convexHull} N-body operator class.
     !!}
     module procedure convexHullConstructorParameters
  end interface nbodyOperatorConvexHull

contains

  function convexHullConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily convexHull} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nbodyOperatorConvexHull)                :: self
    type(inputParameters        ), intent(inout) :: parameters
    
    self=nbodyOperatorConvexHull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function convexHullConstructorParameters

  subroutine convexHullOperate(self,simulations)
    !!{
    Construct the convex hull of the points.
    !!}
    use :: Display           , only : displayIndent, displayUnindent, verbosityLevelStandard
    use :: Points_Convex_Hull, only : convexHull
    implicit none
    class           (nbodyOperatorConvexHull), intent(inout)                 :: self
    type            (nBodyData              ), intent(inout), dimension(  :) :: simulations
    double precision                         , pointer      , dimension(:,:) :: position
    type            (convexHull             ), pointer                       :: hull
    integer                                                                  :: i
    double precision :: volume

    call displayIndent('construct convex hull',verbosityLevelStandard)
    do i=1,size(simulations)
       position => simulations(i)%propertiesRealRank1%value('position')
       allocate(hull)
       hull=convexHull(position)
       call simulations(i)%attributesGeneric%set('convexHull',hull)
       volume=hull%volume()
       nullify(hull    )
       nullify(position)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine convexHullOperate

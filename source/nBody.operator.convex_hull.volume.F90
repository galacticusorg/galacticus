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
Implements an N-body data operator which computes the convex hull volume of the particles.
!!}
  
  !![
  <nbodyOperator name="nbodyOperatorConvexHullVolume">
   <description>An N-body data operator which computes the convex hull volume of the particles.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorConvexHullVolume
     !!{
     An N-body data operator which computes the convex hull volume of the particles.
     !!}
     private
   contains
     procedure :: operate => convexHullVolumeOperate
  end type nbodyOperatorConvexHullVolume

  interface nbodyOperatorConvexHullVolume
     !!{
     Constructors for the {\normalfont \ttfamily convexHullVolume} N-body operator class.
     !!}
     module procedure convexHullVolumeConstructorParameters
  end interface nbodyOperatorConvexHullVolume

contains

  function convexHullVolumeConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily convexHullVolume} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nbodyOperatorConvexHullVolume)                :: self
    type(inputParameters              ), intent(inout) :: parameters
    
    self=nbodyOperatorConvexHullVolume()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function convexHullVolumeConstructorParameters

  subroutine convexHullVolumeOperate(self,simulations)
    !!{
    Compute the convex hull volume of the points.
    !!}
    use :: Error             , only : Error_Report
    use :: Display           , only : displayIndent, displayUnindent, verbosityLevelStandard
    use :: Points_Convex_Hull, only : convexHull
    implicit none
    class  (nbodyOperatorConvexHullVolume), intent(inout)                 :: self
    type   (nBodyData                    ), intent(inout), dimension(  :) :: simulations
    class  (*                            ), pointer                       :: hull
    integer                                                               :: i

    call displayIndent('compute convex hull volume',verbosityLevelStandard)
    do i=1,size(simulations)
       hull => simulations(i)%attributesGeneric%value('convexHull')
       select type (hull)
       class is (convexHull)
          call simulations(i)%attributesReal%set           (keyCH        ='convexHullVolume',value         =hull%volume())
          call simulations(i)%analysis      %writeAttribute(attributeName='convexHullVolume',attributeValue=hull%volume())
       class default
          call Error_Report('incorrect class'//{introspection:location})
       end select
       nullify(hull)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine convexHullVolumeOperate

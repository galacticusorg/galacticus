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
Implements a survey geometry which combines multiple other surveys.
!!}

  type, public :: surveyGeometryList
     class(surveyGeometryClass), pointer :: surveyGeometry_ => null()
     type (surveyGeometryList ), pointer :: next            => null()
  end type surveyGeometryList

  !![
  <surveyGeometry name="surveyGeometryCombined">
   <description>Implements a survey geometry which combines multiple other surveys.</description>
   <linkedList type="surveyGeometryList" variable="surveyGeometries" next="next" object="surveyGeometry_" objectType="surveyGeometryClass"/>
  </surveyGeometry>
  !!]
  type, extends(surveyGeometryClass) :: surveyGeometryCombined
     private
     type(surveyGeometryList), pointer :: surveyGeometries => null()
   contains
     final      ::                            combinedDestructor
     procedure  :: windowFunctionAvailable => combinedWindowFunctionAvailable
     procedure  :: angularPowerAvailable   => combinedAngularPowerAvailable
     procedure  :: solidAngle              => combinedSolidAngle
     procedure  :: windowFunctions         => combinedWindowFunctions
     procedure  :: angularPower            => combinedAngularPower
     procedure  :: pointIncluded           => combinedPointIncluded
  end type surveyGeometryCombined

  interface surveyGeometryCombined
     module procedure combinedConstructorParameters
     module procedure combinedConstructorInternal
  end interface surveyGeometryCombined

contains

  function combinedConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{surveyGeometryCombined} survey geometry class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (surveyGeometryCombined)                :: self
    type   (inputParameters       ), intent(inout) :: parameters
    type   (surveyGeometryList    ), pointer       :: surveyGeometry_
    integer                                        :: i

    self           %surveyGeometries => null()
    surveyGeometry_                  => null()
    do i=1,parameters%copiesCount('surveyGeometry',zeroIfNotPresent=.true.)
       if (associated(surveyGeometry_)) then
          allocate(surveyGeometry_%next)
          surveyGeometry_ => surveyGeometry_%next
       else
          allocate(self%surveyGeometries)
          surveyGeometry_ => self%surveyGeometries
       end if
       !![
       <objectBuilder class="surveyGeometry" name="surveyGeometry_%surveyGeometry_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="surveyGeometry"/>
    !!]
    return
  end function combinedConstructorParameters

  function combinedConstructorInternal(surveyGeometries) result (self)
    !!{
    Internal constructor for the combined output analysis weight operator class.
    !!}
    implicit none
    type(surveyGeometryCombined)                        :: self
    type(surveyGeometryList    ), target, intent(in   ) :: surveyGeometries
    type(surveyGeometryList    ), pointer               :: surveyGeometry_

    self           %surveyGeometries => surveyGeometries
    surveyGeometry_                  => surveyGeometries
    do while (associated(surveyGeometry_))
       !![
       <referenceCountIncrement owner="surveyGeometry_" object="surveyGeometry_"/>
       !!]
       surveyGeometry_ => surveyGeometry_%next
    end do
    return
  end function combinedConstructorInternal

  subroutine combinedDestructor(self)
    !!{
    Destructor for the combined weight operator class.
    !!}
    implicit none
    type(surveyGeometryCombined), intent(inout) :: self
    type(surveyGeometryList    ), pointer       :: surveyGeometry_, surveyGeometryNext

    if (associated(self%surveyGeometries)) then
       surveyGeometry_ => self%surveyGeometries
       do while (associated(surveyGeometry_))
          surveyGeometryNext => surveyGeometry_%next
          !![
          <objectDestructor name="surveyGeometry_%surveyGeometry_"/>
          !!]
          deallocate(surveyGeometry_)
          surveyGeometry_ => surveyGeometryNext
       end do
    end if
    return
  end subroutine combinedDestructor

  logical function combinedWindowFunctionAvailable(self)
    !!{
    Return false to indicate that survey window function is not available.
    !!}
    implicit none
    class(surveyGeometryCombined), intent(inout) :: self
    !$GLC attributes unused :: self

    combinedWindowFunctionAvailable=.false.
    return
  end function combinedWindowFunctionAvailable

  logical function combinedAngularPowerAvailable(self)
    !!{
    Return false to indicate that survey angular power is not available.
    !!}
    implicit none
    class(surveyGeometryCombined), intent(inout) :: self
    !$GLC attributes unused :: self

    combinedAngularPowerAvailable=.false.
    return
  end function combinedAngularPowerAvailable

  double precision function combinedSolidAngle(self,field)
    !!{
    Return the survey solid angle.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (surveyGeometryCombined), intent(inout)               :: self
    integer                        , intent(in   ), optional     :: field
    !$GLC attributes unused :: self, field

    combinedSolidAngle=0.0d0
    call Error_Report('solid angle is not supported'//{introspection:location})
    return
  end function combinedSolidAngle

  subroutine combinedWindowFunctions(self,mass1,mass2,gridCount,boxLength,windowFunction1,windowFunction2)
    !!{
    Provides window functions for combined survey geometries.
    !!}
    use            :: Error        , only : Error_Report
    use, intrinsic :: ISO_C_Binding, only : c_double_complex
    implicit none
    class           (surveyGeometryCombined), intent(inout)                                           :: self
    double precision                        , intent(in   )                                           :: mass1,mass2
    integer                                 , intent(in   )                                           :: gridCount
    double precision                        , intent(  out)                                           :: boxLength
    complex         (c_double_complex      ), intent(  out), dimension(gridCount,gridCount,gridCount) :: windowFunction1,windowFunction2
    !$GLC attributes unused :: self, mass1, mass2, gridCount, boxLength, windowFunction1, windowFunction2

    call Error_Report('window function construction is not supported'//{introspection:location})
    return
  end subroutine combinedWindowFunctions

  double precision function combinedAngularPower(self,i,j,l)
    !!{
    Return the survey angular power $C^{ij}_\ell$ from \gls{mangle} polygons.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (surveyGeometryCombined), intent(inout):: self
    integer                                 , intent(in   ):: i          , j, &
         &                                                    l
    !$GLC attributes unused ::self, i, j, l

    combinedAngularPower=0.0d0
    call Error_Report('angular power is not supported'//{introspection:location})
    return
  end function combinedAngularPower

  logical function combinedPointIncluded(self,point,mass)
    !!{
    Return true if a point is included in the combined survey geometry.
    !!}
    implicit none
    class           (surveyGeometryCombined), intent(inout)               :: self
    double precision                        , intent(in   ), dimension(3) :: point
    double precision                        , intent(in   )               :: mass
    type            (surveyGeometryList    ), pointer                     :: surveyGeometry_

    combinedPointIncluded =  .false.
    surveyGeometry_       => self%surveyGeometries
    do while (associated(surveyGeometry_))
       combinedPointIncluded =  surveyGeometry_%surveyGeometry_%pointIncluded(point,mass)
       if (combinedPointIncluded) exit
       surveyGeometry_       => surveyGeometry_%next
    end do
    return
  end function combinedPointIncluded

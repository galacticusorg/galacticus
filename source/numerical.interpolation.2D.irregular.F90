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
Contains a module which implements a convenient interface to the {\normalfont \ttfamily BIVAR} 2D interpolation on irregularly spaced points
package.
!!}

module Numerical_Interpolation_2D_Irregular
  !!{
  Implements a convenient interface to the {\normalfont \ttfamily BIVAR} 2D interpolation on irregularly spaced points package.
  !!}
  implicit none
  private
  public :: Interpolate_2D_Irregular, interp2dIrregularObject

  ! A derived type used for storing workspace for the interpolation.
  type interp2dIrregularObject
     integer         , allocatable, dimension(:), private :: integerWork
     double precision, allocatable, dimension(:), private :: realWork
  end type interp2dIrregularObject

  interface Interpolate_2D_Irregular
     module procedure Interpolate_2D_Irregular_Array
     module procedure Interpolate_2D_Irregular_Scalar
  end interface Interpolate_2D_Irregular

contains

  function Interpolate_2D_Irregular_Array(dataX,dataY,dataZ,interpolateX,interpolateY,workspace,numberComputePoints,reset)
    !!{
    Perform interpolation on a set of points irregularly spaced on a 2D surface.
    !!}
    use :: Bivar, only : idbvip
    implicit none
    type            (interp2dIrregularObject)                               , intent(inout)           :: workspace
    double precision                         , dimension(:)                 , intent(in   )           :: dataX                         , dataY                    , &
         &                                                                                               dataZ                         , interpolateX             , &
         &                                                                                               interpolateY
    integer                                                                 , intent(in   ), optional :: numberComputePoints
    logical                                                                 , intent(inout), optional :: reset
    double precision                         , dimension(size(interpolateX))                          :: Interpolate_2D_Irregular_Array
    integer                                                                                           :: dataPointCount                , integerWorkspaceSize     , &
         &                                                                                               interpolatedPointCount        , numberComputePointsActual, &
         &                                                                                               realWorkspaceSize             , resetFlag
    logical                                                                                           :: resetActual

    ! Determine reset status.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    else
       resetActual=.true.
    end if
    if (resetActual) then
       resetFlag=1
    else
       resetFlag=2
    end if
    ! Decide how many points to use for computing partial derivatives.
    if (present(numberComputePoints)) then
       ! Use the specified number.
       numberComputePointsActual=numberComputePoints
    else
       ! Use our default of 5.
       numberComputePointsActual=5
    end if

    ! Get number of data points specified.
    dataPointCount=size(dataX)

    ! Get number of interpolation points specified.
    interpolatedPointCount=size(interpolateX)

    ! Ensure that workspace is sufficient.
    integerWorkspaceSize=max(31,27+numberComputePointsActual)*dataPointCount+interpolatedPointCount
    realWorkspaceSize   =8                                   *dataPointCount
    if (allocated(workspace%integerWork)) then
       if (size(workspace%integerWork) < integerWorkspaceSize) then
          deallocate(workspace%integerWork                      )
          allocate  (workspace%integerWork(integerWorkspaceSize))
          workspace%integerWork=0
       end if
       if (size(workspace%realWork   ) < realWorkspaceSize   ) then
          deallocate(workspace%realWork                         )
          allocate  (workspace%realWork   (realWorkspaceSize   ))
          workspace%realWork=0.0d0
       end if
    else
       allocate(workspace%integerWork(integerWorkspaceSize))
       allocate(workspace%realWork   (realWorkspaceSize   ))
       workspace%integerWork=0
       workspace%realWork   =0.0d0
    end if

    ! Call the subroutine that does the interpolation. This call is wrapped inside an OpenMP critical section as the 2D
    ! interpolation code is not thread parallel - I don't understand why, but it's such old and ugly code that I can't figure it
    ! out. It should be replaced.
    !$omp critical (idbvip)
    call idbvip(resetFlag,numberComputePointsActual,dataPointCount,dataX,dataY,dataZ,interpolatedPointCount,interpolateX&
         &,interpolateY,Interpolate_2D_Irregular_Array,workspace%integerWork,workspace%realWork)
    !$omp end critical (idbvip)
    return
  end function Interpolate_2D_Irregular_Array

  double precision function Interpolate_2D_Irregular_Scalar(dataX,dataY,dataZ,interpolateX,interpolateY,workspace,numberComputePoints,reset)
    !!{
    Perform interpolation on a set of points irregularly spaced on a 2D surface. This version is simply a wrapper that does
    look up for a scalar point by calling the array-based version.
    !!}
    implicit none
    type            (interp2dIrregularObject)              , intent(inout)           :: workspace
    double precision                         , dimension(:), intent(in   )           :: dataX              , dataY            , dataZ
    double precision                                       , intent(in   )           :: interpolateX       , interpolateY
    integer                                                , intent(in   ), optional :: numberComputePoints
    logical                                                , intent(inout), optional :: reset
    double precision                         , dimension(1)                          :: interpolateXArray  , interpolateYArray, interpolateZArray

    interpolateXArray(1)=interpolateX
    interpolateYArray(1)=interpolateY
    interpolateZArray=Interpolate_2D_Irregular_Array(dataX,dataY,dataZ,interpolateXArray,interpolateYArray,workspace,numberComputePoints,reset)
    Interpolate_2D_Irregular_Scalar=interpolateZArray(1)
    return
  end function Interpolate_2D_Irregular_Scalar

end module Numerical_Interpolation_2D_Irregular

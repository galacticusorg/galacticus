!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements a convenient interface to the {\tt BIVAR} 2D interpolation on irregularly spaced points
!% package.

module Numerical_Interpolation_2D_Irregular
  !% Implements a convenient interface to the {\tt BIVAR} 2D interpolation on irregularly spaced points package.
  implicit none
  private
  public :: Interpolate_2D_Irregular, interp2dIrregularObject

  ! Specify an explicit dependence on the bivar.o object file.
  !: ./work/build/Bivar/bivar.o

  ! A derived type used for storing workspace for the interpolation.
  type interp2dIrregularObject
     integer,          allocatable, dimension(:), private :: integerWork
     double precision, allocatable, dimension(:), private :: realWork
  end type interp2dIrregularObject

  interface Interpolate_2D_Irregular
     module procedure Interpolate_2D_Irregular_Array
     module procedure Interpolate_2D_Irregular_Scalar
  end interface

contains

  function Interpolate_2D_Irregular_Array(dataX,dataY,dataZ,interpolateX,interpolateY,workspace,numberComputePoints,reset)
    !% Perform interpolation on a set of points irregularly spaced on a 2D surface.
    use Memory_Management
    implicit none
    type(interp2dIrregularObject), intent(inout)                                :: workspace
    double precision,              intent(in),    dimension(:)                  :: dataX,dataY,dataZ,interpolateX,interpolateY
    integer,                       intent(in),    optional                      :: numberComputePoints
    logical,                       intent(inout), optional                      :: reset
    double precision,                             dimension(size(interpolateX)) :: Interpolate_2D_Irregular_Array
    integer                                                                     :: dataPointCount,interpolatedPointCount&
         &,integerWorkspaceSize,realWorkspaceSize,resetFlag,numberComputePointsActual
    logical                                                                     :: resetActual

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
       resetFlag=0
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
          call Memory_Usage_Record(sizeof(workspace%integerWork),addRemove=-1)
          deallocate(workspace%integerWork                      )
          allocate  (workspace%integerWork(integerWorkspaceSize))
          call Memory_Usage_Record(sizeof(workspace%integerWork),addRemove=+1)
       end if
       if (size(workspace%realWork   ) < realWorkspaceSize   ) then
          call Memory_Usage_Record(sizeof(workspace%realWork),addRemove=-1)
          deallocate(workspace%realWork                         )
          allocate  (workspace%realWork   (realWorkspaceSize   ))
          call Memory_Usage_Record(sizeof(workspace%realWork),addRemove=+1)
       end if
    else
       allocate(workspace%integerWork(integerWorkspaceSize))
       allocate(workspace%realWork   (realWorkspaceSize   ))
       call Memory_Usage_Record(sizeof(workspace%integerWork)+sizeof(workspace%realWork),blockCount=2)
    end if

    ! Call the subroutine that does the interpolation.
    !$omp critical(TwoD_Irregular_Interpolation)
    call idbvip(resetFlag,numberComputePointsActual,dataPointCount,dataX,dataY,dataZ,interpolatedPointCount,interpolateX&
         &,interpolateY,Interpolate_2D_Irregular_Array,workspace%integerWork,workspace%realWork)
    !$omp end critical(TwoD_Irregular_Interpolation)
    
    return
  end function Interpolate_2D_Irregular_Array

  double precision function Interpolate_2D_Irregular_Scalar(dataX,dataY,dataZ,interpolateX,interpolateY,workspace,numberComputePoints,reset)
    !% Perform interpolation on a set of points irregularly spaced on a 2D surface. This version is simply a wrapper that does
    !% look up for a scalar point by calling the array-based version.
    implicit none
    type(interp2dIrregularObject), intent(inout)                             :: workspace
    double precision,              intent(in), dimension(:)                  :: dataX,dataY,dataZ
    double precision,              intent(in)                                :: interpolateX,interpolateY
    integer,                       intent(in), optional                      :: numberComputePoints
    logical,                       intent(in), optional                      :: reset
    double precision,                          dimension(1)                  :: interpolateXArray,interpolateYArray,interpolateZArray
    logical                                                                  :: resetActual
    integer                                                                  :: numberComputePointsActual

    ! Determine reset status.
    if (present(reset)) then
       resetActual=reset
    else
       resetActual=.true.
    end if
  
    ! Decide how many points to use for computing partial derivatives.
    if (present(numberComputePoints)) then
       ! Use the specified number.
       numberComputePointsActual=numberComputePoints
    else
       ! Use our default of 5.
       numberComputePointsActual=5
    end if

    interpolateXArray(1)=interpolateX
    interpolateYArray(1)=interpolateY
    interpolateZArray=Interpolate_2D_Irregular_Array(dataX,dataY,dataZ,interpolateXArray,interpolateYArray,workspace,numberComputePointsActual,resetActual)
    Interpolate_2D_Irregular_Scalar=interpolateZArray(1)

    return
  end function Interpolate_2D_Irregular_Scalar

end module Numerical_Interpolation_2D_Irregular

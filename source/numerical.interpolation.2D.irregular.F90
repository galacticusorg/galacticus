!% Contains a module which implements a convenient interface to the {\tt BIVAR} 2D interpolation on irregularly spaced points
!% package.

module Numerical_Interpolation_2D_Irregular
  !% Implements a convenient interface to the {\tt BIVAR} 2D interpolation on irregularly spaced points package.
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
       !!! AJB: At present, we still reset even if told not to, as not resetting seems to cause BIVAR to crash.
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
          deallocate(workspace%integerWork                      )
          allocate  (workspace%integerWork(integerWorkspaceSize))
       end if
       if (size(workspace%realWork   ) < realWorkspaceSize   ) then
          deallocate(workspace%realWork                         )
          allocate  (workspace%realWork   (realWorkspaceSize   ))
       end if
    else
       allocate(workspace%integerWork(integerWorkspaceSize))
       allocate(workspace%realWork   (realWorkspaceSize   ))
    end if

    ! Call the subroutine that does the interpolation.
    !$omp critical (2D_Irregular_Interpolation)
    call idbvip(resetFlag,numberComputePointsActual,dataPointCount,dataX,dataY,dataZ,interpolatedPointCount,interpolateX&
         &,interpolateY,Interpolate_2D_Irregular_Array,workspace%integerWork,workspace%realWork)
    !$omp end critical (2D_Irregular_Interpolation)
    
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

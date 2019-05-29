!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Contains a module which implements the effects of inclination on disk size in an output analysis distribution operator class.

  use Tables
 
  !# <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorDiskSizeInclntn">
  !#  <description>An output analysis distribution operator class which implements the effects of inclination on disk size.</description>
  !# </outputAnalysisDistributionOperator>
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorDiskSizeInclntn
     !% An output distribution operator class which implements the effects of inclination on disk size.
     private
     type (table1DLinearLinear)              :: inclinationTable
     class(table1D            ), allocatable :: sizeTable
  contains
     procedure :: operateScalar       => diskSizeInclinationOperateScalar
     procedure :: operateDistribution => diskSizeInclinationOperateDistribution
  end type outputAnalysisDistributionOperatorDiskSizeInclntn

  interface outputAnalysisDistributionOperatorDiskSizeInclntn
     !% Constructors for the ``diskSizeInclination'' output distribution operator class.
     module procedure diskSizeInclinationConstructorParameters
     module procedure diskSizeInclinationConstructorInternal
  end interface outputAnalysisDistributionOperatorDiskSizeInclntn

  ! Module-scope variables used in root finding and integrations.
  double precision :: diskSizeInclntnXIntegrate, diskSizeInclntnAngle
  !$omp threadprivate(diskSizeInclntnXIntegrate,diskSizeInclntnAngle)
  
contains

  function diskSizeInclinationConstructorParameters(parameters) result(self)
    !% Constructor for the ``diskSizeInclination'' output analysis distribution operator operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type(outputAnalysisDistributionOperatorDiskSizeInclntn)                :: self
    type(inputParameters                                  ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=outputAnalysisDistributionOperatorDiskSizeInclntn()
    return
  end function diskSizeInclinationConstructorParameters

  function diskSizeInclinationConstructorInternal() result(self)
    !% Internal constructor for the ``diskSizeInclination'' output analysis distribution operator class.
    use Root_Finder
    use Numerical_Constants_Math
    use Table_Labels
    use File_Utilities
    use IO_HDF5
    use ISO_Varying_String
    use Galacticus_Paths
    implicit none
    type            (outputAnalysisDistributionOperatorDiskSizeInclntn)                            :: self
    double precision                                                   , parameter                 :: inclinationAngleEpsilon=1.0d-3
    integer                                                            , parameter                 :: inclinationAngleCount  =100
    type            (rootFinder                                       ), save                      :: finder
    !$omp threadprivate (finder)
    double precision                                                                               :: halfLightRadiusFaceOn         , halfLightRadius
    integer                                                                                        :: i
    double precision                                                   , dimension(:), allocatable :: halfMassRadii
    type            (varying_string                                   )                            :: fileName
    type            (hdf5Object                                       )                            :: file
    type            (lockDescriptor                                   )                            :: lockFileDescriptor

    call self%inclinationTable%create      (                                                                           &
         &                                                                1.0d0                                      , &
         &                                                                inclinationAngleEpsilon                    , &
         &                                  inclinationAngleCount                                                    , &
         &                                  extrapolationType            =[extrapolationTypeFix,extrapolationTypeFix]  &
         &                                 )
    call File_Lock_Initialize(lockFileDescriptor)
    fileName=galacticusPath(pathTypeDataStatic)//"/galacticStructure/diskExponentialInclinedHalfMassRadii.hdf5"
    if (File_Exists(fileName)) then
       call hdf5Access%set()
       call File_Lock(char(fileName),lockFileDescriptor,lockIsShared=.true.)
       call file%openFile(char(fileName),readOnly=.true.)
       call file%readDataset('halfMassRadii',halfMassRadii)
       call file%close()
       call File_Unlock(lockFileDescriptor)
       call hdf5Access%unset()
       call self%inclinationTable%populate(halfMassRadii)
    else
       ! Tabulate dependence of projected half-light radius on disk inclination angle.
       call finder               %tolerance   (                                                                           &
            &                                  toleranceRelative            =1.0d-6                                     , &
            &                                  toleranceAbsolute            =1.0d-6                                       &
            &                                 )
       call finder               %rootFunction(                                                                           &
            &                                                               diskSizeInclntnRoot                           &
            &                                 )
       call finder               %rangeExpand (                                                                           &
            &                                  rangeExpandUpward            =2.0d0                                      , &
            &                                  rangeExpandDownward          =0.5d0                                      , &
            &                                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive              , &
            &                                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative              , &
            &                                  rangeExpandType              =rangeExpandMultiplicative                    &
            &                                 )
       diskSizeInclntnAngle =acos(self%inclinationTable%x(1))
       halfLightRadiusFaceOn=finder%find(rootGuess=1.0d0)
       call self%inclinationTable%populate(0.0d0,1)
       !$omp parallel do private (i,halfLightRadius) copyin (finder)
       do i=2,inclinationAngleCount
          diskSizeInclntnAngle=acos(self%inclinationTable%x(i))
          halfLightRadius =finder%find(rootGuess=1.0d0)
          call self%inclinationTable%populate(log10(halfLightRadius/halfLightRadiusFaceOn),i)
       end do
       !$omp end parallel do
       halfMassRadii=reshape(self%inclinationTable%ys(),[inclinationAngleCount])
       call hdf5Access%set()
       call File_Lock(char(fileName),lockFileDescriptor,lockIsShared=.false.)
       call file%openFile(char(fileName))
       call file%writeDataset(halfMassRadii,'halfMassRadii')
       call file%close()
       call File_Unlock(lockFileDescriptor)
       call hdf5Access%unset()
    end if
    ! Reverse the table
    call self%inclinationTable%reverse(self%sizeTable)
    return
  end function diskSizeInclinationConstructorInternal
  
  double precision function diskSizeInclntnRoot(xHalf)
    !% Function used in solving for the half-light radii of inclined disks.
    use FGSL                    , only : fgsl_function, fgsl_integration_workspace
    use Numerical_Integration
    use Numerical_Constants_Math
    implicit none
    double precision                            , intent(in   ) :: xHalf
    double precision                                            :: integralHalf
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    logical                                                     :: integrationReset

    integrationReset=.true.
    integralHalf=Integrate(0.0d0,xHalf,diskSizeInclntnIntegrandX,integrandFunction,integrationWorkspace,toleranceRelative=1.0d-6,reset=integrationReset)
    diskSizeInclntnRoot=integralHalf/2.0d0/Pi/cos(diskSizeInclntnAngle)-0.5d0
    return
  end function diskSizeInclntnRoot

  double precision function diskSizeInclntnIntegrandX(x)
    !% Integral for half-light radius.
    use Numerical_Integration
    use Numerical_Constants_Math
    use FGSL                    , only : fgsl_function, fgsl_integration_workspace
    implicit none
    double precision                            , intent(in   ) :: x
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    logical                                                     :: integrationReset

    integrationReset         =.true.
    diskSizeInclntnXIntegrate=x
    diskSizeInclntnIntegrandX=Integrate(0.0d0,2.0d0*Pi,diskSizeInclntnIntegrandPhi,integrandFunction,integrationWorkspace,toleranceRelative=1.0d-6,reset=integrationReset)*x
    return
  end function diskSizeInclntnIntegrandX

  double precision function diskSizeInclntnIntegrandPhi(phi)
    !% Integral for half-light radius.
    implicit none
    double precision, intent(in   ) :: phi

    diskSizeInclntnIntegrandPhi=exp(-diskSizeInclntnXIntegrate*sqrt((sin(phi)/cos(diskSizeInclntnAngle))**2+cos(phi)**2))  
    return
  end function diskSizeInclntnIntegrandPhi

  function diskSizeInclinationOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !% Implement a disk size inclination output analysis distribution operator.
    implicit none
    class           (outputAnalysisDistributionOperatorDiskSizeInclntn), intent(inout)                                        :: self
    double precision                                                   , intent(in   )                                        :: propertyValue
    integer                                                            , intent(in   )                                        :: propertyType
    double precision                                                   , intent(in   ), dimension(:)                          :: propertyValueMinimum            , propertyValueMaximum
    integer         (c_size_t                                         ), intent(in   )                                        :: outputIndex
    type            (treeNode                                         ), intent(inout)                                        :: node
    double precision                                                                  , dimension(size(propertyValueMinimum)) :: diskSizeInclinationOperateScalar
    double precision                                                                                                          :: ratioLogarithmicMinimum         , ratioLogarithmicMaximum
    integer                                                                                                                   :: i
    !GCC$ attributes unused :: outputIndex, propertyType, node
    
    do i=1,size(propertyValueMinimum)
       ratioLogarithmicMinimum            =min(0.0d0,propertyValueMinimum(i)-propertyValue)
       ratioLogarithmicMaximum            =min(0.0d0,propertyValueMaximum(i)-propertyValue)
       diskSizeInclinationOperateScalar(i)=+self%sizeTable%interpolate(ratioLogarithmicMaximum) &
            &                              -self%sizeTable%interpolate(ratioLogarithmicMinimum)
    end do
    return
  end function diskSizeInclinationOperateScalar

  function diskSizeInclinationOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !% Implement a disk size inclination output analysis distribution operator.
    use Galacticus_Error
    implicit none
    class           (outputAnalysisDistributionOperatorDiskSizeInclntn), intent(inout)                                        :: self
    double precision                                                   , intent(in   ), dimension(:)                          :: distribution
    integer                                                            , intent(in   )                                        :: propertyType
    double precision                                                   , intent(in   ), dimension(:)                          :: propertyValueMinimum                  , propertyValueMaximum
    integer         (c_size_t                                         ), intent(in   )                                        :: outputIndex
    type            (treeNode                                         ), intent(inout)                                        :: node
    double precision                                                                  , dimension(size(propertyValueMinimum)) :: diskSizeInclinationOperateDistribution
    !GCC$ attributes unused :: self, distribution, propertyValueMinimum, propertyValueMaximum, outputIndex, propertyType, node

    diskSizeInclinationOperateDistribution=0.0d0
    call Galacticus_Error_Report('not implemented'//{introspection:location})
    return
  end function diskSizeInclinationOperateDistribution

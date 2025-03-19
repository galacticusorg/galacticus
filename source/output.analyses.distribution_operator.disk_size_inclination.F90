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
  Implements the effects of inclination on disk size in an output analysis distribution operator class.
  !!}

  use :: Tables, only : table1D, table1DLinearLinear

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorDiskSizeInclntn">
   <description>An output analysis distribution operator class which implements the effects of inclination on disk size.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorDiskSizeInclntn
     !!{
     An output distribution operator class which implements the effects of inclination on disk size.
     !!}
     private
     type (table1DLinearLinear)              :: inclinationTable
     class(table1D            ), allocatable :: sizeTable
   contains
     final     ::                        diskSizeInclinationDestructor
     procedure :: operateScalar       => diskSizeInclinationOperateScalar
     procedure :: operateDistribution => diskSizeInclinationOperateDistribution
  end type outputAnalysisDistributionOperatorDiskSizeInclntn

  interface outputAnalysisDistributionOperatorDiskSizeInclntn
     !!{
     Constructors for the ``diskSizeInclination'' output distribution operator class.
     !!}
     module procedure diskSizeInclinationConstructorParameters
     module procedure diskSizeInclinationConstructorInternal
  end interface outputAnalysisDistributionOperatorDiskSizeInclntn

  ! Module-scope variables used in root finding and integrations.
  double precision :: xIntegrate, angle
  !$omp threadprivate(xIntegrate,angle)

contains

  function diskSizeInclinationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``diskSizeInclination'' output analysis distribution operator operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisDistributionOperatorDiskSizeInclntn)                :: self
    type(inputParameters                                  ), intent(inout) :: parameters

    self=outputAnalysisDistributionOperatorDiskSizeInclntn()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function diskSizeInclinationConstructorParameters

  function diskSizeInclinationConstructorInternal() result(self)
    !!{
    Internal constructor for the ``diskSizeInclination'' output analysis distribution operator class.
    !!}
    use :: File_Utilities    , only : Directory_Make           , File_Exists                  , File_Lock                    , File_Unlock, &
          &                           lockDescriptor
    use :: Input_Paths       , only : inputPath                , pathTypeDataDynamic
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : varying_string
    use :: Root_Finder       , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: Table_Labels      , only : extrapolationTypeFix
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
    fileName=inputPath(pathTypeDataDynamic)//"/galacticStructure/diskExponentialInclinedHalfMassRadii.hdf5"
    if (File_Exists(fileName)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),lockFileDescriptor,lockIsShared=.true.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName),readOnly=.true.)
       call file%readDataset('halfMassRadii',halfMassRadii)
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(lockFileDescriptor)
       call self%inclinationTable%populate(halfMassRadii)
    else
       ! Tabulate dependence of projected half-light radius on disk inclination angle.
       finder=rootFinder(                                                             &
            &            rootFunction                 =diskSizeInclntnRoot          , &
            &            toleranceRelative            =1.0d-6                       , &
            &            toleranceAbsolute            =1.0d-6                       , &
            &            rangeExpandUpward            =2.0d0                        , &
            &            rangeExpandDownward          =0.5d0                        , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &            rangeExpandType              =rangeExpandMultiplicative      &
            &           )
       angle =acos(self%inclinationTable%x(1))
       halfLightRadiusFaceOn=finder%find(rootGuess=1.0d0)
       call self%inclinationTable%populate(0.0d0,1)
       !$omp parallel do private (i,halfLightRadius) copyin (finder)
       do i=2,inclinationAngleCount
          angle=acos(self%inclinationTable%x(i))
          halfLightRadius =finder%find(rootGuess=1.0d0)
          call self%inclinationTable%populate(log10(halfLightRadius/halfLightRadiusFaceOn),i)
       end do
       !$omp end parallel do
       halfMassRadii=reshape(self%inclinationTable%ys(),[inclinationAngleCount])
       call Directory_Make(inputPath(pathTypeDataDynamic)//"/galacticStructure")
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),lockFileDescriptor,lockIsShared=.false.)
       !$ call hdf5Access%set()
       call file%openFile(char(fileName),objectsOverwritable=.true.)
       call file%writeDataset(halfMassRadii,'halfMassRadii')
       call file%close()
       !$ call hdf5Access%unset()
       call File_Unlock(lockFileDescriptor)
    end if
    ! Reverse the table
    call self%inclinationTable%reverse(self%sizeTable)
    return
  end function diskSizeInclinationConstructorInternal

  subroutine diskSizeInclinationDestructor(self)
    !!{
    Destructor for the ``diskSizeInclination'' output analysis distribution operator operator class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorDiskSizeInclntn), intent(inout) :: self

    call self%inclinationTable%destroy()
    if (allocated(self%sizeTable)) then
       call self%sizeTable%destroy()
       deallocate(self%sizeTable)
    end if
    return
  end subroutine diskSizeInclinationDestructor

  double precision function diskSizeInclntnRoot(xHalf)
    !!{
    Function used in solving for the half-light radii of inclined disks.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    double precision            , intent(in   ) :: xHalf
    double precision                            :: integralHalf
    type            (integrator)                :: integrator_

    integrator_        =integrator           (diskSizeInclntnIntegrandX,toleranceRelative=1.0d-6)
    integralHalf       =integrator_%integrate(0.0d0                    ,                  xHalf )
    diskSizeInclntnRoot=integralHalf/2.0d0/Pi/cos(angle)-0.5d0
    return
  end function diskSizeInclntnRoot

  double precision function diskSizeInclntnIntegrandX(x)
    !!{
    Integral for half-light radius.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    implicit none
    double precision            , intent(in   ) :: x
    type            (integrator)                :: integrator_

    xIntegrate               =+x
    integrator_              = integrator           (diskSizeInclntnIntegrandPhi,toleranceRelative=1.0d-6   )
    diskSizeInclntnIntegrandX=+integrator_%integrate(0.0d0                      ,                  2.0d+0*Pi) &
         &                    *x
    return
  end function diskSizeInclntnIntegrandX

  double precision function diskSizeInclntnIntegrandPhi(phi)
    !!{
    Integral for half-light radius.
    !!}
    implicit none
    double precision, intent(in   ) :: phi

    diskSizeInclntnIntegrandPhi=exp(-xIntegrate*sqrt((sin(phi)/cos(angle))**2+cos(phi)**2))
    return
  end function diskSizeInclntnIntegrandPhi

  function diskSizeInclinationOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a disk size inclination output analysis distribution operator.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorDiskSizeInclntn), intent(inout)                                        :: self
    double precision                                                   , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType        ), intent(in   )                                        :: propertyType
    double precision                                                   , intent(in   ), dimension(:)                          :: propertyValueMinimum            , propertyValueMaximum
    integer         (c_size_t                                         ), intent(in   )                                        :: outputIndex
    type            (treeNode                                         ), intent(inout)                                        :: node
    double precision                                                                  , dimension(size(propertyValueMinimum)) :: diskSizeInclinationOperateScalar
    double precision                                                                                                          :: ratioLogarithmicMinimum         , ratioLogarithmicMaximum
    integer                                                                                                                   :: i
    !$GLC attributes unused :: outputIndex, propertyType, node

    do i=1,size(propertyValueMinimum)
       ratioLogarithmicMinimum            =min(0.0d0,propertyValueMinimum(i)-propertyValue)
       ratioLogarithmicMaximum            =min(0.0d0,propertyValueMaximum(i)-propertyValue)
       diskSizeInclinationOperateScalar(i)=+self%sizeTable%interpolate(ratioLogarithmicMaximum) &
            &                              -self%sizeTable%interpolate(ratioLogarithmicMinimum)
    end do
    return
  end function diskSizeInclinationOperateScalar

  function diskSizeInclinationOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement a disk size inclination output analysis distribution operator.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisDistributionOperatorDiskSizeInclntn), intent(inout)                                        :: self
    double precision                                                   , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType        ), intent(in   )                                        :: propertyType
    double precision                                                   , intent(in   ), dimension(:)                          :: propertyValueMinimum                  , propertyValueMaximum
    integer         (c_size_t                                         ), intent(in   )                                        :: outputIndex
    type            (treeNode                                         ), intent(inout)                                        :: node
    double precision                                                                  , dimension(size(propertyValueMinimum)) :: diskSizeInclinationOperateDistribution
    !$GLC attributes unused :: self, distribution, propertyValueMinimum, propertyValueMaximum, outputIndex, propertyType, node

    diskSizeInclinationOperateDistribution=0.0d0
    call Error_Report('not implemented'//{introspection:location})
    return
  end function diskSizeInclinationOperateDistribution

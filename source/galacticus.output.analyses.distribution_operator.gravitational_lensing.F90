!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a gravitational lensing output analysis distribution operator class. 

  !$use OMP_Lib
  use Gravitational_Lensing
  use Locks
  
  type :: grvtnlLnsngTransferMatrix
     double precision, allocatable, dimension(:,:) :: matrix
  end type grvtnlLnsngTransferMatrix

  !# <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorGrvtnlLnsng" defaultThreadPrivate="yes">
  !#  <description>A gravitational lensing output analysis distribution operator class.</description>
  !# </outputAnalysisDistributionOperator>
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorGrvtnlLnsng
     !% A gravitational lensing output distribution operator class.
     private
     class           (gravitationalLensingClass), pointer                   :: gravitationalLensing_
     type            (grvtnlLnsngTransferMatrix), allocatable, dimension(:) :: transfer_
     double precision                                                       :: sizeSource
     !$ type         (ompReadWriteLock         )                            :: tabulateLock
   contains
     final     ::                        grvtnlLnsngDestructor
     procedure :: operateScalar       => grvtnlLnsngOperateScalar
     procedure :: operateDistribution => grvtnlLnsngOperateDistribution
  end type outputAnalysisDistributionOperatorGrvtnlLnsng

  interface outputAnalysisDistributionOperatorGrvtnlLnsng
     !% Constructors for the ``gravitational lensing'' output analysis distribution operator class.
     module procedure grvtnlLnsngConstructorParameters
     module procedure grvtnlLnsngConstructorInternal
  end interface outputAnalysisDistributionOperatorGrvtnlLnsng

contains

  function grvtnlLnsngConstructorParameters(parameters) result(self)
    !% Constructor for the ``gravitational lensing'' output analysis distribution operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (gravitationalLensingClass                    ), pointer       :: gravitationalLensing_
    double precision                                                               :: sizeSource
    
    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>sizeSource</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.001d0</defaultValue>
    !#   <description>The source size to assume for gravitational lensing calculations.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    ! Construct the object.
    self=outputAnalysisDistributionOperatorGrvtnlLnsng(gravitationalLensing_,sizeSource)
    !# <inputParametersValidate source="parameters"/>
    return
  end function grvtnlLnsngConstructorParameters

  function grvtnlLnsngConstructorInternal(gravitationalLensing_,sizeSource) result(self)
    !% Internal constructor for the ``gravitational lensing'' output analysis distribution operator class.
    use Galacticus_Output_Times
    implicit none
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng)                        :: self
    class           (gravitationalLensingClass                    ), intent(in   ), target :: gravitationalLensing_
    double precision                                               , intent(in   )         :: sizeSource
    !# <constructorAssign variables="*gravitationalLensing_, sizeSource"/>

    ! Allocate transfer matrices for all outputs.
    allocate(self%transfer_(Galacticus_Output_Time_Count()))
    !$ self%tabulateLock=ompReadWriteLock()
    return
  end function grvtnlLnsngConstructorInternal

  subroutine grvtnlLnsngDestructor(self)
    !% Destructor for the ``gravitational lensing'' output analysis distribution operator class.
    implicit none
    type(outputAnalysisDistributionOperatorGrvtnlLnsng), intent(inout) :: self

    !# <objectDestructor name="self%gravitationalLensing_"/>
    return
  end subroutine grvtnlLnsngDestructor
  
  function grvtnlLnsngOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node) result(distributionNew)
    !% Implement a gravitational lensing output analysis distribution operator.
    use Galacticus_Error
    implicit none
    class           (outputAnalysisDistributionOperatorGrvtnlLnsng), intent(inout)                                        :: self
    double precision                                               , intent(in   )                                        :: propertyValue
    integer                                                        , intent(in   )                                        :: propertyType
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum, propertyValueMaximum
    integer         (c_size_t                                     ), intent(in   )                                        :: outputIndex
    type            (treeNode                                     ), intent(inout)                                        :: node
    double precision                                                              , dimension(size(propertyValueMinimum)) :: distributionNew
    !GCC$ attributes unused :: self, propertyValue, propertyType, propertyValueMinimum, propertyValueMaximum, outputIndex, node

    distributionNew=0.0d0
    call Galacticus_Error_Report('grvtnlLnsngOperateScalar','not implemented')
    return
  end function grvtnlLnsngOperateScalar
  
  function grvtnlLnsngOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node) result(distributionNew)
    !% Implement a gravitational lensing output analysis distribution operator.
    use Memory_Management
    use FGSL
    use Numerical_Integration
    use Galacticus_Output_Times
    use Output_Analyses_Options
    implicit none
    class           (outputAnalysisDistributionOperatorGrvtnlLnsng), intent(inout)                                        :: self
    double precision                                               , intent(in   ), dimension(:)                          :: distribution
    integer                                                        , intent(in   )                                        :: propertyType
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum, propertyValueMaximum
    integer         (c_size_t                                     ), intent(in   )                                        :: outputIndex
    type            (treeNode                                     ), intent(inout)                                        :: node
    double precision                                                              , dimension(size(propertyValueMinimum)) :: distributionNew
    double precision                                                                                                      :: redshift
    integer         (c_size_t                                     )                                                       :: j                   , k
    type            (fgsl_function                                )                                                       :: integrandFunction
    type            (fgsl_integration_workspace                   )                                                       :: integrationWorkspace
    logical                                                                                                               :: integrationReset
    !GCC$ attributes unused :: node

    ! Construct the lensing transfer matrix if not already done.
    !$ call self%tabulateLock%setRead()
    if (.not.allocated(self%transfer_(outputIndex)%matrix)) then
       !$ call self%tabulateLock%setWrite(haveReadLock=.true.)
       if (.not.allocated(self%transfer_(outputIndex)%matrix)) then
          call allocateArray(self%transfer_(outputIndex)%matrix,[size(propertyValueMinimum),size(propertyValueMinimum)])
          redshift=Galacticus_Output_Redshift(outputIndex)
          do j=1,size(propertyValueMinimum)
             do k=1,size(propertyValueMinimum)
                if (j > 1 .and. k > 1) then
                   ! Transfer matrix elements are identical along diagonals of the matrix.
                   self       %transfer_(outputIndex)%matrix(j  ,k  )=                &
                        & self%transfer_(outputIndex)%matrix(j-1,k-1)
                else
                   integrationReset=.true.                               
                   self       %transfer_(outputIndex)%matrix(j  ,k  )=                &
                        & +Integrate(                                                 &
                        &                               propertyValueMinimum     (j), &
                        &                               propertyValueMaximum     (j), &
                        &                               magnificationCDFIntegrand   , &
                        &                               integrandFunction           , &
                        &                               integrationWorkspace        , &
                        &             toleranceRelative=1.0d-3                      , &
                        &             reset            =integrationReset              &
                        &            )                                                &
                        & /(                                                          &
                        &   +propertyValueMaximum(j)                                  &
                        &   -propertyValueMinimum(j)                                  &
                        &  )
                   call Integrate_Done(integrandFunction,integrationWorkspace)
                end if
             end do
          end do
       end if
       !$ call self%tabulateLock%unsetWrite(haveReadLock=.true.)
    end if
    ! Apply the lensing transfer matrix.
    distributionNew=matmul(distribution,self%transfer_(outputIndex)%matrix)
    !$ call self%tabulateLock%unsetRead()
    return

  contains

    double precision function magnificationCDFIntegrand(propertyValue)
      !% Integrand over the gravitational lensing magnification cumulative distribution.
      use Output_Analyses_Options
      use Galacticus_Error
      implicit none
      double precision, intent(in   ) :: propertyValue
      double precision                :: ratioMinimum , ratioMaximum

      ! Find the minimum and maximum magnification ratios.
      select case (propertyType)
      case (outputAnalysisPropertyTypeLinear   )
         ratioMinimum=                 propertyValueMaximum(k)/propertyValue
         ratioMaximum=                 propertyValueMinimum(k)/propertyValue
      case (outputAnalysisPropertyTypeLog10    )
         ratioMinimum=10.0d0**(        propertyValueMinimum(k)-propertyValue)
         ratioMaximum=10.0d0**(        propertyValueMaximum(k)-propertyValue)
      case (outputAnalysisPropertyTypeMagnitude)
         ratioMinimum=10.0d0**(-0.4d0*(propertyValueMinimum(k)-propertyValue))
         ratioMaximum=10.0d0**(-0.4d0*(propertyValueMaximum(k)-propertyValue))
      case default
         call Galacticus_Error_Report('magnificationCDFIntegrand','unknown property type')
      end select
      ! Compute the lensing CDF across this range of magnifications.
      magnificationCDFIntegrand=                                               &
           & max(                                                              &
           &     +0.0d0                                                      , &
           &     +self%gravitationalLensing_%magnificationCDF(                 &
           &                                                  ratioMaximum   , &
           &                                                  redshift       , &
           &                                                  self%sizeSource  &
           &                                                 )                 &
           &     -self%gravitationalLensing_%magnificationCDF(                 &
           &                                                  ratioMinimum   , &
           &                                                  redshift       , &
           &                                                  self%sizeSource  &
           &                                                 )                 &
           &    )
      return
    end function magnificationCDFIntegrand

  end function grvtnlLnsngOperateDistribution

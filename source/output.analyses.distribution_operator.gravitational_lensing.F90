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
Implements a gravitational lensing output analysis distribution operator class.
!!}

  use :: Gravitational_Lensing, only : gravitationalLensingClass
  use :: Locks                , only : ompReadWriteLock
  use :: Output_Times         , only : outputTimesClass

  type :: grvtnlLnsngTransferMatrix
     double precision, allocatable, dimension(:,:) :: matrix
  end type grvtnlLnsngTransferMatrix

  !![
  <enumeration>
   <name>lensedProperty</name>
   <description>Enumeration of properties affected by gravitational lensing.</description>
   <visibility>public</visibility>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="luminosity"/>
   <entry label="size"      />
  </enumeration>
  !!]

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorGrvtnlLnsng">
   <description>A gravitational lensing output analysis distribution operator class.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorGrvtnlLnsng
     !!{
     A gravitational lensing output distribution operator class.
     !!}
     private
     class           (gravitationalLensingClass    ), pointer                   :: gravitationalLensing_ => null()
     class           (outputTimesClass             ), pointer                   :: outputTimes_          => null()
     type            (grvtnlLnsngTransferMatrix    ), allocatable, dimension(:) :: transfer_
     type            (enumerationLensedPropertyType)                            :: lensedProperty
     double precision                                                           :: sizeSource
     !$ type         (ompReadWriteLock             ), allocatable, dimension(:) :: tabulateLock
   contains
     final     ::                        grvtnlLnsngDestructor
     procedure :: operateScalar       => grvtnlLnsngOperateScalar
     procedure :: operateDistribution => grvtnlLnsngOperateDistribution
  end type outputAnalysisDistributionOperatorGrvtnlLnsng

  interface outputAnalysisDistributionOperatorGrvtnlLnsng
     !!{
     Constructors for the \refClass{outputAnalysisDistributionOperatorGrvtnlLnsng} output analysis distribution operator class.
     !!}
     module procedure grvtnlLnsngConstructorParameters
     module procedure grvtnlLnsngConstructorInternal
  end interface outputAnalysisDistributionOperatorGrvtnlLnsng

  ! Module scope lensing object used in parallel evaluation of the lensing matrix.
  integer(c_size_t                 )          :: k_
  class  (gravitationalLensingClass), pointer :: gravitationalLensing_
  !$omp threadprivate(k_,gravitationalLensing_)

contains

  function grvtnlLnsngConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionOperatorGrvtnlLnsng} output analysis distribution operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (gravitationalLensingClass                    ), pointer       :: gravitationalLensing_
    class           (outputTimesClass                             ), pointer       :: outputTimes_
    type            (varying_string                               )                :: lensedProperty
    double precision                                                               :: sizeSource

    !![
    <inputParameter>
      <name>lensedProperty</name>
      <source>parameters</source>
      <defaultValue>var_str('luminosity')</defaultValue>
      <description>The property (luminosity, or size) to be affected by gravitational lensing.</description>
    </inputParameter>
    <inputParameter>
      <name>sizeSource</name>
      <source>parameters</source>
      <defaultValue>0.001d0</defaultValue>
      <description>The source size to assume for gravitational lensing calculations.</description>
    </inputParameter>
    <objectBuilder class="gravitationalLensing" name="gravitationalLensing_" source="parameters"/>
    <objectBuilder class="outputTimes"          name="outputTimes_"          source="parameters"/>
    !!]
    ! Construct the object.
    self=outputAnalysisDistributionOperatorGrvtnlLnsng(gravitationalLensing_,outputTimes_,sizeSource,enumerationLensedPropertyEncode(char(lensedProperty),includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="gravitationalLensing_"/>
    <objectDestructor name="outputTimes_"         />
    !!]
    return
  end function grvtnlLnsngConstructorParameters

  function grvtnlLnsngConstructorInternal(gravitationalLensing_,outputTimes_,sizeSource,lensedProperty) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisDistributionOperatorGrvtnlLnsng} output analysis distribution operator class.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Error        , only : Error_Report
    implicit none
    type            (outputAnalysisDistributionOperatorGrvtnlLnsng)                          :: self
    class           (gravitationalLensingClass                    ), intent(in   ), target   :: gravitationalLensing_
    class           (outputTimesClass                             ), intent(in   ), target   :: outputTimes_
    type            (enumerationLensedPropertyType                ), intent(in   ), optional :: lensedProperty
    double precision                                               , intent(in   )           :: sizeSource
    !$ integer      (c_size_t                                     )                          :: i
    !![
    <optionalArgument name="lensedProperty" defaultsTo="lensedPropertyLuminosity" />
    <constructorAssign variables="*gravitationalLensing_, *outputTimes_, sizeSource"/>
    !!]
    
    if (.not.enumerationLensedPropertyIsValid(lensedProperty_)) call Error_Report('invalid lensedProperty'//{introspection:location})
    self%lensedProperty=lensedProperty_
    ! Allocate transfer matrices for all outputs.
    allocate   (self%transfer_   (self%outputTimes_%count()))
    !$ allocate(self%tabulateLock(self%outputTimes_%count()))
    !$ do i=1,self%outputTimes_%count()
    !$    self%tabulateLock(i)=ompReadWriteLock()
    !$ end do
    return
  end function grvtnlLnsngConstructorInternal

  subroutine grvtnlLnsngDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisDistributionOperatorGrvtnlLnsng} output analysis distribution operator class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorGrvtnlLnsng), intent(inout) :: self

    !![
    <objectDestructor name="self%gravitationalLensing_"/>
    <objectDestructor name="self%outputTimes_"         />
    !!]
    return
  end subroutine grvtnlLnsngDestructor

  function grvtnlLnsngOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node) result(distributionNew)
    !!{
    Implement a gravitational lensing output analysis distribution operator.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisDistributionOperatorGrvtnlLnsng), intent(inout)                                        :: self
    double precision                                               , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   )                                        :: propertyType
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum, propertyValueMaximum
    integer         (c_size_t                                     ), intent(in   )                                        :: outputIndex
    type            (treeNode                                     ), intent(inout)                                        :: node
    double precision                                                              , dimension(size(propertyValueMinimum)) :: distributionNew
    !$GLC attributes unused :: self, propertyValue, propertyType, propertyValueMinimum, propertyValueMaximum, outputIndex, node

    distributionNew=0.0d0
    call Error_Report('not implemented'//{introspection:location})
    return
  end function grvtnlLnsngOperateScalar

  function grvtnlLnsngOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node) result(distributionNew)
    !!{
    Implement a gravitational lensing output analysis distribution operator.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (outputAnalysisDistributionOperatorGrvtnlLnsng), intent(inout)                                        :: self
    double precision                                               , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType    ), intent(in   )                                        :: propertyType
    double precision                                               , intent(in   ), dimension(:)                          :: propertyValueMinimum , propertyValueMaximum
    integer         (c_size_t                                     ), intent(in   )                                        :: outputIndex
    type            (treeNode                                     ), intent(inout)                                        :: node
    double precision                                                              , dimension(size(propertyValueMinimum)) :: distributionNew
    double precision                                                                                                      :: redshift
    integer         (c_size_t                                     )                                                       :: j                    , k                   , &
         &                                                                                                                   i                    , l
    type            (integrator                                   )               , allocatable                           :: integrator_
    !$GLC attributes unused :: node

    ! Construct the lensing transfer matrix if not already done.
    !$ call self%tabulateLock(outputIndex)%setRead()
    if (.not.allocated(self%transfer_(outputIndex)%matrix)) then
       !$ call self%tabulateLock(outputIndex)%setWrite(haveReadLock=.true.)
       if (.not.allocated(self%transfer_(outputIndex)%matrix)) then
          redshift=self%outputTimes_%redshift(outputIndex)
          allocate(self%transfer_(outputIndex)%matrix(size(propertyValueMinimum),size(propertyValueMinimum)))
          !$omp parallel private (i,j,k,l,integrator_)
          allocate(gravitationalLensing_,mold=self%gravitationalLensing_)
          !$omp critical(analysesGravitationalLensingDeepCopy)
          !![
          <deepCopyReset variables="self%gravitationalLensing_"/>
          <deepCopy source="self%gravitationalLensing_" destination="gravitationalLensing_"/>
          <deepCopyFinalize variables="gravitationalLensing_"/>
          !!]
          !$omp end critical(analysesGravitationalLensingDeepCopy)
          allocate(integrator_)
          integrator_=integrator(magnificationCDFIntegrand,toleranceRelative=1.0d-3)
          !$omp do schedule(dynamic)
          do l=1,size(propertyValueMinimum)
             do i=1,2
                if (i == 1) then
                   j=l; k_=1
                else
                   j=1; k_=l
                end if
                self%transfer_(outputIndex)%matrix(j,k_)=+integrator_%integrate( propertyValueMinimum(j),propertyValueMaximum(j)) &
                     &                                   /                     (+propertyValueMaximum(j)-propertyValueMinimum(j))
             end do
          end do
          !$omp end do
          !![
          <objectDestructor name="gravitationalLensing_"/>
          !!]
          deallocate(integrator_)
          !$omp end parallel
          do j=2,size(propertyValueMinimum)
             do k=2,size(propertyValueMinimum)
                ! Transfer matrix elements are identical along diagonals of the matrix.
                self       %transfer_(outputIndex)%matrix(j  ,k  )= &
                     & self%transfer_(outputIndex)%matrix(j-1,k-1)
             end do
          end do
       end if
       !$ call self%tabulateLock(outputIndex)%unsetWrite(haveReadLock=.true.)
    end if
    ! Apply the lensing transfer matrix.
    distributionNew=matmul(distribution,self%transfer_(outputIndex)%matrix)
    !$ call self%tabulateLock(outputIndex)%unsetRead()
    return

  contains

    double precision function magnificationCDFIntegrand(propertyValue)
      !!{
      Integrand over the gravitational lensing magnification cumulative distribution.
      !!}
      use :: Error                  , only : Error_Report
      use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear, outputAnalysisPropertyTypeLog10, outputAnalysisPropertyTypeMagnitude
      implicit none
      double precision, intent(in   ) :: propertyValue
      double precision                :: ratioMinimum , ratioMaximum

      ! Find the minimum and maximum magnification ratios.
      select case (propertyType%ID)
      case (outputAnalysisPropertyTypeLinear   %ID)
         ratioMinimum=                 propertyValueMaximum(k_)/propertyValue
         ratioMaximum=                 propertyValueMinimum(k_)/propertyValue
      case (outputAnalysisPropertyTypeLog10    %ID)
         ratioMinimum=10.0d0**(        propertyValueMinimum(k_)-propertyValue)
         ratioMaximum=10.0d0**(        propertyValueMaximum(k_)-propertyValue)
      case (outputAnalysisPropertyTypeMagnitude%ID)
         ! Note that ratio min/max is related to property value max/min because magnitudes are brighter when more negative.
         ratioMinimum=10.0d0**(-0.4d0*(propertyValueMaximum(k_)-propertyValue))
         ratioMaximum=10.0d0**(-0.4d0*(propertyValueMinimum(k_)-propertyValue))
      case default
         call Error_Report('unknown property type'//{introspection:location})
      end select

      ! Scale ratios depending on the property being lensed.
      select case (self%lensedProperty%ID)
      case (lensedPropertyLuminosity%ID)
         ! Luminosity scales linearly with magnification - no need to modify ratios.
      case (lensedPropertySize      %ID)
         ! Size scales as the square-root of magnification. Our ratios as computed so far are ratios of size. To find the
         ! corresponding ratios of magnification we must therefore square the ratios.
         ratioMinimum=ratioMinimum**2
         ratioMaximum=ratioMaximum**2
      case default
         call Error_Report('unknown lensedProperty'//{introspection:location})
      end select
      ! Compute the lensing CDF across this range of magnifications.
      magnificationCDFIntegrand=                                          &
           & max(                                                         &
           &     +0.0d0                                                 , &
           &     +gravitationalLensing_%magnificationCDF(                 &
           &                                             ratioMaximum   , &
           &                                             redshift       , &
           &                                             self%sizeSource  &
           &                                            )                 &
           &     -gravitationalLensing_%magnificationCDF(                 &
           &                                             ratioMinimum   , &
           &                                             redshift       , &
           &                                             self%sizeSource  &
           &                                            )                 &
           &    )
      return
    end function magnificationCDFIntegrand

  end function grvtnlLnsngOperateDistribution

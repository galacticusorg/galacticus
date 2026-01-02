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
  Implements an output analysis distribution operator class to account for errors on N-body measurements of halo spin.
  !!}

  use :: Halo_Spin_Distributions, only : haloSpinDistributionClass

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorSpinNBodyErrors">
   <description>An output analysis distribution operator class to account for errors on N-body measurements of halo spin.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorSpinNBodyErrors
     !!{
     An output distribution operator class to account for errors on N-body measurements of halo spin.
     !!}
     private
     class  (haloSpinDistributionClass), pointer :: haloSpinDistribution_ => null()
     logical                                     :: errorTolerant
   contains
     final     ::                        spinNBodyErrorsDestructor
     procedure :: operateScalar       => spinNBodyErrorsOperateScalar
     procedure :: operateDistribution => spinNBodyErrorsOperateDistribution
  end type outputAnalysisDistributionOperatorSpinNBodyErrors

  interface outputAnalysisDistributionOperatorSpinNBodyErrors
     !!{
     Constructors for the \refClass{outputAnalysisDistributionOperatorSpinNBodyErrors} output analysis distribution operator class.
     !!}
     module procedure spinNBodyErrorsConstructorParameters
     module procedure spinNBodyErrorsConstructorInternal
  end interface outputAnalysisDistributionOperatorSpinNBodyErrors

contains

  function spinNBodyErrorsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionOperatorSpinNBodyErrors} output analysis distribution operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (outputAnalysisDistributionOperatorSpinNBodyErrors)                :: self
    type   (inputParameters                                  ), intent(inout) :: parameters
    logical                                                                   :: errorTolerant
    class  (haloSpinDistributionClass                        ), pointer       :: haloSpinDistribution_

    !![
    <inputParameter>
      <name>errorTolerant</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, integration tolerance failures are tolerated (a warning is issued but calculations will continue).</description>
    </inputParameter>
    <objectBuilder class="haloSpinDistribution" name="haloSpinDistribution_" source="parameters"/>
    !!]
    self=outputAnalysisDistributionOperatorSpinNBodyErrors(errorTolerant,haloSpinDistribution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="haloSpinDistribution_"/>
    !!]
    return
  end function spinNBodyErrorsConstructorParameters

  function spinNBodyErrorsConstructorInternal(errorTolerant,haloSpinDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisDistributionOperatorSpinNBodyErrors} output analysis distribution operator class.
    !!}
    use :: Error                  , only : Error_Report
    use :: Halo_Spin_Distributions, only : haloSpinDistributionClass, haloSpinDistributionNbodyErrors
    implicit none
    type   (outputAnalysisDistributionOperatorSpinNBodyErrors)                        :: self
    logical                                                   , intent(in   )         :: errorTolerant
    class  (haloSpinDistributionClass                        ), intent(in   ), target :: haloSpinDistribution_
    !![
    <constructorAssign variables="errorTolerant, *haloSpinDistribution_"/>
    !!]

    ! Check that the spin distribution is of the N-body errors class.
    select type (haloSpinDistribution_)
    class is (haloSpinDistributionNbodyErrors)
       ! This is OK, do nothing.
    class default
       ! This is not OK.
       call Error_Report('must provide an N-body errors halo spin distribution class'//{introspection:location})
    end select
    return
  end function spinNBodyErrorsConstructorInternal

  subroutine spinNBodyErrorsDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisDistributionOperatorSpinNBodyErrors} output analysis distribution operator class.
    !!}
    implicit none
    type (outputAnalysisDistributionOperatorSpinNBodyErrors), intent(inout) :: self

    !![
    <objectDestructor name="self%haloSpinDistribution_" />
    !!]
    return
  end subroutine spinNBodyErrorsDestructor

  function spinNBodyErrorsOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement an output analysis distribution operator which accounts for errors in N-body measurements of halo spin.
    !!}
    use :: Error                  , only : Error_Report                    , Warn                           , errorStatusSuccess
    use :: Numerical_Integration  , only : integrator
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear, outputAnalysisPropertyTypeLog10
    implicit none
    class           (outputAnalysisDistributionOperatorSpinNBodyErrors), intent(inout)                                        :: self
    double precision                                                   , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType        ), intent(in   )                                        :: propertyType
    double precision                                                   , intent(in   ), dimension(:)                          :: propertyValueMinimum        , propertyValueMaximum
    integer         (c_size_t                                         ), intent(in   )                                        :: outputIndex
    type            (treeNode                                         ), intent(inout)                                        :: node
    double precision                                                                  , dimension(size(propertyValueMinimum)) :: spinNBodyErrorsOperateScalar
    integer         (c_size_t                                         )                                                       :: i
    integer                                                                                                                   :: status
    double precision                                                                                                          :: spinMeasuredMinimum         , spinMeasuredMaximum     , &
         &                                                                                                                       spinMeasuredRangeMinimum    , spinMeasuredRangeMaximum
    type            (integrator                                       )                                                       :: integrator_
    !$GLC attributes unused :: outputIndex, propertyValue

    select case (propertyType%ID)
    case (outputAnalysisPropertyTypeLinear%ID)
       spinMeasuredRangeMinimum=        propertyValueMinimum(                        1 )
       spinMeasuredRangeMaximum=        propertyValueMaximum(size(propertyValueMaximum))
    case (outputAnalysisPropertyTypeLog10%ID)
       spinMeasuredRangeMinimum=10.0d0**propertyValueMinimum(                        1 )
       spinMeasuredRangeMaximum=10.0d0**propertyValueMaximum(size(propertyValueMaximum))
    case default
       call Error_Report('unhandled property type'//{introspection:location})
    end select
    integrator_=integrator(spinDistributionIntegrate,toleranceRelative=1.0d-6)
    do i=1,size(propertyValueMinimum)
       select case (propertyType%ID)
       case (outputAnalysisPropertyTypeLinear%ID)
          spinMeasuredMinimum=        propertyValueMinimum(i)
          spinMeasuredMaximum=        propertyValueMaximum(i)
       case (outputAnalysisPropertyTypeLog10%ID)
          spinMeasuredMinimum=10.0d0**propertyValueMinimum(i)
          spinMeasuredMaximum=10.0d0**propertyValueMaximum(i)
       case default
          call Error_Report('unhandled property type'//{introspection:location})
       end select
       spinNBodyErrorsOperateScalar(i)=+integrator_%integrate( spinMeasuredMinimum,spinMeasuredMaximum,status=status) &
            &                          /                     (+spinMeasuredMaximum-spinMeasuredMinimum              )
       if (status /= errorStatusSuccess) then
          if (self%errorTolerant) then
             call Warn        ('integration of N-body spin distribution failed'                          )
          else
             call Error_Report('integration of N-body spin distribution failed'//{introspection:location})
          end if
       end if
    end do
    return

  contains

    double precision function spinDistributionIntegrate(spinMeasured)
      !!{
      Integrand function used to find cumulative spin distribution over a bin.
      !!}
      use :: Halo_Spin_Distributions, only : haloSpinDistributionNbodyErrors
     implicit none
      double precision, intent(in   ) :: spinMeasured

      select type (haloSpinDistribution_ => self%haloSpinDistribution_)
      class is (haloSpinDistributionNbodyErrors)
         spinDistributionIntegrate=+                                                  spinMeasured                                                    &
              &                    *haloSpinDistribution_%distributionFixedPoint(node,spinMeasured,spinMeasuredRangeMinimum,spinMeasuredRangeMaximum)
      class default
         spinDistributionIntegrate=0.0d0
      end select
      return
    end function spinDistributionIntegrate

  end function spinNBodyErrorsOperateScalar

  function spinNBodyErrorsOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement an output analysis distribution operator which accounts for errors in N-body measurements of halo spin.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisDistributionOperatorSpinNBodyErrors), intent(inout)                                        :: self
    double precision                                                   , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType        ), intent(in   )                                        :: propertyType
    double precision                                                   , intent(in   ), dimension(:)                          :: propertyValueMinimum              , propertyValueMaximum
    integer         (c_size_t                                         ), intent(in   )                                        :: outputIndex
    type            (treeNode                                         ), intent(inout)                                        :: node
    double precision                                                                  , dimension(size(propertyValueMinimum)) :: spinNBodyErrorsOperateDistribution
    !$GLC attributes unused :: self, distribution, propertyValueMinimum, propertyValueMaximum, outputIndex, propertyType, node

    spinNBodyErrorsOperateDistribution=0.0d0
    call Error_Report('not implemented'//{introspection:location})
    return
  end function spinNBodyErrorsOperateDistribution

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
  Implements an N-body halo mass ratio output analysis distribution operator class.
  !!}
  
  use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass
  use :: Galactic_Filters                 , only : galacticFilterClass

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorMassRatioNBody">
   <description>A random error output analysis distribution operator class.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorMassRatioNBody
     !!{
     An N-body halo mass ratio output analysis distribution operator class.
     !!}
     private
     class           (nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_ => null()
     class           (galacticFilterClass    ), pointer :: galacticFilter_     => null()
     double precision                                   :: massParentMinimum            , massParentMaximum, &
          &                                                timeParent     
   contains
     final     ::                        massRatioNBodyDestructor
     procedure :: operateScalar       => massRatioNBodyOperateScalar
     procedure :: operateDistribution => massRatioNBodyOperateDistribution
  end type outputAnalysisDistributionOperatorMassRatioNBody

  interface outputAnalysisDistributionOperatorMassRatioNBody
     module procedure massRatioNBodyConstructorParameters
     module procedure massRatioNBodyConstructorInternal
  end interface outputAnalysisDistributionOperatorMassRatioNBody
  
contains

  function massRatioNBodyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionOperatorMassRatioNBody} output analysis weight operator class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (outputAnalysisDistributionOperatorMassRatioNBody)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (nbodyHaloMassErrorClass                         ), pointer       :: nbodyHaloMassError_
    class           (galacticFilterClass                             ), pointer       :: galacticFilter_
    double precision                                                                  :: massParentMinimum  , massParentMaximum, &
         &                                                                               timeParent

    !![
    <inputParameter>
      <name>massParentMinimum</name>
      <source>parameters</source>
      <description>Minimum mass of the parent halo over which to integrate.</description>
    </inputParameter>
    <inputParameter>
      <name>massParentMaximum</name>
      <source>parameters</source>
      <description>Maximum mass of the parent halo over which to integrate.</description>
    </inputParameter>
    <inputParameter>
      <name>timeParent</name>
      <source>parameters</source>
      <description>The time at which the parent halo is defined.</description>
    </inputParameter>
    <objectBuilder class="nbodyHaloMassError" name="nbodyHaloMassError_" source="parameters"/>
    <objectBuilder class="galacticFilter"     name="galacticFilter_"     source="parameters"/>
    !!]
    self=outputAnalysisDistributionOperatorMassRatioNBody(massParentMinimum,massParentMaximum,timeParent,nbodyHaloMassError_,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nbodyHaloMassError_"/>
    !!]
    return
  end function massRatioNBodyConstructorParameters

  function massRatioNBodyConstructorInternal(massParentMinimum,massParentMaximum,timeParent,nbodyHaloMassError_,galacticFilter_) result (self)
    !!{
    Internal constructor for the \refClass{outputAnalysisDistributionOperatorMassRatioNBody} output analysis distribution operator class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Node_Property_Extractors, only : nodePropertyExtractorClass, nodePropertyExtractorScalar
    implicit none
    type            (outputAnalysisDistributionOperatorMassRatioNBody)                        :: self
    double precision                                                  , intent(in   )         :: massParentMinimum  , massParentMaximum, &
         &                                                                                       timeParent
    class(nbodyHaloMassErrorClass                                    ), intent(in   ), target :: nbodyHaloMassError_
    class           (galacticFilterClass                             ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="massParentMinimum, massParentMaximum, timeParent, *nbodyHaloMassError_, *galacticFilter_"/>
    !!]

    return
  end function massRatioNBodyConstructorInternal

  subroutine massRatioNBodyDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisDistributionOperatorMassRatioNBody} output analysis weight operator class.
    !!}
    type(outputAnalysisDistributionOperatorMassRatioNBody), intent(inout) :: self

    !![
    <objectDestructor name="self%nbodyHaloMassError_"/>
    <objectDestructor name="self%galacticFilter_"    />
    !!]
    return
  end subroutine massRatioNBodyDestructor

  function massRatioNBodyOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement an N-body mass ratio output analysis distribution operator.
    !!}
    use :: Arrays_Search          , only : searchArray
    use :: Error                  , only : Error_Report
    use :: Galacticus_Nodes       , only : nodeComponentBasic
    use :: Numerical_Comparison   , only : Values_Agree
    use :: Numerical_Integration  , only : integrator
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear, outputAnalysisPropertyTypeLog10
    implicit none
    class           (outputAnalysisDistributionOperatorMassRatioNBody), intent(inout)                                        :: self
    double precision                                                  , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType       ), intent(in   )                                        :: propertyType
    double precision                                                  , intent(in   ), dimension(:)                          :: propertyValueMinimum              , propertyValueMaximum
    integer         (c_size_t                                        ), intent(in   )                                        :: outputIndex
    type            (treeNode                                        ), intent(inout)                                        :: node
    double precision                                                                 , dimension(size(propertyValueMinimum)) :: massRatioNBodyOperateScalar
    double precision                                                  , parameter                                            :: timeTolerance              =1.0d-4
    double precision                                                  , parameter                                            :: integrationExtent          =1.0d+1
    type            (treeNode                                        ), pointer                                              :: nodeParent
    class           (nodeComponentBasic                              ), pointer                                              :: basic
    double precision                                                                                                         :: massUncertaintyRatio              , massUncertaintyParent      , &
         &                                                                                                                      massUncertaintyCorrelation        , massParent                 , &
         &                                                                                                                      massRatio                         , &
         &                                                                                                                      massRatioMinimum                  , massRatioMaximum           , &
         &                                                                                                                      massParentLimitLower              , massParentLimitUpper       , &
         &                                                                                                                      massRatioLimitLower               , massRatioLimitUpper
    integer         (c_size_t                                        )                                                       :: binIndex
    type            (integrator                                      )                                                       :: integrator_
    !$GLC attributes unused :: outputIndex

    ! Get the parent halo mass.
    nodeParent => node%parent
    do while (associated(nodeParent))
       basic => nodeParent%basic()
       if (Values_Agree(basic%time(),self%timeParent,absTol=timeTolerance) .and. self%galacticFilter_%passes(nodeParent)) exit
       nodeParent => nodeParent%parent
    end do
    if (.not.associated(nodeParent)) call Error_Report('unable to find parent halo'//{introspection:location})
    massParent=basic%mass()
    ! Determine progenitor mass.
    select case (propertyType%ID)
    case (outputAnalysisPropertyTypeLinear%ID)
       massRatio=        propertyValue
    case (outputAnalysisPropertyTypeLog10 %ID)
       massRatio=10.0d0**propertyValue
    case default
       massRatio= 0.0d0
       call Error_Report('unsupported property type'//{introspection:location})
    end select
    ! Get the uncertainties and correlation of the masses.
    massUncertaintyRatio       =+self%nbodyHaloMassError_%errorFractional(node           )*massRatio
    massUncertaintyParent      =+self%nbodyHaloMassError_%errorFractional(     nodeParent)*massParent
    massUncertaintyCorrelation =+self%nbodyHaloMassError_%correlation    (node,nodeParent)
    ! Check for zero uncertainties.
    if (massUncertaintyRatio <= 0.0d0 .or. massUncertaintyParent <= 0.0d0) then
       ! Zero uncertainties - add the full weight to the bin if the progenitor and parent are within range.
       massRatioNBodyOperateScalar=0.0d0
       if (massParent <  self%massParentMinimum .or. massParent >=      self%massParentMaximum    ) return
       binIndex=searchArray(propertyValueMinimum,propertyValue)
       !! Capture the final bin.
       if (binIndex == size(propertyValueMinimum)-1_c_size_t .and. propertyValue > propertyValueMaximum(binIndex)) binIndex=size(propertyValueMinimum)
       if (binIndex   <= 0                      .or. binIndex   >  size(     propertyValueMinimum)) return
       if     (                                                 &
            &   propertyValue >= propertyValueMinimum(binIndex) &
            &  .and.                                            &
            &   propertyValue <  propertyValueMaximum(binIndex) &
            & ) massRatioNBodyOperateScalar(binIndex)=1.0d0
    else
       ! Integrate over the bivariate normal distribution to find the contribution to each bin in mass ratio.
       massParentLimitLower=max(                         &
            &                   +self%massParentMinimum, &
            &                   +massParent              &
            &                   -integrationExtent       &
            &                   *massUncertaintyParent/sqrt(1.0d0-massUncertaintyCorrelation**2)   &
            &                  )
       massParentLimitUpper=min(                         &
            &                   +self%massParentMaximum, &
            &                   +massParent              &
            &                   +integrationExtent       &
            &                   *massUncertaintyParent/sqrt(1.0d0-massUncertaintyCorrelation**2)   &
            &                  )
       if (massParentLimitUpper > massParentLimitLower) then
          integrator_=integrator(massRatioBivariateNormalIntegrand,toleranceAbsolute=1.0d-10,toleranceRelative=1.0d-03)
          do binIndex=1,size(propertyValueMinimum)
             select case (propertyType%ID)
             case (outputAnalysisPropertyTypeLinear%ID)
                massRatioMinimum=        propertyValueMinimum(binIndex)
                massRatioMaximum=        propertyValueMaximum(binIndex)
             case (outputAnalysisPropertyTypeLog10 %ID)
                massRatioMinimum=10.0d0**propertyValueMinimum(binIndex)
                massRatioMaximum=10.0d0**propertyValueMaximum(binIndex)
             case default
                massRatioMinimum= 0.0d0
                massRatioMaximum= 0.0d0
                call Error_Report('unsupported property type'//{introspection:location})
             end select
             massRatioLimitLower=+massRatioMinimum           &
                  &              -integrationExtent          &
                  &              *massUncertaintyRatio       &
                  &              +massUncertaintyCorrelation &
                  &              *(                          &
                  &                +massUncertaintyRatio     &
                  &                /massUncertaintyParent    &
                  &               )                          &
                  &              *(                          &
                  &                +massParentLimitLower     &
                  &                -massParent               &
                  &               )
             massRatioLimitUpper=+massRatioMaximum           &
                  &              +integrationExtent          &
                  &              *massUncertaintyRatio       &
                  &              +massUncertaintyCorrelation &
                  &              *(                          &
                  &                +massUncertaintyRatio     &
                  &                /massUncertaintyParent    &
                  &               )                          &
                  &              *(                          &
                  &                +massParentLimitUpper     &
                  &                -massParent               &
                  &               )
             if     (                                 &
                  &   massRatioLimitLower > massRatio &
                  &  .or.                             &
                  &   massRatioLimitUpper < massRatio &
                  & ) then
                massRatioNBodyOperateScalar(binIndex)=0.0d0
             else if (                                                                               &
                  &    +massParent-integrationExtent*massUncertaintyParent >= self%massParentMinimum &
                  &   .and.                                                                          &
                  &    +massParent+integrationExtent*massUncertaintyParent <= self%massParentMaximum &
                  &  ) then
                ! Full parent mass range is included, so use analytic integral over the mass ratio.
                massRatioNBodyOperateScalar(binIndex)=+0.5d0                                                                      &
                     &                                *(                                                                          &
                     &                                  +erf((massRatioMaximum    -massRatio)/massUncertaintyRatio  /sqrt(2.0d0)) &
                     &                                  -erf((massRatioMinimum    -massRatio)/massUncertaintyRatio  /sqrt(2.0d0)) &
                     &                                 )
             else if (                                                                               &
                  &    +massRatio -integrationExtent*massUncertaintyRatio  >= massRatioMinimum       &
                  &   .and.                                                                          &
                  &    +massRatio +integrationExtent*massUncertaintyRatio  <= massRatioMaximum       &
                  &) then
                ! Full mass ratio range is included, so use analytic integral over the parent mass.
                massRatioNBodyOperateScalar(binIndex)=+0.5d0                                                                      &
                     &                                *(                                                                          &
                     &                                  +erf((massParentLimitUpper-massParent)/massUncertaintyParent/sqrt(2.0d0)) &
                     &                                  -erf((massParentLimitLower-massParent)/massUncertaintyParent/sqrt(2.0d0)) &
                     &                                 )
             else
                massRatioNBodyOperateScalar(binIndex)=max(                                                                  &
                     &                                    integrator_%integrate(massParentLimitLower,massParentLimitUpper), &
                     &                                    0.0d0                                                             &
                     &                                   )
             end if
          end do
       else
          massRatioNBodyOperateScalar=0.0d0
       end if
    end if
    return

  contains

    double precision function massRatioBivariateNormalIntegrand(massParentPrimed)
      !!{
      Integrand used in finding the weight given to a bin in the space of parent mass vs. progenitor mass ratio. This is just
      the bivariate normal distribution of parent mass and mass ratio integrated over the required range of mass ratio. The
      remaining integral over the parent mass must then be found numerically.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, intent(in   ) :: massParentPrimed
      double precision                :: xParent         , xRatioMinimum, &
           &                             xRatioMaximum

      xParent                          =+(massParentPrimed-massParent)/massUncertaintyParent
      xRatioMinimum                    =-(massRatioMinimum-massRatio )/massUncertaintyRatio
      xRatioMaximum                    =-(massRatioMaximum-massRatio )/massUncertaintyRatio
      massRatioBivariateNormalIntegrand=+(                                                                                                                                       &
           &                              +exp(-0.5d0*xParent**2)*erf((-xRatioMaximum+xParent*massUncertaintyCorrelation)/sqrt(2.0d0)/sqrt(1.0d0-massUncertaintyCorrelation**2)) &
           &                              -exp(-0.5d0*xParent**2)*erf((-xRatioMinimum+xParent*massUncertaintyCorrelation)/sqrt(2.0d0)/sqrt(1.0d0-massUncertaintyCorrelation**2)) &
           &                             )                                                                                                                                       &
           &                            /     2.0d0                                                                                                                              &
           &                            /sqrt(2.0d0*Pi)                                                                                                                          &
           &                            /massUncertaintyParent
      return
    end function massRatioBivariateNormalIntegrand

  end function massRatioNBodyOperateScalar

  function massRatioNBodyOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node)
    !!{
    Implement an N-body mass ratio output analysis distribution operator.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisDistributionOperatorMassRatioNBody), intent(inout)                                        :: self
    double precision                                                  , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType       ), intent(in   )                                        :: propertyType
    double precision                                                  , intent(in   ), dimension(:)                          :: propertyValueMinimum          , propertyValueMaximum
    integer         (c_size_t                                        ), intent(in   )                                        :: outputIndex
    type            (treeNode                                        ), intent(inout)                                        :: node
    double precision                                                                 , dimension(size(propertyValueMinimum)) :: massRatioNBodyOperateDistribution
    !$GLC attributes unused :: self, distribution, propertyValueMinimum, propertyValueMaximum, outputIndex, propertyType, node

    massRatioNBodyOperateDistribution=0.0d0
    call Error_Report('not implemented'//{introspection:location})
    return
  end function massRatioNBodyOperateDistribution

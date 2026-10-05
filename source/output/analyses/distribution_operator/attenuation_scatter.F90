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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements an output analysis distribution operator which applies scatter in dust attenuation.
  !!}

  use :: Dust_Attenuations, only : dustAttenuationClass, dustAttenuationStellarMassRedshift

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorAttenuationScatter" docformat="rst">
   <description>
   An output analysis distribution operator which applies a truncated-normal scatter in the dust attenuation of a luminosity.
   The attenuation of each galaxy is taken to be :math:`A = \max(\bar{A}+\epsilon,0)` magnitudes, where :math:`\bar{A}` is
   the mean attenuation given by a :galacticus-class:`dustAttenuationStellarMassRedshift` object and :math:`\epsilon` is
   drawn from a normal distribution of root variance ``[rootVarianceAttenuation]``. Rather than drawing :math:`\epsilon`,
   the weight of each galaxy is distributed over the bins of the analysis in proportion to the probability of its attenuated
   luminosity falling in each: the probability, :math:`\Phi(-\bar{A}/\sigma_A)`, that :math:`A \le 0` is placed in the bin
   of its unattenuated luminosity, and the remainder is spread over the bins of lower luminosity.

   The property must be :math:`\log_{10}` of a luminosity which has already been attenuated by :math:`\max(\bar{A},0)`
   magnitudes by the same ``dustAttenuation`` object---as in the H\ :math:`\alpha` luminosity function analyses---so that
   the unattenuated luminosity is recovered by undoing that mean attenuation. Any property operators applied after the
   attenuation (e.g. a correction for cosmological luminosity distance) should therefore be simple shifts in
   :math:`\log_{10}` luminosity. The operator acts on a single value, so must be the first in any sequence of distribution
   operators.
   </description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorClass) :: outputAnalysisDistributionOperatorAttenuationScatter
     !!{RST
     An output analysis distribution operator which applies scatter in dust attenuation.
     !!}
     private
     class           (dustAttenuationStellarMassRedshift), pointer :: dustAttenuation_        => null()
     double precision                                              :: rootVarianceAttenuation
   contains
     final     ::                        attenuationScatterDestructor
     procedure :: operateScalar       => attenuationScatterOperateScalar
     procedure :: operateDistribution => attenuationScatterOperateDistribution
  end type outputAnalysisDistributionOperatorAttenuationScatter

  interface outputAnalysisDistributionOperatorAttenuationScatter
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisDistributionOperatorAttenuationScatter` output analysis distribution
     operator class.
     !!}
     module procedure attenuationScatterConstructorParameters
     module procedure attenuationScatterConstructorInternal
  end interface outputAnalysisDistributionOperatorAttenuationScatter

contains

  function attenuationScatterConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisDistributionOperatorAttenuationScatter` output analysis distribution
    operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisDistributionOperatorAttenuationScatter)                :: self
    type            (inputParameters                                     ), intent(inout) :: parameters
    class           (dustAttenuationClass                                ), pointer       :: dustAttenuation_
    double precision                                                                      :: rootVarianceAttenuation

    !![
    <inputParameter docformat="rst">
      <name>rootVarianceAttenuation</name>
      <description>
      The root variance, :math:`\sigma_A`, of the scatter in attenuation (in magnitudes).
      </description>
      <source>parameters</source>
      <minimum>0.0d0</minimum>
    </inputParameter>
    <objectBuilder class="dustAttenuation" name="dustAttenuation_" source="parameters"/>
    !!]
    self=outputAnalysisDistributionOperatorAttenuationScatter(rootVarianceAttenuation,dustAttenuation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustAttenuation_"/>
    !!]
    return
  end function attenuationScatterConstructorParameters

  function attenuationScatterConstructorInternal(rootVarianceAttenuation,dustAttenuation_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`outputAnalysisDistributionOperatorAttenuationScatter` output analysis
    distribution operator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (outputAnalysisDistributionOperatorAttenuationScatter)                        :: self
    double precision                                                      , intent(in   )         :: rootVarianceAttenuation
    class           (dustAttenuationClass                                ), intent(in   ), target :: dustAttenuation_
    !![
    <constructorAssign variables="rootVarianceAttenuation"/>
    !!]

    select type (dustAttenuation_)
    class is (dustAttenuationStellarMassRedshift)
       self%dustAttenuation_ => dustAttenuation_
       !![
       <referenceCountIncrement owner="self" object="dustAttenuation_"/>
       !!]
    class default
       call Error_Report('attenuation scatter requires a `stellarMassRedshift` dust attenuation'//{introspection:location})
    end select
    return
  end function attenuationScatterConstructorInternal

  subroutine attenuationScatterDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`outputAnalysisDistributionOperatorAttenuationScatter` output analysis distribution
    operator class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorAttenuationScatter), intent(inout) :: self

    !![
    <objectDestructor name="self%dustAttenuation_"/>
    !!]
    return
  end subroutine attenuationScatterDestructor

  function attenuationScatterOperateScalar(self,propertyValue,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node) result(distribution)
    !!{RST
    Distribute the weight of a galaxy over the bins of the analysis, given the probability distribution of its attenuation.
    !!}
    use :: Error                  , only : Error_Report
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLog10
    implicit none
    class           (outputAnalysisDistributionOperatorAttenuationScatter), intent(inout)                                        :: self
    double precision                                                      , intent(in   )                                        :: propertyValue
    type            (enumerationOutputAnalysisPropertyTypeType           ), intent(in   )                                        :: propertyType
    double precision                                                      , intent(in   ), dimension(:)                          :: propertyValueMinimum, propertyValueMaximum
    integer         (c_size_t                                            ), intent(in   )                                        :: outputIndex
    type            (treeNode                                            ), intent(inout)                                        :: node
    double precision                                                                     , dimension(size(propertyValueMinimum)) :: distribution
    double precision                                                                                                             :: attenuationMean     , luminosityUnattenuated, &
         &                                                                                                                          attenuationLower    , attenuationUpper
    integer                                                                                                                      :: i
    !$GLC attributes unused :: outputIndex

    if (propertyType /= outputAnalysisPropertyTypeLog10) call Error_Report('property must be the logarithm of a luminosity'//{introspection:location})
    distribution=0.0d0
    ! A galaxy with no (or infinite) luminosity is not placed in any bin.
    if (propertyValue == +huge(0.0d0) .or. propertyValue == -huge(0.0d0)) return
    ! Recover the unattenuated luminosity, by undoing the mean attenuation (limited to be non-negative) already applied.
    attenuationMean       =self%dustAttenuation_%attenuationMean(node)
    luminosityUnattenuated=propertyValue+0.4d0*max(attenuationMean,0.0d0)
    if (self%rootVarianceAttenuation <= 0.0d0) then
       ! No scatter - the attenuation is exactly max(Ā,0), and the galaxy lies in the bin containing its attenuated luminosity.
       where (propertyValue >= propertyValueMinimum .and. propertyValue < propertyValueMaximum)
          distribution=1.0d0
       end where
       return
    end if
    do i=1,size(propertyValueMinimum)
       ! The range of (positive) attenuation which places the attenuated luminosity in this bin.
       attenuationLower=max((luminosityUnattenuated-propertyValueMaximum(i))/0.4d0,0.0d0)
       attenuationUpper=max((luminosityUnattenuated-propertyValueMinimum(i))/0.4d0,0.0d0)
       if (attenuationUpper > attenuationLower)                                                                  &
            & distribution(i)=+cumulativeNormal((attenuationUpper-attenuationMean)/self%rootVarianceAttenuation) &
            &                 -cumulativeNormal((attenuationLower-attenuationMean)/self%rootVarianceAttenuation)
       ! The probability of zero attenuation is placed in the bin containing the unattenuated luminosity.
       if (luminosityUnattenuated >= propertyValueMinimum(i) .and. luminosityUnattenuated < propertyValueMaximum(i)) &
            & distribution(i)=+distribution(i)                                                                       &
            &                 +cumulativeNormal(-attenuationMean/self%rootVarianceAttenuation)
    end do
    return

  contains

    double precision function cumulativeNormal(x)
      !!{RST
      The cumulative distribution function of the standard normal distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: x

      cumulativeNormal=0.5d0*erfc(-x/sqrt(2.0d0))
      return
    end function cumulativeNormal

  end function attenuationScatterOperateScalar

  function attenuationScatterOperateDistribution(self,distribution,propertyType,propertyValueMinimum,propertyValueMaximum,outputIndex,node) result(distributionNew)
    !!{RST
    Attenuation scatter can be applied only to a single value, not to a distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (outputAnalysisDistributionOperatorAttenuationScatter), intent(inout)                                        :: self
    double precision                                                      , intent(in   ), dimension(:)                          :: distribution
    type            (enumerationOutputAnalysisPropertyTypeType           ), intent(in   )                                        :: propertyType
    double precision                                                      , intent(in   ), dimension(:)                          :: propertyValueMinimum, propertyValueMaximum
    integer         (c_size_t                                            ), intent(in   )                                        :: outputIndex
    type            (treeNode                                            ), intent(inout)                                        :: node
    double precision                                                                     , dimension(size(propertyValueMinimum)) :: distributionNew
    !$GLC attributes unused :: self, distribution, propertyType, propertyValueMaximum, outputIndex, node

    distributionNew=0.0d0
    call Error_Report('attenuation scatter must be the first operator in a sequence, as it can not act on a distribution'//{introspection:location})
    return
  end function attenuationScatterOperateDistribution

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
Implements a merger tree branching probability class using the algorithm of \cite{parkinson_generating_2008} plus an additional term.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Tables                    , only : table1DLogarithmicLinear

  !![
  <mergerTreeBranchingProbability name="mergerTreeBranchingProbabilityPCHPlus">
   <description>Merger tree branching probabilities using the algorithm of \cite{parkinson_generating_2008} plus an additional term.</description>
  </mergerTreeBranchingProbability>
  !!]
  type, extends(mergerTreeBranchingProbabilityParkinsonColeHelly) :: mergerTreeBranchingProbabilityPCHPlus
     !!{
     A merger tree branching probability class using the algorithm of \cite{parkinson_generating_2008} plus an additional term.
     !!}
     private
     double precision :: gamma3
   contains
     procedure :: V               => pchPlusV
     procedure :: modifier        => pchPlusModifier
     procedure :: hypergeometricA => pchPlusHypergeometricA
  end type mergerTreeBranchingProbabilityPCHPlus

  interface mergerTreeBranchingProbabilityPCHPlus
     !!{
     Constructors for the {\normalfont \ttfamily pchPlus} merger tree builder class.
     !!}
     module procedure pchPlusConstructorParameters
     module procedure pchPlusConstructorInternal
  end interface mergerTreeBranchingProbabilityPCHPlus

contains

  function pchPlusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily pchPlus} merger tree branching probability class which reads parameters from a provided
    parameter list.
    !!}
    implicit none
    type            (mergerTreeBranchingProbabilityPCHPlus)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologicalMassVarianceClass        ), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass             ), pointer       :: criticalOverdensity_
    double precision                                                       :: gamma1                   , gamma2            , &
         &                                                                    G0                       , accuracyFirstOrder, &
         &                                                                    precisionHypergeometric  , gamma3
    logical                                                                :: hypergeometricTabulate   , cdmAssumptions    , &
         &                                                                    tolerateRoundOffErrors

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>G0</name>
      <defaultValue>0.57d0</defaultValue>
      <description>The parameter $G_0$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gamma1</name>
      <defaultValue>0.38d0</defaultValue>
      <description>The parameter $\gamma_1$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gamma2</name>
      <defaultValue>-0.01d0</defaultValue>
      <description>The parameter $\gamma_2$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>gamma3</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $\gamma_32$ appearing in the modified merger rate expression of \cite{parkinson_generating_2008}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>accuracyFirstOrder</name>
      <defaultValue>0.1d0</defaultValue>
      <description>Limits the step in $\delta_\mathrm{crit}$ when constructing merger trees using the \cite{parkinson_generating_2008}
         algorithm, so that it never exceeds {\normalfont \ttfamily accuracyFirstOrder}$\sqrt{2[\sigma^2(M_2/2)-\sigma^2(M_2)]}$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>precisionHypergeometric</name>
      <defaultValue>1.0d-6</defaultValue>
      <description>The fractional precision required in evaluates of hypergeometric functions in the modified Press-Schechter tree branching calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>hypergeometricTabulate</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether hypergeometric factors should be precomputed and tabulated in modified Press-Schechter tree branching functions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>cdmAssumptions</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, assume that $\alpha(=-\mathrm{d}\log \sigma/\mathrm{d}\log M)&gt;0$ and $\mathrm{d}\alpha/\mathrm{d}M&gt;0$ (as is true in the case of \gls{cdm}) when constructing merger trees using the \cite{parkinson_generating_2008}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>tolerateRoundOffErrors</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, round-off errors in integrations of branching probability will be tolerated. This may degrade the accuracy of solutions, but can be unavoidable in models with cut-offs in their power spectra.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    !!]
    self=mergerTreeBranchingProbabilityPCHPlus(G0,gamma1,gamma2,gamma3,accuracyFirstOrder,precisionHypergeometric,hypergeometricTabulate,cdmAssumptions,tolerateRoundOffErrors,cosmologicalMassVariance_,criticalOverdensity_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    !!]
    return
  end function pchPlusConstructorParameters

  function pchPlusConstructorInternal(G0,gamma1,gamma2,gamma3,accuracyFirstOrder,precisionHypergeometric,hypergeometricTabulate,cdmAssumptions,tolerateRoundOffErrors,cosmologicalMassVariance_,criticalOverdensity_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily pchPlus} merger tree branching probability class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (mergerTreeBranchingProbabilityPCHPlus)                        :: self
    double precision                                       , intent(in   )         :: gamma1                   , gamma2            , &
         &                                                                            G0                       , accuracyFirstOrder, &
         &                                                                            precisionHypergeometric  , gamma3
    logical                                                , intent(in   )         :: hypergeometricTabulate   , cdmAssumptions    , &
         &                                                                            tolerateRoundOffErrors
    class           (cosmologicalMassVarianceClass        ), intent(in   ), target :: cosmologicalMassVariance_
    class           (criticalOverdensityClass             ), intent(in   ), target :: criticalOverdensity_
    !![
    <constructorAssign variables=" gamma3"/>
    !!]

    ! Validate.
    if (cdmAssumptions .and. gamma3 > 1.5d0) call Error_Report('γ₃>³/₂ violates CDM assumptions'//{introspection:location})
    ! Initialize.
    self%mergerTreeBranchingProbabilityParkinsonColeHelly=mergerTreeBranchingProbabilityParkinsonColeHelly(G0,gamma1,gamma2,accuracyFirstOrder,precisionHypergeometric,hypergeometricTabulate,cdmAssumptions,tolerateRoundOffErrors,cosmologicalMassVariance_,criticalOverdensity_)
    return
  end function pchPlusConstructorInternal

  double precision function pchPlusV(self,massFraction,haloMass)
    !!{
    The function $V(q)$ from \cite[][eqn. A4]{parkinson_generating_2008}.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityPCHPlus), intent(inout) :: self
    double precision                                       , intent(in   ) :: massFraction     , haloMass
    double precision                                                       :: childSigmaSquared

    childSigmaSquared=+self%cosmologicalMassVariance_%rootVariance(massFraction*haloMass,self%timeParent)**2
    pchPlusV         =+       childSigmaSquared  &
         &            /(                         &
         &              +     childSigmaSquared  &
         &              -self%sigmaParentSquared &
         &             )**1.5d0                  &
         &            *(                         &
         &              +1.0d0                   &
         &              -self%sigmaParentSquared &
         &              /     childSigmaSquared  &
         &             )**self%gamma3
    return
  end function pchPlusV

  double precision function pchPlusModifier(self,childSigma)
    !!{
    Empirical modification of the progenitor mass function from \cite{parkinson_generating_2008}. The constant factors of
    $G_0 (\delta_\mathrm{p}/\sigma_\mathrm{p})^{\gamma_2}$ and $1/\sigma_\mathrm{p}^{\gamma_1}$ are not included
    here---instead they are included in a multiplicative prefactor by which integrals over this function are multiplied.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityPCHPlus), intent(inout) :: self
    double precision                                       , intent(in   ) :: childSigma

    pchPlusModifier=childSigma**self%gamma1*(1.0d0-self%sigmaParentSquared/childSigma**2)**self%gamma3
    return
  end function pchPlusModifier

  function pchPlusHypergeometricA(self,gamma) result(a)
    !!{
    Compute the $a$ parameter of the hypergeometric function.
    !!}
    implicit none
    double precision                                       , dimension(2)  :: a
    class           (mergerTreeBranchingProbabilityPCHPlus), intent(inout) :: self
    double precision                                       , intent(in   ) :: gamma

    a=[0.5d0-0.5d0*gamma,1.5d0-self%gamma3]
    return
  end function pchPlusHypergeometricA

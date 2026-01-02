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
  Implements a merger tree branching probability rate modifier which uses the model of \cite{parkinson_generating_2008} plus an additional term.
  !!}

  !![
  <mergerTreeBranchingProbabilityModifier name="mergerTreeBranchingProbabilityModifierPCHPlus">
   <description>
    Provides a merger tree branching probability rate modifier which uses the model of \cite{parkinson_generating_2008} plus an
    additional term. Specifically, the modifier becomes
    \begin{equation}
     G\left( {\sigma_1 \over \sigma_2} , {\delta_2 \over \sigma_2} \right) =
     G_0
     \left({\sigma_1\over\sigma_2}\right)^{\gamma_1}
     \left({\delta_2\over\sigma_2}\right)^{\gamma_2}
     \left(1 - {\sigma_2^2 \over \sigma_1^2}\right)^{\gamma_3},
    \end{equation}
    where $\sigma_i=\sigma(M_i)$, $\sigma(M)$ is the usual present-day, linear-theory mass-variance in spheres enclosing an average
    mass $M$, $M_2$ is the mass of the parent halo, $M_1$ is the mass of the child halo, and $\delta_2$ is the critical overdensity
    for collapse at the epoch of the parent.
   </description>
  </mergerTreeBranchingProbabilityModifier>
  !!]
  type, extends(mergerTreeBranchingProbabilityModifierParkinson2008) :: mergerTreeBranchingProbabilityModifierPCHPlus
     !!{
     A merger tree branching probability rate modifier which uses the model of \cite{parkinson_generating_2008} plus an additional term.
     !!}
     private
     double precision :: gamma3
   contains
     procedure :: rateModifier => pchPlusRateModifier
  end type mergerTreeBranchingProbabilityModifierPCHPlus

  interface mergerTreeBranchingProbabilityModifierPCHPlus
     !!{
     Constructors for the \refClass{mergerTreeBranchingProbabilityModifierPCHPlus} merger tree branching probability rate class.
     !!}
     module procedure pchPlusConstructorParameters
     module procedure pchPlusConstructorInternal
  end interface mergerTreeBranchingProbabilityModifierPCHPlus

contains

  function pchPlusConstructorParameters(parameters) result(self)
    !!{
    A constructor for the {\normalfont \ttfamily pchPlus} merger tree branching probability rate class which builds the
    object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBranchingProbabilityModifierPCHPlus)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (criticalOverdensityClass                     ), pointer       :: criticalOverdensity_
    double precision                                                               :: G0                  , gamma1, &
         &                                                                            gamma2              , gamma3

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
      <description>The parameter $\gamma_3$ appearing in the extension of the modified merger rate expression of \cite{parkinson_generating_2008} as defined in this class.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="criticalOverdensity" name="criticalOverdensity_" source="parameters"/>
    !!]
    self=mergerTreeBranchingProbabilityModifierPCHPlus(G0,gamma1,gamma2,gamma3,criticalOverdensity_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"/>
    !!]
    return
  end function pchPlusConstructorParameters

  function pchPlusConstructorInternal(G0,gamma1,gamma2,gamma3,criticalOverdensity_) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily pchPlus} merger tree branching probability rate class.
    !!}
    implicit none
    type            (mergerTreeBranchingProbabilityModifierPCHPlus)                        :: self
    class           (criticalOverdensityClass                     ), intent(in   ), target :: criticalOverdensity_
    double precision                                               , intent(in   )         :: G0                  , gamma1, &
         &                                                                                    gamma2              , gamma3
    !![
    <constructorAssign variables="gamma3"/>
    !!]

    self%mergerTreeBranchingProbabilityModifierParkinson2008=mergerTreeBranchingProbabilityModifierParkinson2008(G0,gamma1,gamma2,criticalOverdensity_)
    return
  end function pchPlusConstructorInternal

  double precision function pchPlusRateModifier(self,nodeParent,massParent,sigmaParent,sigmaChild,timeParent)
    !!{
    Returns a modifier for merger tree branching rates using the \cite{parkinson_generating_2008} algorithm plus an additional term.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityModifierPCHPlus), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: nodeParent
    double precision                                               , intent(in   ) :: sigmaChild , timeParent, &
         &                                                                            sigmaParent, massParent
    
    ! Compute the modifier.
    pchPlusRateModifier=+self%mergerTreeBranchingProbabilityModifierParkinson2008%rateModifier(nodeParent,massParent,sigmaParent,sigmaChild,timeParent) &
         &              *(                                                                                                                              &
         &                +1.0d0                                                                                                                        &
         &                -(                                                                                                                            &
         &                  +sigmaParent                                                                                                                &
         &                  /sigmaChild                                                                                                                 &
         &                 )**2                                                                                                                         &
         &               )**self%gamma3
    return
  end function pchPlusRateModifier

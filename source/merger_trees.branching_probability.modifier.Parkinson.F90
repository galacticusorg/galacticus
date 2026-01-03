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
Implements a merger tree branching probability rate modifier which uses the model of \cite{parkinson_generating_2008}.
!!}

  use :: Cosmological_Density_Field, only : criticalOverdensity, criticalOverdensityClass

  !![
  <mergerTreeBranchingProbabilityModifier name="mergerTreeBranchingProbabilityModifierParkinson2008">
   <description>
    A merger tree branching probability modifier class which uses the model of \cite{parkinson_generating_2008}. Specifically,
    the branching probability is multiplied by:
    \begin{equation}
     G(\delta_\mathrm{p},\sigma_\mathrm{c},\sigma_\mathrm{p}) = G_0
     \left({\sigma_\mathrm{c}\over\sigma_\mathrm{p}}\right)^{\gamma_1}
     \left({\delta_\mathrm{p}\over\sigma_\mathrm{p}}\right)^{\gamma_2}
    \end{equation}
    where $\delta_\mathrm{p}$ is the current critical overdensity for collapse for the parent halo, and $\sigma_\mathrm{c}$ and
    $\sigma_\mathrm{p}$ are the root-variance of the smooth mass-density field on scales corresponding to the masses of child
    and parent halos respectively. The parameters of the fit can be adjusted via input parameters: $G_0=${\normalfont \ttfamily
    [G0]}, $\gamma_1=${\normalfont \ttfamily [gamma1]}, and $\gamma_2=${\normalfont \ttfamily [gamma2]}.
   </description>
  </mergerTreeBranchingProbabilityModifier>
  !!]
  type, extends(mergerTreeBranchingProbabilityModifierClass) :: mergerTreeBranchingProbabilityModifierParkinson2008
     !!{
     A merger tree branching probability rate modifier which uses the model of \cite{parkinson_generating_2008}.
     !!}
     private
     class           (criticalOverdensityClass), pointer :: criticalOverdensity_ => null()
     double precision                                    :: G0                            , gamma1             , &
          &                                                 gamma2                        , parentTerm         , &
          &                                                 timeParentPrevious            , sigmaParentPrevious
   contains
     final     ::                 parkinson2008Destructor
     procedure :: rateModifier => parkinson2008RateModifier
  end type mergerTreeBranchingProbabilityModifierParkinson2008

  interface mergerTreeBranchingProbabilityModifierParkinson2008
     !!{
     Constructors for the \refClass{mergerTreeBranchingProbabilityModifierParkinson2008} merger tree branching probability rate class.
     !!}
     module procedure parkinson2008ConstructorParameters
     module procedure parkinson2008ConstructorInternal
  end interface mergerTreeBranchingProbabilityModifierParkinson2008

contains

  function parkinson2008ConstructorParameters(parameters) result(self)
    !!{
    A constructor for the {\normalfont \ttfamily parkinson2008} merger tree branching probability rate class which builds the
    object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBranchingProbabilityModifierParkinson2008)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    class           (criticalOverdensityClass                           ), pointer       :: criticalOverdensity_
    double precision                                                                     :: G0                  , gamma1, &
         &                                                                                  gamma2

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
    <objectBuilder class="criticalOverdensity" name="criticalOverdensity_" source="parameters"/>
    !!]
    self=mergerTreeBranchingProbabilityModifierParkinson2008(G0,gamma1,gamma2,criticalOverdensity_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"/>
    !!]
    return
  end function parkinson2008ConstructorParameters

  function parkinson2008ConstructorInternal(G0,gamma1,gamma2,criticalOverdensity_) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily parkinson2008} merger tree branching probability rate class.
    !!}
    implicit none
    type            (mergerTreeBranchingProbabilityModifierParkinson2008)                        :: self
    class           (criticalOverdensityClass                           ), intent(in   ), target :: criticalOverdensity_
    double precision                                                     , intent(in   )         :: G0                  , gamma1, &
         &                                                                                          gamma2
    !![
    <constructorAssign variables="G0, gamma1, gamma2, *criticalOverdensity_"/>
    !!]

    self%timeParentPrevious =-1.0d0
    self%sigmaParentPrevious=-1.0d0
    self%parentTerm         =+0.0d0
    return
  end function parkinson2008ConstructorInternal

  subroutine parkinson2008Destructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBranchingProbabilityModifierParkinson2008} merger tree branching probability rate class.
    !!}
    implicit none
    type(mergerTreeBranchingProbabilityModifierParkinson2008), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"/>
    !!]
    return
  end subroutine parkinson2008Destructor

  double precision function parkinson2008RateModifier(self,nodeParent,massParent,sigmaParent,sigmaChild,timeParent)
    !!{
    Returns a modifier for merger tree branching rates using the \cite{parkinson_generating_2008} algorithm.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityModifierParkinson2008), intent(inout) :: self
    type            (treeNode                                           ), intent(inout) :: nodeParent
    double precision                                                     , intent(in   ) :: sigmaChild , timeParent, &
         &                                                                                  sigmaParent, massParent

    ! Check if we need to update the "parent" term.
    if     (                                         &
         &   timeParent  /= self%timeParentPrevious  &
         &  .or.                                     &
         &   sigmaParent /= self%sigmaParentPrevious &
         & ) then
       ! "Parent" term must be updated. Compute and store it for future re-use.
       self%timeParentPrevious =+timeParent
       self%sigmaParentPrevious=+sigmaParent
       self%parentTerm         =+self%G0                                                                                                       &
            &                   *((self%criticalOverdensity_%value(time=timeParent,mass=massParent,node=nodeParent)/sigmaParent)**self%gamma2) &
            &                   /(                                                                                  sigmaParent **self%gamma1)
    end if
    ! Compute the modifier.
    parkinson2008RateModifier=self%parentTerm*(sigmaChild**self%gamma1)
    return
  end function parkinson2008RateModifier

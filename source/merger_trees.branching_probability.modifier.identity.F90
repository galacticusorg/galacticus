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
Implements a merger tree branching probability rate modifier which always returns the identity modifier.
!!}

  !![
  <mergerTreeBranchingProbabilityModifier name="mergerTreeBranchingProbabilityModifierIdentity">
   <description>
    A merger tree branching probability modifier class which always applies the identity modifier.
   </description>
  </mergerTreeBranchingProbabilityModifier>
  !!]
  type, extends(mergerTreeBranchingProbabilityModifierClass) :: mergerTreeBranchingProbabilityModifierIdentity
     !!{
     A merger tree branching probability rate modifier which always returns the identity modifier.
     !!}
     private
   contains
     procedure :: rateModifier => identityRateModifier
  end type mergerTreeBranchingProbabilityModifierIdentity

  interface mergerTreeBranchingProbabilityModifierIdentity
     !!{
     Constructors for the \refClass{mergerTreeBranchingProbabilityModifierIdentity} merger tree branching probability rate modifier class.
     !!}
     module procedure identityConstructorParameters
  end interface mergerTreeBranchingProbabilityModifierIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{
    A constructor for the {\normalfont \ttfamily identity} merger tree branching probability rate class which builds the
    object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeBranchingProbabilityModifierIdentity)                :: self
    type(inputParameters                               ), intent(inout) :: parameters

    self=mergerTreeBranchingProbabilityModifierIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

  double precision function identityRateModifier(self,nodeParent,massParent,sigmaParent,sigmaChild,timeParent)
    !!{
    Returns a modifier for merger tree branching rates using the \cite{parkinson_generating_2008} algorithm.
    Return the core radius of the hot halo mass distribution.
    !!}
    implicit none
    class           (mergerTreeBranchingProbabilityModifierIdentity), intent(inout) :: self
    type            (treeNode                                      ), intent(inout) :: nodeParent
    double precision                                                , intent(in   ) :: sigmaChild , timeParent, &
         &                                                                             sigmaParent, massParent
    !$GLC attributes unused :: self, nodeParent, massParent, sigmaParent, sigmaChild, timeParent

    identityRateModifier=1.0d0
    return
  end function identityRateModifier

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
  Implementation of a posterior sampling differential evolution proposal size class in which the proposal size is fixed.
  !!}

  !![
  <posteriorSampleDffrntlEvltnProposalSize name="posteriorSampleDffrntlEvltnProposalSizeFixed">
   <description>
    A posterior sampling differential evolution proposal size class in which the proposal size is a fixed value
    $\gamma=${\normalfont \ttfamily [gamma]}.
   </description>
  </posteriorSampleDffrntlEvltnProposalSize>
  !!]
  type, extends(posteriorSampleDffrntlEvltnProposalSizeClass) :: posteriorSampleDffrntlEvltnProposalSizeFixed
     !!{
     Implementation of a posterior sampling differential evolution proposal size class in which the proposal size is fixed.
     !!}
     private
     double precision :: proposalSize
   contains
     procedure :: gamma => fixedGamma
  end type posteriorSampleDffrntlEvltnProposalSizeFixed

  interface posteriorSampleDffrntlEvltnProposalSizeFixed
     !!{
     Constructors for the \refClass{posteriorSampleDffrntlEvltnProposalSizeFixed} posterior sampling differential evolution random jump class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface posteriorSampleDffrntlEvltnProposalSizeFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleDffrntlEvltnProposalSizeFixed} posterior sampling differential evolution random jump class which builds
    the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleDffrntlEvltnProposalSizeFixed)                 :: self
    type            (inputParameters                             ), intent(inout)  :: parameters
    double precision                                                               :: proposalSize

    !![
    <inputParameter>
      <name>proposalSize</name>
      <description>The proposal size, $\gamma$.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleDffrntlEvltnProposalSizeFixed(proposalSize)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(proposalSize) result(self)
    !!{
    Internal constructor for the \refClass{posteriorSampleDffrntlEvltnProposalSizeFixed} posterior sampling differential evolution random jump class.
    !!}
    implicit none
    type            (posteriorSampleDffrntlEvltnProposalSizeFixed)                :: self
    double precision                                              , intent(in   ) :: proposalSize
    !![
    <constructorAssign variables="proposalSize"/>
    !!]

    return
  end function fixedConstructorInternal

  double precision function fixedGamma(self,simulationState,simulationConvergence)
    !!{
    Return the current state.
    !!}
    implicit none
    class(posteriorSampleDffrntlEvltnProposalSizeFixed), intent(inout) :: self
    class(posteriorSampleStateClass                   ), intent(inout) :: simulationState
    class(posteriorSampleConvergenceClass             ), intent(inout) :: simulationConvergence
    !$GLC attributes unused :: simulationState, simulationConvergence

    fixedGamma=self%proposalSize
    return
  end function fixedGamma

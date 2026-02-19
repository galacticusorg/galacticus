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
  Implements a merger trees analysis class which wraps another class and ``heats'' its likelihood to a given temperature.
  !!}

  !![
  <outputAnalysis name="outputAnalysisHeatedLikelihood">
   <description>A merger tree analysis class which wraps another class and ``heats'' its likelihood to a given temperature.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisHeatedLikelihood
     !!{
     Implementation of a merger tree analysis class which wraps another class and ``heats'' its likelihood to a given temperature.
     !!}
     private
     class           (outputAnalysisClass), pointer :: outputAnalysis_ => null()
     double precision                               :: temperature
   contains
     final     ::                  heatedLikelihoodDestructor
     procedure :: newTree       => heatedLikelihoodNewTree
     procedure :: analyze       => heatedLikelihoodAnalyze
     procedure :: finalize      => heatedLikelihoodFinalize
     procedure :: reduce        => heatedLikelihoodReduce
     procedure :: logLikelihood => heatedLikelihoodLoglikelihood
  end type outputAnalysisHeatedLikelihood

  interface outputAnalysisHeatedLikelihood
     !!{
     Constructors for the \refClass{outputAnalysisHeatedLikelihood} merger tree analysis.
     !!}
     module procedure heatedLikelihoodConstructorParameters
     module procedure heatedLikelihoodConstructorInternal
  end interface outputAnalysisHeatedLikelihood

contains

  function heatedLikelihoodConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisHeatedLikelihood} merger tree analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisHeatedLikelihood)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (outputAnalysisClass           ), pointer       :: outputAnalysis_
    double precision                                                :: temperature

    !![
    <inputParameter>
      <name>temperature</name>
      <description>The temperature to which to heat the likelihood.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="outputAnalysis" name="outputAnalysis_" source="parameters"/>
    !!]
    self=outputAnalysisHeatedLikelihood(temperature,outputAnalysis_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputAnalysis_"     />
    !!]
    return
  end function heatedLikelihoodConstructorParameters

  function heatedLikelihoodConstructorInternal(temperature,outputAnalysis_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisHeatedLikelihood} analysis class.
    !!}
    implicit none
    type            (outputAnalysisHeatedLikelihood)                        :: self
    class           (outputAnalysisClass           ), intent(in   ), target :: outputAnalysis_
    double precision                                , intent(in   )         :: temperature
    !![
    <constructorAssign variables="temperature, *outputAnalysis_"/>
    !!]
    return
  end function heatedLikelihoodConstructorInternal

  subroutine heatedLikelihoodDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisHeatedLikelihood} analysis class.
    !!}
    implicit none
    type(outputAnalysisHeatedLikelihood), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"/>
    !!]
    return
  end subroutine heatedLikelihoodDestructor

  subroutine heatedLikelihoodNewTree(self,tree,iOutput)
    !!{
    Output from all analyses.
    !!}
    implicit none
    class  (outputAnalysisHeatedLikelihood), intent(inout) :: self
    type   (mergerTree                    ), intent(inout) :: tree
    integer(c_size_t                      ), intent(in   ) :: iOutput

    call self%outputAnalysis_%newTree(tree,iOutput)
    return
  end subroutine heatedLikelihoodNewTree

  subroutine heatedLikelihoodAnalyze(self,node,iOutput)
    !!{
    Output from all analyses.
    !!}
    implicit none
    class  (outputAnalysisHeatedLikelihood), intent(inout) :: self
    type   (treeNode                      ), intent(inout) :: node
    integer(c_size_t                      ), intent(in   ) :: iOutput

    call self%outputAnalysis_%analyze(node,iOutput)
    return
  end subroutine heatedLikelihoodAnalyze

  subroutine heatedLikelihoodFinalize(self,groupName)
    !!{
    Finalize all analyses.
    !!}
    implicit none
    class(outputAnalysisHeatedLikelihood), intent(inout)           :: self
    type (varying_string                ), intent(in   ), optional :: groupName

    call self%outputAnalysis_%finalize(groupName)
    return
  end subroutine heatedLikelihoodFinalize

  double precision function heatedLikelihoodLogLikelihood(self)
    !!{
    Find the log-likelihood over all analyses.
    !!}
    implicit none
    class(outputAnalysisHeatedLikelihood), intent(inout) :: self

    heatedLikelihoodLogLikelihood=+self%outputAnalysis_%logLikelihood() &
         &                        /self%temperature
    return
  end function heatedLikelihoodLogLikelihood

  subroutine heatedLikelihoodReduce(self,reduced)
    !!{
    Reduce over the analysis.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisHeatedLikelihood), intent(inout) :: self
    class(outputAnalysisClass           ), intent(inout) :: reduced

    select type (reduced)
    type is (outputAnalysisHeatedLikelihood)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine heatedLikelihoodReduce

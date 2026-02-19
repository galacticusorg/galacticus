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
  Implements a merger trees analysis class which combines multiple other analyses.
  !!}

  type, public :: multiAnalysisList
     class(outputAnalysisClass), pointer :: analysis_ => null()
     type (multiAnalysisList  ), pointer :: next      => null()
  end type multiAnalysisList

  !![
  <outputAnalysis name="outputAnalysisMulti">
   <description>A merger tree analysis class which combines multiple other analyses.</description>
   <linkedList type="multiAnalysisList" variable="analyses" next="next" object="analysis_" objectType="outputAnalysisClass"/>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisMulti
     !!{
     Implementation of a merger tree analysis class which combines multiple other analyses.
     !!}
     private
     type(multiAnalysisList), pointer :: analyses => null()
   contains
     final     ::                  multiDestructor
     procedure :: newTree       => multiNewTree
     procedure :: analyze       => multiAnalyze
     procedure :: finalize      => multiFinalize
     procedure :: reduce        => multiReduce
     procedure :: logLikelihood => multiLoglikelihood
  end type outputAnalysisMulti

  interface outputAnalysisMulti
     !!{
     Constructors for the \refClass{outputAnalysisMulti} merger tree analysis.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface outputAnalysisMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisMulti} merger tree analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (outputAnalysisMulti)                :: self
    type   (inputParameters    ), intent(inout) :: parameters
    type   (multiAnalysisList  ), pointer       :: analysis_
    integer                                     :: i

    self      %analyses => null()
    analysis_           => null()
    do i=1,parameters%copiesCount('outputAnalysis',zeroIfNotPresent=.true.)
       if (associated(analysis_)) then
          allocate(analysis_%next)
          analysis_ => analysis_%next
       else
          allocate(self%analyses)
          analysis_ => self%analyses
       end if
       !![
       <objectBuilder class="outputAnalysis" name="analysis_%analysis_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="outputAnalysis"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(analyses) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisMulti} analysis class.
    !!}
    implicit none
    type(outputAnalysisMulti)                        :: self
    type(multiAnalysisList  ), target, intent(in   ) :: analyses
    type(multiAnalysisList  ), pointer               :: analysis_

    self   %analyses => analyses
    analysis_        => analyses
    do while (associated(analysis_))
       !![
       <referenceCountIncrement owner="analysis_" object="analysis_"/>
       !!]
       analysis_ => analysis_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisMulti} analysis class.
    !!}
    implicit none
    type(outputAnalysisMulti), intent(inout) :: self
    type(multiAnalysisList  ), pointer       :: analysis_, analysisNext

    if (associated(self%analyses)) then
       analysis_ => self%analyses
       do while (associated(analysis_))
          analysisNext => analysis_%next
          !![
          <objectDestructor name="analysis_%analysis_"/>
          !!]
          deallocate(analysis_)
          analysis_ => analysisNext
       end do
    end if
    return
  end subroutine multiDestructor

  subroutine multiNewTree(self,tree,iOutput)
    !!{
    Output from all analyses.
    !!}
    implicit none
    class  (outputAnalysisMulti), intent(inout) :: self
    type   (mergerTree         ), intent(inout) :: tree
    integer(c_size_t           ), intent(in   ) :: iOutput
    type   (multiAnalysisList  ), pointer       :: analysis_

    analysis_ => self%analyses
    do while (associated(analysis_))
       call analysis_%analysis_%newTree(tree,iOutput)
       analysis_ => analysis_%next
    end do
    return
  end subroutine multiNewTree

  subroutine multiAnalyze(self,node,iOutput)
    !!{
    Output from all analyses.
    !!}
    implicit none
    class  (outputAnalysisMulti), intent(inout) :: self
    type   (treeNode           ), intent(inout) :: node
    integer(c_size_t           ), intent(in   ) :: iOutput
    type   (multiAnalysisList  ), pointer       :: analysis_

    analysis_ => self%analyses
    do while (associated(analysis_))
       call analysis_%analysis_%analyze(node,iOutput)
       analysis_ => analysis_%next
    end do
    return
  end subroutine multiAnalyze

  subroutine multiFinalize(self,groupName)
    !!{
    Finalize all analyses.
    !!}
    implicit none
    class(outputAnalysisMulti), intent(inout)           :: self
    type (varying_string     ), intent(in   ), optional :: groupName
    type (multiAnalysisList  ), pointer                 :: analysis_

    analysis_ => self%analyses
    do while (associated(analysis_))
       call analysis_%analysis_%finalize(groupName)
       analysis_ => analysis_%next
    end do
    return
  end subroutine multiFinalize

  double precision function multiLogLikelihood(self)
    !!{
    Find the log-likelihood over all analyses. This assumes that the analyses are independent.
    !!}
    implicit none
    class(outputAnalysisMulti), intent(inout) :: self
    type (multiAnalysisList  ), pointer       :: analysis_

    multiLogLikelihood =  0.0d0
    analysis_          => self%analyses
    do while (associated(analysis_))
       multiLogLikelihood=multiLogLikelihood+analysis_%analysis_%logLikelihood()
       analysis_ => analysis_%next
    end do
    return
  end function multiLogLikelihood

  subroutine multiReduce(self,reduced)
    !!{
    Reduce over the analysis.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(outputAnalysisMulti), intent(inout) :: self
    class(outputAnalysisClass), intent(inout) :: reduced
    type (multiAnalysisList  ), pointer       :: analysis_, analysisReduced_

    select type (reduced)
    type is (outputAnalysisMulti)
       analysis_        => self   %analyses
       analysisReduced_ => reduced%analyses
       do while (associated(analysis_))
          call analysis_%analysis_%reduce(analysisReduced_%analysis_)
          analysis_        => analysis_%next
          analysisReduced_ => analysisReduced_%next
       end do
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine multiReduce

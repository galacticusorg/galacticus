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
Implements a null output analysis class.
!!}

  !![
  <outputAnalysis name="outputAnalysisNull">
   <description>A null output analysis class.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisNull
     !!{
     A null output analysis class.
     !!}
     private
   contains
     procedure :: analyze       => nullAnalyze
     procedure :: finalize      => nullFinalize
     procedure :: reduce        => nullReduce
     procedure :: logLikelihood => nullLogLikelihood
  end type outputAnalysisNull

  interface outputAnalysisNull
     !!{
     Constructors for the \refClass{outputAnalysisNull} output analysis class.
     !!}
     module procedure nullConstructorParameters
  end interface outputAnalysisNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisNull} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(outputAnalysisNull)                :: self
    type(inputParameters   ), intent(inout) :: parameters

    self=outputAnalysisNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  subroutine nullAnalyze(self,node,iOutput)
    !!{
    Implement a null output analysis.
    !!}
    implicit none
    class  (outputAnalysisNull), intent(inout) :: self
    type   (treeNode          ), intent(inout) :: node
    integer(c_size_t          ), intent(in   ) :: iOutput
    !$GLC attributes unused :: self, node, iOutput

    return
  end subroutine nullAnalyze

  subroutine nullFinalize(self,groupName)
    !!{
    Implement a null output analysis finalization.
    !!}
    implicit none
    class(outputAnalysisNull), intent(inout)           :: self
    type (varying_string    ), intent(in   ), optional :: groupName
    !$GLC attributes unused :: self, groupName

    return
  end subroutine nullFinalize

  subroutine nullReduce(self,reduced)
    !!{
    Implement a null output analysis reduction.
    !!}
    implicit none
    class(outputAnalysisNull ), intent(inout) :: self
    class(outputAnalysisClass), intent(inout) :: reduced
    !$GLC attributes unused :: self, reduced

    return
  end subroutine nullReduce

  double precision function nullLogLikelihood(self)
    !!{
    Return the log-likelihood of a null output analysis.
    !!}
    implicit none
    class(outputAnalysisNull), intent(inout) :: self
    !$GLC attributes unused :: self

    nullLogLikelihood=0.0d0
    return
  end function nullLogLikelihood

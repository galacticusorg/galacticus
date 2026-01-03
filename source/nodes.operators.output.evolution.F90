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
  Implements a node operator class that outputs complete data on node evolution.
  !!}

  use :: Galactic_Filters, only : galacticFilterClass
  
  !![
  <nodeOperator name="nodeOperatorEvolutionOutput">
   <description>A node operator class that outputs data on mergers between galaxies.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorEvolutionOutput
     !!{
     A node operator class that shifts node indices at node promotion.
     !!}
     private
     type   (varying_string     )          :: outputFileName
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
     integer                               :: outputFile
   contains
     final     ::                              evolutionOutputDestructor
     procedure :: differentialEvolutionPost => evolutionOutputDifferentialEvolutionPost
  end type nodeOperatorEvolutionOutput
  
  interface nodeOperatorEvolutionOutput
     !!{
     Constructors for the \refClass{nodeOperatorEvolutionOutput} node operator class.
     !!}
     module procedure evolutionOutputConstructorParameters
     module procedure evolutionOutputConstructorInternal
  end interface nodeOperatorEvolutionOutput
  
contains
  
  function evolutionOutputConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorEvolutionOutput} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorEvolutionOutput)                :: self
    type (inputParameters            ), intent(inout) :: parameters
    class(galacticFilterClass        ), pointer       :: galacticFilter_
    type (varying_string             )                :: outputFileName

    !![
    <inputParameter>
      <name>outputFileName</name>
      <defaultValue>var_str('mergerTreeEvolution.xml')</defaultValue>
      <description>The name of the file to which merger tree evolution should be output.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="galacticFilter" parameterName="galacticFilter" name="galacticFilter_" source="parameters">
     <default>
      <galacticFilter value="always"/>
     </default>
    </objectBuilder>
    !!]
    self=nodeOperatorEvolutionOutput(outputFileName,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function evolutionOutputConstructorParameters

  function evolutionOutputConstructorInternal(outputFileName,galacticFilter_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorEvolutionOutput} node operator class.
    !!}
    implicit none
    type (nodeOperatorEvolutionOutput)                        :: self
    class(galacticFilterClass        ), intent(in   ), target :: galacticFilter_
    type (varying_string             ), intent(in   )         :: outputFileName
    !![
    <constructorAssign variables="outputFileName, *galacticFilter_"/>
    !!]
    
    open(newUnit=self%outputFile,file=char(outputFileName),status='unknown',form='formatted')
    write (self%outputFile,'(a)') '<evolution>'
    return
  end function evolutionOutputConstructorInternal
  
  subroutine evolutionOutputDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorEvolutionOutput} node operator class.
    !!}
    implicit none
    type(nodeOperatorEvolutionOutput), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    write (self%outputFile,'(a)') '</evolution>'
    close(self%outputFile)
    return
  end subroutine evolutionOutputDestructor
  
  subroutine evolutionOutputDifferentialEvolutionPost(self,node)
    !!{
    Operate on the node after differential evolution
    !!}
    implicit none
    class(nodeOperatorEvolutionOutput), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node

    if (self%galacticFilter_%passes(node)) call node%serializeXML(self%outputFile)    
    return
  end subroutine evolutionOutputDifferentialEvolutionPost
  

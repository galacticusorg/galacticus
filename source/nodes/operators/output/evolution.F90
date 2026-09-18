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
  Implements a node operator class that outputs complete data on node evolution.
  !!}

  use :: Galactic_Filters, only : galacticFilterClass
  
  !![
  <nodeOperator name="nodeOperatorEvolutionOutput" docformat="rst">
   <description>
   A node operator class that records the complete evolutionary trajectory of each node to an XML file at each ODE post-step, for debugging and detailed analysis. ``outputFileName`` specifies the output file (default: ``mergerTreeEvolution.xml``); a ``galacticFilter`` selects which nodes to record (default: all nodes). Useful for inspecting how galaxy properties evolve between output snapshots.
   </description>
   <deepCopy>
    <setTo variables="isOwner" value=".false."/>
   </deepCopy>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorEvolutionOutput
     !!{RST
     A node operator class that outputs complete data on node evolution.
     !!}
     private
     type   (varying_string     )          :: outputFileName
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
     logical                               :: isOwner         =  .false.
   contains
     final     ::                              evolutionOutputDestructor
     procedure :: differentialEvolutionPost => evolutionOutputDifferentialEvolutionPost
  end type nodeOperatorEvolutionOutput
  
  interface nodeOperatorEvolutionOutput
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorEvolutionOutput` node operator class.
     !!}
     module procedure evolutionOutputConstructorParameters
     module procedure evolutionOutputConstructorInternal
  end interface nodeOperatorEvolutionOutput

  ! The output file. This is shared by all instances of the class - including the copies made for each OpenMP thread, which
  ! must all write to the same file - so it is opened on first use (rather than by the constructor), and closed only by the
  ! instance which was constructed (and so "owns" the file). This also ensures that finalization of any temporary or
  ! default-initialized instances (which occurs during construction) can not close the file. Whether the file is open is recorded
  ! separately, as unit numbers returned by `open(newUnit=...)` are negative and so can not serve as an indicator.
  integer                 :: evolutionOutputUnit
  logical                 :: evolutionOutputIsOpen  =.false.
  type   (varying_string) :: evolutionOutputFileName
  
contains
  
  function evolutionOutputConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorEvolutionOutput` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorEvolutionOutput)                :: self
    type (inputParameters            ), intent(inout) :: parameters
    class(galacticFilterClass        ), pointer       :: galacticFilter_
    type (varying_string             )                :: outputFileName

    !![
    <inputParameter docformat="rst">
      <name>outputFileName</name>
      <defaultValue>var_str('mergerTreeEvolution.xml')</defaultValue>
      <description>
      The name of the file to which merger tree evolution should be output.
      </description>
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
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorEvolutionOutput` node operator class.
    !!}
    implicit none
    type (nodeOperatorEvolutionOutput)                        :: self
    class(galacticFilterClass        ), intent(in   ), target :: galacticFilter_
    type (varying_string             ), intent(in   )         :: outputFileName
    !![
    <constructorAssign variables="outputFileName, *galacticFilter_"/>
    !!]

    self%isOwner=.true.
    return
  end function evolutionOutputConstructorInternal
  
  subroutine evolutionOutputDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorEvolutionOutput` node operator class.
    !!}
    implicit none
    type(nodeOperatorEvolutionOutput), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    if (self%isOwner) then
       !$omp critical(evolutionOutputFile)
       if (evolutionOutputIsOpen) then
          write (evolutionOutputUnit,'(a)') '</evolution>'
          close(evolutionOutputUnit)
          evolutionOutputIsOpen=.false.
       end if
       !$omp end critical(evolutionOutputFile)
    end if
    return
  end subroutine evolutionOutputDestructor
  
  subroutine evolutionOutputDifferentialEvolutionPost(self,node)
    !!{RST
    Operate on the node after differential evolution
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : operator(/=)
    implicit none
    class(nodeOperatorEvolutionOutput), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node

    if (.not.self%galacticFilter_%passes(node)) return
    !$omp critical(evolutionOutputFile)
    if (.not.evolutionOutputIsOpen) then
       open(newUnit=evolutionOutputUnit,file=char(self%outputFileName),status='replace',form='formatted')
       write (evolutionOutputUnit,'(a)') '<evolution>'
       evolutionOutputIsOpen  =.true.
       evolutionOutputFileName=self%outputFileName
    else if (self%outputFileName /= evolutionOutputFileName) then
       call Error_Report('only a single evolution output file is supported'//{introspection:location})
    end if
    call node%serializeXML(evolutionOutputUnit)
    !$omp end critical(evolutionOutputFile)
    return
  end subroutine evolutionOutputDifferentialEvolutionPost
  

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
  Implements a merger tree outputter class that outputs the full state of merger trees for later postprocessing.
  !!}

  !![
  <mergerTreeOutputter name="mergerTreeOutputterFullState">
   <description>
    A merger tree outputter class that outputs the full state of merger trees to allow later postprocessing. Complete tree data in
    output in raw binary format to the file specified as {\normalfont \ttfamily [fileName]}. This can be later re-read and
    post-processed using the \refClass{taskPostprocessForests} class.
   </description>
  </mergerTreeOutputter>
  !!]
  type, extends(mergerTreeOutputterClass) :: mergerTreeOutputterFullState
     !!{
     Implementation of a merger tree outputter class that outputs the full state of merger trees to allow later postprocessing.
     !!}
     private
     type(varying_string) :: fileName
   contains
     procedure :: outputTree => fullStateOutputTree
     procedure :: outputNode => fullStateOutputNode
     procedure :: finalize   => fullStateFinalize
  end type mergerTreeOutputterFullState

  interface mergerTreeOutputterFullState
     !!{
     Constructors for the \refClass{mergerTreeOutputterFullState} merger tree outputter.
     !!}
     module procedure fullStateConstructorParameters
     module procedure fullStateConstructorInternal
  end interface mergerTreeOutputterFullState
  
contains
  
  function fullStateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOutputterFullState} merger tree outputter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeOutputterFullState)                :: self
    type(inputParameters             ), intent(inout) :: parameters
    type(varying_string              )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file to which state should be stored.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeOutputterFullState(fileName)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fullStateConstructorParameters

  function fullStateConstructorInternal(fileName) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeOutputterFullState} merger tree outputter class.
    !!}
    implicit none
    type(mergerTreeOutputterFullState)                :: self
    type(varying_string              ), intent(in   ) :: fileName
    !![
    <constructorAssign variables="fileName"/>
    !!]

    return
  end function fullStateConstructorInternal
  
  subroutine fullStateOutputTree(self,tree,indexOutput,time)
    !!{
    Write properties of nodes in {\normalfont \ttfamily tree} to the \glc\ output file.
    !!}
    use :: File_Utilities          , only : File_Lock           , File_Unlock, lockDescriptor
    use :: Merger_Tree_Construction, only : mergerTreeStateStore
    implicit none
    class           (mergerTreeOutputterFullState), intent(inout)         :: self
    type            (mergerTree                  ), intent(inout), target :: tree
    integer         (c_size_t                    ), intent(in   )         :: indexOutput
    double precision                              , intent(in   )         :: time
    type            (lockDescriptor              )                        :: fileLock
    !$GLC attributes unused :: time
    
    call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
    call mergerTreeStateStore(tree,char(self%fileName),indexOutput=indexOutput,snapshot=.false.,append=.true.)
    call File_Unlock(fileLock)
    return
  end subroutine fullStateOutputTree

  subroutine fullStateOutputNode(self,node,indexOutput)
    !!{
    Perform no output.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (mergerTreeOutputterFullState), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer(c_size_t                    ), intent(in   ) :: indexOutput
    !$GLC attributes unused :: self, node, indexOutput

    call Error_Report('output of single nodes is not supported'//{introspection:location})
    return
  end subroutine fullStateOutputNode

  subroutine fullStateFinalize(self)
    !!{
    Finalize full state output---nothing to do here.
    !!}
    implicit none
    class(mergerTreeOutputterFullState), intent(inout) :: self

    return
  end subroutine fullStateFinalize
  

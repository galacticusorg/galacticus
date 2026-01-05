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
Contains a module which implements a class for merger tree evolution timestepping.
!!}

module Merger_Tree_Timesteps
  !!{
  Implements a class for merger tree evolution timestepping.
  !!}
  use :: Galacticus_Nodes                   , only : mergerTree                   , treeNode
  use :: Merger_Trees_Evolve_Deadlock_Status, only : enumerationDeadlockStatusType
  implicit none
  private
  public :: timestepTask

  !![
  <functionClass>
   <name>mergerTreeEvolveTimestep</name>
   <descriptiveName>Merger Tree Evolution Timesteps</descriptiveName>
   <description>Class providing timestep control for merger tree evolution.</description>
   <default>standard</default>
   <method name="timeEvolveTo">
    <description>
     Return the time to which the {\normalfont \ttfamily node} can be evolved. The current limiting time is provided as
     {\normalfont \ttfamily timeEnd}. Optionally, the procedure pointer {\normalfont \ttfamily task} can be set to point to a
     subroutine which will be called after the node is evolved to the end of the timestep. It is acceptable for this pointer to be
     null. The {\normalfont \ttfamily taskSelf} pointer may be set to point to the timestep object and will be made available to
     the timestep task subroutine. Note that the {\normalfont \ttfamily task} will only be called for the task which provided the
     shortest timestep---other tasks can always request to be called again when the next timestep is determined. The subroutine to
     be called at the end of the timestep must have the form:
     \begin{verbatim}
       subroutine timestepTask(self,tree,node,deadlockStatus)
         implicit none
         class(*                            ), intent(inout)          :: self
         type (mergerTree                   ), intent(inout)          :: tree
         type (treeNode                     ), intent(inout), pointer :: node
         type (enumerationDeadlockStatusType), intent(inout)          :: deadlockStatus
         .
         .
         .
         return
       end subroutine timestepTask
     \end{verbatim}
  
     The {\normalfont \ttfamily deadlockStatus} argument should be set to {\normalfont \ttfamily isNotDeadlocked} (provided by
     the
     \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Source.pdf#source.merger_trees_evolve_deadlock_options_F90:merger_trees_evolve_deadlock_status}{\normalfont
     \ttfamily Merger\_Trees\_Evolve\_Deadlock\_Status} module) if, and only if, the end of timestep task makes some change to
     the state of the tree (e.g. merging a node), to indicate that the tree was not deadlocked in this pass (i.e. something
     actually changed in the tree).
  
     If the {\normalfont \ttfamily report} argument is {\normalfont \ttfamily true} then the function should report the value of
     {\normalfont \ttfamily timestep} prior to exiting. (This is used in reporting on timestepping criteri in deadlocked trees.)
     It is recommended that the report be made using the
     \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Source.pdf#source.merger_trees_evolve_timesteps_report_F90:evolve_to_time_reports:evolve_to_time_report}{\normalfont
     \ttfamily Evolve\_To\_Time\_Report()} function. Additionally, if the optional {\normalfont \ttfamily lockNode} and
     {\normalfont \ttfamily lockType} arguments are present then additional information can be supplied to aid in diagnosing
     deadlock conditions. If the current task is limiting the timestep then the {\normalfont \ttfamily lockNode} pointer should
     be set to point to whichever node is causing the limit (which may be {\normalfont \ttfamily node} or some other node,
     e.g. a satellite of {\normalfont \ttfamily node}, etc.), and {\normalfont \ttfamily lockType} should be set to a short
     description label identifying the type of limit.
    </description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision                , intent(in   )                    :: timeEnd </argument>
    <argument>type            (treeNode      ), intent(inout)          , target  :: node    </argument>
    <argument>procedure       (timestepTask  ), intent(  out)          , pointer :: task    </argument>
    <argument>class           (*             ), intent(  out)          , pointer :: taskSelf</argument>
    <argument>logical                         , intent(in   )                    :: report  </argument>
    <argument>type            (treeNode      ), intent(  out), optional, pointer :: lockNode</argument>
    <argument>type            (varying_string), intent(  out), optional          :: lockType</argument>
   </method>
   <method name="refuseToEvolve">
    <description>Return true if evolution should be refused.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
    <code>
      !$GLC attributes unused :: self, node
      mergerTreeEvolveTimestepRefuseToEvolve=.false.
    </code>
   </method>
  </functionClass>
  !!]

  abstract interface
     subroutine timestepTask(self,tree,node,deadlockStatus)
       import mergerTree, treeNode, enumerationDeadlockStatusType
       class(*                            ), intent(inout)          :: self
       type (mergerTree                   ), intent(in   )          :: tree
       type (treeNode                     ), intent(inout), pointer :: node
       type (enumerationDeadlockStatusType), intent(inout)          :: deadlockStatus
     end subroutine timestepTask
  end interface

end module Merger_Tree_Timesteps

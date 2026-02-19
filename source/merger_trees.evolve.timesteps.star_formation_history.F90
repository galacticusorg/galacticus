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

  use :: Star_Formation_Histories, only : starFormationHistoryClass

  !![
  <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepStarFormationHistory">
   <description>
     A merger tree evolution timestepping class that limits the timestep to the next bin in the star formation history.
   </description>
  </mergerTreeEvolveTimestep>
  !!]
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepStarFormationHistory
     !!{
     A merger tree evolution timestepping class that limits the timestep to the next bin in the star formation history.
     !!}
     private
     class(starFormationHistoryClass), pointer :: starFormationHistory_ => null()
   contains
     procedure :: timeEvolveTo => starFormationHistoryTimeEvolveTo
  end type mergerTreeEvolveTimestepStarFormationHistory

  interface mergerTreeEvolveTimestepStarFormationHistory
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepStarFormationHistory} merger tree evolution timestep class.
     !!}
     module procedure starFormationHistoryConstructorParameters
     module procedure starFormationHistoryConstructorInternal
  end interface mergerTreeEvolveTimestepStarFormationHistory

contains

  function starFormationHistoryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepStarFormationHistory} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeEvolveTimestepStarFormationHistory)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(starFormationHistoryClass                   ), pointer       :: starFormationHistory_

    !![
    <objectBuilder class="starFormationHistory" name="starFormationHistory_" source="parameters"/>
    !!]
    self=mergerTreeEvolveTimestepStarFormationHistory(starFormationHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"/>
    !!]
    return
  end function starFormationHistoryConstructorParameters

  function starFormationHistoryConstructorInternal(starFormationHistory_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeEvolveTimestepStarFormationHistory} merger tree timestepping class.
    !!}
    implicit none
    type (mergerTreeEvolveTimestepStarFormationHistory)                        :: self
    class(starFormationHistoryClass                  ), intent(in   ), target :: starFormationHistory_
    !![
    <constructorAssign variables="*starFormationHistory_"/>
    !!]

    return
  end function starFormationHistoryConstructorInternal

  subroutine starFormationHistoryDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolveTimestepStarFormationHistory} merger tree timestepping class.
    !!}
    implicit none
    type(mergerTreeEvolveTimestepStarFormationHistory), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationHistory_"/>
    !!]
    return
  end subroutine starFormationHistoryDestructor

  double precision function starFormationHistoryTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType) result(timeEvolveTo)
    !!{
    Determine a suitable timestep for {\normalfont \ttfamily node} using the starFormationHistory method. This simply selects the smaller of {\normalfont \ttfamily
    timeStepAbsolute} and {\normalfont \ttfamily timeStepRelative}$H^{-1}(t)$.
    !!}
    use :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use :: Galacticus_Nodes      , only : nodeComponentBasic   , nodeComponentDisk, nodeComponentSpheroid
    use :: ISO_Varying_String    , only : varying_string
    use :: Histories             , only : history
    implicit none
    class           (mergerTreeEvolveTimestepStarFormationHistory), intent(inout), target            :: self
    double precision                                              , intent(in   )                    :: timeEnd
    type            (treeNode                                    ), intent(inout), target            :: node
    procedure       (timestepTask                                ), intent(  out), pointer           :: task
    class           (*                                           ), intent(  out), pointer           :: taskSelf
    logical                                                       , intent(in   )                    :: report
    type            (treeNode                                    ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                              ), intent(  out)         , optional :: lockType
    class           (nodeComponentBasic                          )               , pointer           :: basic
    class           (nodeComponentDisk                           )               , pointer           :: disk
    class           (nodeComponentSpheroid                       )               , pointer           :: spheroid
    type            (history                                     )                                   :: historyStarFormationDisk    , historyStarFormationSpheroid
    double precision                                                                                 :: timeStarFormationHistoryNext
    !$GLC attributes unused :: timeEnd

    basic                        => node    %basic               ()
    disk                         => node    %disk                ()
    spheroid                     => node    %spheroid            ()
    historyStarFormationDisk     =  disk    %starFormationHistory()
    historyStarFormationSpheroid =  spheroid%starFormationHistory()
    timeEvolveTo                 =  huge(0.0d0)
    if (historyStarFormationDisk    %exists()) then
       timeStarFormationHistoryNext=self%starFormationHistory_%timeNext(node,historyStarFormationDisk)
       if (basic%time() < timeStarFormationHistoryNext) &
            & timeEvolveTo=min(timeEvolveTo,timeStarFormationHistoryNext)
    end if
    if (historyStarFormationSpheroid%exists()) then
       timeStarFormationHistoryNext=self%starFormationHistory_%timeNext(node,historyStarFormationSpheroid)
       if (basic%time() < timeStarFormationHistoryNext) &
            & timeEvolveTo=min(timeEvolveTo,timeStarFormationHistoryNext)
    end if
    task                            => null()
    taskSelf                        => null()
    if (present(lockNode)) lockNode => node
    if (present(lockType)) lockType =  "starFormationHistory"
    if (        report   ) call Evolve_To_Time_Report("starFormationHistory: ",timeEvolveTo)
    return
  end function starFormationHistoryTimeEvolveTo

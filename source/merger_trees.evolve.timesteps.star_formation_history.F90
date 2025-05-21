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
   contains
     procedure :: timeEvolveTo => starFormationHistoryTimeEvolveTo
  end type mergerTreeEvolveTimestepStarFormationHistory

  interface mergerTreeEvolveTimestepStarFormationHistory
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepStarFormationHistory} merger tree evolution timestep class.
     !!}
     module procedure starFormationHistoryConstructorParameters
  end interface mergerTreeEvolveTimestepStarFormationHistory

contains

  function starFormationHistoryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepStarFormationHistory} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(mergerTreeEvolveTimestepStarFormationHistory)                :: self
    type(inputParameters                             ), intent(inout) :: parameters

    self=mergerTreeEvolveTimestepStarFormationHistory()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function starFormationHistoryConstructorParameters

  double precision function starFormationHistoryTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType) result(timeEvolveTo)
    !!{
    Determine a suitable timestep for {\normalfont \ttfamily node} using the starFormationHistory method. This simply selects the smaller of {\normalfont \ttfamily
    timeStepAbsolute} and {\normalfont \ttfamily timeStepRelative}$H^{-1}(t)$.
    !!}
    use, intrinsic :: ISO_C_Binding         , only : c_size_t
    use            :: Arrays_Search         , only : searchArray
    use            :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use            :: Galacticus_Nodes      , only : nodeComponentBasic   , nodeComponentDisk, nodeComponentSpheroid
    use            :: ISO_Varying_String    , only : varying_string
    use            :: Histories             , only : history
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
    type            (history                                     )                                   :: historyStarFormationDisk, historyStarFormationSpheroid
    integer         (c_size_t                                    )                                   :: i
    !$GLC attributes unused :: timeEnd


    basic                        => node    %basic               ()
    disk                         => node    %disk                ()
    spheroid                     => node    %spheroid            ()
    historyStarFormationDisk     =  disk    %starFormationHistory()
    historyStarFormationSpheroid =  spheroid%starFormationHistory()
    timeEvolveTo                 =  huge(0.0d0)
    if (historyStarFormationDisk    %exists()) then
       i=searchArray(historyStarFormationDisk%time,basic%time())+1
       if (basic%time() < historyStarFormationDisk    %time(i))                   &
            & timeEvolveTo=min(timeEvolveTo,historyStarFormationDisk    %time(i))       
    end if
    if (historyStarFormationSpheroid%exists()) then
       i=searchArray(historyStarFormationSpheroid%time,basic%time())+1
       if (basic%time() < historyStarFormationSpheroid%time(i))                   &
            & timeEvolveTo=min(timeEvolveTo,historyStarFormationSpheroid%time(i))       
    end if
    task                            => null()
    taskSelf                        => null()
    if (present(lockNode)) lockNode => node
    if (present(lockType)) lockType =  "starFormationHistory"
    if (        report   ) call Evolve_To_Time_Report("starFormationHistory: ",timeEvolveTo)
    return
  end function starFormationHistoryTimeEvolveTo

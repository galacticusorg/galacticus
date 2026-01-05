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
  Contains a merger tree evolution timestep class which limits the step to the next lightcone crossing.
  !!}

  use :: Geometry_Lightcones   , only : geometryLightconeClass
  use :: Merger_Tree_Outputters, only : mergerTreeOutputterClass

  !![
  <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepLightconeCrossing">
   <description>  
    A merger tree evolution timestepping class which limits the step to the next lightcone crossing.
   </description>
  </mergerTreeEvolveTimestep>
  !!]
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepLightconeCrossing
     !!{
     Implementation of a merger tree evolution timestep class which limits the step to the next lightcone crossing.
     !!}
     private
     class           (geometryLightconeClass  ), pointer :: geometryLightcone_   => null()
     class           (mergerTreeOutputterClass), pointer :: mergerTreeOutputter_ => null()
     double precision                                    :: timeMinimum                   , timeMaximum
     integer                                             :: timesCrossingID
    contains
     final     ::                 lightconeCrossingDestructor
     procedure :: timeEvolveTo => lightconeCrossingTimeEvolveTo
  end type mergerTreeEvolveTimestepLightconeCrossing

  interface mergerTreeEvolveTimestepLightconeCrossing
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepLightconeCrossing} merger tree evolution timestep class.
     !!}
     module procedure lightconeCrossingConstructorParameters
     module procedure lightconeCrossingConstructorInternal
  end interface mergerTreeEvolveTimestepLightconeCrossing

contains
  
  function lightconeCrossingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepLightconeCrossing} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeEvolveTimestepLightconeCrossing)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (geometryLightconeClass                   ), pointer       :: geometryLightcone_
    class           (mergerTreeOutputterClass                 ), pointer       :: mergerTreeOutputter_
    
    !![
    <objectBuilder class="geometryLightcone"   name="geometryLightcone_"   source="parameters"/>
    <objectBuilder class="mergerTreeOutputter" name="mergerTreeOutputter_" source="parameters"/>
    !!]
    self=mergerTreeEvolveTimestepLightconeCrossing(geometryLightcone_,mergerTreeOutputter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="geometryLightcone_"  />
    <objectDestructor name="mergerTreeOutputter_"/>
    !!]
    return
  end function lightconeCrossingConstructorParameters

  function lightconeCrossingConstructorInternal(geometryLightcone_,mergerTreeOutputter_) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepLightconeCrossing} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    implicit none
    type            (mergerTreeEvolveTimestepLightconeCrossing)                        :: self
    class           (geometryLightconeClass                   ), intent(in   ), target :: geometryLightcone_
    class           (mergerTreeOutputterClass                 ), intent(in   ), target :: mergerTreeOutputter_
    !![
    <constructorAssign variables="*geometryLightcone_, *mergerTreeOutputter_"/>
    !!]

    self%timeMinimum=self%geometryLightcone_%timeMinimum()
    self%timeMaximum=self%geometryLightcone_%timeMaximum()
    !![
    <addMetaProperty component="basic" name="timeStepLightconeCrossingTimes" id="self%timesCrossingID" rank="1" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function lightconeCrossingConstructorInternal

  subroutine lightconeCrossingDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolveTimestepLightconeCrossing} merger tree evolution timestep class.
    !!}
    implicit none
    type(mergerTreeEvolveTimestepLightconeCrossing), intent(inout) :: self

    !![
    <objectDestructor name="self%geometryLightcone_"  />
    <objectDestructor name="self%mergerTreeOutputter_"/>
    !!]
    return
  end subroutine lightconeCrossingDestructor

  double precision function lightconeCrossingTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType)
    !!{
    Determine a suitable timestep for {\normalfont \ttfamily node} such that it does not exceed the time of the next lightconeCrossing merger.
    !!}
    use :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use :: Galacticus_Nodes      , only : nodeComponentBasic
    implicit none
    class           (mergerTreeEvolveTimestepLightconeCrossing), intent(inout), target                 :: self
    double precision                                           , intent(in   )                         :: timeEnd
    type            (treeNode                                 ), intent(inout), target                 :: node
    procedure       (timestepTask                             ), intent(  out), pointer                :: task
    class           (*                                        ), intent(  out), pointer                :: taskSelf
    logical                                                    , intent(in   )                         :: report
    type            (treeNode                                 ), intent(  out), pointer     , optional :: lockNode
    type            (varying_string                           ), intent(  out)              , optional :: lockType
    class           (nodeComponentBasic                       )               , pointer                :: basic
    double precision                                                                                   :: timeCrossing
    double precision                                           , allocatable  , dimension(:)           :: timesCrossing, timesCrossing_

    ! Find the next crossing time for this node.
    lightconeCrossingTimeEvolveTo=huge(0.0d0)
    ! Consider only times after the earliest time specified.
    basic => node%basic()
    if (timeEnd >= max(self%timeMinimum,basic%time())) then
       ! Find the current set of crossing times associated with this node.
       timesCrossing=basic%floatRank1MetaPropertyGet(self%timesCrossingID)
       ! If crossing times have not yet been determined, find them now.
       if (size(timesCrossing) == 0) then
          ! Find the time (if any) of lightcone crossing.
          deallocate(timesCrossing)
          timeCrossing=self%geometryLightcone_%timeLightconeCrossing(node,max(self%timeMinimum,basic%time()),self%timeMaximum,timesCrossing)
          allocate(timesCrossing_(size(timesCrossing)+1))
          timesCrossing_(1:size(timesCrossing ))=timesCrossing
          timesCrossing_(  size(timesCrossing_))=huge(0.0d0)
          call basic%floatRank1MetaPropertySet(self%timesCrossingID,timesCrossing_)
          deallocate(timesCrossing)
          call move_alloc(timesCrossing_,timesCrossing)
       end if
       ! Set the crossing time to the next time.
       if (timesCrossing(1) <= timeEnd) lightconeCrossingTimeEvolveTo=timesCrossing(1)
       ! If a crossing occurs, set the relevant task.
       if (lightconeCrossingTimeEvolveTo < huge(0.0d0)) then
          task                            => lightconeCrossingProcess
          taskSelf                        => self
          if (present(lockNode)) lockNode => node
          if (present(lockType)) lockType =  "lightcone crossing"
          if (        report   ) call Evolve_To_Time_Report("lightcone crossing: ",lightconeCrossingTimeEvolveTo)
       else
          task                            => null()
          taskSelf                        => null()
          if (present(lockNode)) lockNode => null()
          if (present(lockType)) lockType =  ""
       end if
    end if
    return
  end function lightconeCrossingTimeEvolveTo
  
  subroutine lightconeCrossingProcess(self,tree,node,deadlockStatus)
    !!{
    Process a lightconeCrossing node which has undergone a merger with its host node.
    !!}
    use :: Error                              , only : Error_Report
    use :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsNotDeadlocked
    use :: Galacticus_Nodes                   , only : nodeComponentBasic
    implicit none
    class           (*                            ), intent(inout)               :: self
    type            (mergerTree                   ), intent(in   )               :: tree
    type            (treeNode                     ), intent(inout), pointer      :: node
    type            (enumerationDeadlockStatusType), intent(inout)               :: deadlockStatus
    class           (nodeComponentBasic           )               , pointer      :: basic
    double precision                               , allocatable  , dimension(:) :: timesCrossing
    !$GLC attributes unused :: tree

    select type (self)
    class is (mergerTreeEvolveTimestepLightconeCrossing)
       basic         => node %basic                    (                    )     
       timesCrossing =  basic%floatRank1MetaPropertyGet(self%timesCrossingID)
       if (basic%time() /= timesCrossing(1)) return
       call self%mergerTreeOutputter_%outputNode(node,1_c_size_t)
       ! Remove the crossing time for this node.
       timesCrossing=basic%floatRank1MetaPropertyGet(self%timesCrossingID)
       call basic%floatRank1MetaPropertySet(self%timesCrossingID,timesCrossing(2:size(timesCrossing)))
       ! The tree was changed, so mark that it is not deadlocked.
       deadlockStatus=deadlockStatusIsNotDeadlocked
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine lightconeCrossingProcess

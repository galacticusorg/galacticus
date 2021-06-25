!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
     double precision                                    :: timeMinimum
    contains
     final     ::                 lightconeCrossingDestructor
     procedure :: timeEvolveTo => lightconeCrossingTimeEvolveTo
  end type mergerTreeEvolveTimestepLightconeCrossing

  interface mergerTreeEvolveTimestepLightconeCrossing
     !!{
     Constructors for the {\normalfont \ttfamily lightconeCrossing} merger tree evolution timestep class.
     !!}
     module procedure lightconeCrossingConstructorParameters
     module procedure lightconeCrossingConstructorInternal
  end interface mergerTreeEvolveTimestepLightconeCrossing

contains
  
  function lightconeCrossingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily lightconeCrossing} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctionsClass
    use :: Input_Parameters   , only : inputParameter         , inputParameters
    implicit none
    type            (mergerTreeEvolveTimestepLightconeCrossing)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (geometryLightconeClass                   ), pointer       :: geometryLightcone_
    class           (mergerTreeOutputterClass                 ), pointer       :: mergerTreeOutputter_
    class           (cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    double precision                                                           :: redshiftMaximum
    
    !![
    <inputParameter>
      <name>redshiftMaximum</name>
      <description>The maximum redshift at which to limit timesteps to the next lightcone crossing.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="geometryLightcone"   name="geometryLightcone_"   source="parameters"/>
    <objectBuilder class="mergerTreeOutputter" name="mergerTreeOutputter_" source="parameters"/>
    !!]
    self=mergerTreeEvolveTimestepLightconeCrossing(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum)),geometryLightcone_,mergerTreeOutputter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="geometryLightcone_" />
    !!]
    return
  end function lightconeCrossingConstructorParameters

  function lightconeCrossingConstructorInternal(timeMinimum,geometryLightcone_,mergerTreeOutputter_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily lightconeCrossing} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    implicit none
    type            (mergerTreeEvolveTimestepLightconeCrossing)                        :: self
    class           (geometryLightconeClass                   ), intent(in   ), target :: geometryLightcone_
    class           (mergerTreeOutputterClass                 ), intent(in   ), target :: mergerTreeOutputter_
    double precision                                           , intent(in   )         :: timeMinimum
    !![
    <constructorAssign variables="timeMinimum, *geometryLightcone_, *mergerTreeOutputter_"/>
    !!]
  
    return
  end function lightconeCrossingConstructorInternal

  subroutine lightconeCrossingDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily lightconeCrossing} erger tree evolution timestep class.
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
    use :: Galacticus_Nodes      , only : nodeComponentPosition, nodeComponentBasic
    implicit none
    class           (mergerTreeEvolveTimestepLightconeCrossing), intent(inout), target            :: self
    double precision                                           , intent(in   )                    :: timeEnd
    type            (treeNode                                 ), intent(inout), target            :: node
    procedure       (timestepTask                             ), intent(  out), pointer           :: task
    class           (*                                        ), intent(  out), pointer           :: taskSelf
    logical                                                    , intent(in   )                    :: report
    type            (treeNode                                 ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                           ), intent(  out)         , optional :: lockType
    class           (nodeComponentPosition                    )               , pointer           :: position         , positionParent
    class           (nodeComponentBasic                       )               , pointer           :: basic
    double precision                                                                              :: timeCrossing     , timeMaximum   , &
         &                                                                                           timeMaximumParent
    
    ! Find the next crossing time for this node.
    lightconeCrossingTimeEvolveTo=huge(0.0d0)
    ! Consider only times after the earliest time specified.
    if (timeEnd >= self%timeMinimum) then
       ! Limit the maximum time to the latest time for which the position of the node is known.
       position    => node    %position                ()
       timeMaximum =  position%interpolationTimeMaximum()
       if (timeMaximum > 0.0d0) then
          timeMaximum=min(timeMaximum,timeEnd)
       else
          timeMaximum=                timeEnd
       end if
       ! For satellite nodes also limit the maximum time to the latest time for which the position of the parent node is known.
       if (node%isSatellite()) then
          positionParent    => node          %parent%position                ()
          timeMaximumParent =  positionParent       %interpolationTimeMaximum()
          if (timeMaximumParent > 0.0d0) timeMaximum=min(timeMaximum,timeMaximumParent)
       end if
       ! If the maximum time is after the current time, find the time (if any) of lightcone crossing.
       basic => node%basic()
       if (timeMaximum > basic%time()) then       
          timeCrossing=self%geometryLightcone_%timeLightconeCrossing(node,timeMaximum)
          if (timeCrossing <= timeEnd) lightconeCrossingTimeEvolveTo=timeCrossing
       end if
    end if
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
    return
  end function lightconeCrossingTimeEvolveTo
  
  subroutine lightconeCrossingProcess(self,tree,node,deadlockStatus)
    !!{
    Process a lightconeCrossing node which has undergone a merger with its host node.
    !!}
    use :: Galacticus_Error                   , only : Galacticus_Error_Report
    use :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsNotDeadlocked
    use mpi_utilities
    implicit none
    class  (*         ), intent(inout)          :: self
    type   (mergerTree), intent(in   )          :: tree
    type   (treeNode  ), intent(inout), pointer :: node
    integer            , intent(inout)          :: deadlockStatus
    !$GLC attributes unused :: tree

    select type (self)
    class is (mergerTreeEvolveTimestepLightconeCrossing)
       call self%mergerTreeOutputter_%outputNode(node,1_c_size_t)
       ! The tree was changed, so mark that it is not deadlocked.
       deadlockStatus=deadlockStatusIsNotDeadlocked
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine lightconeCrossingProcess

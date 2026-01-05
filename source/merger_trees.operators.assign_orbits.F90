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
Implements a merger tree operator which assigns orbits to non-primary progenitor nodes.
!!}

  use :: Satellite_Merging_Timescales, only : satelliteMergingTimescales, satelliteMergingTimescalesClass
  use :: Virial_Orbits               , only : virialOrbit               , virialOrbitClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorAssignOrbits">
   <description>Provides a merger tree operator which assigns orbits to non-primary progenitor nodes.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorAssignOrbits
     !!{
     An orbit assigning merger tree operator class.
     !!}
     private
     class(virialOrbitClass               ), pointer :: virialOrbit_                => null()
     class(satelliteMergingTimescalesClass), pointer :: satelliteMergingTimescales_ => null()
   contains
     final     ::                        assignOrbitsDestructor
     procedure :: operatePreEvolution => assignOrbitsOperatePreEvolution
  end type mergerTreeOperatorAssignOrbits

  interface mergerTreeOperatorAssignOrbits
     !!{
     Constructors for the orbit assigning merger tree operator class.
     !!}
     module procedure assignOrbitsConstructorParameters
     module procedure assignOrbitsConstructorInternal
  end interface mergerTreeOperatorAssignOrbits

contains

  function assignOrbitsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the orbit assigning merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeOperatorAssignOrbits )                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(virialOrbitClass               ), pointer       :: virialOrbit_
    class(satelliteMergingTimescalesClass), pointer       :: satelliteMergingTimescales_

    !![
    <objectBuilder class="satelliteMergingTimescales" name="satelliteMergingTimescales_" source="parameters"/>
    <objectBuilder class="virialOrbit"                name="virialOrbit_"                source="parameters"/>
    !!]
    self=mergerTreeOperatorAssignOrbits(satelliteMergingTimescales_,virialOrbit_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteMergingTimescales_"/>
    <objectDestructor name="virialOrbit_"               />
    !!]
    return
  end function assignOrbitsConstructorParameters

  function assignOrbitsConstructorInternal(satelliteMergingTimescales_,virialOrbit_) result(self)
    !!{
    Constructor for the orbit assigning merger tree operator class which takes a parameter set as input.
    !!}
    implicit none
    type (mergerTreeOperatorAssignOrbits )                        :: self
    class(virialOrbitClass               ), intent(in   ), target :: virialOrbit_
    class(satelliteMergingTimescalesClass), intent(in   ), target :: satelliteMergingTimescales_
    !![
    <constructorAssign variables="*satelliteMergingTimescales_, *virialOrbit_"/>
    !!]

    return
  end function assignOrbitsConstructorInternal

  subroutine assignOrbitsDestructor(self)
    !!{
    Destructor for the merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorAssignOrbits), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteMergingTimescales_"/>
    <objectDestructor name="self%virialOrbit_"               />
    !!]
    return
  end subroutine assignOrbitsDestructor

  subroutine assignOrbitsOperatePreEvolution(self,tree)
    !!{
    Perform a orbit assigning operation on a merger tree.
    !!}
    use :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, nodeComponentSatellite, treeNode
    use :: Kepler_Orbits      , only : keplerOrbit
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    class  (mergerTreeOperatorAssignOrbits), intent(inout), target :: self
    type   (mergerTree                    ), intent(inout), target :: tree
    type   (treeNode                      ), pointer               :: node                    , nodeProgenitor       , &
         &                                                            mergee
    class  (nodeComponentSatellite        ), pointer               :: satellite               , satelliteProgenitor
    class  (nodeComponentBasic            ), pointer               :: basic                   , basicProgenitor
    type   (mergerTreeWalkerIsolatedNodes )                        :: treeWalker
    type   (keplerOrbit                   )                        :: virialOrbitNode         , virialOrbitProgenitor
    logical                                                        :: satelliteProgenitorFound
    !$GLC attributes unused :: self

    ! Iterate over nodes in forest.
    treeWalker=mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       ! An orbit must be assigned to non-primary progenitors which do not have an orbit already assigned.
       if     (                                             &
            &        associated(node%parent               ) &
            &  .and.                                        &
            &   .not.           node%isPrimaryProgenitor()  &
            & ) then
          satellite       => node     %satellite  (autoCreate=.true.)
          virialOrbitNode =  satellite%virialOrbit(                 )
          if (.not.virialOrbitNode%isDefined()) then
             ! Check for a primary progenitor with a pre-existing satellite.
             satelliteProgenitorFound =  .false.
             nodeProgenitor           => node
             do while (associated(nodeProgenitor))
                satelliteProgenitor   => nodeProgenitor     %satellite  ()
                virialOrbitProgenitor =  satelliteProgenitor%virialOrbit()
                if (virialOrbitProgenitor%isDefined()) then
                   ! A satellite exists in this progenitor.
                   satellite =  satelliteProgenitor
                   basic     => node               %basic()
                   if (satellite%timeOfMerging() < basic%time()) call satellite%timeOfMergingSet(basic%time())
                   if (associated(nodeProgenitor%mergeTarget)) then
                      mergee => nodeProgenitor%mergeTarget%firstMergee
                      if (associated(mergee,nodeProgenitor)) then
                         nodeProgenitor%mergeTarget%firstMergee => node
                      else
                         do while (.not.associated(mergee%siblingMergee,nodeProgenitor))
                            mergee => mergee%siblingMergee
                         end do
                         mergee%siblingMergee => node
                      end if
                      node          %siblingMergee => nodeProgenitor%siblingMergee
                      node          %mergeTarget   => nodeProgenitor%mergeTarget
                      nodeProgenitor%mergeTarget   => null()
                      nodeProgenitor%siblingMergee => null()
                      ! The merge target must be reachable at the merge time. If it is not, find a progenitor of it which is
                      ! reachable at that time.
                      if (associated(node%mergeTarget%firstChild)) then
                         nodeProgenitor  => node          %mergeTarget
                         basicProgenitor => nodeProgenitor%firstChild %basic()
                         if (basicProgenitor%time() > satellite%timeOfMerging()) then
                            ! Shift to an earlier progenitor as merge target.
                            do while (basicProgenitor%time() > satellite%timeOfMerging())
                               nodeProgenitor  => nodeProgenitor%firstChild
                               basicProgenitor => nodeProgenitor%firstChild%basic()
                            end do
                            ! Remove our mergee from its merge target.
                            call node%removeFromMergee()
                            if (associated(nodeProgenitor%firstMergee)) then
                               mergee => nodeProgenitor%firstMergee
                               do while (associated(mergee%siblingMergee))
                                  mergee => mergee%siblingMergee
                               end do
                               mergee%siblingMergee => node
                            else
                               nodeProgenitor%firstMergee => node
                            end if
                            node%mergeTarget   => nodeProgenitor
                            node%siblingMergee => null()
                         end if
                      end if
                      ! The mergee must be a satellite by the time it is due to merge. If it is not (which implies it is merging with a
                      ! halo in a different branch of the tree), reset it to be merging with the appropriate halo in its own branch of
                      ! the tree.
                      if (.not.node%isProgenitorOf(node%mergeTarget)) then
                         nodeProgenitor  => node          %parent
                         basicProgenitor => nodeProgenitor%basic ()
                         if (basicProgenitor%time() > satellite%timeOfMerging()) then
                            ! Adjust the time of merging.
                            call satellite%timeOfMergingSet(basicProgenitor%time())
                            ! Remove our mergee from its merge target.
                            call node%removeFromMergee()
                            if (associated(nodeProgenitor%firstMergee)) then
                               mergee => nodeProgenitor%firstMergee
                               do while (associated(mergee%siblingMergee))
                                  mergee => mergee%siblingMergee
                               end do
                               mergee%siblingMergee => node
                            else
                               nodeProgenitor%firstMergee => node
                            end if
                            node%mergeTarget   => nodeProgenitor
                            node%siblingMergee => null()
                         end if
                      end if
                   end if
                   satelliteProgenitorFound=.true.
                   exit
                end if
                nodeProgenitor => nodeProgenitor%firstChild
             end do
             if (.not.satelliteProgenitorFound) then
                basic           => node             %basic(                        )
                virialOrbitNode =  self%virialOrbit_%orbit(node,node%parent,.false.)
                call satellite%timeOfMergingSet(self%satelliteMergingTimescales_%timeUntilMerging(node,virialOrbitNode)+basic%time())
                call satellite%  virialOrbitSet(                                                       virialOrbitNode              )
             end if
          else
             ! The merge target must be reachable at the merge time. If it is not, find a
             ! progenitor of it which is reachable at that time.
             if (associated(node%mergeTarget).and.associated(node%mergeTarget%firstChild)) then
                nodeProgenitor  => node          %mergeTarget
                basicProgenitor => nodeProgenitor%firstChild %basic()
                if (basicProgenitor%time() > satellite%timeOfMerging()) then
                   ! Shift to an earlier progenitor as merge target.
                   do while (basicProgenitor%time() > satellite%timeOfMerging())
                      nodeProgenitor  => nodeProgenitor%firstChild
                      basicProgenitor => nodeProgenitor%firstChild%basic()
                   end do
                   ! Remove our mergee from its merge target.
                   call node%removeFromMergee()
                   if (associated(nodeProgenitor%firstMergee)) then
                      mergee => nodeProgenitor%firstMergee
                      do while (associated(mergee%siblingMergee))
                         mergee => mergee%siblingMergee
                      end do
                      mergee%siblingMergee => node
                   else
                      nodeProgenitor%firstMergee => node
                   end if
                   node%mergeTarget   => nodeProgenitor
                   node%siblingMergee => null()
                end if
             end if
             ! The mergee must be a satellite by the time it is due to merge. If it is not (which implies it is merging with a
             ! halo in a different branch of the tree), reset it to be merging with the appropriate halo in its own branch of
             ! the tree.
             if (associated(node%mergeTarget).and..not.node%isProgenitorOf(node%mergeTarget)) then
                nodeProgenitor  => node          %parent
                basicProgenitor => nodeProgenitor%basic ()
                if (basicProgenitor%time() > satellite%timeOfMerging()) then
                   ! Adjust the time of merging.
                   call satellite%timeOfMergingSet(basicProgenitor%time())
                   ! Remove our mergee from its merge target.
                   call node%removeFromMergee()
                   if (associated(nodeProgenitor%firstMergee)) then
                      mergee => nodeProgenitor%firstMergee
                      do while (associated(mergee%siblingMergee))
                         mergee => mergee%siblingMergee
                      end do
                      mergee%siblingMergee => node
                   else
                      nodeProgenitor%firstMergee => node
                   end if
                   node%mergeTarget   => nodeProgenitor
                   node%siblingMergee => null()
                end if
             end if
          end if
       end if
    end do
    return
  end subroutine assignOrbitsOperatePreEvolution

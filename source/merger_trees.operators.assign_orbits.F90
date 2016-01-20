!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a merger tree operator which assigns orbits to non-primary progenitor nodes.

  !# <mergerTreeOperator name="mergerTreeOperatorAssignOrbits">
  !#  <description>Provides a merger tree operator which assigns orbits to non-primary progenitor nodes.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorAssignOrbits
     !% An orbit assigning merger tree operator class.
     private
   contains
     final     ::            assignOrbitsDestructor
     procedure :: operate => assignOrbitsOperate
  end type mergerTreeOperatorAssignOrbits

  interface mergerTreeOperatorAssignOrbits
     !% Constructors for the orbit assigning merger tree operator class.
     module procedure assignOrbitsConstructorParameters
     module procedure assignOrbitsConstructorInternal
  end interface mergerTreeOperatorAssignOrbits

contains

  function assignOrbitsConstructorParameters(parameters)
    !% Constructor for the orbit assigning merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorAssignOrbits)                :: assignOrbitsConstructorParameters
    type(inputParameters               ), intent(in   ) :: parameters

    return
  end function assignOrbitsConstructorParameters

  function assignOrbitsConstructorInternal()
    !% Internal constructor for the orbit assigning merger tree operator class.
    implicit none
    type(mergerTreeOperatorAssignOrbits) :: assignOrbitsConstructorInternal
    
    return
  end function assignOrbitsConstructorInternal

  elemental subroutine assignOrbitsDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorAssignOrbits), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine assignOrbitsDestructor

  subroutine assignOrbitsOperate(self,tree)
    !% Perform a orbit assigning operation on a merger tree.
    use Virial_Orbits
    use Kepler_Orbits
    use Satellite_Merging_Timescales
    implicit none
    class(mergerTreeOperatorAssignOrbits ), intent(inout)         :: self
    type (mergerTree                     ), intent(inout), target :: tree
    type (treeNode                       ), pointer               :: node                       , nodeProgenitor, mergee
    type (mergerTree                     ), pointer               :: currentTree
    class(nodeComponentSatellite         ), pointer               :: satellite                  , satelliteProgenitor
    class(satelliteMergingTimescalesClass), pointer               :: satelliteMergingTimescales_
    type (keplerOrbit                    )                        :: virialOrbit, virialOrbitProgenitor

    logical :: satelliteProgenitorFound
    class(nodeComponentBasic), pointer :: basic, basicProgenitor

    
    ! Get required objects.
    satelliteMergingTimescales_ => satelliteMergingTimescales()
    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))   
       ! Walk the tree.
       node => currentTree%baseNode
       do while (associated(node))
          ! An orbit must be assigned to non-primary progenitors which do not have an orbit already assigned.
          if     (                                             &
               &        associated(node%parent               ) &
               &  .and.                                        &
               &   .not.           node%isPrimaryProgenitor()  &               
               & ) then
             satellite   => node     %satellite  (autoCreate=.true.)
             virialOrbit =  satellite%virialOrbit(                 )
             if (.not.virialOrbit%isDefined()) then
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
                               node%mergeTarget => nodeProgenitor
                            end if
                         end if
                      end if
                      satelliteProgenitorFound=.true.
                      exit
                   end if
                   nodeProgenitor => nodeProgenitor%firstChild
                end do
                if (.not.satelliteProgenitorFound) then
                   virialOrbit=Virial_Orbital_Parameters(node,node%parent,.false.)
                   call satellite%  mergeTimeSet(satelliteMergingTimescales_%timeUntilMerging(node,virialOrbit))
                   call satellite%virialOrbitSet(                                                  virialOrbit )
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
                      node%mergeTarget => nodeProgenitor
                   end if
                end if
             end if
          end if
          ! Walk to the next node.
          node => node%walkTree()
       end do
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine assignOrbitsOperate

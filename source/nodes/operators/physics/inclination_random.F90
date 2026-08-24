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
  Implements a node operator which assigns each galaxy a random orientation.
  !!}

  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass

  !![
  <nodeOperator name="nodeOperatorInclinationRandom" docformat="rst">
   <description>
   A node operator which assigns each galaxy an isotropically distributed orientation---uniform in :math:`\cos i`,
   so that :math:`i` is the inclination of the galaxy to the line of sight in radians, zero being face-on. The angle
   is stored in a meta-property of the ``basic`` component, from which
   :galacticus-class:`galacticInclinationRandom` reads it.

   The angle is drawn when the node is initialized, and, if ``resetOnMajorMerger`` is true, drawn again after a major
   galaxy--galaxy merger, on the grounds that such a merger destroys and rebuilds the disk and there is no reason for
   the remnant to remember the orientation of its progenitor. Whether it should is a genuine modeling choice, which
   is why it is a parameter, and why an orientation which is *stored* is worth having at all: a stateless rule
   derived from node identity could not express it.

   Storing the angle rather than drawing it where it is used is what makes it stable---see the discussion in
   :galacticus-class:`galacticInclinationRandom`. Because meta-properties move with their component, the angle also
   survives the promotion of a node up its branch without further work.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorInclinationRandom
     !!{RST
     A node operator which assigns each galaxy a random orientation.
     !!}
     private
     class  (mergerMassMovementsClass), pointer :: mergerMassMovements_ => null()
     integer                                    :: inclinationID
     logical                                    :: resetOnMajorMerger
   contains
     final     ::                    inclinationRandomDestructor
     procedure :: autoHook        => inclinationRandomAutoHook
     procedure :: nodeInitialize  => inclinationRandomNodeInitialize
  end type nodeOperatorInclinationRandom

  interface nodeOperatorInclinationRandom
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorInclinationRandom` node operator class.
     !!}
     module procedure inclinationRandomConstructorParameters
     module procedure inclinationRandomConstructorInternal
  end interface nodeOperatorInclinationRandom

contains

  function inclinationRandomConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorInclinationRandom` node operator class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodeOperatorInclinationRandom)                :: self
    type   (inputParameters              ), intent(inout) :: parameters
    class  (mergerMassMovementsClass     ), pointer       :: mergerMassMovements_
    logical                                               :: resetOnMajorMerger

    !![
    <inputParameter docformat="rst">
      <name>resetOnMajorMerger</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a new orientation is drawn after each major galaxy--galaxy merger.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerMassMovements" name="mergerMassMovements_" source="parameters"/>
    !!]
    self=nodeOperatorInclinationRandom(resetOnMajorMerger,mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerMassMovements_"/>
    !!]
    return
  end function inclinationRandomConstructorParameters

  function inclinationRandomConstructorInternal(resetOnMajorMerger,mergerMassMovements_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorInclinationRandom` node operator class.
    !!}
    implicit none
    type   (nodeOperatorInclinationRandom)                        :: self
    logical                               , intent(in   )         :: resetOnMajorMerger
    class  (mergerMassMovementsClass     ), intent(in   ), target :: mergerMassMovements_
    !![
    <constructorAssign variables="resetOnMajorMerger, *mergerMassMovements_"/>
    <addMetaProperty component="basic" name="inclination" id="self%inclinationID" isCreator="yes" isEvolvable="no"/>
    !!]

    return
  end function inclinationRandomConstructorInternal

  subroutine inclinationRandomAutoHook(self)
    !!{RST
    Attach to the satellite merger event, so that orientations can be re-drawn after a major merger.
    !!}
    use :: Events_Hooks, only : dependencyDirectionAfter, dependencyRegEx, openMPThreadBindingAtLevel, satelliteMergerEvent
    implicit none
    class(nodeOperatorInclinationRandom), intent(inout) :: self
    type (dependencyRegEx              ), dimension(1)  :: dependenciesSatelliteMerger

    if (.not.self%resetOnMajorMerger) return
    ! Attach after the remnant structure is determined, so that the merger has been classified.
    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='inclinationRandom',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine inclinationRandomAutoHook

  subroutine inclinationRandomDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorInclinationRandom` node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorInclinationRandom), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    !![
    <objectDestructor name="self%mergerMassMovements_"/>
    !!]
    return
  end subroutine inclinationRandomDestructor

  double precision function inclinationDraw(node)
    !!{RST
    Draw an inclination from an isotropic distribution of orientations, in radians.

    Orientations uniform on the sphere are uniform in :math:`\cos i`, so drawing :math:`\cos i` uniformly from
    :math:`[0,1]` and taking the arc cosine gives the required distribution over :math:`[0,\pi/2]`. Only the angle
    to the line of sight matters, not which face is presented, so the half range suffices.
    !!}
    implicit none
    type(treeNode), intent(inout), target :: node

    inclinationDraw=acos(node%hostTree%randomNumberGenerator_%uniformSample())
    return
  end function inclinationDraw

  subroutine inclinationRandomNodeInitialize(self,node)
    !!{RST
    Assign an orientation to a newly initialized node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorInclinationRandom), intent(inout), target  :: self
    type (treeNode                     ), intent(inout), target  :: node
    class(nodeComponentBasic           )               , pointer :: basic

    basic => node%basic()
    call basic%floatRank0MetaPropertySet(self%inclinationID,inclinationDraw(node))
    return
  end subroutine inclinationRandomNodeInitialize

  subroutine satelliteMerger(self,node)
    !!{RST
    Draw a new orientation for the remnant of a major galaxy--galaxy merger.
    !!}
    use :: Error                           , only : Error_Report
    use :: Function_Classes                , only : functionClass
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: ISO_Varying_String              , only : char
    use :: Satellite_Merging_Mass_Movements, only : enumerationDestinationMergerType
    implicit none
    class(*                               ), intent(inout)          :: self
    type (treeNode                        ), intent(inout), target  :: node
    type (treeNode                        )               , pointer :: nodeHost
    class(nodeComponentBasic              )               , pointer :: basicHost
    type (enumerationDestinationMergerType)                         :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                             destinationGasHost     , destinationStarsHost
    logical                                                         :: mergerIsMajor

    select type (self)
    class is (nodeOperatorInclinationRandom)
       call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       if (mergerIsMajor) then
          ! The remnant is the host: it is the galaxy which is rebuilt, so it is the one re-oriented.
          nodeHost  => node    %mergesWith()
          basicHost => nodeHost%basic     ()
          call basicHost%floatRank0MetaPropertySet(self%inclinationID,inclinationDraw(nodeHost))
       end if
    class is (functionClass)
       call Error_Report('object is not of [nodeOperatorInclinationRandom] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is of unrecognized class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger

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
  Implements a node operator which assigns each galaxy a random misalignment of the axis of a component from its own.
  !!}

  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass

  !![
  <nodeOperator name="nodeOperatorMisalignmentRandom" docformat="rst">
   <description>
   A node operator which assigns each galaxy a random misalignment of the axis of a component of it---the dusty torus
   around an AGN, for example---from the axis of the galaxy, for use by
   :galacticus-class:`galacticInclinationMisaligned`. Two numbers are stored in meta-properties of the ``basic``
   component: a deviate uniform on :math:`[0,1)`, from which that class finds the angle :math:`\psi` between the
   axes (so that the distribution of :math:`\psi` is set in one place, by that class), and the azimuth
   :math:`\phi` of the axis of the component about that of the galaxy, uniform on :math:`[0,2\pi)`.

   As for :galacticus-class:`nodeOperatorInclinationRandom`, the misalignment is drawn when the node is initialized,
   and, if ``resetOnMajorMerger`` is true, drawn again after a major galaxy--galaxy merger. Storing it, rather than
   drawing it where it is used, keeps it the same for every consumer and at every output, and it moves with the
   ``basic`` component when a node is promoted up its branch.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorMisalignmentRandom
     !!{RST
     A node operator which assigns each galaxy a random misalignment of a component from its axis.
     !!}
     private
     class  (mergerMassMovementsClass), pointer :: mergerMassMovements_  => null()
     integer                                    :: misalignmentDeviateID          , misalignmentAzimuthID
     logical                                    :: resetOnMajorMerger
   contains
     final     ::                    misalignmentRandomDestructor
     procedure :: autoHook        => misalignmentRandomAutoHook
     procedure :: nodeInitialize  => misalignmentRandomNodeInitialize
  end type nodeOperatorMisalignmentRandom

  interface nodeOperatorMisalignmentRandom
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorMisalignmentRandom` node operator class.
     !!}
     module procedure misalignmentRandomConstructorParameters
     module procedure misalignmentRandomConstructorInternal
  end interface nodeOperatorMisalignmentRandom

contains

  function misalignmentRandomConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorMisalignmentRandom` node operator class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodeOperatorMisalignmentRandom)                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    class  (mergerMassMovementsClass      ), pointer       :: mergerMassMovements_
    logical                                                :: resetOnMajorMerger

    !![
    <inputParameter docformat="rst">
      <name>resetOnMajorMerger</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a new misalignment is drawn after each major galaxy--galaxy merger.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerMassMovements" name="mergerMassMovements_" source="parameters"/>
    !!]
    self=nodeOperatorMisalignmentRandom(resetOnMajorMerger,mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerMassMovements_"/>
    !!]
    return
  end function misalignmentRandomConstructorParameters

  function misalignmentRandomConstructorInternal(resetOnMajorMerger,mergerMassMovements_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorMisalignmentRandom` node operator class.
    !!}
    implicit none
    type   (nodeOperatorMisalignmentRandom)                        :: self
    logical                                , intent(in   )         :: resetOnMajorMerger
    class  (mergerMassMovementsClass      ), intent(in   ), target :: mergerMassMovements_
    !![
    <constructorAssign variables="resetOnMajorMerger, *mergerMassMovements_"/>
    <addMetaProperty component="basic" name="misalignmentDeviate" id="self%misalignmentDeviateID" isCreator="yes" isEvolvable="no"/>
    <addMetaProperty component="basic" name="misalignmentAzimuth" id="self%misalignmentAzimuthID" isCreator="yes" isEvolvable="no"/>
    !!]

    return
  end function misalignmentRandomConstructorInternal

  subroutine misalignmentRandomAutoHook(self)
    !!{RST
    Attach to the satellite merger event, so that misalignments can be re-drawn after a major merger.
    !!}
    use :: Events_Hooks, only : dependencyDirectionAfter, dependencyRegEx, openMPThreadBindingAtLevel, satelliteMergerEvent
    implicit none
    class(nodeOperatorMisalignmentRandom), intent(inout) :: self
    type (dependencyRegEx               ), dimension(1)  :: dependenciesSatelliteMerger

    if (.not.self%resetOnMajorMerger) return
    ! Attach after the remnant structure is determined, so that the merger has been classified.
    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='misalignmentRandom',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine misalignmentRandomAutoHook

  subroutine misalignmentRandomDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorMisalignmentRandom` node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorMisalignmentRandom), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    !![
    <objectDestructor name="self%mergerMassMovements_"/>
    !!]
    return
  end subroutine misalignmentRandomDestructor

  subroutine misalignmentDraw(self,node)
    !!{RST
    Draw a misalignment for the galaxy in ``node``, and store it.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(nodeOperatorMisalignmentRandom), intent(inout)          :: self
    type (treeNode                      ), intent(inout), target  :: node
    class(nodeComponentBasic            )               , pointer :: basic

    basic => node%basic()
    call basic%floatRank0MetaPropertySet(self%misalignmentDeviateID,         node%hostTree%randomNumberGenerator_%uniformSample())
    call basic%floatRank0MetaPropertySet(self%misalignmentAzimuthID,2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample())
    return
  end subroutine misalignmentDraw

  subroutine misalignmentRandomNodeInitialize(self,node)
    !!{RST
    Assign a misalignment to a newly initialized node.
    !!}
    implicit none
    class(nodeOperatorMisalignmentRandom), intent(inout), target :: self
    type (treeNode                      ), intent(inout), target :: node

    call misalignmentDraw(self,node)
    return
  end subroutine misalignmentRandomNodeInitialize

  subroutine satelliteMerger(self,node)
    !!{RST
    Draw a new misalignment for the remnant of a major galaxy--galaxy merger.
    !!}
    use :: Error                           , only : Error_Report
    use :: Function_Classes                , only : functionClass
    use :: ISO_Varying_String              , only : char
    use :: Satellite_Merging_Mass_Movements, only : enumerationDestinationMergerType
    implicit none
    class(*                               ), intent(inout)          :: self
    type (treeNode                        ), intent(inout), target  :: node
    type (treeNode                        )               , pointer :: nodeHost
    type (enumerationDestinationMergerType)                         :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                             destinationGasHost     , destinationStarsHost
    logical                                                         :: mergerIsMajor

    select type (self)
    class is (nodeOperatorMisalignmentRandom)
       call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       if (mergerIsMajor) then
          ! The remnant is the host: it is the galaxy which is rebuilt, so it is the one re-oriented.
          nodeHost => node%mergesWith()
          call misalignmentDraw(self,nodeHost)
       end if
    class is (functionClass)
       call Error_Report('object is not of [nodeOperatorMisalignmentRandom] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is of unrecognized class'                                                                      //{introspection:location})
    end select
    return
  end subroutine satelliteMerger

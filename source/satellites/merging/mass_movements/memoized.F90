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

  !+    Contributions to this file made by: Claude.

  !!{RST
  Implements an abstract merger mass movements class which computes the movements once per merger and stores them.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <mergerMassMovements name="mergerMassMovementsMemoized" abstract="yes" docformat="rst">
   <description>
   An abstract merger mass movements class which computes the mass movements for a merger when they are first needed---at the satellite merger event, before any node component responds to it---and stores them, so that the movements are determined by the properties of the merging galaxies prior to any modification by the merger. The stored movements are discarded by the calculation reset event. This is an abstract class---the movements themselves must be computed by a concrete class.
   </description>
  </mergerMassMovements>
  !!]
  type, abstract, extends(mergerMassMovementsClass) :: mergerMassMovementsMemoized
     !!{RST
     An abstract merger mass movements class which computes the movements once per merger and stores them.
     !!}
     private
     integer(kind=kind_int8                  ) :: lastUniqueID           =-huge(0_kind_int8)
     logical                                   :: movementsCalculated    =.false.
     type   (enumerationDestinationMergerType) :: destinationGasSatellite                   , destinationStarsSatellite, &
          &                                       destinationGasHost                        , destinationStarsHost
     logical                                   :: mergerIsMajor
   contains
     !![
     <methods docformat="rst">
       <method description="Compute the movements of stellar and gaseous mass components during a galaxy merger event." method="calculate"  />
       <method description="Detach from the node events. Must be called by the destructor of each concrete class."      method="detachHooks"/>
     </methods>
     !!]
     procedure                                  :: autoHook    => memoizedAutoHook
     procedure                                  :: detachHooks => memoizedDetachHooks
     procedure                                  :: get         => memoizedGet
     procedure(memoizedCalculate    ), deferred :: calculate
  end type mergerMassMovementsMemoized

  abstract interface
     subroutine memoizedCalculate(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       !!{RST
       Interface for the calculation of mass movements by memoized merger mass movements classes.
       !!}
       import mergerMassMovementsMemoized, treeNode, enumerationDestinationMergerType
       class  (mergerMassMovementsMemoized     ), intent(inout)         :: self
       type   (treeNode                        ), intent(inout), target :: node
       type   (enumerationDestinationMergerType), intent(  out)         :: destinationGasSatellite, destinationStarsSatellite, &
            &                                                              destinationGasHost     , destinationStarsHost
       logical                                  , intent(  out)         :: mergerIsMajor
     end subroutine memoizedCalculate
  end interface

contains

  subroutine memoizedAutoHook(self)
    !!{RST
    Attach to the calculation reset and satellite merger events.
    !!}
    use :: Events_Hooks      , only : calculationResetEvent, openMPThreadBindingAllLevels, satelliteMergerEvent
    use :: ISO_Varying_String, only : char
    implicit none
    class(mergerMassMovementsMemoized), intent(inout) :: self

    call calculationResetEvent%attach(self,memoizedCalculationReset,openMPThreadBindingAllLevels,label='remnantStructure:'//char(self%objectType()))
    call satelliteMergerEvent %attach(self,memoizedGetHook         ,openMPThreadBindingAllLevels,label='remnantStructure:'//char(self%objectType()))
    return
  end subroutine memoizedAutoHook

  subroutine memoizedDetachHooks(self)
    !!{RST
    Detach from the calculation reset and satellite merger events.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent
    implicit none
    class(mergerMassMovementsMemoized), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,memoizedCalculationReset)) call calculationResetEvent%detach(self,memoizedCalculationReset)
    if (satelliteMergerEvent %isAttached(self,memoizedGetHook         )) call satelliteMergerEvent %detach(self,memoizedGetHook         )
    return
  end subroutine memoizedDetachHooks

  subroutine memoizedCalculationReset(self,node,uniqueID)
    !!{RST
    Reset the stored mass movements.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    use :: Function_Classes  , only : functionClass
    implicit none
    class  (*        ), intent(inout) :: self
    type   (treeNode ), intent(inout) :: node
    integer(kind_int8), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    select type (self)
    class is (mergerMassMovementsMemoized)
       self%movementsCalculated=.false.
       self%lastUniqueID       =uniqueID
    class is (functionClass)
       call Error_Report('object is not of [mergerMassMovementsMemoized] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is not of [mergerMassMovementsMemoized] class'//{introspection:location})
    end select
    return
  end subroutine memoizedCalculationReset

  subroutine memoizedGetHook(self,node)
    !!{RST
    Hookable wrapper around the get function.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    use :: Function_Classes  , only : functionClass
    implicit none
    class  (*                               ), intent(inout)         :: self
    type   (treeNode                        ), intent(inout), target :: node
    type   (enumerationDestinationMergerType)                        :: destinationGasSatellite, destinationGasHost       , &
         &                                                              destinationStarsHost   , destinationStarsSatellite
    logical                                                          :: mergerIsMajor

    select type (self)
    class is (mergerMassMovementsMemoized)
       call self%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    class is (functionClass)
       call Error_Report('object is not of [mergerMassMovementsMemoized] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is not of [mergerMassMovementsMemoized] class'//{introspection:location})
    end select
    return
  end subroutine memoizedGetHook

  subroutine memoizedGet(self,node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    !!{RST
    Determine where stars and gas move as the result of a merger event, computing the movements when first needed and storing them.
    !!}
    implicit none
    class  (mergerMassMovementsMemoized     ), intent(inout)         :: self
    type   (treeNode                        ), intent(inout), target :: node
    type   (enumerationDestinationMergerType), intent(  out)         :: destinationGasSatellite, destinationGasHost       , &
         &                                                              destinationStarsHost   , destinationStarsSatellite
    logical                                  , intent(  out)         :: mergerIsMajor

    ! The calculation of how mass moves as a result of the merger is computed when first needed and then stored. This ensures that
    ! the results are determined by the properties of the merge target prior to any modification that will occur as node
    ! components are modified in response to the merger.
    if (node%uniqueID() /= self%lastUniqueID) call memoizedCalculationReset(self,node,node%uniqueID())
    if (.not.self%movementsCalculated) then
       self%movementsCalculated=.true.
       call self%calculate(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       self%destinationGasSatellite  =destinationGasSatellite
       self%destinationStarsSatellite=destinationStarsSatellite
       self%destinationGasHost       =destinationGasHost
       self%destinationStarsHost     =destinationStarsHost
       self%mergerIsMajor            =mergerIsMajor
    else
       destinationGasSatellite       =self%destinationGasSatellite
       destinationStarsSatellite     =self%destinationStarsSatellite
       destinationGasHost            =self%destinationGasHost
       destinationStarsHost          =self%destinationStarsHost
       mergerIsMajor                 =self%mergerIsMajor
    end if
    return
  end subroutine memoizedGet

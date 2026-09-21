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
  Implements an abstract merger remnant size class which computes the remnant properties once per merger and stores them.
  !!}

  use :: Kind_Numbers, only : kind_int8

  !![
  <mergerRemnantSize name="mergerRemnantSizeMemoized" abstract="yes" docformat="rst">
   <description>
   An abstract merger remnant size class which computes the size, circular velocity, and specific angular momentum of a merger remnant when they are first needed---at the satellite merger event, before any node component responds to it---and stores them, so that they are determined by the properties of the merging galaxies prior to any modification by the merger. The stored properties are discarded by the calculation reset event. This is an abstract class---the properties themselves must be computed by a concrete class.
   </description>
  </mergerRemnantSize>
  !!]
  type, abstract, extends(mergerRemnantSizeClass) :: mergerRemnantSizeMemoized
     !!{RST
     An abstract merger remnant size class which computes the remnant properties once per merger and stores them.
     !!}
     private
     integer         (kind=kind_int8) :: lastUniqueID        =-huge(0_kind_int8)
     logical                          :: propertiesCalculated=.false.
     double precision                 :: radius                                 , velocityCircular, &
          &                              angularMomentumSpecific
   contains
     !![
     <methods docformat="rst">
       <method description="Compute the size, circular velocity, and specific angular momentum of the merger remnant." method="calculate"   />
       <method description="Detach from the node events. Must be called by the destructor of each concrete class."    method="detachHooks" />
     </methods>
     !!]
     procedure                              :: autoHook    => memoizedAutoHook
     procedure                              :: detachHooks => memoizedDetachHooks
     procedure                              :: get         => memoizedGet
     procedure(memoizedCalculate), deferred :: calculate
  end type mergerRemnantSizeMemoized

  abstract interface
     subroutine memoizedCalculate(self,node,radius,velocityCircular,angularMomentumSpecific)
       !!{RST
       Interface for the calculation of remnant properties by memoized merger remnant size classes.
       !!}
       import mergerRemnantSizeMemoized, treeNode
       class           (mergerRemnantSizeMemoized), intent(inout) :: self
       type            (treeNode                 ), intent(inout) :: node
       double precision                           , intent(  out) :: radius, velocityCircular, angularMomentumSpecific
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
    class(mergerRemnantSizeMemoized), intent(inout) :: self

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
    class(mergerRemnantSizeMemoized), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,memoizedCalculationReset)) call calculationResetEvent%detach(self,memoizedCalculationReset)
    if (satelliteMergerEvent %isAttached(self,memoizedGetHook         )) call satelliteMergerEvent %detach(self,memoizedGetHook         )
    return
  end subroutine memoizedDetachHooks

  subroutine memoizedCalculationReset(self,node,uniqueID)
    !!{RST
    Reset the stored remnant properties.
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
    class is (mergerRemnantSizeMemoized)
       self%propertiesCalculated=.false.
       self%lastUniqueID        =uniqueID
    class is (functionClass)
       call Error_Report('object is not of [mergerRemnantSizeMemoized] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is not of [mergerRemnantSizeMemoized] class'//{introspection:location})
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
    class           (*       ), intent(inout)         :: self
    type            (treeNode), intent(inout), target :: node
    double precision                                  :: radius                 , velocityCircular, &
         &                                               angularMomentumSpecific

    select type (self)
    class is (mergerRemnantSizeMemoized)
       call self%get(node,radius,velocityCircular,angularMomentumSpecific)
    class is (functionClass)
       call Error_Report('object is not of [mergerRemnantSizeMemoized] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is not of [mergerRemnantSizeMemoized] class'//{introspection:location})
    end select
    return
  end subroutine memoizedGetHook

  subroutine memoizedGet(self,node,radius,velocityCircular,angularMomentumSpecific)
    !!{RST
    Compute the size of the merger remnant for ``node``, computing the remnant properties when first needed and storing them.
    !!}
    implicit none
    class           (mergerRemnantSizeMemoized), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: node
    double precision                           , intent(  out) :: radius, velocityCircular, angularMomentumSpecific

    ! The remnant properties are computed when first needed and then stored. This ensures that they are determined by the
    ! properties of the merging galaxies prior to any modification that will occur as node components are modified in response
    ! to the merger.
    if (node%uniqueID() /= self%lastUniqueID) call memoizedCalculationReset(self,node,node%uniqueID())
    if (.not.self%propertiesCalculated) then
       self%propertiesCalculated=.true.
       call self%calculate(node,radius,velocityCircular,angularMomentumSpecific)
       self%radius                 =radius
       self%velocityCircular       =velocityCircular
       self%angularMomentumSpecific=angularMomentumSpecific
    else
       radius                 =self%radius
       velocityCircular       =self%velocityCircular
       angularMomentumSpecific=self%angularMomentumSpecific
    end if
    return
  end subroutine memoizedGet

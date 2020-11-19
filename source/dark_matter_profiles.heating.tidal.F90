!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% A dark matter halo profile heating class which accounts for heating from tidal shocking.

  use :: Kind_Numbers, only : kind_int8

  !# <darkMatterProfileHeating name="darkMatterProfileHeatingTidal">
  !#  <description>A dark matter profile heating model which accounts for heating due to tidal shocking.</description>
  !# </darkMatterProfileHeating>
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingTidal
     !% A dark matter profile heating class which accounts for heating due to tidal shocking.
     private
    double precision            :: specificEnergyOverRadiusSquared_, specificEnergyOverRadiusSquaredParent_
    integer         (kind_int8) :: lastUniqueID                    , parentUniqueID
  contains
     !# <methods>
     !#   <method description="Reset memoized calculations." method="calculationReset" />
     !#   <method description="Compute $Q = E / r^2$." method="specificEnergyOverRadiusSquared" />
     !# </methods>
     final     ::                                    tidalDestructor
     procedure :: autoHook                        => tidalAutoHook
     procedure :: calculationReset                => tidalCalculationReset
     procedure :: specificEnergy                  => tidalSpecificEnergy
     procedure :: specificEnergyGradient          => tidalSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero  => tidalSpecificEnergyIsEverywhereZero
     procedure :: specificEnergyOverRadiusSquared => tidalSpecificEnergyOverRadiusSquared
  end type darkMatterProfileHeatingTidal

  interface darkMatterProfileHeatingTidal
     !% Constructors for the {\normalfont \ttfamily tidal} dark matter profile heating class.
     module procedure tidalConstructorParameters
     module procedure tidalConstructorInternal
  end interface darkMatterProfileHeatingTidal

contains

  function tidalConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily tidal} dark matter profile heating scales class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterProfileHeatingTidal), target        :: self
    type(inputParameters              ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=darkMatterProfileHeatingTidal()
    return
  end function tidalConstructorParameters

  function tidalConstructorInternal() result(self)
    !% Internal constructor for the {\normalfont \ttfamily tidal} dark matter profile heating scales class.
    implicit none
    type(darkMatterProfileHeatingTidal) :: self

    self%specificEnergyOverRadiusSquared_      =-1.0d0
    self%specificEnergyOverRadiusSquaredParent_=-1.0d0
    self%lastUniqueID                          =-1_kind_int8
    self%parentUniqueID                        =-1_kind_int8
    return
  end function tidalConstructorInternal

  subroutine tidalAutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self

    call calculationResetEvent%attach(self,tidalCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine tidalAutoHook

  subroutine tidalCalculationReset(self,node)
    !% Reset the stored tidal radii.
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node

    self   %specificEnergyOverRadiusSquared_      =-1.0d0
    self   %specificEnergyOverRadiusSquaredParent_=-1.0d0
    self   %lastUniqueID                          =node       %uniqueID()
    if (associated(node%parent)) then
       self%parentUniqueID                        =node%parent%uniqueID()
    else
       self%parentUniqueID                        =-1_kind_int8
    end if
    return
  end subroutine tidalCalculationReset

  subroutine tidalDestructor(self)
    !% Destructor for the {\normalfont \ttfamily tidal} dark matter profile heating class.
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileHeatingTidal), intent(inout) :: self

    call calculationResetEvent%detach(self,tidalCalculationReset)
    return
  end subroutine tidalDestructor

  double precision function tidalSpecificEnergy(self,node,darkMatterProfileDMO_,radius)
    !% Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    implicit none
    class           (darkMatterProfileHeatingTidal), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    double precision                               , intent(in   ) :: radius
    !$GLC attributes unused :: darkMatterProfileDMO_

    tidalSpecificEnergy=+self%specificEnergyOverRadiusSquared(node)    &
         &              *radius                                    **2
    return
  end function tidalSpecificEnergy

  double precision function tidalSpecificEnergyGradient(self,node,darkMatterProfileDMO_,radius)
    !% Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    implicit none
    class           (darkMatterProfileHeatingTidal), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    double precision                               , intent(in   ) :: radius
    !$GLC attributes unused :: darkMatterProfileDMO_

    tidalSpecificEnergyGradient=+2.0d0                                      &
         &                      *self%specificEnergyOverRadiusSquared(node) &
         &                      *radius
    return
  end function tidalSpecificEnergyGradient

  logical function tidalSpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !% Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    !$GLC attributes unused :: darkMatterProfileDMO_

    tidalSpecificEnergyIsEverywhereZero=self%specificEnergyOverRadiusSquared(node) <= 0.0d0
    return
  end function tidalSpecificEnergyIsEverywhereZero

  double precision function tidalSpecificEnergyOverRadiusSquared(self,node)
    !% Compute $Q = E / r^2$.
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    class(darkMatterProfileHeatingTidal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(nodeComponentSatellite       ), pointer       :: satellite

    if     (                                        &
         &   node%uniqueID() /= self%parentUniqueID &
         &  .and.                                   &
         &   node%uniqueID() /= self%lastUniqueID   &
         & ) call self%calculationReset(node)
    if (node%uniqueID() == self%parentUniqueID) then
       if (self%specificEnergyOverRadiusSquaredParent_ < 0.0d0) then
          satellite                                   =>      node     %satellite             ()
          self%specificEnergyOverRadiusSquaredParent_ =  max(                                     &
               &                                             +0.0d0                             , &
               &                                             +satellite%tidalHeatingNormalized()  &
               &                                            )
       end if
       tidalSpecificEnergyOverRadiusSquared=self%specificEnergyOverRadiusSquaredParent_
    else
       if (self%specificEnergyOverRadiusSquared_       < 0.0d0) then
          satellite                                   =>      node     %satellite             ()
          self%specificEnergyOverRadiusSquared_       =  max(                                     &
               &                                             +0.0d0                             , &
               &                                             +satellite%tidalHeatingNormalized()  &
               &                                            )
       end if
       tidalSpecificEnergyOverRadiusSquared=self%specificEnergyOverRadiusSquared_
    end if
    return
  end function tidalSpecificEnergyOverRadiusSquared

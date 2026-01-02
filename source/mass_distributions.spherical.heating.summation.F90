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
  Implements a mass distribution heating class that sums heating over other classes.
  !!}

  !![
  <massDistributionHeating name="massDistributionHeatingSummation">
     <description>A mass distribution heating class that sums heating over other classes.</description>
     <linkedList type="massDistributionHeatingList" variable="massDistributionHeatings" next="next" object="massDistributionHeating_" objectType="massDistributionHeatingClass"/>
  </massDistributionHeating>
  !!]
 
  type, public :: massDistributionHeatingList
     class(massDistributionHeatingClass), pointer :: massDistributionHeating_ => null()
     type (massDistributionHeatingList ), pointer :: next                     => null()
  end type massDistributionHeatingList
  
  type, extends(massDistributionHeatingClass) :: massDistributionHeatingSummation
     !!{
     Implementation of a mass distribution heating class that sums heating over other classes.
     !!}
     private
      type(massDistributionHeatingList), pointer :: massDistributionHeatings => null()
   contains
     procedure :: specificEnergy                 => summationSpecificEnergy
     procedure :: specificEnergyGradient         => summationSpecificEnergyGradient
     procedure :: specificEnergyIsEveryWhereZero => summationSpecificEnergyIsEverywhereZero
  end type massDistributionHeatingSummation

  interface massDistributionHeatingSummation
     !!{
     Constructors for the \refClass{massDistributionHeatingSummation} mass distribution class.
     !!}
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface massDistributionHeatingSummation

contains

  function summationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingSummation} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (massDistributionHeatingSummation)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    type   (massDistributionHeatingList     ), pointer       :: massDistributionHeating_
    integer                                                  :: i

    massDistributionHeating_ => null()
    do i=1,parameters%copiesCount('massDistributionHeating',zeroIfNotPresent=.true.)
       if (associated(massDistributionHeating_)) then
          allocate(massDistributionHeating_%next)
          massDistributionHeating_ => massDistributionHeating_%next
       else
          allocate(self%massDistributionHeatings)
          massDistributionHeating_ => self                    %massDistributionHeatings
       end if
       !![
       <objectBuilder class="massDistributionHeating" name="massDistributionHeating_%massDistributionHeating_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="mmassDistributionHeating"/>
    !!]
    return
  end function summationConstructorParameters
  
  function summationConstructorInternal(massDistributionHeatings) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingSummation} dark matter profile heating class.
    !!}
    implicit none
    type(massDistributionHeatingSummation)                         :: self
    type(massDistributionHeatingList     ), pointer, intent(in   ) :: massDistributionHeatings
    type(massDistributionHeatingList     ), pointer                :: massDistributionHeating_

    self             %massDistributionHeatings => massDistributionHeatings
    massDistributionHeating_                   => massDistributionHeatings
    do while (associated(massDistributionHeating_))
       !![
       <referenceCountIncrement owner="massDistributionHeating_" object="massDistributionHeating_"/>
       !!]
       massDistributionHeating_ => massDistributionHeating_%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !!{
    Destructor for composite mass distributions.
    !!}
    implicit none
    type(massDistributionHeatingSummation), intent(inout) :: self
    type(massDistributionHeatingList     ), pointer       :: massDistributionHeating_, massDistributionHeatingNext

    if (associated(self%massDistributionHeatings)) then
       massDistributionHeating_ => self%massDistributionHeatings
       do while (associated(massDistributionHeating_))
          massDistributionHeatingNext => massDistributionHeating_%next
          !![
          <objectDestructor name="massDistributionHeating_%massDistributionHeating_"/>
          !!]
          deallocate(massDistributionHeating_)
          massDistributionHeating_ => massDistributionHeatingNext
       end do
    end if
    return
  end subroutine summationDestructor

  double precision function summationSpecificEnergy(self,radius,massDistribution_) result(energySpecific)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class           (massDistributionHeatingSummation), intent(inout) :: self
    double precision                                  , intent(in   ) :: radius
    class           (massDistributionClass           ), intent(inout) :: massDistribution_
    type            (massDistributionHeatingList     ), pointer       :: massDistributionHeating_

    energySpecific           =  0.0d0
    massDistributionHeating_ => self%massDistributionHeatings
    do while (associated(massDistributionHeating_))
       energySpecific           =  +                                                  energySpecific                          &
            &                      +massDistributionHeating_%massDistributionHeating_%specificEnergy(radius,massDistribution_)
       massDistributionHeating_ =>  massDistributionHeating_%next
    end do
    return
  end function summationSpecificEnergy

  double precision function summationSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating.
    !!}
    implicit none
    class           (massDistributionHeatingSummation), intent(inout) :: self
    double precision                                  , intent(in   ) :: radius
    class           (massDistributionClass           ), intent(inout) :: massDistribution_
    type            (massDistributionHeatingList     ), pointer       :: massDistributionHeating_

    energySpecificGradient   =  0.0d0
    massDistributionHeating_ => self%massDistributionHeatings
    do while (associated(massDistributionHeating_))
       energySpecificGradient   =  +                                                  energySpecificGradient                          &
            &                      +massDistributionHeating_%massDistributionHeating_%specificEnergyGradient(radius,massDistribution_)
       massDistributionHeating_ =>  massDistributionHeating_%next
    end do
    return
  end function summationSpecificEnergyGradient

  logical function summationSpecificEnergyIsEverywhereZero(self) result(energySpecificIsEverywhereZero)
    !!{
    Returns true if the specific energy is everywhere zero.
    !!}
    implicit none
    class(massDistributionHeatingSummation), intent(inout) :: self
    type (massDistributionHeatingList     ), pointer       :: massDistributionHeating_

    energySpecificIsEverywhereZero =  .true.
    massDistributionHeating_       => self%massDistributionHeatings
    do while (associated(massDistributionHeating_))
       energySpecificIsEverywhereZero=massDistributionHeating_%massDistributionHeating_%specificEnergyIsEverywhereZero()
       if (.not.energySpecificIsEverywhereZero) return
       massDistributionHeating_ =>  massDistributionHeating_%next
    end do
    return
  end function summationSpecificEnergyIsEverywhereZero

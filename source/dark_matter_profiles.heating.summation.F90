!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% A dark matter halo profile heating class which sums over other heat sources.

  !# <darkMatterProfileHeating name="darkMatterProfileHeatingSummation">
  !#  <description>A dark matter profile heating model which sums over other heat sources.</description>
  !# </darkMatterProfileHeating>

  type, public :: heatSourceList
     class(darkMatterProfileHeatingClass), pointer :: heatSource
     type (heatSourceList               ), pointer :: next       => null()
  end type heatSourceList

  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingSummation
     !% A dark matter profile heating class which sums over other heat sources.
     private
     type(heatSourceList), pointer :: heatSources
   contains
     final     ::                                   summationDestructor
     procedure :: specificEnergy                 => summationSpecificEnergy
     procedure :: specificEnergyGradient         => summationSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero => summationSpecificEnergyIsEverywhereZero
     procedure :: deepCopy                       => summationDeepCopy
  end type darkMatterProfileHeatingSummation

  interface darkMatterProfileHeatingSummation
     !% Constructors for the {\normalfont \ttfamily summation} dark matter profile heating class.
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface darkMatterProfileHeatingSummation

contains

  function summationConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily summation} dark matter profile heating class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type   (darkMatterProfileHeatingSummation), target        :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    type   (heatSourceList                   ), pointer       :: heatSource
    integer                                                   :: i

    heatSource => null()
    do i=1,parameters%copiesCount('darkMatterProfileHeatingMethod',zeroIfNotPresent=.true.)
       if (associated(heatSource)) then
          allocate(heatSource%next)
          heatSource => heatSource%next
       else
          allocate(self%heatSources)
          heatSource => self%heatSources
       end if
       !# <objectBuilder class="darkMatterProfileHeating" name="heatSource%heatSource" source="parameters" copy="i" />
    end do
    return
  end function summationConstructorParameters

  function summationConstructorInternal(heatSources) result(self)
    !% Internal constructor for the ``summation'' dark matter profile heating class.
    implicit none
    type(darkMatterProfileHeatingSummation)                        :: self
    type(heatSourceList                   ), target, intent(in   ) :: heatSources
    type(heatSourceList                   ), pointer               :: heatSource_

    self       %heatSources => heatSources
    heatSource_             => heatSources
    do while (associated(heatSource_))
       !# <referenceCountIncrement owner="heatSource_" object="heatSource"/>
       heatSource_ => heatSource_%next
    end do
    return
  end function summationConstructorInternal
  
  subroutine summationDestructor(self)
    !% Destructor for the ``summation'' dark matter profile heating class.
    implicit none
    type(darkMatterProfileHeatingSummation), intent(inout) :: self
    type(heatSourceList                   ), pointer       :: heatSource, heatSourceNext

    heatSource => self%heatSources
    do while (associated(heatSource))
       heatSourceNext => heatSource    %next
       !# <objectDestructor name="heatSource%heatSource"/>
       deallocate(heatSource)
       heatSource => heatSourceNext
    end do
    return
  end subroutine summationDestructor

  double precision function summationSpecificEnergy(self,node,darkMatterProfile_,radius)
    !% Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    implicit none
    class           (darkMatterProfileHeatingSummation), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    class           (darkMatterProfileClass           ), intent(inout) :: darkMatterProfile_
    double precision                                   , intent(in   ) :: radius
    type            (heatSourceList                   ), pointer       :: heatSource

    summationSpecificEnergy =  0.0d0
    heatSource              => self%heatSources
    do while (associated(heatSource))
       summationSpecificEnergy=+summationSpecificEnergy                                  &
            &                  +heatSource%heatSource%specificEnergy(                    &
            &                                                        node              , &
            &                                                        darkMatterProfile_, &
            &                                                        radius              &
            &                                                       )
       heatSource => heatSource%next
    end do
    return
  end function summationSpecificEnergy

  double precision function summationSpecificEnergyGradient(self,node,darkMatterProfile_,radius)
    !% Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    implicit none
    class           (darkMatterProfileHeatingSummation), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    class           (darkMatterProfileClass           ), intent(inout) :: darkMatterProfile_
    double precision                                   , intent(in   ) :: radius
    type            (heatSourceList                   ), pointer       :: heatSource

    summationSpecificEnergyGradient =  0.0d0
    heatSource                      => self%heatSources
    do while (associated(heatSource))
       summationSpecificEnergyGradient=+summationSpecificEnergyGradient                                  &
            &                          +heatSource%heatSource%specificEnergyGradient(                    &
            &                                                                        node              , &
            &                                                                        darkMatterProfile_, &
            &                                                                        radius              &
            &                                                                       )
       heatSource => heatSource%next
    end do
    return
  end function summationSpecificEnergyGradient

  logical function summationSpecificEnergyIsEverywhereZero(self,node,darkMatterProfile_)
    !% Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileHeatingSummation), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    class(darkMatterProfileClass           ), intent(inout) :: darkMatterProfile_
    type (heatSourceList                   ), pointer       :: heatSource

    summationSpecificEnergyIsEverywhereZero =  .true.
    heatSource                              => self%heatSources
    do while (associated(heatSource))
       if (heatSource%heatSource%specificEnergyIsEverywhereZero(node,darkMatterProfile_)) then
          summationSpecificEnergyIsEverywhereZero=.false.
          exit
       end if
       heatSource => heatSource%next
    end do
    return
  end function summationSpecificEnergyIsEverywhereZero

  subroutine summationDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily summation} dark matter halo profile heating class.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileHeatingSummation), intent(inout) :: self
    class(darkMatterProfileHeatingClass    ), intent(  out) :: destination
    type (heatSourceList                   ), pointer       :: heatSource_   , heatSourceDestination_, &
         &                                                     heatSourceNew_

    call self%darkMatterProfileHeatingClass%deepCopy(destination)
    select type (destination)
    type is (darkMatterProfileHeatingSummation)
       destination%heatSources => null          ()
       heatSourceDestination_  => null          ()
       heatSource_             => self%heatSources
       do while (associated(heatSource_))
          allocate(heatSourceNew_)
          if (associated(heatSourceDestination_)) then
             heatSourceDestination_%next       => heatSourceNew_
             heatSourceDestination_            => heatSourceNew_             
          else
             destination          %heatSources => heatSourceNew_
             heatSourceDestination_            => heatSourceNew_
          end if
          allocate(heatSourceNew_%heatSource,mold=heatSource_%heatSource)
          !# <deepCopy source="heatSource_%heatSource" destination="heatSourceNew_%heatSource"/>
          heatSource_ => heatSource_%next
       end do       
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine summationDeepCopy

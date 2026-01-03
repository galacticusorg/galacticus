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
  A dark matter halo profile heating class which sums over other heat sources.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingSummation">
   <description>A dark matter profile heating model which sums over other heat sources.</description>
   <linkedList type="heatSourceList" variable="heatSources" next="next" object="heatSource" objectType="darkMatterProfileHeatingClass"/>
  </darkMatterProfileHeating>
  !!]

  type, public :: heatSourceList
     class(darkMatterProfileHeatingClass), pointer :: heatSource => null()
     type (heatSourceList               ), pointer :: next       => null()
  end type heatSourceList

  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingSummation
     !!{
     A dark matter profile heating class which sums over other heat sources.
     !!}
     private
     type(heatSourceList), pointer :: heatSources => null()
   contains
     final     ::        summationDestructor
     procedure :: get => summationGet
  end type darkMatterProfileHeatingSummation

  interface darkMatterProfileHeatingSummation
     !!{
     Constructors for the \refClass{darkMatterProfileHeatingSummation} dark matter profile heating class.
     !!}
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface darkMatterProfileHeatingSummation

contains

  function summationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileHeatingSummation} dark matter profile heating class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (darkMatterProfileHeatingSummation), target        :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    type   (heatSourceList                   ), pointer       :: heatSource
    integer                                                   :: i

    heatSource => null()
    do i=1,parameters%copiesCount('darkMatterProfileHeating',zeroIfNotPresent=.true.)
       if (associated(heatSource)) then
          allocate(heatSource%next)
          heatSource => heatSource%next
       else
          allocate(self%heatSources)
          heatSource => self%heatSources
       end if
       !![
       <objectBuilder class="darkMatterProfileHeating" name="heatSource%heatSource" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="darkMatterProfileHeating"/>
    !!]
    return
  end function summationConstructorParameters

  function summationConstructorInternal(heatSources) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileHeatingSummation} dark matter profile heating class.
    !!}
    implicit none
    type(darkMatterProfileHeatingSummation)                        :: self
    type(heatSourceList                   ), target, intent(in   ) :: heatSources
    type(heatSourceList                   ), pointer               :: heatSource_

    self       %heatSources => heatSources
    heatSource_             => heatSources
    do while (associated(heatSource_))
       !![
       <referenceCountIncrement owner="heatSource_" object="heatSource"/>
       !!]
       heatSource_ => heatSource_%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileHeatingSummation} dark matter profile heating class.
    !!}
    implicit none
    type(darkMatterProfileHeatingSummation), intent(inout) :: self
    type(heatSourceList                   ), pointer       :: heatSource, heatSourceNext

    heatSource => self%heatSources
    do while (associated(heatSource))
       heatSourceNext => heatSource    %next
       !![
       <objectDestructor name="heatSource%heatSource"/>
       !!]
       deallocate(heatSource)
       heatSource => heatSourceNext
    end do
    return
  end subroutine summationDestructor

  function summationGet(self,node) result(massDistributionHeating_)
    !!{
    Return the dark matter mass distribution heating for the given {\normalfont \ttfamily node}.
    !!}
    use :: Mass_Distributions, only : massDistributionHeatingSummation, massDistributionHeatingList
    implicit none
    class(massDistributionHeatingClass     ), pointer       :: massDistributionHeating_
    class(darkMatterProfileHeatingSummation), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    type (massDistributionHeatingList      ), pointer       :: massDistributionHeatings, massDistributionHeating__
    type (heatSourceList                   ), pointer       :: heatSource

    ! Create the mass distribution.
    allocate(massDistributionHeatingSummation :: massDistributionHeating_)
    select type(massDistributionHeating_)
    type is (massDistributionHeatingSummation)
       heatSource                => self%heatSources
       massDistributionHeating__ => null()
       massDistributionHeatings  => null()
       do while (associated(heatSource))
          if (associated(massDistributionHeatings)) then
             allocate(massDistributionHeating__%next)
             massDistributionHeating__ => massDistributionHeating__%next
          else
             allocate(massDistributionHeatings      )
             massDistributionHeating__ => massDistributionHeatings
          end if
          massDistributionHeating__%massDistributionHeating_ => heatSource%heatSource%get(node)
          heatSource                                         => heatSource%next
       end do
       !![
       <referenceConstruct object="massDistributionHeating_">
	 <constructor>
           massDistributionHeatingSummation(                                                  &amp;
            &amp;                           massDistributionHeatings=massDistributionHeatings &amp;
            &amp;                          )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function summationGet

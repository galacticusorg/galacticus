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
  Implementation of a star formation rate in galactic disks which computes the rate by dividing the available gas mass by a
  timescale.
  !!}

  use :: Star_Formation_Timescales   , only : starFormationTimescaleClass
  use :: Star_Formation_Active_Masses, only : starFormationActiveMassClass
  
  !![
  <starFormationRateDisks name="starFormationRateDisksTimescale">
   <description>A star formation rate in galactic disks which computes the rate by integrating a star formation rate over the disk.</description>
  </starFormationRateDisks>
  !!]
  type, extends(starFormationRateDisksClass) :: starFormationRateDisksTimescale
     !!{
     Implementation of a rate for star formation in galactic disks which computes the rate by integrating a star formation rate
     over the disk.
     !!}
     private
     class(starFormationTimescaleClass ), pointer :: starFormationTimescale_  => null()
     class(starFormationActiveMassClass), pointer :: starFormationActiveMass_ => null()
   contains
     final     ::         timescaleDestructor
     procedure :: rate => timescaleRate
  end type starFormationRateDisksTimescale

  interface starFormationRateDisksTimescale
     !!{
     Constructors for the \refClass{starFormationRateDisksTimescale} star formation rate in disks class.
     !!}
     module procedure timescaleConstructorParameters
     module procedure timescaleConstructorInternal
  end interface starFormationRateDisksTimescale

contains

  function timescaleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationRateDisksTimescale} star formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (starFormationRateDisksTimescale)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(starFormationTimescaleClass    ), pointer       :: starFormationTimescale_
    class(starFormationActiveMassClass   ), pointer       :: starFormationActiveMass_

    !![
    <objectBuilder class="starFormationTimescale"  name="starFormationTimescale_"  source="parameters"/>
    <objectBuilder class="starFormationActiveMass" name="starFormationActiveMass_" source="parameters"/>
    !!]
    self=starFormationRateDisksTimescale(starFormationActiveMass_,starFormationTimescale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationTimescale_" />
    <objectDestructor name="starFormationActiveMass_"/>
    !!]
    return
  end function timescaleConstructorParameters

  function timescaleConstructorInternal(starFormationActiveMass_,starFormationTimescale_) result(self)
    !!{
    Internal constructor for the \refClass{starFormationRateDisksTimescale} star formation rate in disks class.
    !!}
    implicit none
    type (starFormationRateDisksTimescale)                        :: self
    class(starFormationTimescaleClass    ), intent(in   ), target :: starFormationTimescale_
    class(starFormationActiveMassClass   ), intent(in   ), target :: starFormationActiveMass_
    !![
    <constructorAssign variables="*starFormationActiveMass_, *starFormationTimescale_"/>
    !!]

    return
  end function timescaleConstructorInternal

  subroutine timescaleDestructor(self)
    !!{
    Destructor for the \refClass{starFormationRateDisksTimescale} star formation rate in disks class.
    !!}
    implicit none
    type(starFormationRateDisksTimescale), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationTimescale_" />
    <objectDestructor name="self%starFormationActiveMass_"/>
    !!]
    return
  end subroutine timescaleDestructor

  double precision function timescaleRate(self,node)
    !!{
    Returns the star formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily node}, by
    dividing the available gas mass by a timescale.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class           (starFormationRateDisksTimescale), intent(inout), target  :: self
    type            (treeNode                       ), intent(inout), target  :: node
    class           (nodeComponentDisk              )               , pointer :: disk
    double precision                                                          :: timescale

    disk      => node%disk                              (    )
    timescale =  self%starFormationTimescale_ %timescale(disk)
    if (timescale > 0.0d0) then
       timescaleRate=+self%starFormationActiveMass_%massActive(disk) &
            &        /                              timescale
    else
       ! Timescale is unphysical (presumably due to an unphysical condition in the disk), assume a zero rate.
       timescaleRate=+0.0d0
    end if
    return
  end function timescaleRate

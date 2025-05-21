!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implementation of a star formation rate in galactic disks which sets rates in satellites to zero.
  !!}

  !![
  <starFormationRateDisks name="starFormationRateDisksCentralsOnly">
   <description>A star formation rate in galactic disks which sets rates in satellites to zero.</description>
  </starFormationRateDisks>
  !!]
  type, extends(starFormationRateDisksClass) :: starFormationRateDisksCentralsOnly
     !!{
     Implementation of a rate for star formation in galactic disks which computes the rate by integrating a star formation rate
     over the disk.
     !!}
     private
     class(starFormationRateDisksClass), pointer :: starFormationRateDisks_ => null()
   contains
     final     ::         centralsOnlyDestructor
     procedure :: rate => centralsOnlyRate
  end type starFormationRateDisksCentralsOnly

  interface starFormationRateDisksCentralsOnly
     !!{
     Constructors for the \refClass{starFormationRateDisksCentralsOnly} star formation rate in disks class.
     !!}
     module procedure centralsOnlyConstructorParameters
     module procedure centralsOnlyConstructorInternal
  end interface starFormationRateDisksCentralsOnly

contains

  function centralsOnlyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationRateDisksCentralsOnly} star formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (starFormationRateDisksCentralsOnly)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(starFormationRateDisksClass       ), pointer       :: starFormationRateDisks_

    !![
    <objectBuilder class="starFormationRateDisks" name="starFormationRateDisks_" source="parameters"/>
    !!]
    self=starFormationRateDisksCentralsOnly(starFormationRateDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"/>
    !!]
    return
  end function centralsOnlyConstructorParameters

  function centralsOnlyConstructorInternal(starFormationRateDisks_) result(self)
    !!{
    Internal constructor for the \refClass{starFormationRateDisksCentralsOnly} star formation rate in disks class.
    !!}
    implicit none
    type (starFormationRateDisksCentralsOnly)                        :: self
    class(starFormationRateDisksClass       ), intent(in   ), target :: starFormationRateDisks_
    !![
    <constructorAssign variables="*starFormationRateDisks_"/>
    !!]

    return
  end function centralsOnlyConstructorInternal

  subroutine centralsOnlyDestructor(self)
    !!{
    Destructor for the \refClass{starFormationRateDisksCentralsOnly} star formation rate in disks class.
    !!}
    implicit none
    type(starFormationRateDisksCentralsOnly), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_" />
    !!]
    return
  end subroutine centralsOnlyDestructor

  double precision function centralsOnlyRate(self,node)
    !!{
    Returns the star formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily
    node}. Assumes zero rate for satellites, falling through to another class for centrals.
    !!}
    implicit none
    class(starFormationRateDisksCentralsOnly), intent(inout), target :: self
    type (treeNode                          ), intent(inout), target :: node
    
    if (node%isSatellite()) then
       centralsOnlyRate=0.0d0
    else
       centralsOnlyRate=self%starFormationRateDisks_%rate(node)
    end if
    return
  end function centralsOnlyRate

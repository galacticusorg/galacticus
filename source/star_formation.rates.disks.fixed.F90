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
  Implementation of a star formation rate in galactic disks which assumes a constant star formation rate.
  !!}

  !![
  <starFormationRateDisks name="starFormationRateDisksFixed">
   <description>A star formation rate in galactic disks which assumes a constant star formation rate.</description>
  </starFormationRateDisks>
  !!]
  type, extends(starFormationRateDisksClass) :: starFormationRateDisksFixed
     !!{
     Implementation of a rate for star formation in galactic disks which assumes a constant star formation rate.
     !!}
     private
     double precision :: rateStarFormation
   contains
     procedure :: rate => fixedRate
  end type starFormationRateDisksFixed

  interface starFormationRateDisksFixed
     !!{
     Constructors for the \refClass{starFormationRateDisksFixed} star formation rate in disks class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface starFormationRateDisksFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationRateDisksFixed} star formation rate in disks class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationRateDisksFixed)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: rateStarFormation
    
    !![
    <inputParameter>
      <name>rateStarFormation</name>
      <defaultValue>1.0d9</defaultValue>
      <description>The rate of star formation in units of $\mathrm{M}_\odot \hbox{Gyr}^{-1}$.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationRateDisksFixed(rateStarFormation)
    !![
    <inputParametersValidate source="parameters"/>
     !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(rateStarFormation) result(self)
    !!{
    Internal constructor for the \refClass{starFormationRateDisksFixed} star formation rate in disks class.
    !!}
    implicit none
    type            (starFormationRateDisksFixed)                :: self
    double precision                             , intent(in   ) :: rateStarFormation
    !![
    <constructorAssign variables="rateStarFormation"/>
    !!]

    return
  end function fixedConstructorInternal

  double precision function fixedRate(self,node)
    !!{
    Returns the star formation rate (in $\mathrm{M}_\odot$ Gyr$^{-1}$) in the galactic disk of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class(starFormationRateDisksFixed), intent(inout), target  :: self
    type (treeNode                   ), intent(inout), target  :: node

    fixedRate=self%rateStarFormation
    return
  end function fixedRate

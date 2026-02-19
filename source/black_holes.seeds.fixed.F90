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
  Implements fixed mass and spin black hole seeds.
  !!}

  !![
  <blackHoleSeeds name="blackHoleSeedsFixed">
   <description>
    A model of black hole seeds in which seeds have fixed mass and spin, independent of the halo in which they form.
   </description>
  </blackHoleSeeds>
  !!]
  type, extends(blackHoleSeedsClass) :: blackHoleSeedsFixed
     !!{
     A model of black hole seeds in which seeds have fixed mass and spin, independent of the halo in which they form.
     !!}
     private
     double precision :: mass_, spin_   
   contains
     procedure :: mass             => fixedMass
     procedure :: spin             => fixedSpin
     procedure :: formationChannel => fixedFormationChannel
  end type blackHoleSeedsFixed
  
  interface blackHoleSeedsFixed
     !!{
     Constructors for the \refClass{blackHoleSeedsFixed} black hole seeds class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface blackHoleSeedsFixed

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleSeedsFixed} black hole seeds class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleSeedsFixed)                :: self
    type            (inputParameters    ), intent(inout) :: parameters
    double precision                                     :: mass      , spin
    
    !![
    <inputParameter>
      <name>mass</name>
      <defaultValue>100.0d0</defaultValue>
      <description>The mass of seed black holes.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>spin</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The spin of seed black holes.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=blackHoleSeedsFixed(mass,spin)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(mass_,spin_) result(self)
    !!{
    Internal constructor for the \refClass{blackHoleSeedsFixed} black hole seed class.
    !!}
    implicit none
    type            (blackHoleSeedsFixed)                :: self
    double precision                     , intent(in   ) :: mass_, spin_
    !![
    <constructorAssign variables="mass_, spin_"/>
    !!]

    return
  end function standardConstructorInternal

  double precision function fixedMass(self,node) result(mass)
    !!{
    Compute the mass of the seed black hole.
    !!}
    implicit none
    class(blackHoleSeedsFixed), intent(inout) :: self
    type (treeNode           ), intent(inout) :: node
    !$GLC attributes unused :: node
    
    mass=self%mass_
    return
  end function fixedMass

  double precision function fixedSpin(self,node) result(spin)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    class(blackHoleSeedsFixed), intent(inout) :: self
    type (treeNode           ), intent(inout) :: node
    !$GLC attributes unused :: node
    
    spin=self%spin_
    return
  end function fixedSpin

  function fixedFormationChannel(self,node) result(channel)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    type (enumerationBlackHoleFormationChannelType)                :: channel
    class(blackHoleSeedsFixed                     ), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    channel=blackHoleFormationChannelUndetermined
    return
  end function fixedFormationChannel

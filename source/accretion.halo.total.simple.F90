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
  An implementation of the intergalactic medium state class for a simplistic model of instantaneous and full reionization.
  !!}

  !![
  <accretionHaloTotal name="accretionHaloTotalSimple">
   <description>
    A halo total accretion class which assumes that the accretion rate equals the growth rate of the basic mass.
   </description>
  </accretionHaloTotal>
  !!]
  type, extends(accretionHaloTotalClass) :: accretionHaloTotalSimple
     !!{
     A halo total accretion class which assumes the accretion corresponds to the basic mass.
     !!}
     private
   contains
     procedure :: accretionRate => simpleAccretionRate
     procedure :: accretedMass  => simpleAccretedMass
  end type accretionHaloTotalSimple

  interface accretionHaloTotalSimple
     !!{
     Constructors for the simple total halo accretion class.
     !!}
     module procedure simpleConstructorParameters
  end interface accretionHaloTotalSimple

contains

  function simpleConstructorParameters(parameters) result (self)
    !!{
    Constructor for the simple total halo accretion state class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(accretionHaloTotalSimple)                :: self
    type(inputParameters         ), intent(inout) :: parameters

    self=accretionHaloTotalSimple()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function simpleConstructorParameters

  double precision function simpleAccretionRate(self,node)
    !!{
    Return the accretion rate onto a halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(accretionHaloTotalSimple), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node
    class(nodeComponentBasic      ), pointer       :: basic
    !$GLC attributes unused :: self

    basic               => node %basic        ()
    simpleAccretionRate =  basic%accretionRate()
    return
  end function simpleAccretionRate

  double precision function simpleAccretedMass(self,node)
    !!{
    Return the mass accreted onto a halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(accretionHaloTotalSimple), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node
    class(nodeComponentBasic      ), pointer       :: basic
    !$GLC attributes unused :: self

    basic              => node %basic()
    simpleAccretedMass =  basic%mass ()
    return
  end function simpleAccretedMass

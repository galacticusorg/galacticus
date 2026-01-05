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
  <accretionHaloTotal name="accretionHaloTotalBertschinger">
   <description>
    A halo total accretion class that assumes that the accretion rate is equal to the growth rate of the basic component
    Bertschinger mass property.
   </description>
  </accretionHaloTotal>
  !!]
  type, extends(accretionHaloTotalClass) :: accretionHaloTotalBertschinger
     !!{
     A halo total accretion class which assumes the accretion corresponds to the Bertschinger mass.
     !!}
     private
     integer :: massBertschingerID, accretionRateBertschingerID
   contains
     procedure :: accretionRate => bertschingerAccretionRate
     procedure :: accretedMass  => bertschingerAccretedMass
  end type accretionHaloTotalBertschinger

  interface accretionHaloTotalBertschinger
     !!{
     Constructors for the bertschinger total halo accretion class.
     !!}
     module procedure bertschingerConstructorParameters
     module procedure bertschingerConstructorInternal
  end interface accretionHaloTotalBertschinger

contains

  function bertschingerConstructorParameters(parameters) result (self)
    !!{
    Constructor for the bertschinger total halo accretion state class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(accretionHaloTotalBertschinger)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=accretionHaloTotalBertschinger()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function bertschingerConstructorParameters

  function bertschingerConstructorInternal() result (self)
    !!{
    Internal constructor for the Bertschinger total halo accretion state class
    !!}
    implicit none
    type(accretionHaloTotalBertschinger) :: self

    !![
    <addMetaProperty component="basic" name="massBertschinger"          id="self%massBertschingerID"          isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="basic" name="accretionRateBertschinger" id="self%accretionRateBertschingerID" isEvolvable="no"  isCreator="no"/>
    !!]
    return
  end function bertschingerConstructorInternal

  double precision function bertschingerAccretionRate(self,node)
    !!{
    Return the accretion rate onto a halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(accretionHaloTotalBertschinger), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentBasic            ), pointer       :: basic
    !$GLC attributes unused :: self

    basic                     => node %basic                    (                                )
    bertschingerAccretionRate =  basic%floatRank0MetaPropertyGet(self%accretionRateBertschingerID)
    return
  end function bertschingerAccretionRate

  double precision function bertschingerAccretedMass(self,node)
    !!{
    Return the mass accreted onto a halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(accretionHaloTotalBertschinger), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentBasic            ), pointer       :: basic
    !$GLC attributes unused :: self

    basic                    => node %basic                    (                       )
    bertschingerAccretedMass =  basic%floatRank0MetaPropertyGet(self%massBertschingerID)
    return
  end function bertschingerAccretedMass

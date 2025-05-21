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
  Implements a node operator class that implements heating of the \gls{cgm} driven by accretion onto black holes.
  !!}

  use :: Black_Hole_CGM_Heating, only : blackHoleCGMHeatingClass

  !![
  <nodeOperator name="nodeOperatorBlackHolesCGMHeating">
   <description>A node operator class that implements heating of the \gls{cgm} driven by accretion onto black holes.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBlackHolesCGMHeating
     !!{
     A node operator class that implements heating of the \gls{cgm} driven by accretion onto black holes.
     !!}
     private
     class(blackHoleCGMHeatingClass), pointer :: blackHoleCGMHeating_ => null()
   contains
     final     ::                          blackHolesCGMHeatingDestructor
     procedure :: differentialEvolution => blackHolesCGMHeatingDifferentialEvolution
  end type nodeOperatorBlackHolesCGMHeating
  
  interface nodeOperatorBlackHolesCGMHeating
     !!{
     Constructors for the \refClass{nodeOperatorBlackHolesCGMHeating} node operator class.
     !!}
     module procedure blackHolesCGMHeatingConstructorParameters
     module procedure blackHolesCGMHeatingConstructorInternal
  end interface nodeOperatorBlackHolesCGMHeating
  
contains

  function blackHolesCGMHeatingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorBlackHolesCGMHeating} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBlackHolesCGMHeating)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(blackHoleCGMHeatingClass        ), pointer       :: blackHoleCGMHeating_
    
    !![
    <objectBuilder class="blackHoleCGMHeating" name="blackHoleCGMHeating_" source="parameters"/>
    !!]
    self=nodeOperatorBlackHolesCGMHeating(blackHoleCGMHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleCGMHeating_"/>
    !!]
    return
  end function blackHolesCGMHeatingConstructorParameters

  function blackHolesCGMHeatingConstructorInternal(blackHoleCGMHeating_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorBlackHolesCGMHeating} node operator class.
    !!}
    implicit none
    type (nodeOperatorBlackHolesCGMHeating)                        :: self
    class(blackHoleCGMHeatingClass        ), intent(in   ), target :: blackHoleCGMHeating_
    !![
    <constructorAssign variables="*blackHoleCGMHeating_"/>
    !!]
    
    return
  end function blackHolesCGMHeatingConstructorInternal

  subroutine blackHolesCGMHeatingDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorBlackHolesCGMHeating} node operator class.
    !!}
    implicit none
    type(nodeOperatorBlackHolesCGMHeating), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleCGMHeating_"/>
    !!]
    return
  end subroutine blackHolesCGMHeatingDestructor

  subroutine blackHolesCGMHeatingDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Account for winds driven by accretion onto black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, nodeComponentHotHalo
    implicit none
    class           (nodeOperatorBlackHolesCGMHeating), intent(inout), target  :: self
    type            (treeNode                        ), intent(inout), target  :: node
    logical                                           , intent(inout)          :: interrupt
    procedure       (interruptTask                   ), intent(inout), pointer :: functionInterrupt
    integer                                           , intent(in   )          :: propertyType
    class           (nodeComponentBlackHole          )               , pointer :: blackHole
    class           (nodeComponentHotHalo            )               , pointer :: hotHalo
    integer                                                                    :: countBlackHole   , indexBlackHole
    double precision                                                           :: rateHeating
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! If there are no black holes in this node, we have nothing to do - return immediately.
    countBlackHole=node%blackHoleCount()
    if (countBlackHole < 1) return
    ! Get the hot halo component so that heating energy can be injected into it.
    hotHalo => node%hotHalo()
    ! Iterate over all black holes in the node.
    do indexBlackHole=1,countBlackHole
       ! Find the wind power from this black hole.
       blackHole   => node%blackHole                       (instance=indexBlackHole)
       rateHeating =  self%blackHoleCGMHeating_%heatingRate(              blackHole)
       ! Inject this energy into the CGM.
       call hotHalo%heatSourceRate(rateHeating,interrupt,functionInterrupt)
    end do
    return
  end subroutine blackHolesCGMHeatingDifferentialEvolution
  

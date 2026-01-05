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
  Implements a node operator class that implements winds driven by accretion onto black holes.
  !!}

  use :: Black_Hole_Winds, only : blackHoleWindClass

  !![
  <nodeOperator name="nodeOperatorBlackHolesWinds">
   <description>A node operator class that implements winds driven by accretion onto black holes.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBlackHolesWinds
     !!{
     A node operator class that implements winds driven by accretion onto black holes.
     !!}
     private
     class(blackHoleWindClass), pointer :: blackHoleWind_ => null()
   contains
     final     ::                          blackHolesWindsDestructor
     procedure :: differentialEvolution => blackHolesWindsDifferentialEvolution
  end type nodeOperatorBlackHolesWinds
  
  interface nodeOperatorBlackHolesWinds
     !!{
     Constructors for the \refClass{nodeOperatorBlackHolesWinds} node operator class.
     !!}
     module procedure blackHolesWindsConstructorParameters
     module procedure blackHolesWindsConstructorInternal
  end interface nodeOperatorBlackHolesWinds
  
contains

  function blackHolesWindsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorBlackHolesWinds} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBlackHolesWinds)                :: self
    type (inputParameters            ), intent(inout) :: parameters
    class(blackHoleWindClass         ), pointer       :: blackHoleWind_
    
    !![
    <objectBuilder class="blackHoleWind" name="blackHoleWind_" source="parameters"/>
    !!]
    self=nodeOperatorBlackHolesWinds(blackHoleWind_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleWind_"/>
    !!]
    return
  end function blackHolesWindsConstructorParameters

  function blackHolesWindsConstructorInternal(blackHoleWind_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorBlackHolesWinds} node operator class.
    !!}
    implicit none
    type (nodeOperatorBlackHolesWinds)                        :: self
    class(blackHoleWindClass         ), intent(in   ), target :: blackHoleWind_
    !![
    <constructorAssign variables="*blackHoleWind_"/>
    !!]
    
    return
  end function blackHolesWindsConstructorInternal

  subroutine blackHolesWindsDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorBlackHolesWinds} node operator class.
    !!}
    implicit none
    type(nodeOperatorBlackHolesWinds), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleWind_"/>
    !!]
    return
  end subroutine blackHolesWindsDestructor

  subroutine blackHolesWindsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Account for winds driven by accretion onto black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorBlackHolesWinds), intent(inout), target  :: self
    type            (treeNode                   ), intent(inout), target  :: node
    logical                                      , intent(inout)          :: interrupt
    procedure       (interruptTask              ), intent(inout), pointer :: functionInterrupt
    integer                                      , intent(in   )          :: propertyType
    class           (nodeComponentBlackHole     )               , pointer :: blackHole
    class           (nodeComponentSpheroid      )               , pointer :: spheroid
    integer                                                               :: countBlackHole   , indexBlackHole
    double precision                                                      :: powerWind
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! If there are no black holes in this node, we have nothing to do - return immediately.
    countBlackHole=node%blackHoleCount()
    if (countBlackHole < 1) return
    ! Get the spheroid halo component so that wind energy can be injected into it.
    spheroid => node%spheroid()
    ! Iterate over all black holes in the node.
    do indexBlackHole=1,countBlackHole
       ! Find the wind power from this black hole.
       blackHole => node%blackHole           (instance=indexBlackHole)
       powerWind =  self%blackHoleWind_%power(              blackHole)
       ! Inject this energy into the spheroid.
       call spheroid%energyGasInputRate(+powerWind)
    end do
    return
  end subroutine blackHolesWindsDifferentialEvolution
  

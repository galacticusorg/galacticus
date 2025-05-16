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
  Implements a node operator class that computes the mass of gas which has cooled out of the \gls{cgm}.
  !!}
  
  use :: Cooling_Rates, only : coolingRateClass

  !![
  <nodeOperator name="nodeOperatorMassCooled">
    <description>
      A node operator class that computes the mass of gas which has cooled out of the \gls{cgm}.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorMassCooled
     !!{
     A node operator class that computes the mass of gas which has cooled out of the \gls{cgm}.
     !!}
     private
     class  (coolingRateClass), pointer :: coolingRate_ => null()
     integer                            :: massCooledID
   contains
     final     ::                                   massCooledDestructor
     procedure :: differentialEvolutionScales    => massCooledDifferentialEvolutionScales
     procedure :: differentialEvolutionInactives => massCooledDifferentialEvolutionInactives
     procedure :: differentialEvolution          => massCooledDifferentialEvolution
  end type nodeOperatorMassCooled
  
  interface nodeOperatorMassCooled
     !!{
     Constructors for the \refClass{nodeOperatorMassCooled} node operator class.
     !!}
     module procedure massCooledConstructorParameters
     module procedure massCooledConstructorInternal
  end interface nodeOperatorMassCooled
  
contains

  function massCooledConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorMassCooled} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorMassCooled)                :: self
    type (inputParameters       ), intent(inout) :: parameters
    class(coolingRateClass      ), pointer       :: coolingRate_
    
    !![
    <objectBuilder class="coolingRate" name="coolingRate_" source="parameters"/>
    !!]
    self=nodeOperatorMassCooled(coolingRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRate_"/>
    !!]
    return
  end function massCooledConstructorParameters

  function massCooledConstructorInternal(coolingRate_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorMassCooled} node operator class.
    !!}
    implicit none
    type (nodeOperatorMassCooled)                        :: self
    class(coolingRateClass      ), intent(in   ), target :: coolingRate_
    !![
    <constructorAssign variables="*coolingRate_"/>
    !!]
    
    !![
    <addMetaProperty component="hotHalo" name="massCooled" id="self%massCooledID" isEvolvable="yes" isCreator="yes"/>
    !!]
    return
  end function massCooledConstructorInternal

  subroutine massCooledDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorMassCooled} node operator class.
    !!}
    implicit none
    type(nodeOperatorMassCooled), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRate_"/>
    !!]
    return
  end subroutine massCooledDestructor

  subroutine massCooledDifferentialEvolutionInactives(self,node)
    !!{
    Mark \gls{cgm} mass cooled integrals as inactive for ODE solving.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class(nodeOperatorMassCooled), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node
    class(nodeComponentHotHalo  ), pointer       :: hotHalo
    
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Hot halo does not yet exist - nothing to do here.
    class default
       call hotHalo%floatRank0MetaPropertyInactive(self%massCooledID)
    end select
    return
  end subroutine massCooledDifferentialEvolutionInactives
  
  subroutine massCooledDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the mass cooled out of the \gls{cgm}.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, nodeComponentBasic
    implicit none
    class           (nodeOperatorMassCooled), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    class           (nodeComponentHotHalo  ), pointer       :: hotHalo
    class           (nodeComponentBasic    ), pointer       :: basic
    double precision                        , parameter     :: scaleMassRelative=1.0d-6

    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Hot halo does not yet exist - nothing to do here.
       class default
       basic => node%basic()
       call hotHalo%floatRank0MetaPropertyScale(self%massCooledID,scaleMassRelative*basic%mass())
    end select
    return
  end subroutine massCooledDifferentialEvolutionScales
  
  subroutine massCooledDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Integrates the mass cooling out of the \gls{cgm}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, propertyActive
    implicit none
    class    (nodeOperatorMassCooled), intent(inout), target  :: self
    type     (treeNode              ), intent(inout), target  :: node
    logical                          , intent(inout)          :: interrupt
    procedure(interruptTask         ), intent(inout), pointer :: functionInterrupt
    integer                          , intent(in   )          :: propertyType
    class    (nodeComponentHotHalo  )               , pointer :: hotHalo
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! Return immediately if active variables are requested.
    if (propertyActive(propertyType)) return
    ! Compute the cooling rate.
    hotHalo => node%hotHalo()
    ! Accumulate rates.
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Hot halo does not yet exist - nothing to do here.
    class default
       call hotHalo%floatRank0MetaPropertyRate(self%massCooledID,self%coolingRate_%rate(node))
    end select
    return
  end subroutine massCooledDifferentialEvolution

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
  Implements a node operator class that implements starvation of gas from the \gls{cgm} of satellite halos.
  !!}

  !![
  <nodeOperator name="nodeOperatorCGMStarvation">
   <description>
    A node operator class that implements starvation of gas from the \gls{cgm} of satellite halos.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMStarvation
     !!{
     A node operator class that implements starvation of gas from the \gls{cgm} of satellite halos.
     !!}
     private
   contains
     procedure :: differentialEvolutionPost => cgmStarvationDifferentialEvolutionPost
  end type nodeOperatorCGMStarvation
  
  interface nodeOperatorCGMStarvation
     !!{
     Constructors for the \refClass{nodeOperatorCGMStarvation} node operator class.
     !!}
     module procedure cgmStarvationConstructorParameters
  end interface nodeOperatorCGMStarvation
  
contains
  
  function cgmStarvationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMStarvation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorCGMStarvation)                :: self
    type(inputParameters          ), intent(inout) :: parameters
    
    self=nodeOperatorCGMStarvation()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cgmStarvationConstructorParameters

  subroutine cgmStarvationDifferentialEvolutionPost(self,node)
    !!{
    Perform starvation of \gls{cgm} gas in satellites.
    !!}
    use :: Abundances_Structure, only : abundances          , zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo
    implicit none
    class(nodeOperatorCGMStarvation), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    type (treeNode                 ), pointer       :: nodeParent
    class(nodeComponentHotHalo     ), pointer       :: hotHalo   , hotHaloParent

    if (.not.node%isSatellite()) return
    ! Transfer any CGM gas to the CGM of the parent node.
    nodeParent => node%parent
    do while (nodeParent%isSatellite())
       nodeParent => nodeParent%parent
    end do
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Generic type - nothing to do.
    class default
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       call    hotHaloParent%               massSet(hotHaloParent%         mass      ()+hotHalo%         mass      ())
       call    hotHaloParent%         abundancesSet(hotHaloParent%         abundances()+hotHalo%         abundances())
       call    hotHalo      %               massSet(                                                            0.0d0)
       call    hotHalo      %         abundancesSet(                                                   zeroAbundances)
       if (hotHalo%outflowedMassIsSettable()) then
          call hotHaloParent%outflowedMassSet      (hotHaloParent%outflowedMass      ()+hotHalo%outflowedMass      ())
          call hotHaloParent%outflowedAbundancesSet(hotHaloParent%outflowedAbundances()+hotHalo%outflowedAbundances())
          call hotHalo      %outflowedMassSet      (                                                            0.0d0)
          call hotHalo      %outflowedAbundancesSet(                                                   zeroAbundances)
       end if
    end select
    return
  end subroutine cgmStarvationDifferentialEvolutionPost
  

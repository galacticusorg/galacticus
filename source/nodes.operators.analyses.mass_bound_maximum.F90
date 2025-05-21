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
  Implements a node operator class that tracks the maximum bound mass achieved by a node.
  !!}
  
  !![
  <nodeOperator name="nodeOperatorMassBoundMaximum">
    <description>
      A node operator class that tracks the maximum bound mass achieved by a node. Intended to be paired with the
      \refClass{nodePropertyExtractorMassBoundMaximum} class to extract these masses for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorMassBoundMaximum
     !!{
     A node operator class that tracks the maximum host halo mass which a node has occupied.
     !!}
     private
     integer :: massBoundMaximumID
   contains
     !![
     <methods>
       <method method="update" description="Update the maximum bound mass of this node."/>
     </methods>
     !!]
     procedure :: nodeInitialize            => massBoundMaximumNodeInitialize
     procedure :: update                    => massBoundMaximumUpdate
     procedure :: nodesMerge                => massBoundMaximumNodesMerge
     procedure :: differentialEvolutionPre  => massBoundMaximumDifferentialEvolutionPre
     procedure :: differentialEvolutionPost => massBoundMaximumDifferentialEvolutionPost
  end type nodeOperatorMassBoundMaximum
  
  interface nodeOperatorMassBoundMaximum
     !!{
     Constructors for the \refClass{nodeOperatorMassBoundMaximum} node operator class.
     !!}
     module procedure massBoundMaximumConstructorParameters
     module procedure massBoundMaximumConstructorInternal
  end interface nodeOperatorMassBoundMaximum
  
contains

  function massBoundMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorMassBoundMaximum} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorMassBoundMaximum)                :: self
    type(inputParameters            ), intent(inout) :: parameters
    
    self=nodeOperatorMassBoundMaximum()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massBoundMaximumConstructorParameters

  function massBoundMaximumConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorMassBoundMaximum} node operator class.
    !!}
    implicit none
    type(nodeOperatorMassBoundMaximum) :: self
    
    !![
    <addMetaProperty component="satellite" name="massBoundMaximum" id="self%massBoundMaximumID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function massBoundMaximumConstructorInternal

  subroutine massBoundMaximumUpdate(self,node)
    !!{
    Update the maximum bound mass of this node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(nodeOperatorMassBoundMaximum), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    class(nodeComponentSatellite      ), pointer       :: satellite
    
    satellite => node%satellite()
    select type (satellite)
    type is (nodeComponentSatellite)
       ! No satellite exists, so nothing to do.
    class default
       call satellite%floatRank0MetaPropertySet(                                                                   &
            &                                                                            self%massBoundMaximumID , &
            &                                    max(                                                              &
            &                                        satellite%floatRank0MetaPropertyGet(self%massBoundMaximumID), &
            &                                        satellite%boundMass                (                       )  &
            &                                       )                                                              &
            &                                   )
    end select
    return
  end subroutine massBoundMaximumUpdate

  subroutine massBoundMaximumNodeInitialize(self,node)
    !!{
    Initialize the maximum bound mass of this node.
    !!}
    implicit none
    class(nodeOperatorMassBoundMaximum), intent(inout), target :: self
    type (treeNode                    ), intent(inout), target :: node
    
    call self%update(node)
    return
  end subroutine massBoundMaximumNodeInitialize

  subroutine massBoundMaximumDifferentialEvolutionPre(self,node)
    !!{
    Update the maximum bound mass of this node.
    !!}
    implicit none
    class(nodeOperatorMassBoundMaximum), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    
    call self%update(node) 
    return
  end subroutine massBoundMaximumDifferentialEvolutionPre

  subroutine massBoundMaximumDifferentialEvolutionPost(self,node)
    !!{
    Update the maximum bound mass of this node.
    !!}
    implicit none
    class(nodeOperatorMassBoundMaximum), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    
    call self%update(node) 
    return
  end subroutine massBoundMaximumDifferentialEvolutionPost

  subroutine massBoundMaximumNodesMerge(self,node)
    !!{
    Update the maximum bound mass of this node.
    !!}
    implicit none
    class(nodeOperatorMassBoundMaximum), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    
    call self%update(node) 
    return
  end subroutine massBoundMaximumNodesMerge

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

  !!{RST
  Implements a node operator class that performs tidal mass loss in spheroids.
  !!}

  use :: Tidal_Stripping_Mass_Loss_Rate, only : tidalStrippingClass
  
  !![
  <nodeOperator name="nodeOperatorTidalMassLossSpheroids" docformat="rst">
   <description>
   Computes and applies gravitational tidal stripping of stellar mass from satellite galaxy spheroids as they orbit within the host halo tidal field, progressively removing outer stars from the bulge component at each orbit.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTidalMassLossSpheroids
     !!{RST
     A node operator class that performs tidal mass loss in spheroids.
     !!}
     private
     class(tidalStrippingClass), pointer :: tidalStripping_ => null()
   contains
     final     ::                          tidalMassLossSpheroidsDestructor
     procedure :: differentialEvolution => tidalMassLossSpheroidsDifferentialEvolution
  end type nodeOperatorTidalMassLossSpheroids
  
  interface nodeOperatorTidalMassLossSpheroids
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorTidalMassLossSpheroids` node operator class.
     !!}
     module procedure tidalMassLossSpheroidsConstructorParameters
     module procedure tidalMassLossSpheroidsConstructorInternal
  end interface nodeOperatorTidalMassLossSpheroids
  
contains

  function tidalMassLossSpheroidsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorTidalMassLossSpheroids` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorTidalMassLossSpheroids)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(tidalStrippingClass               ), pointer       :: tidalStripping_
    
    !![
    <objectBuilder class="tidalStripping" name="tidalStripping_" source="parameters"/>
    !!]
    self=nodeOperatorTidalMassLossSpheroids(tidalStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="tidalStripping_"/>
    !!]
    return
  end function tidalMassLossSpheroidsConstructorParameters

  function tidalMassLossSpheroidsConstructorInternal(tidalStripping_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorTidalMassLossSpheroids` node operator class.
    !!}
    implicit none
    type (nodeOperatorTidalMassLossSpheroids)                        :: self
    class(tidalStrippingClass               ), intent(in   ), target :: tidalStripping_
    !![
    <constructorAssign variables="*tidalStripping_"/>
    !!]

    return
  end function tidalMassLossSpheroidsConstructorInternal

  subroutine tidalMassLossSpheroidsDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorTidalMassLossSpheroids` node operator class.
    !!}
    implicit none
    type(nodeOperatorTidalMassLossSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%tidalStripping_"/>
    !!]
    return
  end subroutine tidalMassLossSpheroidsDestructor
  
  subroutine tidalMassLossSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Apply tidal mass loss to the spheroid.
    !!}
    use :: Galacticus_Nodes         , only : propertyInactive
    use :: Tidal_Mass_Loss_Utilities, only : Tidal_Mass_Loss_Apply_Spheroid
    implicit none
    class    (nodeOperatorTidalMassLossSpheroids), intent(inout), target  :: self
    type     (treeNode                          ), intent(inout), target  :: node
    logical                                      , intent(inout)          :: interrupt
    procedure(interruptTask                     ), intent(inout), pointer :: functionInterrupt
    integer                                      , intent(in   )          :: propertyType
    !$GLC attributes unused :: interrupt, functionInterrupt

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)) return
    call Tidal_Mass_Loss_Apply_Spheroid(node,self%tidalStripping_)
    return
  end subroutine tidalMassLossSpheroidsDifferentialEvolution


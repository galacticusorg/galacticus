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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a node operator class that performs tidal mass loss in nuclear star clusters.
  !!}

  use :: Tidal_Stripping_Mass_Loss_Rate, only : tidalStrippingClass
  
  !![
  <nodeOperator name="nodeOperatorTidalMassLossNuclearStarClusters" docformat="rst">
   <description>
   Computes and applies gravitational tidal stripping of stellar and gaseous mass from the nuclear star clusters of satellite galaxies as they orbit within the host halo tidal field (see :galacticus-class:`tidalStripping`). Stripped gas is moved to the outflowing reservoir of the satellite's hot halo, while stripped stars are removed.
   The degree of tidal truncation of the nuclear star cluster, :math:`x=-\log f` where :math:`f` is the fraction of the mass of its untruncated mass distribution which it retains, is tracked in its ``tidalTruncation`` meta-property (see :galacticus-class:`galacticComponentTruncation`).
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTidalMassLossNuclearStarClusters
     !!{RST
     A node operator class that performs tidal mass loss in nuclear star clusters.
     !!}
     private
     class  (tidalStrippingClass), pointer :: tidalStripping_   => null()
     integer                             :: tidalTruncationID
   contains
     final     ::                                tidalMassLossNuclearStarClustersDestructor
     procedure :: differentialEvolution       => tidalMassLossNuclearStarClustersDifferentialEvolution
     procedure :: differentialEvolutionScales => tidalMassLossNuclearStarClustersDifferentialEvolutionScales
  end type nodeOperatorTidalMassLossNuclearStarClusters
  
  interface nodeOperatorTidalMassLossNuclearStarClusters
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorTidalMassLossNuclearStarClusters` node operator class.
     !!}
     module procedure tidalMassLossNuclearStarClustersConstructorParameters
     module procedure tidalMassLossNuclearStarClustersConstructorInternal
  end interface nodeOperatorTidalMassLossNuclearStarClusters
  
contains

  function tidalMassLossNuclearStarClustersConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorTidalMassLossNuclearStarClusters` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorTidalMassLossNuclearStarClusters)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(tidalStrippingClass                         ), pointer       :: tidalStripping_
    
    !![
    <objectBuilder class="tidalStripping" name="tidalStripping_" source="parameters"/>
    !!]
    self=nodeOperatorTidalMassLossNuclearStarClusters(tidalStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="tidalStripping_"/>
    !!]
    return
  end function tidalMassLossNuclearStarClustersConstructorParameters

  function tidalMassLossNuclearStarClustersConstructorInternal(tidalStripping_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorTidalMassLossNuclearStarClusters` node operator class.
    !!}
    implicit none
    type (nodeOperatorTidalMassLossNuclearStarClusters)                        :: self
    class(tidalStrippingClass                         ), intent(in   ), target :: tidalStripping_
    !![
    <constructorAssign variables="*tidalStripping_"/>
    <addMetaProperty component="NSC" name="tidalTruncation" id="self%tidalTruncationID" isEvolvable="yes" isCreator="yes"/>
    !!]

    return
  end function tidalMassLossNuclearStarClustersConstructorInternal

  subroutine tidalMassLossNuclearStarClustersDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorTidalMassLossNuclearStarClusters` node operator class.
    !!}
    implicit none
    type(nodeOperatorTidalMassLossNuclearStarClusters), intent(inout) :: self

    !![
    <objectDestructor name="self%tidalStripping_"/>
    !!]
    return
  end subroutine tidalMassLossNuclearStarClustersDestructor
  
  subroutine tidalMassLossNuclearStarClustersDifferentialEvolutionScales(self,node)
    !!{RST
    Set absolute ODE solver scale for the degree of tidal truncation of the nuclear star cluster.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class           (nodeOperatorTidalMassLossNuclearStarClusters), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    class           (nodeComponentNSC                            ), pointer       :: component
    double precision                                               , parameter     :: scaleTruncation=1.0d-2

    component => node%NSC()
    select type (component)
    type is (nodeComponentNSC)
       ! Component does not yet exist - nothing to do here.
    class default
       call component%floatRank0MetaPropertyScale(self%tidalTruncationID,scaleTruncation)
    end select
    return
  end subroutine tidalMassLossNuclearStarClustersDifferentialEvolutionScales

  subroutine tidalMassLossNuclearStarClustersDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Apply tidal mass loss to the nuclear star cluster.
    !!}
    use :: Galacticus_Nodes         , only : propertyInactive
    use :: Tidal_Mass_Loss_Utilities, only : Tidal_Mass_Loss_Apply_NSC
    implicit none
    class    (nodeOperatorTidalMassLossNuclearStarClusters), intent(inout), target  :: self
    type     (treeNode                                    ), intent(inout), target  :: node
    logical                                                , intent(inout)          :: interrupt
    procedure(interruptTask                               ), intent(inout), pointer :: functionInterrupt
    integer                                                , intent(in   )          :: propertyType
    !$GLC attributes unused :: interrupt, functionInterrupt

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)) return
    call Tidal_Mass_Loss_Apply_NSC(node,self%tidalStripping_,self%tidalTruncationID)
    return
  end subroutine tidalMassLossNuclearStarClustersDifferentialEvolution


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
  Implements a node operator class that performs tidal mass loss in disks.
  !!}

  use :: Galactic_Component_Truncations, only : galacticComponentTruncationClass
  use :: Tidal_Stripping_Mass_Loss_Rate, only : tidalStrippingClass
  
  !![
  <nodeOperator name="nodeOperatorTidalMassLossDisks" docformat="rst">
   <description>
   Computes and applies gravitational tidal stripping of stellar and gaseous mass from satellite galaxy disks as they orbit within the host halo tidal field, removing material outside the satellite's tidal radius at each pericentric passage.
   The degree of tidal truncation of the disk, :math:`x=-\log f` where :math:`f` is the fraction of the mass of its untruncated mass distribution which it retains, is tracked in its ``tidalTruncation`` meta-property (see :galacticus-class:`galacticComponentTruncation`).

   If the mass distribution of the disk is truncated (i.e. if ``[galacticComponentTruncation]`` is not ``none``), the stripped material carries the specific angular momentum implied by the change in the truncated mass distribution of the disk (see :galacticus-class:`galacticComponentTruncationExponential`), rather than that given by ``[tidalStripping]``. The specific angular momentum of the disk then remains consistent with its truncated mass distribution. (Otherwise, removing material from across the region beyond the tidal radius, while truncating the mass distribution from its edge, leaves the disk with more angular momentum than its mass distribution implies.) The ``[galacticComponentTruncation]`` used here must be the same as that used by the disk component, which is ensured by specifying it at the top level of the parameter file. (A fatal error occurs if one, but not the other, truncates the disk.)
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTidalMassLossDisks
     !!{RST
     A node operator class that performs tidal mass loss in disks.
     !!}
     private
     class  (tidalStrippingClass             ), pointer :: tidalStripping_              => null()
     class  (galacticComponentTruncationClass), pointer :: galacticComponentTruncation_ => null()
     integer                                            :: tidalTruncationID
   contains
     final     ::                                tidalMassLossDisksDestructor
     procedure :: differentialEvolution       => tidalMassLossDisksDifferentialEvolution
     procedure :: differentialEvolutionScales => tidalMassLossDisksDifferentialEvolutionScales
  end type nodeOperatorTidalMassLossDisks
  
  interface nodeOperatorTidalMassLossDisks
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorTidalMassLossDisks` node operator class.
     !!}
     module procedure tidalMassLossDisksConstructorParameters
     module procedure tidalMassLossDisksConstructorInternal
  end interface nodeOperatorTidalMassLossDisks
  
contains

  function tidalMassLossDisksConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorTidalMassLossDisks` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorTidalMassLossDisks  )                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(tidalStrippingClass             ), pointer       :: tidalStripping_
    class(galacticComponentTruncationClass), pointer       :: galacticComponentTruncation_
    
    !![
    <objectBuilder class="tidalStripping"              name="tidalStripping_"              source="parameters"/>
    <objectBuilder class="galacticComponentTruncation" name="galacticComponentTruncation_" source="parameters"/>
    !!]
    self=nodeOperatorTidalMassLossDisks(tidalStripping_,galacticComponentTruncation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="tidalStripping_"             />
    <objectDestructor name="galacticComponentTruncation_"/>
    !!]
    return
  end function tidalMassLossDisksConstructorParameters

  function tidalMassLossDisksConstructorInternal(tidalStripping_,galacticComponentTruncation_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorTidalMassLossDisks` node operator class.
    !!}
    implicit none
    type (nodeOperatorTidalMassLossDisks  )                        :: self
    class(tidalStrippingClass             ), intent(in   ), target :: tidalStripping_
    class(galacticComponentTruncationClass), intent(in   ), target :: galacticComponentTruncation_
    !![
    <constructorAssign variables="*tidalStripping_, *galacticComponentTruncation_"/>
    <addMetaProperty component="disk" name="tidalTruncation" id="self%tidalTruncationID" isEvolvable="yes" isCreator="yes"/>
    !!]

    return
  end function tidalMassLossDisksConstructorInternal

  subroutine tidalMassLossDisksDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorTidalMassLossDisks` node operator class.
    !!}
    implicit none
    type(nodeOperatorTidalMassLossDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%tidalStripping_"             />
    <objectDestructor name="self%galacticComponentTruncation_"/>
    !!]
    return
  end subroutine tidalMassLossDisksDestructor
  
  subroutine tidalMassLossDisksDifferentialEvolutionScales(self,node)
    !!{RST
    Set absolute ODE solver scale for the degree of tidal truncation of the disk.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class           (nodeOperatorTidalMassLossDisks), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    class           (nodeComponentDisk             ), pointer       :: component
    double precision                                , parameter     :: scaleTruncation=1.0d-2

    component => node%disk()
    select type (component)
    type is (nodeComponentDisk)
       ! Component does not yet exist - nothing to do here.
    class default
       call component%floatRank0MetaPropertyScale(self%tidalTruncationID,scaleTruncation)
    end select
    return
  end subroutine tidalMassLossDisksDifferentialEvolutionScales

  subroutine tidalMassLossDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Apply tidal mass loss to the disk. If the mass distribution of the disk is truncated, the stripped material carries the
    specific angular momentum implied by the change in the truncated mass distribution.
    !!}
    use :: Error                                    , only : Error_Report
    use :: Galactic_Component_Truncation_Disk_Status, only : diskTruncationIsActive
    use :: Galacticus_Nodes                         , only : propertyInactive          , nodeComponentDisk
    use :: Tidal_Mass_Loss_Utilities                , only : Tidal_Mass_Loss_Apply_Disk
    implicit none
    class           (nodeOperatorTidalMassLossDisks), intent(inout), target  :: self
    type            (treeNode                      ), intent(inout), target  :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (interruptTask                 ), intent(inout), pointer :: functionInterrupt
    integer                                         , intent(in   )          :: propertyType
    class           (nodeComponentDisk             )               , pointer :: disk
    double precision                                                         :: ratioSpecificAngularMomentum
    !$GLC attributes unused :: interrupt, functionInterrupt

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)) return
    ! The truncation used here must match that of the disk, as otherwise the angular momentum of the disk becomes inconsistent with
    ! its mass distribution.
    if (self%galacticComponentTruncation_%isActive() .neqv. diskTruncationIsActive)                                       &
         & call Error_Report(                                                                                             &
         &                   'the disk is truncated by the [galacticComponentTruncation] used by only one of the disk '// &
         &                   'component and this operator - declare it at the top level of the parameter file'         // &
         &                   {introspection:location}                                                                     &
         &                  )
    if (self%galacticComponentTruncation_%isActive()) then
       disk => node%disk()
       select type (disk)
       type is (nodeComponentDisk)
          ! Disk does not yet exist - nothing to do here.
          return
       end select
       ratioSpecificAngularMomentum=self%galacticComponentTruncation_%ratioAngularMomentumCylindrical(exp(-disk%floatRank0MetaPropertyGet(self%tidalTruncationID)))
       call Tidal_Mass_Loss_Apply_Disk(node,self%tidalStripping_,self%tidalTruncationID,ratioSpecificAngularMomentum)
    else
       call Tidal_Mass_Loss_Apply_Disk(node,self%tidalStripping_,self%tidalTruncationID                             )
    end if
    return
  end subroutine tidalMassLossDisksDifferentialEvolution


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
  Implements a node operator class that performs ram pressure mass loss in disks.
  !!}

  use :: Galactic_Component_Truncations       , only : galacticComponentTruncationClass
  use :: Ram_Pressure_Stripping_Mass_Loss_Rate, only : ramPressureStrippingClass
  
  !![
  <nodeOperator name="nodeOperatorRamPressureMassLossDisks" docformat="rst">
   <description>
   Applies ram pressure stripping of gas from galactic disks. The mass loss rate is computed by the selected
   :galacticus-class:`ramPressureStrippingClass` object; this operator removes that gas from the disk---together with angular
   momentum and abundances---and pipes it into the :term:`CGM` as an outflow.

   If the mass distribution of the disk is not truncated (i.e. if ``[galacticComponentTruncation]`` is ``none``), angular momentum
   is removed in proportion to mass, such that the shape of the disk is unchanged. Otherwise, ram pressure stripping truncates the
   gas of the disk (but not its stars): the degree of truncation of the gas, :math:`x_\mathrm{gas}=-\log f_\mathrm{gas}` where
   :math:`f_\mathrm{gas}` is the fraction of the mass of its untruncated mass distribution which it retains (tracked in the
   ``truncationGas`` meta-property of the disk, see :galacticus-class:`galacticComponentTruncation`), increases at the rate
   :math:`\dot{x}_\mathrm{gas} = \dot{M}_\mathrm{strip}/M_\mathrm{gas}`, and the stripped gas carries the specific angular
   momentum implied by the change in its truncated mass distribution (see
   :galacticus-class:`galacticComponentTruncationExponential`). Assuming a flat rotation curve, the mean specific angular momentum
   of the gas is that of the disk scaled by the ratio of the mean radius of the gas to the mass-weighted mean radius of the stars
   and gas. The ``[galacticComponentTruncation]`` used here must be the same as that used by the disk component, which is ensured
   by specifying it at the top level of the parameter file. (A fatal error occurs if one, but not the other, truncates the disk.)
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorRamPressureMassLossDisks
     !!{RST
     A node operator class that performs ram pressure mass loss in disks.
     !!}
     private
     class  (ramPressureStrippingClass       ), pointer :: ramPressureStripping_        => null()
     class  (galacticComponentTruncationClass), pointer :: galacticComponentTruncation_ => null()
     integer                                            :: tidalTruncationID                     , truncationGasID
   contains
     final     ::                                ramPressureMassLossDisksDestructor
     procedure :: differentialEvolution       => ramPressureMassLossDisksDifferentialEvolution
     procedure :: differentialEvolutionScales => ramPressureMassLossDisksDifferentialEvolutionScales
  end type nodeOperatorRamPressureMassLossDisks
  
  interface nodeOperatorRamPressureMassLossDisks
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorRamPressureMassLossDisks` node operator class.
     !!}
     module procedure ramPressureMassLossDisksConstructorParameters
     module procedure ramPressureMassLossDisksConstructorInternal
  end interface nodeOperatorRamPressureMassLossDisks
  
contains

  function ramPressureMassLossDisksConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorRamPressureMassLossDisks` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorRamPressureMassLossDisks)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(ramPressureStrippingClass           ), pointer       :: ramPressureStripping_
    class(galacticComponentTruncationClass    ), pointer       :: galacticComponentTruncation_
    
    !![
    <objectBuilder class="ramPressureStripping"        name="ramPressureStripping_"        source="parameters"/>
    <objectBuilder class="galacticComponentTruncation" name="galacticComponentTruncation_" source="parameters"/>
    !!]
    self=nodeOperatorRamPressureMassLossDisks(ramPressureStripping_,galacticComponentTruncation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="ramPressureStripping_"       />
    <objectDestructor name="galacticComponentTruncation_"/>
    !!]
    return
  end function ramPressureMassLossDisksConstructorParameters

  function ramPressureMassLossDisksConstructorInternal(ramPressureStripping_,galacticComponentTruncation_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorRamPressureMassLossDisks` node operator class.
    !!}
    implicit none
    type (nodeOperatorRamPressureMassLossDisks)                        :: self
    class(ramPressureStrippingClass           ), intent(in   ), target :: ramPressureStripping_
    class(galacticComponentTruncationClass    ), intent(in   ), target :: galacticComponentTruncation_
    !![
    <constructorAssign variables="*ramPressureStripping_, *galacticComponentTruncation_"/>
    !!]
    ! The degree of truncation of the gas is needed only if the disk is truncated. It is not registered otherwise since every
    ! registered evolvable meta-property adds a variable to the ODE system of every node. The degree of tidal truncation (of the
    ! stars) is needed to find the mean radius of the stars, but is not created by this class.
    if (self%galacticComponentTruncation_%isActive()) then
       !![
       <addMetaProperty component="disk" name="truncationGas"   id="self%truncationGasID"   isEvolvable="yes" isCreator="yes"/>
       <addMetaProperty component="disk" name="tidalTruncation" id="self%tidalTruncationID" isEvolvable="yes" isCreator="no" />
       !!]
    else
       self%truncationGasID  =-1
       self%tidalTruncationID=-1
    end if

    return
  end function ramPressureMassLossDisksConstructorInternal

  subroutine ramPressureMassLossDisksDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorRamPressureMassLossDisks` node operator class.
    !!}
    implicit none
    type(nodeOperatorRamPressureMassLossDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%ramPressureStripping_"       />
    <objectDestructor name="self%galacticComponentTruncation_"/>
    !!]
    return
  end subroutine ramPressureMassLossDisksDestructor

  subroutine ramPressureMassLossDisksDifferentialEvolutionScales(self,node)
    !!{RST
    Set absolute ODE solver scale for the degree of truncation of the gas of the disk.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class           (nodeOperatorRamPressureMassLossDisks), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    class           (nodeComponentDisk                   ), pointer       :: disk
    double precision                                      , parameter     :: scaleTruncation=1.0d-2

    if (.not.self%galacticComponentTruncation_%isActive()) return
    disk => node%disk()
    select type (disk)
    type is (nodeComponentDisk)
    class default
       call disk%floatRank0MetaPropertyScale(self%truncationGasID,scaleTruncation)
    end select
    return
  end subroutine ramPressureMassLossDisksDifferentialEvolutionScales
  
  subroutine ramPressureMassLossDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Perform ram pressure stripping-induced mass loss in a disk.
    !!}
    use :: Abundances_Structure                     , only : operator(*)
    use :: Error                                    , only : Error_Report
    use :: Galactic_Component_Truncation_Disk_Status, only : diskTruncationIsActive
    use :: Galacticus_Nodes                         , only : propertyInactive      , nodeComponentDisk, nodeComponentHotHalo
    implicit none
    class           (nodeOperatorRamPressureMassLossDisks), intent(inout), target  :: self
    type            (treeNode                            ), intent(inout), target  :: node
    logical                                               , intent(inout)          :: interrupt
    procedure       (interruptTask                       ), intent(inout), pointer :: functionInterrupt
    integer                                               , intent(in   )          :: propertyType
    class           (nodeComponentDisk                   )               , pointer :: disk
    class           (nodeComponentHotHalo                )               , pointer :: hotHalo
    double precision                                                               :: massLossRate           , ratioSpecificAngularMomentum, &
         &                                                                            fractionMassRetainedGas, fractionMassRetainedStellar , &
         &                                                                            ratioRadiusMeanGas     , ratioRadiusMeanStellar      , &
         &                                                                            massGas                , massStellar

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)        ) return
    if (self%galacticComponentTruncation_%isActive() .neqv. diskTruncationIsActive)                                       &
         & call Error_Report(                                                                                             &
         &                   'the disk is truncated by the [galacticComponentTruncation] used by only one of the disk '// &
         &                   'component and this operator - declare it at the top level of the parameter file'         // &
         &                   {introspection:location}                                                                     &
         &                  )
    disk => node%disk()
    if (disk%massGas() <= 0.0d0               ) return
    massLossRate=self%ramPressureStripping_%rateMassLoss(disk)
    if (massLossRate   <= 0.0d0               ) return
    hotHalo => node%hotHalo()
    if (self%galacticComponentTruncation_%isActive()) then
       ! The gas is truncated. Find the specific angular momentum of the stripped gas in units of the mean specific angular
       ! momentum of the disk.
       massGas                    =max(0.0d0,disk%massGas    ())
       massStellar                =max(0.0d0,disk%massStellar())
       fractionMassRetainedGas    =exp(-max(0.0d0,disk%floatRank0MetaPropertyGet(self%truncationGasID)))
       fractionMassRetainedStellar=1.0d0
       if (disk%floatRank0MetaPropertyIsCreated(self%tidalTruncationID)) &
            & fractionMassRetainedStellar=exp(-disk%floatRank0MetaPropertyGet(self%tidalTruncationID))
       ratioRadiusMeanGas          =+self%galacticComponentTruncation_%ratioRadiusMeanCylindrical     (fractionMassRetainedGas    )
       ratioRadiusMeanStellar      =+self%galacticComponentTruncation_%ratioRadiusMeanCylindrical     (fractionMassRetainedStellar)
       ratioSpecificAngularMomentum=+self%galacticComponentTruncation_%ratioAngularMomentumCylindrical(fractionMassRetainedGas    ) &
            &                       *ratioRadiusMeanGas                                                                             &
            &                       *(massStellar           +massGas                   )                                            &
            &                       /(massStellar*ratioRadiusMeanStellar+massGas*ratioRadiusMeanGas)
       ! The degree of truncation of the gas increases in proportion to the fractional rate of mass loss.
       call disk%floatRank0MetaPropertyRate(self%truncationGasID,+massLossRate/disk%massGas())
    else
       ! The shape of the disk is unchanged.
       ratioSpecificAngularMomentum=1.0d0
    end if
    call    disk%                  massGasRate(-massLossRate                                                                                        )
    call    disk%          angularMomentumRate(-massLossRate*disk%angularMomentum()/(disk%massGas()+disk%massStellar())*ratioSpecificAngularMomentum)
    call    disk%            abundancesGasRate(-massLossRate*disk%abundancesGas  ()/ disk%massGas()                                                 )
    call hotHalo%           outflowingMassRate(+massLossRate                                                                                        )
    call hotHalo%outflowingAngularMomentumRate(+massLossRate*disk%angularMomentum()/(disk%massGas()+disk%massStellar())*ratioSpecificAngularMomentum)
    call hotHalo%outflowingAbundancesRate     (+massLossRate*disk%abundancesGas  ()/ disk%massGas()                                                 )
    return
  end subroutine ramPressureMassLossDisksDifferentialEvolution


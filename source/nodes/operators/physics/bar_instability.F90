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
  Implements a node operator class that implements bar instabilities in disks.
  !!}

  use :: Bar_Instability_Spheroid_Angular_Momenta, only : barInstabilitySpheroidAngularMomentumClass
  use :: Dark_Matter_Halo_Scales                 , only : darkMatterHaloScaleClass
  use :: Galactic_Dynamics_Bar_Instabilities     , only : galacticDynamicsBarInstabilityClass
  
  !![
  <nodeOperator name="nodeOperatorBarInstability" docformat="rst">
   <description>
   A node operator class that implements bar instabilities in disks. A fraction of the angular momentum of the material transferred from the disk to the spheroid is retained in the disk as suggested by numerical experiments :cite:p:`klypin_cdm-based_2002`.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBarInstability
     !!{RST
     A node operator class that implements bar instabilities in disks.
     !!}
     private
     class  (darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_            => null()
     class  (galacticDynamicsBarInstabilityClass), pointer :: galacticDynamicsBarInstability_ => null()
     class  (barInstabilitySpheroidAngularMomentumClass), pointer :: barInstabilitySpheroidAngularMomentum_ => null()
     logical                                               :: luminositiesStellarInactive
   contains
     final     ::                                   barInstabilityDestructor
     procedure :: differentialEvolutionAnalytics => barInstabilityDifferentialEvolutionAnalytics
     procedure :: differentialEvolution          => barInstabilityDifferentialEvolution
     procedure :: galacticStructureSolverSet      => barInstabilityGalacticStructureSolverSet
  end type nodeOperatorBarInstability

  interface nodeOperatorBarInstability
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorBarInstability` node operator class.
     !!}
     module procedure barInstabilityConstructorParameters
     module procedure barInstabilityConstructorInternal
  end interface nodeOperatorBarInstability
  
contains

  function barInstabilityConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorBarInstability` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorBarInstability         )                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    class  (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class  (galacticDynamicsBarInstabilityClass), pointer       :: galacticDynamicsBarInstability_
    class  (barInstabilitySpheroidAngularMomentumClass), pointer :: barInstabilitySpheroidAngularMomentum_
    logical                                                     :: luminositiesStellarInactive

    !![
    <inputParameter docformat="rst">
      <name>luminositiesStellarInactive</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>
      If true, stellar luminosities will be treated as inactive properties.
      </description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
    <objectBuilder class="galacticDynamicsBarInstability" name="galacticDynamicsBarInstability_" source="parameters"/>
    <objectBuilder class="barInstabilitySpheroidAngularMomentum" name="barInstabilitySpheroidAngularMomentum_" source="parameters"/>
    !!]
    self=nodeOperatorBarInstability(luminositiesStellarInactive,darkMatterHaloScale_,galacticDynamicsBarInstability_,barInstabilitySpheroidAngularMomentum_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticDynamicsBarInstability_"/>
    <objectDestructor name="barInstabilitySpheroidAngularMomentum_"/>
    <objectDestructor name="darkMatterHaloScale_"           />
    !!]
    return
  end function barInstabilityConstructorParameters

  function barInstabilityConstructorInternal(luminositiesStellarInactive,darkMatterHaloScale_,galacticDynamicsBarInstability_,barInstabilitySpheroidAngularMomentum_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorBarInstability` node operator class.
    !!}
    implicit none
    type   (nodeOperatorBarInstability         )                        :: self
    class  (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class  (galacticDynamicsBarInstabilityClass), intent(in   ), target :: galacticDynamicsBarInstability_
    class  (barInstabilitySpheroidAngularMomentumClass), intent(in   ), target :: barInstabilitySpheroidAngularMomentum_
    logical                                     , intent(in   )         :: luminositiesStellarInactive
    !![
    <constructorAssign variables="luminositiesStellarInactive, *darkMatterHaloScale_, *galacticDynamicsBarInstability_, *barInstabilitySpheroidAngularMomentum_"/>
    !!]

    return
  end function barInstabilityConstructorInternal

  subroutine barInstabilityDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorBarInstability` node operator class.
    !!}
    implicit none
    type(nodeOperatorBarInstability), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"           />
    <objectDestructor name="self%galacticDynamicsBarInstability_"/>
    <objectDestructor name="self%barInstabilitySpheroidAngularMomentum_"/>
    !!]
    return
  end subroutine barInstabilityDestructor

  subroutine barInstabilityGalacticStructureSolverSet(self,galacticStructureSolver)
    !!{RST
    Provide the active galactic structure solver to the spheroid pseudo-angular momentum model.
    !!}
    use :: Galactic_Structure_Solvers, only : galacticStructureSolverClass
    implicit none
    class(nodeOperatorBarInstability  ), intent(inout)         :: self
    class(galacticStructureSolverClass), intent(inout), target :: galacticStructureSolver

    call self%barInstabilitySpheroidAngularMomentum_%galacticStructureSolverSet(galacticStructureSolver)
    return
  end subroutine barInstabilityGalacticStructureSolverSet
  
  subroutine barInstabilityDifferentialEvolutionAnalytics(self,node)
    !!{RST
    Initialize the mass transfer fraction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class(nodeOperatorBarInstability), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentDisk         ), pointer       :: disk

    ! Mark the retained stellar mass fraction as analytically-solvable (it is always zero) if we are not solving for luminosities
    ! as inactive properties.
    if (.not.self%luminositiesStellarInactive) then
       disk => node%disk()
       select type (disk)
       type is (nodeComponentDisk)
          ! Disk does not yet exist.
       class default
          call disk%fractionMassRetainedAnalytic()
       end select
    end if
    return
  end subroutine barInstabilityDifferentialEvolutionAnalytics

  subroutine barInstabilityDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Implement the effects of global bar instability on the galactic disk.
    !!}
    use :: Galacticus_Nodes              , only : propertyInactive, nodeComponentDisk  , nodeComponentSpheroid
    use :: Abundances_Structure          , only : operator(*)     , abundances         , zeroAbundances         , max
    use :: Histories                     , only : operator(*)     , history
    use :: Stellar_Luminosities_Structure, only : operator(*)     , stellarLuminosities, zeroStellarLuminosities, max
    implicit none
    class           (nodeOperatorBarInstability), intent(inout), target  :: self
    type            (treeNode                  ), intent(inout), target  :: node
    logical                                     , intent(inout)          :: interrupt
    procedure       (interruptTask             ), intent(inout), pointer :: functionInterrupt
    integer                                     , intent(in   )          :: propertyType
    class           (nodeComponentDisk         )               , pointer :: disk
    class           (nodeComponentSpheroid     )               , pointer :: spheroid
    type            (abundances                ), save                   :: fuelAbundancesRates                , stellarAbundancesRates
    !$omp threadprivate(stellarAbundancesRates,fuelAbundancesRates)
    type            (stellarLuminosities       ), save                   :: luminositiesTransferRate
    !$omp threadprivate(luminositiesTransferRate)
    double precision                                                     :: barInstabilityTimescale                  , barInstabilitySpecificTorque           , &
         &                                                                  fractionAngularMomentumRetainedDisk      , fractionAngularMomentumRetainedSpheroid, &
         &                                                                  rateAngularMomentumDisk                  , rateAngularMomentumSpheroid            , &
         &                                                                  rateAngularMomentumTransfer              , rateMassGasTransfer                    , &
         &                                                                  rateMassStellarTransfer                  , transferRate
    type            (history                   )                         :: historyTransferRate

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)) return
    ! Check for a realistic disk, return immediately if disk is unphysical.
    disk => node%disk()
    if     (     disk%angularMomentum() < 0.0d0 &
         &  .or. disk%radius         () < 0.0d0 &
         &  .or. disk%massGas        () < 0.0d0 &
         & ) return
    ! Determine if the disk is bar unstable and, if so, the rate at which material is moved to the pseudo-bulge.
    if (node%isPhysicallyPlausible) then
       ! Disk has positive angular momentum, so compute an instability timescale.
       call self%galacticDynamicsBarInstability_%timescale(node,barInstabilityTimescale,barInstabilitySpecificTorque,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid)
    else
       ! Disk has non-positive angular momentum, therefore it is unphysical. Do not compute an instability timescale in this
       ! case as the disk radius may be unphysical also.
       barInstabilityTimescale=-1.0d0
    end if
    ! Negative timescale indicates no bar instability.
    if (barInstabilityTimescale < 0.0d0) return
    ! Disk is unstable, so compute rates at which material is transferred to the spheroid.
    spheroid => node%spheroid()
    ! Gas mass.
    rateMassGasTransfer        =max(         0.0d0         ,disk    %massGas             ())/barInstabilityTimescale
    call                                      disk    %massGasRate             (-                                                    rateMassGasTransfer                            )
    call                                      spheroid%massGasRate             (+                                                    rateMassGasTransfer,interrupt,functionInterrupt)
    ! Fraction of stellar mass transferred.
    if (self%luminositiesStellarInactive) then
       transferRate            =max(         0.0d0         ,disk    %fractionMassRetained())/barInstabilityTimescale
       call                                   disk    %fractionMassRetainedRate(-                                                            transferRate                            )
    end if
    ! Stellar mass.
    rateMassStellarTransfer    =max(         0.0d0         ,disk    %massStellar         ())/barInstabilityTimescale
    call                                      disk    %massStellarRate         (-                                                rateMassStellarTransfer                            )
    call                                      spheroid%massStellarRate         (+                                                rateMassStellarTransfer,interrupt,functionInterrupt)
    ! Angular momentum. We remove from the disk only its non-retained fraction. The selected spheroid model determines its
    ! pseudo-angular-momentum rate; the default retained model reproduces the historical retained-fraction prescription.
    rateAngularMomentumTransfer=max(         0.0d0         ,disk    %angularMomentum     ())/barInstabilityTimescale
    rateAngularMomentumDisk    =-(1.0d0-fractionAngularMomentumRetainedDisk)*rateAngularMomentumTransfer
    rateAngularMomentumSpheroid=self%barInstabilitySpheroidAngularMomentum_%rate(                  &
         &                                                                       node            , &
         &                                                                       rateMassGasTransfer, &
         &                                                                       rateMassStellarTransfer, &
         &                                                                       rateAngularMomentumDisk, &
         &                                                                       rateAngularMomentumTransfer, &
         &                                                                       fractionAngularMomentumRetainedSpheroid &
         &                                                                      )
    call                                      disk    %angularMomentumRate     (rateAngularMomentumDisk                                              )
    call                                      spheroid%angularMomentumRate     (rateAngularMomentumSpheroid,interrupt,functionInterrupt)
    ! Gas abundances.
    fuelAbundancesRates        =max(zeroAbundances         ,disk    %abundancesGas       ())/barInstabilityTimescale
    call                                      disk    %abundancesGasRate       (-                                                     fuelAbundancesRates                            )
    call                                      spheroid%abundancesGasRate       (+                                                     fuelAbundancesRates,interrupt,functionInterrupt)
    ! Stellar abundances.
    stellarAbundancesRates     =max(zeroAbundances         ,disk    %abundancesStellar   ())/barInstabilityTimescale
    call                                      disk    %abundancesStellarRate   (-                                                  stellarAbundancesRates                            )
    call                                      spheroid%abundancesStellarRate   (+                                                  stellarAbundancesRates,interrupt,functionInterrupt)
    ! Stellar luminosities.
    if (.not.self%luminositiesStellarInactive) then
       luminositiesTransferRate=max(zeroStellarLuminosities,disk    %luminositiesStellar ())/barInstabilityTimescale
       call                                   disk    %luminositiesStellarRate (-                                                luminositiesTransferRate                            )
       call                                   spheroid%luminositiesStellarRate (+                                                luminositiesTransferRate,interrupt,functionInterrupt)
    end if
    ! Stellar properties history.
    historyTransferRate=disk%stellarPropertiesHistory()
    if (historyTransferRate%exists()) then
       historyTransferRate=historyTransferRate/barInstabilityTimescale
       call disk    %stellarPropertiesHistoryRate(-historyTransferRate                             )
       call spheroid%stellarPropertiesHistoryRate(+historyTransferRate,interrupt,functionInterrupt)
    end if
    call historyTransferRate%destroy()
    ! Star formation history.
    historyTransferRate=disk%starFormationHistory()
    if (historyTransferRate%exists()) then
       historyTransferRate=historyTransferRate/barInstabilityTimescale
       call disk    %starFormationHistoryRate(-historyTransferRate                             )
       call spheroid%starFormationHistoryRate(+historyTransferRate,interrupt,functionInterrupt)
    end if
    call historyTransferRate%destroy()
    ! Additional external torque.
    if     (                                                                                                                                                                      &
         &   spheroid%angularMomentum() < (spheroid%massGas()+spheroid%massStellar())*self%darkMatterHaloScale_%radiusVirial(node)*self%darkMatterHaloScale_%velocityVirial(node) &
         &  .and.                                                                                                                                                                 &
         &   spheroid%radius         () <                                             self%darkMatterHaloScale_%radiusVirial(node)                                                &
         & ) then
       call spheroid%angularMomentumRate(+barInstabilitySpecificTorque*(spheroid%massGas()+spheroid%massStellar()),interrupt,functionInterrupt)
    end if
    return
  end subroutine barInstabilityDifferentialEvolution

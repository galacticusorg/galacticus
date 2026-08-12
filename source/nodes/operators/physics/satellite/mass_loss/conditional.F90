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
  Implements a node operator class that applies conditional mass loss to orbiting satellite halos.
  !!}

  use :: Dark_Matter_Profiles_Soliton_Status, only : enumerationSolitonStatusType      , solitonStatusSolitonNFW, solitonStatusSolitonOnly
  use :: Satellite_Tidal_Stripping          , only : satelliteTidalStrippingClass
  use :: Satellite_Tidal_Stripping_Radii    , only : satelliteTidalStrippingRadiusClass

  !![
  <nodeOperator name="nodeOperatorSatelliteConditionalMassLoss" docformat="rst">
   <description>
   A node operator class that applies conditional tidal mass loss to orbiting satellite halos at each ODE timestep for use with the :galacticus-class:`darkMatterProfileDMOSolitonNFWHeated` dark matter profile. Depending on the halo state, it evolves either the outer halo or the solitonic core using a :galacticus-class:`satelliteTidalStrippingClass`. If no solitonic solution exists, the halo is treated as an NFW halo and tidal stripping is applied using the NFW model throughout its evolution. The satellite bound mass and, where applicable, the solitonic core mass are evolved as ODE variables.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteConditionalMassLoss
     !!{RST
     A node operator class that applies tidal mass loss to orbiting satellite halos.
     !!}
     private
     integer                                              :: massCoreID                              , solitonStatusID                       , &
          &                                                  massCoreNormalID                        , radiusSolitonID
     class  (satelliteTidalStrippingClass      ), pointer :: satelliteTidalStrippingOuter_  => null(), satelliteTidalStrippingCore_ => null()
     class  (satelliteTidalStrippingRadiusClass), pointer :: satelliteTidalStrippingRadius_ => null()
   contains
     final     ::                          satelliteConditionalStrippingDestructor
     procedure :: differentialEvolution => satelliteConditionalStrippingDifferentialEvolution
  end type nodeOperatorSatelliteConditionalMassLoss

  interface nodeOperatorSatelliteConditionalMassLoss
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorSatelliteConditionalMassLoss` node operator class.
     !!}
     module procedure satelliteConditionalStrippingConstructorParameters
     module procedure satelliteConditionalStrippingConstructorInternal
  end interface nodeOperatorSatelliteConditionalMassLoss

  ! Submodule-scope pointer to self, used in callback functions.
  class          (nodeOperatorSatelliteConditionalMassLoss), pointer :: self_
  !$omp threadprivate(self_)

contains

  function satelliteConditionalStrippingConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorSatelliteConditionalMassLoss` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorSatelliteConditionalMassLoss)                :: self
    type   (inputParameters                         ), intent(inout) :: parameters
    class  (satelliteTidalStrippingClass            ), pointer       :: satelliteTidalStrippingOuter_ , satelliteTidalStrippingCore_
    class  (satelliteTidalStrippingRadiusClass      ), pointer       :: satelliteTidalStrippingRadius_

    ! The outer-halo stripping model is built from the unqualified `satelliteTidalStripping` parameter, so that it is inherited
    ! from the enclosing scope and therefore always matches the tidal stripping model used elsewhere in the model. Only the core
    ! stripping model, which is specific to the solitonic core, is named explicitly.
    !
    ! The tidal radius model is used only to detect the transition to the soliton-only state. It must be the *unlimited* model:
    ! if it were wrapped in `satelliteTidalStrippingRadiusLimited` the radius could never fall below the soliton radius, and the
    ! transition would never be detected.
    !![
    <objectBuilder class="satelliteTidalStripping"       name="satelliteTidalStrippingOuter_"                                              source="parameters"/>
    <objectBuilder class="satelliteTidalStripping"       name="satelliteTidalStrippingCore_"   parameterName="satelliteTidalStrippingCore" source="parameters"/>
    <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_"                                             source="parameters"/>
    !!]
    self=nodeOperatorSatelliteConditionalMassLoss(satelliteTidalStrippingOuter_,satelliteTidalStrippingCore_,satelliteTidalStrippingRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingOuter_" />
    <objectDestructor name="satelliteTidalStrippingCore_"  />
    <objectDestructor name="satelliteTidalStrippingRadius_"/>
    !!]
    return
  end function satelliteConditionalStrippingConstructorParameters

  function satelliteConditionalStrippingConstructorInternal(satelliteTidalStrippingOuter_,satelliteTidalStrippingCore_,satelliteTidalStrippingRadius_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorSatelliteConditionalMassLoss` node operator class.
    !!}
    implicit none
    type (nodeOperatorSatelliteConditionalMassLoss)                        :: self
    class(satelliteTidalStrippingClass            ), intent(in   ), target :: satelliteTidalStrippingOuter_ , satelliteTidalStrippingCore_
    class(satelliteTidalStrippingRadiusClass      ), intent(in   ), target :: satelliteTidalStrippingRadius_
    !![
    <constructorAssign variables="*satelliteTidalStrippingOuter_, *satelliteTidalStrippingCore_, *satelliteTidalStrippingRadius_"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton"  id="self%radiusSolitonID"  isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal" id="self%massCoreNormalID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonStatus"         id="self%solitonStatusID"  type="integer"    isCreator="no"/>
    !!]

    return
  end function satelliteConditionalStrippingConstructorInternal

  subroutine satelliteConditionalStrippingDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorSatelliteConditionalMassLoss` node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteConditionalMassLoss), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalStrippingOuter_" />
    <objectDestructor name="self%satelliteTidalStrippingCore_"  />
    <objectDestructor name="self%satelliteTidalStrippingRadius_"/>
    !!]
    return
  end subroutine satelliteConditionalStrippingDestructor

  subroutine satelliteConditionalStrippingDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Perform mass loss from a satellite due to tidal stripping.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentSatellite
    implicit none
    class    (nodeOperatorSatelliteConditionalMassLoss), intent(inout), target  :: self
    type     (treeNode                                ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure(interruptTask                           ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class    (nodeComponentSatellite                  )               , pointer :: satellite
    class    (nodeComponentDarkMatterProfile          ), pointer                :: darkMatterProfile
    double precision                                                            :: massLossRateOuter, massLossRateCore, &
              &                                                                    massCore         , radiusTidal     , &
              &                                                                    radiusSoliton
    type     (enumerationSolitonStatusType            )                         :: solitonStatus
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (.not.node%isSatellite()) return
    ! Get required quantities from the satellite node.
    satellite         => node%satellite        ()
    darkMatterProfile => node%darkMatterProfile()
    massCore          =  darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreID)
    solitonStatus     =  enumerationSolitonStatusType(darkMatterProfile%integerRank0MetaPropertyGet(self%solitonStatusID))
    ! Apply tidal mass loss according to the halo state:
    !  - solitonNFW : strip only the outer NFW component, until the tidal radius reaches the soliton radius.
    !  - solitonOnly: evolve both the bound mass and the solitonic core mass, and destroy the satellite if the core mass becomes non-positive.
    !  - otherwise  : no solitonic solution exists, so treat the halo as a pure NFW halo and evolve the NFW mass loss.
    if (solitonStatus == solitonStatusSolitonNFW) then
        ! Test whether tides have stripped the halo down to the solitonic core. The tidal radius used here is the *unlimited*
        ! one: `satelliteTidalStrippingRadiusLimited` clamps the radius used for stripping to be no smaller than the soliton
        ! radius, so a limited radius could never fall below it and this test would never be satisfied.
        radiusSoliton=darkMatterProfile%floatRank0MetaPropertyGet            (self%radiusSolitonID)
        radiusTidal  =self             %satelliteTidalStrippingRadius_%radius(node                )
        if (radiusSoliton > 0.0d0 .and. radiusTidal <= radiusSoliton) then
            ! The NFW envelope has been stripped away - trigger an interrupt to transition to the soliton-only state.
            interrupt         =  .true.
            functionInterrupt => solitonPhaseTransition
            self_             => self
            return
        else
            ! The halo is in the soliton+NFW phase. Strip only the outer NFW component.
            massLossRateOuter=+self%satelliteTidalStrippingOuter_%massLossRate(node)
            call satellite%boundMassRate(massLossRateOuter)
        end if
    else if (solitonStatus == solitonStatusSolitonOnly) then
        ! The halo is soliton-only. Evolve the bound mass and core mass together.
        massLossRateCore =+self%satelliteTidalStrippingCore_%massLossRate(node)
        call satellite        %boundMassRate             (massLossRateCore    )
        call darkMatterProfile%floatRank0MetaPropertyRate(                       &
                &                                         self%massCoreNormalID, &
                &                                         +0.25d0                &
                &                                         *massLossRateCore      &
                &                                        )
        ! Check if the evolved core mass has become negative.
        massCore = darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreID)
        if (massCore <= 0.0d0) then
            ! Destruction criterion met - trigger an interrupt.
            interrupt         =  .true.
            functionInterrupt => solitonPhaseTransition
            self_             => self
            return
        end if
    else
        ! No solitonic solution exists. Treat the halo as a pure NFW halo.
        massLossRateOuter=+self%satelliteTidalStrippingOuter_%massLossRate(node)
        call satellite%boundMassRate(massLossRateOuter)
    end if
    return
  end subroutine satelliteConditionalStrippingDifferentialEvolution

  subroutine solitonPhaseTransition(node,timeEnd)
    !!{RST
    Advance a satellite to the next stage of its solitonic evolution: transition a soliton+NFW halo whose NFW envelope has been
    tidally stripped away to the soliton-only state, and destroy a soliton-only halo whose core has been fully dissolved.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentSatellite
    implicit none
    type            (treeNode                      ), intent(inout), target   :: node
    double precision                                , intent(in   ), optional :: timeEnd
    class           (nodeComponentSatellite        ), pointer                 :: satellite
    class           (nodeComponentDarkMatterProfile), pointer                 :: darkMatterProfile
    double precision                                                          :: massCore         , radiusTidal, &
         &                                                                       radiusSoliton
    !$GLC attributes unused :: timeEnd

    satellite         => node             %satellite                (                     )
    darkMatterProfile => node             %darkMatterProfile        (                     )
    massCore          =  darkMatterProfile%floatRank0MetaPropertyGet(self_%massCoreID     )
    radiusSoliton     =  darkMatterProfile%floatRank0MetaPropertyGet(self_%radiusSolitonID)
    radiusTidal       =  self_%satelliteTidalStrippingRadius_%radius(      node           )
    ! Transition to the soliton-only state once the tidal radius has reached the soliton radius, so that no NFW envelope remains.
    if     (                                &
         &   radiusSoliton >  0.0d0         &
         &  .and.                           &
         &   radiusTidal   <= radiusSoliton &
         & ) call darkMatterProfile%integerRank0MetaPropertySet(self_%solitonStatusID,solitonStatusSolitonOnly%ID)
    ! Destroy the satellite once the solitonic core has fully dissolved.
    if (massCore <= 0.0d0) call satellite%destructionTimeSet(0.0d0)
    return
  end subroutine solitonPhaseTransition

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

  !+    Contributions to this file made by: Yu Zhao

  !!{RST
  Implements a node operator class that applies conditional mass loss to orbiting satellite halos.
  !!}

  use :: Dark_Matter_Profiles_Soliton_Status , only : enumerationSolitonStatusType      , solitonStatusSolitonNFW, solitonStatusSolitonOnly
  use :: Mass_Distribution_Soliton_Schive2014, only : massTotalPerMassCore
  use :: Satellite_Tidal_Stripping           , only : satelliteTidalStrippingClass
  use :: Satellite_Tidal_Stripping_Radii     , only : satelliteTidalStrippingRadiusClass

  !![
  <nodeOperator name="nodeOperatorSatelliteConditionalMassLoss" docformat="rst">
   <description>
   A node operator class that applies conditional tidal mass loss to orbiting satellite halos at each ODE timestep for use with
   the :galacticus-class:`darkMatterProfileDMOSolitonNFWHeated` dark matter profile. Depending on the halo state, it evolves
   either the outer halo or the solitonic core using a :galacticus-class:`satelliteTidalStrippingClass`. If no solitonic solution
   exists, the halo is treated as an NFW halo and tidal stripping is applied using the NFW model throughout its evolution. The
   satellite bound mass and, where applicable, the solitonic core mass are evolved as ODE variables.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteConditionalMassLoss
     !!{RST
     A node operator class that applies tidal mass loss to orbiting satellite halos.
     !!}
     private
     integer                                              :: massCoreID                              , solitonStatusID                       , &
          &                                                  massCoreNormalID                        , radiusSolitonID                       , &
          &                                                  massCoreTransitionID                    , randomOffsetID
     double precision                                     :: fractionMassCoreDestruction
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
  class(nodeOperatorSatelliteConditionalMassLoss), pointer :: self_ => null()
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
    double precision                                                 :: fractionMassCoreDestruction

    ! The outer-halo stripping model is built from the unqualified `satelliteTidalStripping` parameter, so that it is inherited
    ! from the enclosing scope and therefore always matches the tidal stripping model used elsewhere in the model. Only the core
    ! stripping model, which is specific to the solitonic core, is named explicitly.
    !
    ! The tidal radius model is used only to detect the transition to the soliton-only state. Either the limited or the
    ! unlimited model may be given: since the limited model returns max(radiusTidal,radiusSoliton), the test
    ! `radius <= radiusSoliton` made below is satisfied under exactly the same conditions for both. The unlimited model is the
    ! more natural choice, as the transition is a statement about where tides would strip to in the absence of the limit.
    !![
    <objectBuilder class="satelliteTidalStripping"       name="satelliteTidalStrippingOuter_"                                              source="parameters"/>
    <objectBuilder class="satelliteTidalStripping"       name="satelliteTidalStrippingCore_"   parameterName="satelliteTidalStrippingCore" source="parameters"/>
    <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_"                                             source="parameters"/>
    <inputParameter docformat="rst">
      <name>fractionMassCoreDestruction</name>
      <defaultValue>0.01d0</defaultValue>
      <source>parameters</source>
      <description>
      The fraction of the solitonic core mass, measured at the point at which the halo entered the soliton-only state, below which
      the satellite is destroyed. The core mass decays exponentially under the :cite:t:`du_tidal_2018` model, and so approaches
      zero without ever reaching it, so a threshold of this kind is required if a satellite is ever to be destroyed.
      </description>
    </inputParameter>
    !!]
    self=nodeOperatorSatelliteConditionalMassLoss(satelliteTidalStrippingOuter_,satelliteTidalStrippingCore_,satelliteTidalStrippingRadius_,fractionMassCoreDestruction)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingOuter_" />
    <objectDestructor name="satelliteTidalStrippingCore_"  />
    <objectDestructor name="satelliteTidalStrippingRadius_"/>
    !!]
    return
  end function satelliteConditionalStrippingConstructorParameters

  function satelliteConditionalStrippingConstructorInternal(satelliteTidalStrippingOuter_,satelliteTidalStrippingCore_,satelliteTidalStrippingRadius_,fractionMassCoreDestruction) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorSatelliteConditionalMassLoss` node operator class.
    !!}
    implicit none
    type            (nodeOperatorSatelliteConditionalMassLoss)                        :: self
    class           (satelliteTidalStrippingClass            ), intent(in   ), target :: satelliteTidalStrippingOuter_ , satelliteTidalStrippingCore_
    class           (satelliteTidalStrippingRadiusClass      ), intent(in   ), target :: satelliteTidalStrippingRadius_
    double precision                                          , intent(in   )         :: fractionMassCoreDestruction
    !![
    <constructorAssign variables="*satelliteTidalStrippingOuter_, *satelliteTidalStrippingCore_, *satelliteTidalStrippingRadius_, fractionMassCoreDestruction"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreTransition" id="self%massCoreTransitionID" isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRandomOffset"       id="self%randomOffsetID"       isEvolvable="no"  isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton"      id="self%radiusSolitonID"      isEvolvable="no"  isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal"     id="self%massCoreNormalID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"           id="self%massCoreID"           isEvolvable="no"  isCreator="no" />
    <addMetaProperty component="darkMatterProfile" name="solitonStatus"             id="self%solitonStatusID"      type="integer"    isCreator="no" />
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
    double precision                                                            :: massLossRateOuter, massLossRateCore  , &
              &                                                                    massCore         , radiusTidal       , &
              &                                                                    radiusSoliton    , massCoreTransition, &
              &                                                                    randomOffset
    type     (enumerationSolitonStatusType            )                         :: solitonStatus
    !$GLC attributes unused :: propertyType

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
        ! Test whether tides have stripped the halo down to the solitonic core. Note that this test gives the same result
        ! whether the tidal radius model used here is limited or unlimited, since `max(radiusTidal,radiusSoliton)` is less than
        ! or equal to `radiusSoliton` under precisely the same conditions as `radiusTidal` itself.
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
        call satellite%boundMassRate(massLossRateCore)
        ! Convert the rate of loss of bound mass into a rate of change of the *unscattered* core mass, which is the evolved
        ! property. Dividing by massTotalPerMassCore converts from the bound mass to the core mass, since the bound mass in this
        ! state is taken to be the total mass of the soliton. A further factor of
        ! 10^(-randomOffset) converts from the core mass to the unscattered core mass, as the dark matter profile applies the
        ! log-normal scatter of the core-halo mass relation via massCore = massCoreNormal * 10^randomOffset. Without it the
        ! *fractional* rate of decay of the core mass would be wrong by a factor of 10^randomOffset, and so would differ from
        ! halo to halo for no physical reason.
        randomOffset=+darkMatterProfile%floatRank0MetaPropertyGet(self%randomOffsetID)
        call darkMatterProfile%floatRank0MetaPropertyRate(                       &
                &                                         self%massCoreNormalID, &
                &                                         +massLossRateCore      &
                &                                         /massTotalPerMassCore  &
                &                                         /10.0d0**randomOffset  &
                &                                        )
        ! Test whether the solitonic core has been stripped down to the destruction threshold. Under the Du et al. (2018) model
        ! the core mass decays exponentially, approaching zero without ever reaching it, so destruction must be triggered at some
        ! finite fraction of the core mass at the point at which stripping of the core began.
        massCoreTransition=darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreTransitionID)
        if     (                                                                          &
             &   massCoreTransition > 0.0d0                                               &
             &  .and.                                                                     &
             &   massCore           < self%fractionMassCoreDestruction*massCoreTransition &
             & ) then
            ! Destruction criterion met - trigger an interrupt.
            interrupt         =  .true.
            functionInterrupt => solitonDestruction
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
    Transition a soliton+NFW halo whose NFW envelope has been tidally stripped away to the soliton-only state, recording the
    solitonic core mass at that point so that a destruction threshold can later be measured relative to it.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    type            (treeNode                      ), intent(inout), target   :: node
    double precision                                , intent(in   ), optional :: timeEnd
    class           (nodeComponentDarkMatterProfile), pointer                 :: darkMatterProfile
    double precision                                                          :: radiusTidal      , radiusSoliton
    !$GLC attributes unused :: timeEnd

    darkMatterProfile => node             %darkMatterProfile                    (                     )
    radiusSoliton     =  darkMatterProfile%floatRank0MetaPropertyGet            (self_%radiusSolitonID)
    radiusTidal       =  self_            %satelliteTidalStrippingRadius_%radius(      node           )
    ! Transition to the soliton-only state once the tidal radius has reached the soliton radius, so that no NFW envelope remains.
    if     (                                &
         &   radiusSoliton >  0.0d0         &
         &  .and.                           &
         &   radiusTidal   <= radiusSoliton &
         & ) then
       call darkMatterProfile%integerRank0MetaPropertySet(self_%solitonStatusID     ,solitonStatusSolitonOnly%ID                                  )
       call darkMatterProfile%floatRank0MetaPropertySet  (self_%massCoreTransitionID,darkMatterProfile%floatRank0MetaPropertyGet(self_%massCoreID))
    end if
    return
  end subroutine solitonPhaseTransition

  subroutine solitonDestruction(node,timeEnd)
    !!{RST
    Destroy a satellite whose solitonic core has been stripped below the destruction threshold.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    type            (treeNode              ), intent(inout), target   :: node
    double precision                        , intent(in   ), optional :: timeEnd
    class           (nodeComponentSatellite), pointer                 :: satellite
    !$GLC attributes unused :: timeEnd

    satellite => node%satellite()
    call satellite%destructionTimeSet(0.0d0)
    return
  end subroutine solitonDestruction

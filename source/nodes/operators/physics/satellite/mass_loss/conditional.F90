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

  use :: Satellite_Tidal_Stripping, only : satelliteTidalStrippingClass

  !![
  <nodeOperator name="nodeOperatorSatelliteConditionalMassLoss" docformat="rst">
   <description>
   A node operator class that strips dark matter from orbiting satellite halos at each ODE timestep, using a :galacticus-class:`satelliteTidalStrippingClass` to compute the tidal mass loss rate. The bound mass of the satellite is evolved as an ODE variable.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteConditionalMassLoss
     !!{RST
     A node operator class that applies tidal mass loss to orbiting satellite halos.
     !!}
     private
     integer :: massCoreID, solitonStatusID, massCoreNormalID
     class  (satelliteTidalStrippingClass), pointer :: satelliteTidalStrippingOuter_ => null(), satelliteTidalStrippingCore_ => null(), satelliteTidalStrippingNFW_ => null()
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
    class  (satelliteTidalStrippingClass            ), pointer       :: satelliteTidalStrippingOuter_, satelliteTidalStrippingCore_
    
    !![
    <objectBuilder   class="satelliteTidalStripping" name="satelliteTidalStrippingOuter_" source="parameters"        parameterName="satelliteTidalStrippingOuter"  />
    <objectBuilder   class="satelliteTidalStripping" name="satelliteTidalStrippingCore_"  source="parameters"        parameterName="satelliteTidalStrippingCore"  />
    !!]
    self=nodeOperatorSatelliteConditionalMassLoss(satelliteTidalStrippingOuter_,satelliteTidalStrippingCore_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingOuter_"/>
    <objectDestructor name="satelliteTidalStrippingCore_"/>
    !!]
    return
  end function satelliteConditionalStrippingConstructorParameters

  function satelliteConditionalStrippingConstructorInternal(satelliteTidalStrippingOuter_,satelliteTidalStrippingCore_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorSatelliteConditionalMassLoss` node operator class.
    !!}
    implicit none
    type   (nodeOperatorSatelliteConditionalMassLoss)                        :: self
    class  (satelliteTidalStrippingClass            ), intent(in   ), target :: satelliteTidalStrippingOuter_, satelliteTidalStrippingCore_
    
    !![ 
    <constructorAssign variables="*satelliteTidalStrippingOuter_, *satelliteTidalStrippingCore_"/>
    <addMetaProperty   component="darkMatterProfile"              name="solitonMassCore"           id="self%massCoreID"      isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonStatus"         id="self%solitonStatusID"         type="integer"    isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal" id="self%massCoreNormalID" isEvolvable="yes" isCreator="no"/>
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
    <objectDestructor name="self%satelliteTidalStrippingOuter_"/>
    <objectDestructor name="self%satelliteTidalStrippingCore_"/>
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
              &                                                                    solitonStatus
    double precision                                                            :: massSatellite, massCore
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (.not.node%isSatellite()) return
    ! Get required quantities from the satellite node.
    satellite         => node             %satellite        ()
    massSatellite     =  satellite        %boundMass        ()
    darkMatterProfile => node             %darkMatterProfile()
    massCore          =  darkMatterProfile%floatRank0MetaPropertyGet  (self%massCoreID)
    solitonStatus     =  darkMatterProfile%integerRank0MetaPropertyGet(self%solitonStatusID)
    
    if (solitonStatus == 1) then
        if (massSatellite > 0.0d0 .and. massSatellite < 4.0d0*massCore) then
            ! Destruction criterion met - trigger an interrupt.
            interrupt         =  .true.
            functionInterrupt => destructionTrigger
            self_             => self
            return
        else
            massLossRateOuter=+self%satelliteTidalStrippingOuter_%massLossRate(node)
            call satellite%boundMassRate(massLossRateOuter)
        end if
    else if (solitonStatus == 2) then
        massLossRateCore =+self%satelliteTidalStrippingCore_%massLossRate(node)
        call satellite        %boundMassRate             (massLossRateCore    )
        call darkMatterProfile%floatRank0MetaPropertyRate(                       &
                &                                         self%massCoreNormalID, &
                &                                         +0.25d0                &
                &                                         *massLossRateCore      &
                &                                        )
    else
        massLossRateOuter=+self%satelliteTidalStrippingOuter_%massLossRate(node)
        call satellite%boundMassRate(massLossRateOuter)
    end if
    return
  end subroutine satelliteConditionalStrippingDifferentialEvolution

  subroutine destructionTrigger(node, timeEnd)
    !!{RST
    Trigger destruction of the satellite by setting the time until destruction to zero.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentSatellite
    implicit none
    type            (treeNode                      ), intent(inout), target   :: node
    double precision                                , intent(in   ), optional :: timeEnd
    class           (nodeComponentSatellite        ), pointer                 :: satellite
    class           (nodeComponentDarkMatterProfile), pointer                 :: darkMatterProfile
    double precision                                                          :: massSatellite, massCore
    !$GLC attributes unused :: timeEnd
    
    satellite         => node     %satellite()
    massSatellite     =  satellite%boundMass()
    darkMatterProfile => node     %darkMatterProfile()
    massCore          = darkMatterProfile%floatRank0MetaPropertyGet(self_%massCoreID)
    if     (                                  &
         &   massSatellite > 0.0d0            &
         &  .and.                             &
         &   massSatellite < 4.0d0*massCore   &
         & ) then
        call darkMatterProfile%integerRank0MetaPropertySet(self_%solitonStatusID,2)
    end if
    return
  end subroutine destructionTrigger
  

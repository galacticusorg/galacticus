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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module which provides procedures shared by the galactic (disk, spheroid, and nuclear star cluster) node components.
!!}

module Node_Components_Galactic_Shared
  !!{RST
  Provides procedures shared by the galactic (disk, spheroid, and nuclear star cluster) node components. Since these components
  share no parent type below :galacticus-type:`nodeComponent`, each procedure is written once as a template and instantiated for
  each component which uses it.
  !!}
  implicit none
  private
  public :: Node_Component_Disk_Standard_Post_Evolve       , Node_Component_Spheroid_Standard_Post_Evolve, &
       &    Node_Component_NSC_Standard_Post_Evolve        , Node_Component_Disk_Very_Simple_Post_Evolve , &
       &    Node_Component_Spheroid_Very_Simple_Post_Evolve, Node_Component_Disk_Very_Simple_Post_Step   , &
       &    Node_Component_Spheroid_Very_Simple_Post_Step

  !![
  <generic identifier="historytrim">
   <instance label="Disk_Standard"        intrinsic="nodeComponentDisk"     implementation="nodeComponentDiskStandard"        accessor="disk"     description="disk"                />
   <instance label="Spheroid_Standard"    intrinsic="nodeComponentSpheroid" implementation="nodeComponentSpheroidStandard"    accessor="spheroid" description="spheroid"            />
   <instance label="NSC_Standard"         intrinsic="nodeComponentNSC"      implementation="nodeComponentNSCStandard"         accessor="NSC"      description="nuclear star cluster"/>
   <instance label="Disk_Very_Simple"     intrinsic="nodeComponentDisk"     implementation="nodeComponentDiskVerySimple"      accessor="disk"     description="disk"                />
   <instance label="Spheroid_Very_Simple" intrinsic="nodeComponentSpheroid" implementation="nodeComponentSpheroidVerySimple"  accessor="spheroid" description="spheroid"            />
  </generic>
  <generic identifier="negativegas">
   <instance label="Disk"     intrinsic="nodeComponentDisk"     implementation="nodeComponentDiskVerySimple"     accessor="disk"     default="defaultDiskComponent"     name="disk"     nameCapitalized="Disk"    />
   <instance label="Spheroid" intrinsic="nodeComponentSpheroid" implementation="nodeComponentSpheroidVerySimple" accessor="spheroid" default="defaultSpheroidComponent" name="spheroid" nameCapitalized="Spheroid"/>
  </generic>
  !!]

contains

  subroutine Node_Component_{historytrim¦label}_Post_Evolve(self,node)
    !!{RST
    Trim the future stellar properties history of the {historytrim¦description} component after evolution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, {historytrim¦intrinsic}, {historytrim¦implementation}, treeNode
    use :: Histories       , only : history
    implicit none
    class(*                       ), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node
    class({historytrim¦intrinsic} ), pointer       :: component
    class(nodeComponentBasic      ), pointer       :: basic
    type (history                 )                :: stellarPropertiesHistory
    !$GLC attributes unused :: self

    ! Get the {historytrim¦description} component.
    component => node%{historytrim¦accessor}()
    ! Check if a component of the expected implementation exists.
    select type (component)
    class is ({historytrim¦implementation})
       ! Trim the stellar populations properties future history.
       basic => node%basic()
       stellarPropertiesHistory=component%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(basic%time())
       call component%stellarPropertiesHistorySet(stellarPropertiesHistory)
    end select
    return
  end subroutine Node_Component_{historytrim¦label}_Post_Evolve

  subroutine Node_Component_{negativegas¦label}_Very_Simple_Post_Step(node,status)
    !!{RST
    Catch rounding errors in the very simple {negativegas¦name} gas evolution.
    !!}
    use :: Abundances_Structure          , only : abs                  , zeroAbundances
    use :: Display                       , only : displayMessage       , verbosityLevelWarn
    use :: Galacticus_Nodes              , only : {negativegas¦default}, {negativegas¦intrinsic}, {negativegas¦implementation}, treeNode
    use :: Interface_GSL                 , only : GSL_Success          , GSL_Continue
    use :: ISO_Varying_String            , only : assignment(=)        , operator(//)           , varying_string
    use :: Stellar_Luminosities_Structure, only : abs                  , zeroStellarLuminosities
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode                 ), intent(inout), pointer :: node
    integer                                    , intent(inout)          :: status
    class           ({negativegas¦intrinsic}  )               , pointer :: component
    double precision                           , save                   :: fractionalErrorMaximum=0.0d0
    double precision                                                    :: massComponent                , fractionalError
    character       (len=20                   )                         :: valueString
    type            (varying_string           ), save                   :: message
    !$omp threadprivate(message)

    if (.not.{negativegas¦default}%verySimpleIsActive()) return
    ! Get the {negativegas¦name} component.
    component => node%{negativegas¦accessor}()
    ! Check if a very simple {negativegas¦name} component exists.
    select type (component)
    class is ({negativegas¦implementation})
       ! Note that "status" is not set to failure as these changes in state of the {negativegas¦name} should not change any calculation of
       ! differential evolution rates as a negative gas mass was unphysical anyway.
       !
       ! Trap negative gas masses.
       if (component%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(component%massGas    ()) &
               &          /(                         &
               &                 component%massStellar()  &
               &            +abs(component%massGas    ()) &
               &           )
          !$omp critical (Very_Simple_{negativegas¦nameCapitalized}_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: {negativegas¦name} has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//node%index()//char(10)
             write (valueString,'(e12.6)') component%massGas    ()
             message=message//'  {negativegas¦nameCapitalized} gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') component%massStellar()
             message=message//'  {negativegas¦nameCapitalized} stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure     = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Gas mass will be reset to zero (in future cases also).'                                 //char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerical inaccuracy in the ODE solutions.'                  //char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Very_Simple_{negativegas¦nameCapitalized}_Post_Evolve_Check)
          ! Get the total mass of the {negativegas¦name} material
          massComponent= component%massGas    () &
               &        +component%massStellar()
          if (massComponent == 0.0d0) then
             call component%        massStellarSet(                  0.0d0)
             call component%  abundancesStellarSet(         zeroAbundances)
             call component%luminositiesStellarSet(zeroStellarLuminosities)
          end if
          ! Reset the gas mass of the {negativegas¦name}.
          call component%      massGasSet(         0.0d0)
          call component%abundancesGasSet(zeroAbundances)
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
    end select
    return
  end subroutine Node_Component_{negativegas¦label}_Very_Simple_Post_Step

end module Node_Components_Galactic_Shared

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
Contains a program which tests the stellar mass- and redshift-dependent dust attenuation class.
!!}

program Test_Dust_Attenuation_Stellar_Mass_Redshift
  !!{RST
  Tests the :galacticus-class:`dustAttenuationStellarMassRedshift` dust attenuation class: its mean attenuation, the limits
  on stellar mass, its transmission at the reference wavelength and (through the extinction curve) at other wavelengths, and
  that a negative mean attenuation gives no attenuation.
  !!}
  use :: Cosmology_Functions        , only : cosmologyFunctionsClass
  use :: Display                    , only : displayVerbositySet               , verbosityLevelStandard
  use :: Dust_Attenuation_Descriptors, only : emissionDescriptor
  use :: Dust_Attenuations          , only : dustAttenuationClass              , dustAttenuationStellarMassRedshift
  use :: Dust_Extinction_Curves     , only : dustExtinctionCurveClass
  use :: Error                      , only : Error_Handler_Register
  use :: Events_Hooks               , only : eventsHooksInitialize
  use :: Functions_Global_Utilities , only : Functions_Global_Set
  use :: Galactic_Structure_Options , only : componentTypeAll                  , componentTypeDisk
  use :: Galacticus_Nodes           , only : mergerTree                        , nodeClassHierarchyInitialize, nodeComponentBasic, nodeComponentDisk, &
       &                                     treeNode
  use :: Input_Parameters           , only : inputParameters
  use :: Node_Components            , only : Node_Components_Initialize        , Node_Components_Thread_Initialize
  use :: Unit_Tests                 , only : Assert                            , Unit_Tests_Begin_Group      , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  ! The parameters of the attenuation, as set in the parameter file.
  double precision                                    , parameter                 :: delta0            =+0.10d0   , deltaMass       =+0.20d0, &
       &                                                                             deltaRedshift     =+0.30d0   , deltaMassRedshift=+0.40d0, &
       &                                                                             redshiftPivot     = 1.00d0   , wavelengthHalpha =6564.61d0
  double precision                                    , dimension(3), parameter   :: wavelengths       =[4862.68d0,5008.24d0,3728.49d0]
  type            (inputParameters                   )                            :: parameters
  class           (dustAttenuationClass              ), pointer                   :: dustAttenuation_
  class           (cosmologyFunctionsClass           ), pointer                   :: cosmologyFunctions_
  class           (dustExtinctionCurveClass          ), pointer                   :: dustExtinctionCurve_
  type            (mergerTree                        ), target                    :: tree
  type            (treeNode                          ), pointer                   :: node
  type            (emissionDescriptor                ), dimension(4)              :: descriptors
  double precision                                    , dimension(4)              :: transmission
  double precision                                                                :: time              , redshift         , &
       &                                                                             massLogarithmic   , redshiftTerm     , &
       &                                                                             attenuationExpected, attenuationClipped
  integer                                                                         :: i

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call Unit_Tests_Begin_Group("Dust attenuation: stellar mass and redshift dependent")
  parameters=inputParameters('testSuite/parameters/dustAttenuationStellarMassRedshift.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="dustAttenuation"     name="dustAttenuation_"     source="parameters"/>
  <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
  <objectBuilder class="dustExtinctionCurve" name="dustExtinctionCurve_" source="parameters"/>
  !!]
  ! Emission at H-alpha, and at H-beta, [OIII] and [OII].
  descriptors%componentType=componentTypeDisk
  descriptors(1)%wavelength=wavelengthHalpha
  do i=1,3
     descriptors(i+1)%wavelength=wavelengths(i)
  end do
  select type (dustAttenuation_)
  class is (dustAttenuationStellarMassRedshift)
     call Assert('applies to all components combined',dustAttenuation_%supportsComponent(componentTypeAll),.true.)
     ! At the pivot redshift, and a stellar mass of 10^10 Msun, the attenuation is that of Garn & Best (2010) plus delta0.
     time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftPivot))
     call buildNode(1.0d10,time)
     call Assert('mean attenuation at pivot'            ,dustAttenuation_%attenuationMean(node),0.91d0+delta0,relTol=1.0d-9)
     ! At another mass and redshift, every term contributes.
     time           =cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(2.0d0))
     redshift       =cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(time))
     massLogarithmic=log10(3.0d10/1.0d10)
     redshiftTerm   =log((1.0d0+redshift)/(1.0d0+redshiftPivot))
     attenuationExpected=+0.91d0+0.77d0*massLogarithmic+0.11d0*massLogarithmic**2-0.09d0*massLogarithmic**3 &
          &              +delta0+deltaMass*massLogarithmic+deltaRedshift*redshiftTerm                     &
          &              +deltaMassRedshift*massLogarithmic*redshiftTerm
     call buildNode(3.0d10,time)
     call Assert('mean attenuation at z=2, 3e10 Msun'   ,dustAttenuation_%attenuationMean(node),attenuationExpected,relTol=1.0d-9)
     ! Transmission at H-alpha is exactly that of the mean attenuation; at other wavelengths it scales with the extinction curve.
     transmission=dustAttenuation_%transmission(node,descriptors)
     call Assert('transmission at H-alpha'              ,transmission(1),10.0d0**(-0.4d0*attenuationExpected),relTol=1.0d-9)
     do i=1,3
        call Assert('transmission scales with extinction curve',                                                                           &
             &      transmission(i+1)                                                                                                     , &
             &      10.0d0**(-0.4d0*attenuationExpected*dustExtinctionCurve_%attenuationRelative(wavelengths(i))/dustExtinctionCurve_%attenuationRelative(wavelengthHalpha)), &
             &      relTol=1.0d-9                                                                                                            &
             &     )
     end do
     ! Stellar masses outside the range 10^8 to 10^11 Msun are treated as being at the limit of the range.
     call buildNode(1.0d8,time)
     attenuationClipped=dustAttenuation_%attenuationMean(node)
     call buildNode(1.0d6,time)
     call Assert('stellar mass limited below'           ,dustAttenuation_%attenuationMean(node),attenuationClipped,relTol=1.0d-12)
     call buildNode(1.0d11,time)
     attenuationClipped=dustAttenuation_%attenuationMean(node)
     call buildNode(1.0d13,time)
     call Assert('stellar mass limited above'           ,dustAttenuation_%attenuationMean(node),attenuationClipped,relTol=1.0d-12)
     ! At the lower mass limit and z=3 the mean attenuation (with these parameters) is negative, and no attenuation is applied.
     time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(3.0d0))
     call buildNode(1.0d8,time)
     call Assert('mean attenuation is negative here'    ,dustAttenuation_%attenuationMean(node) < 0.0d0,.true.)
     transmission=dustAttenuation_%transmission(node,descriptors)
     call Assert('negative mean attenuation is not applied',transmission,[1.0d0,1.0d0,1.0d0,1.0d0],relTol=1.0d-12)
  class default
     call Assert('attenuation is of the expected class',.false.,.true.)
  end select
  !![
  <objectDestructor name="dustAttenuation_"    />
  <objectDestructor name="cosmologyFunctions_" />
  <objectDestructor name="dustExtinctionCurve_"/>
  !!]
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

contains

  subroutine buildNode(massStellar,time)
    !!{RST
    Build a node with the given stellar mass (in a disk) at the given time. A new node is built each time, since a node
    memoizes the mass distributions built from its components.
    !!}
    implicit none
    double precision                    , intent(in   ) :: massStellar, time
    class           (nodeComponentBasic), pointer       :: basic
    class           (nodeComponentDisk ), pointer       :: disk

    node  => treeNode(hostTree=tree)
    basic => node%basic(autoCreate=.true.)
    disk  => node%disk (autoCreate=.true.)
    call basic%massSet       (1.0d12     )
    call basic%timeSet       (time       )
    call disk %massStellarSet(massStellar)
    call disk %radiusSet     (3.0d-3     )
    return
  end subroutine buildNode

end program Test_Dust_Attenuation_Stellar_Mass_Redshift

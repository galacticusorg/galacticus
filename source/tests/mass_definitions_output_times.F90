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
Contains a program which tests halo mass definitions and output time conversions against independently computed reference values.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Mass_Definitions_Output_Times
  !!{RST
  Tests virial density contrasts, conversions of halo masses between density contrast definitions, and output time conversions,
  against reference values computed independently using `colossus <https://bdiemer.bitbucket.io/colossus/>`_ and `astropy
  <https://www.astropy.org/>`_ by the ``massDefinitionsOutputTimes.py`` script in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository.

  Note that Galacticus' :galacticus-class:`virialDensityContrastClass` implementations return contrasts relative to the *mean*
  matter density, while the fitting function of :cite:t:`bryan_statistical_1998`, and colossus' ``deltaVir()``, are relative to
  the *critical* density; the two differ by a factor :math:`\Omega_\mathrm{M}(z)`. Note also that the expansion rate of
  :galacticus-class:`cosmologyFunctionsMatterLambda` contains no radiation term, so radiation is excluded from both reference
  codes; retaining it would shift :math:`\Omega_\mathrm{M}(z)` by :math:`1.3\times 10^{-3}` at :math:`z=3`.

  Tolerances are:

  * :math:`\Omega_\mathrm{M}(z)` and virial density contrasts: :math:`10^{-9}` relative, as these involve no physical constants
    and are computed from closed-form expressions on both sides;
  * the virial radius: :math:`10^{-4}` relative. This carries the critical density, and so the difference between the GSL
    (Galacticus) and CODATA 2018/IAU 2015 (reference) values of :math:`\mathrm{G}` and :math:`\mathrm{M}_\odot`, which together
    shift :math:`\rho_\mathrm{crit}` by :math:`4.5\times 10^{-4}` and therefore the radius by :math:`2\times 10^{-5}`;
  * masses converted between density contrast definitions: :math:`10^{-3}` relative. Note that these are *insensitive* to the
    choice of physical constants---the critical density cancels between the target density and the profile
    normalization, so repeating the reference calculation with the GSL constants changes the converted masses by only
    :math:`10^{-16}`. The tolerance is instead set by :galacticus-class:`massDistributionNFW`, whose
    ``radiusEnclosingDensity`` inverts a tabulation of the scale-free enclosed density built with 30 points per octave in
    radius, rather than solving for the radius directly. The resulting interpolation error is amplified by the logarithmic
    slope of the enclosed mass, and so grows for contrasts enclosed well inside the virial radius: the differences from the
    reference are :math:`1.2\times 10^{-4}`, :math:`3.1\times 10^{-4}` and :math:`4.8\times 10^{-4}` for the 200m, virial and
    500c definitions respectively;
  * cosmic times: :math:`10^{-4}` relative, since the definitions of the year differ by :math:`2\times 10^{-5}` and the two
    codes integrate the expansion history independently;
  * redshifts inferred from cosmic times: an absolute tolerance of :math:`10^{-4}`. A relative tolerance is inappropriate here
    because the inference is ill-conditioned at low redshift: :math:`|\mathrm{d}t/\mathrm{d}z|\approx 12.9` Gyr at
    :math:`z\approx 0.06`, so the :math:`2\times 10^{-5}` relative difference in the definition of the year becomes a
    :math:`3\times 10^{-4}` *relative* difference in redshift, while remaining a :math:`2\times 10^{-5}` absolute one.
  !!}
  use, intrinsic :: ISO_C_Binding                       , only : c_size_t
  use            :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use            :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use            :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use            :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
  use            :: Dark_Matter_Profiles_DMO            , only : darkMatterProfileDMONFW
  use            :: Display                             , only : displayVerbositySet                               , &
  &                                                              verbosityLevelStandard
  use            :: Events_Hooks                        , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities          , only : Functions_Global_Set
  use            :: Galacticus_Nodes                    , only : nodeClassHierarchyInitialize                      , &
  &                                                              nodeComponentBasic                                , &
  &                                                              nodeComponentDarkMatterProfile                    , &
  &                                                              treeNode
  use            :: ISO_Varying_String                  , only : assignment(=)                                     , &
  &                                                              varying_string
  use            :: Input_Parameters                    , only : inputParameters
  use            :: Node_Components                     , only : Node_Components_Initialize                        , &
  &                                                              Node_Components_Thread_Initialize                 , &
  &                                                              Node_Components_Thread_Uninitialize               , &
  &                                                              Node_Components_Uninitialize
  use            :: Output_Times                        , only : outputTimesList                                   , &
  &                                                              outputTimesUniformSpacingInRedshift
  use            :: Unit_Tests                          , only : Assert                                            , &
  &                                                              Unit_Tests_Begin_Group                            , &
  &                                                              Unit_Tests_End_Group                              , &
  &                                                              Unit_Tests_Finish
  use            :: Virial_Density_Contrast             , only : fixedDensityTypeCritical                          , &
  &                                                              virialDensityContrastBryanNorman1998              , &
  &                                                              virialDensityContrastFixed
  implicit none
  type            (cosmologyParametersSimple                         )               :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                    )               :: cosmologyFunctions_
  type            (virialDensityContrastFixed                        )               :: virialDensityContrast200Critical_
  type            (virialDensityContrastBryanNorman1998              )               :: virialDensityContrastBryanNorman_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition)               :: darkMatterHaloScale_
  type            (darkMatterProfileDMONFW                           )               :: darkMatterProfileDMO_
  type            (outputTimesList                                   )               :: outputTimesList_
  type            (outputTimesUniformSpacingInRedshift               )               :: outputTimesUniformSpacingInRedshift_
  type            (treeNode                                          ), pointer      :: node
  class           (nodeComponentBasic                                ), pointer      :: basic
  class           (nodeComponentDarkMatterProfile                    ), pointer      :: darkMatterProfile
  type            (varying_string                                    )               :: parameterFile
  type            (inputParameters                                   )               :: parameters
  ! Properties of the test halo: a mass of 10¹²M☉ defined at a contrast of 200 times the critical density, with a concentration
  ! of 5, at z=0.
  double precision                                                    , parameter    :: massHalo200Critical                 =1.0000000000d+12                                    , &
       &                                                                                concentration                       =5.0000000000d+00
  double precision                                                    , parameter    :: radius200CriticalReference          =2.0766574918d-01
  ! Redshifts at which density contrasts are tested, and the corresponding values of Omega_M(z) and of the Bryan & Norman (1998)
  ! contrast relative to the mean matter density.
  double precision                                                    , dimension(3) :: redshiftContrast                    =[0.0000000000d+00,1.0000000000d+00,3.0000000000d+00]
  double precision                                                    , dimension(3) :: omegaMatterReference                =[2.8150000000d-01,7.5812152840d-01,9.6164829590d-01]
  double precision                                                    , dimension(3) :: contrastBryanNormanReference        =[3.5027506739d+02,2.0516122563d+02,1.8140798143d+02]
  double precision                                                    , dimension(3) :: omegaMatter                                                                              , &
       &                                                                                contrastBryanNorman
  ! Contrasts, relative to the mean matter density, of the 200m, virial (Bryan & Norman) and 500c definitions at z=0, and the
  ! corresponding halo masses.
  double precision                                                    , dimension(3) :: contrastTarget                      =[2.0000000000d+02,3.5027506739d+02,1.7761989343d+03]
  double precision                                                    , dimension(3) :: massDefinitionReference             =[1.4234267649d+12,1.2318376262d+12,7.2211395549d+11]
  double precision                                                    , dimension(3) :: massDefinition
  ! Output times: redshifts and the corresponding cosmic times.
  double precision                                                    , dimension(5) :: redshiftOutput                      =[0.0000000000d+00,5.0000000000d-01,1.0000000000d+00,3.0000000000d+00,9.0000000000d+00]
  double precision                                                    , dimension(5) :: timeOutputReference                 =[1.3846077111d+01,8.7237051267d+00,5.9751879363d+00,2.2016460844d+00,5.6040023511d-01]
  double precision                                                    , dimension(3) :: timeInput                           =[1.0000000000d+00,5.0000000000d+00,1.3000000000d+01]
  double precision                                                    , dimension(3) :: redshiftInputReference              =[5.7930295998d+00,1.2737585883d+00,6.2656567102d-02]
  double precision                                                    , dimension(3) :: redshiftFromTime
  double precision                                                    , dimension(5) :: timeFromRedshift
  double precision                                                                   :: radiusVirial                                                                             , &
       &                                                                                massIdentity
  integer                                                                            :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Mass definitions and output times")
  ! Read in controlling parameters and initialize the node component hierarchy.
  parameterFile='testSuite/parameters/basicDarkMatterProfileComponents.xml'
  parameters=inputParameters(parameterFile)
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <referenceConstruct object="cosmologyParameters_"            >
   <constructor>
    cosmologyParametersSimple                         (OmegaMatter=0.2815d0,OmegaBaryon=0.0465d0,OmegaDarkEnergy=0.7185d0,temperatureCMB=2.78d0,HubbleConstant=69.3d0)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"             >
   <constructor>
    cosmologyFunctionsMatterLambda                    (cosmologyParameters_=cosmologyParameters_)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast200Critical_">
   <constructor>
    virialDensityContrastFixed                        (densityContrastValue=200.0d0,densityType=fixedDensityTypeCritical,turnAroundOverVirialRadius=2.0d0,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrastBryanNorman_">
   <constructor>
    virialDensityContrastBryanNorman1998              (allowUnsupportedCosmology=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"            >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition(cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=virialDensityContrast200Critical_)
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileDMO_"           >
   <constructor>
    darkMatterProfileDMONFW                           (velocityDispersionUseSeriesExpansion=.false.,darkMatterHaloScale_=darkMatterHaloScale_)
   </constructor>
  </referenceConstruct>
  !!]

  ! Virial density contrasts. Galacticus returns contrasts relative to the mean matter density, so the Bryan & Norman (1998)
  ! fitting function, which is relative to the critical density, is divided by Omega_M(z) in the reference values.
  call Unit_Tests_Begin_Group("Virial density contrasts")
  do i=1,size(redshiftContrast)
     omegaMatter        (i)=cosmologyFunctions_              %omegaMatterEpochal(expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshiftContrast(i)))
     contrastBryanNorman(i)=virialDensityContrastBryanNorman_%densityContrast   (mass           =massHalo200Critical                                                  , &
          &                                                                     expansionFactor=cosmologyFunctions_%expansionFactorFromRedshift(redshiftContrast(i)) )
  end do
  call Assert("Omega_M(z)"                                  ,omegaMatter        ,omegaMatterReference        ,relTol=1.0d-9)
  call Assert("Bryan & Norman (1998) contrast {mean}"       ,contrastBryanNorman,contrastBryanNormanReference,relTol=1.0d-9)
  call Unit_Tests_End_Group()

  ! Halo mass definitions. The halo is defined by its mass at a contrast of 200 times the critical density, and its scale radius
  ! is set so that its concentration is 5.
  call Unit_Tests_Begin_Group("Halo mass definitions")
  node              => treeNode                  (                 )
  basic             => node    %basic            (autoCreate=.true.)
  darkMatterProfile => node    %darkMatterProfile(autoCreate=.true.)
  call basic%massSet            (massHalo200Critical                 )
  call basic%timeSet            (cosmologyFunctions_%cosmicTime(1.0d0))
  call basic%timeLastIsolatedSet(cosmologyFunctions_%cosmicTime(1.0d0))
  radiusVirial=darkMatterHaloScale_%radiusVirial(node)
  call darkMatterProfile%scaleSet(radiusVirial/concentration)
  call Assert("R200c"                                       ,radiusVirial       ,radius200CriticalReference  ,relTol=1.0d-4)
  ! Requesting the halo's own density contrast must return its own mass.
  massIdentity=Dark_Matter_Profile_Mass_Definition(                                                                                                           &
       &                                           node                 =node                                                                               , &
       &                                           densityContrast      =virialDensityContrast200Critical_%densityContrast(massHalo200Critical,basic%time()), &
       &                                           cosmologyParameters_ =cosmologyParameters_                                                               , &
       &                                           cosmologyFunctions_  =cosmologyFunctions_                                                                , &
       &                                           virialDensityContrast_=virialDensityContrast200Critical_                                                 , &
       &                                           darkMatterProfileDMO_=darkMatterProfileDMO_                                                                &
       &                                          )
  call Assert("mass at the halo's own contrast"             ,massIdentity       ,massHalo200Critical         ,relTol=1.0d-6)
  do i=1,size(contrastTarget)
     massDefinition(i)=Dark_Matter_Profile_Mass_Definition(                                                             &
          &                                                node                  =node                                , &
          &                                                densityContrast       =contrastTarget                   (i), &
          &                                                cosmologyParameters_  =cosmologyParameters_                , &
          &                                                cosmologyFunctions_   =cosmologyFunctions_                 , &
          &                                                virialDensityContrast_=virialDensityContrast200Critical_   , &
          &                                                darkMatterProfileDMO_ =darkMatterProfileDMO_                 &
          &                                               )
  end do
  call Assert("masses at 200m, virial, and 500c definitions",massDefinition     ,massDefinitionReference     ,relTol=1.0d-3)
  call Unit_Tests_End_Group()

  ! Output time conversions.
  call Unit_Tests_Begin_Group("Output times")
  outputTimesUniformSpacingInRedshift_=outputTimesUniformSpacingInRedshift(redshiftMinimum=0.0d0,redshiftMaximum=9.0d0,countRedshifts=5_c_size_t,cosmologyFunctions_=cosmologyFunctions_)
  do i=1,size(redshiftOutput)
     timeFromRedshift(i)=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftOutput(i)))
  end do
  call Assert("cosmic time from redshift"                   ,timeFromRedshift   ,timeOutputReference         ,relTol=1.0d-4)
  outputTimesList_=outputTimesList(times=timeInput,cosmologyFunctions_=cosmologyFunctions_)
  do i=1,size(timeInput)
     redshiftFromTime(i)=outputTimesList_%redshift(int(i,kind=c_size_t))
  end do
  call Assert("redshift from cosmic time"                   ,redshiftFromTime   ,redshiftInputReference      ,absTol=1.0d-4)
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
end program Test_Mass_Definitions_Output_Times

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
  Contains a program to test virial orbit distributions against independently computed reference values.
  !!}

program Test_Satellite_Virial_Orbits
  !!{RST
  Test the :galacticus-class:`virialOrbitBenson2005` and :galacticus-class:`virialOrbitJiang2014` classes against reference values
  computed independently of Galacticus by the ``referenceVirialOrbits.py`` script in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository, which writes each distribution from its paper.

  Both classes report the mean magnitude of the tangential velocity and the root mean squared total velocity at virial crossing,
  in each case as a dimensionless moment of the distribution multiplied by the virial velocity of the host under that class' own
  density contrast definition. Dividing by that virial velocity - obtained here from the same
  ``Dark_Matter_Profile_Mass_Definition`` call the classes themselves make - recovers the moment, which no cosmology enters, so
  the reference values are exact.

  For :cite:t:`benson_orbital_2005` the two moments are of the joint distribution of radial and tangential velocity of their
  eqns. (2)-(4), integrated over the square :math:`[0,3]^2` on which the class samples it. The class stores them as constants, so
  this is a direct check of those constants.

  For :cite:t:`jiang_orbital_2014` the total velocity follows a Voigt profile (their eqns. 9-11) whose parameters are tabulated in
  three host mass by three mass ratio bins in their Table 2, and the radial velocity the distribution
  :math:`P(V_\mathrm{r}/V) = A[\exp(B V_\mathrm{r}/V)-1]` of their eqn. (12), from which the mean tangential velocity at fixed
  total velocity follows in closed form. All nine bins are tested, which exercises the Table 2 transcription, the truncation of
  the Voigt profile, and the bin selection. A second object built from the class defaults is required to agree, which guards those
  defaults.

  The samplers themselves are not tested here - the orbits they return are propagated to the virial radius under the default
  density contrast, so their moments are not those of the distributions above.
  !!}
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda                    , cosmologyFunctionsClass
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
  use :: Dark_Matter_Profiles_DMO            , only : darkMatterProfileDMONFW
  use :: Display                             , only : displayVerbositySet                               , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Functions_Global_Utilities          , only : Functions_Global_Set
  use :: Galacticus_Nodes                    , only : mergerTree                                        , nodeClassHierarchyFinalize       , nodeClassHierarchyInitialize      , nodeComponentBasic          , &
       &                                              nodeComponentDarkMatterProfile                    , treeNode
  use :: Error                               , only : Error_Handler_Register
  use :: IO_HDF5                             , only : ioHDF5AccessInitialize
  use :: Input_Parameters                    , only : inputParameters
  use :: Node_Components                     , only : Node_Components_Initialize                        , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Numerical_Constants_Astronomical    , only : gravitationalConstant_internal
  use :: Unit_Tests                          , only : Assert                                            , Unit_Tests_Begin_Group           , Unit_Tests_End_Group              , Unit_Tests_Finish
  use :: Virial_Density_Contrast             , only : virialDensityContrastBryanNorman1998              , virialDensityContrastClass
  use :: Virial_Orbits                       , only : virialOrbitBenson2005                             , virialOrbitJiang2014
  implicit none
  type            (inputParameters                                   )                 :: parameters                       , parametersDefault
  type            (cosmologyParametersSimple                         ), pointer        :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                    ), pointer        :: cosmologyFunctions_
  type            (virialDensityContrastBryanNorman1998              ), pointer        :: virialDensityContrast_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition), pointer        :: darkMatterHaloScale_
  type            (darkMatterProfileDMONFW                           ), pointer        :: darkMatterProfileDMO_
  type            (virialOrbitBenson2005                             ), pointer        :: benson2005_
  type            (virialOrbitJiang2014                              ), pointer        :: jiang2014_                       , jiang2014Default_
  class           (virialDensityContrastClass                        ), pointer        :: densityContrastDefinition_
  type            (mergerTree                                        )                 :: tree
  type            (treeNode                                          ), pointer        :: hostBenson                       , satelliteBenson
  type            (treeNode                                          ), pointer        :: host                             , satellite
  ! Masses [M☉] of the Benson (2005) host and satellite, and the epoch [Gyr] at which every halo here is placed.
  double precision                                                    , parameter      :: massHostBenson            =1.0d13, massSatelliteBenson       =1.0d11, &
       &                                                                                  time                      =1.38d1
  ! Host masses and satellite-to-host mass ratios spanning the nine bins of Table 2 of Jiang et al. (2015). The bin edges are at
  ! host masses of 10¹²·⁵ and 10¹³·⁵ M☉ and mass ratios of 0.005 and 0.05; each value here sits well inside its bin, so that the
  ! conversion of these masses to the 200ρ_c definition used by the class can not move them across an edge.
  integer                                                             , parameter      :: countBins                 =3
  double precision                                                    , dimension(3)   :: massesHost                =[1.0d12,1.0d13,1.0d15]
  double precision                                                    , dimension(3)   :: ratiosMass                =[1.0d-3,1.0d-2,2.0d-1]
  ! Reference values computed by referenceVirialOrbits.py, in units of the host virial velocity.
  double precision                                                    , parameter      :: velocityTangentialBensonReference  =7.4926511681036d-01, &
       &                                                                                  velocityTotalRootMeanSquaredBensonReference=1.2544764469751d+00
  double precision                                                    , dimension(3,3) :: velocityTangentialJiangReference   =reshape(         &
       &                                                                                   [8.1177938799294d-01,8.0230942507903d-01,7.9165718528094d-01, &
       &                                                                                    7.5236159934832d-01,7.4643627159184d-01,6.8880085529828d-01, &
       &                                                                                    6.3312798321474d-01,5.9181435793411d-01,6.1994195023418d-01],[3,3])
  double precision                                                    , dimension(3,3) :: velocityTotalRootMeanSquaredJiangReference=reshape(  &
       &                                                                                   [1.2421287065605d+00,1.2578080464951d+00,1.2748533838986d+00, &
       &                                                                                    1.2001977653422d+00,1.2190723425164d+00,1.2512894577518d+00, &
       &                                                                                    1.1172623270511d+00,1.1158581114582d+00,1.0964966285461d+00],[3,3])
  ! Tolerances.
  !! The Benson (2005) constants are stored to six decimal digits, and the reference values differ from them in the seventh.
  double precision                                                    , parameter      :: toleranceBenson           =1.0d-6
  !! The Jiang et al. (2015) moments are evaluated by an adaptive integrator at a relative tolerance of 10⁻⁶, and their
  !! normalization by the class' own Voigt cumulative distribution rather than by quadrature of its density.
  double precision                                                    , parameter      :: toleranceJiang            =1.0d-4
  !! Relations which follow algebraically from quantities the classes have already returned.
  double precision                                                    , parameter      :: toleranceAlgebraic        =1.0d-9
  double precision                                                    , dimension(3,3) :: velocityTangentialJiang          , velocityTotalRootMeanSquaredJiang, &
       &                                                                                  ratioDefault                     , ratioExplicit
  double precision                                                                     :: massHost                         , radiusHost                       , &
       &                                                                                  velocityHost
  integer                                                                              :: i                                , j

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call ioHDF5AccessInitialize()
  call Unit_Tests_Begin_Group("Virial orbits")
  parameters=inputParameters('testSuite/parameters/satellites/virialOrbits.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build the objects shared by the orbit classes. These must match those built from the parameter file by the dark matter profile
  ! component.
  allocate(cosmologyParameters_  )
  allocate(cosmologyFunctions_   )
  allocate(virialDensityContrast_)
  allocate(darkMatterHaloScale_  )
  allocate(darkMatterProfileDMO_ )
  allocate(benson2005_           )
  allocate(jiang2014_            )
  allocate(jiang2014Default_     )
  cosmologyParameters_  =cosmologyParametersSimple                         (OmegaMatter=0.3153d0,OmegaBaryon=0.0493d0,OmegaDarkEnergy=0.6847d0,temperatureCMB=2.72548d0,HubbleConstant=67.36d0)
  cosmologyFunctions_   =cosmologyFunctionsMatterLambda                    (cosmologyParameters_=cosmologyParameters_)
  virialDensityContrast_=virialDensityContrastBryanNorman1998              (allowUnsupportedCosmology=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  darkMatterHaloScale_  =darkMatterHaloScaleVirialDensityContrastDefinition(cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=virialDensityContrast_)
  darkMatterProfileDMO_ =darkMatterProfileDMONFW                           (velocityDispersionUseSeriesExpansion=.false.,darkMatterHaloScale_=darkMatterHaloScale_)
  benson2005_           =virialOrbitBenson2005                             (                                                    &
       &                                                                    darkMatterHaloScale_  =darkMatterHaloScale_        , &
       &                                                                    cosmologyFunctions_   =cosmologyFunctions_         , &
       &                                                                    cosmologyParameters_  =cosmologyParameters_        , &
       &                                                                    virialDensityContrast_=virialDensityContrast_      , &
       &                                                                    darkMatterProfileDMO_ =darkMatterProfileDMO_         &
       &                                                                   )
  ! The Jiang et al. (2015) class is built twice: once with the Table 2 parameters written out here, and once from an empty
  ! parameter set so that it takes the class defaults. The two are required to agree below, which guards those defaults.
  jiang2014_            =virialOrbitJiang2014                              (                                                                          &
       &                                                                    bRatioLow             =[+0.049d0,+0.548d0,+1.229d0], &
       &                                                                    bRatioIntermediate    =[+1.044d0,+1.535d0,+3.396d0], &
       &                                                                    bRatioHigh            =[+2.878d0,+3.946d0,+2.982d0], &
       &                                                                    gammaRatioLow         =[+0.109d0,+0.114d0,+0.110d0], &
       &                                                                    gammaRatioIntermediate=[+0.098d0,+0.087d0,+0.050d0], &
       &                                                                    gammaRatioHigh        =[+0.071d0,+0.030d0,-0.012d0], &
       &                                                                    sigmaRatioLow         =[+0.077d0,+0.094d0,+0.072d0], &
       &                                                                    sigmaRatioIntermediate=[+0.073d0,+0.083d0,+0.118d0], &
       &                                                                    sigmaRatioHigh        =[+0.091d0,+0.139d0,+0.187d0], &
       &                                                                    muRatioLow            =[+1.220d0,+1.231d0,+1.254d0], &
       &                                                                    muRatioIntermediate   =[+1.181d0,+1.201d0,+1.236d0], &
       &                                                                    muRatioHigh           =[+1.100d0,+1.100d0,+1.084d0], &
       &                                                                    darkMatterHaloScale_  =darkMatterHaloScale_        , &
       &                                                                    cosmologyParameters_  =cosmologyParameters_        , &
       &                                                                    cosmologyFunctions_   =cosmologyFunctions_         , &
       &                                                                    virialDensityContrast_=virialDensityContrast_      , &
       &                                                                    darkMatterProfileDMO_ =darkMatterProfileDMO_         &
       &                                                                   )
  parametersDefault     =inputParameters                                   (                                                    )
  jiang2014Default_     =virialOrbitJiang2014                              (parametersDefault                                   )
  ! Test the Benson (2005) distribution.
  call Unit_Tests_Begin_Group("Benson (2005)")
  call buildNode(hostBenson     ,massHostBenson     )
  call buildNode(satelliteBenson,massSatelliteBenson)
  densityContrastDefinition_ => benson2005_%densityContrastDefinition()
  massHost                   =  Dark_Matter_Profile_Mass_Definition(hostBenson,densityContrastDefinition_%densityContrast(massHostBenson,time),radiusHost,velocityHost,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_)
  call Assert('mean tangential velocity'        ,benson2005_%velocityTangentialMagnitudeMean(satelliteBenson,hostBenson)/velocityHost,velocityTangentialBensonReference          ,relTol=toleranceBenson   )
  call Assert('root mean squared total velocity',benson2005_%velocityTotalRootMeanSquared   (satelliteBenson,hostBenson)/velocityHost,velocityTotalRootMeanSquaredBensonReference,relTol=toleranceBenson   )
  ! The mean angular momentum and mean energy follow from those two moments, the virial radius under this class' own density
  ! contrast definition, and the reduced mass.
  call Assert('mean angular momentum'           ,benson2005_%angularMomentumMagnitudeMean   (satelliteBenson,hostBenson)                                                                                      &
       &                                        ,+benson2005_%velocityTangentialMagnitudeMean(satelliteBenson,hostBenson)                                                                                     &
       &                                         *radiusHost                                                                                                                                                 &
       &                                         /(1.0d0+massSatelliteBenson/massHostBenson)                                                                                                                 &
       &                                        ,relTol=toleranceAlgebraic)
  call Assert('mean energy'                     ,benson2005_%energyMean                     (satelliteBenson,hostBenson)                                                                                      &
       &                                        ,+0.5d0                                                                                                                                                      &
       &                                         *benson2005_%velocityTotalRootMeanSquared  (satelliteBenson,hostBenson)**2                                                                                   &
       &                                         /(1.0d0+massSatelliteBenson/massHostBenson)                                                                                                                 &
       &                                         -gravitationalConstant_internal*massHost/radiusHost                                                                                                          &
       &                                        ,relTol=toleranceAlgebraic)
  call Unit_Tests_End_Group()
  ! Test the Jiang et al. (2015) distribution, in each of its nine bins.
  call Unit_Tests_Begin_Group("Jiang et al. (2015)")
  densityContrastDefinition_ => jiang2014_%densityContrastDefinition()
  do i=1,countBins
     do j=1,countBins
        call buildNode(host     ,massesHost(i)              )
        call buildNode(satellite,massesHost(i)*ratiosMass(j))
        massHost                              =Dark_Matter_Profile_Mass_Definition(host,densityContrastDefinition_%densityContrast(massesHost(i),time),radiusHost,velocityHost,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterProfileDMO_)
        velocityTangentialJiang        (i,j)  =+jiang2014_       %velocityTangentialMagnitudeMean(satellite,host) &
             &                                 /velocityHost
        velocityTotalRootMeanSquaredJiang(i,j)=+jiang2014_       %velocityTotalRootMeanSquared   (satellite,host) &
             &                                 /velocityHost
        ! The ratio of the two moments is independent of the host virial velocity, so the object built from the class defaults can
        ! be compared with the one built from the Table 2 values written out above without evaluating its own virial velocity.
        ratioExplicit                    (i,j)=+jiang2014_       %velocityTangentialMagnitudeMean(satellite,host) &
             &                                 /jiang2014_       %velocityTotalRootMeanSquared   (satellite,host)
        ratioDefault                     (i,j)=+jiang2014Default_%velocityTangentialMagnitudeMean(satellite,host) &
             &                                 /jiang2014Default_%velocityTotalRootMeanSquared   (satellite,host)
        call host     %destroy()
        call satellite%destroy()
     end do
  end do
  call Assert('mean tangential velocity'        ,velocityTangentialJiang          ,velocityTangentialJiangReference          ,relTol=toleranceJiang    )
  call Assert('root mean squared total velocity',velocityTotalRootMeanSquaredJiang,velocityTotalRootMeanSquaredJiangReference,relTol=toleranceJiang    )
  call Assert('default parameters are those of Table 2',ratioDefault              ,ratioExplicit                             ,relTol=toleranceAlgebraic)
  call Unit_Tests_End_Group()
  ! Clean up.
  call hostBenson     %destroy()
  call satelliteBenson%destroy()
  deallocate(benson2005_      )
  deallocate(jiang2014_       )
  deallocate(jiang2014Default_)
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()

contains

  subroutine buildNode(node_,mass_)
    !!{RST
    Build a node of the given mass, with a dark matter profile scale length following a plausible mass-concentration relation -
    the precise value is immaterial here, entering only through the conversion of the mass to each class' own density contrast
    definition.
    !!}
    implicit none
    type            (treeNode                      ), intent(  out), pointer :: node_
    double precision                                , intent(in   )          :: mass_
    class           (nodeComponentBasic             )              , pointer :: basic_
    class           (nodeComponentDarkMatterProfile )              , pointer :: darkMatterProfile_

    node_             => treeNode(hostTree=tree)
    basic_            => node_   %basic            (autoCreate=.true.)
    darkMatterProfile_=> node_   %darkMatterProfile(autoCreate=.true.)
    call basic_            %massSet            (mass_                                     )
    call basic_            %timeSet            (time                                      )
    call basic_            %timeLastIsolatedSet(time                                      )
    call darkMatterProfile_%scaleSet           (2.0d-2*(mass_/1.0d12)**(1.0d0/3.0d0)       )
    return
  end subroutine buildNode

end program Test_Satellite_Virial_Orbits

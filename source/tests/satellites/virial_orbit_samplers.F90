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
  Contains a program to test that virial orbit samplers reproduce the distributions they are drawn from.
  !!}

program Test_Satellite_Virial_Orbit_Samplers
  !!{RST
  Test that the orbit samplers of the :galacticus-class:`virialOrbitBenson2005`, :galacticus-class:`virialOrbitJiang2014` and
  :galacticus-class:`virialOrbitLi2020` classes draw from the distributions whose moments those classes report.

  ``tests.satellites.virial_orbits`` checks those reported moments against the papers. This test closes the other half of the
  chain: that the orbits actually returned by ``orbit()`` have those moments. Each sampler is asked for a large number of orbits,
  and the sample mean of the tangential velocity and the sample root mean squared total velocity are required to agree with
  :math:`\langle V_\theta \rangle` and :math:`\sqrt{\langle V^2 \rangle}` as reported by the same object.

  Each class samples its distribution at the virial radius under *its own* density contrast definition, and then propagates the
  resulting orbit to the virial radius under the density contrast of the ``darkMatterHaloScale`` it was given - which in general
  changes the velocities, so that the moments of the returned orbits are not those of the distribution. Here each class is
  therefore built with a ``darkMatterHaloScale`` whose density contrast is precisely that class' own definition: spherical
  collapse in a universe of collisionless matter and a cosmological constant for :cite:t:`benson_orbital_2005`, a fixed contrast
  of 200 times the critical density for :cite:t:`jiang_orbital_2014`, and :cite:t:`bryan_statistical_1998` for
  :cite:t:`li_orbital_2020`. The two radii then coincide and the propagation is a no-op, so the returned velocities are those
  drawn. Orbits are requested with ``acceptUnboundOrbits`` set, since the reported moments are of the unrestricted distribution.

  The agreement required is statistical: each sample mean must lie within five standard errors of the reported moment, with the
  standard error estimated from the sample itself. So that this can not pass merely by being imprecise, each standard error is
  also required to be a small fraction of the moment it tests, which is what makes the comparison sharp.
  !!}
  use, intrinsic :: ISO_C_Binding                     , only : c_long
  use            :: Cosmology_Functions               , only : cosmologyFunctionsMatterLambda
  use            :: Cosmology_Parameters              , only : cosmologyParametersSimple
  use            :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use            :: Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMONFW
  use            :: Display                           , only : displayVerbositySet                              , verbosityLevelStandard           , displayMessage
  use            :: Error                             , only : Error_Handler_Register
  use            :: Events_Hooks                      , only : eventsHooksInitialize
  use            :: Functions_Global_Utilities        , only : Functions_Global_Set
  use            :: Galacticus_Nodes                  , only : mergerTree                                       , nodeClassHierarchyFinalize       , nodeClassHierarchyInitialize      , nodeComponentBasic          , &
       &                                                       nodeComponentDarkMatterProfile                   , treeNode
  use            :: IO_HDF5                           , only : ioHDF5AccessInitialize
  use            :: Input_Parameters                  , only : inputParameters
  use            :: Kepler_Orbits                     , only : keplerOrbit
  use            :: Node_Components                   , only : Node_Components_Initialize                       , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use            :: Numerical_Random_Numbers          , only : randomNumberGeneratorGSL
  use            :: Unit_Tests                        , only : Assert                                           , Unit_Tests_Begin_Group           , Unit_Tests_End_Group              , Unit_Tests_Finish
  use            :: Virial_Density_Contrast           , only : virialDensityContrastFixed                       , fixedDensityTypeCritical         , virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use            :: Virial_Orbits                     , only : virialOrbitBenson2005                            , virialOrbitJiang2014             , virialOrbitLi2020                 , virialOrbitClass
  implicit none
  type            (inputParameters                                             )          :: parameters
  type            (cosmologyParametersSimple                                   ), pointer :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              ), pointer :: cosmologyFunctions_
  type            (darkMatterProfileDMONFW                                     ), pointer :: darkMatterProfileDMO_
  ! One density contrast per class, each precisely that class' own definition, with a dark matter halo scale built upon it.
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: contrastBenson_
  type            (virialDensityContrastFixed                                  ), pointer :: contrastJiang_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition           ), pointer :: scaleBenson_                    , scaleJiang_
  type            (virialOrbitBenson2005                                       ), pointer :: benson2005_
  type            (virialOrbitJiang2014                                        ), pointer :: jiang2014_
  class           (virialOrbitClass                                            ), pointer :: virialOrbit_
  type            (mergerTree                                                  )          :: tree
  type            (treeNode                                                    ), pointer :: host                            , satellite
  ! Halo masses [M☉] and the epoch [Gyr] at which they are placed.
  double precision                                                              , parameter :: massHost                =1.0d13, massSatellite             =1.0d11, &
       &                                                                                       time                    =1.38d1
  ! The number of orbits drawn from each sampler. At this size the standard error on the mean tangential velocity is a few parts
  ! in ten thousand of the mean, so a five standard error band is a fractional agreement of order 10⁻³.
  integer                                                                       , parameter :: countOrbits             =2000000
  ! The width of the band, in standard errors, and the largest standard error - as a fraction of the moment it tests - for which
  ! the comparison is considered sharp enough to be meaningful.
  double precision                                                              , parameter :: widthBand               =5.0d0 , standardErrorMaximum      =2.0d-3
  double precision                                                                          :: velocityTangentialMean          , velocityTotalRootMeanSquared     , &
       &                                                                                       errorVelocityTangential         , errorVelocityTotalRootMeanSquared

  call displayVerbositySet(verbosityLevelStandard)
  call Error_Handler_Register()
  call ioHDF5AccessInitialize()
  call Unit_Tests_Begin_Group("Virial orbit samplers")
  parameters=inputParameters('testSuite/parameters/satellites/virialOrbitSamplers.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build the objects shared by all three classes.
  allocate(cosmologyParameters_ )
  allocate(cosmologyFunctions_  )
  allocate(contrastBenson_      )
  allocate(contrastJiang_       )
  allocate(scaleBenson_         )
  allocate(scaleJiang_          )
  allocate(darkMatterProfileDMO_)
  allocate(benson2005_          )
  allocate(jiang2014_           )
  cosmologyParameters_ =cosmologyParametersSimple     (OmegaMatter=0.3153d0,OmegaBaryon=0.0493d0,OmegaDarkEnergy=0.6847d0,temperatureCMB=2.72548d0,HubbleConstant=67.36d0)
  cosmologyFunctions_  =cosmologyFunctionsMatterLambda(cosmologyParameters_=cosmologyParameters_)
  ! The density contrast definitions of the three classes, reproduced here so that each class can be given a dark matter halo
  ! scale which uses its own definition.
  contrastBenson_      =virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(tableStore=.true.,cosmologyFunctions_=cosmologyFunctions_)
  contrastJiang_       =virialDensityContrastFixed                                    (densityContrastValue=200.0d0,densityType=fixedDensityTypeCritical,turnAroundOverVirialRadius=2.0d0,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  scaleBenson_         =darkMatterHaloScaleVirialDensityContrastDefinition            (cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=contrastBenson_)
  scaleJiang_          =darkMatterHaloScaleVirialDensityContrastDefinition            (cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,virialDensityContrast_=contrastJiang_ )
  darkMatterProfileDMO_=darkMatterProfileDMONFW                                       (velocityDispersionUseSeriesExpansion=.false.,darkMatterHaloScale_=scaleBenson_)
  benson2005_          =virialOrbitBenson2005                                         (                                             &
       &                                                                               darkMatterHaloScale_  =scaleBenson_        , &
       &                                                                               cosmologyFunctions_   =cosmologyFunctions_ , &
       &                                                                               cosmologyParameters_  =cosmologyParameters_, &
       &                                                                               virialDensityContrast_=contrastBenson_     , &
       &                                                                               darkMatterProfileDMO_ =darkMatterProfileDMO_ &
       &                                                                              )
  jiang2014_           =virialOrbitJiang2014                                          (                                                                          &
       &                                                                               bRatioLow             =[+0.049d0,+0.548d0,+1.229d0], &
       &                                                                               bRatioIntermediate    =[+1.044d0,+1.535d0,+3.396d0], &
       &                                                                               bRatioHigh            =[+2.878d0,+3.946d0,+2.982d0], &
       &                                                                               gammaRatioLow         =[+0.109d0,+0.114d0,+0.110d0], &
       &                                                                               gammaRatioIntermediate=[+0.098d0,+0.087d0,+0.050d0], &
       &                                                                               gammaRatioHigh        =[+0.071d0,+0.030d0,-0.012d0], &
       &                                                                               sigmaRatioLow         =[+0.077d0,+0.094d0,+0.072d0], &
       &                                                                               sigmaRatioIntermediate=[+0.073d0,+0.083d0,+0.118d0], &
       &                                                                               sigmaRatioHigh        =[+0.091d0,+0.139d0,+0.187d0], &
       &                                                                               muRatioLow            =[+1.220d0,+1.231d0,+1.254d0], &
       &                                                                               muRatioIntermediate   =[+1.181d0,+1.201d0,+1.236d0], &
       &                                                                               muRatioHigh           =[+1.100d0,+1.100d0,+1.084d0], &
       &                                                                               darkMatterHaloScale_  =scaleJiang_                 , &
       &                                                                               cosmologyParameters_  =cosmologyParameters_        , &
       &                                                                               cosmologyFunctions_   =cosmologyFunctions_         , &
       &                                                                               virialDensityContrast_=contrastJiang_              , &
       &                                                                               darkMatterProfileDMO_ =darkMatterProfileDMO_         &
       &                                                                              )
  ! Build the nodes, and a random number generator for the tree from which the samplers draw.
  call buildNode(host     ,massHost     )
  call buildNode(satellite,massSatellite)
  tree%nodeBase => host
  allocate(randomNumberGeneratorGSL :: tree%randomNumberGenerator_)
  select type (randomNumberGenerator_ => tree%randomNumberGenerator_)
  type is (randomNumberGeneratorGSL)
     randomNumberGenerator_=randomNumberGeneratorGSL(seed_=1379_c_long)
  end select
  call tree%properties%initialize()
  ! Test each sampler in turn.
  call Unit_Tests_Begin_Group("Benson (2005)")
  virialOrbit_ => benson2005_
  call sampleOrbits(velocityTangentialMean,errorVelocityTangential,velocityTotalRootMeanSquared,errorVelocityTotalRootMeanSquared)
  call assertSample(velocityTangentialMean,errorVelocityTangential,velocityTotalRootMeanSquared,errorVelocityTotalRootMeanSquared)
  call Unit_Tests_End_Group()
  call Unit_Tests_Begin_Group("Jiang et al. (2015)")
  virialOrbit_ => jiang2014_
  call sampleOrbits(velocityTangentialMean,errorVelocityTangential,velocityTotalRootMeanSquared,errorVelocityTotalRootMeanSquared)
  call assertSample(velocityTangentialMean,errorVelocityTangential,velocityTotalRootMeanSquared,errorVelocityTotalRootMeanSquared)
  call Unit_Tests_End_Group()
  ! The Li et al. (2020) class is built from the parameter file, so that the critical overdensity and cosmological mass variance
  ! objects it requires are built too, but with the dark matter halo scale built here so that its own density contrast definition
  ! is used.
  call Unit_Tests_Begin_Group("Li et al. (2020)")
  !![
  <objectBuilder class="virialOrbit" name="virialOrbit_" source="parameters"/>
  !!]
  select type (orbit_ => virialOrbit_)
  type is (virialOrbitLi2020)
     call sampleOrbits(velocityTangentialMean,errorVelocityTangential,velocityTotalRootMeanSquared,errorVelocityTotalRootMeanSquared)
     call assertSample(velocityTangentialMean,errorVelocityTangential,velocityTotalRootMeanSquared,errorVelocityTotalRootMeanSquared)
  class default
     call Assert('the virial orbit class built is li2020',.false.,.true.)
  end select
  call Unit_Tests_End_Group()
  ! Clean up.
  call host     %destroy()
  call satellite%destroy()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()

contains

  subroutine sampleOrbits(velocityTangentialMean_,errorVelocityTangential_,velocityTotalRootMeanSquared_,errorVelocityTotalRootMeanSquared_)
    !!{RST
    Draw orbits from the current sampler, returning the sample mean of the tangential velocity and the sample root mean squared
    total velocity, each with its standard error. The standard error of the root mean squared velocity follows from that of the
    mean squared velocity by propagation through the square root.
    !!}
    implicit none
    double precision            , intent(  out) :: velocityTangentialMean_      , errorVelocityTangential_          , &
         &                                         velocityTotalRootMeanSquared_, errorVelocityTotalRootMeanSquared_
    type            (keplerOrbit)                :: orbit
    double precision                             :: velocityTangential           , velocityTotalSquared             , &
         &                                          sumTangential                , sumTangentialSquared             , &
         &                                          sumTotalSquared              , sumTotalFourth                   , &
         &                                          velocityTotalMeanSquared     , varianceTangential               , &
         &                                          varianceTotalSquared
    integer                                      :: i

    sumTangential       =0.0d0
    sumTangentialSquared=0.0d0
    sumTotalSquared     =0.0d0
    sumTotalFourth      =0.0d0
    do i=1,countOrbits
       orbit               = virialOrbit_%orbit(satellite,host,acceptUnboundOrbits=.true.)
       velocityTangential  = orbit%velocityTangential()
       velocityTotalSquared=+orbit%velocityTangential()**2 &
            &               +orbit%velocityRadial    ()**2
       sumTangential       =sumTangential       +velocityTangential
       sumTangentialSquared=sumTangentialSquared+velocityTangential  **2
       sumTotalSquared     =sumTotalSquared     +velocityTotalSquared
       sumTotalFourth      =sumTotalFourth      +velocityTotalSquared**2
    end do
    velocityTangentialMean_      =+sumTangential  /dble(countOrbits)
    velocityTotalMeanSquared     =+sumTotalSquared/dble(countOrbits)
    velocityTotalRootMeanSquared_=sqrt(velocityTotalMeanSquared)
    varianceTangential           =+sumTangentialSquared/dble(countOrbits)-velocityTangentialMean_  **2
    varianceTotalSquared         =+sumTotalFourth      /dble(countOrbits)-velocityTotalMeanSquared **2
    errorVelocityTangential_     =sqrt(varianceTangential  /dble(countOrbits))
    ! d(√x)/dx = 1/(2√x), so the standard error of the root mean squared velocity is that of the mean squared velocity divided by
    ! twice the root mean squared velocity.
    errorVelocityTotalRootMeanSquared_=+sqrt(varianceTotalSquared/dble(countOrbits)) &
         &                             /2.0d0                                       &
         &                             /velocityTotalRootMeanSquared_
    return
  end subroutine sampleOrbits

  subroutine assertSample(velocityTangentialMean_,errorVelocityTangential_,velocityTotalRootMeanSquared_,errorVelocityTotalRootMeanSquared_)
    !!{RST
    Assert that the sample moments agree with those reported by the current class, to within the chosen band of standard errors,
    and that those standard errors are small enough for the comparison to be meaningful.
    !!}
    implicit none
    double precision, intent(in   ) :: velocityTangentialMean_      , errorVelocityTangential_          , &
         &                             velocityTotalRootMeanSquared_, errorVelocityTotalRootMeanSquared_
    double precision                :: velocityTangentialReported   , velocityTotalRootMeanSquaredReported
    character       (len=256      ) :: message

    velocityTangentialReported          =virialOrbit_%velocityTangentialMagnitudeMean(satellite,host)
    velocityTotalRootMeanSquaredReported=virialOrbit_%velocityTotalRootMeanSquared   (satellite,host)
    ! Report the comparison, so that the margin by which it passes is visible and a statistical test can not quietly become
    ! uninformative.
    write (message,'(a,f9.4,a,f9.4,a,f7.4,a,f6.2,a)') '  <V_θ>: sampled ',velocityTangentialMean_      ,' km/s, distribution ',velocityTangentialReported          ,' km/s, standard error ',errorVelocityTangential_          ,' km/s, difference ',abs(velocityTangentialMean_      -velocityTangentialReported          )/errorVelocityTangential_          ,'σ'
    call displayMessage(trim(message))
    write (message,'(a,f9.4,a,f9.4,a,f7.4,a,f6.2,a)') '  √<V²>: sampled ',velocityTotalRootMeanSquared_,' km/s, distribution ',velocityTotalRootMeanSquaredReported,' km/s, standard error ',errorVelocityTotalRootMeanSquared_,' km/s, difference ',abs(velocityTotalRootMeanSquared_-velocityTotalRootMeanSquaredReported)/errorVelocityTotalRootMeanSquared_,'σ'
    call displayMessage(trim(message))
    ! Require the comparison to be sharp before making it, so that it can not pass by being imprecise.
    call Assert('standard error of the mean tangential velocity is small'   ,errorVelocityTangential_          /velocityTangentialReported           < standardErrorMaximum,.true.)
    call Assert('standard error of the root mean squared velocity is small' ,errorVelocityTotalRootMeanSquared_/velocityTotalRootMeanSquaredReported < standardErrorMaximum,.true.)
    call Assert('sampled mean tangential velocity matches the distribution' ,abs(velocityTangentialMean_      -velocityTangentialReported          ) < widthBand*errorVelocityTangential_          ,.true.)
    call Assert('sampled root mean squared velocity matches the distribution',abs(velocityTotalRootMeanSquared_-velocityTotalRootMeanSquaredReported) < widthBand*errorVelocityTotalRootMeanSquared_,.true.)
    return
  end subroutine assertSample

  subroutine buildNode(node_,mass_)
    !!{RST
    Build a node of the given mass, with a dark matter profile scale length following a plausible mass-concentration relation.
    !!}
    implicit none
    type            (treeNode                      ), intent(  out), pointer :: node_
    double precision                                , intent(in   )          :: mass_
    class           (nodeComponentBasic            )               , pointer :: basic_
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile_

    node_             => treeNode(hostTree=tree)
    basic_            => node_   %basic            (autoCreate=.true.)
    darkMatterProfile_=> node_   %darkMatterProfile(autoCreate=.true.)
    call basic_            %massSet            (mass_                               )
    call basic_            %timeSet            (time                                )
    call basic_            %timeLastIsolatedSet(time                                )
    call darkMatterProfile_%scaleSet           (2.0d-2*(mass_/1.0d12)**(1.0d0/3.0d0))
    return
  end subroutine buildNode

end program Test_Satellite_Virial_Orbit_Samplers

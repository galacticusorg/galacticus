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
  Contains a program to test the molecular fraction star formation rate surface density laws against independently computed
  reference values.
  !!}

program Test_Star_Formation_Molecular_Fractions
  !!{RST
  Test the :galacticus-class:`starFormationRateSurfaceDensityDisksBlitz2006` and
  :galacticus-class:`starFormationRateSurfaceDensityDisksKrumholz2009` classes against reference values computed independently of
  Galacticus by the ``referenceMolecularFractions.py`` script in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository, which writes both laws directly from their papers. Both
  classes are built from their own default parameters, so that those defaults are tested along with the formulas.

  For :cite:t:`blitz_role_2006` the midplane pressure of a disk of locally isothermal gas and stars,

  .. math::

     P_\mathrm{ext} = {\pi \over 2} \G \Sigma_\mathrm{gas} \left[ \Sigma_\mathrm{gas} + {\sigma_\mathrm{gas} \over \sigma_\star}
     \Sigma_\star \right], \,\,\, \sigma_\star = \sqrt{\pi \G h_\star \Sigma_\star},

  is compared with the reference in physical units, which fixes the normalization :math:`P_0/k_\mathrm{B} = 3.5 \times 10^4`
  K cm\ :math:`^{-3}` of Table 2 (and eqn. 12) of that paper. The ratio of molecular to atomic hydrogen, :math:`R_\mathrm{mol} =
  (P_\mathrm{ext}/P_0)^\alpha` (their eqn. 11), is checked to scale as :math:`P_\mathrm{ext}^{0.92}`, and the molecular fraction
  recovered from the star formation rate is checked against :math:`f_\mathrm{H_2} = R_\mathrm{mol}/(1+R_\mathrm{mol})` (their
  eqn. 21) rather than the :math:`\min(R_\mathrm{mol},1)` which preceded it. The fully molecular limit and the full rate surface
  density are compared with the reference also.

  For :cite:t:`krumholz_star_2009` the rate surface density is compared with the reference at five radii, for both the full
  fitting function of their eqn. (2) and the faster :cite:p:`mckee_atomic--molecular_2010` form, and the limits of zero
  metallicity, of vanishing molecular fraction at low surface density, and of unit molecular fraction at high surface density are
  checked. Their eqn. (10) uses the *total* gas surface density, including helium, so the test also confirms that the star
  formation rate reduces to :math:`\nu_\mathrm{SF} f_\mathrm{H_2} \Sigma_\mathrm{gas}` precisely where the *total* gas surface
  density equals the transition surface density of 85 :math:`\mathrm{M}_\odot` pc\ :math:`^{-2}`.
  !!}
  use :: Abundances_Structure                     , only : abundances
  use :: Display                                  , only : displayVerbositySet                          , verbosityLevelStandard
  use :: Events_Hooks                             , only : eventsHooksInitialize
  use :: Functions_Global_Utilities               , only : Functions_Global_Set
  use :: Galacticus_Nodes                         , only : mergerTree                                   , nodeClassHierarchyFinalize                     , nodeClassHierarchyInitialize      , nodeComponentBasic          , &
       &                                                   nodeComponentDisk                            , treeNode
  use :: Input_Parameters                         , only : inputParameters
  use :: Node_Components                          , only : Node_Components_Initialize                    , Node_Components_Thread_Initialize              , Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Numerical_Constants_Astronomical         , only : hydrogenByMassSolar                          , metallicitySolar
  use :: Numerical_Constants_Math                 , only : Pi
  use :: Numerical_Constants_Prefixes             , only : giga                                         , mega
  use :: Star_Formation_Rate_Surface_Density_Disks, only : starFormationRateSurfaceDensityDisksBlitz2006, starFormationRateSurfaceDensityDisksKrumholz2009
  use :: Unit_Tests                               , only : Assert                                       , Unit_Tests_Begin_Group                         , Unit_Tests_End_Group              , Unit_Tests_Finish
  implicit none
  type            (inputParameters                                 )                 :: parameters                            , parametersDefault
  type            (starFormationRateSurfaceDensityDisksBlitz2006   ), pointer        :: blitz2006_
  type            (starFormationRateSurfaceDensityDisksKrumholz2009), pointer        :: krumholz2009_                         , krumholz2009Fast_
  type            (mergerTree                                      )                 :: tree
  type            (treeNode                                        ), pointer        :: node                                  , nodeCompressedTen                    , &
       &                                                                                nodeCompressedTwenty
  class           (nodeComponentBasic                              ), pointer        :: basic
  class           (nodeComponentDisk                               ), pointer        :: disk
  type            (abundances                                      )                 :: abundancesGas
  ! Disk properties. These, and the Solar metallicity of the gas, must match those assumed by the reference script.
  double precision                                                  , parameter      :: massGas                   =1.0d10     , massStellar               =5.0d10 , & ! [M☉]
       &                                                                                radiusScale               =3.0d-3                                            ! [Mpc]
  ! Parameters of the two star formation laws. Each is the default of the corresponding class, and is used here only to recover
  ! the molecular fraction from the star formation rate - the classes themselves are built from their defaults.
  double precision                                                  , parameter      :: surfaceDensityCritical    =2.0d+02    , surfaceDensityExponent    =0.4d0  , & ! [M☉/pc²]
       &                                                                                frequencyNormalization    =5.25d-10                                        , & ! [yr⁻¹]
       &                                                                                pressureCharacteristic    =3.5d+04    , pressureExponent          =0.92d0 , & ! [K/cm³]
       &                                                                                frequencyStarFormation    =0.385d0    , clumpingFactor            =5.0d0      ! [Gyr⁻¹]
  ! The transition surface density of eqn. (10) of Krumholz et al. (2009), in Galacticus' internal units.
  double precision                                                  , parameter      :: surfaceDensityTransition  =85.0d12                                           ! [M☉/Mpc²]
  ! Radii at which each law is tested, in units of the disk scale length. For the Krumholz et al. (2009) law the second is the
  ! radius at which the gas surface density equals the transition surface density.
  integer                                                           , parameter      :: countRadiiBlitz           =6          , countRadiiKrumholz        =5
  double precision                                                  , dimension(6)   :: radiiBlitz                =[0.5d0,1.0d0,2.0d0,3.0d0,4.0d0,6.0d0]
  double precision                                                  , dimension(5)   :: radiiKrumholz             =[0.5d0,0.7325874717403015d0,1.0d0,2.0d0,3.0d0]
  ! Reference values computed by referenceMolecularFractions.py.
  !! Gas surface densities [M☉/pc²] at the Blitz & Rosolowsky (2006) radii.
  double precision                                                  , dimension(6)   :: surfaceDensityGasReference=[1.0725816958895d+02,6.5055368360355d+01,2.3932532557610d+01,8.8042867031108d+00,3.2389160722535d+00,4.3833962401805d-01]
  !! Midplane pressures, P_ext/k_B [K/cm³].
  double precision                                                  , dimension(6)   :: pressureReference         =[7.3027497000209d+05,3.0513637366315d+05,5.5776059586966d+04,1.0779469858510d+04,2.1797763851529d+03,9.7587490195865d+01]
  !! Molecular fractions.
  double precision                                                  , dimension(6)   :: fractionMolecularReference=[9.4240651254240d-01,8.7997335645664d-01,6.0556724446457d-01,2.5284645579126d-01,7.2155830553603d-02,4.4439150381712d-03]
  !! Star formation rate surface densities [M☉/Gyr/Mpc²].
  double precision                                                  , dimension(6)   :: rateBlitzReference        =[6.2973793711964d+13,3.3051869160264d+13,7.3823495252989d+12,1.0325201952966d+12,1.0125935377225d+11,7.7738260031637d+08]
  !! Krumholz et al. (2009) molecular fractions, for the full and fast fitting functions, and rate surface densities
  !! [M☉/Gyr/Mpc²].
  double precision                                                  , dimension(5)   :: fractionMolecularKrumholzReference    =[9.6184243210179d-01,9.5198615456030d-01,9.3757069811107d-01,8.3682408513234d-01,5.8319551403055d-01], &
       &                                                                                fractionMolecularKrumholzFastReference=[9.6214522812381d-01,9.5238986128348d-01,9.3809477201003d-01,8.3748636933749d-01,5.9586006438119d-01], &
       &                                                                                rateKrumholzReference                 =[4.2887311608501d+13,3.1153746907986d+13,2.5649133906021d+13,1.1714548304555d+13,4.1776134042290d+12]
  ! Tolerances.
  !! The pressure is compared in physical units, so this comparison is sensitive to the values of the physical constants - GSL,
  !! and therefore Galacticus, has used both 6.673×10⁻¹¹ and 6.6743×10⁻¹¹ N m² kg⁻² for the gravitational constant. Since
  !! P ∝ G^¾ here that is a 5×10⁻⁴ effect, so we allow 2×10⁻³. The error which this test exists to catch - the characteristic
  !! pressure being used as though it were a logarithm - is a factor of 7,700.
  double precision                                                  , parameter      :: tolerancePressure         =2.0d-3
  !! The molecular fraction of the Blitz & Rosolowsky (2006) law is tabulated by a `fastExponentiator` over a pressure ratio of
  !! [0,1] at a spacing of 10⁻³, which is exact above unity and accurate to 2×10⁻⁶ down to a ratio of 0.06 - but only to 9×10⁻⁴
  !! at the ratio of 0.003 reached at six scale lengths. Relative to the reference the rate is subject to the pressure
  !! uncertainty above also, amplified by α(1-f_H₂) ≤ 0.92.
  double precision                                                  , parameter      :: toleranceBlitz            =3.0d-3
  !! The logarithmic slope of R_mol with pressure is formed from the class' own pressures, so only the tabulation enters.
  double precision                                                  , parameter      :: toleranceSlope            =1.0d-5
  !! The Krumholz et al. (2009) law involves no physical constants beyond the definition of a parsec, and its molecular fraction
  !! is tabulated in s over [0,10] with 1000 points, which is accurate to 6×10⁻⁶ at the values of s reached here.
  double precision                                                  , parameter      :: toleranceKrumholz         =1.0d-4
  !! The gas surface density of an exponential disk is analytic, so agreement here is limited only by round-off.
  double precision                                                  , parameter      :: toleranceSurfaceDensity   =1.0d-9
  !! In the limiting cases the molecular fraction differs from unity by 10⁻⁴ (Blitz & Rosolowsky) and 10⁻⁴ (Krumholz et al.).
  double precision                                                  , parameter      :: toleranceLimit            =1.0d-3
  double precision                                                  , dimension(6)   :: surfaceDensityGas                     , pressureRatio                        , &
       &                                                                                rateBlitz                             , fractionMolecular
  double precision                                                  , dimension(5)   :: rateKrumholz                          , rateKrumholzFast                     , &
       &                                                                                fractionMolecularKrumholz             , fractionMolecularKrumholzFast
  double precision                                                                   :: frequencyNormalizationBlitz           , surfaceDensityHydrogen               , &
       &                                                                                radius                                , rateFullyMolecular                   , &
       &                                                                                slopePressure                         , surfaceDensityGas_
  integer                                                                            :: i

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Molecular fraction star formation rate surface densities")
  parameters=inputParameters('testSuite/parameters/starFormation/molecularFractions.xml')
  call eventsHooksInitialize            (          )
  call Functions_Global_Set             (          )
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  ! Build the star formation rate surface density objects. Both are built from an empty parameter set, so that they take the
  ! default parameters of their classes - which are therefore tested here along with the formulas. The exception is the fast
  ! molecular fraction fitting function, which is not the default and so is selected explicitly.
  parametersDefault=inputParameters()
  allocate(blitz2006_       )
  allocate(krumholz2009_    )
  allocate(krumholz2009Fast_)
  blitz2006_       =starFormationRateSurfaceDensityDisksBlitz2006   (parametersDefault)
  krumholz2009_    =starFormationRateSurfaceDensityDisksKrumholz2009(parametersDefault)
  krumholz2009Fast_=starFormationRateSurfaceDensityDisksKrumholz2009(                                                &
       &                                                             frequencyStarFormation        =frequencyStarFormation, &
       &                                                             clumpingFactorMolecularComplex=clumpingFactor        , &
       &                                                             molecularFractionFast         =.true.                  &
       &                                                            )
  ! Build a node with an exponential disk of gas and stars at Solar metallicity. At Solar metallicity Galacticus' hydrogen mass
  ! fraction is exactly its Solar value, `hydrogenByMassSolar`, which is what the reference script assumes. Two further nodes
  ! carry the same disk compressed by factors of ten and twenty, used for the fully molecular limits. These must be separate
  ! nodes: a node memoizes the mass distributions built from its components, so simply resetting the scale length of this disk
  ! would leave the original surface density profile in place.
  call buildNode(node                ,radiusScale        )
  call buildNode(nodeCompressedTen   ,radiusScale/10.0d0 )
  call buildNode(nodeCompressedTwenty,radiusScale/20.0d0 )
  tree%nodeBase => node
  basic         => node%basic(autoCreate=.true.)
  disk          => node%disk (autoCreate=.true.)
  ! The star formation frequency normalization, converted from yr⁻¹ to Gyr⁻¹.
  frequencyNormalizationBlitz=frequencyNormalization*giga
  ! Evaluate the Blitz & Rosolowsky (2006) law.
  call Unit_Tests_Begin_Group("Blitz & Rosolowsky (2006)")
  do i=1,countRadiiBlitz
     radius                =+radiiBlitz(i)                                                  &
          &                 *radiusScale
     rateBlitz        (i)  = blitz2006_%rate         (node,radius                          )
     ! The pressure ratio, P_ext/P₀, is available directly from the class, along with the gas surface density.
     pressureRatio    (i)  = blitz2006_%pressureRatio(node,radius,surfaceDensityGas=surfaceDensityGas_)
     surfaceDensityGas(i)  = surfaceDensityGas_
     ! Recover the molecular fraction from the rate surface density.
     surfaceDensityHydrogen=+hydrogenByMassSolar                                            &
          &                 *surfaceDensityGas(i)
     fractionMolecular(i)  =+rateBlitz        (i)                                           &
          &                 /surfaceDensityHydrogen                                         &
          &                 /frequencyNormalizationBlitz                                    &
          &                 /(                                                              &
          &                   +1.0d0                                                        &
          &                   +(                                                            &
          &                     +surfaceDensityHydrogen                                     &
          &                     /surfaceDensityCritical                                     &
          &                     /mega                  **2                                  &
          &                    )                       **surfaceDensityExponent             &
          &                  )
  end do
  ! The gas surface density must be that of an exponential disk - this simply confirms that the disk built here is the disk
  ! assumed by the reference script.
  call Assert('gas surface density is exponential',surfaceDensityGas/mega**2            ,surfaceDensityGasReference,relTol=toleranceSurfaceDensity)
  ! The midplane pressure, in physical units, fixes the characteristic pressure of Table 2 of Blitz & Rosolowsky (2006).
  call Assert('midplane pressure [K/cm³]'         ,pressureRatio*pressureCharacteristic  ,pressureReference         ,relTol=tolerancePressure      )
  ! The molecular fraction must be R_mol/(1+R_mol), not min(R_mol,1).
  call Assert('molecular fraction'                ,fractionMolecular                     ,fractionMolecularReference,relTol=toleranceBlitz         )
  ! The full star formation rate surface density.
  call Assert('rate surface density'              ,rateBlitz                             ,rateBlitzReference        ,relTol=toleranceBlitz         )
  ! The ratio of molecular to atomic hydrogen must scale as P_ext^α. This is measured between the fourth and fifth radii, where
  ! the molecular fraction is small enough for R_mol to be well determined, but the pressure ratio is still large enough that the
  ! tabulated exponentiation is accurate.
  if (fractionMolecular(4) < 1.0d0 .and. fractionMolecular(5) < 1.0d0) then
     slopePressure=+log(                                                                     &
          &             +(fractionMolecular(4)/(1.0d0-fractionMolecular(4)))                 &
          &             /(fractionMolecular(5)/(1.0d0-fractionMolecular(5)))                 &
          &            )                                                                     &
          &        /log(                                                                     &
          &             +pressureRatio     (4)                                               &
          &             /pressureRatio     (5)                                               &
          &            )
  else
     ! A molecular fraction of precisely unity means that R_mol can not be recovered - which is itself a failure of the model, so
     ! report a slope which can not match the exponent rather than dividing by zero.
     slopePressure=-1.0d0
  end if
  call Assert('R_mol ∝ P_ext^α'                   ,slopePressure                         ,pressureExponent          ,relTol=toleranceSlope         )
  ! In the high pressure limit the disk must be fully molecular. Compressing the disk by a factor of ten raises the midplane
  ! pressure by four orders of magnitude, leaving a molecular fraction of 1-1.8×10⁻⁵.
  radius                =+radiiBlitz(1)                                                      &
       &                 *radiusScale                                                        &
       &                 /10.0d0
  rateBlitz         (1) = blitz2006_%rate         (nodeCompressedTen,radius                  )
  pressureRatio     (1) = blitz2006_%pressureRatio(nodeCompressedTen,radius,surfaceDensityGas=surfaceDensityGas_)
  surfaceDensityHydrogen=+hydrogenByMassSolar                                                &
       &                 *surfaceDensityGas_
  rateFullyMolecular    =+surfaceDensityHydrogen                                             &
       &                 *frequencyNormalizationBlitz                                        &
       &                 *(                                                                  &
       &                   +1.0d0                                                            &
       &                   +(                                                                &
       &                     +surfaceDensityHydrogen                                         &
       &                     /surfaceDensityCritical                                         &
       &                     /mega                  **2                                      &
       &                    )                       **surfaceDensityExponent                 &
       &                  )
  call Assert('fully molecular at high pressure'  ,rateBlitz(1)                           ,rateFullyMolecular        ,relTol=toleranceLimit         )
  call Unit_Tests_End_Group()
  ! Evaluate the Krumholz, McKee & Tumlinson (2009) law.
  call Unit_Tests_Begin_Group("Krumholz, McKee & Tumlinson (2009)")
  do i=1,countRadiiKrumholz
     radius                          =+radiiKrumholz(i)                                                       &
          &                           *radiusScale
     rateKrumholz                 (i)= krumholz2009_    %rate(node,radius)
     rateKrumholzFast             (i)= krumholz2009Fast_%rate(node,radius)
     ! Recover the molecular fraction from the rate surface density. The surface density factor of eqn. (10) of Krumholz et al.
     ! (2009) is evaluated here from the exponential disk directly.
     surfaceDensityGas_              =+massGas                                                                &
          &                           /2.0d0                                                                  &
          &                           /Pi                                                                     &
          &                           /radiusScale                 **2                                        &
          &                           *exp(-radiiKrumholz(i))
     fractionMolecularKrumholz    (i)=+rateKrumholz    (i)                                                    &
          &                           /frequencyStarFormation                                                 &
          &                           /surfaceDensityGas_                                                     &
          &                           /(                                                                      &
          &                             +surfaceDensityGas_                                                   &
          &                             /surfaceDensityTransition                                             &
          &                            )**sign(0.33d0,surfaceDensityGas_-surfaceDensityTransition)
     fractionMolecularKrumholzFast(i)=+rateKrumholzFast(i)                                                    &
          &                           /frequencyStarFormation                                                 &
          &                           /surfaceDensityGas_                                                     &
          &                           /(                                                                      &
          &                             +surfaceDensityGas_                                                   &
          &                             /surfaceDensityTransition                                             &
          &                            )**sign(0.33d0,surfaceDensityGas_-surfaceDensityTransition)
  end do
  call Assert('molecular fraction'                    ,fractionMolecularKrumholz    ,fractionMolecularKrumholzReference    ,relTol=toleranceKrumholz)
  call Assert('molecular fraction (fast fit)'         ,fractionMolecularKrumholzFast,fractionMolecularKrumholzFastReference,relTol=toleranceKrumholz)
  call Assert('rate surface density'                  ,rateKrumholz                 ,rateKrumholzReference                 ,relTol=toleranceKrumholz)
  ! Where the total gas surface density equals the transition surface density of 85 M☉/pc² the surface density factor of eqn. (10)
  ! is unity, so that the rate reduces to ν f_H₂ Σ_gas. This holds for the *total* gas surface density, not for the surface
  ! density of hydrogen alone.
  call Assert('rate is ν f Σ_gas at Σ_gas = 85 M☉/pc²',rateKrumholz(2),frequencyStarFormation*fractionMolecularKrumholzReference(2)*surfaceDensityTransition,relTol=toleranceKrumholz)
  ! At zero metallicity there is no dust on which molecules can form, and so no star formation.
  call abundancesGas %metallicitySet  (0.0d0                )
  call disk          %abundancesGasSet(abundancesGas        )
  call krumholz2009_ %calculationReset(node,node%uniqueID() )
  call Assert('no star formation at zero metallicity'  ,krumholz2009_%rate(node,radiusScale),0.0d0)
  call abundancesGas %metallicitySet  (metallicitySolar*massGas)
  call disk          %abundancesGasSet(abundancesGas           )
  call krumholz2009_ %calculationReset(node,node%uniqueID()    )
  ! In the high surface density limit the gas is fully molecular, and the rate reduces to ν Σ_gas (Σ_gas/85 M☉ pc⁻²)^0.33. Here
  ! the disk is compressed by a factor of twenty, leaving a molecular fraction of 1-9.6×10⁻⁵.
  surfaceDensityGas_=+massGas                                                                                 &
       &             /2.0d0                                                                                   &
       &             /Pi                                                                                      &
       &             /(radiusScale/20.0d0)     **2                                                            &
       &             *exp(-radiiKrumholz(1))
  call Assert('fully molecular at high surface density',krumholz2009_%rate(nodeCompressedTwenty,radiiKrumholz(1)*radiusScale/20.0d0),frequencyStarFormation*surfaceDensityGas_*(surfaceDensityGas_/surfaceDensityTransition)**0.33d0,relTol=toleranceLimit)
  ! In the low surface density limit the molecular fraction, and so the rate, vanishes.
  call Assert('vanishing molecular fraction at low surface density',krumholz2009_%rate(node,20.0d0*radiusScale),0.0d0,absTol=1.0d-3)
  call Unit_Tests_End_Group()
  ! Clean up.
  call node                %destroy()
  call nodeCompressedTen   %destroy()
  call nodeCompressedTwenty%destroy()
  deallocate(blitz2006_       )
  deallocate(krumholz2009_    )
  deallocate(krumholz2009Fast_)
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()

contains

  subroutine buildNode(node_,radiusDisk_)
    !!{RST
    Build a node containing an exponential disk of gas and stars at Solar metallicity, with the given disk scale length.
    !!}
    implicit none
    type            (treeNode         ), intent(  out), pointer :: node_
    double precision                   , intent(in   )          :: radiusDisk_
    class           (nodeComponentBasic)               , pointer :: basic_
    class           (nodeComponentDisk )               , pointer :: disk_
    type            (abundances       )                         :: abundancesGas_

    node_  => treeNode(hostTree=tree                )
    basic_ => node_   %basic   (autoCreate=.true.)
    disk_  => node_   %disk    (autoCreate=.true.)
    call abundancesGas_%metallicitySet  (metallicitySolar*massGas)
    call basic_        %massSet         (massGas+massStellar     )
    call basic_        %timeSet         (13.8d0                  )
    call disk_         %massGasSet      (massGas                 )
    call disk_         %massStellarSet  (massStellar             )
    call disk_         %radiusSet       (radiusDisk_             )
    call disk_         %abundancesGasSet(abundancesGas_          )
    return
  end subroutine buildNode

end program Test_Star_Formation_Molecular_Fractions

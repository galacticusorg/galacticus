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

program Tests_Merger_Tree_Branching_CDM
  !!{RST
  Tests of merger tree branching probability bounds for a :term:`CDM` power spectrum.

  The companion tests in ``tests.merger_tree_branching.exe`` use a power-law power spectrum, for which the logarithmic gradient of
  the mass variance, :math:`\alpha`, is constant. That makes the bounds on the branching probability exact, which is what gives
  those tests their analytic targets---but it also means they never exercise the case the tree builder actually meets, in which
  :math:`\alpha` varies with mass and the bounds are genuine inequalities. These tests use a :term:`CDM` power spectrum instead,
  and check the property on which the tree building algorithm depends: candidate branching events are drawn at the rate given by
  the upper bound and then rejected using the exact rate, so an upper bound which is not one biases every tree built, silently.

  They also probe the poles of the hypergeometric functions used to evaluate the bounds. The third parameter of those functions
  vanishes whenever the effective exponent :math:`\gamma_1-1/\alpha` reaches an odd integer :math:`\ge 3`, and for a :term:`CDM`
  power spectrum that happens at halo masses which are entirely ordinary.
  !!}
  use            :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower          , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use            :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use            :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use            :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use            :: Display                             , only : displayVerbositySet                            , verbosityLevelWorking
  use            :: Error                               , only : Error_Handler_Register
  use            :: Events_Hooks                        , only : eventsHooksInitialize
  use            :: Galacticus_Nodes                    , only : treeNode
  use            :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use            :: Merger_Tree_Branching               , only : mergerTreeBranchingProbabilityParkinsonColeHelly, mergerTreeBranchingBoundLower                              , mergerTreeBranchingBoundUpper
  use            :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use            :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use            :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionTopHat
  use            :: Transfer_Functions                  , only : transferFunctionEisensteinHu1999
  use            :: Unit_Tests                          , only : Assert                                         , Unit_Tests_Begin_Group                                     , Unit_Tests_End_Group         , Unit_Tests_Finish, &
       &                                                         compareLessThanOrEqual
  implicit none
  type            (cosmologyParametersSimple                                   )                :: cosmologyParametersSimple_
  type            (cosmologyFunctionsMatterLambda                              )                :: cosmologyFunctionsMatterLambda_
  type            (cosmologicalMassVarianceFilteredPower                       )                :: cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionTopHat                           )                :: powerSpectrumWindowFunctionTopHat_
  type            (powerSpectrumPrimordialPowerLaw                             )                :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionEisensteinHu1999                            )                :: transferFunctionEisensteinHu1999_
  type            (powerSpectrumPrimordialTransferredSimple                    )                :: powerSpectrumPrimordialTransferredSimple_
  type            (linearGrowthCollisionlessMatter                             )                :: linearGrowthCollisionlessMatter_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                :: criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_
  type            (darkMatterParticleCDM                                       )                :: darkMatterParticleCDM_
  type            (treeNode                                                    )                :: node
  ! Configurations over which the bounds are checked. All four combinations of the two options which change how the bounds are
  ! evaluated are covered: whether CDM assumptions are made, and whether the hypergeometric functions are tabulated.
  integer                                                                       , parameter     :: countConfiguration=4                 , countMass         =13, &
       &                                                                                           countRedshift     =3                 , countPole         = 3, &
       &                                                                                           countOffset       =9
  logical                                                                       , dimension(countConfiguration), parameter :: cdmAssumptionsConfiguration=[.true. ,.true. ,.false.,.false.], &
       &                                                                                                                     tabulateConfiguration      =[.true. ,.false.,.true. ,.false.]
  double precision                                                              , dimension(countRedshift      ), parameter :: redshift                   =[0.0d0,1.0d0,3.0d0]
  ! Poles of the hypergeometric functions, and the halo masses at which each is reached for a value of γ₁ close to the value
  ! favored by Parkinson et al. (2008). The masses are chosen only so that the required γ₁ comes out physically plausible - the
  ! test below asserts that it does, which is what establishes that these poles are reachable in practice rather than only in
  ! principle. The exact value of γ₁ which places the effective exponent on the pole is computed from the mass variance itself, so
  ! these need be no more than approximately right.
  double precision                                                              , dimension(countPole          ), parameter :: poleTarget                 =[3.0d0 ,5.0d0 ,7.0d0 ]                , &
       &                                                                                                                       massPole                   =[1.0d16,1.0d14,1.0d12]
  ! Offsets applied to the effective exponent about each pole. The innermost of these lie deep within the neighborhood in which the
  ! CDM assumptions are declined, the outermost well outside it, and one pair just outside its edge. The innermost pair carry most
  ! of the discriminating power of the precision-independence test below: the precision lost to the cancellation grows as the
  ! reciprocal of the offset, so it is there that a failure to decline the CDM assumptions shows up most strongly.
  double precision                                                              , dimension(countOffset        ), parameter :: offsetPole                 =[-1.0d-1,-4.0d-3,-1.0d-4,-1.0d-5,0.0d0,+1.0d-5,+1.0d-4,+4.0d-3,+1.0d-1]
  double precision                                                              , parameter     :: massResolutionFraction=1.0d-3
  type            (mergerTreeBranchingProbabilityParkinsonColeHelly            ), dimension(countConfiguration        ) :: branchingConfiguration_
  type            (mergerTreeBranchingProbabilityParkinsonColeHelly            ), dimension(countPole*countOffset     ) :: branchingPole_
  type            (mergerTreeBranchingProbabilityParkinsonColeHelly            ), dimension(countPole*countOffset     ) :: branchingPolePrecise_
  double precision                                                              , dimension(countMass                 ) :: mass                                  , probabilityMass    , &
       &                                                                                                                  boundLowerMass                        , boundUpperMass     , &
       &                                                                                                                  boundUpperDirect
  double precision                                                              , dimension(countOffset               ) :: probabilityPole                       , boundUpperPole     , &
       &                                                                                                                  boundUpperPolePrecise
  double precision                                                              , dimension(countPole*countOffset     ) :: gamma1Pole
  double precision                                                              , dimension(countPole                 ) :: gamma1PoleRequired
  double precision                                                                               :: time                                  , expansionFactor    , &
       &                                                                                            deltaCritical_                        , massResolution_    , &
       &                                                                                            sigmaScratch                          , alphaHalf
  integer                                                                                        :: i                                     , j                  , &
       &                                                                                            k                                     , iConfiguration
  ! Long enough for the composed group labels below, whose subscripted characters occupy several bytes each.
  character       (len=64                                                      )                  :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelWorking)
  ! Register error handlers, so that Galacticus' own GSL error handler is installed. Without it the default GSL handler aborts
  ! irrespective of any `status` argument requesting that a failure be returned instead, and the recovery from a failed evaluation
  ! of a hypergeometric function which these tests are meant to exercise cannot be reached.
  call Error_Handler_Register()
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Build a ΛCDM cosmology with an Eisenstein-Hu transfer function and a top-hat filter, so that α = dlnσ/dlnM varies with mass as
  ! it does in any real model.
  darkMatterParticleCDM_                                       =darkMatterParticleCDM                                        (                                                                                          &
       &                                                                                                                     )
  cosmologyParametersSimple_                                   =cosmologyParametersSimple                                    (                                                                                          &
       &                                                                                                                      OmegaMatter                             = 0.3000d0                                      , &
       &                                                                                                                      OmegaBaryon                             = 0.0450d0                                      , &
       &                                                                                                                      OmegaDarkEnergy                         = 0.7000d0                                      , &
       &                                                                                                                      temperatureCMB                          = 2.7250d0                                      , &
       &                                                                                                                      HubbleConstant                          =70.0000d0                                        &
       &                                                                                                                     )
  cosmologyFunctionsMatterLambda_                              =cosmologyFunctionsMatterLambda                               (                                                                                          &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                       &
       &                                                                                                                     )
  linearGrowthCollisionlessMatter_                             =linearGrowthCollisionlessMatter                              (                                                                                          &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                  &
       &                                                                                                                     )
  powerSpectrumPrimordialPowerLaw_                             =powerSpectrumPrimordialPowerLaw                              (                                                                                          &
       &                                                                                                                      index_                                  =+0.96d0                                        , &
       &                                                                                                                      running                                 =+0.00d0                                        , &
       &                                                                                                                      runningRunning                          =+0.00d0                                        , &
       &                                                                                                                      wavenumberReference                     =+1.00d0                                        , &
       &                                                                                                                      runningSmallScalesOnly                  =.false.                                          &
       &                                                                                                                     )
  transferFunctionEisensteinHu1999_                            =transferFunctionEisensteinHu1999                             (                                                                                          &
       &                                                                                                                      neutrinoNumberEffective                 =3.046d0                                        , &
       &                                                                                                                      neutrinoMassSummed                      =0.000d0                                        , &
       &                                                                                                                      darkMatterParticle_                     =darkMatterParticleCDM_                         , &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                  &
       &                                                                                                                     )
  powerSpectrumPrimordialTransferredSimple_                    =powerSpectrumPrimordialTransferredSimple                     (                                                                                          &
       &                                                                                                                      powerSpectrumPrimordial_                =powerSpectrumPrimordialPowerLaw_               , &
       &                                                                                                                      transferFunction_                       =transferFunctionEisensteinHu1999_              , &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_                 &
       &                                                                                                                     )
  powerSpectrumWindowFunctionTopHat_                           =powerSpectrumWindowFunctionTopHat                            (                                                                                          &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                       &
       &                                                                                                                     )
  cosmologicalMassVarianceFilteredPower_                       =cosmologicalMassVarianceFilteredPower                        (                                                                                          &
       &                                                                                                                      sigma8                                  =0.8d+0                                         , &
       &                                                                                                                      tolerance                               =1.0d-4                                         , &
       &                                                                                                                      toleranceTopHat                         =1.0d-4                                         , &
       &                                                                                                                      rootVarianceLogarithmicGradientTolerance=1.0d-6                                         , &
       &                                                                                                                      integrationFailureIsFatal               =.true.                                         , &
       &                                                                                                                      storeTabulations                        =.true.                                         , &
       &                                                                                                                      nonMonotonicIsFatal                     =.true.                                         , &
       &                                                                                                                      monotonicInterpolation                  =.false.                                        , &
       &                                                                                                                      truncateAtParticleHorizon               =.false.                                        , &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                , &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_               , &
       &                                                                                                                      powerSpectrumPrimordialTransferred_     =powerSpectrumPrimordialTransferredSimple_      , &
       &                                                                                                                      powerSpectrumWindowFunction_            =powerSpectrumWindowFunctionTopHat_               &
       &                                                                                                                     )
  criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_=criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt (                                                                                          &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_               , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                , &
       &                                                                                                                      cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_         , &
       &                                                                                                                      darkMatterParticle_                     =darkMatterParticleCDM_                         , &
       &                                                                                                                      tableStore                              =.true.                                           &
       &                                                                                                                     )
  ! Build one branching probability object per configuration, using the parameter values favored by Parkinson et al. (2008).
  do iConfiguration=1,countConfiguration
     branchingConfiguration_(iConfiguration)=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                                  &
          &                                                                                   G0                       =0.57d0                                                               , &
          &                                                                                   gamma1                   =0.38d0                                                               , &
          &                                                                                   gamma2                   =-0.01d0                                                              , &
          &                                                                                   accuracyFirstOrder       =1.0d-1                                                               , &
          &                                                                                   precisionHypergeometric  =1.0d-6                                                               , &
          &                                                                                   hypergeometricTabulate   =tabulateConfiguration      (iConfiguration)                          , &
          &                                                                                   cdmAssumptions           =cdmAssumptionsConfiguration(iConfiguration)                          , &
          &                                                                                   tolerateRoundOffErrors   =.false.                                                              , &
          &                                                                                   cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                               , &
          &                                                                                   criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_          &
          &                                                                                  )
  end do
  ! Masses at which the bounds are checked, spanning the range over which merger trees are built.
  do i=1,countMass
     mass(i)=10.0d0**(9.0d0+6.0d0*dble(i-1)/dble(countMass-1))
  end do

  call Unit_Tests_Begin_Group("Merger tree branching (CDM power spectrum)")

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! The bounds must bound.
  !
  ! With α varying with mass the bounds are genuine inequalities, and this is the property the tree builder relies upon: an upper
  ! bound which falls below the true probability does not fail, it biases.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Bounds bracket the probability")
  do i=1,countRedshift
     expansionFactor=cosmologyFunctionsMatterLambda_%expansionFactorFromRedshift(redshift(i))
     time           =cosmologyFunctionsMatterLambda_%cosmicTime                 (expansionFactor)
     write (label,'(a,f3.1)') "z = ",redshift(i)
     call Unit_Tests_Begin_Group(trim(label))
     do iConfiguration=1,countConfiguration
        do j=1,countMass
           massResolution_   =mass(j)*massResolutionFraction
           deltaCritical_    =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_%value(time=time,mass=mass(j),node=node)
           probabilityMass(j)=branchingConfiguration_(iConfiguration)%probability     (mass(j),deltaCritical_,time,massResolution_                              ,node)
           boundLowerMass (j)=branchingConfiguration_(iConfiguration)%probabilityBound(mass(j),deltaCritical_,time,massResolution_,mergerTreeBranchingBoundLower,node)
           boundUpperMass (j)=branchingConfiguration_(iConfiguration)%probabilityBound(mass(j),deltaCritical_,time,massResolution_,mergerTreeBranchingBoundUpper,node)
        end do
        if (cdmAssumptionsConfiguration(iConfiguration)) then
           label="CDM assumptions"
        else
           label="no CDM assumptions"
        end if
        if (tabulateConfiguration(iConfiguration)) then
           label=trim(label)//"; tabulated ₂F₁"
        else
           label=trim(label)//"; computed ₂F₁"
        end if
        call Unit_Tests_Begin_Group(trim(label))
        call Assert('lower bound ≤ probability',boundLowerMass ,probabilityMass,compareLessThanOrEqual)
        call Assert('probability ≤ upper bound',probabilityMass,boundUpperMass ,compareLessThanOrEqual)
        call Unit_Tests_End_Group  (                          )
     end do
     call Unit_Tests_End_Group  (         )
  end do
  call Unit_Tests_End_Group  (            )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! Tabulating the hypergeometric functions must not change the bound.
  !
  ! The upper bound is tabulated against mass under the CDM assumptions, and each entry of that tabulation is evaluated at its own
  ! effective exponent. A single entry spoiled by the cancellation which afflicts those functions near their poles would corrupt
  ! the interpolation over a whole interval in mass, which comparing against the directly computed bound will reveal.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Tabulated and computed bounds agree")
  time=cosmologyFunctionsMatterLambda_%cosmicTime(1.0d0)
  do j=1,countMass
     massResolution_  =mass(j)*massResolutionFraction
     deltaCritical_   =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_%value(time=time,mass=mass(j),node=node)
     boundUpperMass(j)=branchingConfiguration_(1)%probabilityBound(mass(j),deltaCritical_,time,massResolution_,mergerTreeBranchingBoundUpper,node)
     boundUpperDirect(j)=branchingConfiguration_(2)%probabilityBound(mass(j),deltaCritical_,time,massResolution_,mergerTreeBranchingBoundUpper,node)
  end do
  ! The tabulation is interpolated at thirty points per decade in mass, which is what limits the agreement here.
  call Assert('upper bound, tabulated against computed',boundUpperMass,boundUpperDirect,relTol=3.0d-2)
  call Unit_Tests_End_Group  (                                        )

  ! ----------------------------------------------------------------------------------------------------------------------------------
  ! Behavior at the poles of the hypergeometric functions.
  !
  ! The effective exponent γ₁-1/α reaches an odd integer at halo masses which are entirely ordinary, and there the two
  ! hypergeometric factors from which the bound is formed both diverge while their difference stays finite. For each pole the value
  ! of γ₁ which places the exponent exactly upon it is computed from the mass variance, the exponent is then displaced about it,
  ! and the bound is required to remain a bound throughout - and, outside the neighborhood in which the CDM assumptions are
  ! declined, to be independent of the precision to which the hypergeometric functions are evaluated. That last is the sensitive
  ! test: it is precisely the cancellation which makes the result depend on that precision.
  ! ----------------------------------------------------------------------------------------------------------------------------------
  call Unit_Tests_Begin_Group("Hypergeometric poles")
  do i=1,countPole
     ! The upper bound is evaluated at the effective exponent formed from α at half the parent mass, so that is the α which must
     ! be used to place the exponent on the pole.
     call cosmologicalMassVarianceFilteredPower_%rootVarianceAndLogarithmicGradient(0.5d0*massPole(i),time,sigmaScratch,alphaHalf)
     gamma1PoleRequired(i)=poleTarget(i)+1.0d0/alphaHalf
     do j=1,countOffset
        k            =(i-1)*countOffset+j
        gamma1Pole(k)=gamma1PoleRequired(i)+offsetPole(j)
     end do
  end do
  ! Establish that these poles are reachable with a physically plausible γ₁ - that is, that this is not merely a mathematical
  ! curiosity but a configuration a user could arrive at. Parkinson et al. (2008) favor γ₁ = 0.38.
  call Assert('γ₁ placing the exponent on a pole is physically plausible',abs(gamma1PoleRequired),spread(1.5d0,1,countPole),compareLessThanOrEqual)
  do k=1,countPole*countOffset
     branchingPole_       (k)=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                          &
          &                                                                    G0                       =0.57d0                                                        , &
          &                                                                    gamma1                   =gamma1Pole(k)                                                 , &
          &                                                                    gamma2                   =-0.01d0                                                       , &
          &                                                                    accuracyFirstOrder       =1.0d-1                                                        , &
          &                                                                    precisionHypergeometric  =1.0d-6                                                        , &
          &                                                                    hypergeometricTabulate   =.false.                                                       , &
          &                                                                    cdmAssumptions           =.true.                                                        , &
          &                                                                    tolerateRoundOffErrors   =.false.                                                       , &
          &                                                                    cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                        , &
          &                                                                    criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_   &
          &                                                                   )
     branchingPolePrecise_(k)=mergerTreeBranchingProbabilityParkinsonColeHelly(                                                                                          &
          &                                                                    G0                       =0.57d0                                                        , &
          &                                                                    gamma1                   =gamma1Pole(k)                                                 , &
          &                                                                    gamma2                   =-0.01d0                                                       , &
          &                                                                    accuracyFirstOrder       =1.0d-1                                                        , &
          &                                                                    precisionHypergeometric  =1.0d-12                                                       , &
          &                                                                    hypergeometricTabulate   =.false.                                                       , &
          &                                                                    cdmAssumptions           =.true.                                                        , &
          &                                                                    tolerateRoundOffErrors   =.false.                                                       , &
          &                                                                    cosmologicalMassVariance_=cosmologicalMassVarianceFilteredPower_                        , &
          &                                                                    criticalOverdensity_     =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_   &
          &                                                                   )
  end do
  do i=1,countPole
     massResolution_=massPole(i)*massResolutionFraction
     deltaCritical_ =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_%value(time=time,mass=massPole(i),node=node)
     do j=1,countOffset
        k                       =(i-1)*countOffset+j
        probabilityPole      (j)=branchingPole_       (k)%probability     (massPole(i),deltaCritical_,time,massResolution_                              ,node)
        boundUpperPole       (j)=branchingPole_       (k)%probabilityBound(massPole(i),deltaCritical_,time,massResolution_,mergerTreeBranchingBoundUpper,node)
        boundUpperPolePrecise(j)=branchingPolePrecise_(k)%probabilityBound(massPole(i),deltaCritical_,time,massResolution_,mergerTreeBranchingBoundUpper,node)
     end do
     write (label,'(a,f3.1,a,es7.1e2)') "γ₁-1/α = ",poleTarget(i)," at M = ",massPole(i)
     call Unit_Tests_Begin_Group(trim(label))
     call Assert('probability ≤ upper bound'                 ,probabilityPole,boundUpperPole       ,compareLessThanOrEqual)
     ! The two objects compared here differ only in the precision to which they evaluate the hypergeometric functions, so their
     ! results always differ by of order the looser of those precisions, 10⁻⁶ - that is the floor on this comparison, and no
     ! tolerance below a few times it is meaningful. Against that, declining to guard the poles was measured to move these bounds
     ! by 1.3-1.7·10⁻⁴. The tolerance is set between the two: some four times the largest difference measured with the guard in
     ! place, and four times smaller than the smallest difference measured without it.
     call Assert('upper bound independent of ₂F₁ precision'  ,boundUpperPole ,boundUpperPolePrecise,relTol=3.0d-5         )
     call Unit_Tests_End_Group  (                                            )
  end do
  call Unit_Tests_End_Group  (                                               )

  call Unit_Tests_End_Group  (                                               )
  call Unit_Tests_Finish     (                                               )
end program Tests_Merger_Tree_Branching_CDM

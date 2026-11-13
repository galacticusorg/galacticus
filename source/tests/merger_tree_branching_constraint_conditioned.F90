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

program Tests_Merger_Tree_Branching_Constraint_Conditioned
  !!{RST
  Tests of the merger tree branching probability class which conditions a branch on reaching a specified progenitor. These use a
  critical density universe and an :math:`n=-1` power spectrum filtered with a sharp :math:`k`-space window function, for which
  :math:`S(M) \propto M^{-2/3}`, together with the analytic first crossing rate for a constant barrier. The reference rates are
  computed here directly from those analytic forms, per unit critical overdensity scaled to the present day, and so do not rely on
  the conversions made within the class being tested.
  !!}
  use, intrinsic :: ISO_C_Binding                       , only : c_long
  use            :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower                       , &
  &                                                              criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use            :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use            :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use            :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use            :: Display                             , only : displayVerbositySet                                         , &
  &                                                              verbosityLevelWorking
  use            :: Error                               , only : Error_Handler_Register
  use            :: Events_Hooks                        , only : eventsHooksInitialize
  use            :: Excursion_Sets_Barriers             , only : excursionSetBarrierCriticalOverdensity
  use            :: Excursion_Sets_First_Crossings      , only : excursionSetFirstCrossingLinearBarrier
  use            :: Functions_Global_Utilities          , only : Functions_Global_Set
  use            :: Galacticus_Nodes                    , only : mergerTree                                                  , &
  &                                                              treeNode
  use            :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use            :: Merger_Tree_Branching               , only : mergerTreeBranchingProbabilityConstraintConditioned         , &
  &                                                              mergerTreeBranchingProbabilityGeneralizedPressSchechter
  use            :: Merger_Tree_Branching_Modifiers     , only : mergerTreeBranchingProbabilityModifierIdentity
  use            :: Merger_Trees_Build_Mass_Resolution  , only : mergerTreeMassResolutionFixed
  use            :: Numerical_Constants_Math            , only : Pi
  use            :: Numerical_Random_Numbers            , only : randomNumberGeneratorGSL
  use            :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use            :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use            :: Power_Spectrum_Normalizers          , only : powerSpectrumNormalizerSigma8
  use            :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionSharpKSpace
  use            :: Transfer_Functions                  , only : transferFunctionIdentity
  use            :: Unit_Tests                          , only : Assert                                                      , &
  &                                                              Unit_Tests_Begin_Group                                      , &
  &                                                              Unit_Tests_End_Group                                        , &
  &                                                              Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple                                   )                                                :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda                              )                                                :: cosmologyFunctions_
  type            (linearGrowthCollisionlessMatter                             )                                                :: linearGrowth_
  type            (powerSpectrumPrimordialPowerLaw                             )                                                :: powerSpectrumPrimordial_
  type            (transferFunctionIdentity                                    )                                                :: transferFunction_
  type            (powerSpectrumNormalizerSigma8                               )                                                :: powerSpectrumNormalizer_
  type            (powerSpectrumPrimordialTransferredSimple                    )                                                :: powerSpectrumPrimordialTransferred_
  type            (powerSpectrumWindowFunctionSharpKSpace                      )                                                :: powerSpectrumWindowFunction_
  type            (cosmologicalMassVarianceFilteredPower                       ), target                                        :: cosmologicalMassVariance_
  type            (darkMatterParticleCDM                                       )                                                :: darkMatterParticle_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                                                :: criticalOverdensity_
  type            (excursionSetBarrierCriticalOverdensity                      )                                                :: excursionSetBarrier_
  type            (excursionSetFirstCrossingLinearBarrier                      )                                                :: excursionSetFirstCrossing_
  type            (mergerTreeBranchingProbabilityModifierIdentity              )                                                :: modifier_
  type            (mergerTreeBranchingProbabilityGeneralizedPressSchechter     )                                                :: wrapped
  type            (mergerTreeMassResolutionFixed                               )                                                :: massResolution_
  type            (randomNumberGeneratorGSL                                    )                                                :: randomNumberGenerator_
  type            (mergerTree                                                  ), target                                        :: tree
  type            (treeNode                                                    ), target                                        :: node
  ! Cases use a wide window at two redshifts, and the default window, which needs a finer tabulation.
  integer                                                                       , parameter                                     :: countCase                          =3                                   , &
       &                                                                                                                           countFine                          =20000                               , &
       &                                                                                                                           countParent                        =3                                   , &
       &                                                                                                                           countFraction                      =5                                   , &
       &                                                                                                                           countOffset                        =2
  double precision                                                              , parameter                                     :: massConstraint                     =1.0d12                              , &
       &                                                                                                                           massResolutionFixed                =1.0d10                              , &
       &                                                                                                                           stepDelta                          =0.01d0                              , &
       &                                                                                                                           massMinimum                        =1.0d6
  double precision                                                              , dimension(              countCase), parameter :: redshiftConstraint                 =[1.00d0,3.00d0,1.000d0]             , &
       &                                                                                                                           windowsLogMass                     =[0.30d0,0.30d0,0.050d0]             , &
       &                                                                                                                           stepsLogMass                       =[0.02d0,0.02d0,0.005d0]
  ! Offsets of the critical overdensity from the constraint at which the tabulation is checked. The largest comes first, so that
  ! the table is built only once.
  double precision                                                              , dimension(countOffset            ), parameter :: deltaOffsets                       =[1.00d0,0.30d0]
  ! For the default window the parents include masses just above the window, where the expectation is hardest to converge.
  double precision                                                              , dimension(countParent  ,countCase), parameter :: massesParent                       =reshape([3.0d12,1.0d13,3.0d13,3.0d12,1.0d13,3.0d13,1.1d12,2.0d12,1.0d13],[countParent,countCase])
  ! Each redshift has its own objects, as assigning a constructor result onto an object already in use risks a double finalization
  ! (GCC PR 110626).
  type            (mergerTreeBranchingProbabilityConstraintConditioned         ), dimension(              countCase)            :: conditioned                                                             , &
       &                                                                                                                           conditionedFine
  double precision                                                              , dimension(countFraction          ), parameter :: fractionsSought                    =[1.0d-3,1.0d-2,0.25d0,0.50d0,0.75d0]
  ! Tolerances. The largest deviations attained are 0.16% (first step), 0.30% (backward equation), 0.39% (grid convergence), and
  ! 0.005% (inversion of the cumulative rate). A quadrature that does not resolve the default window fails the grid convergence
  ! assertions by up to 7%. Omitting the growth factor from the tabulated rates fails 36 of these assertions, and omitting it
  ! from only the sub-resolution mass loss fails 5 of them.
  double precision                                                              , parameter                                     :: toleranceFirstStep                 =1.0d-2                              , &
       &                                                                                                                           toleranceResidual                  =2.0d-2                              , &
       &                                                                                                                           toleranceConvergence               =1.0d-2                              , &
       &                                                                                                                           toleranceInverse                   =2.0d-2                              , &
       &                                                                                                                           toleranceExact                     =1.0d-12
  double precision                                                                                                              :: timePresent                                                             , &
       &                                                                                                                           timeConstraint                                                          , &
       &                                                                                                                           deltaConstraint                                                         , &
       &                                                                                                                           varianceConstraint                                                      , &
       &                                                                                                                           growthFactor                                                            , &
       &                                                                                                                           delta                                                                   , &
       &                                                                                                                           time                                                                    , &
       &                                                                                                                           massParent                                                              , &
       &                                                                                                                           expected                                                                , &
       &                                                                                                                           actual                                                                  , &
       &                                                                                                                           H                                                                       , &
       &                                                                                                                           HFine                                                                   , &
       &                                                                                                                           probability                                                             , &
       &                                                                                                                           massBranch                                                              , &
       &                                                                                                                           rateExpected                                                            , &
       &                                                                                                                           massDrawn                                                               , &
       &                                                                                                                           drift                                                                   , &
       &                                                                                                                           gradient                                                                , &
       &                                                                                                                           derivative                                                              , &
       &                                                                                                                           windowLogMass                                                           , &
       &                                                                                                                           stepLogMass
  double precision                                                              , dimension(countParent            )            :: massParents
  integer                                                                                                                       :: i                                                                       , &
       &                                                                                                                           j                                                                       , &
       &                                                                                                                           k
  character       (len=16                                                      )                                                :: label

  call displayVerbositySet     (verbosityLevelWorking)
  call Error_Handler_Register  (                     )
  call eventsHooksInitialize   (                     )
  call Functions_Global_Set    (                     )
  darkMatterParticle_                =darkMatterParticleCDM                         ()
  cosmologyParameters_               =cosmologyParametersSimple                     (OmegaMatter=1.0d0,OmegaBaryon=0.0d0,OmegaDarkEnergy=0.0d0,temperatureCMB=2.78d0,HubbleConstant=70.0d0)
  cosmologyFunctions_                =cosmologyFunctionsMatterLambda                (cosmologyParameters_=cosmologyParameters_)
  linearGrowth_                      =linearGrowthCollisionlessMatter               (cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_)
  powerSpectrumPrimordial_           =powerSpectrumPrimordialPowerLaw               (index_=-1.0d0,running=0.0d0,runningRunning=0.0d0,wavenumberReference=1.0d0,runningSmallScalesOnly=.false.)
  transferFunction_                  =transferFunctionIdentity                      (cosmologyParameters_=cosmologyParameters_,time=13.8d0)
  powerSpectrumNormalizer_           =powerSpectrumNormalizerSigma8                 (sigma8=1.0d0,tolerance=1.0d-4,integrationFailureIsFatal=.true.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,transferFunction_=transferFunction_)
  powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferredSimple      (powerSpectrumPrimordial_=powerSpectrumPrimordial_,transferFunction_=transferFunction_,linearGrowth_=linearGrowth_,powerSpectrumNormalizer_=powerSpectrumNormalizer_)
  powerSpectrumWindowFunction_       =powerSpectrumWindowFunctionSharpKSpace        (cosmologyParameters_=cosmologyParameters_,normalization=0.0d0)
  cosmologicalMassVariance_          =cosmologicalMassVarianceFilteredPower         (tolerance=1.0d-4,toleranceTopHat=1.0d-4,rootVarianceLogarithmicGradientTolerance=1.0d-9,integrationFailureIsFatal=.true.,storeTabulations=.true.,nonMonotonicIsFatal=.true.,monotonicInterpolation=.false.,truncateAtParticleHorizon=.false.,cosmologyParameters_=cosmologyParameters_,cosmologyFunctions_=cosmologyFunctions_,linearGrowth_=linearGrowth_,powerSpectrumPrimordialTransferred_=powerSpectrumPrimordialTransferred_,powerSpectrumWindowFunction_=powerSpectrumWindowFunction_)
  criticalOverdensity_               =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt(linearGrowth_=linearGrowth_,cosmologyFunctions_=cosmologyFunctions_,cosmologicalMassVariance_=cosmologicalMassVariance_,darkMatterParticle_=darkMatterParticle_,tableStore=.true.)
  excursionSetBarrier_               =excursionSetBarrierCriticalOverdensity        (criticalOverdensity_=criticalOverdensity_,cosmologicalMassVariance_=cosmologicalMassVariance_)
  excursionSetFirstCrossing_         =excursionSetFirstCrossingLinearBarrier        (fractionalTimeStep=1.0d-3,excursionSetBarrier_=excursionSetBarrier_,cosmologicalMassVariance_=cosmologicalMassVariance_)
  modifier_                          =mergerTreeBranchingProbabilityModifierIdentity()
  wrapped                            =mergerTreeBranchingProbabilityGeneralizedPressSchechter(deltaStepMaximum=1.0d-1,massMinimum=massMinimum,smoothAccretion=.true.,distributionFunctionLowerHalfOnly=.false.,distributionFunctionNormalize=.true.,cosmologyFunctions_=cosmologyFunctions_,criticalOverdensity_=criticalOverdensity_,cosmologicalMassVariance_=cosmologicalMassVariance_,excursionSetFirstCrossing_=excursionSetFirstCrossing_,mergerTreeBranchingProbabilityModifier_=modifier_)
  massResolution_                    =mergerTreeMassResolutionFixed                 (massResolution=massResolutionFixed)
  randomNumberGenerator_             =randomNumberGeneratorGSL                      (8213_c_long)
  node%hostTree => tree
  timePresent=cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
  call Unit_Tests_Begin_Group("Constraint-conditioned merger tree branching")
  do i=1,countCase
     windowLogMass=windowsLogMass(  i)
     stepLogMass  =stepsLogMass  (  i)
     massParents  =massesParent  (:,i)
     write (label,'(f4.1)') redshiftConstraint(i)
     call Unit_Tests_Begin_Group("z₁ = "//trim(adjustl(label)))
     write (label,'(f5.3)') windowLogMass
     call Unit_Tests_Begin_Group("window = ± "//trim(adjustl(label)))
     timeConstraint    =+cosmologyFunctions_      %cosmicTime(expansionFactor=1.0d0/(1.0d0+redshiftConstraint(i)))
     deltaConstraint   =+criticalOverdensity_     %value     (time=timeConstraint,mass=massConstraint)             &
          &             /linearGrowth_            %value     (time=timeConstraint                    )
     varianceConstraint=+cosmologicalMassVariance_%rootVariance(time=timePresent ,mass=massConstraint)**2
     conditioned    (i)=mergerTreeBranchingProbabilityConstraintConditioned(deltaConstraint,varianceConstraint,windowLogMass,      stepLogMass,      stepDelta,wrapped,massResolution_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_)
     conditionedFine(i)=mergerTreeBranchingProbabilityConstraintConditioned(deltaConstraint,varianceConstraint,windowLogMass,0.5d0*stepLogMass,0.5d0*stepDelta,wrapped,massResolution_,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_)
     ! One step away from the constraint the expectation is the step multiplied by the rate of producing a progenitor in the
     ! window, as nothing else has yet contributed. This fixes the normalization of the rates used in the tabulation, including
     ! the conversion to the present-day-scaled critical overdensity. Limited by the quadrature across the edges of the window.
     delta=deltaConstraint-stepDelta
     do j=1,countParent
        if (massParents(j) < 2.0d0*massConstraint*exp(windowLogMass)) cycle
        expected=stepDelta*gainWindow(massParents(j))
        actual  =conditioned(i)%progenitorExpectation(massParents(j),delta,massResolutionFixed,node)
        write (label,'(e8.1)') massParents(j)
        call Assert('first step from the constraint, M = '//trim(adjustl(label)),actual,expected,relTol=toleranceFirstStep)
     end do
     ! Away from the constraint, the tabulated solution must satisfy the backward equation, and be converged.
     do k=1,countOffset
        delta=deltaConstraint-deltaOffsets(k)
        time =criticalOverdensity_%timeOfCollapse(criticalOverdensity=delta,mass=massConstraint,node=node)
        write (label,'(f4.2)') deltaOffsets(k)
        call Unit_Tests_Begin_Group("δ₁ - δ = "//trim(adjustl(label)))
        do j=1,countParent
           massParent=massParents(j)
           H         =conditioned(i)%progenitorExpectation(massParent,delta,massResolutionFixed,node)
           derivative=-(                                                                                                   &
                &       +conditioned(i)%progenitorExpectation(massParent,delta+stepDelta,massResolutionFixed,node)         &
                &       -conditioned(i)%progenitorExpectation(massParent,delta-stepDelta,massResolutionFixed,node)         &
                &      )                                                                                                   &
                &     /2.0d0                                                                                               &
                &     /stepDelta
           gradient  =+(                                                                                                   &
                &       +conditioned(i)%progenitorExpectation(massParent*exp(+stepLogMass),delta,massResolutionFixed,node) &
                &       -conditioned(i)%progenitorExpectation(massParent*exp(-stepLogMass),delta,massResolutionFixed,node) &
                &      )                                                                                                   &
                &     /2.0d0                                                                                               &
                &     /stepLogMass
           growthFactor=+cosmologicalMassVariance_%rootVariance(time=time       ,mass=massParent) &
                &       /cosmologicalMassVariance_%rootVariance(time=timePresent,mass=massParent)
           drift     =+wrapped%fractionSubresolution(massParent,delta,time,massResolutionFixed,node) &
                &     *growthFactor
           expected  =-drift*gradient+integralBackward(massParent,delta)
           write (label,'(e8.1)') massParent
           call Assert('backward equation, M = '//trim(adjustl(label)),derivative,expected,relTol=toleranceResidual)
           HFine=conditionedFine(i)%progenitorExpectation(massParent,delta,massResolutionFixed,node)
           call Assert('grid convergence, M = '//trim(adjustl(label)),H,HFine,relTol=toleranceConvergence)
        end do
        call Unit_Tests_End_Group()
     end do
     ! The remaining checks are made well away from the constraint.
     delta=deltaConstraint-deltaOffsets(1)
     time =criticalOverdensity_%timeOfCollapse(criticalOverdensity=delta,mass=massConstraint,node=node)
     ! Halos can not gain mass going back in time, so masses below the window can never contribute.
     massParent=massConstraint*exp(-windowLogMass-3.0d0*stepLogMass)
     call Assert('expectation vanishes below the window',conditioned(i)%progenitorExpectation(massParent,deltaConstraint-deltaOffsets(1),massResolutionFixed,node),0.0d0)
     ! The conditioned rate must be the symmetrized rate multiplied by the ratio of expectations.
     massParent=1.0d13
     do j=1,3
        select case (j)
        case (1)
           massDrawn=2.0d11                               ! In the lower half.
        case (2)
           massDrawn=7.0d12                               ! In the upper half.
        case (3)
           massDrawn=massParent-0.5d0*massResolutionFixed ! Unresolved complement.
        end select
        rateExpected=+wrapped%rate(massParent,delta,time,massDrawn           ,node)
        if (massParent-massDrawn >= massResolutionFixed) rateExpected=rateExpected+wrapped%rate(massParent,delta,time,massParent-massDrawn,node)
        rateExpected=+0.5d0                                                                  &
             &       *rateExpected                                                           &
             &       *conditioned(i)%progenitorExpectation(massDrawn ,delta,massResolutionFixed,node) &
             &       /conditioned(i)%progenitorExpectation(massParent,delta,massResolutionFixed,node)
        call Assert('conditioned rate',conditioned(i)%rate(massParent,delta,time,massDrawn,node),rateExpected,relTol=toleranceExact)
     end do
     ! Quantities unaffected by the conditioning are forwarded unchanged.
     call Assert('sub-resolution fraction forwarded',conditioned(i)%fractionSubresolution(massParent,delta,time,massResolutionFixed,node),wrapped%fractionSubresolution(massParent,delta,time,massResolutionFixed,node),relTol=toleranceExact)
     call Assert('maximum step forwarded'           ,conditioned(i)%stepMaximum          (massParent,delta,time,massResolutionFixed     ),wrapped%stepMaximum          (massParent,delta,time,massResolutionFixed     ),relTol=toleranceExact)
     ! The branch mass returned must invert the cumulative rate. The accuracy is limited by that of the integrals in the class.
     probability=conditioned(i)%probability(massParent,delta,time,massResolutionFixed,node)
     do j=1,countFraction
        massBranch=conditioned(i)%massBranch(massParent,delta,time,massResolutionFixed,fractionsSought(j)*probability,randomNumberGenerator_,node)
        call Assert('branch mass inverts the cumulative rate',integralRate(massParent,delta,time,massBranch),fractionsSought(j)*probability,relTol=toleranceInverse)
     end do
     call Unit_Tests_End_Group()
     call Unit_Tests_End_Group()
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  double precision function variance(mass)
    !!{RST
    The variance at the present day.
    !!}
    implicit none
    double precision, intent(in   ) :: mass

    variance=cosmologicalMassVariance_%rootVariance(time=timePresent,mass=mass)**2
    return
  end function variance

  double precision function rateDraw(massHalo,massDrawn)
    !!{RST
    The analytic rate per unit drawn mass and per unit present-day-scaled critical overdensity for a halo to draw the given mass.
    !!}
    implicit none
    double precision, intent(in   ) :: massHalo, massDrawn

    rateDraw=+(massHalo/massDrawn)                                       &
         &   *(2.0d0/3.0d0)*variance(massDrawn)/massDrawn                &
         &   /sqrt(2.0d0*Pi*(variance(massDrawn)-variance(massHalo))**3)
    return
  end function rateDraw

  double precision function gainWindow(massHalo)
    !!{RST
    The rate of producing a progenitor in the window, by drawing it directly or as the complement of the drawn mass.
    !!}
    implicit none
    double precision, intent(in   ) :: massHalo
    double precision                :: x
    integer                         :: k

    gainWindow=0.0d0
    do k=1,countFine
       x         =massConstraint*exp(windowLogMass*(2.0d0*(dble(k)-0.5d0)/dble(countFine)-1.0d0))
       gainWindow=+gainWindow                                                                     &
            &     +0.5d0                                                                          &
            &     *(rateDraw(massHalo,x)+rateDraw(massHalo,massHalo-x))                           &
            &     *x                                                                              &
            &     *2.0d0*windowLogMass/dble(countFine)
    end do
    return
  end function gainWindow

  double precision function expectation(mass,delta)
    !!{RST
    The tabulated expectation, or zero below the mass resolution.
    !!}
    implicit none
    double precision, intent(in   ) :: mass, delta

    expectation=0.0d0
    if (mass >= massResolutionFixed) expectation=conditioned(i)%progenitorExpectation(mass,delta,massResolutionFixed,node)
    return
  end function expectation

  double precision function integralBackward(massHalo,delta)
    !!{RST
    The jump term of the backward equation: draws of mass in the lower half, integrated over the logarithm of the drawn mass, and in the upper half, integrated over the logarithm of the complement.
    !!}
    implicit none
    double precision, intent(in   ) :: massHalo, delta
    double precision                :: logLower, logUpper, &
         &                             x       , H_
    integer                         :: k

    integralBackward=0.0d0
    H_              =expectation(massHalo,delta)
    if (massHalo > 2.0d0*massResolutionFixed) then
       logLower=log(      massResolutionFixed)
       logUpper=log(0.5d0*massHalo           )
       do k=1,countFine
          x               =exp(logLower+(logUpper-logLower)*(dble(k)-0.5d0)/dble(countFine))
          integralBackward=+integralBackward                                                 &
               &           +0.5d0                                                            &
               &           *rateDraw(massHalo,x)                                             &
               &           *(expectation(x,delta)+expectation(massHalo-x,delta)-H_)          &
               &           *x                                                                &
               &           *(logUpper-logLower)/dble(countFine)
       end do
    end if
    logLower=log(1.0d-10*massHalo)
    logUpper=log(0.5d0  *massHalo)
    do k=1,countFine
       x               =exp(logLower+(logUpper-logLower)*(dble(k)-0.5d0)/dble(countFine))
       integralBackward=+integralBackward                                                 &
            &           +0.5d0                                                            &
            &           *rateDraw(massHalo,massHalo-x)                                    &
            &           *(expectation(massHalo-x,delta)+expectation(x,delta)-H_)          &
            &           *x                                                                &
            &           *(logUpper-logLower)/dble(countFine)
    end do
    return
  end function integralBackward

  double precision function integralRate(massHalo,delta,time,massBranch)
    !!{RST
    The integral of the conditioned rate from the mass resolution to the given branch mass, over the logarithm of the complement.
    !!}
    implicit none
    double precision, intent(in   ) :: massHalo, delta     , &
         &                             time    , massBranch
    double precision                :: logLower, logUpper  , &
         &                             x
    integer                         :: k

    integralRate=0.0d0
    logLower    =log(massHalo-massBranch         )
    logUpper    =log(massHalo-massResolutionFixed)
    do k=1,countFine
       x           =exp(logLower+(logUpper-logLower)*(dble(k)-0.5d0)/dble(countFine))
       integralRate=+integralRate                                             &
            &       +conditioned(i)%rate(massHalo,delta,time,massHalo-x,node) &
            &       *x                                                        &
            &       *(logUpper-logLower)/dble(countFine)
    end do
    return
  end function integralRate

end program Tests_Merger_Tree_Branching_Constraint_Conditioned

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

program Tests_Merger_Tree_Branching_Modifiers
  !!{RST
  Tests of merger tree branching probability rate modifiers.

  Each modifier is a small, closed-form multiplicative factor applied to the branching rate, so each can be checked directly
  against the expression it implements. The modifiers are checked over a grid of child root-variances spanning the full range
  encountered when building a tree, from a progenitor barely above the resolution limit to one at half the parent mass, and at
  more than one parent mass and epoch---the last of these because two of the modifiers memoize the factors which depend only on
  the parent, and a memoization which failed to notice that the parent had changed would otherwise go undetected.
  !!}
  use :: Cosmological_Density_Field          , only : cosmologicalMassVarianceFilteredPower              , criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Cosmology_Functions                 , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters                , only : cosmologyParametersSimple
  use :: Dark_Matter_Particles               , only : darkMatterParticleCDM
  use :: Display                             , only : displayVerbositySet                               , verbosityLevelWorking
  use :: Error                               , only : Error_Handler_Register
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Galacticus_Nodes                    , only : treeNode
  use :: Linear_Growth                       , only : linearGrowthCollisionlessMatter
  use :: Merger_Tree_Branching_Modifiers     , only : mergerTreeBranchingProbabilityModifierIdentity     , mergerTreeBranchingProbabilityModifierParkinson2008        , mergerTreeBranchingProbabilityModifierPCHPlus, mergerTreeBranchingProbabilityModifierMulti, &
       &                                              mergerTreeBranchingProbabilityModifierClass        , multiModifierList
  use :: Power_Spectra_Primordial            , only : powerSpectrumPrimordialPowerLaw
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredSimple
  use :: Power_Spectrum_Window_Functions     , only : powerSpectrumWindowFunctionSharpKSpace
  use :: Transfer_Functions                  , only : transferFunctionIdentity
  use :: Unit_Tests                          , only : Assert                                            , Unit_Tests_Begin_Group                                     , Unit_Tests_End_Group                         , Unit_Tests_Finish
  implicit none
  type            (cosmologyParametersSimple                                   )                :: cosmologyParametersSimple_
  type            (cosmologyFunctionsMatterLambda                              )                :: cosmologyFunctionsMatterLambda_
  type            (cosmologicalMassVarianceFilteredPower                       )                :: cosmologicalMassVarianceFilteredPower_
  type            (powerSpectrumWindowFunctionSharpKSpace                      )                :: powerSpectrumWindowFunctionSharpKSpace_
  type            (powerSpectrumPrimordialPowerLaw                             )                :: powerSpectrumPrimordialPowerLaw_
  type            (transferFunctionIdentity                                    )                :: transferFunctionIdentity_
  type            (powerSpectrumPrimordialTransferredSimple                    )                :: powerSpectrumPrimordialTransferredSimple_
  type            (linearGrowthCollisionlessMatter                             )                :: linearGrowthCollisionlessMatter_
  type            (criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt)                :: criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_
  type            (darkMatterParticleCDM                                       )                :: darkMatterParticleCDM_
  type            (treeNode                                                    )                :: node
  type            (mergerTreeBranchingProbabilityModifierIdentity              )                :: modifierIdentity_
  type            (mergerTreeBranchingProbabilityModifierParkinson2008         )                :: modifierParkinson2008_
  type            (mergerTreeBranchingProbabilityModifierPCHPlus               )                :: modifierPCHPlus_
  type            (mergerTreeBranchingProbabilityModifierMulti                 )                :: modifierMulti_
  ! The list chained by the `multi` modifier must be allocated on the heap: its destructor deallocates each element of the list.
  type            (multiModifierList                                           ), pointer      :: modifierListHead     , modifierListTail
  ! Parameters of the modifiers. G₀, γ₁ and γ₂ take the values favored by Parkinson et al. (2008); γ₃ is given a value which is
  ! neither zero nor unity, so that the term it controls is neither absent nor linear.
  double precision                                                              , parameter     :: G0    =0.57d0, gamma1=0.38d0, &
       &                                                                                           gamma2=-0.01d0, gamma3=0.30d0
  integer                                                                       , parameter     :: countChild=9  , countParent=3
  ! Ratios of the child to the parent root-variance at which the modifiers are evaluated. Unity corresponds to a progenitor of the
  ! parent's own mass and large values to one far below it, so this spans the full range met when building a tree. The value just
  ! above unity probes the limit in which the term controlled by γ₃ vanishes.
  double precision                                                              , dimension(countChild ), parameter :: ratioSigma =[1.000001d0,1.0001d0,1.01d0,1.1d0,1.5d0,2.0d0,5.0d0,2.0d1,1.0d3]
  double precision                                                              , dimension(countParent), parameter :: massParent =[1.0d10,1.0d12,1.0d14]
  double precision                                                              , dimension(countParent), parameter :: redshift   =[0.0d0 ,1.0d0 ,3.0d0 ]
  double precision                                                              , dimension(countChild ) :: modifierActual        , modifierExpected     , &
       &                                                                                                   modifierParkinsonValue, modifierPCHPlusValue , &
       &                                                                                                   sigmaChild
  double precision                                                                               :: sigmaParent           , timeParent           , &
       &                                                                                            expansionFactor       , deltaCritical_
  integer                                                                                        :: i                     , j
  character       (len=32                                                      )                  :: label

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelWorking)
  ! Register error handlers.
  call Error_Handler_Register()
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Build the objects needed. Only the critical overdensity is actually required by the modifiers - the root-variances on which
  ! they depend are passed to them directly - but it in turn requires the rest of the chain.
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
       &                                                                                                                      index_                                  =-1.0d0                                         , &
       &                                                                                                                      running                                 =+0.0d0                                         , &
       &                                                                                                                      runningRunning                          =+0.0d0                                         , &
       &                                                                                                                      wavenumberReference                     =+1.0d0                                         , &
       &                                                                                                                      runningSmallScalesOnly                  =.false.                                          &
       &                                                                                                                     )
  transferFunctionIdentity_                                    =transferFunctionIdentity                                     (                                                                                          &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      time                                    =13.8d0                                           &
       &                                                                                                                     )
  powerSpectrumPrimordialTransferredSimple_                    =powerSpectrumPrimordialTransferredSimple                     (                                                                                          &
       &                                                                                                                      powerSpectrumPrimordial_                =powerSpectrumPrimordialPowerLaw_               , &
       &                                                                                                                      transferFunction_                       =transferFunctionIdentity_                      , &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_                 &
       &                                                                                                                     )
  powerSpectrumWindowFunctionSharpKSpace_                      =powerSpectrumWindowFunctionSharpKSpace                       (                                                                                          &
       &                                                                                                                      cosmologyParameters_                    =cosmologyParametersSimple_                     , &
       &                                                                                                                      normalization                           =0.0d0                                            &
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
       &                                                                                                                      powerSpectrumWindowFunction_            =powerSpectrumWindowFunctionSharpKSpace_          &
       &                                                                                                                     )
  criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_=criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt (                                                                                          &
       &                                                                                                                      linearGrowth_                           =linearGrowthCollisionlessMatter_               , &
       &                                                                                                                      cosmologyFunctions_                     =cosmologyFunctionsMatterLambda_                , &
       &                                                                                                                      cosmologicalMassVariance_               =cosmologicalMassVarianceFilteredPower_         , &
       &                                                                                                                      darkMatterParticle_                     =darkMatterParticleCDM_                         , &
       &                                                                                                                      tableStore                              =.true.                                           &
       &                                                                                                                     )
  modifierIdentity_     =mergerTreeBranchingProbabilityModifierIdentity     (                                                                                                  )
  modifierParkinson2008_=mergerTreeBranchingProbabilityModifierParkinson2008(G0,gamma1,gamma2       ,criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_)
  modifierPCHPlus_      =mergerTreeBranchingProbabilityModifierPCHPlus      (G0,gamma1,gamma2,gamma3,criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_)
  ! Chain the two non-trivial modifiers together, so that the product taken by the `multi` modifier can be compared against the
  ! product of the two evaluated separately.
  allocate(modifierListHead)
  allocate(modifierListTail)
  allocate(mergerTreeBranchingProbabilityModifierParkinson2008 :: modifierListHead%modifier_)
  allocate(mergerTreeBranchingProbabilityModifierPCHPlus       :: modifierListTail%modifier_)
  select type (modifier_ => modifierListHead%modifier_)
  type is (mergerTreeBranchingProbabilityModifierParkinson2008)
     modifier_=mergerTreeBranchingProbabilityModifierParkinson2008(G0,gamma1,gamma2       ,criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_)
  end select
  select type (modifier_ => modifierListTail%modifier_)
  type is (mergerTreeBranchingProbabilityModifierPCHPlus       )
     modifier_=mergerTreeBranchingProbabilityModifierPCHPlus      (G0,gamma1,gamma2,gamma3,criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_)
  end select
  modifierListHead%next => modifierListTail
  modifierMulti_        =  mergerTreeBranchingProbabilityModifierMulti(modifierListHead)

  call Unit_Tests_Begin_Group("Merger tree branching probability modifiers")
  do i=1,countParent
     expansionFactor=cosmologyFunctionsMatterLambda_                              %expansionFactorFromRedshift(redshift(i)                                        )
     timeParent     =cosmologyFunctionsMatterLambda_                              %cosmicTime                 (expansionFactor                                    )
     sigmaParent    =cosmologicalMassVarianceFilteredPower_                       %rootVariance               (massParent(i),timeParent                           )
     deltaCritical_ =criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt_%value                      (time=timeParent,mass=massParent(i),node=node       )
     sigmaChild     =sigmaParent*ratioSigma
     write (label,'(a,es8.2e2,a,f3.1)') "M = ",massParent(i)," , z = ",redshift(i)
     call Unit_Tests_Begin_Group(trim(label))
     ! The identity modifier leaves the rate untouched.
     do j=1,countChild
        modifierActual(j)=modifierIdentity_%rateModifier(node,massParent(i),sigmaParent,sigmaChild(j),timeParent)
     end do
     call Assert('identity'      ,modifierActual,spread(1.0d0,1,countChild),relTol=0.0d0)
     ! G(δ,σ_c,σ_p) = G₀ (σ_c/σ_p)^γ₁ (δ_p/σ_p)^γ₂.
     do j=1,countChild
        modifierParkinsonValue(j)=modifierParkinson2008_%rateModifier(node,massParent(i),sigmaParent,sigmaChild(j),timeParent)
        modifierExpected      (j)=+G0                                          &
             &                    *(sigmaChild(j)/sigmaParent)**gamma1         &
             &                    *(deltaCritical_/sigmaParent)**gamma2
     end do
     call Assert('Parkinson (2008)',modifierParkinsonValue,modifierExpected,relTol=1.0d-12)
     ! The same, multiplied by [1-(σ_p/σ_c)²]^γ₃.
     do j=1,countChild
        modifierPCHPlusValue(j)=modifierPCHPlus_%rateModifier(node,massParent(i),sigmaParent,sigmaChild(j),timeParent)
        modifierExpected    (j)=+modifierParkinsonValue(j)                                  &
             &                  *(1.0d0-(sigmaParent/sigmaChild(j))**2)**gamma3
     end do
     call Assert('PCH+'          ,modifierPCHPlusValue,modifierExpected,relTol=1.0d-12)
     ! The multi modifier takes the product over those it chains.
     do j=1,countChild
        modifierActual  (j)=modifierMulti_%rateModifier(node,massParent(i),sigmaParent,sigmaChild(j),timeParent)
        modifierExpected(j)=modifierParkinsonValue(j)*modifierPCHPlusValue(j)
     end do
     call Assert('multi'         ,modifierActual,modifierExpected,relTol=1.0d-12)
     call Unit_Tests_End_Group  (              )
  end do
  call Unit_Tests_End_Group  (                 )
  call Unit_Tests_Finish     (                 )
end program Tests_Merger_Tree_Branching_Modifiers

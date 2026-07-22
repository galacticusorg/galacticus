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
Contains a program to test the revised spherical collapse model for decaying dark matter (DDM) of
:cite:t:`montandon_decaying_2026`, as implemented in the
``Decaying_Dark_Matter_Spherical_Collapse`` module. Target values are taken
from an independent (Python) reimplementation of the paper's equations and cross-checked against the
paper's figures (their Figs. 5 and 8).
!!}

program Test_Decaying_Dark_Matter_Spherical_Collapse
  !!{RST
  Tests of the revised spherical collapse model for decaying dark matter.
  !!}
  use :: Display                                , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Unit_Tests                             , only : Assert                                  , Unit_Tests_Begin_Group       , Unit_Tests_End_Group           , &
       &                                                 Unit_Tests_Finish                       , compareLessThan
  use :: Decaying_Dark_Matter_Spherical_Collapse, only : decayingDarkMatterCriticalOverdensityEdS, decayingDarkMatterGammaTilde , decayingDarkMatterEpsilon      , &
       &                                                 decayingDarkMatterJIntegral             , decayingDarkMatterDeltaCLarge, decayingDarkMatterDeltaCSmall  , &
       &                                                 decayingDarkMatterMassScale1            , decayingDarkMatterDeltaCFit  , decayingDarkMatterMassCollapsed
  implicit none
  double precision, parameter :: timeCollapse=13.8d0 ! Gyr, corresponding to z=0.
  double precision            :: deltaCEdS

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Decaying dark matter: revised spherical collapse")

  ! Einstein-de Sitter critical overdensity [Montandon et al. (2026), their eq. 32].
  call Unit_Tests_Begin_Group("Critical overdensity: limits and fitting function")
  deltaCEdS=decayingDarkMatterCriticalOverdensityEdS()
  call Assert("Einstein-de Sitter critical overdensity",deltaCEdS,1.68647020d0,relTol=1.0d-5)

  ! The J integral [their eq. 38] is negative (collapse is delayed) and grows in magnitude with the
  ! dimensionless decay rate. These are constant-independent (analytic + quadrature).
  call Assert("J(Γ̃     ) is negative",decayingDarkMatterJIntegral(0.69d0),+0.00000d0,compareLessThan)
  call Assert("J(Γ̃=0.69)"            ,decayingDarkMatterJIntegral(0.69d0),-7.60137d0,relTol=1.0d-3  )

  ! Small-mass plateau [their eq. 42]: depends only on the lifetime (via Γ̃), and is larger for
  ! shorter lifetimes. Constant-independent.
  call Assert("small-mass plateau, lifetime= 5 Gyr",deltaCSmallExcess( 5.0d0),2.700269d0,relTol=1.0d-3)
  call Assert("small-mass plateau, lifetime=10 Gyr",deltaCSmallExcess(10.0d0),1.449455d0,relTol=1.0d-3)
  call Assert("small-mass plateau, lifetime=20 Gyr",deltaCSmallExcess(20.0d0),0.755020d0,relTol=1.0d-3)

  ! The transition mass scale M₁ [their eq. 44] scales as vₖ³. This ratio is independent of the
  ! normalization.
  call Assert("transition mass scale M₁ ~ vₖ³",massScale1Ratio(625.0d0,100.0d0),244.140625d0,relTol=1.0d-6)

  ! Absolute normalization of M₁ [their eq. 44]. The target values were obtained by digitizing the
  ! δ_c(M₀) curves of their Fig. 5 (which is drawn at z=1.083, i.e. t_coll=5.558 Gyr for their
  ! cosmology) and fitting their eq. 43 to each curve with M₁ free. This pins the one otherwise
  ! uncertain constant in the model; note that their eq. 45 is only an order-of-magnitude analytic
  ! rationale for the vₖ³ scaling and exceeds the calibrated value by a factor of ~200.
  call Assert("transition mass scale M₁, lifetime= 5 Gyr, vₖ= 100 km/s",massScale1At(5.558d0, 5.0d0, 100.0d0),3.741d9 ,relTol=1.5d-1)
  call Assert("transition mass scale M₁, lifetime= 5 Gyr, vₖ= 707 km/s",massScale1At(5.558d0, 5.0d0, 707.0d0),1.291d12,relTol=1.5d-1)
  call Assert("transition mass scale M₁, lifetime=20 Gyr, vₖ=1880 km/s",massScale1At(5.558d0,20.0d0,1880.0d0),4.072d13,relTol=3.0d-1)

  ! Full fitting function [their eq. 43] recovers the two plateaux and is monotonically decreasing in
  ! mass (delta_c larger at small mass).
  call Assert("δ_c(M₀) decreases with mass",deltaCFitFull(1.0d18,timeCollapse,20.0d0,2250.0d0),deltaCFitFull(1.0d12,timeCollapse,20.0d0,2250.0d0),compareLessThan)
  call Unit_Tests_End_Group()

  ! Mass mapping M_coll/M₀ [their eq. 46]. These values match the theory curves of their Fig. 8 (four
  ! N-body models, z=0). They depend only weakly on the physical constants (through the turnaround
  ! radius), so a modest tolerance is used.
  call Unit_Tests_Begin_Group("Mass mapping: collapsed vs. Lagrangian mass")
  call Assert("Γ⁻¹=10 Gyr, vₖ=1250 km/s, M₀=1e14",massCollapsedFraction(1.0d14,10.0d0,1250.0d0),0.299075d0,relTol=5.0d-3)
  call Assert("Γ⁻¹=10 Gyr, vₖ=1250 km/s, M₀=1e15",massCollapsedFraction(1.0d15,10.0d0,1250.0d0),0.838012d0,relTol=5.0d-3)
  call Assert("Γ⁻¹= 5 Gyr, vₖ= 625 km/s, M₀=1e14",massCollapsedFraction(1.0d14, 5.0d0, 625.0d0),0.685252d0,relTol=5.0d-3)
  call Assert("Γ⁻¹= 5 Gyr, vₖ= 625 km/s, M₀=1e15",massCollapsedFraction(1.0d15, 5.0d0, 625.0d0),0.990827d0,relTol=5.0d-3)
  call Assert("Γ⁻¹=20 Gyr, vₖ= 625 km/s, M₀=1e14",massCollapsedFraction(1.0d14,20.0d0, 625.0d0),0.824764d0,relTol=5.0d-3)
  call Assert("Γ⁻¹=20 Gyr, vₖ=2250 km/s, M₀=1e14",massCollapsedFraction(1.0d14,20.0d0,2250.0d0),0.506446d0,relTol=5.0d-3)
  call Unit_Tests_End_Group()

  ! Limits recovering LCDM: as vₖ → 0 (transition mass → 0, all masses in the large-mass regime) or
  ! as the lifetime → infinity (Γ̃ → 0), the critical overdensity returns to its EdS value and
  ! the collapsed mass returns to the Lagrangian mass. These are structural and constant-independent.
  call Unit_Tests_Begin_Group("Cold dark matter limits")
  call Assert("vₖ → 0: δ_c → EdS"        ,deltaCFitFull        (1.0d14,timeCollapse,10.0d0,   1.0d-1)-deltaCEdS,1.0d-4,compareLessThan)
  call Assert("vₖ → 0: M_coll → M₀"      ,massCollapsedFraction(1.0d14,             10.0d0,   1.0d-1)          ,1.0d+0,relTol=1.0d-4  )
  call Assert("lifetime → ∞: δ_c → EdS"  ,deltaCFitFull        (1.0d14,timeCollapse, 1.0d6,1250.0d+0)-deltaCEdS,1.0d-3,compareLessThan)
  call Assert("lifetime → ∞: M_coll → M₀",massCollapsedFraction(1.0d14,              1.0d6,1250.0d+0)          ,1.0d+0,relTol=1.0d-3  )
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish  ()

contains

  double precision function deltaCSmallExcess(lifetime)
    !!{RST
    Return the small-mass plateau excess :math:`\delta_\mathrm{c}^\mathrm{small}-\delta_\mathrm{c}^\mathrm{EdS}` for the given lifetime (Gyr).
    !!}
    implicit none
    double precision, intent(in   ) :: lifetime

    deltaCSmallExcess=+decayingDarkMatterDeltaCSmall           (decayingDarkMatterGammaTilde(timeCollapse,lifetime)) &
         &            -decayingDarkMatterCriticalOverdensityEdS(                                                   )
    return
  end function deltaCSmallExcess

  double precision function massScale1Ratio(velocityKick1,velocityKick2)
    !!{RST
    Return the ratio of transition mass scales :math:`M_1(v_{k,1})/M_1(v_{k,2})` (which should equal :math:`(v_{k,1}/v_{k,2})^3`).
    !!}
    implicit none
    double precision, intent(in   ) :: velocityKick1, velocityKick2
    double precision                :: gammaTilde

    gammaTilde     =+decayingDarkMatterGammaTilde(timeCollapse ,10.0d0                 )
    massScale1Ratio=+decayingDarkMatterMassScale1(velocityKick1,gammaTilde,timeCollapse) &
         &          /decayingDarkMatterMassScale1(velocityKick2,gammaTilde,timeCollapse)
    return
  end function massScale1Ratio

  double precision function massScale1At(timeCollapse_,lifetime,velocityKick)
    !!{RST
    Return the transition mass scale :math:`M_1` (in :math:`\mathrm{M}_\odot`) at the given collapse
    time (Gyr), lifetime (Gyr), and velocity kick (km/s).
    !!}
    implicit none
    double precision, intent(in   ) :: timeCollapse_, lifetime, velocityKick

    massScale1At=decayingDarkMatterMassScale1(velocityKick,decayingDarkMatterGammaTilde(timeCollapse_,lifetime),timeCollapse_)
    return
  end function massScale1At

  double precision function deltaCFitFull(mass0,timeCollapse_,lifetime,velocityKick)
    !!{RST
    Return the full mass-dependent critical overdensity :math:`\delta_\mathrm{c}^\mathrm{fit}(M_0)`, assembled
    from the module's individual ingredients.
    !!}
    implicit none
    double precision, intent(in   ) :: mass0      , timeCollapse_, lifetime, velocityKick
    double precision                :: gammaTilde , epsilon      , jIntegral, deltaCLarge, &
         &                             deltaCSmall, massScale1

    gammaTilde   =decayingDarkMatterGammaTilde (timeCollapse_,lifetime                            )
    epsilon      =decayingDarkMatterEpsilon    (velocityKick                                      )
    jIntegral    =decayingDarkMatterJIntegral  (gammaTilde                                        )
    deltaCLarge  =decayingDarkMatterDeltaCLarge(gammaTilde   ,epsilon    ,jIntegral               )
    deltaCSmall  =decayingDarkMatterDeltaCSmall(gammaTilde                                        )
    massScale1   =decayingDarkMatterMassScale1 (velocityKick ,gammaTilde ,timeCollapse_           )
    deltaCFitFull=decayingDarkMatterDeltaCFit  (mass0        ,deltaCLarge,deltaCSmall  ,massScale1)
    return
  end function deltaCFitFull

  double precision function massCollapsedFraction(mass0,lifetime,velocityKick)
    !!{RST
    Return the ratio of collapsed mass to Lagrangian mass, :math:`M_\mathrm{coll}/M_0`.
    !!}
    implicit none
    double precision, intent(in   ) :: mass0, lifetime, velocityKick

    massCollapsedFraction=+decayingDarkMatterMassCollapsed(mass0,timeCollapse,lifetime,velocityKick) &
         &                /mass0
    return
  end function massCollapsedFraction

end program Test_Decaying_Dark_Matter_Spherical_Collapse

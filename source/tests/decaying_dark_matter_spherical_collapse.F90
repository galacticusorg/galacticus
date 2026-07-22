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
{\normalfont \ttfamily Decaying\_Dark\_Matter\_Spherical\_Collapse} module. Target values are taken
from an independent (Python) reimplementation of the paper's equations and cross-checked against the
paper's figures (their Figs. 5 and 8).
!!}

program Test_Decaying_Dark_Matter_Spherical_Collapse
  !!{RST
  Tests of the revised spherical collapse model for decaying dark matter.
  !!}
  use :: Display                                , only : displayVerbositySet             , verbosityLevelStandard
  use :: Unit_Tests                             , only : Assert                          , Unit_Tests_Begin_Group       , Unit_Tests_End_Group    , Unit_Tests_Finish, &
       &                                                 compareLessThan
  use :: Decaying_Dark_Matter_Spherical_Collapse, only : decayingDarkMatterCriticalOverdensityEdS, decayingDarkMatterGammaTilde , decayingDarkMatterEpsilon      , &
       &                                                 decayingDarkMatterJIntegral     , decayingDarkMatterDeltaCLarge, decayingDarkMatterDeltaCSmall , &
       &                                                 decayingDarkMatterMassScale1    , decayingDarkMatterDeltaCFit  , decayingDarkMatterMassCollapsed
  implicit none
  double precision, parameter :: timeCollapse=13.8d0 ! Gyr, corresponding to z=0.
  double precision            :: deltaCEdS

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Decaying dark matter: revised spherical collapse")

  deltaCEdS=decayingDarkMatterCriticalOverdensityEdS()

  ! Einstein-de Sitter critical overdensity [Montandon et al. (2026), their eq. 32].
  call Unit_Tests_Begin_Group("Critical overdensity: limits and fitting function")
  call Assert("Einstein-de Sitter critical overdensity",deltaCEdS,1.68647020d0,relTol=1.0d-5)

  ! The J integral [their eq. 38] is negative (collapse is delayed) and grows in magnitude with the
  ! dimensionless decay rate. These are constant-independent (analytic + quadrature).
  call Assert("J(gammaTilde) is negative",decayingDarkMatterJIntegral(0.69d0),0.0d0,compareLessThan)
  call Assert("J(gammaTilde=0.69)"       ,decayingDarkMatterJIntegral(0.69d0),-7.60137d0,relTol=1.0d-3)

  ! Small-mass plateau [their eq. 42]: depends only on the lifetime (via gammaTilde), and is larger for
  ! shorter lifetimes. Constant-independent.
  call Assert("small-mass plateau, lifetime= 5 Gyr",deltaCSmallExcess( 5.0d0),2.700269d0,relTol=1.0d-3)
  call Assert("small-mass plateau, lifetime=10 Gyr",deltaCSmallExcess(10.0d0),1.449455d0,relTol=1.0d-3)
  call Assert("small-mass plateau, lifetime=20 Gyr",deltaCSmallExcess(20.0d0),0.755020d0,relTol=1.0d-3)

  ! The transition mass scale M1 [their eqs. 44,45] scales as v_k^3. This ratio is independent of the
  ! physical constants.
  call Assert("transition mass scale M1 ~ v_k^3",massScale1Ratio(625.0d0,100.0d0),244.140625d0,relTol=1.0d-6)

  ! Full fitting function [their eq. 43] recovers the two plateaux and is monotonically decreasing in
  ! mass (delta_c larger at small mass).
  call Assert("delta_c(M0) decreases with mass",deltaCFitFull(1.0d18,timeCollapse,20.0d0,2250.0d0),deltaCFitFull(1.0d12,timeCollapse,20.0d0,2250.0d0),compareLessThan)
  call Unit_Tests_End_Group()

  ! Mass mapping M_coll/M0 [their eq. 46]. These values match the theory curves of their Fig. 8 (four
  ! N-body models, z=0). They depend only weakly on the physical constants (through the turnaround
  ! radius), so a modest tolerance is used.
  call Unit_Tests_Begin_Group("Mass mapping: collapsed vs. Lagrangian mass")
  call Assert("Gamma^-1=10 Gyr, v_k=1250 km/s, M0=1e14",massCollapsedFraction(1.0d14,10.0d0,1250.0d0),0.299075d0,relTol=5.0d-3)
  call Assert("Gamma^-1=10 Gyr, v_k=1250 km/s, M0=1e15",massCollapsedFraction(1.0d15,10.0d0,1250.0d0),0.838012d0,relTol=5.0d-3)
  call Assert("Gamma^-1= 5 Gyr, v_k= 625 km/s, M0=1e14",massCollapsedFraction(1.0d14, 5.0d0, 625.0d0),0.685252d0,relTol=5.0d-3)
  call Assert("Gamma^-1= 5 Gyr, v_k= 625 km/s, M0=1e15",massCollapsedFraction(1.0d15, 5.0d0, 625.0d0),0.990827d0,relTol=5.0d-3)
  call Assert("Gamma^-1=20 Gyr, v_k= 625 km/s, M0=1e14",massCollapsedFraction(1.0d14,20.0d0, 625.0d0),0.824764d0,relTol=5.0d-3)
  call Assert("Gamma^-1=20 Gyr, v_k=2250 km/s, M0=1e14",massCollapsedFraction(1.0d14,20.0d0,2250.0d0),0.506446d0,relTol=5.0d-3)
  call Unit_Tests_End_Group()

  ! Limits recovering LCDM: as v_k -> 0 (transition mass -> 0, all masses in the large-mass regime) or
  ! as the lifetime -> infinity (gammaTilde -> 0), the critical overdensity returns to its EdS value and
  ! the collapsed mass returns to the Lagrangian mass. These are structural and constant-independent.
  call Unit_Tests_Begin_Group("Cold dark matter limits")
  call Assert("v_k -> 0: delta_c -> EdS"      ,deltaCFitFull        (1.0d14,timeCollapse,10.0d0,1.0d-1)-deltaCEdS,1.0d-4,compareLessThan)
  call Assert("v_k -> 0: M_coll -> M0"        ,massCollapsedFraction(1.0d14,             10.0d0,1.0d-1)          ,1.0d0 ,relTol=1.0d-4  )
  call Assert("lifetime -> inf: delta_c -> EdS",deltaCFitFull       (1.0d14,timeCollapse,1.0d6 ,1250.0d0)-deltaCEdS,1.0d-3,compareLessThan)
  call Assert("lifetime -> inf: M_coll -> M0"  ,massCollapsedFraction(1.0d14,            1.0d6 ,1250.0d0)         ,1.0d0 ,relTol=1.0d-3  )
  call Unit_Tests_End_Group()

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish  ()

contains

  double precision function deltaCSmallExcess(lifetime)
    !!{RST
    Return the small-mass plateau excess $\delta_\mathrm{c}^\mathrm{small}-\delta_\mathrm{c}^\mathrm{EdS}$ for the given lifetime (Gyr).
    !!}
    implicit none
    double precision, intent(in   ) :: lifetime

    deltaCSmallExcess=+decayingDarkMatterDeltaCSmall(decayingDarkMatterGammaTilde(timeCollapse,lifetime)) &
         &            -decayingDarkMatterCriticalOverdensityEdS()
    return
  end function deltaCSmallExcess

  double precision function massScale1Ratio(velocityKick1,velocityKick2)
    !!{RST
    Return the ratio of transition mass scales $M_1(v_{k,1})/M_1(v_{k,2})$ (which should equal $(v_{k,1}/v_{k,2})^3$).
    !!}
    implicit none
    double precision, intent(in   ) :: velocityKick1, velocityKick2
    double precision                :: gammaTilde

    gammaTilde     =decayingDarkMatterGammaTilde(timeCollapse,10.0d0)
    massScale1Ratio=+decayingDarkMatterMassScale1(velocityKick1,gammaTilde,timeCollapse) &
         &          /decayingDarkMatterMassScale1(velocityKick2,gammaTilde,timeCollapse)
    return
  end function massScale1Ratio

  double precision function deltaCFitFull(mass0,timeCollapse_,lifetime,velocityKick)
    !!{RST
    Return the full mass-dependent critical overdensity $\delta_\mathrm{c}^\mathrm{fit}(M_0)$, assembled
    from the module's individual ingredients.
    !!}
    implicit none
    double precision, intent(in   ) :: mass0      , timeCollapse_, lifetime, velocityKick
    double precision                :: gammaTilde , epsilon      , jIntegral, deltaCLarge, &
         &                             deltaCSmall, massScale1

    gammaTilde   =decayingDarkMatterGammaTilde (timeCollapse_,lifetime               )
    epsilon      =decayingDarkMatterEpsilon    (velocityKick                          )
    jIntegral    =decayingDarkMatterJIntegral  (gammaTilde                            )
    deltaCLarge  =decayingDarkMatterDeltaCLarge(gammaTilde   ,epsilon      ,jIntegral )
    deltaCSmall  =decayingDarkMatterDeltaCSmall(gammaTilde                            )
    massScale1   =decayingDarkMatterMassScale1 (velocityKick ,gammaTilde   ,timeCollapse_)
    deltaCFitFull=decayingDarkMatterDeltaCFit  (mass0        ,deltaCLarge  ,deltaCSmall,massScale1)
    return
  end function deltaCFitFull

  double precision function massCollapsedFraction(mass0,lifetime,velocityKick)
    !!{RST
    Return the ratio of collapsed mass to Lagrangian mass, $M_\mathrm{coll}/M_0$.
    !!}
    implicit none
    double precision, intent(in   ) :: mass0, lifetime, velocityKick

    massCollapsedFraction=+decayingDarkMatterMassCollapsed(mass0,timeCollapse,lifetime,velocityKick) &
         &                /mass0
    return
  end function massCollapsedFraction

end program Test_Decaying_Dark_Matter_Spherical_Collapse

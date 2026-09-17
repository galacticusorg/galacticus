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
Contains a program which tests the Einasto and :cite:t:`burkert_structure_1995` mass distributions against numerical
integration of their own density profiles.
!!}

program Test_Mass_Distributions_Einasto_Burkert
  !!{RST
  Tests :galacticus-class:`massDistributionEinasto` and :galacticus-class:`massDistributionBurkert` against an independent
  reference (``massDistributionProfilesCheck.py`` in `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_).

  Both classes are defined by their density alone, and evaluate everything else from closed forms - incomplete :math:`\Gamma`
  functions for Einasto, logarithms and arctangents for Burkert. The reference instead integrates the density numerically, so
  the closed forms are checked against something which shares none of their algebra. This is what
  ``tests.dark_matter_profiles.generic`` cannot do: it compares each profile's analytic results against *its own* numerical
  integrals, which establishes that the closed forms match the density as coded, but not that either is right.

  The densities are

  .. math::

     \rho_\mathrm{Einasto}(r) = \rho_{-2} \exp \left( - {2 \over \alpha} \left[ \left( {r \over r_{-2}} \right)^\alpha - 1 \right] \right), \quad \rho_\mathrm{Burkert}(r) = \rho_0 \left( 1 + {r \over r_\mathrm{s}} \right)^{-1} \left( 1 + \left[ {r \over r_\mathrm{s}} \right]^2 \right)^{-1},

  and the potential is taken with its zero point at infinity, as both classes do. The density is compared directly as well as
  through its integrals: the mass and potential comparisons test the closed forms, and would not notice a density which was
  wrong in shape but right at the radius where the normalization is defined.

  The profiles are built dimensionally, with a scale radius and density normalization which are neither of them unity, and the
  scale-free reference values are multiplied by :math:`\rho_0 r_\mathrm{s}^3` for masses and :math:`\mathrm{G} \rho_0
  r_\mathrm{s}^2` for potentials. A comparison made in scale-free units would not detect a wrong power of either.

  The Burkert radii reach to :math:`10^{-6} r_\mathrm{s}`. The enclosed mass there is the difference of terms which cancel to
  leading order, so the class evaluates it from a series expansion below :math:`10^{-4} r_\mathrm{s}`; these are the radii at
  which an error in that expansion shows. The Burkert profile has :math:`\rho \propto r^{-3}` at large radii, so its total mass
  diverges logarithmically and the class returns ``huge(0.0d0)``; there is correspondingly no half-mass radius to check, unlike
  the Einasto case.
  !!}
  use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
  use :: Display                         , only : displayVerbositySet           , verbosityLevelStandard
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Mass_Distributions              , only : massDistributionEinasto       , massDistributionBurkert
  use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
  use :: Unit_Tests                      , only : Assert                        , Unit_Tests_Begin_Group , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (massDistributionEinasto)                              :: massDistributionEinasto_
  type            (massDistributionBurkert)                              :: massDistributionBurkert_
  type            (coordinateSpherical    )                              :: coordinates
  double precision                         , parameter                   :: radiusScale              =2.0d-2, densityNormalization=1.0d15 ! [Mpc], [M☉/Mpc³]
  ! Einasto shape parameters. The range spans the values found for simulated halos.
  integer                                  , parameter                   :: countShapes              =3
  double precision                         , dimension(countShapes)      :: shapeParameters          =[1.20d-1,1.70d-1,3.00d-1]
  integer                                  , parameter                   :: countRadiiEinasto        =5
  double precision                         , dimension(countRadiiEinasto):: radiiEinasto             =[1.00000000000000d-02,1.00000000000000d-01,1.00000000000000d+00,3.00000000000000d+00,1.00000000000000d+01]
  ! Reference densities [ρ₋₂], masses [ρ₀r_s³] and potentials [Gρ₀r_s²], one row per shape parameter.
  double precision                         , dimension(countShapes,countRadiiEinasto) :: densityEinastoReference  =reshape([                                                                                                                       &
       &                                                                                                                   1.18322737126341d+03,5.94189189347290d+02,1.47241738004144d+02,                                                        &
       &                                                                                                                   5.59079776493434d+01,4.51860551474204d+01,2.78106333587981d+01,                                                        &
       &                                                                                                                   1.00000000000000d+00,1.00000000000000d+00,1.00000000000000d+00,                                                        &
       &                                                                                                                   9.54990931919721d-02,8.92959808518340d-02,7.40811274966597d-02,                                                        &
       &                                                                                                                   4.97028062870189d-03,3.56507826169525d-03,1.31347090504436d-03                                                         &
       &                                                                                                                  ],[countShapes,countRadiiEinasto])
  double precision                         , dimension(countShapes,countRadiiEinasto) :: massEinastoReference     =reshape([                                                                                                                       &
       &                                                                                                                   7.76378727443195d-03,3.47093026305449d-03,7.25675904884169d-04,                                                        &
       &                                                                                                                   4.43334377528103d-01,3.22338261211660d-01,1.65036143560088d-01,                                                        &
       &                                                                                                                   1.07429803394901d+01,1.02961974842208d+01,9.46144616737493d+00,                                                        &
       &                                                                                                                   3.43446111240311d+01,3.32733258551800d+01,3.08574252540443d+01,                                                        &
       &                                                                                                                   9.14771030171143d+01,8.04434056933233d+01,5.87645659940589d+01                                                         &
       &                                                                                                                  ],[countShapes,countRadiiEinasto])
  double precision                         , dimension(countShapes,countRadiiEinasto) :: potentialEinastoReference=reshape([                                                                                                                       &
       &                                                                                                                  -6.38583959514343d+01,-5.41261631819720d+01,-4.11315818624678d+01,                                                      &
       &                                                                                                                  -5.87202474703593d+01,-5.09369952107571d+01,-3.98787432338112d+01,                                                      &
       &                                                                                                                  -4.09473650775716d+01,-3.53884276629913d+01,-2.79277826540795d+01,                                                      &
       &                                                                                                                  -2.85059885420603d+01,-2.33064095616353d+01,-1.65607092299740d+01,                                                      &
       &                                                                                                                  -1.58311687982197d+01,-1.14706522675515d+01,-6.54900477888132d+00                                                       &
       &                                                                                                                  ],[countShapes,countRadiiEinasto])
  double precision                         , dimension(countShapes)      :: massTotalEinastoReference     =[3.19710509633007d+02,1.58264937211993d+02,6.88748678780043d+01]
  double precision                         , dimension(countShapes)      :: radiusHalfMassEinastoReference=[2.62400859193370d+01,9.71142629454652d+00,3.45302722317064d+00]
  ! Burkert.
  integer                                  , parameter                   :: countRadiiBurkert        =9
  double precision                         , dimension(countRadiiBurkert):: radiiBurkert             =[1.00000000000000d-06,1.00000000000000d-05,1.00000000000000d-04,1.00000000000000d-03,1.00000000000000d-02,1.00000000000000d-01,1.00000000000000d+00,1.00000000000000d+01,1.00000000000000d+02]
  double precision                         , dimension(countRadiiBurkert):: densityBurkertReference  =[9.99999000000000d-01,9.99990000000000d-01,9.99900000000000d-01,9.99000000000999d-01,9.90000009900000d-01,9.00090009000900d-01,2.50000000000000d-01,9.00090009000900d-04,9.90000009900000d-07]
  double precision                         , dimension(countRadiiBurkert):: massBurkertReference     =[4.18878706319374d-18,4.18875878885986d-15,4.18847604552103d-12,4.18564861213460d-09,4.15737429604537d-06,3.87479476142188d-03,1.59795607036613d+00,2.03218812691716d+01,4.81263345732133d+01]
  double precision                         , dimension(countRadiiBurkert):: potentialBurkertReference=[-9.86960440108391d+00,-9.86960440143550d+00,-9.86960438014645d+00,-9.86960230774145d+00,-9.86939600877637d+00,-9.84970737084679d+00,-8.71034436121441d+00,-3.22601637489810d+00,-6.06298733594240d-01]
  ! The peak of the Burkert rotation curve, in units of r_s and √(Gρ₀r_s²).
  double precision                         , parameter                   :: radiusRotationCurveMaximumBurkertReference  =3.2446257246042634d0, &
       &                                                                    velocityRotationCurveMaximumBurkertReference=1.6442977505324985d0
  double precision                         , dimension(countRadiiEinasto):: massEinasto                        , potentialEinasto     , &
       &                                                                    densityEinasto
  double precision                         , dimension(countRadiiBurkert):: massBurkert                        , potentialBurkert     , &
       &                                                                    densityBurkert
  double precision                         , dimension(countShapes)      :: massTotalEinasto                   , radiusHalfMassEinasto, &
       &                                                                    massAtRadiusHalfMassEinasto
  ! Tolerances. The reference integrals are converged to better than 10⁻¹², and the closed forms are evaluated through GSL's Γ
  ! functions, so the comparison is limited by the quadrature rather than by either implementation.
  double precision                         , parameter                   :: tolerance                =1.0d-9
  ! `radiusEnclosingMass` is not a closed form: it inverts a tabulation of the scale-free mass built at 30 points per octave,
  ! so its accuracy is set by that interpolation and not by the profile. It is checked here only to confirm that the inversion
  ! finds the right radius; the closed form itself is checked to `tolerance` by evaluating the mass at the reference radius.
  double precision                         , parameter                   :: toleranceTabulated       =1.0d-3
  integer                                                                :: i                                  , j

  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call Unit_Tests_Begin_Group("Einasto and Burkert mass distributions")
  ! Einasto.
  call Unit_Tests_Begin_Group("Einasto profile"      )
  do j=1,countShapes
     massDistributionEinasto_=massDistributionEinasto(shapeParameter=shapeParameters(j),scaleLength=radiusScale,densityNormalization=densityNormalization)
     do i=1,countRadiiEinasto
        coordinates        =[radiiEinasto(i)*radiusScale,0.0d0,0.0d0]
        densityEinasto  (i)=massDistributionEinasto_%density             (coordinates                )
        massEinasto     (i)=massDistributionEinasto_%massEnclosedBySphere(radiiEinasto(i)*radiusScale)
        potentialEinasto(i)=massDistributionEinasto_%potential           (coordinates                )
     end do
     massTotalEinasto           (j)=massDistributionEinasto_%massTotal           (                                             )
     radiusHalfMassEinasto      (j)=massDistributionEinasto_%radiusEnclosingMass (0.5d0*massDistributionEinasto_%massTotal()   )
     massAtRadiusHalfMassEinasto(j)=massDistributionEinasto_%massEnclosedBySphere(radiusHalfMassEinastoReference(j)*radiusScale)
     ! The density at the scale radius is the density normalization by construction, since the exponent vanishes there. This
     ! fixes the meaning of both, which the mass and potential comparisons below then rest on.
     coordinates=[radiusScale,0.0d0,0.0d0]
     call Assert('ρ(r₋₂) = ρ₋₂'       ,massDistributionEinasto_%density(coordinates),densityNormalization                                                                             ,relTol=tolerance)
     call Assert('density, ρ(r)'      ,densityEinasto                               ,densityEinastoReference  (j,:)*densityNormalization                                              ,relTol=tolerance)
     call Assert('enclosed mass, M(r)',massEinasto                                  ,massEinastoReference     (j,:)*densityNormalization*radiusScale**3                               ,relTol=tolerance)
     call Assert('potential, Φ(r)'    ,potentialEinasto                             ,potentialEinastoReference(j,:)*densityNormalization*radiusScale**2*gravitationalConstant_internal,relTol=tolerance)
  end do
  call Assert('total mass, M(∞)',massTotalEinasto,massTotalEinastoReference*densityNormalization*radiusScale**3,relTol=tolerance)
  ! The closed form and the total mass together: half of the total must be enclosed by the independently-located half-mass radius.
  call Assert('half the mass lies inside r_½' ,massAtRadiusHalfMassEinasto,0.5d0*massTotalEinasto                    ,relTol=tolerance         )
  call Assert('radius enclosing half the mass',radiusHalfMassEinasto      ,radiusHalfMassEinastoReference*radiusScale,relTol=toleranceTabulated)
  call Unit_Tests_End_Group  (                       )
  ! Burkert.
  call Unit_Tests_Begin_Group("Burkert profile"      )
  massDistributionBurkert_=massDistributionBurkert(scaleLength=radiusScale,densityNormalization=densityNormalization)
  do i=1,countRadiiBurkert
     coordinates        =[radiiBurkert(i)*radiusScale,0.0d0,0.0d0]
     densityBurkert  (i)=massDistributionBurkert_%density             (coordinates                )
     massBurkert     (i)=massDistributionBurkert_%massEnclosedBySphere(radiiBurkert(i)*radiusScale)
     potentialBurkert(i)=massDistributionBurkert_%potential           (coordinates                )
  end do
  ! The central density is the density normalization by construction, as both factors of the profile are unity at the origin.
  coordinates=[0.0d0,0.0d0,0.0d0]
  call Assert('ρ(0) = ρ₀'                    ,massDistributionBurkert_%density(coordinates)                         ,densityNormalization                                                                                              ,relTol=tolerance)
  call Assert('density, ρ(r)'                ,densityBurkert                                                        ,densityBurkertReference  *densityNormalization                                                                    ,relTol=tolerance)
  call Assert('enclosed mass, M(r)'          ,massBurkert                                                           ,massBurkertReference     *densityNormalization*radiusScale**3                                                     ,relTol=tolerance)
  call Assert('potential, Φ(r)'              ,potentialBurkert                                                      ,potentialBurkertReference*densityNormalization*radiusScale**2*gravitationalConstant_internal                      ,relTol=tolerance)
  ! The total mass diverges logarithmically, so the class reports it as infinite.
  call Assert('total mass diverges'          ,massDistributionBurkert_%massTotal                   () >= huge(0.0d0) ,.true.                                                                                                                            )
  call Assert('radius of peak rotation curve',massDistributionBurkert_%radiusRotationCurveMaximum  ()               ,radiusRotationCurveMaximumBurkertReference  *radiusScale                                                          ,relTol=tolerance)
  call Assert('peak of the rotation curve'   ,massDistributionBurkert_%velocityRotationCurveMaximum()               ,velocityRotationCurveMaximumBurkertReference*radiusScale*sqrt(densityNormalization*gravitationalConstant_internal),relTol=tolerance)
  call Unit_Tests_End_Group()
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Mass_Distributions_Einasto_Burkert

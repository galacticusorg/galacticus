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
Contains a program which tests that the cusp-NFW and soliton-NFW mass distributions reduce to the NFW mass distribution in the
appropriate limits.
!!}

program Test_Mass_Distributions_Cusp_Soliton_NFW_Limits
  !!{RST
  Tests that :galacticus-class:`massDistributionCuspNFW` and :galacticus-class:`massDistributionSolitonNFW` reduce to the
  independently-implemented :galacticus-class:`massDistributionNFW` where they should, and checks their enclosed masses against
  an independent reference (``massDistributionProfilesCheck.py`` in `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_, which integrates each density numerically).

  The cusp-NFW density of :cite:t:`delos_cusp-halo_2025` is

  .. math::

     \rho(r) = \rho_\mathrm{s} \left( y^2 + {r \over r_\mathrm{s}} \right)^{1/2} \left( {r \over r_\mathrm{s}} \right)^{-3/2} \left[ 1 + {r \over r_\mathrm{s}} \right]^{-2},

  so it is *exactly* :math:`\sqrt{1 + y^2 r_\mathrm{s}/r}` times the NFW density of the same :math:`\rho_\mathrm{s}` and
  :math:`r_\mathrm{s}`. That factor is asserted directly, which is sharper than checking that the two agree at large radius: it
  holds at every radius, and it reduces to equality for :math:`y=0`. The cusp dominates within :math:`r \sim y^2 r_\mathrm{s}`,
  where the logarithmic slope tends to :math:`-3/2` rather than the NFW :math:`-1`.

  The soliton-NFW distribution of :cite:t:`schive_understanding_2014` is a soliton core inside ``radiusSoliton`` and an NFW
  profile outside it, so beyond that radius its density must equal NFW exactly, and the mass it adds beyond that radius must
  equal the mass NFW adds. The two regimes are joined by construction rather than by matching, so the profile need not be
  continuous there and no continuity is asserted; what is asserted is that each regime is what it claims to be.

  The enclosed mass of the cusp-NFW profile is evaluated through four different branches depending on how :math:`r/r_\mathrm{s}`
  compares with :math:`y^2` and with unity, so the radii below are chosen to exercise each of them.
  !!}
  use :: Coordinates             , only : coordinateSpherical  , assignment(=)
  use :: Display                 , only : displayVerbositySet  , verbosityLevelStandard
  use :: Events_Hooks            , only : eventsHooksInitialize
  use :: Mass_Distributions      , only : massDistributionNFW  , massDistributionCuspNFW, massDistributionSolitonNFW
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert               , Unit_Tests_Begin_Group , Unit_Tests_End_Group      , Unit_Tests_Finish
  implicit none
  type            (massDistributionNFW       )                            :: massDistributionNFW_
  type            (massDistributionCuspNFW   )                            :: massDistributionCuspNFW_          , massDistributionCuspNFWNoCusp_
  type            (massDistributionSolitonNFW)                            :: massDistributionSolitonNFW_
  type            (coordinateSpherical       )                            :: coordinates
  double precision                            , parameter                 :: radiusScale             =2.0d-2   , densityNormalization=1.0d15 ! [Mpc], [M☉/Mpc³]
  ! The soliton core. `radiusSoliton` is four core radii, so the soliton radii below all lie strictly inside it.
  double precision                            , parameter                 :: radiusCore              =2.0d-3   , radiusSoliton       =8.0d-3, &
       &                                                                     densitySolitonCentral   =1.0d17                                    ! [Mpc], [Mpc], [M☉/Mpc³]
  ! Cusp amplitudes. 0.225 is the value found by `tests.prompt_cusps` for a representative halo.
  double precision                            , parameter                 :: yCusp                   =1.00d-1  , yCuspLarge          =2.25d-1
  integer                                     , parameter                 :: countRadiiCusp          =10
  double precision                            , dimension(countRadiiCusp) :: radiiCusp               =[1.00000000000000d-09,1.00000000000000d-07,1.00000000000000d-06,1.00000000000000d-05,1.00000000000000d-04,1.00000000000000d-03,1.00000000000000d-02,1.00000000000000d-01,1.00000000000000d+00,1.00000000000000d+01]
  double precision                            , dimension(countRadiiCusp) :: massCuspReference       =[2.64922361384326d-14,2.64923116729412d-11,8.37782167905250d-10,2.64998636104605d-08,8.40165965958284d-07,2.72404659142722d-05,1.04297949607350d-03,6.05340260429755d-02,2.45789232707557d+00,1.87652750600798d+01]
  double precision                            , dimension(countRadiiCusp) :: massCuspLargeReference  =[5.96075298764774d-14,5.96075577648368d-11,1.88496450026144d-09,5.96103464350899d-08,1.88584587923866d-06,5.98875683325889d-05,1.96909105327401d-03,7.74854246921795d-02,2.57476621868235d+00,1.89858522250268d+01]
  ! The Schive et al. (2014) soliton core, in units of ρ_c r_c³.
  integer                                     , parameter                 :: countRadiiSoliton       =6
  double precision                            , dimension(countRadiiSoliton) :: radiiSoliton         =[1.00000000000000d-02,1.00000000000000d-01,3.00000000000000d-01,1.00000000000000d+00,2.00000000000000d+00,3.00000000000000d+00]
  double precision                            , dimension(countRadiiSoliton) :: massSolitonReference =[4.18860724378189d-06,4.17054696082707d-03,1.08765856026096d-01,2.78981859039829d+00,8.77699680095753d+00,1.10658801542291d+01]
  ! Radii outside the soliton, in units of the NFW scale radius.
  integer                                     , parameter                 :: countRadiiOuter         =4
  double precision                            , dimension(countRadiiOuter):: radiiOuter              =[1.0d0,3.0d0,1.0d1,1.0d2]
  double precision                            , dimension(countRadiiCusp) :: densityCusp                     , densityCuspNoCusp          , &
       &                                                                     densityNFWCusp                  , densityCuspExpected        , &
       &                                                                     massCusp                        , massCuspLarge              , &
       &                                                                     massCuspNoCusp                  , massNFWCusp
  double precision                            , dimension(countRadiiSoliton) :: massSoliton
  double precision                            , dimension(countRadiiOuter):: densitySolitonNFWOuter          , densityNFWOuter            , &
       &                                                                     massSolitonNFWIncrement         , massNFWIncrement
  ! Tolerances. The limits are exact relations between closed forms, so they are checked at close to machine precision. The
  ! reference integrals are converged to better than 10⁻¹², and the enclosed-mass branches agree with them to better than
  ! 3×10⁻⁹ for these cusp amplitudes.
  double precision                            , parameter                 :: toleranceExact          =1.0d-12, toleranceReference  =1.0d-7
  ! With no cusp the profile is exactly NFW, so its enclosed mass must be too. At the smallest radii the cusp-NFW value comes
  ! from a simplified solution which drops the (1+r/r_s)⁻² factor of the density, costing O(r/r_s) - at most 1.3×10⁻⁶ over the
  ! radii used here.
  double precision                            , parameter                 :: toleranceNoCusp         =1.0d-5
  integer                                                                 :: i

  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call Unit_Tests_Begin_Group("Cusp-NFW and soliton-NFW limits")
  massDistributionNFW_          =massDistributionNFW       (scaleLength=radiusScale,densityNormalization=densityNormalization)
  massDistributionCuspNFW_      =massDistributionCuspNFW   (radiusScale=radiusScale,densityNormalization=densityNormalization,y=yCusp)
  massDistributionCuspNFWNoCusp_=massDistributionCuspNFW   (radiusScale=radiusScale,densityNormalization=densityNormalization,y=0.0d0)
  massDistributionSolitonNFW_   =massDistributionSolitonNFW(radiusScale=radiusScale,densityNormalizationNFW=densityNormalization,radiusCore=radiusCore,radiusSoliton=radiusSoliton,densitySolitonCentral=densitySolitonCentral)
  ! Cusp-NFW.
  call Unit_Tests_Begin_Group("cusp-NFW → NFW")
  do i=1,countRadiiCusp
     coordinates           =[radiiCusp(i)*radiusScale,0.0d0,0.0d0]
     densityCusp        (i)=massDistributionCuspNFW_      %density             (coordinates             )
     densityCuspNoCusp  (i)=massDistributionCuspNFWNoCusp_%density             (coordinates             )
     densityNFWCusp     (i)=massDistributionNFW_          %density             (coordinates             )
     massCusp           (i)=massDistributionCuspNFW_      %massEnclosedBySphere(radiiCusp(i)*radiusScale)
     massCuspNoCusp     (i)=massDistributionCuspNFWNoCusp_%massEnclosedBySphere(radiiCusp(i)*radiusScale)
     massNFWCusp        (i)=massDistributionNFW_          %massEnclosedBySphere(radiiCusp(i)*radiusScale)
     ! The cusp multiplies the NFW density by exactly √(1+y²r_s/r).
     densityCuspExpected(i)=densityNFWCusp(i)*sqrt(1.0d0+yCusp**2/radiiCusp(i))
  end do
  call Assert('y=0 gives the NFW density'                    ,densityCuspNoCusp                                              ,densityNFWCusp                                       ,relTol=toleranceExact    )
  call Assert('y=0 gives the NFW enclosed mass'              ,massCuspNoCusp                                                 ,massNFWCusp                                          ,relTol=toleranceNoCusp   )
  call Assert('the cusp multiplies NFW by √(1+y²r_s/r)'      ,densityCusp                                                    ,densityCuspExpected                                  ,relTol=toleranceExact    )
  call Assert('enclosed mass against numerical integration'  ,massCusp                                                       ,massCuspReference*densityNormalization*radiusScale**3,relTol=toleranceReference)
  call Assert('central slope is -3/2 with a cusp'            ,massDistributionCuspNFW_      %densitySlopeLogarithmicCentral(),-1.5d0                                               ,relTol=toleranceExact    )
  call Assert('central slope is -1 without one'              ,massDistributionCuspNFWNoCusp_%densitySlopeLogarithmicCentral(),-1.0d0                                               ,relTol=toleranceExact    )
  call Unit_Tests_End_Group  (                       )
  ! A second cusp amplitude, which routes several radii through different branches of the enclosed mass.
  call Unit_Tests_Begin_Group("cusp-NFW, larger cusp")
  massDistributionCuspNFW_=massDistributionCuspNFW(radiusScale=radiusScale,densityNormalization=densityNormalization,y=yCuspLarge)
  do i=1,countRadiiCusp
     massCuspLarge(i)=massDistributionCuspNFW_%massEnclosedBySphere(radiiCusp(i)*radiusScale)
  end do
  call Assert('enclosed mass against numerical integration'  ,massCuspLarge,massCuspLargeReference*densityNormalization*radiusScale**3,relTol=toleranceReference)
  call Unit_Tests_End_Group  (                       )
  ! Soliton-NFW.
  call Unit_Tests_Begin_Group("soliton-NFW → NFW")
  do i=1,countRadiiSoliton
     massSoliton(i)=massDistributionSolitonNFW_%massEnclosedBySphere(radiiSoliton(i)*radiusCore)
  end do
  do i=1,countRadiiOuter
     coordinates               =[radiiOuter(i)*radiusScale,0.0d0,0.0d0]
     densitySolitonNFWOuter (i)=massDistributionSolitonNFW_%density             (coordinates               )
     densityNFWOuter        (i)=massDistributionNFW_       %density             (coordinates               )
     ! Beyond the soliton the profile is NFW, so the mass it adds must be the mass NFW adds.
     massSolitonNFWIncrement(i)=massDistributionSolitonNFW_%massEnclosedBySphere(radiiOuter(i)*radiusScale) &
          &                    -massDistributionSolitonNFW_%massEnclosedBySphere(radiusSoliton            )
     massNFWIncrement       (i)=massDistributionNFW_       %massEnclosedBySphere(radiiOuter(i)*radiusScale) &
          &                    -massDistributionNFW_       %massEnclosedBySphere(radiusSoliton            )
  end do
  call Assert('density beyond the soliton is NFW'              ,densitySolitonNFWOuter ,densityNFWOuter                                         ,relTol=toleranceExact    )
  call Assert('mass added beyond the soliton is NFW'           ,massSolitonNFWIncrement,massNFWIncrement                                        ,relTol=toleranceExact    )
  call Assert('soliton core mass against numerical integration',massSoliton            ,massSolitonReference*densitySolitonCentral*radiusCore**3,relTol=toleranceReference)
  call Unit_Tests_End_Group  (                       )
  call Unit_Tests_End_Group  (                       )
  call Unit_Tests_Finish     (                       )
end program Test_Mass_Distributions_Cusp_Soliton_NFW_Limits

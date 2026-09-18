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
Contains a program which tests the response of an NFW dark matter profile to tidal heating.
!!}

!+    Contributions to this file made by: Andrew Benson, Claude.

program Test_Dark_Matter_Profiles_Heated_Tidal_NFW
  !!{RST
  Tests the response of an NFW dark matter profile to :cite:t:`gnedin_tidal_1999` tidal heating, against values computed
  independently by the ``satelliteHeatedProfile.py`` script in the `galacticusDevTools
  &lt;https://github.com/galacticusorg/galacticusDevTools&gt;`_ repository.

  This completes the chain for tidal heating: ``tests.satellite_tidal_heating_rate.exe`` verifies the heating *rate*
  pointwise, ``test-satellite-tidal-heating-evolution.py`` verifies its *accumulation* along an orbit, and this test verifies
  what the accumulated heat then does to the satellite's density profile.

  A shell initially at radius :math:`r_\mathrm{i}`, given specific energy :math:`\epsilon(r_\mathrm{i})`, expands to
  :math:`r_\mathrm{f}` where :math:`\epsilon(r_\mathrm{i}) + \frac{1}{2} \mathrm{G} M(&lt;r_\mathrm{i}) (1/r_\mathrm{f} -
  1/r_\mathrm{i}) = 0`; the enclosed mass is carried with the shell, and the density follows from the Jacobian of the
  mapping. All three are compared here.

  The existing :galacticus-class:`massDistributionSphericalHeated` test, ``tests.dark_matter_profiles.heated.exe``, compares
  against an analytic solution available for an *isothermal* profile, and sets every second-order heating coefficient to
  **zero**. So neither the second-order energy perturbation - which brings in the velocity dispersion and the density
  logarithmic slope - nor any profile without a closed-form solution is covered by it. This test uses an NFW profile and
  non-zero second-order coefficients, and asserts against the first-order-only case as well, so that each term is exercised
  separately.
  !!}
  use :: Display              , only : displayMessage             , displayVerbositySet         , verbosityLevelStandard
  use :: Error                , only : Error_Handler_Register
  use :: Events_Hooks         , only : eventsHooksInitialize
  use :: IO_HDF5              , only : ioHDF5AccessInitialize
  use :: Coordinates          , only : coordinateSpherical        , assignment(=)
  use :: Mass_Distributions   , only : massDistributionClass       , massDistributionNFW        , massDistributionSpherical, massDistributionSphericalHeated, &
       &                               massDistributionHeatingTidal, kinematicsDistributionClass, kinematicsDistributionNFW, nonAnalyticSolversNumerical
  use :: Unit_Tests           , only : Assert                     , Unit_Tests_Begin_Group      , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  ! The satellite's unheated profile: an NFW halo holding the dark matter fraction of a 10¹⁰ M☉ halo, at the concentration and
  ! virial radius used by the companion orbit references.
  double precision, parameter :: massProfile        =8.436409768474469d+09, radiusVirial=6.699153395662194d-02
  double precision, parameter :: radiusScale        =4.466102263774796d-03
  ! The second-order heating coefficients. Non-zero, unlike in tests.dark_matter_profiles.heated.exe. They are small because
  ! the second-order term scales as √Q against the first order's Q, and so dominates at weak heating for any coefficient of
  ! order unity; as set, it runs from 1.5% to 57% of the first-order term across this grid.
  double precision, parameter :: coefficientSecondOrder0  = 0.05d0, coefficientSecondOrder1=0.01d0, &
       &                         coefficientSecondOrder2  = 0.005d0
  double precision, parameter :: correlationVelocityRadius=-0.30d0
  ! Options for the heated mass distribution, matching those of tests.dark_matter_profiles.heated.exe.
  logical         , parameter :: tolerateVelocityMaximumFailure        =.false., tolerateEnclosedMassIntegrationFailure=.false., &
       &                         toleratePotentialIntegrationFailure   =.false.
  double precision, parameter :: fractionRadiusFinalSmall              =1.0d-12, toleranceRelativePotential            =1.0d-3
  ! The reference configurations and values, emitted by `satelliteHeatedProfile.py --fortran` and, for the first-order-only
  ! case, `--fortran --first-order`. Heating rates are in (km/s/Mpc)², radii in Mpc, masses in M☉ and densities in M☉/Mpc³.
  integer         , parameter :: countConfigurations=12
  double precision, dimension(countConfigurations), parameter :: heatingNormalized                =[ 7.153163174399995d+04, 7.153163174399995d+04, 7.153163174399995d+04, 7.153163174399995d+04, 2.145948952319999d+05, 2.145948952319999d+05, 2.145948952319999d+05, 2.145948952319999d+05, 5.722530539519996d+05, 5.722530539519996d+05, 5.722530539519996d+05, 5.722530539519996d+05]
  double precision, dimension(countConfigurations), parameter :: radius                           =[ 6.699153395662194d-03, 2.009746018698658d-02, 4.019492037397316d-02, 6.699153395662194d-02, 6.699153395662194d-03, 2.009746018698658d-02, 4.019492037397316d-02, 6.699153395662194d-02, 6.699153395662194d-03, 2.009746018698658d-02, 4.019492037397316d-02, 6.699153395662194d-02]
  double precision, dimension(countConfigurations), parameter :: radiusInitialReference           =[ 6.627845289853967d-03, 1.873522445661988d-02, 3.173047863465878d-02, 4.123526809303012d-02, 6.522962399709909d-03, 1.699580219313952d-02, 2.543405131902638d-02, 3.034791951262045d-02, 6.299505596428228d-03, 1.445315456409896d-02, 1.933864127063351d-02, 2.180553328458031d-02]
  double precision, dimension(countConfigurations), parameter :: massEnclosedReference            =[ 1.436440594300117d+09, 3.862570861125665d+09, 5.589521336036562d+09, 6.543467974567850d+09, 1.410434815002989d+09, 3.576027326162958d+09, 4.830426025872082d+09, 5.433009386863104d+09, 1.354769497217541d+09, 3.124885522104216d+09, 3.958175672398074d+09, 4.330500545916596d+09]
  double precision, dimension(countConfigurations), parameter :: densityReference                 =[ 4.274092047646993d+14, 2.558908198102040d+13, 2.703848714223727d+12, 4.024672242535334d+11, 4.109221109855131d+14, 2.051154637546063d+13, 1.778742149324999d+12, 2.459994438121197d+11, 3.767937311867849d+14, 1.460522743413092d+13, 1.115459513388445d+12, 1.504667588504663d+11]
  double precision, dimension(countConfigurations), parameter :: radiusInitialFirstOrderReference =[ 6.653662202461202d-03, 1.893251305152432d-02, 3.225239606515413d-02, 4.201038654580513d-02, 6.565999292744095d-03, 1.725226747438226d-02, 2.592088085708914d-02, 3.095543205933896d-02, 6.363999727359859d-03, 1.472217870602326d-02, 1.973802837441915d-02, 2.226095409927675d-02]
  double precision, dimension(countConfigurations), parameter :: massEnclosedFirstOrderReference  =[ 1.442829697886611d+09, 3.894036222333392d+09, 5.647273303943746d+09, 6.613293634549033d+09, 1.421115420420325d+09, 3.619341189930073d+09, 4.893673985055480d+09, 5.502426043065147d+09, 1.370871104481437d+09, 3.174581208652310d+09, 4.020435020806918d+09, 4.396194886761092d+09]
  double precision, dimension(countConfigurations), parameter :: densityFirstOrderReference       =[ 4.302051925108132d+14, 2.594952475607273d+13, 2.743755439821815d+12, 4.064325791608717d+11, 4.154190601967658d+14, 2.087117041769065d+13, 1.801381775719251d+12, 2.478503601209171d+11, 3.830153923018826d+14, 1.487039952666410d+13, 1.128024208562542d+12, 1.515404573109971d+11]
  ! Tolerances. The heated profile's root finder works to a nominal relative tolerance of 10⁻⁶ (`toleranceRelative` in
  ! `massDistributionSphericalHeated`), which bounds all three comparisons; in practice it converges far more tightly and the
  ! largest difference measured over this grid is 1.0 x 10⁻¹⁰. The tolerances below keep four orders of margin over that
  ! while remaining far tighter than any error in the physics: the missing coefficient-derivative term in
  ! `tidalSpecificEnergyGradient`, found by this test, displaced the density by 1.8 x 10⁻³ at the coefficients used here and
  ! by 1.7-3.8 x 10⁻² at those of the reference tidal heating model.
  double precision, parameter :: toleranceRadius=1.0d-7, toleranceMass=1.0d-7, toleranceDensity=1.0d-6
  class           (massDistributionClass       ), pointer      :: massDistributionHeated_
  class           (massDistributionSpherical   ), pointer      :: massDistributionNFW_
  class           (kinematicsDistributionClass ), pointer      :: kinematicsDistributionNFW_
  class           (massDistributionHeatingTidal), pointer      :: massDistributionHeating_
  type            (coordinateSpherical         )               :: coordinates
  character       (len=128                     )               :: message
  integer                                                      :: iConfiguration, iOrder
  double precision                                             :: radiusInitial , massEnclosed        , &
       &                                                          density       , differenceMaximum   , &
       &                                                          differenceRadius, differenceMass    , &
       &                                                          differenceDensity

  call displayVerbositySet   (verbosityLevelStandard)
  call Error_Handler_Register(                      )
  call eventsHooksInitialize (                      )
  call ioHDF5AccessInitialize(                      )
  call Unit_Tests_Begin_Group("Tidally heated NFW profile")
  differenceMaximum=0.0d0
  ! Loop over the two orders: iOrder=1 includes the second-order energy perturbation, iOrder=2 switches it off by zeroing its
  ! coefficients. Comparing both separates the two terms - with only the second-order case asserted, a coefficient which was
  ! ignored entirely would still have to be caught by the first-order comparison.
  do iOrder=1,2
     do iConfiguration=1,countConfigurations
        ! Build the unheated profile and its kinematics afresh for every configuration. The heated distribution memoizes the
        ! last (final radius, initial radius) pair it solved for and returns it unchanged when the same final radius is
        ! requested again, so a single object reused across heating rates would silently return the previous rate's answer.
        allocate(massDistributionNFW :: massDistributionNFW_)
        select type (massDistributionNFW_)
        type is (massDistributionNFW)
           !![
           <referenceConstruct object="massDistributionNFW_" constructor="massDistributionNFW(scaleLength=radiusScale,mass=massProfile,radiusVirial=radiusVirial)"/>
           !!]
        end select
        allocate(kinematicsDistributionNFW :: kinematicsDistributionNFW_)
        select type (kinematicsDistributionNFW_)
        type is (kinematicsDistributionNFW)
           !![
           <referenceConstruct object="kinematicsDistributionNFW_" constructor="kinematicsDistributionNFW(useSeriesApproximation=.false.)"/>
           !!]
        end select
        call massDistributionNFW_%setKinematicsDistribution(kinematicsDistributionNFW_)
        allocate(massDistributionHeatingTidal :: massDistributionHeating_)
        select type (massDistributionHeating_)
        type is (massDistributionHeatingTidal)
           if (iOrder == 1) then
              !![
              <referenceConstruct object="massDistributionHeating_" constructor="massDistributionHeatingTidal(heatingNormalized(iConfiguration),coefficientSecondOrder0,coefficientSecondOrder1,coefficientSecondOrder2,correlationVelocityRadius)"/>
              !!]
           else
              !![
              <referenceConstruct object="massDistributionHeating_" constructor="massDistributionHeatingTidal(heatingNormalized(iConfiguration),0.0d0,0.0d0,0.0d0,correlationVelocityRadius)"/>
              !!]
           end if
        end select
        allocate(massDistributionSphericalHeated :: massDistributionHeated_)
        select type (massDistributionHeated_)
        type is (massDistributionSphericalHeated)
           !![
           <referenceConstruct object="massDistributionHeated_" constructor="massDistributionSphericalHeated(nonAnalyticSolversNumerical,tolerateVelocityMaximumFailure,tolerateEnclosedMassIntegrationFailure,toleratePotentialIntegrationFailure,fractionRadiusFinalSmall,toleranceRelativePotential,massDistributionNFW_,massDistributionHeating_)"/>
           !!]
           radiusInitial=massDistributionHeated_%radiusInitial(radius(iConfiguration))
        end select
        coordinates =[radius(iConfiguration),0.0d0,0.0d0]
        massEnclosed=massDistributionHeated_%massEnclosedBySphere(radius(iConfiguration))
        density     =massDistributionHeated_%density             (coordinates           )
        if (iOrder == 1) then
           differenceRadius =abs(radiusInitial/radiusInitialReference          (iConfiguration)-1.0d0)
           differenceMass   =abs(massEnclosed /massEnclosedReference           (iConfiguration)-1.0d0)
           differenceDensity=abs(density      /densityReference                (iConfiguration)-1.0d0)
           write (message,'(a,i0)') 'initial radius, second order, configuration ',iConfiguration
           call Assert(trim(message),radiusInitial,radiusInitialReference(iConfiguration),relTol=toleranceRadius )
           write (message,'(a,i0)') 'enclosed mass, second order, configuration ',iConfiguration
           call Assert(trim(message),massEnclosed ,massEnclosedReference (iConfiguration),relTol=toleranceMass   )
           write (message,'(a,i0)') 'density, second order, configuration ',iConfiguration
           call Assert(trim(message),density      ,densityReference      (iConfiguration),relTol=toleranceDensity)
        else
           differenceRadius =abs(radiusInitial/radiusInitialFirstOrderReference(iConfiguration)-1.0d0)
           differenceMass   =abs(massEnclosed /massEnclosedFirstOrderReference (iConfiguration)-1.0d0)
           differenceDensity=abs(density      /densityFirstOrderReference      (iConfiguration)-1.0d0)
           write (message,'(a,i0)') 'initial radius, first order only, configuration ',iConfiguration
           call Assert(trim(message),radiusInitial,radiusInitialFirstOrderReference(iConfiguration),relTol=toleranceRadius )
           write (message,'(a,i0)') 'enclosed mass, first order only, configuration ',iConfiguration
           call Assert(trim(message),massEnclosed ,massEnclosedFirstOrderReference (iConfiguration),relTol=toleranceMass   )
           write (message,'(a,i0)') 'density, first order only, configuration ',iConfiguration
           call Assert(trim(message),density      ,densityFirstOrderReference      (iConfiguration),relTol=toleranceDensity)
        end if
        differenceMaximum=max(differenceMaximum,differenceRadius,differenceMass,differenceDensity)
        write (message,'(a,i0,a,i0,a,e12.5,a,e12.5,a,e12.5)') 'order ',iOrder,', configuration ',iConfiguration, &
             & ': fractional difference, r_i ',differenceRadius,', M ',differenceMass,', rho ',differenceDensity
        call displayMessage(trim(message))
        !![
        <objectDestructor name="massDistributionHeated_"   />
        <objectDestructor name="massDistributionHeating_"  />
        <objectDestructor name="kinematicsDistributionNFW_"/>
        <objectDestructor name="massDistributionNFW_"      />
        !!]
     end do
  end do
  write (message,'(a,e12.5)') 'largest fractional difference over all configurations: ',differenceMaximum
  call displayMessage(trim(message))
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Dark_Matter_Profiles_Heated_Tidal_NFW

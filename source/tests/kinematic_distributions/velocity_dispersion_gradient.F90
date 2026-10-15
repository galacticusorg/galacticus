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
Contains a program which tests the logarithmic gradient of the 1D velocity dispersion supplied by kinematics distributions.
!!}

program Test_Kinematic_Distributions_Velocity_Dispersion_Gradient
  !!{RST
  Tests ``velocityDispersion1DGradientLogarithmic()`` of kinematics distributions against a central finite difference of
  ``velocityDispersion1D()`` of the same distribution. This checks that each class's gradient is the derivative of the velocity
  dispersion that class actually returns. Cases covered are:

  * an NFW distribution, which uses the default (isotropic Jeans equation) gradient;
  * spherical scaler, local, and truncated distributions, which override the gradient;
  * an SIDM isothermal distribution inside its interaction radius, where the velocity dispersion is constant, but the Jeans
    equation gradient is not zero;
  * a tidally-heated NFW distribution using the approximate velocity dispersion,
    :math:`\sigma_\mathrm{r}^2 \rightarrow \sigma_\mathrm{r}^2 - (2/3) \epsilon`, which is not a solution of the Jeans equation.
  !!}
  use :: Coordinates                     , only : coordinateSpherical            , assignment(=)
  use :: Dark_Matter_Particles           , only : darkMatterParticleCDM          , darkMatterParticleSelfInteractingDarkMatterConstant
  use :: Display                         , only : displayMessage                 , displayVerbositySet                                , verbosityLevelStandard
  use :: Error                           , only : Error_Handler_Register
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: IO_HDF5                         , only : ioHDF5AccessInitialize
  use :: Mass_Distributions              , only : massDistributionClass          , massDistributionNFW                                , massDistributionSpherical            , massDistributionSphericalHeated        , &
       &                                          massDistributionHeatingTidal   , massDistributionSphericalScaler                    , massDistributionSphericalTruncated   , massDistributionSphericalSIDMIsothermal, &
       &                                          kinematicsDistributionClass    , kinematicsDistributionNFW                          , kinematicsDistributionSphericalScaler, kinematicsDistributionLocal            , &
       &                                          kinematicsDistributionTruncated, kinematicsDistributionSIDMIsothermal               , kinematicsDistributionHeated         , nonAnalyticSolversNumerical
  use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
  use :: Unit_Tests                      , only : Assert                         , Unit_Tests_Begin_Group                             , Unit_Tests_End_Group                 , Unit_Tests_Finish
  implicit none
  ! An NFW halo holding the dark matter fraction of a 10¹⁰ M☉ halo, as in tests.dark_matter_profiles.heated_tidal_NFW.exe.
  double precision                                                                   , parameter :: massProfile                =8.436409768474469d+09                        , radiusVirial             =6.699153395662194d-02, &
       &                                                                                            radiusScale                =4.466102263774796d-03
  ! Tidal heating parameters, matching those of tests.dark_matter_profiles.heated_tidal_NFW.exe.
  double precision                                                                   , parameter :: coefficientSecondOrder0    =0.05d0                                       , coefficientSecondOrder1  =0.01d0               , &
       &                                                                                            coefficientSecondOrder2    =0.005d0                                      , correlationVelocityRadius=-0.30d0
  double precision                                                     , dimension(2), parameter :: heatingNormalized          =[7.153163174399995d+04,5.722530539519996d+05]
  ! Step in log(r) used in the finite difference. The truncation error of the central difference is O(step²)≈10⁻⁷, while the
  ! rounding error, dominated by the root finder used in the heated distribution, is far smaller.
  double precision                                                                   , parameter :: stepLogarithmic            =1.0d-3
  double precision                                                                   , parameter :: toleranceGradient          =1.0d-5
  ! Tolerance for comparison of the velocity dispersion with the numerical solution of the isotropic Jeans equation. That solution is
  ! tabulated with a maximum relative tolerance of 10⁻³, and differs from the analytic NFW solution by up to 1.7×10⁻³ at the radii
  ! used here. Passing unscaled coordinates to the kinematics wrapped by the scaler gives errors of 9-25%.
  double precision                                                                   , parameter :: toleranceVelocityDispersion=3.0d-3
  class           (massDistributionSpherical                          ), pointer                 :: massDistributionNFW_
  class           (massDistributionClass                              ), pointer                 :: massDistribution_
  class           (kinematicsDistributionClass                        ), pointer                 :: kinematicsDistributionNFW_                                               , kinematicsDistribution_
  class           (massDistributionHeatingTidal                       ), pointer                 :: massDistributionHeating_
  type            (darkMatterParticleCDM                              ), pointer                 :: darkMatterParticleCDM_
  type            (darkMatterParticleSelfInteractingDarkMatterConstant), pointer                 :: darkMatterParticleSIDM_
  double precision                                                                               :: radiusInteraction
  integer                                                                                        :: i
  character       (len=32                                             )                          :: label

  call displayVerbositySet   (verbosityLevelStandard)
  call Error_Handler_Register(                      )
  call eventsHooksInitialize (                      )
  call ioHDF5AccessInitialize(                      )
  call Unit_Tests_Begin_Group("Velocity dispersion logarithmic gradient")

  ! NFW distribution, using the default gradient.
  call nfwConstruct()
  call checkGradient('NFW',kinematicsDistributionNFW_,massDistributionNFW_,radiusScale*[0.03d0,0.1d0,0.3d0,1.0d0,3.0d0,10.0d0])

  ! Spherical scaler distribution. A length scaling factor other than unity checks that the scaled coordinates are passed to the
  ! unscaled kinematics distribution.
  allocate(massDistributionSphericalScaler :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionSphericalScaler)
     !![
     <referenceConstruct object="massDistribution_" constructor="massDistributionSphericalScaler(factorScalingLength=2.0d0,factorScalingMass=0.8d0,massDistribution_=massDistributionNFW_)"/>
     !!]
  end select
  allocate(kinematicsDistributionSphericalScaler :: kinematicsDistribution_)
  select type (kinematicsDistribution_)
  type is (kinematicsDistributionSphericalScaler)
     !![
     <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionSphericalScaler(factorScalingLength=2.0d0,factorScalingMass=0.8d0,kinematicsDistribution_=kinematicsDistributionNFW_)"/>
     !!]
  end select
  call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
  call checkGradient('scaler',kinematicsDistribution_,massDistribution_,radiusScale*[0.03d0,0.3d0,3.0d0],checkJeans=.true.)
  !![
  <objectDestructor name="kinematicsDistribution_"/>
  <objectDestructor name="massDistribution_"      />
  !!]

  ! Local distribution.
  allocate(kinematicsDistributionLocal :: kinematicsDistribution_)
  select type (kinematicsDistribution_)
  type is (kinematicsDistributionLocal)
     !![
     <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionLocal(alpha=0.7d0)"/>
     !!]
  end select
  call checkGradient('local',kinematicsDistribution_,massDistributionNFW_,radiusScale*[0.03d0,0.3d0,3.0d0])
  !![
  <objectDestructor name="kinematicsDistribution_"/>
  !!]

  ! Truncated distribution, inside the minimum truncation radius where the piecewise solution is used.
  allocate(massDistributionSphericalTruncated :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionSphericalTruncated)
     !![
     <referenceConstruct object="massDistribution_" constructor="massDistributionSphericalTruncated(radiusTruncateMinimum=3.0d0*radiusScale,radiusTruncateMaximum=6.0d0*radiusScale,nonAnalyticSolver=nonAnalyticSolversNumerical,massDistribution_=massDistributionNFW_)"/>
     !!]
  end select
  allocate(kinematicsDistributionTruncated :: kinematicsDistribution_)
  select type (kinematicsDistribution_)
  type is (kinematicsDistributionTruncated)
     !![
     <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionTruncated()"/>
     !!]
  end select
  call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
  call checkGradient('truncated',kinematicsDistribution_,massDistribution_,radiusScale*[0.1d0,0.5d0,1.0d0,2.0d0])
  !![
  <objectDestructor name="kinematicsDistribution_"/>
  <objectDestructor name="massDistribution_"      />
  !!]

  ! SIDM isothermal distribution, inside the interaction radius.
  allocate(darkMatterParticleCDM_ )
  allocate(darkMatterParticleSIDM_)
  !![
  <referenceConstruct object="darkMatterParticleCDM_"  constructor="darkMatterParticleCDM()"/>
  <referenceConstruct object="darkMatterParticleSIDM_" constructor="darkMatterParticleSelfInteractingDarkMatterConstant(crossSectionSelfInteraction=1.0d0,darkMatterParticle_=darkMatterParticleCDM_)"/>
  !!]
  allocate(massDistributionSphericalSIDMIsothermal :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionSphericalSIDMIsothermal)
     !![
     <referenceConstruct object="massDistribution_" constructor="massDistributionSphericalSIDMIsothermal(timeAge=13.8d0,velocityRelativeMean=sqrt(gravitationalConstant_internal*massProfile/radiusVirial),nonAnalyticSolver=nonAnalyticSolversNumerical,massDistribution_=massDistributionNFW_,darkMatterParticle_=darkMatterParticleSIDM_)"/>
     !!]
     radiusInteraction=massDistribution_%radiusInteraction()
  end select
  allocate(kinematicsDistributionSIDMIsothermal :: kinematicsDistribution_)
  select type (kinematicsDistribution_)
  type is (kinematicsDistributionSIDMIsothermal)
     !![
     <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionSIDMIsothermal()"/>
     !!]
  end select
  call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
  call checkGradient('SIDM isothermal',kinematicsDistribution_,massDistribution_,radiusInteraction*[0.1d0,0.5d0,0.9d0])
  !![
  <objectDestructor name="kinematicsDistribution_"/>
  <objectDestructor name="massDistribution_"      />
  <objectDestructor name="darkMatterParticleSIDM_"/>
  <objectDestructor name="darkMatterParticleCDM_" />
  !!]

  ! Tidally-heated NFW distribution, using the approximate velocity dispersion.
  do i=1,size(heatingNormalized)
     allocate(massDistributionHeatingTidal :: massDistributionHeating_)
     select type (massDistributionHeating_)
     type is (massDistributionHeatingTidal)
        !![
        <referenceConstruct object="massDistributionHeating_" constructor="massDistributionHeatingTidal(heatingNormalized(i),coefficientSecondOrder0,coefficientSecondOrder1,coefficientSecondOrder2,correlationVelocityRadius)"/>
        !!]
     end select
     allocate(massDistributionSphericalHeated :: massDistribution_)
     select type (massDistribution_)
     type is (massDistributionSphericalHeated)
        !![
        <referenceConstruct object="massDistribution_" constructor="massDistributionSphericalHeated(nonAnalyticSolversNumerical,.false.,.false.,.false.,1.0d-12,1.0d-3,massDistributionNFW_,massDistributionHeating_)"/>
        !!]
     end select
     allocate(kinematicsDistributionHeated :: kinematicsDistribution_)
     select type (kinematicsDistribution_)
     type is (kinematicsDistributionHeated)
        !![
        <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionHeated(nonAnalyticSolver=nonAnalyticSolversNumerical,velocityDispersionApproximate=.true.,toleranceRelativeVelocityDispersion=1.0d-6,toleranceRelativeVelocityDispersionMaximum=1.0d-3)"/>
        !!]
     end select
     call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
     write (label,'(a,i1)') 'heated (approximate) ',i
     call checkGradient(trim(label),kinematicsDistribution_,massDistribution_,radiusVirial*[0.1d0,0.3d0,0.6d0,1.0d0])
     !![
     <objectDestructor name="kinematicsDistribution_" />
     <objectDestructor name="massDistribution_"       />
     <objectDestructor name="massDistributionHeating_"/>
     !!]
  end do

  !![
  <objectDestructor name="kinematicsDistributionNFW_"/>
  <objectDestructor name="massDistributionNFW_"      />
  !!]
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  subroutine nfwConstruct()
    !!{RST
    Construct the NFW mass distribution and its kinematics distribution.
    !!}
    implicit none

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
    return
  end subroutine nfwConstruct

  subroutine checkGradient(label,kinematicsDistribution_,massDistribution_,radii,checkJeans)
    !!{RST
    Compare the logarithmic gradient of the velocity dispersion with a central finite difference at each of the given ``radii``.
    The gradient found by assuming the isotropic Jeans equation is reported alongside, for information. If ``checkJeans`` is true
    the distribution is expected to obey the isotropic Jeans equation, so the velocity dispersion is also compared with a
    numerical solution of that equation, and the gradient with the isotropic Jeans gradient. The finite difference check alone
    is self-consistent, and so cannot detect an incorrect velocity dispersion.
    !!}
    implicit none
    character       (len=*                      ), intent(in   )               :: label
    logical                                      , intent(in   ), optional     :: checkJeans
    class           (kinematicsDistributionClass), intent(inout)               :: kinematicsDistribution_
    class           (massDistributionClass      ), intent(inout), target       :: massDistribution_
    double precision                             , intent(in   ), dimension(:) :: radii
    type            (coordinateSpherical        )                              :: coordinates
    double precision                                                           :: velocityDispersionPlus , velocityDispersionMinus, &
         &                                                                        velocityDispersion     , gradient               , &
         &                                                                        gradientDifference     , gradientJeans          , &
         &                                                                        velocityDispersionJeans
    integer                                                                    :: j
    character       (len=256                    )                              :: message

    do j=1,size(radii)
       coordinates            =[radii(j)*exp(+stepLogarithmic),0.0d0,0.0d0]
       velocityDispersionPlus =kinematicsDistribution_%velocityDispersion1D                        (coordinates,massDistribution_,massDistribution_                   )
       coordinates            =[radii(j)*exp(-stepLogarithmic),0.0d0,0.0d0]
       velocityDispersionMinus=kinematicsDistribution_%velocityDispersion1D                        (coordinates,massDistribution_,massDistribution_                   )
       coordinates            =[radii(j)                      ,0.0d0,0.0d0]
       velocityDispersion     =kinematicsDistribution_%velocityDispersion1D                        (coordinates,massDistribution_,massDistribution_                   )
       gradient               =kinematicsDistribution_%velocityDispersion1DGradientLogarithmic     (coordinates,massDistribution_,massDistribution_                   )
       gradientJeans          =kinematicsDistribution_%jeansEquationGradientLogarithmic            (coordinates,massDistribution_,massDistribution_,velocityDispersion)
       gradientDifference     =log(velocityDispersionPlus/velocityDispersionMinus)/2.0d0/stepLogarithmic
       write (message,'(a,a,e10.3,a,f12.8,a,f12.8,a,f12.8,a)') label,': r = ',radii(j),' Mpc; dlogσ/dlogr = ',gradient,' (finite difference ',gradientDifference,'; isotropic Jeans ',gradientJeans,')'
       call displayMessage(trim(message))
       write (message,'(a,a,e10.3,a)') label,': r = ',radii(j),' Mpc'
       call Assert(trim(message),gradient,gradientDifference,absTol=toleranceGradient)
       if (present(checkJeans)) then
          if (checkJeans) then
             velocityDispersionJeans=kinematicsDistribution_%velocityDispersion1DNumerical(coordinates,massDistribution_,massDistribution_)
             write (message,'(a,a,e10.3,a)') label,': r = ',radii(j),' Mpc; σ vs. isotropic Jeans solution'
             call Assert(trim(message),velocityDispersion,velocityDispersionJeans,relTol=toleranceVelocityDispersion)
             write (message,'(a,a,e10.3,a)') label,': r = ',radii(j),' Mpc; dlogσ/dlogr vs. isotropic Jeans gradient'
             call Assert(trim(message),gradient,gradientJeans,absTol=toleranceGradient)
          end if
       end if
    end do
    return
  end subroutine checkGradient

end program Test_Kinematic_Distributions_Velocity_Dispersion_Gradient

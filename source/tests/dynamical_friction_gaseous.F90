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
Tests of the dimensionless drag force for dynamical friction in a gaseous medium. Reference values for the
:cite:t:`kim_dynamical_2007` class are evaluated directly from their fitting formulae (their eqn. 12b), and for the
:cite:t:`chandrasekhar_dynamical_1943` class from the collisionless formula. The Chandrasekhar integral of a collisional
spherical mass distribution is also checked against the drag computed directly from the sound speed of the gas, and that of
an isothermal distribution flagged as collisionless against the collisionless formula, whichever gaseous drag is supplied.
For a collisional (gaseous) cylindrical distribution, with a uniform effective sound speed, the Chandrasekhar integral is
checked against the drag computed directly from the velocity of the perturber relative to the rotating gas, while for a
collisionless cylindrical distribution it is checked to be independent of the gaseous drag supplied.
!!}

program Tests_Dynamical_Friction_Gaseous
  use :: Coordinates                         , only : assignment(=)                            , coordinateCartesian
  use :: Display                             , only : displayVerbositySet                      , verbosityLevelStandard
  use :: Events_Hooks                        , only : eventsHooksInitialize
  use :: Ideal_Gases_Thermodynamics          , only : Ideal_Gas_Sound_Speed
  use :: IO_HDF5                             , only : ioHDF5AccessInitialize
  use :: Mass_Distributions                  , only : kinematicsDistributionIsothermal         , massDistributionIsothermal        , &
  &                                                   kinematicsDistributionEffectiveSoundSpeed, kinematicsDistributionLocal       , &
  &                                                   massDistributionExponentialDisk
  use :: Satellite_Dynamical_Friction_Gaseous, only : dynamicalFrictionGaseousChandrasekhar1943, dynamicalFrictionGaseousKimKim2007
  use :: Unit_Tests                          , only : Assert                                   , Unit_Tests_Begin_Group            , &
  &                                                   Unit_Tests_End_Group                     , Unit_Tests_Finish
  implicit none
  type            (dynamicalFrictionGaseousKimKim2007       )               :: kimKim2007Exact                  , kimKim2007
  type            (dynamicalFrictionGaseousChandrasekhar1943)               :: chandrasekhar1943
  type            (massDistributionIsothermal               )               :: massDistributionHost             , massDistributionPerturber                 , &
       &                                                                       massDistributionHostCollisionless
  type            (kinematicsDistributionIsothermal         )               :: kinematicsDistributionGas        , kinematicsDistributionCollisionless
  type            (kinematicsDistributionEffectiveSoundSpeed)               :: kinematicsDistributionISM
  type            (kinematicsDistributionLocal              )               :: kinematicsDistributionStellar
  type            (massDistributionExponentialDisk          )               :: massDistributionDiskGas          , massDistributionDiskStellar
  type            (coordinateCartesian                      )               :: position                         , velocity                                  , &
       &                                                                       integral
  double precision                                           , dimension(3) :: integralVector                   , integralVectorChandrasekhar1943           , &
       &                                                                       velocityDisk                     , velocityRelative
  double precision                                                          :: forceFactorBelow                 , forceFactorAbove                          , &
       &                                                                       soundSpeed                       , density                                   , &
       &                                                                       speed                            , velocityCircular
  double precision                                           , parameter    :: velocityDispersion       =2.0d+2 , offsetMachNumber                   =1.0d-9, &
       &                                                                       logarithmCoulomb         =2.5d+0 , machNumber                         =1.5d+0, &
       &                                                                       massPerturber            =1.0d+9 , soundSpeedISM                      =1.0d+1, &
       &                                                                       massPerturberCompact     =1.0d+5

  call displayVerbositySet  (verbosityLevelStandard)
  call eventsHooksInitialize(                      )
  call ioHDF5AccessInitialize(                     )
  call Unit_Tests_Begin_Group("Gaseous dynamical friction")
  kimKim2007Exact  =dynamicalFrictionGaseousKimKim2007       (widthTransition=0.00d0)
  kimKim2007       =dynamicalFrictionGaseousKimKim2007       (widthTransition=0.05d0)
  chandrasekhar1943=dynamicalFrictionGaseousChandrasekhar1943(                      )
  ! Kim & Kim (2007) fitting formulae, without interpolation, in each regime.
  call Assert("Kim & Kim (2007): subsonic, M=0.5"              ,kimKim2007Exact%forceFactor(0.5d0,3.0d0),5.0495103738251834d-2,relTol=1.0d-9)
  call Assert("Kim & Kim (2007): transonic, M=2.0, lnΛ=3"      ,kimKim2007Exact%forceFactor(2.0d0,3.0d0),3.6152961537931320d+0,relTol=1.0d-9)
  call Assert("Kim & Kim (2007): supersonic, M=6.0, lnΛ=3"     ,kimKim2007Exact%forceFactor(6.0d0,3.0d0),2.1627524754662977d+0,relTol=1.0d-9)
  ! Outside of the transition intervals, interpolation must have no effect.
  call Assert("Kim & Kim (2007): interpolation leaves M=0.5"   ,kimKim2007%forceFactor(0.5d0,3.0d0),kimKim2007Exact%forceFactor(0.5d0,3.0d0),relTol=1.0d-12)
  call Assert("Kim & Kim (2007): interpolation leaves M=2.0"   ,kimKim2007%forceFactor(2.0d0,3.0d0),kimKim2007Exact%forceFactor(2.0d0,3.0d0),relTol=1.0d-12)
  ! Within the transition intervals, the force is interpolated linearly between the fits at the interval boundaries.
  call Assert("Kim & Kim (2007): interpolated, M=1.00, lnΛ=3"  ,kimKim2007%forceFactor(1.00d0,3.0d0),1.4319659526659578d+0,relTol=1.0d-9)
  call Assert("Kim & Kim (2007): interpolated, M=0.97, lnΛ=3"  ,kimKim2007%forceFactor(0.97d0,3.0d0),9.9461055306357540d-1,relTol=1.0d-9)
  call Assert("Kim & Kim (2007): interpolated, M=4.40, lnΛ=3"  ,kimKim2007%forceFactor(4.40d0,3.0d0),2.1721912802761930d+0,relTol=1.0d-9)
  ! The interpolated force must be continuous at the boundaries of each transition interval.
  forceFactorBelow=kimKim2007%forceFactor(0.95d0-offsetMachNumber,3.0d0)
  forceFactorAbove=kimKim2007%forceFactor(0.95d0+offsetMachNumber,3.0d0)
  call Assert("Kim & Kim (2007): continuous at M=0.95"         ,forceFactorAbove,forceFactorBelow,relTol=1.0d-6)
  forceFactorBelow=kimKim2007%forceFactor(1.05d0-offsetMachNumber,3.0d0)
  forceFactorAbove=kimKim2007%forceFactor(1.05d0+offsetMachNumber,3.0d0)
  call Assert("Kim & Kim (2007): continuous at M=1.05"         ,forceFactorAbove,forceFactorBelow,relTol=1.0d-6)
  forceFactorBelow=kimKim2007%forceFactor(4.35d0-offsetMachNumber,3.0d0)
  forceFactorAbove=kimKim2007%forceFactor(4.35d0+offsetMachNumber,3.0d0)
  call Assert("Kim & Kim (2007): continuous at M=4.35"         ,forceFactorAbove,forceFactorBelow,relTol=1.0d-6)
  forceFactorBelow=kimKim2007%forceFactor(4.45d0-offsetMachNumber,3.0d0)
  forceFactorAbove=kimKim2007%forceFactor(4.45d0+offsetMachNumber,3.0d0)
  call Assert("Kim & Kim (2007): continuous at M=4.45"         ,forceFactorAbove,forceFactorBelow,relTol=1.0d-6)
  ! Without interpolation the fits are discontinuous at M=1 for lnΛ=3.
  forceFactorBelow=kimKim2007Exact%forceFactor(1.0d0-offsetMachNumber,3.0d0)
  forceFactorAbove=kimKim2007Exact%forceFactor(1.0d0+offsetMachNumber,3.0d0)
  call Assert("Kim & Kim (2007): discontinuous at M=1 if exact",abs(forceFactorAbove-forceFactorBelow) > 0.5d0,.true.)
  ! The force is limited to be non-negative where the fits become negative.
  call Assert("Kim & Kim (2007): non-negative, M=0.01"         ,kimKim2007%forceFactor(1.0d-2,2.0d0),0.0d0)
  call Assert("Kim & Kim (2007): non-negative, M=1000, lnΛ=2"  ,kimKim2007%forceFactor(1.0d+3,2.0d0),0.0d0)
  ! Chandrasekhar (1943) formula.
  call Assert("Chandrasekhar (1943): M=1, lnΛ=2"               ,chandrasekhar1943%forceFactor(1.0d0,2.0d0),3.9749608619759813d-1,relTol=1.0d-9)
  call Assert("Chandrasekhar (1943): M=100, lnΛ=2"             ,chandrasekhar1943%forceFactor(1.0d2,2.0d0),2.0000000000000000d+0,relTol=1.0d-9)
  ! Chandrasekhar integral of a collisional (gaseous) spherical distribution. The perturber is much more compact than its own
  ! enclosed mass at the position considered, so that no suppression due to its extent is applied. The integral must then be
  ! -ρ v I(M,lnΛ)/V³, with the Mach number found from the sound speed of the gas.
  kinematicsDistributionGas=kinematicsDistributionIsothermal(velocityDispersion_=velocityDispersion             )
  massDistributionHost     =massDistributionIsothermal      (mass               =1.0d12 ,lengthReference=0.30d0)
  massDistributionPerturber=massDistributionIsothermal      (mass               =1.0d10 ,lengthReference=0.01d0)
  call massDistributionHost%setKinematicsDistribution(kinematicsDistributionGas)
  position  =[0.05d0,0.0d0,0.0d0]
  soundSpeed=Ideal_Gas_Sound_Speed(kinematicsDistributionGas%temperature(position))
  speed     =machNumber*soundSpeed
  velocity  =[0.0d0,speed,0.0d0]
  density   =massDistributionHost%density(position)
  integral  =massDistributionHost%chandrasekharIntegral(massDistributionHost,massDistributionPerturber,massPerturber,position,velocity,logarithmCoulomb,kimKim2007)
  integralVector=integral
  call Assert("Chandrasekhar integral, collisional, Kim & Kim (2007)"   ,integralVector,[0.0d0,-density*kimKim2007       %forceFactor(machNumber,logarithmCoulomb)/speed**2,0.0d0],relTol=1.0d-9)
  integral  =massDistributionHost%chandrasekharIntegral(massDistributionHost,massDistributionPerturber,massPerturber,position,velocity,logarithmCoulomb,chandrasekhar1943)
  integralVector=integral
  call Assert("Chandrasekhar integral, collisional, Chandrasekhar (1943)",integralVector,[0.0d0,-density*chandrasekhar1943%forceFactor(machNumber,logarithmCoulomb)/speed**2,0.0d0],relTol=1.0d-9)
  ! Chandrasekhar integral of the same distribution flagged as collisionless (as used for isothermal dark matter halos). The
  ! gaseous drag must then be ignored, and the integral must be the collisionless result, -ρ v lnΛ [erf(X)-2X exp(-X²)/sqrt(π)]/V³
  ! with X=V/sqrt(2)/sigma - which is the Chandrasekhar (1943) drag factor evaluated at V/sigma.
  kinematicsDistributionCollisionless=kinematicsDistributionIsothermal(velocityDispersion_=velocityDispersion,collisional=.false.)
  massDistributionHostCollisionless  =massDistributionIsothermal      (mass               =1.0d12            ,lengthReference=0.30d0)
  call massDistributionHostCollisionless%setKinematicsDistribution(kinematicsDistributionCollisionless)
  integral  =massDistributionHostCollisionless%chandrasekharIntegral(massDistributionHostCollisionless,massDistributionPerturber,massPerturber,position,velocity,logarithmCoulomb,kimKim2007)
  integralVector=integral
  call Assert("Chandrasekhar integral, collisionless, Kim & Kim (2007)"   ,integralVector,[0.0d0,-density*chandrasekhar1943%forceFactor(speed/velocityDispersion,logarithmCoulomb)/speed**2,0.0d0],relTol=1.0d-9)
  integral  =massDistributionHostCollisionless%chandrasekharIntegral(massDistributionHostCollisionless,massDistributionPerturber,massPerturber,position,velocity,logarithmCoulomb,chandrasekhar1943)
  integralVector=integral
  call Assert("Chandrasekhar integral, collisionless, Chandrasekhar (1943)",integralVector,[0.0d0,-density*chandrasekhar1943%forceFactor(speed/velocityDispersion,logarithmCoulomb)/speed**2,0.0d0],relTol=1.0d-9)
  ! Effective sound speed kinematics, as used for the interstellar medium.
  kinematicsDistributionISM=kinematicsDistributionEffectiveSoundSpeed(soundSpeedEffective=soundSpeedISM)
  call Assert("Effective sound speed: collisional"                     ,kinematicsDistributionISM%isCollisional       (                                                  ),.true.                                 )
  call Assert("Effective sound speed: sound speed"                     ,kinematicsDistributionISM%soundSpeed          (position                                          ),soundSpeedISM            ,relTol=1.0d-12)
  call Assert("Effective sound speed: velocity dispersion"             ,kinematicsDistributionISM%velocityDispersion1D(position,massDistributionHost,massDistributionHost),soundSpeedISM            ,relTol=1.0d-12)
  call Assert("Effective sound speed: isothermal temperature"          ,Ideal_Gas_Sound_Speed(kinematicsDistributionISM%temperature(position))                            ,sqrt(5.0d0/3.0d0)*soundSpeedISM,relTol=1.0d-9)
  ! Chandrasekhar integral of a collisional (gaseous) cylindrical distribution. The perturber moves vertically through the disk
  ! plane with a speed relative to the rotating gas corresponding to a Mach number of 1.5. It is much more compact than its own
  ! enclosed mass within the scale height of the disk, so that no suppression due to its extent is applied. The integral must then
  ! be -ρ v I(M,lnΛ)/V³, with v the velocity of the perturber relative to the gas.
  massDistributionDiskGas=massDistributionExponentialDisk(radiusScale=3.0d-3,scaleHeight=3.0d-4,mass=1.0d10)
  call massDistributionDiskGas%setKinematicsDistribution(kinematicsDistributionISM)
  position        =[3.0d-3,0.0d0,0.0d0]
  velocityCircular=massDistributionDiskGas%rotationCurve(3.0d-3)
  velocityDisk    =[0.0d0,-velocityCircular,0.0d0]
  velocityRelative=[0.0d0,0.0d0,machNumber*soundSpeedISM]
  velocity        =velocityDisk+velocityRelative
  density         =massDistributionDiskGas%density(position)
  integral        =massDistributionDiskGas%chandrasekharIntegral(massDistributionDiskGas,massDistributionPerturber,massPerturberCompact,position,velocity,logarithmCoulomb,kimKim2007)
  integralVector  =integral
  call Assert("Chandrasekhar integral, cylindrical, collisional"        ,integralVector,-density*velocityRelative*kimKim2007%forceFactor(machNumber,logarithmCoulomb)/(machNumber*soundSpeedISM)**3,relTol=1.0d-9)
  ! Chandrasekhar integral of a collisionless cylindrical distribution must not depend on the gaseous drag.
  kinematicsDistributionStellar=kinematicsDistributionLocal(alpha=1.0d0/sqrt(2.0d0))
  massDistributionDiskStellar  =massDistributionExponentialDisk(radiusScale=3.0d-3,scaleHeight=3.0d-4,mass=1.0d10)
  call massDistributionDiskStellar%setKinematicsDistribution(kinematicsDistributionStellar)
  integral                       =massDistributionDiskStellar%chandrasekharIntegral(massDistributionDiskStellar,massDistributionPerturber,massPerturberCompact,position,velocity,logarithmCoulomb,kimKim2007       )
  integralVector                 =integral
  integral                       =massDistributionDiskStellar%chandrasekharIntegral(massDistributionDiskStellar,massDistributionPerturber,massPerturberCompact,position,velocity,logarithmCoulomb,chandrasekhar1943)
  integralVectorChandrasekhar1943=integral
  call Assert("Chandrasekhar integral, cylindrical, collisionless"      ,integralVector,integralVectorChandrasekhar1943,relTol=1.0d-12)
  call Assert("Chandrasekhar integral, cylindrical, collisionless, drag",integralVector(3) < 0.0d0,.true.)
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Dynamical_Friction_Gaseous

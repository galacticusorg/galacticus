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

program Tests_Kepler_Orbits
  !!{
  Tests for orbital parameter conversions.
  !!}
  use :: Display                         , only : displayVerbositySet           , verbosityLevelStandard
  use :: Kepler_Orbits                   , only : keplerOrbit
  use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
  use :: Numerical_Constants_Math        , only : Pi
  use :: Unit_Tests                      , only : Assert                        , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
          &                                       compareEquals
  implicit none
  type            (keplerOrbit) :: orbit
  double precision              :: valueActual, valueExpected, velocityScale

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Orbital parameter conversions")

  ! Compute velocity scale for unit mass and radius.
  velocityScale=sqrt(gravitationalConstant_internal)

  ! Create a circular orbit.
  call orbit%reset()
  call orbit%velocityRadialSet    (0.0d0*velocityScale)
  call orbit%velocityTangentialSet(1.0d0*velocityScale)
  call orbit%radiusSet            (1.0d0              )
  call orbit%massesSet            (0.0d0,1.0d0        )

  ! Check values are computed correctly.
  valueActual  =orbit%energy         ()
  valueExpected=-0.50d0*velocityScale**2 ! Energy of a circular orbit at the virial radius.
  call Assert('Energy of circular orbit'          ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%angularMomentum()
  valueExpected=velocityScale ! Angular momentum of a circular orbit at the virial radius.
  call Assert('Angular momentum of circular orbit',valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%eccentricity   ()
  valueExpected=0.0d0
  call Assert('Eccentricity of circular orbit'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%semiMajorAxis  ()
  valueExpected=1.0d0
  call Assert('Semi-major axis of circular orbit' ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)

  ! Create an elliptical orbit.
  call orbit%reset()
  call orbit%velocityRadialSet    (0.5d0*velocityScale)
  call orbit%velocityTangentialSet(0.5d0*velocityScale)
  call orbit%radiusSet            (1.0d0      )
  call orbit%massesSet            (0.0d0,1.0d0)

  ! Check values are computed correctly.
  valueActual  =orbit%energy         ()
  valueExpected=-0.750d0*velocityScale**2
  call Assert('Energy of elliptical orbit'          ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%angularMomentum()
  valueExpected=0.5d0*velocityScale
  call Assert('Angular momentum of elliptical orbit',valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%eccentricity   ()
  valueExpected=sqrt(5.0d0/8.0d0)
  call Assert('Eccentricity of elliptical orbit'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%semiMajorAxis  ()
  valueExpected=(2.0d0/3.0d0)
  call Assert('Semi-major axis of elliptical orbit' ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)

  ! Create a circular orbit with equal mass satellite and host.
  call orbit%reset()
  call orbit%velocityRadialSet    (0.0d0*velocityScale)
  call orbit%velocityTangentialSet(1.0d0*velocityScale)
  call orbit%radiusSet            (1.0d0              )
  call orbit%massesSet            (1.0d0,1.0d0        )

  ! Check values are computed correctly.
  valueActual  =orbit%energy         ()
  valueExpected=-0.75d0*velocityScale**2 ! Energy of a circular orbit at the virial radius.
  call Assert('Energy of "circular" orbit with equal mass objects'          ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%angularMomentum()
  valueExpected=0.5d0*velocityScale ! Angular momentum of a circular orbit at the virial radius.
  call Assert('Angular momentum of "circular" orbit with equal mass objects',valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%eccentricity   ()
  valueExpected=0.5d0
  call Assert('Eccentricity of "circular" orbit with equal mass objects'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%semiMajorAxis  ()
  valueExpected=2.0d0/3.0d0
  call Assert('Semi-major axis of "circular" orbit with equal mass objects' ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)

  ! Create a circular orbit, specifying eccentricity, radius and periapsis.
  call orbit%reset()
  call orbit%eccentricitySet    (0.0d0      )
  call orbit%radiusPericenterSet(1.0d0      )
  call orbit%radiusSet          (1.0d0      )
  call orbit%massesSet          (0.0d0,1.0d0)

  ! Check values are computed correctly.
  valueActual  =orbit%energy            ()
  valueExpected=-0.50d0*velocityScale**2 ! Energy of a circular orbit at the virial radius.
  call Assert('Energy of circular orbit'             ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%angularMomentum   ()
  valueExpected=velocityScale ! Angular momentum of a circular orbit at the virial radius.
  call Assert('Angular momentum of circular orbit'   ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%velocityRadial    ()
  valueExpected=0.0d0
  call Assert('Radial velocity of circular orbit'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%velocityTangential()
  valueExpected=velocityScale
  call Assert('Tangential velocity of circular orbit',valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =orbit%semiMajorAxis     ()
  valueExpected=1.0d0
  call Assert('Semi-major axis of circular orbit'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)

  ! Create an elliptical orbit, starting from pericenter, propagate to apocenter. The orbit is placed in the x-y plane, with the
  ! pericenter located on the +x axis.
  call orbit%reset                (                            )
  call orbit%massesSet            (0.0d0                 ,1.0d0)
  call orbit%velocityRadialSet    (0.0d0                       )
  call orbit%velocityTangentialSet(2.29541910309191000d-4      )
  call orbit%radiusSet            (0.14285714285714285d+0      )
  call orbit%thetaSet             (0.5d0*Pi                    )
  call orbit%phiSet               (0.0d0                       )
  call orbit%epsilonSet           (0.0d0                       )
  ! Propagate to orbit to apocenter (or very close to it, to avoid rounding errors).
  call orbit%propagate            (1.0d0-1.0d-6)
  ! Check that we are at the expected position.
  call Assert('Propagation: θ',orbit%theta             (),0.5d0               *Pi,absTol=1.0d-2              )
  call Assert('Propagation: φ',orbit%phi               (),1.0d0               *Pi,absTol=1.0d-2              )
  call Assert('Propagation: ε',orbit%epsilon           (),0.0d0                  ,absTol=1.0d-2              )
  call Assert('Propagation: vᵣ',orbit%velocityRadial    (),0.0d0                  ,absTol=1.0d-2*velocityScale)
  call Assert('Propagation: vᵩ',orbit%velocityTangential(),3.279170147274157d-5   ,absTol=1.0d-2*velocityScale)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Tests_Kepler_Orbits

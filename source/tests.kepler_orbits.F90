!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


program Tests_Kepler_Orbits
  !% Tests for orbital parameter conversions.
  use Unit_Tests
  use Input_Parameters
  use ISO_Varying_String
  use Memory_Management
  use Kepler_Orbits_Structure
  use Numerical_Constants_Physical
  implicit none
  type(varying_string) :: parameterFile
  type(keplerOrbit)    :: thisOrbit
  double precision     :: valueExpected,valueActual,velocityScale

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.orbits.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Orbital parameter conversions")

  ! Open the parameter file.
  parameterFile='parameters.xml'
  call Input_Parameters_File_Open(parameterFile)

  ! Compute velocity scale for unit mass and radius.
  velocityScale=dsqrt(gravitationalConstantGalacticus)

  ! Create a circular orbit.
  call thisOrbit%reset()
  call thisOrbit%velocityRadialSet    (0.0d0*velocityScale)
  call thisOrbit%velocityTangentialSet(1.0d0*velocityScale)
  call thisOrbit%radiusSet            (1.0d0              )
  call thisOrbit%massesSet            (0.0d0,1.0d0        )
  
  ! Check values are computed correctly.
  valueActual  =thisOrbit%energy         ()
  valueExpected=-0.50d0*velocityScale**2 ! Energy of a circular orbit at the virial radius.
  call Assert('Energy of circular orbit'          ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%angularMomentum()
  valueExpected=velocityScale ! Angular momentum of a circular orbit at the virial radius.
  call Assert('Angular momentum of circular orbit',valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%eccentricity   ()
  valueExpected=0.0d0
  call Assert('Eccentricity of circular orbit'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%semiMajorAxis  ()
  valueExpected=1.0d0
  call Assert('Semi-major axis of circular orbit' ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  
  ! Create an elliptical orbit.
  call thisOrbit%reset()
  call thisOrbit%velocityRadialSet    (0.5d0*velocityScale)
  call thisOrbit%velocityTangentialSet(0.5d0*velocityScale)
  call thisOrbit%radiusSet            (1.0d0      )
  call thisOrbit%massesSet            (0.0d0,1.0d0)

  ! Check values are computed correctly.
  valueActual  =thisOrbit%energy         ()
  valueExpected=-0.750d0*velocityScale**2
  call Assert('Energy of elliptical orbit'          ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%angularMomentum()
  valueExpected=0.5d0*velocityScale
  call Assert('Angular momentum of elliptical orbit',valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%eccentricity   ()
  valueExpected=dsqrt(5.0d0/8.0d0)
  call Assert('Eccentricity of elliptical orbit'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%semiMajorAxis  ()
  valueExpected=(2.0d0/3.0d0)
  call Assert('Semi-major axis of elliptical orbit' ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  
  ! Create a circular orbit with equal mass satellite and host.
  call thisOrbit%reset()
  call thisOrbit%velocityRadialSet    (0.0d0*velocityScale)
  call thisOrbit%velocityTangentialSet(1.0d0*velocityScale)
  call thisOrbit%radiusSet            (1.0d0              )
  call thisOrbit%massesSet            (1.0d0,1.0d0        )
  
  ! Check values are computed correctly.
  valueActual  =thisOrbit%energy         ()
  valueExpected=-0.75d0*velocityScale**2 ! Energy of a circular orbit at the virial radius.
  call Assert('Energy of "circular" orbit with equal mass objects'          ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%angularMomentum()
  valueExpected=0.5d0*velocityScale ! Angular momentum of a circular orbit at the virial radius.
  call Assert('Angular momentum of "circular" orbit with equal mass objects',valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%eccentricity   ()
  valueExpected=0.5d0
  call Assert('Eccentricity of "circular" orbit with equal mass objects'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%semiMajorAxis  ()
  valueExpected=2.0d0/3.0d0
  call Assert('Semi-major axis of "circular" orbit with equal mass objects' ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)

  ! Create a circular orbit, specifying eccentricity, radius and periapsis.
  call thisOrbit%reset()
  call thisOrbit%eccentricitySet    (0.0d0      )
  call thisOrbit%radiusPericenterSet(1.0d0      )
  call thisOrbit%radiusSet          (1.0d0      )
  call thisOrbit%massesSet          (0.0d0,1.0d0)
  
  ! Check values are computed correctly.
  valueActual  =thisOrbit%energy            ()
  valueExpected=-0.50d0*velocityScale**2 ! Energy of a circular orbit at the virial radius.
  call Assert('Energy of circular orbit'             ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%angularMomentum   ()
  valueExpected=velocityScale ! Angular momentum of a circular orbit at the virial radius.
  call Assert('Angular momentum of circular orbit'   ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%velocityRadial    ()
  valueExpected=0.0d0
  call Assert('Radial velocity of circular orbit'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%velocityTangential()
  valueExpected=velocityScale
  call Assert('Tangential velocity of circular orbit',valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  valueActual  =thisOrbit%semiMajorAxis     ()
  valueExpected=1.0d0
  call Assert('Semi-major axis of circular orbit'    ,valueActual,valueExpected,compare=compareEquals,relTol=1.0d-6)
  
  ! Close the parameter file.
  call Input_Parameters_File_Close  

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Kepler_Orbits

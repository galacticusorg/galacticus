!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a node operator class that propagates satellite halos along their orbits.

  !# <nodeOperator name="nodeOperatorSatelliteOrbit">
  !#  <description>A node operator class that propagates satellite halos along their orbits.</description>
  !# </nodeOperator>
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteOrbit
     !% A node operator class that propagates satellite halos along their orbits.
     private
   contains
     procedure :: differentialEvolution => satelliteOrbitDifferentialEvolution
  end type nodeOperatorSatelliteOrbit
  
  interface nodeOperatorSatelliteOrbit
     !% Constructors for the {\normalfont \ttfamily satelliteOrbit} node operator class.
     module procedure satelliteOrbitConstructorParameters
  end interface nodeOperatorSatelliteOrbit
  
contains

  function satelliteOrbitConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily satelliteOrbit} node operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorSatelliteOrbit)                :: self
    type(inputParameters           ), intent(inout) :: parameters
    
    self=nodeOperatorSatelliteOrbit()
    !# <inputParametersValidate source="parameters"/>
    return
  end function satelliteOrbitConstructorParameters
  
  subroutine satelliteOrbitDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !% Perform deceleration of a satellite due to dark matter self-interactions.
    use :: Galacticus_Nodes                  , only : nodeComponentSatellite
    use :: Galactic_Structure_Accelerations  , only : Galactic_Structure_Acceleration
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Numerical_Constants_Astronomical  , only : gigaYear                        , megaParsec
    use :: Numerical_Constants_Astronomical  , only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes      , only : kilo
    use :: Vectors                           , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorSatelliteOrbit), intent(inout)          :: self
    type            (treeNode                  ), intent(inout)          :: node
    logical                                     , intent(inout)          :: interrupt
    procedure       (interruptTask             ), intent(inout), pointer :: functionInterrupt
    integer                                     , intent(in   )          :: propertyType
    type            (treeNode                  ), pointer                :: nodeHost
    class           (nodeComponentSatellite    )               , pointer :: satellite
    double precision                            , dimension(3)           :: position         , velocity, &
         &                                                                  acceleration
    double precision                                                     :: massEnclosed     , radius
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (.not.node%isSatellite()) return
    satellite => node     %satellite (        )
    nodeHost  => node     %mergesWith(        )
    position  =  satellite%position  (        )
    velocity  =  satellite%velocity  (        )
    radius    =  Vector_Magnitude    (position)
    call satellite%positionRate(            &
         &                      +kilo       &
         &                      *gigaYear   &
         &                      /megaParsec &
         &                      *velocity   &
         &                     )
    if (radius <= 0.0d0) return ! If radius is non-positive, assume no acceleration.
    massEnclosed=Galactic_Structure_Enclosed_Mass(nodeHost,radius  )
    acceleration=Galactic_Structure_Acceleration (nodeHost,position)           
    ! Include a factor (1+m_{sat}/m_{host})=m_{sat}/µ (where µ is the reduced mass) to convert from the two-body problem of
    ! satellite and host orbitting their common center of mass to the equivalent one-body problem (since we're solving for the
    ! motion of the satellite relative to the center of the host which is held fixed).
    call satellite%velocityRate(                            &
         &                      +acceleration               &
         &                      *(                          &
         &                        +1.0d0                    &
         &                        +satellite%boundMass   () &
         &                        /          massEnclosed   &
         &                      )                           &
         &                     )
    return
  end subroutine satelliteOrbitDifferentialEvolution
  

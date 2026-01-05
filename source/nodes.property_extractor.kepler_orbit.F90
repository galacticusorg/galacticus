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

  use :: Kepler_Orbits, only : keplerOrbit, enumerationKeplerOrbitType

  !![
  <nodePropertyExtractor name="nodePropertyExtractorKeplerOrbit" abstract="yes">
   <description>A property extractor class for {\normalfont \ttfamily keplerOrbit} objects.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple), abstract :: nodePropertyExtractorKeplerOrbit
     !!{
     A property extractor for {\normalfont \ttfamily keplerOrbit} objects.
     !!}
     private
     type            (varying_string            ), allocatable, dimension(:) :: properties   , names_, &
          &                                                                     descriptions_
     double precision                            , allocatable, dimension(:) :: unitsInSI_
     type            (enumerationKeplerOrbitType), allocatable, dimension(:) :: propertyIDs
     integer                                                                 :: count_
     type            (varying_string            )                            :: prefix
   contains
     !![
     <methods>
       <method description="Initialize the properties to be extracted." method="initialize" />
       <method description="Extract properties from a {\normalfont \ttfamily keplerOrbit} object." method="extractFromOrbit" />
     </methods>
     !!]
     procedure :: initialize       => keplerOrbitInitialize
     procedure :: elementCount     => keplerOrbitElementCount
     procedure :: extractFromOrbit => keplerOrbitExtractFromOrbit
     procedure :: names            => keplerOrbitNames
     procedure :: descriptions     => keplerOrbitDescriptions
     procedure :: unitsInSI        => keplerOrbitUnitsInSI
  end type nodePropertyExtractorKeplerOrbit

contains

  subroutine keplerOrbitInitialize(self,properties,prefix)
    !!{
    Initializer for the {\normalfont \ttfamily keplerOrbit} output extractor property extractor class.
    !!}
    use :: Error                           , only : Error_Report
    use :: ISO_Varying_String              , only : trim
    use :: Numerical_Constants_Astronomical, only : megaParsec                  , massSolar
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Kepler_Orbits                   , only : keplerOrbitMasses           , keplerOrbitSpecificReducedMass, keplerOrbitRadius          , keplerOrbitTheta            , &
         &                                          keplerOrbitPhi              , keplerOrbitEpsilon            , keplerOrbitRadiusPericenter, keplerOrbitRadiusApocenter  , &
         &                                          keplerOrbitVelocityRadial   , keplerOrbitVelocityTangential , keplerOrbitEnergy          , keplerOrbitAngularMomentum  , &
         &                                          keplerOrbitEccentricity     , keplerOrbitSemiMajorAxis      , keplerOrbitMassHost        , keplerOrbitMassSatellite    , &
         &                                          enumerationKeplerOrbitEncode
    implicit none
    class    (nodePropertyExtractorKeplerOrbit), intent(inout)               :: self
    type     (varying_string                  ), intent(in   ), dimension(:) :: properties
    character(len=*                           ), intent(in   )               :: prefix
    integer                                                                  :: i

    self%count_=size(properties)
    allocate(self%properties   (self%count_))
    allocate(self%names_       (self%count_))
    allocate(self%descriptions_(self%count_))
    allocate(self%unitsInSI_   (self%count_))
    allocate(self%propertyIDs  (self%count_))
    self%properties=properties
    self%prefix    =prefix
    do i=1,self%count_
       self%propertyIDs(i)=enumerationKeplerOrbitEncode(char(properties(i)),includesPrefix=.false.)
       select case (self%propertyIDs(i)%ID)
       case (keplerOrbitMasses             %ID)
          call Error_Report('"masses" property is unsupported'//{introspection:location})
       case (keplerOrbitMassHost           %ID)
          self%names_       (i)=prefix//'HostMass'
          self%descriptions_(i)='The mass of the host system [M☉].'
          self%unitsInSI_   (i)=massSolar
       case (keplerOrbitMassSatellite      %ID)
          self%names_       (i)=prefix//'SatelliteMass'
          self%descriptions_(i)='The mass of the satellite system [M☉].'
          self%unitsInSI_   (i)=massSolar
       case (keplerOrbitSpecificReducedMass%ID)
          self%names_       (i)=prefix//'SpecificReducedMass'
          self%descriptions_(i)='The reduced mass divided by the mass of the orbiting object.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitRadius             %ID)
          self%names_       (i)=prefix//'Radius'
          self%descriptions_(i)='The current orbital radius [Mpc].'
          self%unitsInSI_   (i)=megaParsec
       case (keplerOrbitTheta              %ID)
          self%names_       (i)=prefix//'Theta'
          self%descriptions_(i)='Orbital angular position, θ, in spherical coordinates.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitPhi                %ID)
          self%names_       (i)=prefix//'Phi'
          self%descriptions_(i)='Orbital angular position, ɸ, in spherical coordinates.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitEpsilon            %ID)
          self%names_       (i)=prefix//'Epsilon'
          self%descriptions_(i)='Direction of the tangential component of velocity.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitRadiusPericenter   %ID)
          self%names_       (i)=prefix//'RadiusPericenter'
          self%descriptions_(i)='Radius of the orbital pericenter [Mpc].'
          self%unitsInSI_   (i)=megaParsec
       case (keplerOrbitRadiusApocenter    %ID)
          self%names_       (i)=prefix//'RadiusApocenter'
          self%descriptions_(i)='Radius of the orbital apocenter [Mpc].'
          self%unitsInSI_   (i)=megaParsec
       case (keplerOrbitVelocityRadial     %ID)
          self%names_       (i)=prefix//'VelocityRadial'
          self%descriptions_(i)='Radial velocity of the orbit [km/s].'
          self%unitsInSI_   (i)=kilo
       case (keplerOrbitVelocityTangential %ID)
          self%names_       (i)=prefix//'VelocityTangential'
          self%descriptions_(i)='Tangential velocity of the orbit [km/s].'
          self%unitsInSI_   (i)=kilo
       case (keplerOrbitEnergy             %ID)
          self%names_       (i)=prefix//'Energy'
          self%descriptions_(i)='Energy of the orbit [M☉ (km/s)²].'
          self%unitsInSI_   (i)=massSolar*kilo**2
       case (keplerOrbitAngularMomentum    %ID)
          self%names_       (i)=prefix//'AngularMomentum'
          self%descriptions_(i)='Angular momentum of the orbit [M☉ Mpc km/s].'
          self%unitsInSI_   (i)=massSolar*kilo*megaParsec
       case (keplerOrbitEccentricity       %ID)
          self%names_       (i)=prefix//'Eccentricity'
          self%descriptions_(i)='Eccentricity of the orbit.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitSemiMajorAxis      %ID)
          self%names_       (i)=prefix//'SemiMajorAxis'
          self%descriptions_(i)='Semi-major axis of the orbit [Mpc].'
          self%unitsInSI_   (i)=megaParsec
       case default
          call Error_Report('unexpected property "'//trim(properties(i))//'"'//{introspection:location})
       end select
    end do
    return
  end subroutine keplerOrbitInitialize

  integer function keplerOrbitElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily keplerOrbit} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorKeplerOrbit), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    keplerOrbitElementCount=self%count_
    return
  end function keplerOrbitElementCount

  function keplerOrbitExtractFromOrbit(self,orbit)
    !!{
    Extract properties from a {\normalfont \ttfamily keplerOrbit} object.
    !!}
    use :: Kepler_Orbits, only : keplerOrbitMassHost      , keplerOrbitSpecificReducedMass, keplerOrbitRadius          , keplerOrbitTheta          , &
         &                       keplerOrbitPhi           , keplerOrbitEpsilon            , keplerOrbitRadiusPericenter, keplerOrbitRadiusApocenter, &
         &                       keplerOrbitVelocityRadial, keplerOrbitVelocityTangential , keplerOrbitEnergy          , keplerOrbitAngularMomentum, &
         &                       keplerOrbitEccentricity  , keplerOrbitSemiMajorAxis      , keplerOrbitMassSatellite
    implicit none
    double precision                                  , dimension(:) , allocatable :: keplerOrbitExtractFromOrbit
    class           (nodePropertyExtractorKeplerOrbit), intent(inout), target      :: self
    type            (keplerOrbit                     ), intent(inout)              :: orbit
    integer                                                                        :: i

    allocate(keplerOrbitExtractFromOrbit(self%count_))
    if (orbit%isDefined()) then
       ! Orbit is defined - extract required properties.
       do i=1,self%count_
          select case (self%propertyIDs(i)%ID)
          case (keplerOrbitMassHost           %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%massHost           ()
          case (keplerOrbitMassSatellite      %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%massSatellite      ()
          case (keplerOrbitSpecificReducedMass%ID)
             keplerOrbitExtractFromOrbit(i)=orbit%specificReducedMass()
          case (keplerOrbitRadius             %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%radius             ()
          case (keplerOrbitTheta              %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%theta              ()
          case (keplerOrbitPhi                %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%phi                ()
          case (keplerOrbitEpsilon            %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%epsilon            ()
          case (keplerOrbitRadiusPericenter   %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%radiusPericenter   ()
          case (keplerOrbitRadiusApocenter    %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%radiusApocenter    ()
          case (keplerOrbitVelocityRadial     %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%velocityRadial     ()
          case (keplerOrbitVelocityTangential %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%velocityTangential ()
          case (keplerOrbitEnergy             %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%energy             ()
          case (keplerOrbitAngularMomentum    %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%angularMomentum    ()
          case (keplerOrbitEccentricity       %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%eccentricity       ()
          case (keplerOrbitSemiMajorAxis      %ID)
             keplerOrbitExtractFromOrbit(i)=orbit%semiMajorAxis      ()
          end select
       end do
    else
       ! Orbit is undefined - set all values to zero.
       keplerOrbitExtractFromOrbit=0.0d0
    end if
    return
  end function keplerOrbitExtractFromOrbit

  subroutine keplerOrbitNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily keplerOrbit} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorKeplerOrbit), intent(inout)                             :: self
    double precision                                  , intent(in   )                             :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%count_))
    names=self%names_
    return
  end subroutine keplerOrbitNames

  subroutine keplerOrbitDescriptions(self,time,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily virialProperies} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorKeplerOrbit), intent(inout)                             :: self
    double precision                                  , intent(in   )                             :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%count_))
    descriptions=self%descriptions_
    return
  end subroutine keplerOrbitDescriptions

  function keplerOrbitUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily virialProperies} properties in the SI system.
    !!}
    implicit none
    double precision                                  , dimension(:) , allocatable :: keplerOrbitUnitsInSI
    class           (nodePropertyExtractorKeplerOrbit), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
   !$GLC attributes unused :: time

    allocate(keplerOrbitUnitsInSI(self%count_))
    keplerOrbitUnitsInSI=self%unitsInSI_
    return
  end function keplerOrbitUnitsInSI


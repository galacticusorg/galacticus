!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  use :: Kepler_Orbits, only : keplerOrbit

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
     type            (varying_string), allocatable, dimension(:) :: properties   , names_, &
          &                                                         descriptions_
     double precision                , allocatable, dimension(:) :: unitsInSI_
     integer                         , allocatable, dimension(:) :: propertyIDs
     integer                                                     :: count_
     type            (varying_string)                            :: prefix
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
     procedure :: type             => keplerOrbitType
  end type nodePropertyExtractorKeplerOrbit

contains

  subroutine keplerOrbitInitialize(self,properties,prefix)
    !!{
    Initializor for the {\normalfont \ttfamily keplerOrbit} output extractor property extractor class.
    !!}
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: ISO_Varying_String              , only : trim
    use :: Numerical_Constants_Astronomical, only : megaParsec               , massSolar
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Kepler_Orbits                   , only : keplerOrbitMasses        , keplerOrbitSpecificReducedMass, keplerOrbitRadius          , keplerOrbitTheta            , &
         &                                          keplerOrbitPhi           , keplerOrbitEpsilon            , keplerOrbitRadiusPericenter, keplerOrbitRadiusApocenter  , &
         &                                          keplerOrbitVelocityRadial, keplerOrbitVelocityTangential , keplerOrbitEnergy          , keplerOrbitAngularMomentum  , &
         &                                          keplerOrbitEccentricity  , keplerOrbitSemiMajorAxis      , keplerOrbitHostMass        , enumerationKeplerOrbitEncode
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
       select case (self%propertyIDs(i))
       case (keplerOrbitMasses             )
          call Galacticus_Error_Report('"masses" property is unsupported'//{introspection:location})
       case (keplerOrbitHostMass           )
          self%names_       (i)=prefix//'HostMass'
          self%descriptions_(i)='The mass of the host system [M☉].'
          self%unitsInSI_   (i)=massSolar
       case (keplerOrbitSpecificReducedMass)
          self%names_       (i)=prefix//'SpecificReducedMass'
          self%descriptions_(i)='The reduced mass divided by the mass of the orbiting object.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitRadius             )
          self%names_       (i)=prefix//'Radius'
          self%descriptions_(i)='The current orbital radius [Mpc].'
          self%unitsInSI_   (i)=megaParsec
       case (keplerOrbitTheta              )
          self%names_       (i)=prefix//'Theta'
          self%descriptions_(i)='Orbital angular position, θ, in spherical coordinates.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitPhi                )
          self%names_       (i)=prefix//'Phi'
          self%descriptions_(i)='Orbital angular position, ɸ, in spherical coordinates.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitEpsilon            )
          self%names_       (i)=prefix//'Epsilon'
          self%descriptions_(i)='Direction of the tangential component of velocity.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitRadiusPericenter   )
          self%names_       (i)=prefix//'RadiusPericenter'
          self%descriptions_(i)='Radius of the orbital pericenter [Mpc].'
          self%unitsInSI_   (i)=megaParsec
       case (keplerOrbitRadiusApocenter    )
          self%names_       (i)=prefix//'RadiusApocenter'
          self%descriptions_(i)='Radius of the orbital apocenter [Mpc].'
          self%unitsInSI_   (i)=megaParsec
       case (keplerOrbitVelocityRadial     )
          self%names_       (i)=prefix//'VelocityRadial'
          self%descriptions_(i)='Radial velocity of the orbit [km/s].'
          self%unitsInSI_   (i)=kilo
       case (keplerOrbitVelocityTangential )
          self%names_       (i)=prefix//'VelocityTangential'
          self%descriptions_(i)='Tangential velocity of the orbit [km/s].'
          self%unitsInSI_   (i)=kilo
       case (keplerOrbitEnergy             )
          self%names_       (i)=prefix//'Energy'
          self%descriptions_(i)='Energy of the orbit [M☉ (km/s)²].'
          self%unitsInSI_   (i)=massSolar*kilo**2
       case (keplerOrbitAngularMomentum    )
          self%names_       (i)=prefix//'AngularMomentum'
          self%descriptions_(i)='Angular momentum of the orbit [M☉ Mpc km/s].'
          self%unitsInSI_   (i)=massSolar*kilo*megaParsec
       case (keplerOrbitEccentricity       )
          self%names_       (i)=prefix//'Eccentricity'
          self%descriptions_(i)='Eccentricity of the orbit.'
          self%unitsInSI_   (i)=1.0d0
       case (keplerOrbitSemiMajorAxis      )
          self%names_       (i)=prefix//'SemiMajorAxis'
          self%descriptions_(i)='Semi-major axis of the orbit [Mpc].'
          self%unitsInSI_   (i)=megaParsec
       case default
          call Galacticus_Error_Report('unexpected property "'//trim(properties(i))//'"'//{introspection:location})
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
    use :: Kepler_Orbits, only : keplerOrbitHostMass      , keplerOrbitSpecificReducedMass, keplerOrbitRadius          , keplerOrbitTheta          , &
         &                       keplerOrbitPhi           , keplerOrbitEpsilon            , keplerOrbitRadiusPericenter, keplerOrbitRadiusApocenter, &
         &                       keplerOrbitVelocityRadial, keplerOrbitVelocityTangential , keplerOrbitEnergy          , keplerOrbitAngularMomentum, &
         &                       keplerOrbitEccentricity  , keplerOrbitSemiMajorAxis
    implicit none
    double precision                                  , dimension(:) , allocatable :: keplerOrbitExtractFromOrbit
    class           (nodePropertyExtractorKeplerOrbit), intent(inout), target      :: self
    type            (keplerOrbit                     ), intent(inout)              :: orbit
    integer                                                                        :: i

    allocate(keplerOrbitExtractFromOrbit(self%count_))
    if (orbit%isDefined()) then
       ! Orbit is defined - extract required properties.
       do i=1,self%count_
          select case (self%propertyIDs(i))
          case (keplerOrbitHostMass           )
             keplerOrbitExtractFromOrbit(i)=orbit%hostMass           ()
          case (keplerOrbitSpecificReducedMass)
             keplerOrbitExtractFromOrbit(i)=orbit%specificReducedMass()
          case (keplerOrbitRadius             )
             keplerOrbitExtractFromOrbit(i)=orbit%radius             ()
          case (keplerOrbitTheta              )
             keplerOrbitExtractFromOrbit(i)=orbit%theta              ()
          case (keplerOrbitPhi                )
             keplerOrbitExtractFromOrbit(i)=orbit%phi                ()
          case (keplerOrbitEpsilon            )
             keplerOrbitExtractFromOrbit(i)=orbit%epsilon            ()
          case (keplerOrbitRadiusPericenter   )
             keplerOrbitExtractFromOrbit(i)=orbit%radiusPericenter   ()
          case (keplerOrbitRadiusApocenter    )
             keplerOrbitExtractFromOrbit(i)=orbit%radiusApocenter    ()
          case (keplerOrbitVelocityRadial     )
             keplerOrbitExtractFromOrbit(i)=orbit%velocityRadial     ()
          case (keplerOrbitVelocityTangential )
             keplerOrbitExtractFromOrbit(i)=orbit%velocityTangential ()
          case (keplerOrbitEnergy             )
             keplerOrbitExtractFromOrbit(i)=orbit%energy             ()
          case (keplerOrbitAngularMomentum    )
             keplerOrbitExtractFromOrbit(i)=orbit%angularMomentum    ()
          case (keplerOrbitEccentricity       )
             keplerOrbitExtractFromOrbit(i)=orbit%eccentricity       ()
          case (keplerOrbitSemiMajorAxis      )
             keplerOrbitExtractFromOrbit(i)=orbit%semiMajorAxis      ()
          end select
       end do
    else
       ! Orbit is undefined - set all values to zero.
       keplerOrbitExtractFromOrbit=0.0d0
    end if
    return
  end function keplerOrbitExtractFromOrbit

  function keplerOrbitNames(self,time)
    !!{
    Return the names of the {\normalfont \ttfamily keplerOrbit} properties.
    !!}
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: keplerOrbitNames
    class           (nodePropertyExtractorKeplerOrbit), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(keplerOrbitNames(self%count_))
    keplerOrbitNames=self%names_
    return
  end function keplerOrbitNames

  function keplerOrbitDescriptions(self,time)
    !!{
    Return the descriptions of the {\normalfont \ttfamily virialProperies} properties.
    !!}
    implicit none
    type            (varying_string                  ), dimension(:) , allocatable :: keplerOrbitDescriptions
    class           (nodePropertyExtractorKeplerOrbit), intent(inout)              :: self
    double precision                                  , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(keplerOrbitDescriptions(self%count_))
    keplerOrbitDescriptions=self%descriptions_
    return
  end function keplerOrbitDescriptions

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

  integer function keplerOrbitType(self)
    !!{
    Return the type of the keplerOrbit property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorKeplerOrbit), intent(inout) :: self
    !$GLC attributes unused :: self

    keplerOrbitType=outputAnalysisPropertyTypeLinear
    return
  end function keplerOrbitType

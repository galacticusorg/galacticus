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

!!{
Contains a module which defines an orbit structure for use in \glc.
!!}

module Kepler_Orbits
  !!{
  Defines an orbit structure for use in \glc.
  !!}
  implicit none
  private
  public :: keplerOrbit

  ! Effective infinite radius used for apocenters in unbound orbits.
  double precision, parameter :: radiusEffectiveInfinity=huge(1.0d0)

  !![
  <enumeration>
   <name>keplerOrbit</name>
   <description>Properties of Kepler orbit objects.</description>
   <indexing>1</indexing>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <entry label="masses"              description="The masses of the orbiting bodies."                           />
   <entry label="massSatellite"       description="The mass of the satellite body."                              />
   <entry label="massHost"            description="The mass of the host body."                                   />
   <entry label="specificReducedMass" description="The specific reduced mass of the orbiting body."              />
   <entry label="radius"              description="The initial radius of the orbit."                             />
   <entry label="theta"               description="The initial polar angle, theta, of the orbit."                />
   <entry label="phi"                 description="The initial azimuthal angle, phi, of the orbit."              />
   <entry label="epsilon"             description="The initial rotation angle, epsilon, of the orbital velocity."/>
   <entry label="radiusPericenter"    description="The orbital pericenter radius."                               />
   <entry label="radiusApocenter"     description="The orbital apocenter radius."                                />
   <entry label="velocityRadial"      description="The initial radial velocity of the orbit."                    />
   <entry label="velocityTangential"  description="The initial tangential velocity of the orbit."                />
   <entry label="energy"              description="The orbital energy."                                          />
   <entry label="angularMomentum"     description="The orbital angular momentum"                                 />
   <entry label="eccentricity"        description="The orbital eccentricity."                                    />
   <entry label="semiMajorAxis"       description="The semi-major axis of the orbit."                            />
   <entry label="timeInitial"         description="The time at which the orbit was initialized."                 />
   <entry label="timeCurrent"         description="The current time for the orbit."                              />
  </enumeration>
  !!]

  type keplerOrbit
     !!{     
     The structure used for describing orbits in \glc. This object will automatically convert from one set of orbital parameters
     to another where possible. The orbiting bodies (a satellite orbiting around its host) are treated as point masses, and the
     usual ``reduced mass'' framework is used, such that radii and velocities are measured relative to a stationary host. Energy
     and angular momentum are defined per unit satellite mass (not per unit reduced mass). Note that not all interconversions
     between elements are implemented. The object works by attempting to get the radial and tangential velocities and the
     radius. If it can obtain these, any other parameter can be computed. Getting these three parameters relies on having known
     conversions from other possible combinations of parameters. The position of the object is described by $(r,\theta,\phi)$ in
     standard spherical coordinates. The direction of the tangential component is velocity is taken to be the direction of the
     vector $\mathbf{r} \times \mathbf{\hat{e}}_\mathrm{z}$, rotated by an angle $\epsilon$ around the vector $\mathbf{r}$. In
     cases where $\mathbf{r}$ is parallel to $\mathbf{\hat{e}}_\mathrm{z}$ we set $\epsilon=\pi/2$, since at some infinitesimally
     later time, $\delta t$ the position vector will become $\mathbf{r} + \mathbf{v}\delta t$ such that the above cross product is
     $(\mathbf{r} + \mathbf{v}\delta t) \times \mathbf{\hat{e}}_\mathrm{z} \propto \mathbf{v}_\mathrm{t} \times
     \mathbf{\hat{e}}_\mathrm{z}$ where $\mathbf{v}_\mathrm{t}$ is the tangential component of the velocity vector, and, by
     definition, we then must have $\epsilon=\pi/2$ as this vector is normal to $\mathbf{v}_\mathrm{t}$.
     !!}
     private
     double precision :: massHostValue        , specificReducedMassValue
     double precision :: radiusApocenterValue , radiusPericenterValue   , &
          &              radiusValue          , massSatelliteValue
     double precision :: velocityRadialValue  , velocityTangentialValue
     double precision :: thetaValue           , phiValue                , &
          &              epsilonValue
     double precision :: angularMomentumValue
     double precision :: energyValue
     double precision :: eccentricityValue
     double precision :: semimajorAxisValue
     logical          :: massesIsSet
     logical          :: radiusApocenterIsSet , radiusIsSet             , &
          &              radiusPericenterIsSet
     logical          :: thetaIsSet           , phiIsSet                , &
          &              epsilonIsSet
     logical          :: velocityRadialIsSet  , velocityTangentialIsSet
     logical          :: angularMomentumIsSet
     logical          :: energyIsSet
     logical          :: eccentricityIsSet
     logical          :: semimajorAxisIsSet
   contains
     ! Orbit methods.
     !![
     <methods>
       <method description="Build a Kepler orbit from an XML definition." method="builder" />
       <method description="Dump an orbit." method="dump" />
       <method description="Dump an orbit in binary." method="dumpRaw" />
       <method description="Read an orbit in binary." method="readRaw" />
       <method description="Resets orbit properties. If the optional {\normalfont \ttfamily keep} argument is provided and listed properties will \emph{not} be reset." method="reset" />
       <method description="Destroys an orbit." method="destroy" />
       <method description="Returns true if an orbit is fully defined." method="isDefined" />
       <method description="Asserts that an orbit is fully defined." method="assertIsDefined" />
       <method description="Returns true if the orbit is bound." method="isBound" />
       <method description="Propagates an orbit to a new position." method="propagate" />
       <method description="Sets the radial velocity of an orbit." method="velocityRadialSet" />
       <method description="Sets the masses of satellite and host objects." method="massesSet" />
       <method description="Sets the radius of an orbit." method="radiusSet" />
       <method description="Sets the angle $\theta$ of an orbit." method="thetaSet" />
       <method description="Sets the angle $\phi$ of an orbit." method="phiSet" />
       <method description="Sets the angle $\epsilon$ of an orbit." method="epsilonSet" />
       <method description="Sets the pericenter radius of an orbit." method="radiusPericenterSet" />
       <method description="Sets the apocenter radius of an orbit." method="radiusApocenterSet" />
       <method description="Sets the tangential velocity of an orbit." method="velocityTangentialSet" />
       <method description="Sets the energy of an orbit." method="energySet" />
       <method description="Sets the eccentricity of an orbit." method="eccentricitySet" />
       <method description="Sets the angular momentum of an orbit." method="angularMomentumSet" />
       <method description="Sets the semi-major axis of an orbit." method="semiMajorAxisSet" />
       <method description="Returns the host mass of an orbit." method="massHost" />
       <method description="Returns the satellite mass of an orbit." method="massSatellite" />
       <method description="Returns the velocity scale of an orbit." method="velocityScale" />
       <method description="Returns the specific reduced mass (i.e. the reduced mass per unit satellite mass, $\mu_\mathrm{s} = M_\mathrm{host}/(M_\mathrm{satellite}+M_\mathrm{host})$) of the orbit." method="specificReducedMass" />
       <method description="Returns the radius of an orbit." method="radius" />
       <method description="Returns the angle $\theta$ of an orbit." method="theta" />
       <method description="Returns the angle $\phi$ of an orbit." method="phi" />
       <method description="Returns the angle $\epsilon$ of an orbit." method="epsilon" />
       <method description="Returns the pericenter radius of an orbit." method="radiusPericenter" />
       <method description="Returns the apocenter radius of an orbit." method="radiusApocenter" />
       <method description="Returns the radial velocity of an orbit." method="velocityRadial" />
       <method description="Returns the tangential velocity of an orbit." method="velocityTangential" />
       <method description="Returns the energy of an orbit." method="energy" />
       <method description="Returns the eccentricity of an orbit." method="eccentricity" />
       <method description="Returns the angular momentum of an orbit." method="angularMomentum" />
       <method description="Returns the semi-major axis of an orbit." method="semiMajorAxis" />
       <method description="Returns the position coordinates." method="position" />
       <method description="Returns the velocity coordinates." method="velocity" />
       <method description="Returns the size of any non-static components of the type." method="nonStaticSizeOf" />
       <method description="Return true if two orbits are equivalent are equal." method="operator(==)" />
     </methods>
     !!]
     procedure :: builder               => Kepler_Orbits_Builder
     procedure :: dump                  => Kepler_Orbits_Dump
     procedure :: dumpRaw               => Kepler_Orbits_Dump_Raw
     procedure :: readRaw               => Kepler_Orbits_Read_Raw
     procedure :: nonStaticSizeOf       => Kepler_Orbits_Non_Static_Size_Of
     procedure :: reset                 => Kepler_Orbits_Reset
     procedure :: destroy               => Kepler_Orbits_Destroy
     procedure :: isDefined             => Kepler_Orbits_Is_Defined
     procedure :: assertIsDefined       => Kepler_Orbits_Assert_Is_Defined
     procedure :: isBound               => Kepler_Orbits_Is_Bound
     procedure :: propagate             => Kepler_Orbits_Propagate
     procedure :: massesSet             => Kepler_Orbits_Masses_Set
     procedure :: radiusSet             => Kepler_Orbits_Radius_Set
     procedure :: thetaSet              => Kepler_Orbits_Theta_Set
     procedure :: phiSet                => Kepler_Orbits_Phi_Set
     procedure :: epsilonSet            => Kepler_Orbits_Epsilon_Set
     procedure :: radiusPericenterSet   => Kepler_Orbits_Pericenter_Radius_Set
     procedure :: radiusApocenterSet    => Kepler_Orbits_Apocenter_Radius_Set
     procedure :: velocityRadialSet     => Kepler_Orbits_Velocity_Radial_Set
     procedure :: velocityTangentialSet => Kepler_Orbits_Velocity_Tangential_Set
     procedure :: energySet             => Kepler_Orbits_Energy_Set
     procedure :: angularMomentumSet    => Kepler_Orbits_Angular_Momentum_Set
     procedure :: eccentricitySet       => Kepler_Orbits_Eccentricity_Set
     procedure :: semiMajorAxisSet      => Kepler_Orbits_Semi_Major_Axis_Set
     procedure :: specificReducedMass   => Kepler_Orbits_Specific_Reduced_Mass
     procedure :: massHost              => Kepler_Orbits_Host_Mass
     procedure :: massSatellite         => Kepler_Orbits_Satellite_Mass
     procedure :: velocityScale         => Kepler_Orbits_Velocity_Scale
     procedure :: radius                => Kepler_Orbits_Radius
     procedure :: theta                 => Kepler_Orbits_Theta
     procedure :: phi                   => Kepler_Orbits_Phi
     procedure :: epsilon               => Kepler_Orbits_Epsilon
     procedure :: position              => Kepler_Orbits_Position
     procedure :: radiusPericenter      => Kepler_Orbits_Pericenter_Radius
     procedure :: radiusApocenter       => Kepler_Orbits_Apocenter_Radius
     procedure :: velocityRadial        => Kepler_Orbits_Velocity_Radial
     procedure :: velocityTangential    => Kepler_Orbits_Velocity_Tangential
     procedure :: velocity              => Kepler_Orbits_Velocity
     procedure :: energy                => Kepler_Orbits_Energy
     procedure :: angularMomentum       => Kepler_Orbits_Angular_Momentum
     procedure :: eccentricity          => Kepler_Orbits_Eccentricity
     procedure :: semiMajorAxis         => Kepler_Orbits_Semi_Major_Axis
     procedure ::                          Kepler_Orbits_Equivalent
     generic   :: operator(==)          => Kepler_Orbits_Equivalent
  end type keplerOrbit

  interface keplerOrbit
     !!{
     Constructors for Kepler orbits.
     !!}
     module procedure keplerOrbitConstructorNull
  end interface keplerOrbit

  ! A null orbit.
  type(keplerOrbit), public :: zeroKeplerOrbit=keplerOrbit(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.)

contains

  function keplerOrbitConstructorNull() result(self)
    !!{
    Null constructor for Kepler orbit objects.
    !!}
    type(keplerOrbit) :: self

    self=zeroKeplerOrbit
    return
  end function keplerOrbitConstructorNull

  subroutine Kepler_Orbits_Destroy(orbit)
    !!{
    Destroy an orbit.
    !!}
    implicit none
    class(keplerOrbit), intent(inout) :: orbit
    !$GLC attributes unused :: orbit

    ! Nothing to do.
    return
  end subroutine Kepler_Orbits_Destroy

  subroutine Kepler_Orbits_Builder(self,keplerOrbitDefinition)
    !!{
    Build a {\normalfont \ttfamily keplerOrbit} object from the given XML {\normalfont \ttfamily keplerOrbitDefinition}.
    !!}
    use :: FoX_DOM, only : getNodeName                 , node       , extractDataContent
    use :: Error  , only : Error_Report
    use :: IO_XML , only : XML_Get_Elements_By_Tag_Name, xmlNodeList
    implicit none
    class           (keplerOrbit), intent(inout)               :: self
    type            (node       ), pointer                     :: keplerOrbitDefinition
    type            (node       ), pointer                     :: property
    type            (xmlNodeList), allocatable  , dimension(:) :: propertyList
    character       (len=18     ), parameter    , dimension(5) :: propertyNames        =                                         &
         &                                                                              [                                        &
         &                                                                               'massHost          ',                   &
         &                                                                               'massSatellite     ',                   &
         &                                                                               'velocityRadial    ',                   &
         &                                                                               'velocityTangential',                   &
         &                                                                               'radius            '                    &
         &                                                                              ]
    integer                                                    :: i
    double precision                                           :: massHost                                   , massSatellite   , &
         &                                                        propertyValue
    logical                                                    :: massHostSet                                , massSatelliteSet
    character       (len=32     )                              :: nodeName
    
    ! Get the radius.
    do i=1,size(propertyNames)
       !$omp critical (FoX_DOM_Access)
       call XML_Get_Elements_By_Tag_Name(keplerOrbitDefinition,trim(propertyNames(i)),propertyList)
       !$omp end critical (FoX_DOM_Access)
       if (size(propertyList) >  1) call Error_Report('multiple '//trim(propertyNames(i))//' values specified'//{introspection:location})
       if (size(propertyList) == 1) then
          !$omp critical (FoX_DOM_Access)
          property => propertyList(0)%element
          call extractDataContent(property,propertyValue)
          nodeName=getNodeName(property)
          !$omp end critical (FoX_DOM_Access)
          select case (trim(nodeName))
          case ( 'radius'             )
             call self%            radiusSet(propertyValue)
          case ( 'velocityRadial'     )
             call self%    velocityRadialSet(propertyValue)
          case ( 'velocityTangential' )
             call self%velocityTangentialSet(propertyValue)
          case ( 'massHost'           )
             massHost        =propertyValue
             massHostSet     =.true.
          case ( 'massSatellite'      )
             massSatellite   =propertyValue
             massSatelliteSet=.true.
          case default
             call Error_Report('unrecognized property name'//{introspection:location})
          end select
       end if
    end do
    if (.not.(massHostSet.and.massSatelliteSet)) call Error_Report('satellite and host masses must be specified'//{introspection:location})
    call self%massesSet(massSatellite,massHost)
    call self%assertIsDefined()
    return
  end subroutine Kepler_Orbits_Builder

  subroutine Kepler_Orbits_Dump(self,verbosityLevel)
    !!{
    Reset an orbit to a null state.
    !!}
    use :: Display           , only : displayMessage, enumerationVerbosityLevelType
    use :: ISO_Varying_String, only : assignment(=) , varying_string
    implicit none
    class    (keplerOrbit                  ), intent(in   ) :: self
    type     (enumerationVerbosityLevelType), intent(in   ) :: verbosityLevel
    character(len=22                       )                :: label
    type     (varying_string               )                :: message

    if (self%massesIsSet             ) then
       write (label,'(e22.16)') self%massHostValue
       message='host mass:             '//label
       call displayMessage(message,verbosityLevel)
       write (label,'(e22.16)') self%specificReducedMassValue
       message='specific reduced mass: '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%radiusIsSet             ) then
       write (label,'(e22.16)') self%radiusValue
       message='radius:                '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%radiusPericenterIsSet   ) then
       write (label,'(e22.16)') self%radiusPericenterValue
       message='radius pericenter:     '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%radiusApocenterIsSet   ) then
       write (label,'(e22.16)') self%radiusApocenterValue
       message='radius apocenter:      '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%velocityRadialIsSet    ) then
       write (label,'(e22.16)') self%velocityRadialValue
       message='velocity radial:       '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%velocityTangentialIsSet) then
       write (label,'(e22.16)') self%velocityTangentialValue
       message='velocity tangential:   '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%angularMomentumIsSet   ) then
       write (label,'(e22.16)') self%angularMomentumValue
       message='angular momentum:      '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%energyIsSet            ) then
       write (label,'(e22.16)') self%energyValue
       message='energy:                '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%eccentricityIsSet      ) then
       write (label,'(e22.16)') self%eccentricityValue
       message='eccentricity:          '//label
       call displayMessage(message,verbosityLevel)
    end if
    if (self%semimajorAxisIsSet     ) then
       write (label,'(e22.16)') self%semimajorAxisValue
       message='semi-major axis:       '//label
       call displayMessage(message,verbosityLevel)
    end if
    return
  end subroutine Kepler_Orbits_Dump

  subroutine Kepler_Orbits_Dump_Raw(self,fileHandle)
    !!{
    Dump a {\normalfont \ttfamily keplerOrbit} object in binary.
    !!}
    implicit none
    class  (keplerOrbit), intent(in   ) :: self
    integer             , intent(in   ) :: fileHandle

    write (fileHandle) self%massesIsSet,self%massesIsSet,self%radiusIsSet,self%radiusPericenterIsSet,self%radiusApocenterIsSet&
         &,self%velocityRadialIsSet,self%velocityTangentialIsSet,self%angularMomentumIsSet,self%energyIsSet&
         &,self%eccentricityIsSet,self%semimajorAxisIsSet
    if (self%massesIsSet            ) write (fileHandle) self%massHostValue,self%specificReducedMassValue
    if (self%radiusIsSet            ) write (fileHandle) self%radiusValue
    if (self%radiusPericenterIsSet  ) write (fileHandle) self%radiusPericenterValue
    if (self%radiusApocenterIsSet   ) write (fileHandle) self%radiusApocenterValue
    if (self%velocityRadialIsSet    ) write (fileHandle) self%velocityRadialValue
    if (self%velocityTangentialIsSet) write (fileHandle) self%velocityTangentialValue
    if (self%angularMomentumIsSet   ) write (fileHandle) self%angularMomentumValue
    if (self%energyIsSet            ) write (fileHandle) self%energyValue
    if (self%eccentricityIsSet      ) write (fileHandle) self%eccentricityValue
    if (self%semimajorAxisIsSet     ) write (fileHandle) self%semimajorAxisValue
    return
  end subroutine Kepler_Orbits_Dump_Raw

  subroutine Kepler_Orbits_Read_Raw(self,fileHandle)
    !!{
    Read a {\normalfont \ttfamily keplerOrbit} object in binary.
    !!}
    implicit none
    class  (keplerOrbit), intent(inout) :: self
    integer             , intent(in   ) :: fileHandle

    read (fileHandle) self%massesIsSet,self%massesIsSet,self%radiusIsSet,self%radiusPericenterIsSet,self%radiusApocenterIsSet&
         &,self%velocityRadialIsSet,self%velocityTangentialIsSet,self%angularMomentumIsSet,self%energyIsSet&
         &,self%eccentricityIsSet,self%semimajorAxisIsSet
    if (self%massesIsSet            ) read (fileHandle) self%massHostValue,self%specificReducedMassValue
    if (self%radiusIsSet            ) read (fileHandle) self%radiusValue
    if (self%radiusPericenterIsSet  ) read (fileHandle) self%radiusPericenterValue
    if (self%radiusApocenterIsSet   ) read (fileHandle) self%radiusApocenterValue
    if (self%velocityRadialIsSet    ) read (fileHandle) self%velocityRadialValue
    if (self%velocityTangentialIsSet) read (fileHandle) self%velocityTangentialValue
    if (self%angularMomentumIsSet   ) read (fileHandle) self%angularMomentumValue
    if (self%energyIsSet            ) read (fileHandle) self%energyValue
    if (self%eccentricityIsSet      ) read (fileHandle) self%eccentricityValue
    if (self%semimajorAxisIsSet     ) read (fileHandle) self%semimajorAxisValue
    return
  end subroutine Kepler_Orbits_Read_Raw

  subroutine Kepler_Orbits_Reset(orbit,keep)
    !!{
    Reset an orbit to a null state.
    !!}
    implicit none
    class  (keplerOrbit               ), intent(inout)                         :: orbit
    type   (enumerationKeplerOrbitType), intent(in   ), dimension(:), optional :: keep
    type   (enumerationKeplerOrbitType), allocatable  , dimension(:)           :: keep_

    ! Set list of properties to keep.
    if (present(keep)) then
       keep_=keep
    else
       allocate(keep_(0))
    end if
    ! Simply specify that no properties have been set as yet.
    if (.not.any(keep_ == keplerOrbitMasses            )) orbit%massesIsSet            =.false.
    if (.not.any(keep_ == keplerOrbitRadius            )) orbit%radiusIsSet            =.false.
    if (.not.any(keep_ == keplerOrbitRadiusPericenter  )) orbit%radiusPericenterIsSet  =.false.
    if (.not.any(keep_ == keplerOrbitRadiusApocenter   )) orbit%radiusApocenterIsSet   =.false.
    if (.not.any(keep_ == keplerOrbitVelocityRadial    )) orbit%velocityRadialIsSet    =.false.
    if (.not.any(keep_ == keplerOrbitVelocityTangential)) orbit%velocityTangentialIsSet=.false.
    if (.not.any(keep_ == keplerOrbitAngularMomentum   )) orbit%angularMomentumIsSet   =.false.
    if (.not.any(keep_ == keplerOrbitEnergy            )) orbit%energyIsSet            =.false.
    if (.not.any(keep_ == keplerOrbitEccentricity      )) orbit%eccentricityIsSet      =.false.
    if (.not.any(keep_ == keplerOrbitSemiMajorAxis     )) orbit%semimajorAxisIsSet     =.false.
    if (.not.any(keep_ == keplerOrbitTheta             )) orbit%thetaIsSet             =.false.
    if (.not.any(keep_ == keplerOrbitPhi               )) orbit%phiIsSet               =.false.
    return
  end subroutine Kepler_Orbits_Reset

  subroutine Kepler_Orbits_Masses_Set(orbit,massSatellite,massHost)
    !!{
    Sets the masses of the two orbiting objects in a {\normalfont \ttfamily keplerOrbit} object.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: massHost , massSatellite

    ! Set the mass factor and flag that is set.
    orbit%specificReducedMassValue=1.0d0/(1.0d0+massSatellite/massHost)
    orbit%massHostValue           =massHost
    orbit%massSatelliteValue      =massSatellite
    orbit%massesIsSet             =.true.
    return
  end subroutine Kepler_Orbits_Masses_Set

  subroutine Kepler_Orbits_Radius_Set(orbit,radius)
    !!{
    Sets the radius to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: radius

    ! Set the radius and flag that is set.
    orbit%radiusValue=radius
    orbit%radiusIsSet=.true.
    return
  end subroutine Kepler_Orbits_Radius_Set

  subroutine Kepler_Orbits_Pericenter_Radius_Set(orbit,radius)
    !!{
    Sets the pericenter radius to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: radius

    ! Set the pericenter radius and flag that is set.
    orbit%radiusPericenterValue=radius
    orbit%radiusPericenterIsSet=.true.
    return
  end subroutine Kepler_Orbits_Pericenter_Radius_Set

  subroutine Kepler_Orbits_Apocenter_Radius_Set(orbit,radius)
    !!{
    Sets the apocenter radius to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: radius

    ! Set the apocenter radius and flag that is set.
    orbit%radiusApocenterValue=radius
    orbit%radiusApocenterIsSet=.true.
    return
  end subroutine Kepler_Orbits_Apocenter_Radius_Set

  subroutine Kepler_Orbits_Velocity_Radial_Set(orbit,velocityRadial)
    !!{
    Sets the radial velocity to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: velocityRadial

    ! Set the radial velocity and flag that is set.
    orbit%velocityRadialValue=velocityRadial
    orbit%velocityRadialIsSet=.true.
    return
  end subroutine Kepler_Orbits_Velocity_Radial_Set

  subroutine Kepler_Orbits_Velocity_Tangential_Set(orbit,velocityTangential)
    !!{
    Sets the tangential velocity to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: velocityTangential

    ! Set the tangential velocity and flag that is set.
    orbit%velocityTangentialValue=velocityTangential
    orbit%velocityTangentialIsSet=.true.
    return
  end subroutine Kepler_Orbits_Velocity_Tangential_Set

  subroutine Kepler_Orbits_Energy_Set(orbit,energy)
    !!{
    Sets the tangential velocity to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: energy

    ! Set the tangential velocity and flag that is set.
    orbit%energyValue=energy
    orbit%energyIsSet=.true.
    return
  end subroutine Kepler_Orbits_Energy_Set

  subroutine Kepler_Orbits_Eccentricity_Set(orbit,eccentricity)
    !!{
    Sets the tangential velocity to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: eccentricity

    ! Set the tangential velocity and flag that is set.
    orbit%eccentricityValue=eccentricity
    orbit%eccentricityIsSet=.true.
    return
  end subroutine Kepler_Orbits_Eccentricity_Set

  subroutine Kepler_Orbits_Angular_Momentum_Set(orbit,angularMomentum)
    !!{
    Sets the tangential velocity to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: angularMomentum

    ! Set the tangential velocity and flag that is set.
    orbit%angularMomentumValue=angularMomentum
    orbit%angularMomentumIsSet=.true.
    return
  end subroutine Kepler_Orbits_Angular_Momentum_Set

  subroutine Kepler_Orbits_Semi_Major_Axis_Set(orbit,semiMajorAxis)
    !!{
    Sets the semi-major axis to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: semiMajorAxis

    ! Set the tangential velocity and flag that is set.
    orbit%semiMajorAxisValue=semiMajorAxis
    orbit%semiMajorAxisIsSet=.true.
    return
  end subroutine Kepler_Orbits_Semi_Major_Axis_Set

  subroutine Kepler_Orbits_Theta_Set(orbit,theta)
    !!{
    Sets the angle $\theta$ to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: theta

    ! Set the theta and flag that is set.
    orbit%thetaValue=theta
    orbit%thetaIsSet=.true.
    return
  end subroutine Kepler_Orbits_Theta_Set

  subroutine Kepler_Orbits_Phi_Set(orbit,phi)
    !!{
    Sets the angle $\phi$ to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: phi

    ! Set the phi and flag that is set.
    orbit%phiValue=phi
    orbit%phiIsSet=.true.
    return
  end subroutine Kepler_Orbits_Phi_Set

  subroutine Kepler_Orbits_Epsilon_Set(orbit,epsilon)
    !!{
    Sets the $\epsilon$ to the specified value.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision             , intent(in   ) :: epsilon

    ! Set ε and flag that is set.
    orbit%epsilonValue=epsilon
    orbit%epsilonIsSet=.true.
    return
  end subroutine Kepler_Orbits_Epsilon_Set

  double precision function Kepler_Orbits_Theta(orbit)
    !!{
    Return the angle $\theta$ for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%thetaIsSet) call Error_Report('theta has not been set for this orbit'//{introspection:location})
    Kepler_Orbits_Theta=orbit%thetaValue
    return
  end function Kepler_Orbits_Theta

  double precision function Kepler_Orbits_Phi(orbit)
    !!{
    Return the angle $\phi$ for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%phiIsSet) call Error_Report('phi has not been set for this orbit'//{introspection:location})
    Kepler_Orbits_Phi=orbit%phiValue
    return
  end function Kepler_Orbits_Phi

  double precision function Kepler_Orbits_Epsilon(orbit)
    !!{
    Return the angle $\epsilon$ for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%epsilonIsSet) call Error_Report('epsilon has not been set for this orbit'//{introspection:location})
    Kepler_Orbits_Epsilon=orbit%epsilonValue
    return
  end function Kepler_Orbits_Epsilon

  double precision function Kepler_Orbits_Specific_Reduced_Mass(orbit)
    !!{
    Return the specific reduced mass for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%massesIsSet) call Error_Report('mass factor has not been set for this orbit'//{introspection:location})
    Kepler_Orbits_Specific_Reduced_Mass=orbit%specificReducedMassValue
    return
  end function Kepler_Orbits_Specific_Reduced_Mass

  double precision function Kepler_Orbits_Satellite_Mass(orbit)
    !!{
    Return the mass for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%massesIsSet) call Error_Report('satellite mass has not been set for this orbit'//{introspection:location})
    Kepler_Orbits_Satellite_Mass=orbit%massSatelliteValue
    return
  end function Kepler_Orbits_Satellite_Mass

  double precision function Kepler_Orbits_Host_Mass(orbit)
    !!{
    Return the host mass for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%massesIsSet) call Error_Report('host mass has not been set for this orbit'//{introspection:location})
    Kepler_Orbits_Host_Mass=orbit%massHostValue
    return
  end function Kepler_Orbits_Host_Mass

  double precision function Kepler_Orbits_Radius(orbit)
    !!{
    Return the radius for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%radiusIsSet) call Error_Report('radius has not been set for this orbit'//{introspection:location})
    Kepler_Orbits_Radius=orbit%radiusValue
    return
  end function Kepler_Orbits_Radius

  double precision function Kepler_Orbits_Pericenter_Radius(orbit)
    !!{
    Return the pericenter radius for this orbit.
    !!}
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%radiusPericenterIsSet) then
       ! Assert that the orbit is defined.
       call orbit%assertIsDefined()
       ! Compute the pericenter radius.
       orbit%radiusPericenterValue=orbit%semiMajorAxis()*(1.0d0-orbit%eccentricity())
       orbit%radiusPericenterIsSet=.true.
    end if
    Kepler_Orbits_Pericenter_Radius=orbit%radiusPericenterValue
    return
  end function Kepler_Orbits_Pericenter_Radius

  double precision function Kepler_Orbits_Apocenter_Radius(orbit)
    !!{
    Return the apocenter radius for this orbit.
    !!}
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%radiusApocenterIsSet) then
       ! Assert that the orbit is defined.
       call orbit%assertIsDefined()
       ! Compute the pericenter radius.
       if (orbit%isBound()) then
          orbit%radiusApocenterValue=orbit%semiMajorAxis()*(1.0d0+orbit%eccentricity())
       else
          orbit%radiusApocenterValue=radiusEffectiveInfinity
       end if
       orbit%radiusApocenterIsSet=.true.
    end if
    Kepler_Orbits_Apocenter_Radius=orbit%radiusApocenterValue
    return
  end function Kepler_Orbits_Apocenter_Radius

  double precision function Kepler_Orbits_Velocity_Radial(orbit)
    !!{
    Return the radial velocity for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%velocityRadialIsSet) then
       ! Assert that the orbit is defined.
       call orbit%assertIsDefined()
       ! Compute from eccentricity, radius and periapsis if possible.
       if (orbit%eccentricityIsSet.and.orbit%radiusIsSet.and.orbit%radiusPericenterIsSet) then
          orbit%velocityRadialValue=orbit%velocityScale()*sqrt((2.0d0*(1.0d0-orbit%radius()&
               &/orbit%radiusPericenter())+(1.0d0+orbit%eccentricity())*(orbit%radius()/orbit%radiusPericenter()&
               &-orbit%radiusPericenter()/orbit%radius()))/orbit%specificReducedMass())
          orbit%velocityRadialIsSet=.true.
       end if
       ! If we were not able to compute the radial velocity, exit.
       if (.not.orbit%velocityRadialIsSet) call Error_Report('radial velocity has not been set for this orbit and can not be computed'//{introspection:location})
    end if
    Kepler_Orbits_Velocity_Radial=orbit%velocityRadialValue
    return
  end function Kepler_Orbits_Velocity_Radial

  double precision function Kepler_Orbits_Velocity_Tangential(orbit)
    !!{
    Return the tangential velocity for this orbit.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%velocityTangentialIsSet) then
       ! Assert that the orbit is defined.
       call orbit%assertIsDefined()
       ! Compute from angular momentum and radius if possible.
       if (orbit%angularMomentumIsSet.and.orbit%radiusIsSet) then
          orbit%velocityTangentialValue=orbit%angularMomentum()/orbit%radius()/orbit%specificReducedMass()
          orbit%velocityTangentialIsSet=.true.
       end if
       ! Compute from eccentricity, radius and periapsis if possible.
       if (orbit%eccentricityIsSet.and.orbit%radiusIsSet.and.orbit%radiusPericenterIsSet) then
          orbit%velocityTangentialValue=orbit%velocityScale()*sqrt((1.0d0+orbit%eccentricity())&
               &*orbit%radiusPericenter() /orbit%radius() /orbit%specificReducedMass())
          orbit%velocityTangentialIsSet=.true.
       end if
       ! If we were not able to compute the tangential velocity, exit.
       if (.not.orbit%velocityTangentialIsSet) call Error_Report('tangential velocity has not been set for this orbit and can not be computed'//{introspection:location})
    end if
    Kepler_Orbits_Velocity_Tangential=orbit%velocityTangentialValue
    return
  end function Kepler_Orbits_Velocity_Tangential

  double precision function Kepler_Orbits_Energy(orbit)
    !!{
    Return the energy for this orbit.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    ! Check if energy is set.
    if (.not.orbit%energyIsSet) then
       ! Assert that the orbit is defined.
       call orbit%assertIsDefined()
       ! Compute the energy.
       orbit%energyValue=-gravitationalConstant_internal*orbit%massHost()/orbit%radius()+0.5d0&
            &*(orbit%velocityRadial()**2+orbit%velocityTangential()**2)*orbit%specificReducedMass()
       orbit%energyIsSet=.true.
    end if
    Kepler_Orbits_Energy=orbit%energyValue
    return
  end function Kepler_Orbits_Energy

  double precision function Kepler_Orbits_Angular_Momentum(orbit)
    !!{
    Return the angular momentum for this orbit.
    !!}
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    if (.not.orbit%angularMomentumIsSet) then
       orbit%angularMomentumValue=orbit%radius()*orbit%velocityTangential()*orbit%specificReducedMass()
       orbit%angularMomentumIsSet=.true.
    end if
    Kepler_Orbits_Angular_Momentum=orbit%angularMomentumValue
    return
  end function Kepler_Orbits_Angular_Momentum

  double precision function Kepler_Orbits_Eccentricity(orbit)
    !!{
    Return the eccentricity for this orbit.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision                             :: velocityRadial, velocityTangential

    if (.not.orbit%eccentricityIsSet) then
       velocityTangential         =orbit%velocityTangential()/orbit%velocityScale()
       velocityRadial             =orbit%velocityRadial    ()/orbit%velocityScale()
       orbit%eccentricityValue=sqrt(1.0d0+2.0d0*orbit%energy()*orbit%angularMomentum()**2/orbit%radius()**2&
            &/orbit%velocityScale()**4/orbit%specificReducedMass())
       orbit%eccentricityIsSet=.true.
    end if
    Kepler_Orbits_Eccentricity=orbit%eccentricityValue
    return
  end function Kepler_Orbits_Eccentricity

  double precision function Kepler_Orbits_Semi_Major_Axis(orbit)
    !!{
    Return the semi-major axis for this orbit.
    !!}
    implicit none
    class           (keplerOrbit), intent(inout) :: orbit
    double precision                             :: velocityRadial, velocityTangential

    if (.not.orbit%semiMajorAxisIsSet) then
       velocityTangential          =orbit%velocityTangential()/orbit%velocityScale()
       velocityRadial              =orbit%velocityRadial    ()/orbit%velocityScale()
       orbit%semiMajorAxisValue=orbit%radius()/orbit%specificReducedMass()/(2.0d0/orbit%specificReducedMass()&
            &-velocityRadial**2-velocityTangential**2)
       orbit%semiMajorAxisIsSet=.true.
    end if
    Kepler_Orbits_Semi_Major_Axis=orbit%semiMajorAxisValue
    return
  end function Kepler_Orbits_Semi_Major_Axis

  logical function Kepler_Orbits_Is_Defined(orbit)
    !!{
    Returns true if the orbit is fully defined. For the orbits consider here, in which we don't care about the orientation of
    the orbital plane or the argument of pericenter, this requires that three orbital parameter be set (in addition to the
    masses of the orbiting bodies).
    !!}
    implicit none
    class  (keplerOrbit), intent(in   ) :: orbit
    integer                             :: orbitalParameterCount

    ! Assume orbit is not defined by default.
    Kepler_Orbits_Is_Defined=.false.
    ! Ensure that masses are set.
    if (orbit%massesIsSet) then
       ! Count how many parameters are set.
       orbitalParameterCount=0
       if (orbit%radiusIsSet            ) orbitalParameterCount=orbitalParameterCount+1
       if (orbit%radiusPericenterIsSet  ) orbitalParameterCount=orbitalParameterCount+1
       if (orbit%radiusApocenterIsSet   ) orbitalParameterCount=orbitalParameterCount+1
       if (orbit%energyIsSet            ) orbitalParameterCount=orbitalParameterCount+1
       if (orbit%angularMomentumIsSet   ) orbitalParameterCount=orbitalParameterCount+1
       if (orbit%velocityRadialIsSet    ) orbitalParameterCount=orbitalParameterCount+1
       if (orbit%velocityTangentialIsSet) orbitalParameterCount=orbitalParameterCount+1
       if (orbit%eccentricityIsSet      ) orbitalParameterCount=orbitalParameterCount+1
       if (orbit%semiMajorAxisIsSet     ) orbitalParameterCount=orbitalParameterCount+1
       ! Orbit is defined if at least 3 parameters are set.
       Kepler_Orbits_Is_Defined=(orbitalParameterCount >= 3)
    end if
    return
  end function Kepler_Orbits_Is_Defined

  subroutine Kepler_Orbits_Assert_Is_Defined(orbit)
    !!{
    Assert that an orbit is defined - quit with an error if it is not.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(keplerOrbit), intent(in   ) :: orbit

    if (.not.orbit%isDefined()) call Error_Report('orbit is not defined'//{introspection:location})
    return
  end subroutine Kepler_Orbits_Assert_Is_Defined

  logical function Kepler_Orbits_Is_Bound(orbit)
    !!{
    Returns true if the orbit is bound.
    !!}
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    ! Assert that the orbit is defined.
    call orbit%assertIsDefined()
    ! Test if the energy is negative (indicating a bound orbit).
    Kepler_Orbits_Is_Bound=(orbit%energy() < 0.0d0)
    return
  end function Kepler_Orbits_Is_Bound

  double precision function Kepler_Orbits_Velocity_Scale(orbit)
    !!{
    Return the velocity scale for the orbit.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(keplerOrbit), intent(inout) :: orbit

    ! Check that masses and radius have been specified.
    if (.not.(orbit%radiusIsSet.and.orbit%massesIsSet)) call Error_Report('orbit masses and radius must be specified'//{introspection:location})
    ! Compute the velocity scale.
    Kepler_Orbits_Velocity_Scale=sqrt(gravitationalConstant_internal*orbit%massHost()/orbit%radius())
    return
  end function Kepler_Orbits_Velocity_Scale

  subroutine Kepler_Orbits_Propagate(orbit,newRadius,infalling)
    !!{
    Propagate an orbit along its path.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Vectors                         , only : Vector_Magnitude              , Vector_Product
    implicit none
    class           (keplerOrbit), intent(inout)           :: orbit
    double precision             , intent(in   )           :: newRadius
    logical                      , intent(in   ), optional :: infalling
    double precision             , dimension(3)            :: position             , positionNew          , &
         &                                                    velocity             , velocityNew          , &
         &                                                    vectorEpsilon1       , vectorEpsilon2       , &
         &                                                    vectorNormal         , vectorPlane          , &
         &                                                    vectorRadial         , vectorRadialNew      , &
         &                                                    vectorTangential1    , vectorTangential2    , &
         &                                                    velocityTangentialNew
    double precision                                       :: angularMomentum      , energy               , &
         &                                                    newVelocityRadial    , newVelocityTangential, &
         &                                                    eccentricity         , semiMajorAxis        , &
         &                                                    thetaInitial         , thetaFinal           , &
         &                                                    velocityEpsilon1     , velocityEpsilon2

    ! Assert that the orbit is defined.
    call orbit%assertIsDefined()
    ! Check that the radius is within the allowed range.
    if     (                                                                                                    &
         &    newRadius < orbit%radiusPericenter()                                                              &
         &  .or.                                                                                                &
         &   (newRadius > orbit%radiusApocenter () .and. orbit%radiusApocenter() > 0.0d0)                       &
         & ) call Error_Report('radius lies outside of allowed range for this orbit'//{introspection:location})
    ! Get the energy and angular momentum
    energy         =orbit%energy         ()
    angularMomentum=orbit%angularMomentum()
    ! Compute velocity components.
    newVelocityTangential=angularMomentum/newRadius/orbit%specificReducedMass()
    newVelocityRadial    =sqrt(2.0d0*(energy+gravitationalConstant_internal*orbit%massHost()/newRadius)/orbit%specificReducedMass()-newVelocityTangential**2)
    ! Move to the infalling phase of the orbit if requested.
    if (present(infalling)) then
       if (infalling) newVelocityRadial=-newVelocityRadial
    end if
    ! If position information is set, propagate that also.
    if     (                    &
         &   orbit%  thetaIsSet &
         &  .and.               &
         &   orbit%    phiIsSet &
         &  .and.               &
         &   orbit%epsilonIsSet &
         & ) then
       ! Get orbital eccentricity and semi-major axis.
       eccentricity =orbit%eccentricity ()
       semiMajorAxis=orbit%semiMajorAxis()
       ! Get initial position and velocity.
       vectorRadial     = [sin(orbit%theta())*cos(orbit%phi()),sin(orbit%theta())*sin(orbit%phi()),cos(orbit%theta())]
       vectorTangential1= Vector_Product(vectorRadial     ,[0.0d0,0.0d0,1.0d0])
       ! If the radial vector happens to lie precisely along the z-axis, try constructing a tangent vector using a different basis
       ! vector.
       if (Vector_Magnitude(vectorTangential1) <= 0.0d0) &
            & vectorTangential1= Vector_Product(vectorRadial     ,[0.0d0,1.0d0,0.0d0])
        vectorTangential2= Vector_Product(vectorTangential1,vectorRadial       )
       vectorTangential1= vectorTangential1/Vector_Magnitude(vectorTangential1)
       vectorTangential2= vectorTangential2/Vector_Magnitude(vectorTangential2)
       position         =+  orbit%radius            ()             &
            &            *  vectorRadial
       velocity         =+  orbit%velocityRadial    ()             &
            &            *  vectorRadial                           &
            &            +  orbit%velocityTangential()             &
            &            *(                                        &
            &              +vectorTangential1*cos(orbit%epsilon()) &
            &              +vectorTangential2*sin(orbit%epsilon()) &
            &             )
       ! Find the vector defining the orbital plane.
       vectorNormal=+Vector_Product  (position,velocity)
       vectorNormal=+vectorNormal                        &
            &       /Vector_Magnitude(vectorNormal)
       ! Find the orbital angle corresponding to the initial position.
       thetaInitial=acos(max(-1.0d0,min(1.0d0,((1.0d0-eccentricity**2)*semiMajorAxis/orbit%radius    ()-1.0d0)/eccentricity)))
       if (orbit%velocityRadial() < 0.0d0) thetaInitial=-thetaInitial
       ! Find the orbital angle corresponding to the new position.
       thetaFinal  =acos(max(-1.0d0,min(1.0d0,((1.0d0-eccentricity**2)*semiMajorAxis/       newRadius  -1.0d0)/eccentricity)))
       if (newVelocityRadial      < 0.0d0) thetaFinal  =-thetaFinal
       ! Rotate the radial vector around the normal vector by the difference in angle.
       vectorPlane    = Vector_Product(vectorRadial,vectorNormal)
       vectorRadialNew=+vectorRadial*cos(thetaFinal-thetaInitial) &
            &          +vectorPlane* sin(thetaFinal-thetaInitial)
       ! Find the new position.
       positionNew=vectorRadialNew*newRadius
       ! Construct the new velocity.
       velocityTangentialNew=newVelocityTangential*vectorPlane
       velocityNew          =vectorRadialNew*newVelocityRadial+velocityTangentialNew
       ! Construct epsilon vectors.
       vectorEpsilon1  =Vector_Product(vectorRadialNew,[0.0d0,0.0d0,1.0d0])
       vectorEpsilon2  =Vector_Product(vectorEpsilon1 ,vectorRadialNew    )
       vectorEpsilon1  =vectorEpsilon1/Vector_Magnitude(vectorEpsilon1)
       vectorEpsilon2  =vectorEpsilon2/Vector_Magnitude(vectorEpsilon2)
       velocityEpsilon1=Dot_Product(velocityTangentialNew,vectorEpsilon1)
       velocityEpsilon2=Dot_Product(velocityTangentialNew,vectorEpsilon2)
       ! Set new angles.
       call orbit%thetaSet  (acos (positionNew        (3)/newRadius                    ))
       call orbit%phiSet    (atan2(positionNew        (2)          ,positionNew     (1)))
       call orbit%epsilonSet(atan2(velocityEpsilon2                ,velocityEpsilon1   ))
    else
       ! Incomplete position information present - mark position as unset.
       orbit%  thetaIsSet=.false.
       orbit%    phiIsSet=.false.
       orbit%epsilonIsSet=.false.
    end if
    ! Set new values for the state vector.
    call orbit%radiusSet            (newRadius            )
    call orbit%velocityTangentialSet(newVelocityTangential)
    call orbit%velocityRadialSet    (newVelocityRadial    )
    return
  end subroutine Kepler_Orbits_Propagate

  function Kepler_Orbits_Position(orbit) result(position)
    !!{
    Return the position of the orbit in Cartesian coordinates.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCartesian
    implicit none
    type (coordinateCartesian)                :: position
    class(keplerOrbit        ), intent(inout) :: orbit

    position=+      orbit%radius()                    &
         &   *[                                       &
         &     +sin(orbit%theta ())*cos(orbit%phi()), &
         &     +sin(orbit%theta ())*sin(orbit%phi()), &
         &     +cos(orbit%theta ())                   &
         &    ]
    return
  end function Kepler_Orbits_Position

  function Kepler_Orbits_Velocity(orbit) result(velocity)
    !!{
    Return the position of the orbit in Cartesian coordinates.
    !!}
    use :: Coordinates, only : assignment(=) , coordinateCartesian
    use :: Vectors    , only : Vector_Product
    implicit none
    type            (coordinateCartesian)                :: velocity
    class           (keplerOrbit        ), intent(inout) :: orbit
    double precision                     , dimension(3)  :: radialVector             , velocityRadialVector     , &
         &                                                  velocityTangentialVector1, velocityTangentialVector2

    radialVector=[                                       &
         &        +sin(orbit%theta ())*cos(orbit%phi()), &
         &        +sin(orbit%theta ())*sin(orbit%phi()), &
         &        +cos(orbit%theta ())                   &
         &       ]
    velocityRadialVector     =orbit%velocityRadial()*radialVector
    velocityTangentialVector1=Vector_Product(radialVector,[1.0d0,0.0d0,0.0d0]      )
    velocityTangentialVector1=velocityTangentialVector1/sqrt(sum(velocityTangentialVector1**2))
    velocityTangentialVector2=Vector_Product(radialVector,velocityTangentialVector1)
    velocityTangentialVector1=velocityTangentialVector1*orbit%velocityTangential()*cos(orbit%epsilon())
    velocityTangentialVector2=velocityTangentialVector2*orbit%velocityTangential()*sin(orbit%epsilon())
    velocity                 =velocityRadialVector+velocityTangentialVector1+velocityTangentialVector2
    return
  end function Kepler_Orbits_Velocity

  function Kepler_Orbits_Non_Static_Size_Of(self)
    !!{
    Return the size of any non-static components of the object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t   )                :: Kepler_Orbits_Non_Static_Size_Of
    class  (keplerOrbit), intent(in   ) :: self
    !$GLC attributes unused :: self

    Kepler_Orbits_Non_Static_Size_Of=0_c_size_t
    return
  end function Kepler_Orbits_Non_Static_Size_Of

  logical function Kepler_Orbits_Equivalent(orbit1,orbit2) result(equivalent)
    !!{
    Test whether two {\normalfont \ttfamily keplerOrbit} objects are equal.
    !!}
    implicit none
    class(keplerOrbit), intent(in   ) :: orbit1, orbit2

    equivalent= (.not.orbit1%             massesIsSet.or..not.orbit2%           massesIsSet.or.orbit1%               massHostValue == orbit2%               massHostValue) &
         &     .and.                                     &
         &      (.not.orbit1%             massesIsSet.or..not.orbit2%           massesIsSet.or.orbit1%          massSatelliteValue == orbit2%          massSatelliteValue) &
         &     .and.                                     &
         &      (.not.orbit1%   radiusApocenterIsSet.or..not.orbit2%   radiusApocenterIsSet.or.orbit1%        radiusApocenterValue == orbit2%        radiusApocenterValue) &
         &     .and.                                     &
         &      (.not.orbit1%            radiusIsSet.or..not.orbit2%            radiusIsSet.or.orbit1%                 radiusValue == orbit2%                 radiusValue) &
         &     .and.                                     &
         &      (.not.orbit1%  radiusPericenterIsSet.or..not.orbit2%  radiusPericenterIsSet.or.orbit1%       radiusPericenterValue == orbit2%       radiusPericenterValue) &
         &     .and.                                     &
         &      (.not.orbit1%             thetaIsSet.or..not.orbit2%             thetaIsSet.or.orbit1%                  thetaValue == orbit2%                  thetaValue) &
         &     .and.                                     &
         &      (.not.orbit1%               phiIsSet.or..not.orbit2%               phiIsSet.or.orbit1%                    phiValue == orbit2%                    phiValue) &
         &     .and.                                     &
         &      (.not.orbit1%           epsilonIsSet.or..not.orbit2%           epsilonIsSet.or.orbit1%                epsilonValue == orbit2%                epsilonValue) &
         &     .and.                                     &
         &      (.not.orbit1%    velocityRadialIsSet.or..not.orbit2%    velocityRadialIsSet.or.orbit1%         velocityRadialValue == orbit2%         velocityRadialValue) &
         &     .and.                                     &
         &      (.not.orbit1%velocityTangentialIsSet.or..not.orbit2%velocityTangentialIsSet.or.orbit1%     velocityTangentialValue == orbit2%     velocityTangentialValue) &
         &     .and.                                     &
         &      (.not.orbit1%   angularMomentumIsSet.or..not.orbit2%   angularMomentumIsSet.or.orbit1%        angularMomentumValue == orbit2%        angularMomentumValue) &
         &     .and.                                     &
         &      (.not.orbit1%            energyIsSet.or..not.orbit2%            energyIsSet.or.orbit1%                 energyValue == orbit2%                 energyValue) &
         &     .and.                                     &
         &      (.not.orbit1%      eccentricityIsSet.or..not.orbit2%      eccentricityIsSet.or.orbit1%           eccentricityValue == orbit2%           eccentricityValue) &
         &     .and.                                     &
         &      (.not.orbit1%     semimajorAxisIsSet.or..not.orbit2%     semimajorAxisIsSet.or.orbit1%          semimajorAxisValue == orbit2%          semimajorAxisValue)
    return
  end function Kepler_Orbits_Equivalent

end module Kepler_Orbits

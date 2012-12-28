!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which defines an orbit structure for use in \glc.

module Kepler_Orbits
  !% Defines an orbit structure for use in \glc.
  implicit none
  private
  public :: keplerOrbit

  ! Effective infinite radius used for apocenters in unbound orbits.
  double precision, parameter :: radiusEffectiveInfinity=1.0d30

  type keplerOrbit
     !% The structure used for describing orbits in \glc. This object will automatically convert from one set of orbital
     !% parameters to another where possible. The orbitting bodies (a satellite orbitting around its host) are treated as point
     !% masses, and the usual ``reduced mass'' framework is used, such that radii and velocities are measured relative to a
     !% stationary host. Energy and angular momentum are defined per unit satellite mass (not per unit reduced mass). Note that
     !% not all interconversions between elements are implemented. The object works by attempting to get the radial and tangential
     !% velocities and the radius. If it can obtain these, any other parameter can be computed. Getting these three parameters
     !% relies on having known conversions from other possible combinations of parameters. 
     private
     double precision :: specificReducedMassValue,hostMassValue
     double precision :: radiusValue,radiusPericenterValue,radiusApocenterValue
     double precision :: velocityRadialValue,velocityTangentialValue
     double precision :: angularMomentumValue
     double precision :: energyValue
     double precision :: eccentricityValue
     double precision :: semimajorAxisValue
     logical          :: massesIsSet
     logical          :: radiusIsSet,radiusPericenterIsSet,radiusApocenterIsSet
     logical          :: velocityRadialIsSet,velocityTangentialIsSet
     logical          :: angularMomentumIsSet
     logical          :: energyIsSet
     logical          :: eccentricityIsSet
     logical          :: semimajorAxisIsSet
   contains
     ! Orbit methods.
     !@ <objectMethods>
     !@   <object>orbit</object>
     !@   <objectMethod>
     !@     <method>dump</method>
     !@     <description>Dump an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dumpRaw</method>
     !@     <description>Dump an orbit in binary.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>readRaw</method>
     !@     <description>Read an orbit in binary.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>rese</method>
     !@     <description>Resets an orbit to a null state.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <description>Destroys an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isDefined</method>
     !@     <description>Returns true if an orbit is fully defined.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>assertIsDefined</method>
     !@     <description>Asserts that an orbit is fully defined.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isBound</method>
     !@     <description>Returns true if the orbit is bound.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>propagte</method>
     !@     <description>Propagates an orbit to a new position.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>velocityRadialSet</method>
     !@     <description>Sets the radial velocity of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>massesSet</method>
     !@     <description>Sets the masses of satellite and host objects.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusSet</method>
     !@     <description>Sets the radius of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusPericenterSet</method>
     !@     <description>Sets the pericenter radius of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusApocenterSet</method>
     !@     <description>Sets the apocenter radius of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>velocityTangentialSet</method>
     !@     <description>Sets the tangential velocity of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>energySet</method>
     !@     <description>Sets the energy of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>eccentricitySet</method>
     !@     <description>Sets the eccentricity of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>angularMomentumSet</method>
     !@     <description>Sets the angular momentum of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>semiMajorAxisSet</method>
     !@     <description>Sets the semi-major axis of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>hostMass</method>
     !@     <description>Returns the host mass of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>velocityScale</method>
     !@     <description>Returns the velocity scale of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reducedMassSpecific</method>
     !@     <description>Returns the specific reduced mass (i.e. the reduced mass per unit satellite mass, $\mu_{\rm s} = M_{\rm host}/(M_{\rm satellite}+M_{\rm host})$) of the orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radius</method>
     !@     <description>Returns the radius of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusPericenter</method>
     !@     <description>Returns the pericenter radius of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusApocenter</method>
     !@     <description>Returns the apocenter radius of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>velocityRadial</method>
     !@     <description>Returns the radial velocity of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>velocityTangential</method>
     !@     <description>Returns the tangential velocity of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>energy</method>
     !@     <description>Returns the energy of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>eccentricity</method>
     !@     <description>Returns the eccentricity of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>angularMomentum</method>
     !@     <description>Returns the angular momentum of an orbit.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>semiMajorAxis</method>
     !@     <description>Returns the semi-major axis of an orbit.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: dump                  => Kepler_Orbits_Dump
     procedure :: dumpRaw               => Kepler_Orbits_Dump_Raw
     procedure :: readRaw               => Kepler_Orbits_Read_Raw
     procedure :: reset                 => Kepler_Orbits_Reset
     procedure :: destroy               => Kepler_Orbits_Destroy
     procedure :: isDefined             => Kepler_Orbits_Is_Defined
     procedure :: assertIsDefined       => Kepler_Orbits_Assert_Is_Defined
     procedure :: isBound               => Kepler_Orbits_Is_Bound
     procedure :: propagate             => Kepler_Orbits_Propagate
     procedure :: massesSet             => Kepler_Orbits_Masses_Set
     procedure :: radiusSet             => Kepler_Orbits_Radius_Set
     procedure :: radiusPericenterSet   => Kepler_Orbits_Pericenter_Radius_Set
     procedure :: radiusApocenterSet    => Kepler_Orbits_Apocenter_Radius_Set
     procedure :: velocityRadialSet     => Kepler_Orbits_Velocity_Radial_Set
     procedure :: velocityTangentialSet => Kepler_Orbits_Velocity_Tangential_Set
     procedure :: energySet             => Kepler_Orbits_Energy_Set
     procedure :: angularMomentumSet    => Kepler_Orbits_Angular_Momentum_Set
     procedure :: eccentricitySet       => Kepler_Orbits_Eccentricity_Set
     procedure :: semiMajorAxisSet      => Kepler_Orbits_Semi_Major_Axis_Set
     procedure :: specificReducedMass   => Kepler_Orbits_Specific_Reduced_Mass
     procedure :: hostMass              => Kepler_Orbits_Host_Mass
     procedure :: velocityScale         => Kepler_Orbits_Velocity_Scale
     procedure :: radius                => Kepler_Orbits_Radius
     procedure :: radiusPericenter      => Kepler_Orbits_Pericenter_Radius
     procedure :: radiusApocenter       => Kepler_Orbits_Apocenter_Radius
     procedure :: velocityRadial        => Kepler_Orbits_Velocity_Radial
     procedure :: velocityTangential    => Kepler_Orbits_Velocity_Tangential
     procedure :: energy                => Kepler_Orbits_Energy
     procedure :: angularMomentum       => Kepler_Orbits_Angular_Momentum
     procedure :: eccentricity          => Kepler_Orbits_Eccentricity
     procedure :: semiMajorAxis         => Kepler_Orbits_Semi_Major_Axis
  end type keplerOrbit

contains

  subroutine Kepler_Orbits_Destroy(thisOrbit)
    !% Destroy an orbit.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    ! Nothing to do.
    return
  end subroutine Kepler_Orbits_Destroy

  subroutine Kepler_Orbits_Dump(self)
    !% Reset an orbit to a null state.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class    (keplerOrbit   ), intent(in   ) :: self
    character(len=12        )                :: label
    type     (varying_string)                :: message

    if (self%massesIsSet             ) then
       write (label,'(e12.6)') self%hostMassValue
       message='host mass:             '//label
       call Galacticus_Display_Message(message)
       write (label,'(e12.6)') self%specificReducedMassValue
       message='specific reduced mass: '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%radiusIsSet             ) then
       write (label,'(e12.6)') self%radiusValue
       message='radius            :     '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%radiusPericenterIsSet   ) then
       write (label,'(e12.6)') self%radiusPericenterValue
       message='radius pericenter:      '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%radiusApocenterIsSet   ) then
       write (label,'(e12.6)') self%radiusApocenterValue
       message='radius apocenter:       '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%velocityRadialIsSet    ) then
       write (label,'(e12.6)') self%velocityRadialValue
       message='velocity radial:        '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%velocityTangentialIsSet) then
       write (label,'(e12.6)') self%velocityTangentialValue
       message='velocity tangential:    '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%angularMomentumIsSet   ) then
       write (label,'(e12.6)') self%angularMomentumValue
       message='angular momentum:       '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%energyIsSet            ) then
       write (label,'(e12.6)') self%energyValue
       message='energy:                 '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%eccentricityIsSet      ) then
       write (label,'(e12.6)') self%eccentricityValue
       message='eccentricity:           '//label
       call Galacticus_Display_Message(message)
    end if
    if (self%semimajorAxisIsSet     ) then
       write (label,'(e12.6)') self%semimajorAxisValue
       message='semi-major axis:        '//label
       call Galacticus_Display_Message(message)
    end if
    return
  end subroutine Kepler_Orbits_Dump

  subroutine Kepler_Orbits_Dump_Raw(self,fileHandle)
    !% Dump a {\tt keplerOrbit} object in binary.
    implicit none
    class  (keplerOrbit), intent(in   ) :: self
    integer             , intent(in   ) :: fileHandle 

    write (fileHandle) self%massesIsSet,self%massesIsSet,self%radiusIsSet,self%radiusPericenterIsSet,self%radiusApocenterIsSet&
         &,self%velocityRadialIsSet,self%velocityTangentialIsSet,self%angularMomentumIsSet,self%energyIsSet&
         &,self%eccentricityIsSet,self%semimajorAxisIsSet
    if (self%massesIsSet            ) write (fileHandle) self%hostMassValue,self%specificReducedMassValue
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
    !% Read a {\tt keplerOrbit} object in binary.
    implicit none
    class  (keplerOrbit), intent(inout) :: self
    integer             , intent(in   ) :: fileHandle 

    read (fileHandle) self%massesIsSet,self%massesIsSet,self%radiusIsSet,self%radiusPericenterIsSet,self%radiusApocenterIsSet&
         &,self%velocityRadialIsSet,self%velocityTangentialIsSet,self%angularMomentumIsSet,self%energyIsSet&
         &,self%eccentricityIsSet,self%semimajorAxisIsSet
    if (self%massesIsSet            ) read (fileHandle) self%hostMassValue,self%specificReducedMassValue
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

  subroutine Kepler_Orbits_Reset(thisOrbit)
    !% Reset an orbit to a null state.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Simply specify that no properties have been set as yet.
       thisOrbit%massesIsSet            =.false.
       thisOrbit%radiusIsSet            =.false.
       thisOrbit%radiusPericenterIsSet  =.false.
       thisOrbit%radiusApocenterIsSet   =.false.
       thisOrbit%velocityRadialIsSet    =.false.
       thisOrbit%velocityTangentialIsSet=.false.
       thisOrbit%angularMomentumIsSet   =.false.
       thisOrbit%energyIsSet            =.false.
       thisOrbit%eccentricityIsSet      =.false.
       thisOrbit%semimajorAxisIsSet     =.false.
    end select
    return
  end subroutine Kepler_Orbits_Reset

  subroutine Kepler_Orbits_Masses_Set(thisOrbit,satelliteMass,hostMass)
    !% Sets the masses of the two orbitting objects in a {\tt keplerOrbit} object.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: satelliteMass,hostMass
    
    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the mass factor and flag that is set.
       thisOrbit%specificReducedMassValue=1.0d0/(1.0d0+satelliteMass/hostMass)
       thisOrbit%hostMassValue           =hostMass
       thisOrbit%massesIsSet             =.true.
    end select
    return
  end subroutine Kepler_Orbits_Masses_Set

  subroutine Kepler_Orbits_Radius_Set(thisOrbit,radius)
    !% Sets the radius to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: radius

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the radius and flag that is set.
       thisOrbit%radiusValue=radius
       thisOrbit%radiusIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Radius_Set

  subroutine Kepler_Orbits_Pericenter_Radius_Set(thisOrbit,radius)
    !% Sets the pericenter radius to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: radius

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the pericenter radius and flag that is set.
       thisOrbit%radiusPericenterValue=radius
       thisOrbit%radiusPericenterIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Pericenter_Radius_Set

  subroutine Kepler_Orbits_Apocenter_Radius_Set(thisOrbit,radius)
    !% Sets the apocenter radius to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: radius

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the apocenter radius and flag that is set.
       thisOrbit%radiusApocenterValue=radius
       thisOrbit%radiusApocenterIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Apocenter_Radius_Set

  subroutine Kepler_Orbits_Velocity_Radial_Set(thisOrbit,velocityRadial)
    !% Sets the radial velocity to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: velocityRadial

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the radial velocity and flag that is set.
       thisOrbit%velocityRadialValue=velocityRadial
       thisOrbit%velocityRadialIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Velocity_Radial_Set

  subroutine Kepler_Orbits_Velocity_Tangential_Set(thisOrbit,velocityTangential)
    !% Sets the tangential velocity to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: velocityTangential

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the tangential velocity and flag that is set.
       thisOrbit%velocityTangentialValue=velocityTangential
       thisOrbit%velocityTangentialIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Velocity_Tangential_Set

  subroutine Kepler_Orbits_Energy_Set(thisOrbit,energy)
    !% Sets the tangential velocity to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: energy

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the tangential velocity and flag that is set.
       thisOrbit%energyValue=energy
       thisOrbit%energyIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Energy_Set

  subroutine Kepler_Orbits_Eccentricity_Set(thisOrbit,eccentricity)
    !% Sets the tangential velocity to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: eccentricity

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the tangential velocity and flag that is set.
       thisOrbit%eccentricityValue=eccentricity
       thisOrbit%eccentricityIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Eccentricity_Set

  subroutine Kepler_Orbits_Angular_Momentum_Set(thisOrbit,angularMomentum)
    !% Sets the tangential velocity to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: angularMomentum

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the tangential velocity and flag that is set.
       thisOrbit%angularMomentumValue=angularMomentum
       thisOrbit%angularMomentumIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Angular_Momentum_Set

  subroutine Kepler_Orbits_Semi_Major_Axis_Set(thisOrbit,semiMajorAxis)
    !% Sets the semi-major axis to the specified value.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision,   intent(in)    :: semiMajorAxis

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Set the tangential velocity and flag that is set.
       thisOrbit%semiMajorAxisValue=semiMajorAxis
       thisOrbit%semiMajorAxisIsSet=.true.
    end select
    return
  end subroutine Kepler_Orbits_Semi_Major_Axis_Set

  double precision function Kepler_Orbits_Specific_Reduced_Mass(thisOrbit)
    !% Return the specific reduced mass for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%massesIsSet) call Galacticus_Error_Report('Kepler_Orbits_Specific_Reduced_Mass', &
            & 'mass factor has not been set for this orbit')
       Kepler_Orbits_Specific_Reduced_Mass=thisOrbit%specificReducedMassValue
    end select
    return
  end function Kepler_Orbits_Specific_Reduced_Mass

  double precision function Kepler_Orbits_Host_Mass(thisOrbit)
    !% Return the host mass for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%massesIsSet) call Galacticus_Error_Report('Kepler_Orbits_Host_Mass', &
            & 'host mass has not been set for this orbit')
       Kepler_Orbits_Host_Mass=thisOrbit%hostMassValue
    end select
    return
  end function Kepler_Orbits_Host_Mass

  double precision function Kepler_Orbits_Radius(thisOrbit)
    !% Return the radius for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%radiusIsSet) call Galacticus_Error_Report('Kepler_Orbits_Radius', &
            & 'radius has not been set for this orbit')
       Kepler_Orbits_Radius=thisOrbit%radiusValue
    end select
    return
  end function Kepler_Orbits_Radius

  double precision function Kepler_Orbits_Pericenter_Radius(thisOrbit)
    !% Return the pericenter radius for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%radiusPericenterIsSet) then
          ! Assert that the orbit is defined.
          call thisOrbit%assertIsDefined()
          ! Compute the pericenter radius.
          thisOrbit%radiusPericenterValue=thisOrbit%semiMajorAxis()*(1.0d0-thisOrbit%eccentricity())
          thisOrbit%radiusPericenterIsSet=.true.
       end if
       Kepler_Orbits_Pericenter_Radius=thisOrbit%radiusPericenterValue
    end select
    return
  end function Kepler_Orbits_Pericenter_Radius

  double precision function Kepler_Orbits_Apocenter_Radius(thisOrbit)
    !% Return the apocenter radius for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%radiusApocenterIsSet) then
          ! Assert that the orbit is defined.
          call thisOrbit%assertIsDefined()
          ! Compute the pericenter radius.
          if (thisOrbit%isBound()) then
             thisOrbit%radiusApocenterValue=thisOrbit%semiMajorAxis()*(1.0d0+thisOrbit%eccentricity())
          else
             thisOrbit%radiusApocenterValue=radiusEffectiveInfinity
          end if
          thisOrbit%radiusApocenterIsSet=.true.
       end if
       Kepler_Orbits_Apocenter_Radius=thisOrbit%radiusApocenterValue
    end select
    return
  end function Kepler_Orbits_Apocenter_Radius

  double precision function Kepler_Orbits_Velocity_Radial(thisOrbit)
    !% Return the radial velocity for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%velocityRadialIsSet) then
          ! Assert that the orbit is defined.
          call thisOrbit%assertIsDefined()
          ! Compute from eccentricity, radius and periapsis if possible.
          if (thisOrbit%eccentricityIsSet.and.thisOrbit%radiusIsSet.and.thisOrbit%radiusPericenterIsSet) then
             thisOrbit%velocityRadialValue=thisOrbit%velocityScale()*dsqrt((2.0d0*(1.0d0-thisOrbit%radius()&
                  &/thisOrbit%radiusPericenter())+(1.0d0+thisOrbit%eccentricity())*(thisOrbit%radius()/thisOrbit%radiusPericenter()&
                  &-thisOrbit%radiusPericenter()/thisOrbit%radius()))/thisOrbit%specificReducedMass())
             thisOrbit%velocityRadialIsSet=.true.
          end if
          ! If we were not able to compute the radial velocity, exit.
          if (.not.thisOrbit%velocityRadialIsSet) call Galacticus_Error_Report('Kepler_Orbits_Velocity_Radial','radial velocity&
               & has not been set for this orbit and can not be computed')
       end if
       Kepler_Orbits_Velocity_Radial=thisOrbit%velocityRadialValue
    end select
    return
  end function Kepler_Orbits_Velocity_Radial

  double precision function Kepler_Orbits_Velocity_Tangential(thisOrbit)
    !% Return the tangential velocity for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%velocityTangentialIsSet) then
          ! Assert that the orbit is defined.
          call thisOrbit%assertIsDefined()
          ! Compute from angular momentum and radius if possible.
          if (thisOrbit%angularMomentumIsSet.and.thisOrbit%radiusIsSet) then
             thisOrbit%velocityTangentialValue=thisOrbit%angularMomentum()/thisOrbit%radius()/thisOrbit%specificReducedMass()
             thisOrbit%velocityTangentialIsSet=.true.
          end if
          ! Compute from eccentricity, radius and periapsis if possible.
          if (thisOrbit%eccentricityIsSet.and.thisOrbit%radiusIsSet.and.thisOrbit%radiusPericenterIsSet) then
             thisOrbit%velocityTangentialValue=thisOrbit%velocityScale()*dsqrt((1.0d0+thisOrbit%eccentricity())&
                  &*thisOrbit%radiusPericenter() /thisOrbit%radius() /thisOrbit%specificReducedMass())
             thisOrbit%velocityTangentialIsSet=.true.
          end if
          ! If we were not able to compute the tangential velocity, exit.
          if (.not.thisOrbit%velocityTangentialIsSet) call Galacticus_Error_Report('Kepler_Orbits_Velocity_Tangential','tangential velocity&
               & has not been set for this orbit and can not be computed')
       end if
       Kepler_Orbits_Velocity_Tangential=thisOrbit%velocityTangentialValue
    end select
    return
  end function Kepler_Orbits_Velocity_Tangential

  double precision function Kepler_Orbits_Energy(thisOrbit)
    !% Return the energy for this orbit.
    use Galacticus_Error
    use Numerical_Constants_Physical
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Check if energy is set.
       if (.not.thisOrbit%energyIsSet) then
          ! Assert that the orbit is defined.
          call thisOrbit%assertIsDefined()
          ! Compute the energy.
          thisOrbit%energyValue=-gravitationalConstantGalacticus*thisOrbit%hostMass()/thisOrbit%radius()+0.5d0&
               &*(thisOrbit%velocityRadial()**2+thisOrbit%velocityTangential()**2)*thisOrbit%specificReducedMass()
          thisOrbit%energyIsSet=.true.
       end if
       
       Kepler_Orbits_Energy=thisOrbit%energyValue
    end select
    return
  end function Kepler_Orbits_Energy

  double precision function Kepler_Orbits_Angular_Momentum(thisOrbit)
    !% Return the angular momentum for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%angularMomentumIsSet) then
          thisOrbit%angularMomentumValue=thisOrbit%radius()*thisOrbit%velocityTangential()*thisOrbit%specificReducedMass()
          thisOrbit%angularMomentumIsSet=.true.
       end if
       
       Kepler_Orbits_Angular_Momentum=thisOrbit%angularMomentumValue
    end select
    return
  end function Kepler_Orbits_Angular_Momentum

  double precision function Kepler_Orbits_Eccentricity(thisOrbit)
    !% Return the eccentricity for this orbit.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision                 :: velocityTangential,velocityRadial

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%eccentricityIsSet) then
          velocityTangential         =thisOrbit%velocityTangential()/thisOrbit%velocityScale()
          velocityRadial             =thisOrbit%velocityRadial    ()/thisOrbit%velocityScale()
          thisOrbit%eccentricityValue=dsqrt(1.0d0+2.0d0*thisOrbit%energy()*thisOrbit%angularMomentum()**2/thisOrbit%radius()**2&
               &/thisOrbit%velocityScale()**4/thisOrbit%specificReducedMass())
          thisOrbit%eccentricityIsSet=.true.
       end if
       
       Kepler_Orbits_Eccentricity=thisOrbit%eccentricityValue
    end select
    return
  end function Kepler_Orbits_Eccentricity

  double precision function Kepler_Orbits_Semi_Major_Axis(thisOrbit)
    !% Return the semi-major axis for this orbit.
    implicit none
    class(keplerOrbit), intent(inout) :: thisOrbit
    double precision                 :: velocityTangential,velocityRadial

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%semiMajorAxisIsSet) then
          velocityTangential          =thisOrbit%velocityTangential()/thisOrbit%velocityScale()
          velocityRadial              =thisOrbit%velocityRadial    ()/thisOrbit%velocityScale()
          thisOrbit%semiMajorAxisValue=thisOrbit%radius()/thisOrbit%specificReducedMass()/(2.0d0/thisOrbit%specificReducedMass()&
               &-velocityRadial**2-velocityTangential**2)
          thisOrbit%semiMajorAxisIsSet=.true.
       end if
       
       Kepler_Orbits_Semi_Major_Axis=thisOrbit%semiMajorAxisValue
    end select
    return
  end function Kepler_Orbits_Semi_Major_Axis

  logical function Kepler_Orbits_Is_Defined(thisOrbit)
    !% Returns true if the orbit is fully defined. For the orbits consider here, in which we don't care about the orientation of
    !% the orbital plane or the argument of pericenter, this requires that three orbital parameter be set (in addition to the
    !% masses of the orbitting bodies).
    implicit none
    class(keplerOrbit), intent(in) :: thisOrbit
    integer                        :: orbitalParameterCount          

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Assume orbit is not defined by default.
       Kepler_Orbits_Is_Defined=.false. 
       ! Ensure that masses are set.
       if (thisOrbit%massesIsSet) then
          ! Count how many parameters are set.
          orbitalParameterCount=0
          if (thisOrbit%radiusIsSet            ) orbitalParameterCount=orbitalParameterCount+1
          if (thisOrbit%radiusPericenterIsSet  ) orbitalParameterCount=orbitalParameterCount+1
          if (thisOrbit%radiusApocenterIsSet   ) orbitalParameterCount=orbitalParameterCount+1
          if (thisOrbit%energyIsSet            ) orbitalParameterCount=orbitalParameterCount+1
          if (thisOrbit%angularMomentumIsSet   ) orbitalParameterCount=orbitalParameterCount+1
          if (thisOrbit%velocityRadialIsSet    ) orbitalParameterCount=orbitalParameterCount+1
          if (thisOrbit%velocityTangentialIsSet) orbitalParameterCount=orbitalParameterCount+1
          if (thisOrbit%eccentricityIsSet      ) orbitalParameterCount=orbitalParameterCount+1
          if (thisOrbit%semiMajorAxisIsSet     ) orbitalParameterCount=orbitalParameterCount+1
          ! Orbit is defined if at least 3 parameters are set.
          Kepler_Orbits_Is_Defined=(orbitalParameterCount >= 3)
       end if
    end select
    return
  end function Kepler_Orbits_Is_Defined

  subroutine Kepler_Orbits_Assert_Is_Defined(thisOrbit)
    !% Assert that an orbit is defined - quit with an error if it is not.
    use Galacticus_Error
    implicit none
    class(keplerOrbit), intent(in) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       if (.not.thisOrbit%isDefined()) call Galacticus_Error_Report('Kepler_Orbits_Assert_Is_Defined','orbit is not defined')
    end select
    return
  end subroutine Kepler_Orbits_Assert_Is_Defined

  logical function Kepler_Orbits_Is_Bound(thisOrbit)
    !% Returns true if the orbit is bound.
    implicit none
    class(keplerOrbit), intent(in) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Assert that the orbit is defined.
       call thisOrbit%assertIsDefined()
       ! Test if the energy is negative (indicating a bound orbit).
       Kepler_Orbits_Is_Bound=(thisOrbit%energy() < 0.0d0)
    end select
    return
  end function Kepler_Orbits_Is_Bound

  double precision function Kepler_Orbits_Velocity_Scale(thisOrbit)
    !% Return the velocity scale for the orbit.
    use Galacticus_Error
    use Numerical_Constants_Physical
    implicit none
    class(keplerOrbit), intent(in) :: thisOrbit

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Check that masses and radius have been specified.
       if (.not.(thisOrbit%radiusIsSet.and.thisOrbit%massesIsSet)) call Galacticus_Error_Report('Kepler_Orbits_Velocity_Scale','orbit&
            & masses and radius must be specified')
       ! Compute the velocity scale.
       Kepler_Orbits_Velocity_Scale=dsqrt(gravitationalConstantGalacticus*thisOrbit%hostMass()/thisOrbit%radius())
    end select
    return
  end function Kepler_Orbits_Velocity_Scale

  subroutine Kepler_Orbits_Propagate(thisOrbit,newRadius,infalling)
    !% Propagate an orbit along its path.
    use Galacticus_Error
    use Numerical_Constants_Physical
    implicit none
    class(keplerOrbit), intent(inout)           :: thisOrbit
    double precision,   intent(in   )           :: newRadius
    logical,            intent(in   ), optional :: infalling
    double precision                            :: energy,angularMomentum,newVelocityRadial,newVelocityTangential

    select type(thisOrbit)
    type is (keplerOrbit)
       ! Assert that the orbit is defined.
       call thisOrbit%assertIsDefined()
       ! Check that the radius is within the allowed range.
       if     (                                                                                                               &
            &    newRadius < thisOrbit%radiusPericenter()                                                                     &
            &  .or.                                                                                                           &
            &   (newRadius > thisOrbit%radiusApocenter () .and. thisOrbit%radiusApocenter() > 0.0d0)                          &
            & ) call Galacticus_Error_Report('Kepler_Orbits_Propagate','radius lies outside of allowed range for this orbit')
       ! Get the energy and angular momentum
       energy         =thisOrbit%energy         ()
       angularMomentum=thisOrbit%angularMomentum()
       ! Compute velocity components.
       newVelocityTangential=angularMomentum/newRadius
       newVelocityRadial    =sqrt(2.0d0*(energy+gravitationalConstantGalacticus*thisOrbit%hostMass()/newRadius)/thisOrbit%specificReducedMass()-newVelocityTangential**2)
       ! Move to the infalling phase of the orbit if requested.
       if (present(infalling)) then
          if (infalling) newVelocityRadial=-newVelocityRadial
       end if
       ! Set new values for the state vector.
       call thisOrbit%radiusSet            (newRadius            )
       call thisOrbit%velocityTangentialSet(newVelocityTangential)
       call thisOrbit%velocityRadialSet    (newVelocityRadial    )
    end select
    return
  end subroutine Kepler_Orbits_Propagate

end module Kepler_Orbits

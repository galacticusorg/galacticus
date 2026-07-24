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
Contains a module which implements the coordinates class.
!!}

module Coordinates
  !!{RST
  Implements the coordinates class.

  .. _coordinates-array-convention:

  **Convention for bare 3-element arrays.** Not every position or velocity in Galacticus is carried as a
  ``coordinate`` object: bulk and serialized data --- node-component ``position``/``velocity`` properties,
  lightcone positions and velocities, virial-orbit vector means, N-body particle blocks, and merger-tree
  reader records --- are stored and passed as plain ``double precision, dimension(3)`` arrays, for which
  an object per particle would be a regression. **Wherever such a bare 3-element array represents a
  position, velocity, or angular momentum, its components are Cartesian, in the order** :math:`(x,y,z)`.
  The coordinate system is a convention of the representation, not a property carried by the type.

  Code which needs a system-aware position should assign the array into a ``coordinateCartesian`` object,
  after which assignment to a ``coordinateSpherical`` or ``coordinateCylindrical`` object converts
  automatically. Note that the reverse assignment (array ``=`` coordinate) yields the *raw* component triple
  of whatever system the object is in --- ``(r,theta,phi)`` for a spherical object, not ``(x,y,z)`` --- so
  prefer ``%toCartesian()`` when a Cartesian array is what is wanted.
  !!}
  implicit none
  private
  public :: assignment(=), operator(*)

  ! Define assignment interfaces.
  interface assignment(=)
     module procedure Coordinates_Assign_To
     module procedure Coordinates_Assign_From
     module procedure Coordinates_Assign
  end interface assignment(=)

  type, abstract, public :: coordinate
     !!{RST
     The base coordinate object class.
     !!}
     ! The raw component triple. This is `private` because its meaning is specific to the coordinate system of
     ! the extending type ((x,y,z), (r,theta,phi), or (r,phi,z)), so direct access outside this module couples
     ! callers to that per-system interpretation (see issue \#75). Use the named accessors of the concrete
     ! types (`%x()`, `%r()`, `%theta()`, ...) or `%toCartesian()` instead. The three concrete coordinate
     ! types are all defined in this module, so they retain access.
     double precision, private :: position(3)
   contains
     !![
     <methods docformat="rst">
       <method description="Return the coordinates in a Cartesian system as a 3-element array."          method="toCartesian"      />
       <method description="Set the coordinates from a Cartesian system specified as a 3-element array." method="fromCartesian"    />
       <method description="Return the cylindrical radial coordinate."                                   method="rCylindrical"     />
       <method description="Return the spherical radial coordinate."                                     method="rSpherical"       />
       <method description="Return the square of the spherical radial coordinate."                       method="rSphericalSquared"/>
       <method description="Multiply the coordinate by a scalar."                                        method="operator(*)"      />
       <method description="Divide the coordinate by a scalar."                                          method="operator(/)"      />
       <method description="Multiply the coordinate by a scalar."                                        method="scalarMultiply"   />
       <method description="Divide the coordinate by a scalar."                                          method="scalarDivide"     />
       <method description="Scale the coordinates by a scalar."                                          method="scale"            />
     </methods>
     !!]
     procedure                                      :: toCartesian       => Coordinates_Null_To
     procedure                                      :: fromCartesian     => Coordinates_Null_From
     procedure                                      :: rCylindrical      => Coordinates_Radius_Cylindrical
     procedure                                      :: rSpherical        => Coordinates_Radius_Spherical
     procedure(rSphericalSquaredTemplate), deferred :: rSphericalSquared
     procedure(scalarMultiplyTemplate   ), deferred :: scalarMultiply
     procedure(scalarDivideTemplate     ), deferred :: scalarDivide
     generic                                        :: operator(*)       => scalarMultiply
     generic                                        :: operator(/)       => scalarDivide
     procedure(scaleTemplate            ), deferred :: scale
  end type coordinate

  type, public, extends(coordinate) :: coordinateCartesian
     !!{RST
     A Cartesian coordinate object class.
     !!}
   contains
     !![
     <methods docformat="rst">
       <method description="Get the :math:`x`-coordinate." method="x"   />
       <method description="Get the :math:`y`-coordinate." method="y"   />
       <method description="Get the :math:`z`-coordinate." method="z"   />
       <method description="set the :math:`x`-coordinate." method="xSet"/>
       <method description="set the :math:`y`-coordinate." method="ySet"/>
       <method description="set the :math:`z`-coordinate." method="zSet"/>
     </methods>
     !!]
     procedure :: toCartesian       => Coordinates_Cartesian_To_Cartesian
     procedure :: fromCartesian     => Coordinates_Cartesian_From_Cartesian
     procedure :: x                 => Coordinates_Cartesian_X
     procedure :: y                 => Coordinates_Cartesian_Y
     procedure :: z                 => Coordinates_Cartesian_Z
     procedure :: xSet              => Coordinates_Cartesian_Set_X
     procedure :: ySet              => Coordinates_Cartesian_Set_Y
     procedure :: zSet              => Coordinates_Cartesian_Set_Z
     procedure :: rSphericalSquared => Coordinates_Cartesian_R_Spherical_Squared
     procedure :: scalarMultiply    => Coordinates_Cartesian_Scalar_Multiply
     procedure :: scalarDivide      => Coordinates_Cartesian_Scalar_Divide
     procedure :: scale             => Coordinates_Cartesian_Scale
  end type coordinateCartesian

  type, public, extends(coordinate) :: coordinateSpherical
     !!{RST
     A spherical coordinate object class.
     !!}
   contains
     !![
     <methods docformat="rst">
       <method description="Get the :math:`r`-coordinate."      method="r"       />
       <method description="Get the :math:`\theta`-coordinate." method="theta"   />
       <method description="Get the :math:`\phi`-coordinate."   method="phi"     />
       <method description="set the :math:`r`-coordinate."      method="rSet"    />
       <method description="set the :math:`\theta`-coordinate." method="thetaSet"/>
       <method description="set the :math:`\phi`-coordinate."   method="phiSet"  />
     </methods>
     !!]
     procedure :: toCartesian       => Coordinates_Spherical_To_Cartesian
     procedure :: fromCartesian     => Coordinates_Spherical_From_Cartesian
     procedure :: r                 => Coordinates_Spherical_R
     procedure :: theta             => Coordinates_Spherical_Theta
     procedure :: phi               => Coordinates_Spherical_Phi
     procedure :: rSet              => Coordinates_Spherical_Set_R
     procedure :: thetaSet          => Coordinates_Spherical_Set_Theta
     procedure :: phiSet            => Coordinates_Spherical_Set_Phi
     procedure :: rSpherical        => Coordinates_Spherical_R_Spherical
     procedure :: rSphericalSquared => Coordinates_Spherical_R_Spherical_Squared
     procedure :: scalarMultiply    => Coordinates_Spherical_Scalar_Multiply
     procedure :: scalarDivide      => Coordinates_Spherical_Scalar_Divide
     procedure :: scale             => Coordinates_Spherical_Scale
  end type coordinateSpherical

  type, public, extends(coordinate) :: coordinateCylindrical
     !!{RST
     A cylindrical coordinate object class.
     !!}
   contains
     !![
     <methods docformat="rst">
       <method description="Get the :math:`r`-coordinate."    method="r"     />
       <method description="Get the :math:`\phi`-coordinate." method="phi"   />
       <method description="Get the :math:`z`-coordinate."    method="z"     />
       <method description="set the :math:`r`-coordinate."    method="rSet"  />
       <method description="set the :math:`\phi`-coordinate." method="phiSet"/>
       <method description="set the :math:`z`-coordinate."    method="zSet"  />
     </methods>
     !!]
     procedure :: toCartesian       => Coordinates_Cylindrical_To_Cartesian
     procedure :: fromCartesian     => Coordinates_Cylindrical_From_Cartesian
     procedure :: r                 => Coordinates_Cylindrical_R
     procedure :: phi               => Coordinates_Cylindrical_Phi
     procedure :: z                 => Coordinates_Cylindrical_Z
     procedure :: rSet              => Coordinates_Cylindrical_Set_R
     procedure :: phiSet            => Coordinates_Cylindrical_Set_Phi
     procedure :: zSet              => Coordinates_Cylindrical_Set_Z
     procedure :: rSphericalSquared => Coordinates_Cylindrical_R_Spherical_Squared
     procedure :: scalarMultiply    => Coordinates_Cylindrical_Scalar_Multiply
     procedure :: scalarDivide      => Coordinates_Cylindrical_Scalar_Divide
     procedure :: scale             => Coordinates_Cylindrical_Scale
  end type coordinateCylindrical

  ! Explicit constructors for the coordinate types, allowing e.g. `coordinateCartesian(x,y,z)` and
  ! `coordinateSpherical(r,theta,phi)`. These are preferred over array assignment (`coord = [a,b,c]`) in new
  ! code because they make the coordinate system explicit at the construction site (avoiding the footgun where
  ! a bare array is silently interpreted in whatever system the target happens to be -- see issue \#75).
  !
  ! The argument-less `coordinateCartesian()` form is also provided: it is emitted by auto-generated
  ! `functionClass` code for methods that return a `coordinateCartesian` (e.g. the massDistribution
  ! acceleration/positionSample/chandrasekharIntegral). Using an explicit constructor rather than a
  ! default-initializer on the `position` component avoids zero-initializing ordinary coordinate objects on
  ! every construction (which, because coordinate assignment is a defined assignment with an intent(out)
  ! argument, the compiler does not elide -- it adds measurable cost in coordinate-heavy hot paths).
  interface coordinateCartesian
     module procedure coordinatesCartesianConstructorNull
     module procedure coordinatesCartesianConstructor
  end interface coordinateCartesian

  interface coordinateSpherical
     module procedure coordinatesSphericalConstructor
  end interface coordinateSpherical

  interface coordinateCylindrical
     module procedure coordinatesCylindricalConstructor
  end interface coordinateCylindrical

  abstract interface
     double precision function rSphericalSquaredTemplate(self)
       import coordinate
       class(coordinate), intent(in   ) :: self
     end function rSphericalSquaredTemplate
  end interface
  
  abstract interface
     function scalarMultiplyTemplate(self,multiplier)
       import coordinate
       class           (coordinate), allocatable   :: scalarMultiplyTemplate
       class           (coordinate), intent(in   ) :: self
       double precision            , intent(in   ) :: multiplier
     end function scalarMultiplyTemplate
  end interface
  
  abstract interface
     function scalarDivideTemplate(self,divisor)
       import coordinate
       class           (coordinate), allocatable   :: scalarDivideTemplate
       class           (coordinate), intent(in   ) :: self
       double precision            , intent(in   ) :: divisor
     end function scalarDivideTemplate
  end interface
  
  abstract interface
     subroutine scaleTemplate(self,scalar,selfScaled)
       import coordinate
       class           (coordinate), intent(in   )              :: self
       double precision            , intent(in   )              :: scalar
       class           (coordinate), intent(inout), allocatable :: selfScaled
     end subroutine scaleTemplate
  end interface
  
  ! Interface to multiplication operators with coordinate objects as their second argument.
  interface operator(*)
     module procedure Coordinates_Scalar_Multiply_Switched
  end interface operator(*)

contains

  function coordinatesCartesianConstructorNull() result(self)
    !!{RST
    Null constructor for a Cartesian ``coordinate`` object, returning the origin. See the interface block for
    why this exists.
    !!}
    implicit none
    type(coordinateCartesian) :: self

    self%position=0.0d0
    return
  end function coordinatesCartesianConstructorNull

  function coordinatesCartesianConstructor(x,y,z) result(self)
    !!{RST
    Constructor for a Cartesian ``coordinate`` object from its :math:`(x,y,z)` components.
    !!}
    implicit none
    type            (coordinateCartesian)                :: self
    double precision                     , intent(in   ) :: x   , y, z

    self%position=[x,y,z]
    return
  end function coordinatesCartesianConstructor

  function coordinatesSphericalConstructor(r,theta,phi) result(self)
    !!{RST
    Constructor for a spherical ``coordinate`` object from its :math:`(r,\theta,\phi)` components.
    !!}
    implicit none
    type            (coordinateSpherical)                :: self
    double precision                     , intent(in   ) :: r   , theta, phi

    self%position=[r,theta,phi]
    return
  end function coordinatesSphericalConstructor

  function coordinatesCylindricalConstructor(r,phi,z) result(self)
    !!{RST
    Constructor for a cylindrical ``coordinate`` object from its :math:`(r,\phi,z)` components.
    !!}
    implicit none
    type            (coordinateCylindrical)                :: self
    double precision                       , intent(in   ) :: r   , phi, z

    self%position=[r,phi,z]
    return
  end function coordinatesCylindricalConstructor

  subroutine Coordinates_Null_From(self,x)
    !!{RST
    Set generic coordinate object from Cartesian point. Simply quits with an error.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (coordinate)              , intent(  out) :: self
    double precision            , dimension(3), intent(in   ) :: x
    !$GLC attributes unused :: self, x

    call Error_Report('no transformation from cartesian coordinates defined'//{introspection:location})
    return
  end subroutine Coordinates_Null_From

  function Coordinates_Null_To(self)
    !!{RST
    Convert generic coordinate object to Cartesian point. Simply quits with an error.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (coordinate), intent(in   ) :: self
    double precision            , dimension(3)  :: Coordinates_Null_To
    !$GLC attributes unused :: self

    Coordinates_Null_To=0.0d0
    call Error_Report('no transformation to cartesian coordinates defined'//{introspection:location})
    return
  end function Coordinates_Null_To

  subroutine Coordinates_Assign(coordinatesTo,coordinatesFrom)
    !!{RST
    Assign one coordinate object to another, automatically handling the conversion between coordinate systems.
    !!}
    implicit none
    class           (coordinate), intent(  out) :: coordinatesTo
    class           (coordinate), intent(in   ) :: coordinatesFrom
    double precision            , dimension(3)  :: x

    ! Handle special cases of conversion (this is done for optimization).
    if (same_type_as(coordinatesTo,coordinatesFrom)) then
       coordinatesTo%position=coordinatesFrom%position
    else
       ! Assign by transforming through cartesian coordinates.
       x=coordinatesFrom%toCartesian()
       call coordinatesTo%fromCartesian(x)
    end if
    return
  end subroutine Coordinates_Assign

  subroutine Coordinates_Assign_To(coordinates,x)
    !!{RST
    Assign a 3-component vector to a ``coordinate`` object.
    !!}
    implicit none
    class           (coordinate)              , intent(  out) :: coordinates
    double precision            , dimension(3), intent(in   ) :: x

    coordinates%position=x
    return
  end subroutine Coordinates_Assign_To

  subroutine Coordinates_Assign_From(x,coordinates)
    !!{RST
    Return a 3-component vector from a ``coordinate`` object.
    !!}
    implicit none
    class           (coordinate)              , intent(in   ) :: coordinates
    double precision            , dimension(3), intent(  out) :: x

    x=coordinates%position
    return
  end subroutine Coordinates_Assign_From

  ! Cartesian coordinate object.
  subroutine Coordinates_Cartesian_From_Cartesian(self,x)
    !!{RST
    Create a Cartesian ``coordinate`` object from a Cartesian vector.
    !!}
    implicit none
    class           (coordinateCartesian)              , intent(  out) :: self
    double precision                     , dimension(3), intent(in   ) :: x

    self%position=x
    return
  end subroutine Coordinates_Cartesian_From_Cartesian

  function Coordinates_Cartesian_To_Cartesian(self)
    !!{RST
    Return a Cartesian vector from a Cartesian ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateCartesian), intent(in   ) :: self
    double precision                     , dimension(3)  :: Coordinates_Cartesian_To_Cartesian

    Coordinates_Cartesian_To_Cartesian=self%position
    return
  end function Coordinates_Cartesian_To_Cartesian

  double precision function Coordinates_Cartesian_X(self)
    !!{RST
    Return the :math:`x`-component of a Cartesian ``coordinate`` object.
    !!}
    implicit none
    class(coordinateCartesian), intent(in   ) :: self

    Coordinates_Cartesian_X=self%position(1)
    return
  end function Coordinates_Cartesian_X

  double precision function Coordinates_Cartesian_Y(self)
    !!{RST
    Return the :math:`y`-component of a Cartesian ``coordinate`` object.
    !!}
    implicit none
    class(coordinateCartesian), intent(in   ) :: self

    Coordinates_Cartesian_Y=self%position(2)
    return
  end function Coordinates_Cartesian_Y

  double precision function Coordinates_Cartesian_Z(self)
    !!{RST
    Return the :math:`z`-component of a Cartesian ``coordinate`` object.
    !!}
    implicit none
    class(coordinateCartesian), intent(in   ) :: self

    Coordinates_Cartesian_Z=self%position(3)
    return
  end function Coordinates_Cartesian_Z

  subroutine Coordinates_Cartesian_Set_X(self,x)
    !!{RST
    Return the :math:`x`-component of a Cartesian ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: x

    self%position(1)=x
    return
  end subroutine Coordinates_Cartesian_Set_X

  subroutine Coordinates_Cartesian_Set_Y(self,y)
    !!{RST
    Return the :math:`y`-component of a Cartesian ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: y

    self%position(2)=y
    return
  end subroutine Coordinates_Cartesian_Set_Y

  subroutine Coordinates_Cartesian_Set_Z(self,z)
    !!{RST
    Return the :math:`z`-component of a Cartesian ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: z

    self%position(3)=z
    return
  end subroutine Coordinates_Cartesian_Set_Z

  function Coordinates_Cartesian_Scalar_Multiply(self,multiplier) result(scaled)
    !!{RST
    Multiply a Cartesian ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinate         ), allocatable   :: scaled
    class           (coordinateCartesian), intent(in   ) :: self
    double precision                     , intent(in   ) :: multiplier

    allocate(scaled,source=self)
    scaled%position=+self%position   &
         &          *     multiplier
    return
  end function Coordinates_Cartesian_Scalar_Multiply

  function Coordinates_Cartesian_Scalar_Divide(self,divisor) result(scaled)
    !!{RST
    Divide a Cartesian ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinate         ), allocatable   :: scaled
    class           (coordinateCartesian), intent(in   ) :: self
    double precision                     , intent(in   ) :: divisor

    allocate(scaled,source=self)
    scaled%position=+self%position &
         &          /     divisor
    return
  end function Coordinates_Cartesian_Scalar_Divide
  
  ! The scale() method scales a coordinate in-place, populating an allocatable, polymorphic result. Unlike
  ! the overloaded operator(*)/operator(/)---which return function results and so require an
  ! already-allocated (typically concrete-typed) assignment target---scale() can target an unallocated
  ! `class(coordinate), allocatable` argument directly (assigning an operator result to such a target would
  ! instead invoke the defined assignment on the unallocated, abstract-typed LHS and fail). It is therefore
  ! the scaling idiom used by the polymorphic scaler decorators. (Historically this was documented as a
  ! gfortran PR 37336 finalization workaround; that leak is fixed in gfortran >= 13, but scale() is retained
  ! because it remains the only way to scale into an unallocated polymorphic target.)
  subroutine Coordinates_Cartesian_Scale(self,scalar,selfScaled)
    !!{RST
    Scale a Cartesian ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinateCartesian), intent(in   )              :: self
    double precision                     , intent(in   )              :: scalar
    class           (coordinate         ), intent(inout), allocatable :: selfScaled

    allocate(selfScaled,source=self)
    selfScaled%position=+selfScaled%position &
         &              *           scalar
    return
  end subroutine Coordinates_Cartesian_Scale

  double precision function Coordinates_Cartesian_R_Spherical_Squared(self)
    !!{RST
    Return the squared spherical radius, :math:`r^2` of a Cartesian ``coordinate`` object.
    !!}
    implicit none
    class(coordinateCartesian), intent(in   ) :: self

    Coordinates_Cartesian_R_Spherical_Squared=sum(self%position**2)
    return
  end function Coordinates_Cartesian_R_Spherical_Squared

  ! Spherical coordinate object.
  subroutine Coordinates_Spherical_From_Cartesian(self,x)
    !!{RST
    Create a spherical ``coordinate`` object from a Cartesian vector.
    !!}
    implicit none
    class           (coordinateSpherical)              , intent(  out) :: self
    double precision                     , dimension(3), intent(in   ) :: x
    double precision                                                   :: phi , r, theta

    r     =sqrt (sum(x**2))
    if (r > 0.0d0) then
       theta =acos (x(3)/r   )
       phi   =atan2(x(2),x(1))
    else
       ! theta and phi angles are undefined at zero radius.
       theta=0.0d0
       phi  =0.0d0
    end if
    self%position=[r,theta,phi]
    return
  end subroutine Coordinates_Spherical_From_Cartesian

  function Coordinates_Spherical_To_Cartesian(self)
    !!{RST
    Return a Cartesian vector from a spherical ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateSpherical), intent(in   ) :: self
    double precision                     , dimension(3)  :: Coordinates_Spherical_To_Cartesian
    double precision                                     :: phi                               , r, &
         &                                                  theta

    r    =self%position(1)
    theta=self%position(2)
    phi  =self%position(3)
    Coordinates_Spherical_To_Cartesian= &
         & r                            &
         & *[                           &
         &   sin(theta)*cos(phi),       &
         &   sin(theta)*sin(phi),       &
         &   cos(theta)                 &
         &  ]
    return
  end function Coordinates_Spherical_To_Cartesian

  double precision function Coordinates_Spherical_R(self)
    !!{RST
    Return the :math:`r`-component of a Spherical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_R=self%position(1)
    return
  end function Coordinates_Spherical_R

  double precision function Coordinates_Spherical_Theta(self)
    !!{RST
    Return the :math:`\theta`-component of a Spherical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_Theta=self%position(2)
    return
  end function Coordinates_Spherical_Theta

  double precision function Coordinates_Spherical_Phi(self)
    !!{RST
    Return the :math:`\phi`-component of a Spherical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_Phi=self%position(3)
    return
  end function Coordinates_Spherical_Phi

  subroutine Coordinates_Spherical_Set_R(self,r)
    !!{RST
    Return the :math:`r`-component of a Spherical ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: r

    self%position(1)=r
    return
  end subroutine Coordinates_Spherical_Set_R

  subroutine Coordinates_Spherical_Set_Theta(self,theta)
    !!{RST
    Return the :math:`\theta`-component of a Spherical ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: theta

    self%position(2)=theta
    return
  end subroutine Coordinates_Spherical_Set_Theta

  subroutine Coordinates_Spherical_Set_Phi(self,phi)
    !!{RST
    Return the :math:`\phi`-component of a Spherical ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: phi

    self%position(3)=phi
    return
  end subroutine Coordinates_Spherical_Set_Phi

  double precision function Coordinates_Spherical_R_Spherical(self)
    !!{RST
    Return the spherical radius, :math:`r` of a spherical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_R_Spherical=self%position(1)
    return
  end function Coordinates_Spherical_R_Spherical

  double precision function Coordinates_Spherical_R_Spherical_Squared(self)
    !!{RST
    Return the squared spherical radius, :math:`r^2` of a spherical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_R_Spherical_Squared=self%position(1)**2
    return
  end function Coordinates_Spherical_R_Spherical_Squared

  function Coordinates_Spherical_Scalar_Multiply(self,multiplier) result(scaled)
    !!{RST
    Multiply a spherical ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinate         ), allocatable   :: scaled
    class           (coordinateSpherical), intent(in   ) :: self
    double precision                     , intent(in   ) :: multiplier

    allocate(scaled,source=self)
    scaled%position   =+self  %position
    scaled%position(1)=+scaled%position(1) &
         &             *multiplier
    return
  end function Coordinates_Spherical_Scalar_Multiply

  function Coordinates_Spherical_Scalar_Divide(self,divisor) result(scaled)
    !!{RST
    Divide a spherical ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinate         ), allocatable   :: scaled
    class           (coordinateSpherical), intent(in   ) :: self
    double precision                     , intent(in   ) :: divisor

    allocate(scaled,source=self)
    scaled%position(1)=+scaled%position(1) &
         &             /divisor
    return
  end function Coordinates_Spherical_Scalar_Divide

  ! See the note on Coordinates_Cartesian_Scale for why scale() is retained alongside the * and / operators.
  subroutine Coordinates_Spherical_Scale(self,scalar,selfScaled)
    !!{RST
    Scale a spherical ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinateSpherical), intent(in   )              :: self
    double precision                     , intent(in   )              :: scalar
    class           (coordinate         ), intent(inout), allocatable :: selfScaled

    allocate(selfScaled,source=self)
    selfScaled%position(1)=+selfScaled%position(1) &
         &                 *           scalar
    return
  end subroutine Coordinates_Spherical_Scale

  ! Cylindrical coordinate object.
  subroutine Coordinates_Cylindrical_From_Cartesian(self,x)
    !!{RST
    Create a cylindrical ``coordinate`` object from a Cartesian vector.
    !!}
    implicit none
    class           (coordinateCylindrical)              , intent(  out) :: self
    double precision                       , dimension(3), intent(in   ) :: x
    double precision                                                     :: phi , r, z

    r     =sqrt (x(1)**2+x(2)**2)
    phi   =atan2(x(2),x(1))
    z     =x(3)
    self%position=[r,phi,z]
    return
  end subroutine Coordinates_Cylindrical_From_Cartesian

  function Coordinates_Cylindrical_To_Cartesian(self)
    !!{RST
    Return a Cartesian vector from a cylindrical ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateCylindrical), intent(in   ) :: self
    double precision                       , dimension(3)  :: Coordinates_Cylindrical_To_Cartesian
    double precision                                       :: phi                                 , r, &
         &                                                    z

    r  =self%position(1)
    phi=self%position(2)
    z  =self%position(3)
    Coordinates_Cylindrical_To_Cartesian= &
         & [                              &
         &  r*cos(phi),                   &
         &  r*sin(phi),                   &
         &  z                             &
         & ]
    return
  end function Coordinates_Cylindrical_To_Cartesian

  double precision function Coordinates_Cylindrical_R(self)
    !!{RST
    Return the :math:`r`-component of a Cylindrical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self

    Coordinates_Cylindrical_R=self%position(1)
    return
  end function Coordinates_Cylindrical_R

  double precision function Coordinates_Cylindrical_Phi(self)
    !!{RST
    Return the :math:`\phi`-component of a Cylindrical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self

    Coordinates_Cylindrical_Phi=self%position(2)
    return
  end function Coordinates_Cylindrical_Phi

  double precision function Coordinates_Cylindrical_Z(self)
    !!{RST
    Return the :math:`z`-component of a Cylindrical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self

    Coordinates_Cylindrical_Z=self%position(3)
    return
  end function Coordinates_Cylindrical_Z

  subroutine Coordinates_Cylindrical_Set_R(self,r)
    !!{RST
    Return the :math:`r`-component of a Cylindrical ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: r

    self%position(1)=r
    return
  end subroutine Coordinates_Cylindrical_Set_R

  subroutine Coordinates_Cylindrical_Set_Phi(self,phi)
    !!{RST
    Return the :math:`\phi`-component of a Cylindrical ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: phi

    self%position(2)=phi
    return
  end subroutine Coordinates_Cylindrical_Set_Phi

  subroutine Coordinates_Cylindrical_Set_Z(self,z)
    !!{RST
    Return the :math:`z`-component of a Cylindrical ``coordinate`` object.
    !!}
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: z

    self%position(3)=z
    return
  end subroutine Coordinates_Cylindrical_Set_Z

  double precision function Coordinates_Cylindrical_R_Spherical_Squared(self)
    !!{RST
    Return the squared spherical radius, :math:`r^2` of a cylindrical ``coordinate`` object.
    !!}
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self

    Coordinates_Cylindrical_R_Spherical_Squared=self%position(1)**2+self%position(3)**2
    return
  end function Coordinates_Cylindrical_R_Spherical_Squared

  function Coordinates_Cylindrical_Scalar_Multiply(self,multiplier) result(scaled)
    !!{RST
    Multiply a cylindrical ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinate           ), allocatable   :: scaled
    class           (coordinateCylindrical), intent(in   ) :: self
    double precision                       , intent(in   ) :: multiplier

    allocate(scaled,source=self)
    scaled%position(1)=+scaled%position(1) &
         &             *multiplier
    scaled%position(3)=+scaled%position(3) &
         &             *multiplier
    return
  end function Coordinates_Cylindrical_Scalar_Multiply

  function Coordinates_Cylindrical_Scalar_Divide(self,divisor) result(scaled)
    !!{RST
    Divide a cylindrical ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinate           ), allocatable   :: scaled
    class           (coordinateCylindrical), intent(in   ) :: self
    double precision                       , intent(in   ) :: divisor

    allocate(scaled,source=self)
    scaled%position(1)=+scaled%position(1) &
         &             /divisor
    scaled%position(3)=+scaled%position(3) &
         &             /divisor
    return
  end function Coordinates_Cylindrical_Scalar_Divide

  ! See the note on Coordinates_Cartesian_Scale for why scale() is retained alongside the * and / operators.
  subroutine Coordinates_Cylindrical_Scale(self,scalar,selfScaled)
    !!{RST
    Scale a cylindrical ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinateCylindrical), intent(in   )              :: self
    double precision                       , intent(in   )              :: scalar
    class           (coordinate           ), intent(inout), allocatable :: selfScaled

    allocate(selfScaled,source=self)
    selfScaled%position(1)=+selfScaled%position(1) &
         &                 *           scalar
    selfScaled%position(3)=+selfScaled%position(3) &
         &                 *           scalar
    return
  end subroutine Coordinates_Cylindrical_Scale

  ! General functions.
  double precision function Coordinates_Radius_Cylindrical(self)
    implicit none
    class(coordinate           ), intent(in   ) :: self
    type (coordinateCylindrical)                :: coordinateCylindrical_

    select type (self)
    type is (coordinateCylindrical)
       ! Already in cylindrical coordinates - simply return the radial coordinate.
       Coordinates_Radius_Cylindrical=self%position(1)
    class default
       ! Convert to cylindrical coordinates and then return the radial coordinate.
       coordinateCylindrical_        =self
       Coordinates_Radius_Cylindrical=coordinateCylindrical_%position(1)
    end select
   return
  end function Coordinates_Radius_Cylindrical

  double precision function Coordinates_Radius_Spherical(self)
    implicit none
    class(coordinate), intent(in   ) :: self

    Coordinates_Radius_Spherical=sqrt(self%rSphericalSquared())
    return
  end function Coordinates_Radius_Spherical

  function Coordinates_Scalar_Multiply_Switched(multiplier,self) result(scaled)
    !!{RST
    Multiply a Cartesian ``coordinate`` object by a scalar.
    !!}
    implicit none
    class           (coordinate), allocatable   :: scaled
    class           (coordinate), intent(in   ) :: self
    double precision            , intent(in   ) :: multiplier

    scaled=self*multiplier
    return
  end function Coordinates_Scalar_Multiply_Switched

end module Coordinates

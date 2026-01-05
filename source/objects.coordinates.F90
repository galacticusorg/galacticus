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
Contains a module which implements the coordinates class.
!!}

module Coordinates
  !!{
  Implements the coordinates class.
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
     !!{
     The base coordinate object class.
     !!}
     double precision :: position(3)
   contains
     !![
     <methods>
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
     !!{
     A Cartesian coordinate object class.
     !!}
   contains
     !![
     <methods>
       <method description="Get the $x$-coordinate." method="x"   />
       <method description="Get the $y$-coordinate." method="y"   />
       <method description="Get the $z$-coordinate." method="z"   />
       <method description="set the $x$-coordinate." method="xSet"/>
       <method description="set the $y$-coordinate." method="ySet"/>
       <method description="set the $z$-coordinate." method="zSet"/>
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
     !!{
     A spherical coordinate object class.
     !!}
   contains
     !![
     <methods>
       <method description="Get the $r$-coordinate."      method="r"       />
       <method description="Get the $\theta$-coordinate." method="theta"   />
       <method description="Get the $\phi$-coordinate."   method="phi"     />
       <method description="set the $r$-coordinate."      method="rSet"    />
       <method description="set the $\theta$-coordinate." method="thetaSet"/>
       <method description="set the $\phi$-coordinate."   method="phiSet"  />
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
     !!{
     A cylindrical coordinate object class.
     !!}
   contains
     !![
     <methods>
       <method description="Get the $r$-coordinate."    method="r"     />
       <method description="Get the $\phi$-coordinate." method="phi"   />
       <method description="Get the $z$-coordinate."    method="z"     />
       <method description="set the $r$-coordinate."    method="rSet"  />
       <method description="set the $\phi$-coordinate." method="phiSet"/>
       <method description="set the $z$-coordinate."    method="zSet"  />
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

  subroutine Coordinates_Null_From(self,x)
    !!{
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
    !!{
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
    !!{
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
    !!{
    Assign a 3-component vector to a {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinate)              , intent(  out) :: coordinates
    double precision            , dimension(3), intent(in   ) :: x

    coordinates%position=x
    return
  end subroutine Coordinates_Assign_To

  subroutine Coordinates_Assign_From(x,coordinates)
    !!{
    Return a 3-component vector from a {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinate)              , intent(in   ) :: coordinates
    double precision            , dimension(3), intent(  out) :: x

    x=coordinates%position
    return
  end subroutine Coordinates_Assign_From

  ! Cartesian coordinate object.
  subroutine Coordinates_Cartesian_From_Cartesian(self,x)
    !!{
    Create a Cartesian {\normalfont \ttfamily coordinate} object from a Cartesian vector.
    !!}
    implicit none
    class           (coordinateCartesian)              , intent(  out) :: self
    double precision                     , dimension(3), intent(in   ) :: x

    self%position=x
    return
  end subroutine Coordinates_Cartesian_From_Cartesian

  function Coordinates_Cartesian_To_Cartesian(self)
    !!{
    Return a Cartesian vector from a Cartesian {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateCartesian), intent(in   ) :: self
    double precision                     , dimension(3)  :: Coordinates_Cartesian_To_Cartesian

    Coordinates_Cartesian_To_Cartesian=self%position
    return
  end function Coordinates_Cartesian_To_Cartesian

  double precision function Coordinates_Cartesian_X(self)
    !!{
    Return the $x$-component of a Cartesian {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateCartesian), intent(in   ) :: self

    Coordinates_Cartesian_X=self%position(1)
    return
  end function Coordinates_Cartesian_X

  double precision function Coordinates_Cartesian_Y(self)
    !!{
    Return the $y$-component of a Cartesian {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateCartesian), intent(in   ) :: self

    Coordinates_Cartesian_Y=self%position(2)
    return
  end function Coordinates_Cartesian_Y

  double precision function Coordinates_Cartesian_Z(self)
    !!{
    Return the $z$-component of a Cartesian {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateCartesian), intent(in   ) :: self

    Coordinates_Cartesian_Z=self%position(3)
    return
  end function Coordinates_Cartesian_Z

  subroutine Coordinates_Cartesian_Set_X(self,x)
    !!{
    Return the $x$-component of a Cartesian {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: x

    self%position(1)=x
    return
  end subroutine Coordinates_Cartesian_Set_X

  subroutine Coordinates_Cartesian_Set_Y(self,y)
    !!{
    Return the $y$-component of a Cartesian {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: y

    self%position(2)=y
    return
  end subroutine Coordinates_Cartesian_Set_Y

  subroutine Coordinates_Cartesian_Set_Z(self,z)
    !!{
    Return the $z$-component of a Cartesian {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: z

    self%position(3)=z
    return
  end subroutine Coordinates_Cartesian_Set_Z

  function Coordinates_Cartesian_Scalar_Multiply(self,multiplier) result(scaled)
    !!{
    Multiply a Cartesian {\normalfont \ttfamily coordinate} object by a scalar.
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
    !!{
    Divide a Cartesian {\normalfont \ttfamily coordinate} object by a scalar.
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
  
  !![
  <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336">
    <description>
      This function is needed to allow scaling of coordinate objects. It should not be needed as we overload the * and / operators
      for coordinate objects. But, until finalization is completed, function results are not finalized causing the overloaded *
      and / operators to leak memory. This is a workaround to avoid that.
    </description>
  </workaround>
  !!]
  subroutine Coordinates_Cartesian_Scale(self,scalar,selfScaled)
    !!{
    Scale a Cartesian {\normalfont \ttfamily coordinate} object by a scalar.
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
    !!{
    Return the squared spherical radius, $r^2$ of a Cartesian {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateCartesian), intent(in   ) :: self

    Coordinates_Cartesian_R_Spherical_Squared=sum(self%position**2)
    return
  end function Coordinates_Cartesian_R_Spherical_Squared

  ! Spherical coordinate object.
  subroutine Coordinates_Spherical_From_Cartesian(self,x)
    !!{
    Create a spherical {\normalfont \ttfamily coordinate} object from a Cartesian vector.
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
    !!{
    Return a Cartesian vector from a spherical {\normalfont \ttfamily coordinate} object.
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
    !!{
    Return the $r$-component of a Spherical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_R=self%position(1)
    return
  end function Coordinates_Spherical_R

  double precision function Coordinates_Spherical_Theta(self)
    !!{
    Return the $\theta$-component of a Spherical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_Theta=self%position(2)
    return
  end function Coordinates_Spherical_Theta

  double precision function Coordinates_Spherical_Phi(self)
    !!{
    Return the $\phi$-component of a Spherical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_Phi=self%position(3)
    return
  end function Coordinates_Spherical_Phi

  subroutine Coordinates_Spherical_Set_R(self,r)
    !!{
    Return the $r$-component of a Spherical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: r

    self%position(1)=r
    return
  end subroutine Coordinates_Spherical_Set_R

  subroutine Coordinates_Spherical_Set_Theta(self,theta)
    !!{
    Return the $\theta$-component of a Spherical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: theta

    self%position(2)=theta
    return
  end subroutine Coordinates_Spherical_Set_Theta

  subroutine Coordinates_Spherical_Set_Phi(self,phi)
    !!{
    Return the $\phi$-component of a Spherical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: phi

    self%position(3)=phi
    return
  end subroutine Coordinates_Spherical_Set_Phi

  double precision function Coordinates_Spherical_R_Spherical(self)
    !!{
    Return the spherical radius, $r$ of a spherical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_R_Spherical=self%position(1)
    return
  end function Coordinates_Spherical_R_Spherical

  double precision function Coordinates_Spherical_R_Spherical_Squared(self)
    !!{
    Return the squared spherical radius, $r^2$ of a spherical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateSpherical), intent(in   ) :: self

    Coordinates_Spherical_R_Spherical_Squared=self%position(1)**2
    return
  end function Coordinates_Spherical_R_Spherical_Squared

  function Coordinates_Spherical_Scalar_Multiply(self,multiplier) result(scaled)
    !!{
    Multiply a spherical {\normalfont \ttfamily coordinate} object by a scalar.
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
    !!{
    Divide a spherical {\normalfont \ttfamily coordinate} object by a scalar.
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

  !![
  <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336">
    <description>
      This function is needed to allow scaling of coordinate objects. It should not be needed as we overload the * and / operators
      for coordinate objects. But, until finalization is completed, function results are not finalized causing the overloaded *
      and / operators to leak memory. This is a workaround to avoid that.
    </description>
  </workaround>
  !!]
  subroutine Coordinates_Spherical_Scale(self,scalar,selfScaled)
    !!{
    Scale a spherical {\normalfont \ttfamily coordinate} object by a scalar.
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
    !!{
    Create a cylindrical {\normalfont \ttfamily coordinate} object from a Cartesian vector.
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
    !!{
    Return a Cartesian vector from a cylindrical {\normalfont \ttfamily coordinate} object.
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
    !!{
    Return the $r$-component of a Cylindrical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self

    Coordinates_Cylindrical_R=self%position(1)
    return
  end function Coordinates_Cylindrical_R

  double precision function Coordinates_Cylindrical_Phi(self)
    !!{
    Return the $\phi$-component of a Cylindrical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self

    Coordinates_Cylindrical_Phi=self%position(2)
    return
  end function Coordinates_Cylindrical_Phi

  double precision function Coordinates_Cylindrical_Z(self)
    !!{
    Return the $z$-component of a Cylindrical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self

    Coordinates_Cylindrical_Z=self%position(3)
    return
  end function Coordinates_Cylindrical_Z

  subroutine Coordinates_Cylindrical_Set_R(self,r)
    !!{
    Return the $r$-component of a Cylindrical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: r

    self%position(1)=r
    return
  end subroutine Coordinates_Cylindrical_Set_R

  subroutine Coordinates_Cylindrical_Set_Phi(self,phi)
    !!{
    Return the $\phi$-component of a Cylindrical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: phi

    self%position(2)=phi
    return
  end subroutine Coordinates_Cylindrical_Set_Phi

  subroutine Coordinates_Cylindrical_Set_Z(self,z)
    !!{
    Return the $z$-component of a Cylindrical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: z

    self%position(3)=z
    return
  end subroutine Coordinates_Cylindrical_Set_Z

  double precision function Coordinates_Cylindrical_R_Spherical_Squared(self)
    !!{
    Return the squared spherical radius, $r^2$ of a cylindrical {\normalfont \ttfamily coordinate} object.
    !!}
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self

    Coordinates_Cylindrical_R_Spherical_Squared=self%position(1)**2+self%position(3)**2
    return
  end function Coordinates_Cylindrical_R_Spherical_Squared

  function Coordinates_Cylindrical_Scalar_Multiply(self,multiplier) result(scaled)
    !!{
    Multiply a cylindrical {\normalfont \ttfamily coordinate} object by a scalar.
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
    !!{
    Divide a cylindrical {\normalfont \ttfamily coordinate} object by a scalar.
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

  !![
  <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336">
    <description>
      This function is needed to allow scaling of coordinate objects. It should not be needed as we overload the * and / operators
      for coordinate objects. But, until finalization is completed, function results are not finalized causing the overloaded *
      and / operators to leak memory. This is a workaround to avoid that.
    </description>
  </workaround>
  !!]
  subroutine Coordinates_Cylindrical_Scale(self,scalar,selfScaled)
    !!{
    Scale a cylindrical {\normalfont \ttfamily coordinate} object by a scalar.
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
    !!{
    Multiply a Cartesian {\normalfont \ttfamily coordinate} object by a scalar.
    !!}
    implicit none
    class           (coordinate), allocatable   :: scaled
    class           (coordinate), intent(in   ) :: self
    double precision            , intent(in   ) :: multiplier

    scaled=self*multiplier
    return
  end function Coordinates_Scalar_Multiply_Switched

end module Coordinates

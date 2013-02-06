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

!% Contains a module which implements the coordinates class.

module Coordinates
  !% Implements the coordinates class.
  implicit none
  private
  public :: assignment(=)
  
  ! Define assignment interfaces.
  interface assignment(=)
     module procedure Coordinates_Assign_To
     module procedure Coordinates_Assign_From
     module procedure Coordinates_Assign
  end interface assignment(=)

  type, public :: coordinate
     !% The base coordinate object class.
     double precision :: position(3)
   contains
     !@ <objectMethods>
     !@   <object>coordinate</object>
     !@   <objectMethod>
     !@     <method>toCartesian</method>
     !@     <description>Return the coordinates in a Cartesian system as a 3-element array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>fromCartesian</method>
     !@     <description>Set the coordinates from a Cartesian system specified as a 3-element array.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: toCartesian   => Coordinates_Null_To
     procedure :: fromCartesian => Coordinates_Null_From
  end type coordinate

  type, public, extends(coordinate) :: coordinateCartesian
     !% A Cartesian coordinate object class.
   contains
     !@ <objectMethods>
     !@   <object>coordinateCartesian</object>
     !@   <objectMethod>
     !@     <method>toCartesian</method>
     !@     <description>Return the coordinates in a Cartesian system as a 3-element array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>fromCartesian</method>
     !@     <description>Set the coordinates from a Cartesian system specified as a 3-element array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>x</method>
     !@     <description>Get the $x$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>y</method>
     !@     <description>Get the $y$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>z</method>
     !@     <description>Get the $z$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>xSet</method>
     !@     <description>set the $x$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>ySet</method>
     !@     <description>set the $y$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>zSet</method>
     !@     <description>set the $z$-coordinate.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: toCartesian   => Coordinates_Cartesian_To_Cartesian
     procedure :: fromCartesian => Coordinates_Cartesian_From_Cartesian
     procedure :: x             => Coordinates_Cartesian_X
     procedure :: y             => Coordinates_Cartesian_Y
     procedure :: z             => Coordinates_Cartesian_Z
     procedure :: xSet          => Coordinates_Cartesian_Set_X
     procedure :: ySet          => Coordinates_Cartesian_Set_Y
     procedure :: zSet          => Coordinates_Cartesian_Set_Z
  end type coordinateCartesian

  type, public, extends(coordinate) :: coordinateSpherical
     !% A spherical coordinate object class.
   contains
     !@ <objectMethods>
     !@   <object>coordinateSpherical</object>
     !@   <objectMethod>
     !@     <method>toCartesian</method>
     !@     <description>Return the coordinates in a Cartesian system as a 3-element array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>fromCartesian</method>
     !@     <description>Set the coordinates from a Cartesian system specified as a 3-element array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>r</method>
     !@     <description>Get the $r$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>theta</method>
     !@     <description>Get the $\theta$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>phi</method>
     !@     <description>Get the $\phi$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>rSet</method>
     !@     <description>set the $r$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>thetaSet</method>
     !@     <description>set the $\theta$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>phiSet</method>
     !@     <description>set the $\phi$-coordinate.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: toCartesian   => Coordinates_Spherical_To_Cartesian
     procedure :: fromCartesian => Coordinates_Spherical_From_Cartesian
     procedure :: r             => Coordinates_Spherical_R
     procedure :: theta         => Coordinates_Spherical_Theta
     procedure :: phi           => Coordinates_Spherical_Phi
     procedure :: rSet          => Coordinates_Spherical_Set_R
     procedure :: thetaSet      => Coordinates_Spherical_Set_Theta
     procedure :: phiSet        => Coordinates_Spherical_Set_Phi
  end type coordinateSpherical

  type, public, extends(coordinate) :: coordinateCylindrical
     !% A cylindrical coordinate object class.
   contains
     !@ <objectMethods>
     !@   <object>coordinateCylindrical</object>
     !@   <objectMethod>
     !@     <method>toCartesian</method>
     !@     <description>Return the coordinates in a Cartesian system as a 3-element array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>fromCartesian</method>
     !@     <description>Set the coordinates from a Cartesian system specified as a 3-element array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>r</method>
     !@     <description>Get the $r$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>phi</method>
     !@     <description>Get the $\phi$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>z</method>
     !@     <description>Get the $z$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>rSet</method>
     !@     <description>set the $r$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>phiSet</method>
     !@     <description>set the $\phi$-coordinate.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>zSet</method>
     !@     <description>set the $z$-coordinate.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: toCartesian   => Coordinates_Cylindrical_To_Cartesian
     procedure :: fromCartesian => Coordinates_Cylindrical_From_Cartesian
     procedure :: r             => Coordinates_Cylindrical_R
     procedure :: phi           => Coordinates_Cylindrical_Phi
     procedure :: z             => Coordinates_Cylindrical_Z
     procedure :: rSet          => Coordinates_Cylindrical_Set_R
     procedure :: phiSet        => Coordinates_Cylindrical_Set_Phi
     procedure :: zSet          => Coordinates_Cylindrical_Set_Z
  end type coordinateCylindrical

contains

  subroutine Coordinates_Null_From(self,x)
    !% Set generic coordinate object from Cartesian point. Simply quits with an error.
    use Galacticus_Error
    implicit none
    class           (coordinate), intent(  out)               :: self
    double precision            , intent(in   ), dimension(3) :: x

    call Galacticus_Error_Report('Coordinates_Null_To','no transformation from cartesian coordinates defined')
    return
  end subroutine Coordinates_Null_From

  function Coordinates_Null_To(self)
    !% Convert generic coordinate object to Cartesian point. Simply quits with an error.
    use Galacticus_Error
    implicit none
    class           (coordinate), intent(in   ) :: self
    double precision            , dimension(3)  :: Coordinates_Null_To

    call Galacticus_Error_Report('Coordinates_Null_To','no transformation to cartesian coordinates defined')
    return
  end function Coordinates_Null_To

  subroutine Coordinates_Assign(coordinatesTo,coordinatesFrom)
    !% Assign one coordinate object to another, automatically handling the conversion between coordinate systems.
    implicit none
    class           (coordinate), intent(  out) :: coordinatesTo
    class           (coordinate), intent(in   ) :: coordinatesFrom
    double precision            , dimension(3)  :: x

    ! Assign by transforming through cartesian coordinates.
    x=coordinatesFrom%toCartesian()
    call coordinatesTo%fromCartesian(x)
    return
  end subroutine Coordinates_Assign

  subroutine Coordinates_Assign_To(coordinates,x)
    !% Assign a 3-component vector to a {\tt coordinate} object.
    implicit none
    class           (coordinate), intent(  out)               :: coordinates
    double precision            , intent(in   ), dimension(3) :: x

    coordinates%position=x
    return
  end subroutine Coordinates_Assign_To

  subroutine Coordinates_Assign_From(x,coordinates)
    !% Return a 3-component vector from a {\tt coordinate} object.
    implicit none
    class           (coordinate), intent(in   )               :: coordinates
    double precision            , intent(  out), dimension(3) :: x

    x=coordinates%position
    return
  end subroutine Coordinates_Assign_From

  ! Cartesian coordinate object.
  subroutine Coordinates_Cartesian_From_Cartesian(self,x)
    !% Create a Cartesian {\tt coordinate} object from a Cartesian vector.
    implicit none
    class           (coordinateCartesian), intent(  out)               :: self
    double precision                     , intent(in   ), dimension(3) :: x
    
    self%position=x
    return
  end subroutine Coordinates_Cartesian_From_Cartesian

  function Coordinates_Cartesian_To_Cartesian(self)
    !% Return a Cartesian vector from a Cartesian {\tt coordinate} object.
    implicit none
    class           (coordinateCartesian), intent(in   ) :: self
    double precision                     , dimension(3)  :: Coordinates_Cartesian_To_Cartesian

    Coordinates_Cartesian_To_Cartesian=self%position
    return
  end function Coordinates_Cartesian_To_Cartesian

  double precision function Coordinates_Cartesian_X(self)
    !% Return the $x$-component of a Cartesian {\tt coordinate} object.
    implicit none
    class(coordinateCartesian), intent(in   ) :: self
    
    Coordinates_Cartesian_X=self%position(1)
    return
  end function Coordinates_Cartesian_X

  double precision function Coordinates_Cartesian_Y(self)
    !% Return the $y$-component of a Cartesian {\tt coordinate} object.
    implicit none
    class(coordinateCartesian), intent(in   ) :: self
    
    Coordinates_Cartesian_Y=self%position(2)
    return
  end function Coordinates_Cartesian_Y

  double precision function Coordinates_Cartesian_Z(self)
    !% Return the $z$-component of a Cartesian {\tt coordinate} object.
    implicit none
    class(coordinateCartesian), intent(in   ) :: self
    
    Coordinates_Cartesian_Z=self%position(3)
    return
  end function Coordinates_Cartesian_Z

  subroutine Coordinates_Cartesian_Set_X(self,x)
    !% Return the $x$-component of a Cartesian {\tt coordinate} object.
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: x
 
    self%position(1)=x
    return
  end subroutine Coordinates_Cartesian_Set_X

  subroutine Coordinates_Cartesian_Set_Y(self,y)
    !% Return the $y$-component of a Cartesian {\tt coordinate} object.
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: y
 
    self%position(2)=y
    return
  end subroutine Coordinates_Cartesian_Set_Y

  subroutine Coordinates_Cartesian_Set_Z(self,z)
    !% Return the $z$-component of a Cartesian {\tt coordinate} object.
    implicit none
    class           (coordinateCartesian), intent(inout) :: self
    double precision                     , intent(in   ) :: z
 
    self%position(3)=z
    return
  end subroutine Coordinates_Cartesian_Set_Z

  ! Spherical coordinate object.
  subroutine Coordinates_Spherical_From_Cartesian(self,x)
    !% Create a spherical {\tt coordinate} object from a Cartesian vector.
    implicit none
    class           (coordinateSpherical), intent(  out)               :: self
    double precision                     , intent(in   ), dimension(3) :: x
    double precision                                                   :: r,theta,phi
    
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
    !% Return a Cartesian vector from a spherical {\tt coordinate} object.
    implicit none
    class           (coordinateSpherical), intent(in   ) :: self
    double precision                     , dimension(3)  :: Coordinates_Spherical_To_Cartesian
    double precision                                     :: r,theta,phi

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
    !% Return the $r$-component of a Spherical {\tt coordinate} object.
    implicit none
    class(coordinateSpherical), intent(in   ) :: self
    
    Coordinates_Spherical_R=self%position(1)
    return
  end function Coordinates_Spherical_R

  double precision function Coordinates_Spherical_Theta(self)
    !% Return the $\theta$-component of a Spherical {\tt coordinate} object.
    implicit none
    class(coordinateSpherical), intent(in   ) :: self
    
    Coordinates_Spherical_Theta=self%position(2)
    return
  end function Coordinates_Spherical_Theta

  double precision function Coordinates_Spherical_Phi(self)
    !% Return the $\phi$-component of a Spherical {\tt coordinate} object.
    implicit none
    class(coordinateSpherical), intent(in   ) :: self
    
    Coordinates_Spherical_Phi=self%position(3)
    return
  end function Coordinates_Spherical_Phi

  subroutine Coordinates_Spherical_Set_R(self,r)
    !% Return the $r$-component of a Spherical {\tt coordinate} object.
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: r
 
    self%position(1)=r
    return
  end subroutine Coordinates_Spherical_Set_R

  subroutine Coordinates_Spherical_Set_Theta(self,theta)
    !% Return the $\theta$-component of a Spherical {\tt coordinate} object.
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: theta
 
    self%position(2)=theta
    return
  end subroutine Coordinates_Spherical_Set_Theta
  
  subroutine Coordinates_Spherical_Set_Phi(self,phi)
    !% Return the $\phi$-component of a Spherical {\tt coordinate} object.
    implicit none
    class           (coordinateSpherical), intent(inout) :: self
    double precision                     , intent(in   ) :: phi
 
    self%position(3)=phi
    return
  end subroutine Coordinates_Spherical_Set_Phi

  ! Cylindrical coordinate object.
  subroutine Coordinates_Cylindrical_From_Cartesian(self,x)
    !% Create a cylindrical {\tt coordinate} object from a Cartesian vector.
    implicit none
    class           (coordinateCylindrical), intent(  out)               :: self
    double precision                       , intent(in   ), dimension(3) :: x
    double precision                                                     :: r,phi,z
    
    r     =sqrt (x(1)**2+x(2)**2)
    phi   =atan2(x(2),x(1))
    z     =x(3)
    self%position=[r,phi,z]
    return
  end subroutine Coordinates_Cylindrical_From_Cartesian

  function Coordinates_Cylindrical_To_Cartesian(self)
    !% Return a Cartesian vector from a cylindrical {\tt coordinate} object.
    implicit none
    class           (coordinateCylindrical), intent(in   ) :: self
    double precision                       , dimension(3)  :: Coordinates_Cylindrical_To_Cartesian
    double precision                                       :: r,phi,z

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
    !% Return the $r$-component of a Cylindrical {\tt coordinate} object.
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self
    
    Coordinates_Cylindrical_R=self%position(1)
    return
  end function Coordinates_Cylindrical_R

  double precision function Coordinates_Cylindrical_Phi(self)
    !% Return the $\phi$-component of a Cylindrical {\tt coordinate} object.
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self
    
    Coordinates_Cylindrical_Phi=self%position(2)
    return
  end function Coordinates_Cylindrical_Phi

  double precision function Coordinates_Cylindrical_Z(self)
    !% Return the $z$-component of a Cylindrical {\tt coordinate} object.
    implicit none
    class(coordinateCylindrical), intent(in   ) :: self
    
    Coordinates_Cylindrical_Z=self%position(3)
    return
  end function Coordinates_Cylindrical_Z

  subroutine Coordinates_Cylindrical_Set_R(self,r)
    !% Return the $r$-component of a Cylindrical {\tt coordinate} object.
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: r
 
    self%position(1)=r
    return
  end subroutine Coordinates_Cylindrical_Set_R

  subroutine Coordinates_Cylindrical_Set_Phi(self,phi)
    !% Return the $\phi$-component of a Cylindrical {\tt coordinate} object.
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: phi
 
    self%position(2)=phi
    return
  end subroutine Coordinates_Cylindrical_Set_Phi
  
  subroutine Coordinates_Cylindrical_Set_Z(self,z)
    !% Return the $z$-component of a Cylindrical {\tt coordinate} object.
    implicit none
    class           (coordinateCylindrical), intent(inout) :: self
    double precision                       , intent(in   ) :: z
 
    self%position(3)=z
    return
  end subroutine Coordinates_Cylindrical_Set_Z

end module Coordinates

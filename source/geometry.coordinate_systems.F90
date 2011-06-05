!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations related to coordinate systems and transformations.

module Coordinate_Systems
  !% Implements calculations related to coordinate systems and transformations.
  private
  public :: Coordinates_Cylindrical_To_Spherical, Coordinates_Cartesian_To_Spherical, Coordinates_Spherical_To_Cylindrical, Coordinates_Cartesian_To_Cylindrical

contains

  function Coordinates_Cartesian_To_Spherical(cartesianPosition)
    !% Convert $(x,y,z)$ in Cartesian coordinates into $(r,\theta,\phi)$ in spherical coordinates, with $\theta=0$ corresponding
    !% to the $z$-axis and $\phi=0$ corresponding to the $x$-axis.
    implicit none
    double precision, dimension(3)             :: Coordinates_Cartesian_To_Spherical
    double precision, dimension(3), intent(in) :: cartesianPosition

    ! Spherical radius.
    Coordinates_Cartesian_To_Spherical(1)=dsqrt(cartesianPosition(1)**2+cartesianPosition(2)**2+cartesianPosition(3)**2)
    ! Spherical theta.
    Coordinates_Cartesian_To_Spherical(2)=dacos(cartesianPosition(3)/dsqrt(cartesianPosition(1)**2+cartesianPosition(2)**2))
    ! Spherical phi.
    Coordinates_Cartesian_To_Spherical(3)=datan2(cartesianPosition(2),cartesianPosition(1))
    return
  end function Coordinates_Cartesian_To_Spherical

  function Coordinates_Cartesian_To_Cylindrical(cartesianPosition)
    !% Convert $(x,y,z)$ in Cartesian coordinates into $(r,\phi,z)$ in cylindrical coordinates, with $\phi=0$ corresponding to the $x$-axis.
    implicit none
    double precision, dimension(3)             :: Coordinates_Cartesian_To_Cylindrical
    double precision, dimension(3), intent(in) :: cartesianPosition

    ! Cylindrical radius.
    Coordinates_Cartesian_To_Cylindrical(1)=dsqrt(cartesianPosition(1)**2+cartesianPosition(2)**2)
    ! Spherical phi.
    Coordinates_Cartesian_To_Cylindrical(2)=datan2(cartesianPosition(2),cartesianPosition(1))
    ! Spherical z.
    Coordinates_Cartesian_To_Cylindrical(3)=cartesianPosition(3)
    return
  end function Coordinates_Cartesian_To_Cylindrical

  function Coordinates_Cylindrical_To_Spherical(cylindricalPosition)
    !% Convert $(R,\phi,z)$ in cylindrical coordinates into $(r,\theta,\phi)$ in spherical coordinates, with $\phi=0$
    !% corresponding to the $x$-axis.
    implicit none
    double precision, dimension(3)             :: Coordinates_Cylindrical_To_Spherical
    double precision, dimension(3), intent(in) :: cylindricalPosition

    ! Spherical radius.
    Coordinates_Cylindrical_To_Spherical(1)=dsqrt(cylindricalPosition(1)**2+cylindricalPosition(3)**2)
    ! Spherical theta.
    Coordinates_Cylindrical_To_Spherical(2)=dacos(cylindricalPosition(3)/cylindricalPosition(1))
    ! Spherical phi.
    Coordinates_Cylindrical_To_Spherical(3)=cylindricalPosition(2)
    return
  end function Coordinates_Cylindrical_To_Spherical

  function Coordinates_Spherical_To_Cylindrical(sphericalPosition)
    !% Convert $(r,\theta,\phi)$ in spherical coordinates into $(R,\phi,z)$ in cylindrical coordinates, with $\phi=0$
    !% corresponding to the $x$-axis.
    implicit none
    double precision, dimension(3)             :: Coordinates_Spherical_To_Cylindrical
    double precision, dimension(3), intent(in) :: sphericalPosition

    ! Cylindrical radius.
    Coordinates_Spherical_To_Cylindrical(1)=sphericalPosition(1)*dsin(sphericalPosition(2))
    ! Cylinderical phi.
    Coordinates_Spherical_To_Cylindrical(2)=sphericalPosition(3)
    ! Cylindrical z.
    Coordinates_Spherical_To_Cylindrical(3)=sphericalPosition(1)*dcos(sphericalPosition(2))
    return
  end function Coordinates_Spherical_To_Cylindrical

end module Coordinate_Systems

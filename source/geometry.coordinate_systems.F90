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
Contains a module which implements calculations related to coordinate systems and transformations.
!!}

module Coordinate_Systems
  !!{
  Implements calculations related to coordinate systems and transformations.
  !!}
  implicit none
  private
  public :: Coordinates_Cylindrical_To_Spherical, Coordinates_Cartesian_To_Spherical, Coordinates_Spherical_To_Cylindrical, Coordinates_Cartesian_To_Cylindrical

contains

  function Coordinates_Cartesian_To_Spherical(cartesianPosition)
    !!{
    Convert $(x,y,z)$ in Cartesian coordinates into $(r,\theta,\phi)$ in spherical coordinates, with $\theta=0$ corresponding
    to the $z$-axis and $\phi=0$ corresponding to the $x$-axis.
    !!}
    implicit none
    double precision, dimension(3)                :: Coordinates_Cartesian_To_Spherical
    double precision, dimension(3), intent(in   ) :: cartesianPosition

    ! Spherical radius.
    Coordinates_Cartesian_To_Spherical(1)=sqrt(cartesianPosition(1)**2+cartesianPosition(2)**2+cartesianPosition(3)**2)
    ! Check for zero radius.
    if (Coordinates_Cartesian_To_Spherical(1) == 0.0d0) then
       ! Other coordinates are arbitrary - set to zero.
       Coordinates_Cartesian_To_Spherical(2:3)=0.0d0
    else
       ! Spherical theta.
       Coordinates_Cartesian_To_Spherical(2)=acos(cartesianPosition(3)/Coordinates_Cartesian_To_Spherical(1))
       ! Spherical phi.
       Coordinates_Cartesian_To_Spherical(3)=atan2(cartesianPosition(2),cartesianPosition(1))
    end if
    return
  end function Coordinates_Cartesian_To_Spherical

  function Coordinates_Cartesian_To_Cylindrical(cartesianPosition)
    !!{
    Convert $(x,y,z)$ in Cartesian coordinates into $(r,\phi,z)$ in cylindrical coordinates, with $\phi=0$ corresponding to the $x$-axis.
    !!}
    implicit none
    double precision, dimension(3)                :: Coordinates_Cartesian_To_Cylindrical
    double precision, dimension(3), intent(in   ) :: cartesianPosition

    ! Cylindrical radius.
    Coordinates_Cartesian_To_Cylindrical(1)=sqrt(cartesianPosition(1)**2+cartesianPosition(2)**2)
    ! Spherical phi.
    Coordinates_Cartesian_To_Cylindrical(2)=atan2(cartesianPosition(2),cartesianPosition(1))
    ! Spherical z.
    Coordinates_Cartesian_To_Cylindrical(3)=cartesianPosition(3)
    return
  end function Coordinates_Cartesian_To_Cylindrical

  function Coordinates_Cylindrical_To_Spherical(cylindricalPosition)
    !!{
    Convert $(R,\phi,z)$ in cylindrical coordinates into $(r,\theta,\phi)$ in spherical coordinates, with $\phi=0$
    corresponding to the $x$-axis.
    !!}
    implicit none
    double precision, dimension(3)                :: Coordinates_Cylindrical_To_Spherical
    double precision, dimension(3), intent(in   ) :: cylindricalPosition

    ! Spherical radius.
    Coordinates_Cylindrical_To_Spherical(1)=sqrt(cylindricalPosition(1)**2+cylindricalPosition(3)**2)
    ! Check for zero radius.
    if (Coordinates_Cylindrical_To_Spherical(1) == 0.0d0) then
       ! Angular coordinate is undefined - set to zero.
       Coordinates_Cylindrical_To_Spherical(2)=0.0d0
    else
       ! Spherical theta.
       Coordinates_Cylindrical_To_Spherical(2)=acos(cylindricalPosition(3)/Coordinates_Cylindrical_To_Spherical(1))
    end if
    ! Spherical phi.
    Coordinates_Cylindrical_To_Spherical(3)=cylindricalPosition(2)
    return
  end function Coordinates_Cylindrical_To_Spherical

  function Coordinates_Spherical_To_Cylindrical(sphericalPosition)
    !!{
    Convert $(r,\theta,\phi)$ in spherical coordinates into $(R,\phi,z)$ in cylindrical coordinates, with $\phi=0$
    corresponding to the $x$-axis.
    !!}
    implicit none
    double precision, dimension(3)                :: Coordinates_Spherical_To_Cylindrical
    double precision, dimension(3), intent(in   ) :: sphericalPosition

    ! Cylindrical radius.
    Coordinates_Spherical_To_Cylindrical(1)=sphericalPosition(1)*sin(sphericalPosition(2))
    ! Cylindrical phi.
    Coordinates_Spherical_To_Cylindrical(2)=sphericalPosition(3)
    ! Cylindrical z.
    Coordinates_Spherical_To_Cylindrical(3)=sphericalPosition(1)*cos(sphericalPosition(2))
    return
  end function Coordinates_Spherical_To_Cylindrical

end module Coordinate_Systems

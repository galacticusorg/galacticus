!% Contains a module which implements calculations related to coordinate systems and transformations.

module Coordinate_Systems
  !% Implements calculations related to coordinate systems and transformations.
  private
  public :: Coordinates_Cylindrical_To_Spherical, Coordinates_Cartesian_To_Spherical, Coordinates_Spherical_To_Cylindrical

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

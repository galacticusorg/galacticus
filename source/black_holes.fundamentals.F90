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
Contains a module which implements fundamental properties of black holes.
!!}

module Black_Hole_Fundamentals
  !!{
  Implements fundamental properties of black holes.
  !!}
  implicit none
  private
  public :: Black_Hole_ISCO_Radius                   , Black_Hole_ISCO_Specific_Energy       , Black_Hole_Gravitational_Radius    , Black_Hole_Frame_Dragging_Frequency, &
       &    Black_Hole_Metric_A_Factor               , Black_Hole_Metric_D_Factor            , Black_Hole_Horizon_Radius          , Black_Hole_Static_Radius           , &
       &    Black_Hole_ISCO_Specific_Angular_Momentum, Black_Hole_Rotational_Energy_Spin_Down, Black_Hole_Eddington_Accretion_Rate

  ! Identifiers for unit system options.
  integer, parameter, public :: unitsGravitational=0
  integer, parameter, public :: unitsPhysical     =1

  ! Identifiers for orbit orientation options.
  integer, parameter, public :: orbitPrograde     =0
  integer, parameter, public :: orbitRetrograde   =1

  ! Generic interfaces for functions.
  interface Black_Hole_ISCO_Radius
     module procedure Black_Hole_ISCO_Radius_Node
     module procedure Black_Hole_ISCO_Radius_Spin
  end interface Black_Hole_ISCO_Radius
  interface Black_Hole_Rotational_Energy_Spin_Down
     module procedure Black_Hole_Rotational_Energy_Spin_Down_Node
     module procedure Black_Hole_Rotational_Energy_Spin_Down_Spin
  end interface Black_Hole_Rotational_Energy_Spin_Down
  interface Black_Hole_Metric_A_Factor
     module procedure Black_Hole_Metric_A_Factor_Node
     module procedure Black_Hole_Metric_A_Factor_Spin
  end interface Black_Hole_Metric_A_Factor
  interface Black_Hole_Metric_D_Factor
     module procedure Black_Hole_Metric_D_Factor_Node
     module procedure Black_Hole_Metric_D_Factor_Spin
  end interface Black_Hole_Metric_D_Factor
  interface Black_Hole_Frame_Dragging_Frequency
     module procedure Black_Hole_Frame_Dragging_Frequency_Node
     module procedure Black_Hole_Frame_Dragging_Frequency_Spin
  end interface Black_Hole_Frame_Dragging_Frequency
  interface Black_Hole_Horizon_Radius
     module procedure Black_Hole_Horizon_Radius_Node
     module procedure Black_Hole_Horizon_Radius_Spin
  end interface Black_Hole_Horizon_Radius
  interface Black_Hole_Static_Radius
     module procedure Black_Hole_Static_Radius_Node
     module procedure Black_Hole_Static_Radius_Spin
  end interface Black_Hole_Static_Radius
  interface Black_Hole_ISCO_Specific_Energy
     module procedure Black_Hole_ISCO_Specific_Energy_Node
     module procedure Black_Hole_ISCO_Specific_Energy_Spin
  end interface Black_Hole_ISCO_Specific_Energy

contains

  double precision function Black_Hole_ISCO_Radius_Spin(spinBlackHole,orbit)
    !!{
    Returns the radius (in gravitational units and for a prograde or retrograde orbit) of the innermost stable
    circular orbit for a black hole with spin {\normalfont \ttfamily spinBlackHole}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision      , intent(in   )           :: spinBlackHole
    integer               , intent(in   ), optional :: orbit
    integer         , save                          :: orbitPrevious
    double precision, save                          :: spinBlackHolePrevious=2.0d0, radiusISCOPrevious
    !$omp threadprivate(orbitPrevious,spinBlackHolePrevious,radiusISCOPrevious)
    integer                                         :: orbitActual
    double precision                                :: A1Factor                   , A2Factor

    ! Determine what orbit to use.
    if (present(orbit)) then
       orbitActual=orbit
    else
       orbitActual=orbitPrograde
    end if

    ! Check if we're being called with the same spin and orbit as the last time.
    if (orbitActual == orbitPrevious .and. spinBlackHole == spinBlackHolePrevious) then
       ! We are, so just return our previously computed radius.
       Black_Hole_ISCO_Radius_Spin=radiusISCOPrevious
    else
       ! We are not, so compute the radius of the innermost stable circular orbit in units of the gravitational radius.
       A1Factor=A1(spinBlackHole)
       A2Factor=A2(spinBlackHole)
       select case (orbitActual)
       case (orbitPrograde)
          Black_Hole_ISCO_Radius_Spin=3.0d0+A2Factor-sqrt((3.0d0-A1Factor)*(3.0d0+A1Factor+2.0d0*A2Factor))
       case (orbitRetrograde)
          Black_Hole_ISCO_Radius_Spin=3.0d0+A2Factor+sqrt((3.0d0-A1Factor)*(3.0d0+A1Factor+2.0d0*A2Factor))
       case default
          Black_Hole_ISCO_Radius_Spin=0.0d0
          call Error_Report('unrecognized orbit parameter'//{introspection:location})
       end select
       orbitPrevious        =orbitActual
       spinBlackHolePrevious=spinBlackHole
       radiusISCOPrevious   =Black_Hole_ISCO_Radius_Spin
    end if
    return
  end function Black_Hole_ISCO_Radius_Spin

  double precision function Black_Hole_Eddington_Accretion_Rate(blackHole)
    !!{
    Return the Eddington accretion rate (in $M_\odot$ Gyr$^{-1}$) for the black hole in {\normalfont \ttfamily blackHole}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole
    use :: Numerical_Constants_Astronomical, only : gigaYear
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : gravitationalConstant , speedLight, thomsonCrossSection
    implicit none
    class(nodeComponentBlackHole), intent(inout) :: blackHole

    Black_Hole_Eddington_Accretion_Rate=4.0d0*Pi*gravitationalConstant*blackHole%mass()*massHydrogenAtom&
         &*gigaYear/thomsonCrossSection/speedLight
    return
  end function Black_Hole_Eddington_Accretion_Rate

  double precision function Black_Hole_ISCO_Radius_Node(blackHole,units,orbit)
    !!{
    Returns the radius (in physical or gravitational units and for a prograde or retrograde orbit) of the innermost stable
    circular orbit for the black hole in {\normalfont \ttfamily blackHole}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (nodeComponentBlackHole), intent(inout)           :: blackHole
    integer                                 , intent(in   ), optional :: orbit        , units
    integer                                                           :: orbitActual  , unitsActual
    double precision                                                  :: spinBlackHole

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Determine what orbit to use.
    if (present(orbit)) then
       orbitActual=orbit
    else
       orbitActual=orbitPrograde
    end if

    ! Get the black hole spin.
    spinBlackHole=blackHole%spin()

    ! Compute the radius of the innermost stable circular orbit in units of the gravitational radius.
    Black_Hole_ISCO_Radius_Node=Black_Hole_ISCO_Radius_Spin(spinBlackHole,orbitActual)

    ! Convert to physical units if necessary.
    if (unitsActual == unitsPhysical) Black_Hole_ISCO_Radius_Node=Black_Hole_ISCO_Radius_Node*Black_Hole_Gravitational_Radius(blackHole)
    return
  end function Black_Hole_ISCO_Radius_Node

  double precision function Black_Hole_ISCO_Specific_Energy_Node(blackHole,units,orbit)
    !!{
    Returns the specific energy (in physical or gravitational units and for a prograde or retrograde orbit) of the innermost
    stable circular orbit for the given {\normalfont \ttfamily blackHole}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (nodeComponentBlackHole), intent(inout)           :: blackHole
    integer                                 , intent(in   ), optional :: orbit        , units
    integer                                                           :: unitsActual
    double precision                                                  :: spinBlackHole

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get the black hole spin.
    spinBlackHole=blackHole%spin()

    ! Get the dimensionless ISCO energy.
    Black_Hole_ISCO_Specific_Energy_Node=Black_Hole_ISCO_Specific_Energy_Spin(spinBlackHole,orbit)

    ! Convert to physical units if necessary.
    if (unitsActual == unitsPhysical) Black_Hole_ISCO_Specific_Energy_Node=Black_Hole_ISCO_Specific_Energy_Node&
         &*gravitationalConstant_internal*blackHole%mass()**2/Black_Hole_Gravitational_Radius(blackHole)
    return
  end function Black_Hole_ISCO_Specific_Energy_Node

  double precision function Black_Hole_ISCO_Specific_Energy_Spin(spinBlackHole,orbit)
    !!{
    Returns the specific energy (in physical or gravitational units and for a prograde or retrograde orbit) of the innermost
    stable circular orbit for a black hole of given {\normalfont \ttfamily spinBlackHole}.
    !!}
    implicit none
    double precision           , intent(inout)           :: spinBlackHole
    integer                    , intent(in   ), optional :: orbit
    ! Maximum spin above which we use series solution.
    double precision, parameter                          :: maximumSpin        =0.99999d0
    ! Coefficients in the expansion of ISCO energy with spin close to maximal spin.
    double precision, parameter                          :: coefficientFirst   =0.9164864242d0, coefficientZeroth=0.5773502693d0
    double precision                                     :: blackHoleIscoRadius

    ! Get black hole ISCO radius.
    blackHoleIscoRadius=Black_Hole_ISCO_Radius(spinBlackHole,orbit)

    ! Compute the specific energy in gravitational units.
    if (spinBlackHole >= maximumSpin) then
       Black_Hole_ISCO_Specific_Energy_Spin=coefficientZeroth+coefficientFirst*(1.0d0-spinBlackHole)**(1.0d0/3.0d0)
    else
       Black_Hole_ISCO_Specific_Energy_Spin=(blackHoleIscoRadius**2-2.0d0*blackHoleIscoRadius+spinBlackHole &
            &*sqrt(blackHoleIscoRadius))/blackHoleIscoRadius/sqrt(blackHoleIscoRadius**2-3.0d0*blackHoleIscoRadius+2.0d0 &
            &*spinBlackHole*sqrt(blackHoleIscoRadius))
    end if

    return
  end function Black_Hole_ISCO_Specific_Energy_Spin

  double precision function Black_Hole_ISCO_Specific_Angular_Momentum(blackHole,units,orbit)
    !!{
    Returns the specific angular momentum (in physical or gravitational units and for a prograde or retrograde orbit) of the
    innermost stable circular orbit for the black hole in {\normalfont \ttfamily blackHole}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (nodeComponentBlackHole)           , intent(inout)           :: blackHole
    integer                                            , intent(in   ), optional :: orbit                            , units
    ! Maximum spin above which we use series solution.
    double precision                        , parameter                          :: maximumSpin        =0.99999d0
    ! Coefficients in the expansion of ISCO energy with spin close to maximal spin.
    double precision                        , parameter                          :: coefficientFirst   =1.832972849d0, coefficientZeroth=1.154700538d0
    integer                                                                      :: orbitActual                      , unitsActual
    double precision                                                             :: blackHoleIscoRadius              , spinBlackHole

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Determine what orbit to use.
    if (present(orbit)) then
       orbitActual=orbit
    else
       orbitActual=orbitPrograde
    end if

    ! Get black hole spin and ISCO radius.
    spinBlackHole      =blackHole%spin()
    blackHoleIscoRadius=Black_Hole_ISCO_Radius(blackHole,units=unitsGravitational,orbit=orbitActual)

    ! Compute the specific angular momentum in gravitational units.
    if (spinBlackHole > maximumSpin) then
       Black_Hole_ISCO_Specific_Angular_Momentum=coefficientZeroth+coefficientFirst*(1.0d0-spinBlackHole)**(1.0d0/3.0d0)
    else
       Black_Hole_ISCO_Specific_Angular_Momentum=sqrt(blackHoleIscoRadius)*(blackHoleIscoRadius**2-2.0d0*spinBlackHole &
            &*sqrt(blackHoleIscoRadius)+spinBlackHole**2)/blackHoleIscoRadius/sqrt(blackHoleIscoRadius**2-3.0d0 &
            &*blackHoleIscoRadius+2.0d0*spinBlackHole*sqrt(blackHoleIscoRadius))
    end if

    ! Convert to physical units if necessary.
    if (unitsActual == unitsPhysical) Black_Hole_ISCO_Specific_Angular_Momentum=Black_Hole_ISCO_Specific_Angular_Momentum &
         &*sqrt(gravitationalConstant_internal*blackHole%mass()*Black_Hole_Gravitational_Radius(blackHole))
    return
  end function Black_Hole_ISCO_Specific_Angular_Momentum

  double precision function Black_Hole_Gravitational_Radius(blackHole)
    !!{
    Computes the gravitational radius (in Mpc) for the {\normalfont \ttfamily blackHole}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Prefixes    , only : milli
    implicit none
    class(nodeComponentBlackHole), intent(inout) :: blackHole

    Black_Hole_Gravitational_Radius=gravitationalConstant_internal*blackHole%mass()/(milli*speedLight)**2
    return
  end function Black_Hole_Gravitational_Radius

  double precision function Black_Hole_Frame_Dragging_Frequency_Node(blackHole,radius,units)
    !!{
    Returns the frame-dragging angular velocity in the Kerr metric.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (nodeComponentBlackHole), intent(inout), pointer  :: blackHole
    double precision                        , intent(in   )           :: radius
    integer                                 , intent(in   ), optional :: units
    integer                                                           :: unitsActual
    double precision                                                  :: spinBlackHole, radiusDimensionless

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get the dimensionless radius.
    select case (unitsActual)
    case (unitsGravitational)
       radiusDimensionless=radius
    case (unitsPhysical)
       radiusDimensionless=radius/Black_Hole_Gravitational_Radius(blackHole)
    case default
       radiusDimensionless=0.0d0
       call Error_Report('unrecognized units'//{introspection:location})
    end select

    ! Get the black hole spin.
    spinBlackHole=blackHole%spin()

    Black_Hole_Frame_Dragging_Frequency_Node=Black_Hole_Frame_Dragging_Frequency_Spin(spinBlackHole,radiusDimensionless)
    return
  end function Black_Hole_Frame_Dragging_Frequency_Node

  double precision function Black_Hole_Frame_Dragging_Frequency_Spin(spinBlackHole,radius)
    !!{
    Returns the frame-dragging angular velocity in the Kerr metric.
    !!}
    implicit none
    double precision, intent(in   ) :: spinBlackHole, radius

    Black_Hole_Frame_Dragging_Frequency_Spin=2.0d0*spinBlackHole/Black_Hole_Metric_A_Factor_Spin(spinBlackHole,radius)/radius**3
    return
  end function Black_Hole_Frame_Dragging_Frequency_Spin

  double precision function Black_Hole_Metric_A_Factor_Node(blackHole,radius,units)
    !!{
    Returns the $\mathcal{A}$ factor appearing in the Kerr metric for {\normalfont \ttfamily blackHole}.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (nodeComponentBlackHole), intent(inout)           :: blackHole
    double precision                        , intent(in   )           :: radius
    integer                                 , intent(in   ), optional :: units
    integer                                                           :: unitsActual
    double precision                                                  :: spinBlackHole, radiusDimensionless

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get dimensionless radius.
    select case (unitsActual)
    case (unitsGravitational)
       radiusDimensionless=radius
    case (unitsPhysical)
       radiusDimensionless=radius/Black_Hole_Gravitational_Radius(blackHole)
    case default
       radiusDimensionless=0.0d0
       call Error_Report('unrecognized units'//{introspection:location})
    end select

    ! Get the black hole spin.
    spinBlackHole=blackHole%spin()

    Black_Hole_Metric_A_Factor_Node=Black_Hole_Metric_A_Factor_Spin(spinBlackHole,radiusDimensionless)
    return
  end function Black_Hole_Metric_A_Factor_Node

  double precision function Black_Hole_Metric_A_Factor_Spin(spinBlackHole,radius)
    !!{
    Returns the $\mathcal{A}$ factor appearing in the Kerr metric for spin {\normalfont \ttfamily spinBlackHole}.
    !!}
    implicit none
    double precision, intent(in   ) :: spinBlackHole, radius

    Black_Hole_Metric_A_Factor_Spin=1.0d0+spinBlackHole/radius**2+2.0d0*spinBlackHole**2/radius**3
    return
  end function Black_Hole_Metric_A_Factor_Spin

  double precision function Black_Hole_Metric_D_Factor_Node(blackHole,radius,units)
    !!{
    Returns the $\mathcal{D}$ factor appearing in the Kerr metric for {\normalfont \ttfamily blackHole}.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (nodeComponentBlackHole), intent(inout), pointer  :: blackHole
    double precision                        , intent(in   )           :: radius
    integer                                 , intent(in   ), optional :: units
    integer                                                           :: unitsActual
    double precision                                                  :: spinBlackHole, radiusDimensionless

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get dimensionless radius.
    select case (unitsActual)
    case (unitsGravitational)
       radiusDimensionless=radius
    case (unitsPhysical)
       radiusDimensionless=radius/Black_Hole_Gravitational_Radius(blackHole)
    case default
       radiusDimensionless=0.0d0
       call Error_Report('unrecognized units'//{introspection:location})
    end select

    ! Get the black hole spin.
    spinBlackHole=blackHole%spin()

    Black_Hole_Metric_D_Factor_Node=Black_Hole_Metric_D_Factor_Spin(spinBlackHole,radiusDimensionless)
    return
  end function Black_Hole_Metric_D_Factor_Node

  double precision function Black_Hole_Metric_D_Factor_Spin(spinBlackHole,radius)
    !!{
    Returns the $\mathcal{D}$ factor appearing in the Kerr metric for spin {\normalfont \ttfamily spinBlackHole}.
    !!}
    implicit none
    double precision, intent(in   ) :: spinBlackHole, radius

    Black_Hole_Metric_D_Factor_Spin=1.0d0-2.0d0/radius+(spinBlackHole/radius)**2
    return
  end function Black_Hole_Metric_D_Factor_Spin

  double precision function Black_Hole_Horizon_Radius_Node(blackHole,units)
    !!{
    Return the radius of the horizon for a Kerr metric with dimensionless angular momentum {\normalfont \ttfamily j}.
    The radius is in units of the gravitational radius.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (nodeComponentBlackHole), intent(inout), pointer  :: blackHole
    integer                                 , intent(in   ), optional :: units
    double precision                                                  :: spinBlackHole, radiusDimensionless
    integer                                                           :: unitsActual

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get the black hole spin.
    spinBlackHole=blackHole%spin()

    radiusDimensionless=Black_Hole_Horizon_Radius_Spin(spinBlackHole)
    select case (unitsActual)
    case (unitsGravitational)
       Black_Hole_Horizon_Radius_Node=radiusDimensionless
    case (unitsPhysical)
       Black_Hole_Horizon_Radius_Node=radiusDimensionless*Black_Hole_Gravitational_Radius(blackHole)
    case default
       Black_Hole_Horizon_Radius_Node=0.0d0
       call Error_Report('unrecognized units'//{introspection:location})
    end select
    return
  end function Black_Hole_Horizon_Radius_Node

  double precision function Black_Hole_Horizon_Radius_Spin(spinBlackHole)
    !!{
    Return the radius of the horizon for a Kerr metric with dimensionless angular momentum {\normalfont \ttfamily j}.
    The radius is in units of the gravitational radius.
    !!}
    implicit none
    double precision, intent(in   ) :: spinBlackHole

    Black_Hole_Horizon_Radius_Spin=1.0d0+sqrt(1.0d0-spinBlackHole**2)
    return
  end function Black_Hole_Horizon_Radius_Spin

  double precision function Black_Hole_Static_Radius_Node(blackHole,theta,units)
    !!{
    Return the radius of the static limit for a Kerr metric for the black hole in {\normalfont \ttfamily blackHole} and angle {\normalfont \ttfamily theta}.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (nodeComponentBlackHole), intent(inout)           :: blackHole
    integer                                 , intent(in   ), optional :: units
    double precision                        , intent(in   ), optional :: theta
    double precision                                                  :: spinBlackHole, radiusDimensionless
    integer                                                           :: unitsActual

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get the black hole spin.
    spinBlackHole=blackHole%spin()

    ! Get the dimensionless static radius.
    radiusDimensionless=Black_Hole_Static_Radius_Spin(spinBlackHole,theta)

    ! Convert to the appropriate units.
    select case (unitsActual)
    case (unitsGravitational)
       Black_Hole_Static_Radius_Node=radiusDimensionless
    case (unitsPhysical)
       Black_Hole_Static_Radius_Node=radiusDimensionless*Black_Hole_Gravitational_Radius(blackHole)
    case default
       Black_Hole_Static_Radius_Node=0.0d0
       call Error_Report('unrecognized units'//{introspection:location})
    end select
    return
  end function Black_Hole_Static_Radius_Node

  double precision function Black_Hole_Static_Radius_Spin(spinBlackHole,theta)
    !!{
    Return the radius of the static limit for a Kerr metric for a black hole of given {\normalfont \ttfamily spinBlackHole} and angle {\normalfont \ttfamily theta}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   )           :: spinBlackHole
    double precision, intent(in   ), optional :: theta
    double precision                          :: thetaActual

    ! Determine what angle to use.
    if (present(theta)) then
       thetaActual=theta
    else
       ! Assume rotational plane if no angle specified.
       thetaActual=Pi/2.0d0
    end if

    Black_Hole_Static_Radius_Spin=1.0d0+sqrt(1.0d0-(spinBlackHole*cos(thetaActual))**2)
    return
  end function Black_Hole_Static_Radius_Spin

  double precision function A1(spinBlackHole)
    !!{
    Return the function $A_1(j)$ that appears in the Kerr metric with spin {\normalfont \ttfamily spinBlackHole}.
    !!}
    implicit none
    double precision, intent(in   ) :: spinBlackHole

    A1=1.0d0+((1.0d0-spinBlackHole**2)**(1.0d0/3.0d0))*((1.0d0+spinBlackHole)**(1.0d0/3.0d0)+(1.0d0-spinBlackHole)**(1.0d0/3.0d0))
    return
  end function A1

  double precision function A2(spinBlackHole)
    !!{
    Return the function $A_2(j)$ that appears in the Kerr metric with spin {\normalfont \ttfamily spinBlackHole}.
    !!}
    implicit none
    double precision, intent(in   ) :: spinBlackHole

    A2=sqrt(3.0d0*spinBlackHole**2+A1(spinBlackHole)**2)
    return
  end function A2

  double precision function Black_Hole_Rotational_Energy_Spin_Down_Node(blackHole)
    !!{
    Computes the spin down rate of a black hole due to extraction of rotational energy for the primary black hole.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                                        :: spinBlackHole

    ! Get the black hole spin.
    spinBlackHole=blackHole%spin()

    ! Get the spin down function.
    Black_Hole_Rotational_Energy_Spin_Down_Node=Black_Hole_Rotational_Energy_Spin_Down_Spin(spinBlackHole)
    return
  end function Black_Hole_Rotational_Energy_Spin_Down_Node

  double precision function Black_Hole_Rotational_Energy_Spin_Down_Spin(spinBlackHole)
    !!{
    Computes the spin down rate of a black hole due to extraction of rotational energy. Specifically, it returns the factor $S$
    in the relation:
    \begin{equation}
    s = - S {P_\mathrm{rotation} \over \dot M_{\bullet, 0} \clight^2},
    \end{equation}
    where $P_\mathrm{rotation}$ is the power of rotational energy extraction and
    \begin{equation}
    S = [(1+\sqrt{1-j^2})^2+j^2] {\sqrt{1-j^2}\over j},
    \end{equation}
    for black hole spin $j$. This result is derived as follows. Starting from equation~(9) in \cite{benson_maximum_2009}:
    \begin{equation}
     M_{\bullet,\mathrm{irr}} = \frac{1}{2} M_\bullet \left[ (1+\sqrt{1-j^2})^2 + j^2 \right]^{1/2},
     \label{eq:massBlackHoleIrreducible}
    \end{equation}
    which we can rearrange to get
    \begin{equation}
     (1 + \sqrt{1-j^2})^2 + j^2 = 4 (M_{\bullet,\mathrm{irr}}/M_\bullet)^2.
    \end{equation}
    Differentiating this gives
    \begin{equation}
    \left[ -2(1 + \sqrt{1-j^2}) \frac{1}{2} (1-j^2)^{-1/2} 2 j + 2 j \right] \frac{\mathrm{d} j}{\mathrm{d} M_\bullet} = - 8 (M_{\bullet,\mathrm{irr}}/M_\bullet)^2 M^{-1}_\bullet,
    \end{equation}
    which simplifies to
    \begin{equation}
    \frac{j}{(1-j^2)^{1/2}} \frac{\mathrm{d} j}{\mathrm{d} M_\bullet} = 4 (M_{\bullet,\mathrm{irr}}/M_\bullet)^2 M^{-1}_\bullet,
    \end{equation}
    and therefore
    \begin{equation}
     \frac{\mathrm{d} j}{\mathrm{d} M_\bullet} = \frac{4}{M_\bullet} \left(\frac{M_{\bullet,\mathrm{irr}}}{M_\bullet}\right)^2  \frac{(1-j^2)^{1/2}}{j},
    \end{equation}
    which is equation~(10) in \cite{benson_maximum_2009}.
    
    The spin-down rate is then
    \begin{equation}
     \frac{\mathrm{d}j}{\mathrm{d}t} = - \frac{\mathrm{d}j}{\mathrm{d}M_\bullet} \frac{P_\mathrm{jet}}{\mathrm{c}^2} = - 4 \frac{P_\mathrm{jet}}{M_\bullet \mathrm{c}^2} \left(\frac{M_{\bullet,\mathrm{irr}}}{M_\bullet}\right)^2  \frac{(1-j^2)^{1/2}}{j}.
     \label{eq:blackHoleJetSpinDownRate}
    \end{equation}
    
    The spin-down parameter is then
    \begin{equation}
     s_\mathrm{jet} = \frac{M_\bullet}{\dot{M}_{\bullet,0}} \frac{\mathrm{d}j}{\mathrm{d}t}.
    \end{equation}
    Using equations~(\ref{eq:massBlackHoleIrreducible}) and (\ref{eq:blackHoleJetSpinDownRate}) the above becomes
    \begin{equation}
      s_\mathrm{jet} = - \frac{P_\mathrm{jet}}{\dot{M}_{\bullet,0} \mathrm{c}^2} \left[ (1+\sqrt{1-j^2})^2 + j^2 \right] \frac{(1-j^2)^{1/2}}{j},
    \end{equation}
    which is equation~(12)\footnote{Note that equation~(12) of \cite{benson_maximum_2009} is missing the factor of
    $\dot{M}_{\bullet,0} \mathrm{c}^2$ in the denominator} of \cite{benson_maximum_2009}.
    !!}
    implicit none
    double precision, intent(in   ) :: spinBlackHole
    double precision, parameter     :: spinBlackHoleMinimum      =5.0d-8
    double precision, parameter     :: spinBlackHoleSeriesMinimum=1.0d-20

    if (spinBlackHole > spinBlackHoleMinimum) then
       ! Full solution.
       Black_Hole_Rotational_Energy_Spin_Down_Spin=((1.0d0+sqrt(1.0d0-spinBlackHole**2))**2+spinBlackHole**2)*sqrt(1.0d0&
            &-spinBlackHole**2)/spinBlackHole
    else
       if (spinBlackHole > spinBlackHoleSeriesMinimum) then
          ! Series solution.
          Black_Hole_Rotational_Energy_Spin_Down_Spin=4.0d0/spinBlackHole-3.0d0*spinBlackHole-0.25d0*(spinBlackHole**3)
       else
          ! Series diverges as spin tends to zero. Simply set to zero below some low spin. Rotational energy is zero close to zero
          ! spin so this should not cause a problem.
          Black_Hole_Rotational_Energy_Spin_Down_Spin=0.0d0
       end if
    end if
    return
  end function Black_Hole_Rotational_Energy_Spin_Down_Spin

end module Black_Hole_Fundamentals

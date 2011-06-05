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


!% Contains a module which implements fundamental properties of black holes.

module Black_Hole_Fundamentals
  !% Implements fundamental properties of black holes.
  use Tree_Nodes
  private
  public :: Black_Hole_ISCO_Radius, Black_Hole_ISCO_Specific_Energy, Black_Hole_Gravitational_Radius,&
       & Black_Hole_Frame_Dragging_Frequency, Black_Hole_Metric_A_Factor, Black_Hole_Metric_D_Factor, Black_Hole_Horizon_Radius,&
       & Black_Hole_Static_Radius, Black_Hole_ISCO_Specific_Angular_Momentum, Black_Hole_Rotational_Energy_Spin_Down,&
       & Black_Hole_Eddington_Accretion_Rate, Black_Hole_Fundamentals_Unit_Test

  ! Identifiers for unit system options.
  integer, parameter, public :: unitsGravitational=0
  integer, parameter, public :: unitsPhysical     =1

  ! Identifiers for orbit orientation options.
  integer, parameter, public :: orbitPrograde  =0
  integer, parameter, public :: orbitRetrograde=1

  ! Generic interfaces for functions.
  interface Black_Hole_ISCO_Radius
     module procedure Black_Hole_ISCO_Radius_Node
     module procedure Black_Hole_ISCO_Radius_Spin
  end interface
   interface Black_Hole_Rotational_Energy_Spin_Down
     module procedure Black_Hole_Rotational_Energy_Spin_Down_Node
     module procedure Black_Hole_Rotational_Energy_Spin_Down_Spin
  end interface
  interface Black_Hole_Metric_A_Factor
     module procedure Black_Hole_Metric_A_Factor_Node
     module procedure Black_Hole_Metric_A_Factor_Spin
  end interface
  interface Black_Hole_Metric_D_Factor
     module procedure Black_Hole_Metric_D_Factor_Node
     module procedure Black_Hole_Metric_D_Factor_Spin
  end interface
  interface Black_Hole_Frame_Dragging_Frequency
     module procedure Black_Hole_Frame_Dragging_Frequency_Node
     module procedure Black_Hole_Frame_Dragging_Frequency_Spin
  end interface
  interface Black_Hole_Horizon_Radius
     module procedure Black_Hole_Horizon_Radius_Node
     module procedure Black_Hole_Horizon_Radius_Spin
  end interface
  interface Black_Hole_Static_Radius
     module procedure Black_Hole_Static_Radius_Node
     module procedure Black_Hole_Static_Radius_Spin
  end interface
  interface Black_Hole_ISCO_Specific_Energy
     module procedure Black_Hole_ISCO_Specific_Energy_Node
     module procedure Black_Hole_ISCO_Specific_Energy_Spin
  end interface Black_Hole_ISCO_Specific_Energy

contains

  double precision function Black_Hole_ISCO_Radius_Spin(blackHoleSpin,orbit)
    !% Returns the radius (in gravitational units and for a prograde or retorgrade orbit) of the innermost stable
    !% circular orbit for a black hole with spin {\tt blackHoleSpin}.
    use Galacticus_Error
    implicit none
    double precision, intent(in)           :: blackHoleSpin
    integer,          intent(in), optional :: orbit
    integer,          save                 :: orbitPrevious
    double precision, save                 :: blackHoleSpinPrevious=2.0d0, radiusISCOPrevious
    !$omp threadprivate(orbitPrevious,blackHoleSpinPrevious,radiusISCOPrevious)
    integer                                :: orbitActual
    double precision                       :: A1Factor,A2Factor
 
    ! Determine what orbit to use.
    if (present(orbit)) then
       orbitActual=orbit
    else
       orbitActual=orbitPrograde
    end if

    ! Check if we're being called with the same spin and orbit as the last time.
    if (orbitActual == orbitPrevious .and. blackHoleSpin == blackHoleSpinPrevious) then
       ! We are, so just return our previously computed radius.
       Black_Hole_ISCO_Radius_Spin=radiusISCOPrevious
    else
       ! We are not, so compute the radius of the innermost stable circular orbit in units of the gravitational radius.
       A1Factor=A1(blackHoleSpin)
       A2Factor=A2(blackHoleSpin)
       select case (orbitActual)
       case (orbitPrograde)
          Black_Hole_ISCO_Radius_Spin=3.0d0+A2Factor-dsqrt((3.0d0-A1Factor)*(3.0d0+A1Factor+2.0d0*A2Factor))
       case (orbitRetrograde)
          Black_Hole_ISCO_Radius_Spin=3.0d0+A2Factor+dsqrt((3.0d0-A1Factor)*(3.0d0+A1Factor+2.0d0*A2Factor))
       case default
          call Galacticus_Error_Report('Black_Hole_ISCO_Radius_Spin','unrecognized orbit parameter')
       end select
       orbitPrevious        =orbitActual
       blackHoleSpinPrevious=blackHoleSpin
       radiusISCOPrevious   =Black_Hole_ISCO_Radius_Spin
    end if
    return
  end function Black_Hole_ISCO_Radius_Spin

  double precision function Black_Hole_Eddington_Accretion_Rate(thisNode)
    !% Return the Eddington accretion rate (in $M_\odot$ Gyr$^{-1}$) for the black hole in {\tt thisNode}.
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Black_Hole_Eddington_Accretion_Rate=4.0d0*Pi*gravitationalConstant*Tree_Node_Black_Hole_Mass(thisNode)*massHydrogenAtom&
         &*gigaYear/thomsonCrossSection/speedLight

    return
  end function Black_Hole_Eddington_Accretion_Rate

  double precision function Black_Hole_ISCO_Radius_Node(thisNode,units,orbit)
    !% Returns the radius (in physical or gravitational units and for a prograde or retorgrade orbit) of the innermost stable
    !% circular orbit for the black hole in {\tt thisNode}.
    use Galacticus_Error
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(in),    optional :: units,orbit
    integer                                 :: unitsActual,orbitActual
    double precision                        :: blackHoleSpin
 
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
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    ! Compute the radius of the innermost stable circular orbit in units of the gravitational radius.
    Black_Hole_ISCO_Radius_Node=Black_Hole_ISCO_Radius_Spin(blackHoleSpin,orbitActual)

    ! Convert to physical units if necessary.
    if (unitsActual == unitsPhysical) Black_Hole_ISCO_Radius_Node=Black_Hole_ISCO_Radius_Node*Black_Hole_Gravitational_Radius(thisNode)
    return
  end function Black_Hole_ISCO_Radius_Node

  double precision function Black_Hole_ISCO_Specific_Energy_Node(thisNode,units,orbit)
    !% Returns the specific energy (in physical or gravitational units and for a prograde or retorgrade orbit) of the innermost
    !% stable circular orbit for the black hole in {\tt thisNode}.
    use Galacticus_Error
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: units,orbit
    integer                                   :: unitsActual
    double precision                          :: blackHoleSpin,blackHoleIscoRadius
 
    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    ! Get the dimensionless ISCO energy.
    Black_Hole_ISCO_Specific_Energy_Node=Black_Hole_ISCO_Specific_Energy_Spin(blackHoleSpin,orbit)

    ! Convert to physical units if necessary.
    if (unitsActual == unitsPhysical) Black_Hole_ISCO_Specific_Energy_Node=Black_Hole_ISCO_Specific_Energy_Node&
         &*gravitationalConstantGalacticus*Tree_Node_Black_Hole_Mass(thisNode)**2/Black_Hole_Gravitational_Radius(thisNode)
    return
  end function Black_Hole_ISCO_Specific_Energy_Node

  double precision function Black_Hole_ISCO_Specific_Energy_Spin(blackHoleSpin,orbit)
    !% Returns the specific energy (in physical or gravitational units and for a prograde or retorgrade orbit) of the innermost
    !% stable circular orbit for a black hole of given {\tt blackHoleSpin}.
    implicit none
    double precision, intent(inout)          :: blackHoleSpin
    integer,          intent(in),   optional :: orbit
    ! Maximum spin above which we use series solution.
    double precision, parameter              :: maximumSpin=0.99999d0
    ! Coefficients in the expansion of ISCO energy with spin close to maximal spin.
    double precision, parameter              :: coefficientZeroth=0.5773502693d0, coefficientFirst=0.9164864242d0
    double precision                         :: blackHoleIscoRadius
 
    ! Get black hole ISCO radius.
    blackHoleIscoRadius=Black_Hole_ISCO_Radius(blackHoleSpin,orbit)

    ! Compute the specific energy in gravitational units.
    if (blackHoleSpin >= maximumSpin) then
       Black_Hole_ISCO_Specific_Energy_Spin=coefficientZeroth+coefficientFirst*(1.0d0-blackHoleSpin)**(1.0d0/3.0d0)
    else
       Black_Hole_ISCO_Specific_Energy_Spin=(blackHoleIscoRadius**2-2.0d0*blackHoleIscoRadius+blackHoleSpin &
            &*dsqrt(blackHoleIscoRadius))/blackHoleIscoRadius/dsqrt(blackHoleIscoRadius**2-3.0d0*blackHoleIscoRadius+2.0d0 &
            &*blackHoleSpin*dsqrt(blackHoleIscoRadius))
    end if

    return
  end function Black_Hole_ISCO_Specific_Energy_Spin

  double precision function Black_Hole_ISCO_Specific_Angular_Momentum(thisNode,units,orbit)
    !% Returns the specific angular momentum (in physical or gravitational units and for a prograde or retorgrade orbit) of the
    !% innermost stable circular orbit for the black hole in {\tt thisNode}.
    use Galacticus_Error
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: units,orbit
    ! Maximum spin above which we use series solution.
    double precision, parameter               :: maximumSpin=0.99999d0
    ! Coefficients in the expansion of ISCO energy with spin close to maximal spin.
    double precision, parameter               :: coefficientZeroth=1.154700538d0, coefficientFirst=1.832972849d0
    integer                                   :: unitsActual,orbitActual
    double precision                          :: blackHoleSpin,blackHoleIscoRadius
 
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
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)
    blackHoleIscoRadius=Black_Hole_ISCO_Radius(thisNode,units=unitsGravitational,orbit=orbitActual)

    ! Compute the specific angular momentum in gravitational units.
    if (blackHoleSpin > maximumSpin) then
       Black_Hole_ISCO_Specific_Angular_Momentum=coefficientZeroth+coefficientFirst*(1.0d0-blackHoleSpin)**(1.0d0/3.0d0)
    else
       Black_Hole_ISCO_Specific_Angular_Momentum=dsqrt(blackHoleIscoRadius)*(blackHoleIscoRadius**2-2.0d0*blackHoleSpin &
            &*dsqrt(blackHoleIscoRadius)+blackHoleSpin**2)/blackHoleIscoRadius/dsqrt(blackHoleIscoRadius**2-3.0d0 &
            &*blackHoleIscoRadius+2.0d0*blackHoleSpin*dsqrt(blackHoleIscoRadius))
    end if

    ! Convert to physical units if necessary.
    if (unitsActual == unitsPhysical) Black_Hole_ISCO_Specific_Angular_Momentum=Black_Hole_ISCO_Specific_Angular_Momentum &
         &*dsqrt(gravitationalConstantGalacticus*Tree_Node_Black_Hole_Mass(thisNode)*Black_Hole_Gravitational_Radius(thisNode))
    return
  end function Black_Hole_ISCO_Specific_Angular_Momentum

  double precision function Black_Hole_Gravitational_Radius(thisNode)
    !% Computes the gravitational radius (in Mpc) for the black hole in {\tt thisNode}.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Physical
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Black_Hole_Gravitational_Radius=gravitationalConstantGalacticus*Tree_Node_Black_Hole_Mass(thisNode)/(milli*speedLight)**2
    return
  end function Black_Hole_Gravitational_Radius

  double precision function Black_Hole_Frame_Dragging_Frequency_Node(thisNode,radius,units)
    !% Returns the frame-dragging angular velocity in the Kerr metric.
    use Galacticus_Error
    implicit none
    type(treeNode),   intent(inout), pointer    :: thisNode
    double precision, intent(in)                :: radius
    integer,          intent(in),    optional   :: units
    integer                                     :: unitsActual
    double precision                            :: blackHoleSpin,radiusDimensionless

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
       radiusDimensionless=radius/Black_Hole_Gravitational_Radius(thisNode)
    case default
       call Galacticus_Error_Report('Black_Hole_Frame_Dragging_Frequency_Spin','unrecognized units')
    end select

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    Black_Hole_Frame_Dragging_Frequency_Node=Black_Hole_Frame_Dragging_Frequency_Spin(blackHoleSpin,radiusDimensionless)
    return
  end function Black_Hole_Frame_Dragging_Frequency_Node

  double precision function Black_Hole_Frame_Dragging_Frequency_Spin(blackHoleSpin,radius)
    !% Returns the frame-dragging angular velocity in the Kerr metric.
    use Galacticus_Error
    implicit none
    double precision, intent(in) :: radius,blackHoleSpin
 
    Black_Hole_Frame_Dragging_Frequency_Spin=2.0d0*blackHoleSpin/Black_Hole_Metric_A_Factor_Spin(blackHoleSpin,radius)/radius**3
    return
  end function Black_Hole_Frame_Dragging_Frequency_Spin

  double precision function Black_Hole_Metric_A_Factor_Node(thisNode,radius,units)
    !% Returns the $\mathcal{A}$ factor appearing in the Kerr metric for {\tt thisNode}.
    use Galacticus_Error
    implicit none
    type(treeNode),   intent(inout), pointer    :: thisNode
    double precision, intent(in)                :: radius
    integer,          intent(in),    optional   :: units
    integer                                     :: unitsActual
    double precision                            :: blackHoleSpin,radiusDimensionless

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
       radiusDimensionless=radius/Black_Hole_Gravitational_Radius(thisNode)
    case default
       call Galacticus_Error_Report('Black_Hole_Metric_A_Factor_Spin','unrecognized units')
    end select

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    Black_Hole_Metric_A_Factor_Node=Black_Hole_Metric_A_Factor_Spin(blackHoleSpin,radiusDimensionless)
    return
  end function Black_Hole_Metric_A_Factor_Node

  double precision function Black_Hole_Metric_A_Factor_Spin(blackHoleSpin,radius)
    !% Returns the $\mathcal{A}$ factor appearing in the Kerr metric for spin {\tt blackHoleSpin}.
    use Galacticus_Error
    implicit none
    double precision, intent(in)                :: radius,blackHoleSpin

    Black_Hole_Metric_A_Factor_Spin=1.0d0+blackHoleSpin/radius**2+2.0d0*blackHoleSpin**2/radius**3
    return
  end function Black_Hole_Metric_A_Factor_Spin

  double precision function Black_Hole_Metric_D_Factor_Node(thisNode,radius,units)
    !% Returns the $\mathcal{D}$ factor appearing in the Kerr metric for {\tt thisNode}.
    use Galacticus_Error
    implicit none
    type(treeNode),   intent(inout), pointer    :: thisNode
    double precision, intent(in)                :: radius
    integer,          intent(in),    optional   :: units
    integer                                     :: unitsActual
    double precision                            :: blackHoleSpin,radiusDimensionless

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
       radiusDimensionless=radius/Black_Hole_Gravitational_Radius(thisNode)
    case default
       call Galacticus_Error_Report('Black_Hole_Metric_D_Factor_Spin','unrecognized units')
    end select

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    Black_Hole_Metric_D_Factor_Node=Black_Hole_Metric_D_Factor_Spin(blackHoleSpin,radiusDimensionless)
    return
  end function Black_Hole_Metric_D_Factor_Node

  double precision function Black_Hole_Metric_D_Factor_Spin(blackHoleSpin,radius)
    !% Returns the $\mathcal{D}$ factor appearing in the Kerr metric for spin {\tt blackHoleSpin}.
    use Galacticus_Error
    implicit none
    double precision, intent(in)                :: radius,blackHoleSpin

    Black_Hole_Metric_D_Factor_Spin=1.0d0-2.0d0/radius+(blackHoleSpin/radius)**2
    return
  end function Black_Hole_Metric_D_Factor_Spin

  double precision function Black_Hole_Horizon_Radius_Node(thisNode,units)
    !% Return the radius of the horizon for a Kerr metric with dimensionless angular momentum {\tt j}.
    !% The radius is in units of the gravitational radius.
    use Galacticus_Error
    implicit none
    type(treeNode),  intent(inout), pointer  :: thisNode
    integer,         intent(in),    optional :: units
    double precision                         :: blackHoleSpin,radiusDimensionless
    integer                                  :: unitsActual

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    radiusDimensionless=Black_Hole_Horizon_Radius_Spin(blackHoleSpin)
    select case (unitsActual)
    case (unitsGravitational)
       Black_Hole_Horizon_Radius_Node=radiusDimensionless
    case (unitsPhysical)
       Black_Hole_Horizon_Radius_Node=radiusDimensionless*Black_Hole_Gravitational_Radius(thisNode)
    case default
       call Galacticus_Error_Report('Black_Hole_Horizon_Radius_Node','unrecognized units')
    end select
    return
  end function Black_Hole_Horizon_Radius_Node

  double precision function Black_Hole_Horizon_Radius_Spin(blackHoleSpin)
    !% Return the radius of the horizon for a Kerr metric with dimensionless angular momentum {\tt j}.
    !% The radius is in units of the gravitational radius.
    use Galacticus_Error
    implicit none
    double precision, intent(in) :: blackHoleSpin

    Black_Hole_Horizon_Radius_Spin=1.0d0+dsqrt(1.0d0-blackHoleSpin**2)   
    return
  end function Black_Hole_Horizon_Radius_Spin
  
  double precision function Black_Hole_Static_Radius_Node(thisNode,theta,units)
    !% Return the radius of the static limit for a Kerr metric for the black hole in {\tt thisNode} and angle {\tt theta}.
    use Galacticus_Error
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: units
    double precision, intent(in),    optional :: theta
    double precision                          :: blackHoleSpin,radiusDimensionless
    integer                                   :: unitsActual

    ! Determine what system of units to use.
    if (present(units)) then
       unitsActual=units
    else
       unitsActual=unitsPhysical
    end if

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    ! Get the dimensionless static radius.
    radiusDimensionless=Black_Hole_Static_Radius_Spin(blackHoleSpin,theta)
 
    ! Convert to the appropriate units.
    select case (unitsActual)
    case (unitsGravitational)
       Black_Hole_Static_Radius_Node=radiusDimensionless
    case (unitsPhysical)
       Black_Hole_Static_Radius_Node=radiusDimensionless*Black_Hole_Gravitational_Radius(thisNode)
    case default
       call Galacticus_Error_Report('Black_Hole_Static_Radius_Node','unrecognized units')
    end select
    return
  end function Black_Hole_Static_Radius_Node

  double precision function Black_Hole_Static_Radius_Spin(blackHoleSpin,theta)
    !% Return the radius of the static limit for a Kerr metric for a black hole of given {\tt blackHoleSpin} and angle {\tt theta}.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in)           :: blackHoleSpin
    double precision, intent(in), optional :: theta
    double precision                       :: thetaActual

    ! Determine what angle to use.
    if (present(theta)) then
       thetaActual=theta
    else
       ! Assume rotational plane if no angle specified.
       thetaActual=Pi/2.0d0
    end if

    Black_Hole_Static_Radius_Spin=1.0d0+dsqrt(1.0d0-(blackHoleSpin*dcos(thetaActual))**2)
    return
  end function Black_Hole_Static_Radius_Spin

  double precision function A1(blackHoleSpin)
    !% Return the function $A_1(j)$ that appears in the Kerr metric with spin {\tt blackHoleSpin}.
    implicit none
    double precision, intent(in) :: blackHoleSpin

    A1=1.0d0+((1.0d0-blackHoleSpin**2)**(1.0d0/3.0d0))*((1.0d0+blackHoleSpin)**(1.0d0/3.0d0)+(1.0d0-blackHoleSpin)**(1.0d0/3.0d0))
    return
  end function A1

  double precision function A2(blackHoleSpin)
    !% Return the function $A_2(j)$ that appears in the Kerr metric with spin {\tt blackHoleSpin}.
    implicit none
    double precision, intent(in)           :: blackHoleSpin

    A2=dsqrt(3.0d0*blackHoleSpin**2+A1(blackHoleSpin)**2)
    return
  end function A2

  double precision function Black_Hole_Rotational_Energy_Spin_Down_Node(thisNode)
    !% Wrapper function for \href{func:black_hole_rotational_energy_spin_down_spin}{{\tt
    !% Black\_Hole\_Rotational\_Energy\_Spin\_Down\_Node}} which takes a tree node as input.
    implicit none
    type(treeNode),  intent(inout), pointer :: thisNode
    double precision                        :: blackHoleSpin

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    ! Get the spin down function.
    Black_Hole_Rotational_Energy_Spin_Down_Node=Black_Hole_Rotational_Energy_Spin_Down_Spin(blackHoleSpin)
    return
  end function Black_Hole_Rotational_Energy_Spin_Down_Node

  double precision function Black_Hole_Rotational_Energy_Spin_Down_Spin(blackHoleSpin)
    !% Computes the spin down rate of a black hole due to extraction of rotational energy. Specifically, it returns the factor $S$
    !% in the relation:
    !% \begin{equation}
    !% s = - S {P_{\rm rotation} \over \dot M_{\bullet, 0} \clight^2},
    !% \end{equation}
    !% where $P_{\rm rotation}$ is the power of rotational energy extraction and
    !% \begin{equation}
    !% S = [(1+\sqrt{1-j^2})^2+j^2] {\sqrt{1-j^2}\over j},
    !% \end{equation}
    !% for black hole spin $j$.
    implicit none
    double precision, intent(in) :: blackHoleSpin
    double precision, parameter  :: blackHoleSpinMinimum      =5.0d-8
    double precision, parameter  :: blackHoleSpinSeriesMinimum=1.0d-20

    if (blackHoleSpin > blackHoleSpinMinimum) then
       ! Full solution.
       Black_Hole_Rotational_Energy_Spin_Down_Spin=((1.0d0+dsqrt(1.0d0-blackHoleSpin**2))**2+blackHoleSpin**2)*dsqrt(1.0d0&
            &-blackHoleSpin**2)/blackHoleSpin
    else
       if (blackHoleSpin > blackHoleSpinSeriesMinimum) then
          ! Series solution.
          Black_Hole_Rotational_Energy_Spin_Down_Spin=4.0d0/blackHoleSpin-3.0d0*blackHoleSpin-0.25d0*(blackHoleSpin**3)
       else
          ! Series diverges as spin tends to zero. Simply set to zero below some low spin. Rotational energy is zero close to zero
          ! spin so this should not cause a problem.
          Black_Hole_Rotational_Energy_Spin_Down_Spin=0.0d0
       end if
    end if
    return
  end function Black_Hole_Rotational_Energy_Spin_Down_Spin

  !# <galacticusSelfTest>
  !#   <unitName>Black_Hole_Fundamentals_Unit_Test</unitName>
  !# </galacticusSelfTest>
  subroutine Black_Hole_Fundamentals_Unit_Test
    !% Unit tests for the black holes fundamentals module.
    use Unit_Tests
    implicit none

    ! Begin a unit testing group.
    call Unit_Tests_Begin_Group("black hole fundamentals")

    ! ISCO radius for a Schwarzchild black hole should be 6 for prograde orbits.
    call Assert("Schwarzchild metric ISCO radius",Black_Hole_ISCO_Radius_Spin(0.0d0,orbitPrograde),6.0d0,compareEquals,1.0d-6)

    ! ISCO radius for an extreme Kerr black hole should be 1 for prograde orbits.
    call Assert("Extreme Kerr metric ISCO radius",Black_Hole_ISCO_Radius_Spin(1.0d0,orbitPrograde),1.0d0,compareEquals,1.0d-6)

    ! Horizon radius for a Schwarzchild black hole should be 2.
    call Assert("Schwarzchild metric horizon radius",Black_Hole_Horizon_Radius_Spin(0.0d0),2.0d0,compareEquals,1.0d-6)

    ! Horizon radius for an extreme Kerr black hole should be 1.
    call Assert("Extreme Kerr metric horizon radius",Black_Hole_Horizon_Radius_Spin(1.0d0),1.0d0,compareEquals,1.0d-6)

    ! End the unit testing group.
    call Unit_Tests_End_Group
    return
  end subroutine Black_Hole_Fundamentals_Unit_Test

end module Black_Hole_Fundamentals

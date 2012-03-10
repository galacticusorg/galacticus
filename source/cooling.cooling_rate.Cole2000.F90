!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements a \cite{cole_hierarchical_2000} cooling rate calculation.

module Cooling_Rates_Cole2000
  !% Implements a \cite{cole_hierarchical_2000} cooling rate calculation.
  use, intrinsic :: ISO_C_Binding
  use Tree_Nodes
  implicit none
  private
  public :: Cooling_Rate_Cole2000_Initialize

contains

  !# <coolingRateMethod>
  !#  <unitName>Cooling_Rate_Cole2000_Initialize</unitName>
  !# </coolingRateMethod>
  subroutine Cooling_Rate_Cole2000_Initialize(coolingRateMethod,Cooling_Rate_Get)
    !% Initializes the ``Cole et al. (2000)'' cooling rate module.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: coolingRateMethod
    procedure(double precision), pointer, intent(inout) :: Cooling_Rate_Get
    
    ! Return a pointer to our implementation.
    if (coolingRateMethod == 'Cole2000') Cooling_Rate_Get => Cooling_Rate_Cole2000
    return
  end subroutine Cooling_Rate_Cole2000_Initialize

  double precision function Cooling_Rate_Cole2000(thisNode)
    !% Computes the mass cooling rate in a hot gas halo utilizing the \cite{cole_hierarchical_2000} method. This is based on the
    !% properties of the halo at formation time, and gives a zero cooling rate when the cooling radius exceeds the virial radius.
    use Tree_Nodes
    use Cooling_Infall_Radii
    use Numerical_Constants_Math
    use Hot_Halo_Density_Profile
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: infallRadius,infallDensity,outerRadius,infallRadiusGrowthRate
    
    ! Check for empty halos.
    if (Tree_Node_Hot_Halo_Mass(thisNode%formationNode) <= 0.0d0) then
       Cooling_Rate_Cole2000=0.0d0
       return
    end if

    ! Get the outer radius of the hot halo.
    outerRadius=Tree_Node_Hot_Halo_Outer_Radius(thisNode%formationNode)

    ! Get the infall radius.
    infallRadius=Infall_Radius                 (thisNode%formationNode)

    if (infallRadius >= outerRadius) then
       ! Cooling radius exceeds the outer radius - zero infall rate.
       Cooling_Rate_Cole2000=0.0d0
    else
       ! Find the density at the infall radius.
       infallDensity=Hot_Halo_Density                  (thisNode%formationNode,infallRadius)
       ! Find infall radius growth rate.
       infallRadiusGrowthRate=Infall_Radius_Growth_Rate(thisNode%formationNode             )
       ! Compute the infall rate.
       Cooling_Rate_Cole2000=4.0d0*Pi*(infallRadius**2)*infallDensity*infallRadiusGrowthRate
    end if
    return
  end function Cooling_Rate_Cole2000

end module Cooling_Rates_Cole2000

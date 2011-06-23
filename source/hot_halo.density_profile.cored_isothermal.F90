!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements a cored isothermal profile for hot gas halos.

module Hot_Halo_Density_Profile_Cored_Isothermal
  !% Implements a cored isothermal profile for hot gas halos.
  private
  public :: Hot_Halo_Density_Cored_Isothermal

contains

  !# <hotHaloDensityMethod>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal</unitName>
  !# </hotHaloDensityMethod>
  subroutine Hot_Halo_Density_Cored_Isothermal(hotHaloDensityMethod,Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get&
       &,Hot_Halo_Enclosed_Mass_Get,Hot_Halo_Profile_Rotation_Normalization_Get)
    !% Initialize the cored isothermal hot halo density profile module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: hotHaloDensityMethod
    procedure(double precision), pointer, intent(inout) :: Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get&
         &,Hot_Halo_Enclosed_Mass_Get,Hot_Halo_Profile_Rotation_Normalization_Get
    
    if (hotHaloDensityMethod == 'cored isothermal') then
       Hot_Halo_Density_Get                        => Hot_Halo_Density_Cored_Isothermal_Get
       Hot_Halo_Density_Log_Slope_Get              => Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get
       Hot_Halo_Enclosed_Mass_Get                  => Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get
       Hot_Halo_Profile_Rotation_Normalization_Get => Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get
    end if
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal

  double precision function Density_Normalization_Factor(coreRadius,virialRadius)
    !% Computes the density profile normalization factor for a given core radius and virial radius.
    implicit none
    double precision, intent(in) :: coreRadius,virialRadius
    
    Density_Normalization_Factor=1.0d0/(virialRadius/coreRadius-datan(virialRadius/coreRadius))
    return
  end function Density_Normalization_Factor

  double precision function Hot_Halo_Density_Cored_Isothermal_Get(thisNode,radius)
    !% Compute the density at radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Tree_Nodes
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: hotGasMass,virialRadius,coreRadius,densityNormalization

    hotGasMass          =Tree_Node_Hot_Halo_Mass(thisNode)
    virialRadius        =Dark_Matter_Halo_Virial_Radius(thisNode)
    coreRadius          =Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)
    densityNormalization=Density_Normalization_Factor(coreRadius,virialRadius)*hotGasMass/4.0d0/Pi/(coreRadius**3)
    Hot_Halo_Density_Cored_Isothermal_Get=densityNormalization/(1.0d0+(radius/coreRadius)**2)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Get
  
  double precision function Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get(thisNode,radius)
    !% Compute the density at radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: coreRadius,radiusInCoreUnitsSquared

    coreRadius              =Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)
    radiusInCoreUnitsSquared=(radius/coreRadius)**2
    Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get=-2.0d0*radiusInCoreUnitsSquared/(1.0d0+radiusInCoreUnitsSquared)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get
  
  double precision function Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get(thisNode,radius)
    !% Compute the mass enclosed within radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: hotGasMass,virialRadius,coreRadius

    hotGasMass  =Tree_Node_Hot_Halo_Mass                      (thisNode)
    virialRadius=Dark_Matter_Halo_Virial_Radius               (thisNode)
    ! Truncate the profile at the virial radius.
    if (radius >= virialRadius) then
       Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get=hotGasMass
       return
    end if
    coreRadius  =Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)
    Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get=hotGasMass*Density_Normalization_Factor(coreRadius,virialRadius)*(radius&
         &/coreRadius-datan(radius/coreRadius))

    return
  end function Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get

  double precision function Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get(thisNode)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    use Hot_Halo_Density_Cored_Isothermal_Core_Radii
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: radiusVirial,radiusCoreOverRadiusVirial

    ! Get virial radius and ratio of core radius to virial radius.
    radiusVirial              =Dark_Matter_Halo_Virial_Radius(thisNode)
    radiusCoreOverRadiusVirial=Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)/radiusVirial

    ! Compute the normalization.
    Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get=                                                   &
         &                                                       (                                                  &
         &                                                        1.0d0                                             &
         &                                                       - radiusCoreOverRadiusVirial                       &
         &                                                        *datan(1.0d0/radiusCoreOverRadiusVirial)          &
         &                                                      )                                                   &
         &                                                     /(                                                   &
         &                                                        0.5d0                                             &
         &                                                       + radiusCoreOverRadiusVirial**2                    &
         &                                                        *dlog(                                            &
         &                                                               radiusCoreOverRadiusVirial                 &
         &                                                              /dsqrt(1.0d0+radiusCoreOverRadiusVirial**2) &
         &                                                             )                                            &
         &                                                      )                                                   &
         &                                                     /radiusVirial
    return
  end function Hot_Halo_Profile_Rotation_Normalization_Cored_Isothermal_Get
  
end module Hot_Halo_Density_Profile_Cored_Isothermal

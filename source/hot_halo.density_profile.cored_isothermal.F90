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

  ! Parameter which specifies the ratio of the cored isothermal profile core radius to the virial radius.
  double precision :: isothermalCoreRadiusOverVirialRadius

  ! Pre-computed factor that appears in the density normalization.
  double precision :: denistyNormalizationFactor

contains

  !# <hotHaloDensityMethod>
  !#  <unitName>Hot_Halo_Density_Cored_Isothermal</unitName>
  !# </hotHaloDensityMethod>
  subroutine Hot_Halo_Density_Cored_Isothermal(hotHaloDensityMethod,Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get&
       &,Hot_Halo_Enclosed_Mass_Get)
    !% Initialize the cored isothermal hot halo density profile module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: hotHaloDensityMethod
    procedure(double precision), pointer, intent(inout) :: Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get&
         &,Hot_Halo_Enclosed_Mass_Get
    
    if (hotHaloDensityMethod == 'cored isothermal') then
       Hot_Halo_Density_Get           => Hot_Halo_Density_Cored_Isothermal_Get
       Hot_Halo_Density_Log_Slope_Get => Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get
       Hot_Halo_Enclosed_Mass_Get     => Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get
       !@ <inputParameter>
       !@   <name>isothermalCoreRadiusOverVirialRadius</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The core radius in the ``cored isothermal'' hot halo density profile in units of the virial radius.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('isothermalCoreRadiusOverVirialRadius',isothermalCoreRadiusOverVirialRadius,defaultValue=0.1d0)
       denistyNormalizationFactor=1.0d0/(1.0d0/isothermalCoreRadiusOverVirialRadius-datan(1.0d0&
            &/isothermalCoreRadiusOverVirialRadius))
    end if
    return
  end subroutine Hot_Halo_Density_Cored_Isothermal

  double precision function Hot_Halo_Density_Cored_Isothermal_Get(thisNode,radius)
    !% Compute the density at radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: hotGasMass,virialRadius,coreRadius,densityNormalization

    hotGasMass=Tree_Node_Hot_Halo_Mass(thisNode)
    virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
    coreRadius=isothermalCoreRadiusOverVirialRadius*virialRadius
    densityNormalization=denistyNormalizationFactor*hotGasMass/4.0d0/Pi/(coreRadius**3)
    Hot_Halo_Density_Cored_Isothermal_Get=densityNormalization/(1.0d0+(radius/coreRadius)**2)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Get
  
  double precision function Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get(thisNode,radius)
    !% Compute the density at radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: virialRadius,coreRadius,radiusInCoreUnitsSquared

    virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
    coreRadius=isothermalCoreRadiusOverVirialRadius*virialRadius
    radiusInCoreUnitsSquared=(radius/coreRadius)**2
    Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get=-2.0d0*radiusInCoreUnitsSquared/(1.0d0+radiusInCoreUnitsSquared)
    return
  end function Hot_Halo_Density_Cored_Isothermal_Log_Slope_Get
  
  double precision function Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get(thisNode,radius)
    !% Compute the mass enclosed within radius {\tt radius} in a cored isothermal hot halo density profile for {\tt thisNode}.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Math
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: hotGasMass,virialRadius,coreRadius

    hotGasMass=Tree_Node_Hot_Halo_Mass(thisNode)
    virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
    coreRadius=isothermalCoreRadiusOverVirialRadius*virialRadius
    Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get=hotGasMass*denistyNormalizationFactor*(radius/coreRadius-datan(radius&
         &/coreRadius))
    return
  end function Hot_Halo_Density_Cored_Isothermal_Enclosed_Mass_Get
  
end module Hot_Halo_Density_Profile_Cored_Isothermal

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


!% Contains a module which implements a warm dark matter scaling of critical overdensities for collapse based on the work of \cite{barkana_constraints_2001}.

module Critical_Overdensity_Mass_Scalings_WDM
  !% Implements a warm dark matter scaling of critical overdensities for collapse based on the work of \cite{barkana_constraints_2001}.
  use Tree_Nodes
  implicit none
  private
  public :: Critical_Overdensity_Mass_Scaling_WDM_Initialize
  
  ! The Jeans mass used in the warm dark matter fitting formula.
  double precision :: jeansMass

contains

  !# <criticalOverdensityMassScalingMethod>
  !#  <unitName>Critical_Overdensity_Mass_Scaling_WDM_Initialize</unitName>
  !# </criticalOverdensityMassScalingMethod>
  subroutine Critical_Overdensity_Mass_Scaling_WDM_Initialize(criticalOverdensityMassScalingMethod,Critical_Overdensity_Mass_Scaling_Get)
    !% Initializes the ``warmDarkMatter'' critical overdensity mass scaling method.
    use ISO_Varying_String
    use Input_Parameters
    use Cosmological_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: criticalOverdensityMassScalingMethod
    procedure(double precision), pointer, intent(inout) :: Critical_Overdensity_Mass_Scaling_Get
    double precision                                    :: matterRadiationEqualityRedshift,warmDarkMatterCriticalOverdensityGX&
         &,warmDarkMatterCriticalOverdensityMX    

    if (criticalOverdensityMassScalingMethod == 'warmDarkMatter') then
       ! Return a pointer to our implementation of the mass scaling function.
       Critical_Overdensity_Mass_Scaling_Get => Critical_Overdensity_Mass_Scaling_WDM

       ! Get warm dark matter particle properties.
       !@ <inputParameter>
       !@   <name>warmDarkMatterCriticalOverdensityGX</name>
       !@   <defaultValue>1.5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The effective number of degrees of freedom for the warm dark matter particle.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('warmDarkMatterCriticalOverdensityGX',warmDarkMatterCriticalOverdensityGX,defaultValue=1.5d0)
       !@ <inputParameter>
       !@   <name>warmDarkMatterCriticalOverdensityMX</name>
       !@   <defaultValue>1.0 keV</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass (in keV) of the warm dark matter particle.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('warmDarkMatterCriticalOverdensityMX',warmDarkMatterCriticalOverdensityMX,defaultValue=1.0d0)

       ! Compute corresponding Jeans mass.
       matterRadiationEqualityRedshift=3600.0d0*(Omega_Matter()*Little_H_0()**2/0.15d0)-1.0d0
       jeansMass=3.06d8*((1.0d0+matterRadiationEqualityRedshift)/3000.0d0)**1.5d0*dsqrt(Omega_Matter()*Little_H_0()**2/0.15d0)&
            &/(warmDarkMatterCriticalOverdensityGX/1.5d0)/(warmDarkMatterCriticalOverdensityMX/1.0d0)**4

    end if
    return
  end subroutine Critical_Overdensity_Mass_Scaling_WDM_Initialize

  double precision function Critical_Overdensity_Mass_Scaling_WDM(mass)
    !% Returns a mass scaling for critical overdensities based on the results of \cite{barkana_constraints_2001}. This method
    !% assumes that their results for the original collapse barrier (i.e. the critical overdensity, and which they call $B_0$)
    !% scale with the effective Jeans mass of the warm dark matter particle as computed using their eqn.~(10). A fitting function
    !% to $B_0$ with this scaling was found by fitting to the modfied barrier, $B$, as shown in their Fig.~2 (upper panel).
    implicit none
    double precision, intent(in) :: mass
    double precision, parameter  :: jeansMassFraction=0.125d0, exponent1=1.4d0, exponent2=0.45d0

    Critical_Overdensity_Mass_Scaling_WDM=dexp((jeansMassFraction*jeansMass/mass)**exponent1+(jeansMassFraction*jeansMass/mass)**exponent2)
    return
  end function Critical_Overdensity_Mass_Scaling_WDM
  
end module Critical_Overdensity_Mass_Scalings_WDM

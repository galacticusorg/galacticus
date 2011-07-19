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


!% Contains a module which implements calculations of halo bias using the fitting function of \cite{sheth_ellipsoidal_2001}.

module Dark_Matter_Halo_Biases_SMT
  !% Implements calculations of halo bias using the fitting function of \cite{sheth_ellipsoidal_2001}.
  private
  public :: Dark_Matter_Halo_Bias_SMT_Initialize

contains

  !# <darkMatterHaloBiasMethod>
  !#  <unitName>Dark_Matter_Halo_Bias_SMT_Initialize</unitName>
  !# </darkMatterHaloBiasMethod>
  subroutine Dark_Matter_Halo_Bias_SMT_Initialize(darkMatterHaloBiasMethod,Dark_Matter_Halo_Bias_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterHaloBiasMethod
    procedure(double precision), pointer, intent(inout) :: Dark_Matter_Halo_Bias_Get

    if (darkMatterHaloBiasMethod == 'SMT') Dark_Matter_Halo_Bias_Get => Dark_Matter_Halo_Bias_SMT
    return
  end subroutine Dark_Matter_Halo_Bias_SMT_Initialize

  double precision function Dark_Matter_Halo_Bias_SMT(thisNode)
    !% Computes the bias for a dark matter halo using the method of \cite{sheth_ellipsoidal_2001}.
    use Critical_Overdensity
    use CDM_Power_Spectrum
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, parameter              :: a=0.707d0, b=0.5d0, c=0.6d0
    double precision                         :: deltaCritical,sigma,nu

    ! Get critical overdensity for collapse and root-variance, then compute peak height parameter, nu.
    deltaCritical=Critical_Overdensity_for_Collapse(time=Tree_Node_Time(thisNode),mass=Tree_Node_Mass(thisNode))
    sigma        =sigma_CDM(Tree_Node_Mass(thisNode))
    nu           =deltaCritical/sigma
    
    ! Compute halo bias.
    Dark_Matter_Halo_Bias_SMT=1.0d0+(dsqrt(a)*(a*nu**2)+dsqrt(a)*b*(a*nu**2)**(1.0d0-c)-(a*nu**2)**c/((a*nu**2)**c+b*(1.0d0-c)&
         &*(1.0d0-c/2.0d0)))/dsqrt(a)/deltaCritical
    return
  end function Dark_Matter_Halo_Bias_SMT

end module Dark_Matter_Halo_Biases_SMT

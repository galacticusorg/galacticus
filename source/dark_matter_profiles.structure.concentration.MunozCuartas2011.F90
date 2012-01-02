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


!% Contains a module which implements the \cite{munoz-cuartas_redshift_2011} NFW halo concentration algorithm.

module Dark_Matter_Profiles_Concentrations_MunozCuartas2011
  !% Implements the \cite{munoz-cuartas_redshift_2011} NFW halo concentration algorithm.
  use Tree_Nodes
  implicit none
  private
  public :: Dark_Matter_Concentrations_MunozCuartas2011_Initialize

contains

  !# <darkMatterConcentrationMethod>
  !#  <unitName>Dark_Matter_Concentrations_MunozCuartas2011_Initialize</unitName>
  !# </darkMatterConcentrationMethod>
  subroutine Dark_Matter_Concentrations_MunozCuartas2011_Initialize(darkMatterConcentrationMethod,Dark_Matter_Profile_Concentration_Get)
    !% Initializes the ``Munoz-Cuartas2011'' halo concentration module.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterConcentrationMethod
    procedure(double precision), pointer, intent(inout) :: Dark_Matter_Profile_Concentration_Get
    
    if (darkMatterConcentrationMethod == 'Munoz-Cuartas2011') Dark_Matter_Profile_Concentration_Get => Dark_Matter_Profile_Concentration_MunozCuartas2011
  
    return
  end subroutine Dark_Matter_Concentrations_MunozCuartas2011_Initialize

  double precision function Dark_Matter_Profile_Concentration_MunozCuartas2011(thisNode)
    !% Returns the concentration of the dark matter profile of {\tt thisNode} using the method of \cite{munoz-cuartas_redshift_2011}.
    use Tree_Nodes
    use Cosmology_Functions
    use Cosmological_Parameters
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, parameter              :: w=0.029d0, m=0.097d0, alpha=-110.001d0, beta=2469.720d0, gamma=16.885d0
    double precision                         :: a,b,redshift,haloMassLogarithmic,concentrationLogarithmic

    redshift                =Redshift_from_Expansion_Factor(Expansion_Factor(Tree_Node_Time(thisNode)))
    a                       =w*redshift-m
    b                       =alpha/(redshift+gamma)+beta/(redshift+gamma)**2
    haloMassLogarithmic     =dlog10(Tree_Node_Mass(thisNode)*Little_H_0())
    concentrationLogarithmic=a*haloMassLogarithmic+b
    Dark_Matter_Profile_Concentration_MunozCuartas2011=10.0d0**concentrationLogarithmic
    return
  end function Dark_Matter_Profile_Concentration_MunozCuartas2011
  
end module Dark_Matter_Profiles_Concentrations_MunozCuartas2011

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


!% Contains a module which implements the \cite{prada_halo_2011} halo concentration algorithm.

module Dark_Matter_Profiles_Concentrations_Prada2011
  !% Implements the \cite{prada_halo_2011} halo concentration algorithm.
  use Tree_Nodes
  implicit none
  private
  public :: Dark_Matter_Concentrations_Prada2011_Initialize

  ! Parameters of the model.
  double precision :: prada2011ConcentrationA    ,prada2011ConcentrationB            ,prada2011ConcentrationC            , &
       &              prada2011ConcentrationD    ,prada2011ConcentrationC0           ,prada2011ConcentrationC1           , &
       &              prada2011ConcentrationAlpha,prada2011ConcentrationBeta         ,prada2011ConcentrationX0           , &
       &              prada2011ConcentrationX1   ,prada2011ConcentrationInverseSigma0,prada2011ConcentrationInverseSigma1

contains

  !# <darkMatterConcentrationMethod>
  !#  <unitName>Dark_Matter_Concentrations_Prada2011_Initialize</unitName>
  !# </darkMatterConcentrationMethod>
  subroutine Dark_Matter_Concentrations_Prada2011_Initialize(darkMatterConcentrationMethod,Dark_Matter_Profile_Concentration_Get)
    !% Initializes the ``Prada2011'' halo concentration module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterConcentrationMethod
    procedure(double precision), pointer, intent(inout) :: Dark_Matter_Profile_Concentration_Get
    
    if (darkMatterConcentrationMethod == 'Prada2011') then
       ! Return a pointer to our implementation.
       Dark_Matter_Profile_Concentration_Get => Dark_Matter_Profile_Concentration_Prada2011
       ! Get parameters of the model.
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationA</name>
       !@   <defaultValue>2.881 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $A$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationA",prada2011ConcentrationA,defaultValue=2.881d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationB</name>
       !@   <defaultValue>1.257 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $b$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationB",prada2011ConcentrationB,defaultValue=1.257d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationC</name>
       !@   <defaultValue>1.022 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $c$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationC",prada2011ConcentrationC,defaultValue=1.022d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationD</name>
       !@   <defaultValue>0.060 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $d$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationD",prada2011ConcentrationD,defaultValue=0.060d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationC0</name>
       !@   <defaultValue>3.681 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $c_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationC0",prada2011ConcentrationC0,defaultValue=3.681d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationC1</name>
       !@   <defaultValue>5.033 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $c_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationC1",prada2011ConcentrationC1,defaultValue=5.033d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationX0</name>
       !@   <defaultValue>0.424 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $x_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationX0",prada2011ConcentrationX0,defaultValue=0.424d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationX1</name>
       !@   <defaultValue>0.526 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $x_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationX1",prada2011ConcentrationX1,defaultValue=0.526d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationInverseSigma0</name>
       !@   <defaultValue>1.047 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\sigma^{-1}_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationInverseSigma0",prada2011ConcentrationInverseSigma0,defaultValue=1.047d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationInverseSigma1</name>
       !@   <defaultValue>1.646 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\sigma^{-1}_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationInverseSigma1",prada2011ConcentrationInverseSigma1,defaultValue=1.646d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationAlpha</name>
       !@   <defaultValue>6.948 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\alpha$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationAlpha",prada2011ConcentrationAlpha,defaultValue=6.948d0)
       !@ <inputParameter>
       !@   <name>prada2011ConcentrationBeta</name>
       !@   <defaultValue>7.386 \cite{prada_halo_2011}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\beta$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("prada2011ConcentrationBeta",prada2011ConcentrationBeta,defaultValue=7.386d0)
    end if
    return
  end subroutine Dark_Matter_Concentrations_Prada2011_Initialize

  double precision function Dark_Matter_Profile_Concentration_Prada2011(thisNode)
    !% Returns the concentration of the dark matter profile of {\tt thisNode} using the method of \cite{prada_halo_2011}.
    use Tree_Nodes
    use CDM_Power_Spectrum
    use Cosmology_Functions
    use Cosmological_Parameters
    use Linear_Growth
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: massNode,timeNode,x,sigmaPrime

    massNode  =Tree_Node_Mass(thisNode)
    timeNode  =Tree_Node_Time(thisNode)
    x         =(Omega_DE()/Omega_Matter())**(1.0d0/3.0d0)*Expansion_Factor(timeNode)
    sigmaPrime=B1(x)*sigma_CDM(massNode)*Linear_Growth_Factor(timeNode)
    Dark_Matter_Profile_Concentration_Prada2011=B0(x)*C(sigmaPrime)
    return
  end function Dark_Matter_Profile_Concentration_Prada2011
  
  double precision function B0(x)
    !% The function $B_0(x)$ as defined in eqn.~(18) of \cite{prada_halo_2011}.
    implicit none
    double precision, intent(in) :: x
    B0=           cMin(x)/           cMin(1.393d0)
  end function B0
  
  double precision function B1(x)
    !% The function $B_1(x)$ as defined in eqn.~(18) of \cite{prada_halo_2011}.
    implicit none
    double precision, intent(in) :: x
    B1=inverseSigmaMin(x)/inverseSigmaMin(1.393d0)
  end function B1
  
  double precision function cMin           (x)
    !% The function $c_{\rm min}(x)$ as defined in eqn.~(19) of \cite{prada_halo_2011}.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in) :: x
    cMin           =           prada2011ConcentrationC0+(           prada2011ConcentrationC1-           prada2011ConcentrationC0)*(atan(prada2011ConcentrationAlpha*(x-prada2011ConcentrationX0))/Pi+0.5d0)
  end function cMin
  
  double precision function inverseSigmaMin(x)
    !% The function $\sigma^{-1}_{\rm min}(x)$ as defined in eqn.~(20) of \cite{prada_halo_2011}.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in) :: x
    inverseSigmaMin=prada2011ConcentrationInverseSigma0+(prada2011ConcentrationInverseSigma1-prada2011ConcentrationInverseSigma0)*(atan(prada2011ConcentrationBeta *(x-prada2011ConcentrationX1))/Pi+0.5d0)
  end function inverseSigmaMin

  double precision function C(sigmaPrime)
    !% The function $\mathcal{C}(\sigma^\prime)$ as defined in eqn.~(17) of \cite{prada_halo_2011}.
    implicit none
    double precision, intent(in) :: sigmaPrime
    C=prada2011ConcentrationA*((sigmaPrime/prada2011ConcentrationB)**prada2011ConcentrationC+1.0d0)*exp(prada2011ConcentrationD/sigmaPrime**2)
  end function C

end module Dark_Matter_Profiles_Concentrations_Prada2011

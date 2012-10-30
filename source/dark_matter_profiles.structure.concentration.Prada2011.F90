!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the \cite{prada_halo_2011} halo concentration algorithm.

module Dark_Matter_Profiles_Concentrations_Prada2011
  !% Implements the \cite{prada_halo_2011} halo concentration algorithm.
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
    use Galacticus_Nodes
    use CDM_Power_Spectrum
    use Cosmology_Functions
    use Cosmological_Parameters
    use Linear_Growth
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic),                pointer :: thisBasicComponent
    double precision                                  :: massNode,timeNode,x,sigmaPrime

    thisBasicComponent => thisNode%basic()
    massNode  =thisBasicComponent%mass()
    timeNode  =thisBasicComponent%time()
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

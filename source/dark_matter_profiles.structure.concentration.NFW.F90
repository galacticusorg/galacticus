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


!% Contains a module which implements the \cite{navarro_structure_1996} NFW halo concentration algorithm.

module Dark_Matter_Profiles_Concentrations_NFW1996
  !% Implements the \cite{navarro_structure_1996} NFW halo concentration algorithm.
  use Tree_Nodes
  implicit none
  private
  public :: Dark_Matter_Concentrations_NFW1996_Initialize

  ! Module global variable used in root finding.
  double precision :: targetValue
  !$omp threadprivate(targetValue)

  ! Parameters of the fit.
  double precision :: nfw96ConcentrationC,nfw96ConcentrationF

contains

  !# <darkMatterConcentrationMethod>
  !#  <unitName>Dark_Matter_Concentrations_NFW1996_Initialize</unitName>
  !# </darkMatterConcentrationMethod>
  subroutine Dark_Matter_Concentrations_NFW1996_Initialize(darkMatterConcentrationMethod,Dark_Matter_Profile_Concentration_Get)
    !% Initializes the ``NFW 1996'' halo concentration module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterConcentrationMethod
    procedure(double precision), pointer, intent(inout) :: Dark_Matter_Profile_Concentration_Get
    
    if (darkMatterConcentrationMethod == 'NFW 1996') then
       ! Return a pointer to our implementation.
       Dark_Matter_Profile_Concentration_Get => Dark_Matter_Profile_Concentration_NFW1996
       ! Get parameters of the model.
       !@ <inputParameter>
       !@   <name>nfw96ConcentrationF</name>
       !@   <defaultValue>0.01 \cite{navarro_structure_1996}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $C$ appearing in the halo concentration algorithm of \cite{navarro_structure_1996}.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("nfw96ConcentrationF",nfw96ConcentrationF,defaultValue=   0.01d0)
       !@ <inputParameter>
       !@   <name>nfw96ConcentrationC</name>
       !@   <defaultValue>2000 \cite{navarro_structure_1996}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $f$ appearing in the halo concentration algorithm of \cite{navarro_structure_1996}.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("nfw96ConcentrationC",nfw96ConcentrationC,defaultValue=2000.0d0)
    end if

    return
  end subroutine Dark_Matter_Concentrations_NFW1996_Initialize

  double precision function Dark_Matter_Profile_Concentration_NFW1996(thisNode)
    !% Returns the concentration of the dark matter profile of {\tt thisNode} using the method of \cite{navarro_structure_1996}.
    use, intrinsic :: ISO_C_Binding
    use Tree_Nodes
    use CDM_Power_Spectrum
    use Cosmology_Functions
    use Critical_Overdensity
    use Root_Finder
    use FGSL
    use Virial_Density_Contrast
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    double precision,        parameter              :: fitParameterNuHalf=0.47693628d0
    double precision,        parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6
    type(fgsl_function),     save                   :: rootFunction
    type(fgsl_root_fsolver), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision                                :: nodeMass,nodeTime,collapseMass,collapseCriticalOverdensity,collapseTime&
         &,collapseExpansionFactor,expansionFactor,collapseOverdensity,concentrationMinimum,concentrationMaximum
    type(c_ptr)                                     :: parameterPointer

    ! Get the properties of the node.
    nodeMass                   =Tree_Node_Mass(thisNode)
    nodeTime                   =Tree_Node_Time(thisNode)
    expansionFactor            =Expansion_Factor(nodeTime)

    ! Compute the mass of a progenitor as defined by NFW.
    collapseMass               =nfw96ConcentrationF*nodeMass

    ! Find the time of collapse for this progenitor.
    collapseCriticalOverdensity=dsqrt(2.0d0*fitParameterNuHalf**2*(sigma_CDM(collapseMass)**2-sigma_CDM(nodeMass)**2))+Critical_Overdensity_for_Collapse(nodeTime)
    collapseTime               =Time_of_Collapse(collapseCriticalOverdensity)
    collapseExpansionFactor    =Expansion_Factor(collapseTime)

    ! Compute the overdensity of the progenitor at collapse using the scaling given by NFW.
    collapseOverdensity        =nfw96ConcentrationC*(expansionFactor/collapseExpansionFactor)**3

    ! Find the ratio of this overdensity to that at for the present node.
    targetValue                =collapseOverdensity/Halo_Virial_Density_Contrast(nodeTime)

    ! Solve for the corresponding concentration.
    concentrationMinimum= 1.0d0
    concentrationMaximum=20.0d0
    Dark_Matter_Profile_Concentration_NFW1996=Root_Find(concentrationMinimum,concentrationMaximum,NFW_Concentration_Function_Root,parameterPointer &
            &,rootFunction,rootFunctionSolver,toleranceAbsolute,toleranceRelative)

    return
  end function Dark_Matter_Profile_Concentration_NFW1996
  
  function NFW_Concentration_Function_Root(concentration,parameterPointer) bind(c)
    !% Root function used in finding concentrations in the \cite{navarro_structure_1996} method.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double), value   :: concentration
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: NFW_Concentration_Function_Root
    
    NFW_Concentration_Function_Root=concentration**3/(dlog(1.0d0+concentration)-concentration/(1.0d0+concentration))/3.0d0-targetValue
    return
  end function NFW_Concentration_Function_Root

end module Dark_Matter_Profiles_Concentrations_NFW1996

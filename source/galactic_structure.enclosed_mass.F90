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


!% Contains a module which implements calculations of the mass enclosed within a specified radius.

module Galactic_Structure_Enclosed_Masses
  !% Implements calculations of the mass enclosed within a specified radius.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  use Tree_Nodes
  use Galactic_Structure_Options
  private
  public :: Galactic_Structure_Enclosed_Mass, Galactic_Structure_Radius_Enclosing_Mass

  ! Variables used in root finding.
  integer                  :: massTypeRoot,componentTypeRoot,weightByRoot,weightIndexRoot
  double precision         :: massRoot
  type(treeNode),  pointer :: activeNode
  !$omp threadprivate(massRoot,massTypeRoot,componentTypeRoot,weightByRoot,weightIndexRoot,activeNode)

contains

  double precision function Galactic_Structure_Enclosed_Mass(thisNode,radius,massType,componentType,weightBy,weightIndex)
    !% Solve for the mass within a given radius, or the total mass if no radius is specified. Assumes that galactic structure has
    !% already been computed.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="enclosedMassTask" type="moduleUse">
    include 'galactic_structure.enclosed_mass.tasks.modules.inc'
    !# </include>
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: massType,componentType,weightBy,weightIndex
    double precision, intent(in),    optional :: radius
    integer                                   :: massTypeActual,componentTypeActual,weightByActual,weightIndexActual
    double precision                          :: radiusActual,componentMass

    ! Determine which radius to use.
    if (present(radius)) then
       radiusActual=radius
    else
       radiusActual=radiusLarge
    end if

    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeActual=massType
    else
       massTypeActual=massTypeAll
    end if

    ! Determine which weighting to use.
    if (present(weightBy)) then
       weightByActual=weightBy
       select case (weightByActual)
       case (weightByLuminosity)
          if (.not.present(weightIndex)) call Galacticus_Error_Report('Galactic_Structure_Enclosed_Mass','weightIndex should be specified for luminosity weighting')
          weightIndexActual=weightIndex
       end select
    else
       weightByActual=weightByMass
    end if

    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeActual=componentType
    else
       componentTypeActual=componentTypeAll
    end if

    ! Initialize to zero mass.
    Galactic_Structure_Enclosed_Mass=0.0d0

    ! Call routines to supply the masses for all components.
    !# <include directive="enclosedMassTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,radiusActual,massTypeActual,componentTypeActual,weightByActual,weightIndexActual,componentMass</subroutineArgs>
    !#  <subroutineAction>Galactic_Structure_Enclosed_Mass=Galactic_Structure_Enclosed_Mass+componentMass</subroutineAction>
    include 'galactic_structure.enclosed_mass.tasks.inc'
    !# </include>

    return
  end function Galactic_Structure_Enclosed_Mass
  
  double precision function Galactic_Structure_Radius_Enclosing_Mass(thisNode,mass,fractionalMass,massType,componentType,weightBy,weightIndex)
    !% Return the radius enclosing a given mass (or fractional mass) in {\tt thisNode}.
    use Galacticus_Error
    use Root_Finder
    use FGSL
    use Dark_Matter_Halo_Scales
    use Galacticus_Display
    use ISO_Varying_String
    use String_Handling
    implicit none
    type(treeNode),          intent(inout), pointer  :: thisNode
    integer,                 intent(in),    optional :: massType,componentType,weightBy,weightIndex
    double precision,        intent(in),    optional :: mass,fractionalMass
    type(fgsl_function),     save                    :: rootFunction
    type(fgsl_root_fsolver), save                    :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision                                 :: radiusMinimum,radiusMaximum
    type(c_ptr)                                      :: parameterPointer
    type(varying_string)                             :: message
    character(len=11)                                :: massLabel

    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeRoot=massType
    else
       massTypeRoot=massTypeAll
    end if

    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeRoot=componentType
    else
       componentTypeRoot=componentTypeAll
    end if

    ! Determine which weighting to use.
    if (present(weightBy)) then
       weightByRoot=weightBy
       select case (weightByRoot)
       case (weightByLuminosity)
          if (.not.present(weightIndex)) call Galacticus_Error_Report('Galactic_Structure_Radius_Enclosing_Mass','weightIndex should be specified for luminosity weighting')
          weightIndexRoot=weightIndex
       end select
    else
       weightByRoot=weightByMass
    end if

    ! Determine what mass to use.
    if (present(mass)) then
       if (present(fractionalMass)) call Galacticus_Error_Report('Galactic_Structure_Radius_Enclosing_Mass','only one mass or&
            & fractionalMass can be specified')
       massRoot=mass
    else if (present(fractionalMass)) then
       if (fractionalMass >= 1.0d0) then
          Galactic_Structure_Radius_Enclosing_Mass=radiusLarge
          return
       end if
       massRoot=fractionalMass*Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeRoot,weightBy=weightByRoot,weightIndex=weightIndexRoot)
    else
       call Galacticus_Error_Report('Galactic_Structure_Radius_Enclosing_Mass','either mass or fractionalMass must be specified')
    end if
    if (massRoot <= 0.0d0) then
       Galactic_Structure_Radius_Enclosing_Mass=0.0d0
       return
    end if

    ! Solve for the radius.
    activeNode => thisNode
    radiusMinimum=0.0d0
    if (Enclosed_Mass_Root(radiusMinimum,parameterPointer) >= 0.0d0) then
       message='Enclosed mass in galaxy (ID='
       write (massLabel,'(e10.4)') Galactic_Structure_Enclosed_Mass(activeNode,radiusMinimum,massTypeRoot,componentTypeRoot)
       message=message//thisNode%index()//') seems to be finite ('//trim(massLabel)
       write (massLabel,'(e10.4)') massRoot
       message=message//') at zero radius (was seeking '//trim(massLabel)
       message=message//') - expect a crash.'
       call Galacticus_Display_Message(message,verbosityInfo)
    end if
    radiusMaximum=Dark_Matter_Halo_Virial_Radius(thisNode)
    do while (Enclosed_Mass_Root(radiusMaximum,parameterPointer) <= 0.0d0)
       radiusMaximum=radiusMaximum*2.0d0
    end do
    Galactic_Structure_Radius_Enclosing_Mass=Root_Find(radiusMinimum,radiusMaximum,Enclosed_Mass_Root,parameterPointer&
         &,rootFunction,rootFunctionSolver,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    
    return
  end function Galactic_Structure_Radius_Enclosing_Mass

  function Enclosed_Mass_Root(radius,parameterPointer) bind(c)
    !% Root function used in solving for the radius that encloses a given mass.
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: Enclosed_Mass_Root
   
    ! Evaluate the root function.
    Enclosed_Mass_Root=Galactic_Structure_Enclosed_Mass(activeNode,radius,massTypeRoot,componentTypeRoot,weightByRoot,weightIndexRoot)-massRoot
  end function Enclosed_Mass_Root

end module Galactic_Structure_Enclosed_Masses

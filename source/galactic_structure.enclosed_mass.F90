!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  integer                  :: massTypeRoot,componentTypeRoot
  double precision         :: massRoot
  type(treeNode),  pointer :: activeNode
  !$omp threadprivate(massRoot,massTypeRoot,componentTypeRoot,activeNode)

contains

  double precision function Galactic_Structure_Enclosed_Mass(thisNode,radius,massType,componentType)
    !% Solve for the mass within a given radius, or the total mass if no radius is specified. Assumes that galactic structure has
    !% already been computed.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="enclosedMassTask" type="moduleUse">
    include 'galactic_structure.enclosed_mass.tasks.modules.inc'
    !# </include>
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: massType,componentType
    double precision, intent(in),    optional :: radius
    integer                                   :: massTypeActual,componentTypeActual
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
    !#  <subroutineArgs>thisNode,radiusActual,massTypeActual,componentTypeActual,componentMass</subroutineArgs>
    !#  <subroutineAction>Galactic_Structure_Enclosed_Mass=Galactic_Structure_Enclosed_Mass+componentMass</subroutineAction>
    include 'galactic_structure.enclosed_mass.tasks.inc'
    !# </include>

    return
  end function Galactic_Structure_Enclosed_Mass
  
  double precision function Galactic_Structure_Radius_Enclosing_Mass(thisNode,mass,fractionalMass,massType,componentType)
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
    integer,                 intent(in),    optional :: massType,componentType
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
       massRoot=fractionalMass*Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeRoot)
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
       call Galacticus_Display_Message(message,1)
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
    Enclosed_Mass_Root=Galactic_Structure_Enclosed_Mass(activeNode,radius,massTypeRoot,componentTypeRoot)-massRoot
  end function Enclosed_Mass_Root

end module Galactic_Structure_Enclosed_Masses

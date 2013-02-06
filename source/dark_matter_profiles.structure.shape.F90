!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of dark matter halo density profile shapes.

module Dark_Matter_Profiles_Shapes
  !% Implements calculations of dark matter halo density profile shapes.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Dark_Matter_Profile_Shape

  ! Flag to indicate if this module has been initialized.  
  logical              :: darkMatterShapeInitialized=.false.

  ! Name of dark matter profile shape method used.
  type(varying_string) :: darkMatterShapeMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Dark_Matter_Profile_Template), pointer :: Dark_Matter_Profile_Shape_Get => null()
  abstract interface
     double precision function Dark_Matter_Profile_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Dark_Matter_Profile_Template
  end interface

contains

  subroutine Dark_Matter_Shapes_Initialize
    !% Initialize the dark matter profile module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="darkMatterShapeMethod" type="moduleUse">
    include 'dark_matter_profiles.structure.shape.modules.inc'
    !# </include>
    implicit none
    
    ! Initialize if necessary.
    if (.not.darkMatterShapeInitialized) then
       !$omp critical(Dark_Matter_Shapes_Initialization) 
       if (.not.darkMatterShapeInitialized) then
          ! Get the halo spin distribution method parameter.
          !@ <inputParameter>
          !@   <name>darkMatterShapeMethod</name>
          !@   <defaultValue>Gao2008</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of dark matter halo density profile shapes.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('darkMatterShapeMethod',darkMatterShapeMethod,defaultValue='Gao2008')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="darkMatterShapeMethod" type="functionCall" functionType="void">
          !#  <functionArgs>darkMatterShapeMethod,Dark_Matter_Profile_Shape_Get</functionArgs>
          include 'dark_matter_profiles.structure.shape.inc'
          !# </include>
          if (.not.associated(Dark_Matter_Profile_Shape_Get)) &
               & call Galacticus_Error_Report('Dark_Matter_Shapes_Initialize','method ' //char(darkMatterShapeMethod)//' is unrecognized')
          darkMatterShapeInitialized=.true.
       end if
       !$omp end critical(Dark_Matter_Shapes_Initialization) 
    end if
    return
  end subroutine Dark_Matter_Shapes_Initialize

  double precision function Dark_Matter_Profile_Shape(thisNode)
    !% Returns the shape of the dark matter profile of {\tt thisNode}.
    implicit none
    type(treeNode),pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Dark_Matter_Shapes_Initialize

    ! Get the shape using the selected method.
    Dark_Matter_Profile_Shape=Dark_Matter_Profile_Shape_Get(thisNode)

    return
  end function Dark_Matter_Profile_Shape

end module Dark_Matter_Profiles_Shapes

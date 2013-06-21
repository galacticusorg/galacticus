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

!% Contains a module which implements calculations of the initial radius in the dark matter halo for use when solving for galactic
!% structure.

module Galactic_Structure_Initial_Radii
  !% Implements calculations of the initial radius in the dark matter halo for use when solving for galactic structure.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Galactic_Structure_Radius_Initial, Galactic_Structure_Radius_Initial_Derivative

  ! Flag to indicate if this module has been initialized.  
  logical                                                          :: moduleInitialized                               =.false.  
  
  ! Name of cooling rate available method used.                                                                                                                           
  type     (varying_string                              )          :: galacticStructureRadiusSolverInitialRadiusMethod          
  
  ! Pointer to the function that actually does the calculation.                                                                                                                           
  procedure(Galactic_Structure_Radius_Initial           ), pointer :: Galactic_Structure_Radius_Initial_Get           =>null()  
  procedure(Galactic_Structure_Radius_Initial_Derivative), pointer :: Galactic_Structure_Radius_Initial_Derivative_Get=>null()  
                                                                                                                             
contains

  subroutine Galactic_Structure_Radius_Initial_Initialize()
    !% Initialize the initial radius module for the galacticu structure subsystem.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="galacticStructureRadiusSolverInitialRadiusMethod" type="moduleUse">
    include 'galactic_structure.radius_solver.initial_radii.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galactic_Structure_Radius_Initial_Initialization) 
       if (.not.moduleInitialized) then
          ! Get the galactic structure radii solver method parameter.
          !@ <inputParameter>
          !@   <name>galacticStructureRadiusSolverInitialRadiusMethod</name>
          !@   <defaultValue>adiabatic</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Selects the method to be used to determine initial radii in the dark matter halo when solving for galactic structure.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('galacticStructureRadiusSolverInitialRadiusMethod',galacticStructureRadiusSolverInitialRadiusMethod,defaultValue='adiabatic')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="galacticStructureRadiusSolverInitialRadiusMethod" type="functionCall" functionType="void">
          !#  <functionArgs>galacticStructureRadiusSolverInitialRadiusMethod,Galactic_Structure_Radius_Initial_Get,Galactic_Structure_Radius_Initial_Derivative_Get</functionArgs>
          include 'galactic_structure.radius_solver.initial_radii.inc'
          !# </include>
          if     (                                                                                      &
               &  .not.                                                                                 &
               &       (                                                                                &
               &         associated(Galactic_Structure_Radius_Initial_Get           )                   &
               &        .and.                                                                           &
               &         associated(Galactic_Structure_Radius_Initial_Derivative_Get)                   &
               &       )                                                                                &
               & )                                                                                      &
               & call Galacticus_Error_Report(                                                          &
               &                              'Galactic_Structure_Radius_Initial_Initialize'          , &
               &                              'method '                                                 &
               &                              //char(galacticStructureRadiusSolverInitialRadiusMethod)  &
               &                              //' is unrecognized'                                      &
               &                             )
          moduleInitialized=.true.
       end if
       !$omp end critical(Galactic_Structure_Radius_Initial_Initialization) 
    end if

    return
  end subroutine Galactic_Structure_Radius_Initial_Initialize

  double precision function Galactic_Structure_Radius_Initial(thisNode,radius)
    !% Find the initial radius in the dark matter halo of {\tt thisNode} corresponding to the given final {\tt radius}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: radius    
    
    ! Ensure that the module is initialized.                                                             
    call Galactic_Structure_Radius_Initial_Initialize()
    ! Get the initial radius.
    Galactic_Structure_Radius_Initial=Galactic_Structure_Radius_Initial_Get(thisNode,radius)
    return
  end function Galactic_Structure_Radius_Initial

  double precision function Galactic_Structure_Radius_Initial_Derivative(thisNode,radius)
    !% Find the derivative of the initial radius in the dark matter halo of {\tt thisNode} with respect to the final radius
    !% corresponding to the given final {\tt radius}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: radius    
    
    ! Ensure that the module is initialized.                                                             
    call Galactic_Structure_Radius_Initial_Initialize()
    ! Get the initial radius.
    Galactic_Structure_Radius_Initial_Derivative=Galactic_Structure_Radius_Initial_Derivative_Get(thisNode,radius)
    return
  end function Galactic_Structure_Radius_Initial_Derivative
  
end module Galactic_Structure_Initial_Radii

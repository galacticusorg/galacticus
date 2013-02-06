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

!% Contains a module which implements calculations of bar instability in galactic disks.

module Galactic_Dynamics_Bar_Instabilities
  !% Implements calculations of bar instability in galactic disks.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Bar_Instability_Timescale

  ! Flag to indicate if this module has been initialized.  
  logical              :: barInstabilitiesInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: barInstabilityMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Bar_Instability_Template), pointer :: Bar_Instability_Timescale_Get => null()
  abstract interface
     double precision function Bar_Instability_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Bar_Instability_Template
  end interface

contains

  subroutine Galactic_Dynamics_Bar_Instability_Initialize
    !% Initialize the bar instability module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="barInstabilityMethod" type="moduleUse">
    include 'galactic_dynamics.bar_instability.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.barInstabilitiesInitialized) then
       !$omp critical(Galactic_Dynamics_Bar_Instability_Initialize)
       if (.not.barInstabilitiesInitialized) then
          ! Get the halo spin distribution method parameter.
          !@ <inputParameter>
          !@   <name>barInstabilityMethod</name>
          !@   <defaultValue>ELN</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for bar instability calculations.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('barInstabilityMethod',barInstabilityMethod,defaultValue='ELN')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="barInstabilityMethod" type="functionCall" functionType="void">
          !#  <functionArgs>barInstabilityMethod,Bar_Instability_Timescale_Get</functionArgs>
          include 'galactic_dynamics.bar_instability.inc'
          !# </include>
          if (.not.associated(Bar_Instability_Timescale_Get)) &
               & call Galacticus_Error_Report('Galactic_Dynamics_Bar_Instability_Initialize','method ' //char(barInstabilityMethod)//' is unrecognized')
          barInstabilitiesInitialized=.true.
       end if
       !$omp end critical(Galactic_Dynamics_Bar_Instability_Initialize) 
    end if
    return
  end subroutine Galactic_Dynamics_Bar_Instability_Initialize

  double precision function Bar_Instability_Timescale(thisNode)
    !% Returns a timescale on which the bar instability depletes material from a disk into a pseudo-bulge. A negative value
    !% indicates no instability.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Galactic_Dynamics_Bar_Instability_Initialize

    ! Get the timescale using the selected method.
    Bar_Instability_Timescale=Bar_Instability_Timescale_Get(thisNode)

    return
  end function Bar_Instability_Timescale

end module Galactic_Dynamics_Bar_Instabilities

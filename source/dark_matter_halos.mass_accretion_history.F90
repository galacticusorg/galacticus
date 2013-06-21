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

!% Contains a module which implements calculations of dark matter halo mass accretion histories.

module Dark_Matter_Halo_Mass_Accretion_Histories
  !% Implements calculations of dark matter halo mass accretion histories.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Dark_Matter_Halo_Mass_Accretion_Time

  ! Flag to indicate if this module has been initialized.  
  logical                                            :: darkMatterAccretionHistoryInitialized   =.false.  
  
  ! Name of cooling rate available method used.                                                                                                     
  type     (varying_string                )          :: darkMatterAccretionHistoryMethod                  
  
  ! Pointer to the function that actually does the calculation.                                                                                                     
  procedure(Dark_Matter_Accretion_Template), pointer :: Dark_Matter_Halo_Mass_Accretion_Time_Get=>null()  
  abstract interface
     double precision function Dark_Matter_Accretion_Template(baseNode,nodeMass)
       import treeNode
       type            (treeNode), intent(inout), pointer :: baseNode  
       double precision          , intent(in   )          :: nodeMass  
     end function Dark_Matter_Accretion_Template
  end interface

contains

  subroutine Dark_Matter_Mass_Accretion_Initialize
    !% Initialize the dark matter mass accretion history module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="darkMatterAccretionHistoryMethod" type="moduleUse">
    include 'dark_matter_halos.mass_accretion_history.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.darkMatterAccretionHistoryInitialized) then
       !$omp critical(Dark_Matter_Mass_Accretion_Initialization) 
       if (.not.darkMatterAccretionHistoryInitialized) then
          ! Get the mass accretion history method parameter.
          !@ <inputParameter>
          !@   <name>darkMatterAccretionHistoryMethod</name>
          !@   <defaultValue>Wechsler2002</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of dark matter halo mass accretion histories.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('darkMatterAccretionHistoryMethod',darkMatterAccretionHistoryMethod,defaultValue='Wechsler2002')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="darkMatterAccretionHistoryMethod" type="functionCall" functionType="void">
          !#  <functionArgs>darkMatterAccretionHistoryMethod,Dark_Matter_Halo_Mass_Accretion_Time_Get</functionArgs>
          include 'dark_matter_halos.mass_accretion_history.inc'
          !# </include>
          if (.not.associated(Dark_Matter_Halo_Mass_Accretion_Time_Get)) &
               & call Galacticus_Error_Report('Dark_Matter_Mass_Accretion_Initialize','method ' //char(darkMatterAccretionHistoryMethod)//' is unrecognized')
          darkMatterAccretionHistoryInitialized=.true.
       end if
       !$omp end critical(Dark_Matter_Mass_Accretion_Initialization) 
    end if
    return
  end subroutine Dark_Matter_Mass_Accretion_Initialize

  double precision function Dark_Matter_Halo_Mass_Accretion_Time(baseNode,nodeMass)
    !% Returns the time for {\tt thisNode} in {\tt thisTree} according to the mass accretion history.
    implicit none
    type            (treeNode), intent(inout), pointer :: baseNode  
    double precision          , intent(in   )          :: nodeMass  
    
    ! Initialize the module.                                                             
    call Dark_Matter_Mass_Accretion_Initialize

    ! Get the time for the node.
    Dark_Matter_Halo_Mass_Accretion_Time=Dark_Matter_Halo_Mass_Accretion_Time_Get(baseNode,nodeMass)

    return
  end function Dark_Matter_Halo_Mass_Accretion_Time
  
end module Dark_Matter_Halo_Mass_Accretion_Histories

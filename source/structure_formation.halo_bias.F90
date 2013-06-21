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

!% Contains a module which implements calculations of dark matter halo bias.

module Dark_Matter_Halo_Biases
  !% Implements calculations of dark matter halo bias.
  use Galacticus_Nodes
  use ISO_Varying_String
  implicit none
  private
  public :: Dark_Matter_Halo_Bias

  ! Flag to indicate if this module has been initialized.
  logical                                     :: haloBiasInitialized           =.false.

  ! Name of halo bias method used.
  type     (varying_string         )          :: darkMatterHaloBiasMethod

  ! Pointer to the function that returns halo bias.
  procedure(Halo_Bias_Node_Template), pointer :: Dark_Matter_Halo_Bias_Node_Get=>null()
  procedure(Halo_Bias_Template     ), pointer :: Dark_Matter_Halo_Bias_Get     =>null()
  abstract interface
     double precision function Halo_Bias_Node_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Halo_Bias_Node_Template
  end interface
  abstract interface
     double precision function Halo_Bias_Template(mass,time)
       double precision, intent(in   ) :: mass, time
     end function Halo_Bias_Template
  end interface

  interface Dark_Matter_Halo_Bias
     module procedure Dark_Matter_Halo_Bias_By_Node
     module procedure Dark_Matter_Halo_Bias_By_Mass
  end interface Dark_Matter_Halo_Bias

contains

  subroutine Dark_Matter_Halo_Bias_Initialize
    !% Initalize the dark matter halo bias module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="darkMatterHaloBiasMethod" type="moduleUse">
    include 'structure_formation.halo_bias.modules.inc'
    !# </include>
    implicit none

    if (.not.haloBiasInitialized) then
       !$omp critical(haloBiasInitialize)
       if (.not.haloBiasInitialized) then
          ! Get the halo bias method parameter.
          !@ <inputParameter>
          !@   <name>darkMatterHaloBiasMethod</name>
          !@   <defaultValue>Tinker2010</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Selects which dark matter halo bias method to use.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('darkMatterHaloBiasMethod',darkMatterHaloBiasMethod,defaultValue='Tinker2010')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="darkMatterHaloBiasMethod" type="functionCall" functionType="void">
          !#  <functionArgs>darkMatterHaloBiasMethod,Dark_Matter_Halo_Bias_Node_Get,Dark_Matter_Halo_Bias_Get</functionArgs>
          include 'structure_formation.halo_bias.inc'
          !# </include>
          if (.not.(associated(Dark_Matter_Halo_Bias_Get).and.associated(Dark_Matter_Halo_Bias_Node_Get))) call&
               & Galacticus_Error_Report('Dark_Matter_Halo_Bias_Initialize','method '//char(darkMatterHaloBiasMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          haloBiasInitialized=.true.
       end if
       !$omp end critical(haloBiasInitialize)
    end if
    return
  end subroutine Dark_Matter_Halo_Bias_Initialize

  double precision function Dark_Matter_Halo_Bias_By_Node(thisNode)
    !% Computes the bias for a dark matter halo.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Ensure the module is initalized.
    call Dark_Matter_Halo_Bias_Initialize

    ! Get the dark matter halo bias.
    Dark_Matter_Halo_Bias_By_Node=Dark_Matter_Halo_Bias_Node_Get(thisNode)

    return
  end function Dark_Matter_Halo_Bias_By_Node

  double precision function Dark_Matter_Halo_Bias_By_Mass(mass,time)
    !% Computes the bias for a dark matter halo.
    implicit none
    double precision, intent(in   ) :: mass, time

    ! Ensure the module is initalized.
    call Dark_Matter_Halo_Bias_Initialize

    ! Get the dark matter halo bias.
    Dark_Matter_Halo_Bias_By_Mass=Dark_Matter_Halo_Bias_Get(mass,time)

    return
  end function Dark_Matter_Halo_Bias_By_Mass

end module Dark_Matter_Halo_Biases

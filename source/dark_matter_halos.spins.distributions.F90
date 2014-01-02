!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module that implements calculations of dark matter halo spin distributions

module Halo_Spin_Distributions
  !% Implements calculations of dark matter halo spin distributions
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Halo_Spin_Distribution_Sample

  ! Flag to indicate if this module has been initialized.
  logical                                           :: haloSpinDistributionInitialized=.false.

  ! Name of cooling rate available method used.
  type     (varying_string               )          :: haloSpinDistributionMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Halo_Spin_Sample_Get_Template), pointer :: Halo_Spin_Sample_Get           =>null()
  interface Halo_Spin_Sample_Get_Template
     double precision function Halo_Spin_Sample_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Halo_Spin_Sample_Get_Template
  end interface

contains

  double precision function Halo_Spin_Distribution_Sample(thisNode)
    !% Return a halo spin selected randomly from a distribution.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="haloSpinDistributionMethod" type="moduleUse">
    include 'dark_matter_halos.spins.distributions.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize if necessary.
    if (.not.haloSpinDistributionInitialized) then
       !$omp critical(Halo_Spin_Distribution_Initialization)
       if (.not.haloSpinDistributionInitialized) then
          ! Get the halo spin distribution method parameter.
          !@ <inputParameter>
          !@   <name>haloSpinDistributionMethod</name>
          !@   <defaultValue>Bett2007</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be use for computing halo spin distributions.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('haloSpinDistributionMethod',haloSpinDistributionMethod,defaultValue='Bett2007')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="haloSpinDistributionMethod" type="functionCall" functionType="void">
          !#  <functionArgs>haloSpinDistributionMethod,Halo_Spin_Sample_Get</functionArgs>
          include 'dark_matter_halos.spins.distributions.inc'
          !# </include>
          if (.not.associated(Halo_Spin_Sample_Get)) call Galacticus_Error_Report('Halo_Spin_Distribution','method ' &
               &//char(haloSpinDistributionMethod)//' is unrecognized')
          haloSpinDistributionInitialized=.true.
       end if
       !$omp end critical(Halo_Spin_Distribution_Initialization)
    end if

    ! Get the cooling rate using the selected method.
    Halo_Spin_Distribution_Sample=Halo_Spin_Sample_Get(thisNode)
    return
  end function Halo_Spin_Distribution_Sample

end module Halo_Spin_Distributions

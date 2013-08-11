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

!% Contains a module which implements a delta function halo spin distribution.

module Halo_Spin_Distributions_Delta_Function
  !% Implements a delta function halo spin distribution (i.e. all halos have the same spin).
  implicit none
  private
  public :: Halo_Spin_Distribution_Delta_Function_Initialize

  ! Parameters of the spin distribution.
  double precision :: deltaFunctionSpinDistributionSpin

contains

  !# <haloSpinDistributionMethod>
  !#  <unitName>Halo_Spin_Distribution_Delta_Function_Initialize</unitName>
  !# </haloSpinDistributionMethod>
  subroutine Halo_Spin_Distribution_Delta_Function_Initialize(haloSpinDistributionMethod,Halo_Spin_Sample_Get)
    !% Initializes the ``delta function'' halo spin distribution module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                       ), intent(in   )          :: haloSpinDistributionMethod
    procedure(Halo_Spin_Distribution_Delta_Function), intent(inout), pointer :: Halo_Spin_Sample_Get

    if (haloSpinDistributionMethod == 'deltaFunction') then
       Halo_Spin_Sample_Get => Halo_Spin_Distribution_Delta_Function
       !@ <inputParameter>
       !@   <name>deltaFunctionSpinDistributionSpin</name>
       !@   <defaultValue>0.03687 \citep{bett_spin_2007}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The fixed value of spin in a delta function spin distribution.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('deltaFunctionSpinDistributionSpin',deltaFunctionSpinDistributionSpin,defaultValue=0.03687d0)
    end if
    return
  end subroutine Halo_Spin_Distribution_Delta_Function_Initialize

  double precision function Halo_Spin_Distribution_Delta_Function(thisNode)
    !% Return a halo spin from a delta function distribution.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Halo_Spin_Distribution_Delta_Function=deltaFunctionSpinDistributionSpin
    return
  end function Halo_Spin_Distribution_Delta_Function

end module Halo_Spin_Distributions_Delta_Function

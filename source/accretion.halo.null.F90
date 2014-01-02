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

!% Contains a module which implements calculations of null baryonic accretion.

module Accretion_Halos_Null
  !% Implements calculations of null baryonic accretion.
  use Abundances_Structure
  implicit none
  private
  public :: Accretion_Halos_Null_Initialize

contains

  !# <accretionHalosMethod>
  !#  <unitName>Accretion_Halos_Null_Initialize</unitName>
  !# </accretionHalosMethod>
  subroutine Accretion_Halos_Null_Initialize(accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get &
       &,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get &
       &,Halo_Baryonic_Accretion_Rate_Abundances_Get,Halo_Baryonic_Accreted_Abundances_Get&
       &,Halo_Baryonic_Accretion_Rate_Chemicals_Get,Halo_Baryonic_Accreted_Chemicals_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                  ), intent(in   )          :: accretionHalosMethod
    procedure(Halo_Baryonic_Accretion_Rate_Null_Get           ), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Get
    procedure(Halo_Baryonic_Accreted_Mass_Null_Get            ), intent(inout), pointer :: Halo_Baryonic_Accreted_Mass_Get
    procedure(Halo_Baryonic_Failed_Accretion_Rate_Null_Get    ), intent(inout), pointer :: Halo_Baryonic_Failed_Accretion_Rate_Get
    procedure(Halo_Baryonic_Failed_Accreted_Mass_Null_Get     ), intent(inout), pointer :: Halo_Baryonic_Failed_Accreted_Mass_Get
    procedure(Halo_Baryonic_Accretion_Rate_Abundances_Null_Get), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Abundances_Get
    procedure(Halo_Baryonic_Accreted_Abundances_Null_Get      ), intent(inout), pointer :: Halo_Baryonic_Accreted_Abundances_Get
    procedure(Halo_Baryonic_Accretion_Rate_Chemicals_Null_Get ), intent(inout), pointer :: Halo_Baryonic_Accretion_Rate_Chemicals_Get
    procedure(Halo_Baryonic_Accreted_Chemicals_Null_Get       ), intent(inout), pointer :: Halo_Baryonic_Accreted_Chemicals_Get

    if (accretionHalosMethod == 'null') then
       ! Set pointers to our implementations of accretion functions.
       Halo_Baryonic_Accretion_Rate_Get            => Halo_Baryonic_Accretion_Rate_Null_Get
       Halo_Baryonic_Accreted_Mass_Get             => Halo_Baryonic_Accreted_Mass_Null_Get
       Halo_Baryonic_Failed_Accretion_Rate_Get     => Halo_Baryonic_Failed_Accretion_Rate_Null_Get
       Halo_Baryonic_Failed_Accreted_Mass_Get      => Halo_Baryonic_Failed_Accreted_Mass_Null_Get
       Halo_Baryonic_Accretion_Rate_Abundances_Get => Halo_Baryonic_Accretion_Rate_Abundances_Null_Get
       Halo_Baryonic_Accreted_Abundances_Get       => Halo_Baryonic_Accreted_Abundances_Null_Get
       Halo_Baryonic_Accretion_Rate_Chemicals_Get  => Halo_Baryonic_Accretion_Rate_Chemicals_Null_Get
       Halo_Baryonic_Accreted_Chemicals_Get        => Halo_Baryonic_Accreted_Chemicals_Null_Get
    end if
    return
  end subroutine Accretion_Halos_Null_Initialize

  double precision function Halo_Baryonic_Accretion_Rate_Null_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Halo_Baryonic_Accretion_Rate_Null_Get=0.0d0
    return
  end function Halo_Baryonic_Accretion_Rate_Null_Get

  double precision function Halo_Baryonic_Accreted_Mass_Null_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Halo_Baryonic_Accreted_Mass_Null_Get=0.0d0
    return
  end function Halo_Baryonic_Accreted_Mass_Null_Get

  double precision function Halo_Baryonic_Failed_Accretion_Rate_Null_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Halo_Baryonic_Failed_Accretion_Rate_Null_Get=0.0d0
    return
  end function Halo_Baryonic_Failed_Accretion_Rate_Null_Get

  double precision function Halo_Baryonic_Failed_Accreted_Mass_Null_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Halo_Baryonic_Failed_Accreted_Mass_Null_Get=0.0d0
    return
  end function Halo_Baryonic_Failed_Accreted_Mass_Null_Get

  subroutine Halo_Baryonic_Accretion_Rate_Abundances_Null_Get(thisNode,accretionRateAbundances)
    !% Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type(treeNode  ), intent(inout), pointer :: thisNode
    type(abundances), intent(inout)          :: accretionRateAbundances

    accretionRateAbundances=zeroAbundances
    return
  end subroutine Halo_Baryonic_Accretion_Rate_Abundances_Null_Get

  subroutine Halo_Baryonic_Accreted_Abundances_Null_Get(thisNode,accretedAbundances)
    !% Computes the mass of abundances accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type(treeNode  ), intent(inout), pointer :: thisNode
    type(abundances), intent(inout)          :: accretedAbundances

    accretedAbundances=zeroAbundances
    return
  end subroutine Halo_Baryonic_Accreted_Abundances_Null_Get

  subroutine Halo_Baryonic_Accretion_Rate_Chemicals_Null_Get(thisNode,accretionRateChemicals)
    !% Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type(treeNode          ), intent(inout), pointer :: thisNode
    type(chemicalAbundances), intent(inout)          :: accretionRateChemicals

    accretionRateChemicals=zeroChemicals
    return
  end subroutine Halo_Baryonic_Accretion_Rate_Chemicals_Null_Get

  subroutine Halo_Baryonic_Accreted_Chemicals_Null_Get(thisNode,accretedChemicals)
    !% Computes the mass of chemicals accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type(treeNode          ), intent(inout), pointer :: thisNode
    type(chemicalAbundances), intent(inout)          :: accretedChemicals

    accretedChemicals=zeroChemicals
    return
  end subroutine Halo_Baryonic_Accreted_Chemicals_Null_Get

end module Accretion_Halos_Null

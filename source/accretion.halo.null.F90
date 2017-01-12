!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !% An implementation of null accretion from the \gls{igm} onto halos.

  !# <accretionHalo name="accretionHaloNull">
  !#  <description>Accretion onto halos assuming no accretion.</description>
  !# </accretionHalo>

  type, extends(accretionHaloClass) :: accretionHaloNull
     !% A halo accretion class that assumes no accretion.
     private
   contains
     procedure :: branchHasBaryons       => nullBranchHasBaryons
     procedure :: accretionRate          => nullAccretionRate
     procedure :: accretedMass           => nullAccretedMass
     procedure :: failedAccretionRate    => nullFailedAccretionRate
     procedure :: failedAccretedMass     => nullFailedAccretedMass
     procedure :: accretionRateMetals    => nullAccretionRateMetals
     procedure :: accretedMassMetals     => nullAccretedMassMetals
     procedure :: accretionRateChemicals => nullAccretionRateChemicals
     procedure :: accretedMassChemicals  => nullAccretedMassChemicals
  end type accretionHaloNull

contains

  logical function nullBranchHasBaryons(self,node)
    !% Returns true if this branch can accrete any baryons.
    use Galacticus_Nodes
    implicit none
    class(accretionHaloNull), intent(inout)         :: self
    type (treeNode         ), intent(inout), target :: node
    !GCC$ attributes unused :: self, node
    
    nullBranchHasBaryons=.false.
    return
  end function nullBranchHasBaryons

  double precision function nullAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    implicit none
    class  (accretionHaloNull), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    integer                   , intent(in   ) :: accretionMode
    !GCC$ attributes unused :: self, node, accretionMode
    
    nullAccretionRate=0.0d0
    return
  end function nullAccretionRate

  double precision function nullAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    implicit none
    class  (accretionHaloNull), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    integer                   , intent(in   ) :: accretionMode
    !GCC$ attributes unused :: self, node, accretionMode

    nullAccretedMass=0.0d0
    return
  end function nullAccretedMass

  double precision function nullFailedAccretionRate(self,node,accretionMode)
    !% Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    implicit none
    class  (accretionHaloNull), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    integer                   , intent(in   ) :: accretionMode
    !GCC$ attributes unused :: self, node, accretionMode

    nullFailedAccretionRate=0.0d0
    return
  end function nullFailedAccretionRate

  double precision function nullFailedAccretedMass(self,node,accretionMode)
    !% Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    use Galacticus_Nodes
    implicit none
    class  (accretionHaloNull), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    integer                   , intent(in   ) :: accretionMode
    !GCC$ attributes unused :: self, node, accretionMode

    nullFailedAccretedMass=0.0d0
    return
  end function nullFailedAccretedMass

  function nullAccretionRateMetals(self,node,accretionMode)
    !% Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type   (abundances       )                :: nullAccretionRateMetals
    class  (accretionHaloNull), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    integer                   , intent(in   ) :: accretionMode
    !GCC$ attributes unused :: self, node, accretionMode

    nullAccretionRateMetals=zeroAbundances
    return
  end function nullAccretionRateMetals

  function nullAccretedMassMetals(self,node,accretionMode)
    !% Computes the mass of abundances accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    use Galacticus_Nodes
    implicit none
    type   (abundances       )                :: nullAccretedMassMetals
    class  (accretionHaloNull), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    integer                   , intent(in   ) :: accretionMode
    !GCC$ attributes unused :: self, node, accretionMode

    nullAccretedMassMetals=zeroAbundances
    return
  end function nullAccretedMassMetals

  function nullAccretionRateChemicals(self,node,accretionMode)
    !% Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type   (chemicalAbundances)                :: nullAccretionRateChemicals
    class  (accretionHaloNull ), intent(inout) :: self
    type   (treeNode          ), intent(inout) :: node
    integer                    , intent(in   ) :: accretionMode
    !GCC$ attributes unused :: self, node, accretionMode

    nullAccretionRateChemicals=zeroChemicalAbundances
    return
  end function nullAccretionRateChemicals

  function nullAccretedMassChemicals(self,node,accretionMode)
    !% Computes the mass of chemicals accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    use Galacticus_Nodes
    use Chemical_Abundances_Structure
    implicit none
    type   (chemicalAbundances)                :: nullAccretedMassChemicals
    class  (accretionHaloNull ), intent(inout) :: self
    type   (treeNode          ), intent(inout) :: node
    integer                    , intent(in   ) :: accretionMode
    !GCC$ attributes unused :: self, node, accretionMode

    nullAccretedMassChemicals=zeroChemicalAbundances
    return
  end function nullAccretedMassChemicals

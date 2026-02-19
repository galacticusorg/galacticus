!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  An implementation of zero accretion from the \gls{igm} onto halos.
  !!}

  !![
  <accretionHalo name="accretionHaloZero">
   <description>Accretion onto halos assuming no accretion.</description>
  </accretionHalo>
  !!]
  type, extends(accretionHaloClass) :: accretionHaloZero
     !!{
     A halo accretion class that assumes no accretion.
     !!}
     private
   contains
     procedure :: branchHasBaryons          => zeroBranchHasBaryons
     procedure :: accretionRate             => zeroAccretionRate
     procedure :: accretedMass              => zeroAccretedMass
     procedure :: failedAccretionRate       => zeroFailedAccretionRate
     procedure :: failedAccretedMass        => zeroFailedAccretedMass
     procedure :: accretionRateMetals       => zeroAccretionRateMetals
     procedure :: accretedMassMetals        => zeroAccretedMassMetals
     procedure :: failedAccretionRateMetals => zeroFailedAccretionRateMetals
     procedure :: failedAccretedMassMetals  => zeroFailedAccretedMassMetals
     procedure :: accretionRateChemicals    => zeroAccretionRateChemicals
     procedure :: accretedMassChemicals     => zeroAccretedMassChemicals
  end type accretionHaloZero

  interface accretionHaloZero
     !!{
     Constructors for the \refClass{accretionHaloZero} halo accretion class.
     !!}
     module procedure zeroConstructorParameters
  end interface accretionHaloZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{accretionHaloZero} halo accretion class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(accretionHaloZero)                :: self
    type(inputParameters  ), intent(inout) :: parameters

    self=accretionHaloZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  logical function zeroBranchHasBaryons(self,node)
    !!{
    Returns true if this branch can accrete any baryons.
    !!}
    implicit none
    class(accretionHaloZero), intent(inout)         :: self
    type (treeNode         ), intent(inout), target :: node
    !$GLC attributes unused :: self, node

    zeroBranchHasBaryons=.false.
    return
  end function zeroBranchHasBaryons

  double precision function zeroAccretionRate(self,node,accretionMode)
    !!{
    Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(accretionHaloZero           ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroAccretionRate=0.0d0
    return
  end function zeroAccretionRate

  double precision function zeroAccretedMass(self,node,accretionMode)
    !!{
    Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(accretionHaloZero           ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroAccretedMass=0.0d0
    return
  end function zeroAccretedMass

  double precision function zeroFailedAccretionRate(self,node,accretionMode)
    !!{
    Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(accretionHaloZero           ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroFailedAccretionRate=0.0d0
    return
  end function zeroFailedAccretionRate

  double precision function zeroFailedAccretedMass(self,node,accretionMode)
    !!{
    Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(accretionHaloZero           ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroFailedAccretedMass=0.0d0
    return
  end function zeroFailedAccretedMass

  function zeroAccretionRateMetals(self,node,accretionMode)
    !!{
    Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : abundances, zeroAbundances
    implicit none
    type (abundances                  )                :: zeroAccretionRateMetals
    class(accretionHaloZero           ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroAccretionRateMetals=zeroAbundances
    return
  end function zeroAccretionRateMetals

  function zeroAccretedMassMetals(self,node,accretionMode)
    !!{
    Computes the mass of abundances accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : abundances, zeroAbundances
    implicit none
    type (abundances                  )                :: zeroAccretedMassMetals
    class(accretionHaloZero           ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroAccretedMassMetals=zeroAbundances
    return
  end function zeroAccretedMassMetals

  function zeroFailedAccretionRateMetals(self,node,accretionMode)
    !!{
    Computes the rate of failed mass of abundance accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : abundances, zeroAbundances
    implicit none
    type   (abundances       )                :: zeroFailedAccretionRateMetals
    class  (accretionHaloZero), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    type (enumerationAccretionModeType)                   , intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroFailedAccretionRateMetals=zeroAbundances
    return
  end function zeroFailedAccretionRateMetals

  function zeroFailedAccretedMassMetals(self,node,accretionMode)
    !!{
    Computes the mass of abundances that failed to accrete (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : abundances, zeroAbundances
    implicit none
    type   (abundances       )                :: zeroFailedAccretedMassMetals
    class  (accretionHaloZero), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    type (enumerationAccretionModeType)                   , intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroFailedAccretedMassMetals=zeroAbundances
    return
  end function zeroFailedAccretedMassMetals

  function zeroAccretionRateChemicals(self,node,accretionMode)
    !!{
    Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances, zeroChemicalAbundances
    implicit none
    type (chemicalAbundances          )                :: zeroAccretionRateChemicals
    class(accretionHaloZero           ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroAccretionRateChemicals=zeroChemicalAbundances
    return
  end function zeroAccretionRateChemicals

  function zeroAccretedMassChemicals(self,node,accretionMode)
    !!{
    Computes the mass of chemicals accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances, zeroChemicalAbundances
    implicit none
    type (chemicalAbundances          )                :: zeroAccretedMassChemicals
    class(accretionHaloZero           ), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    type (enumerationAccretionModeType), intent(in   ) :: accretionMode
    !$GLC attributes unused :: self, node, accretionMode

    zeroAccretedMassChemicals=zeroChemicalAbundances
    return
  end function zeroAccretedMassChemicals

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which stores data for the standard disk node component.
!!}

module Node_Component_Disk_Standard_Data
  !!{RST
  Stores data for the standard disk node component.
  !!}
  use :: Galactic_Component_Truncations, only : galacticComponentTruncationClass
  use :: Galactic_Disk_Heatings        , only : galacticDiskHeatingClass
  use :: Mass_Distributions            , only : massDistributionClass           , kinematicsDistributionLocal, kinematicsDistributionClass
  use :: Object_Pools                  , only : objectPool
  implicit none
  public

  ! The mass distribution objects.
  class  (massDistributionClass           ), pointer :: massDistributionStellar_    , massDistributionGas_
  type   (kinematicsDistributionLocal     ), pointer :: kinematicDistribution_
  class  (kinematicsDistributionClass     ), pointer :: kinematicDistributionGas_
  !$omp threadprivate(massDistributionStellar_,massDistributionGas_,kinematicDistribution_,kinematicDistributionGas_)

  ! The model of truncation of the mass distribution, and the IDs of the meta-properties giving the degree of tidal truncation (of
  ! the stars) and the degree of truncation of the gas.
  class  (galacticComponentTruncationClass), pointer :: galacticComponentTruncation_
  !$omp threadprivate(galacticComponentTruncation_)
  integer                                            :: tidalTruncationID           , truncationGasID
  ! The model of heating of the mass distribution.
  class  (galacticDiskHeatingClass        ), pointer :: galacticDiskHeating_
  !$omp threadprivate(galacticDiskHeating_)

  ! Per-thread pools of the (lightweight) scaler mass distributions that wrap the dimensionless
  ! stellar and gas distributions above, re-used across nodes to avoid allocating and destroying a
  ! scaler on every call. A separate pool is used per wrapped distribution so the wrapped pointer
  ! never changes on re-use. Released in the component's thread uninitialization task.
  type   (objectPool                      )          :: scalerStellarPool          , scalerGasPool
  !$omp threadprivate(scalerStellarPool,scalerGasPool)

contains

  double precision function Node_Component_Disk_Standard_Truncation_Gas_Arrive(truncationGas,massGas,massGasArriving) result(truncationGasNew)
    !!{RST
    Return the truncation of the gas of a disk, given its current truncation, ``truncationGas``, and gas mass, ``massGas``, after
    a mass ``massGasArriving`` of gas (assumed to have the untruncated distribution) arrives in it, such that the untruncated mass,
    :math:`M_\mathrm{gas}/f_\mathrm{gas}`, grows by the mass of the arriving gas.
    !!}
    implicit none
    double precision, intent(in   ) :: truncationGas       , massGas       , &
         &                             massGasArriving
    double precision                :: fractionMassRetained, massGasCurrent

    truncationGasNew=truncationGas
    if (massGasArriving <= 0.0d0) return
    fractionMassRetained=exp(-max(0.0d0,truncationGas))
    if (fractionMassRetained >= 1.0d0) return
    massGasCurrent  =+max(0.0d0,massGas)
    truncationGasNew=-log(                                                        &
         &                 +(massGasCurrent                     +massGasArriving) &
         &                 /(massGasCurrent/fractionMassRetained+massGasArriving) &
         &               )
    return
  end function Node_Component_Disk_Standard_Truncation_Gas_Arrive

end module Node_Component_Disk_Standard_Data

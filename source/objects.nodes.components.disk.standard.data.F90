!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which stores data for the standard disk node component.
!!}

module Node_Component_Disk_Standard_Data
  !!{
  Stores data for the standard disk node component.
  !!}
  use :: Kind_Numbers      , only : kind_int8
  use :: Mass_Distributions, only : massDistributionClass, massDistributionCylindrical
  implicit none
  public

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8  )          :: lastUniqueID                    =-1
  !$omp threadprivate(lastUniqueID)
  ! Records of previously computed and stored quantities.
  logical                                     :: surfaceDensityCentralGasComputed   , surfaceDensityCentralStellarComputed, &
       &                                         surfaceDensityCentralTotalComputed
  !$omp threadprivate(surfaceDensityCentralGasComputed,surfaceDensityCentralStellarComputed,surfaceDensityCentralTotalComputed)
  double precision                            :: surfaceDensityCentralGas           , surfaceDensityCentralStellar        , &
       &                                         surfaceDensityCentralTotal
  !$omp threadprivate(surfaceDensityCentralGas,surfaceDensityCentralStellar,surfaceDensityCentralTotal)
  logical                                     :: radiusScaleDiskComputed
  !$omp threadprivate(radiusScaleDiskComputed)
  double precision                            :: radiusScaleDisk
  !$omp threadprivate(radiusScaleDisk)

  ! The mass distribution object.
  class           (massDistributionClass      ), pointer :: diskMassDistribution_
  class           (massDistributionCylindrical), pointer :: diskMassDistribution
  !$omp threadprivate(diskMassDistribution_,diskMassDistribution)

contains

  subroutine Node_Component_Disk_Standard_Reset(uniqueID)
    !!{
    Reset calculations for the standard disk component.
    !!}
    implicit none
    integer(kind=kind_int8), intent(in   ) :: uniqueID

    radiusScaleDiskComputed             =.false.
    surfaceDensityCentralGasComputed    =.false.
    surfaceDensityCentralStellarComputed=.false.
    surfaceDensityCentralTotalComputed  =.false.
    lastUniqueID                        =uniqueID
    return
  end subroutine Node_Component_Disk_Standard_Reset

end module Node_Component_Disk_Standard_Data

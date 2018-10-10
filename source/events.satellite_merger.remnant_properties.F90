!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which computes remnant properties for satellite mergers.

module Satellite_Merging_Remnant_Properties
  !% Computes remnant properties for satellite mergers.
  implicit none
  private
  public :: Satellite_Merging_Remnant_Compute

  ! Remnant properties.
  integer         , public :: destinationGasSatellite       , destinationStarsSatellite, &
       &                      destinationGasHost            , destinationStarsHost
  logical         , public :: mergerIsMajor
  double precision, public :: radiusRemnant                 , velocityCircularRemnant  , &
       &                      angularMomentumSpecificRemnant
  !$omp threadprivate(destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor,radiusRemnant,velocityCircularRemnant,angularMomentumSpecificRemnant)
  
contains

  !# <satelliteMergerTask>
  !#  <unitName>Satellite_Merging_Remnant_Compute</unitName>
  !# </satelliteMergerTask>
  subroutine Satellite_Merging_Remnant_Compute(node)
    !% Compute properties for satellite merger remnants.
    use Satellite_Merging_Mass_Movements
    use Satellite_Merging_Remnant_Sizes
    use Galacticus_Nodes
    implicit none
    type (treeNode                ), intent(inout), pointer :: node
    class(mergerMassMovementsClass)               , pointer :: mergerMassMovements_
    class(mergerRemnantSizeClass  )               , pointer :: mergerRemnantSize_

    ! Movement of mass and size calculations for the merger remnant are done once, here, so that they are determined by the
    ! properties of the merge target prior to any modification that will occur as node components are modified in response to the
    ! merger.
    mergerMassMovements_ => mergerMassMovements()
    mergerRemnantSize_   => mergerRemnantSize  ()
    call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    call mergerRemnantSize_  %get(node,radiusRemnant,velocityCircularRemnant,angularMomentumSpecificRemnant)
    return
  end subroutine Satellite_Merging_Remnant_Compute

end module Satellite_Merging_Remnant_Properties

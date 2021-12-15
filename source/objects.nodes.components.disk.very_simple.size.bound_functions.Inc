!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains custom functions for the very simple size disk component.
!!}

double precision function Node_Component_Disk_Very_Simple_Size_Half_Mass_Radius(self)
  !!{
  Return the half-mass radius of the very simple size disk.
  !!}
  implicit none
  class           (nodeComponentDiskVerySimpleSize), intent(inout) :: self

  Node_Component_Disk_Very_Simple_Size_Half_Mass_Radius=self%radius()
  return
end function Node_Component_Disk_Very_Simple_Size_Half_Mass_Radius

double precision function Node_Component_Disk_Very_Simple_Size_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the mass within a given radius for an standard disk.
  !!}
  use :: Galactic_Structure_Options               , only : componentTypeAll    , componentTypeDisk, massTypeAll    , massTypeBaryonic, &
          &                                                massTypeGalactic    , massTypeGaseous  , massTypeStellar, radiusLarge     , &
          &                                                weightByLuminosity  , weightByMass
  use :: Node_Component_Disk_Very_Simple_Size_Data, only : diskMassDistribution
  implicit none
  class           (nodeComponentDiskVerySimpleSize), intent(inout) :: self
  integer                                          , intent(in   ) :: componentType   , massType        , weightBy, weightIndex
  double precision                                 , intent(in   ) :: radius
  double precision                                                 :: radiusDisk      , fractionalRadius
  type            (stellarLuminosities            ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Very_Simple_Size_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Get the total mass.
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        Node_Component_Disk_Very_Simple_Size_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous)
        Node_Component_Disk_Very_Simple_Size_Enclosed_Mass=self%massGas()
     case (massTypeStellar)
        Node_Component_Disk_Very_Simple_Size_Enclosed_Mass=               self%massStellar()
     end select
  case (weightByLuminosity)
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Very_Simple_Size_Enclosed_Mass   =luminositiesDisk%luminosity(weightIndex)
     end select
  end select
  ! Return if no mass.
  if (Node_Component_Disk_Very_Simple_Size_Enclosed_Mass <=       0.0d0) return
  ! Return if the total mass was requested.
  if (radius                                     >= radiusLarge) return
  ! Compute the actual mass.
  radiusDisk=self%radius()
  if (radiusDisk > 0.0d0) then
     fractionalRadius=radius/radiusDisk
     Node_Component_Disk_Very_Simple_Size_Enclosed_Mass=                         &
          & +Node_Component_Disk_Very_Simple_Size_Enclosed_Mass                  &
          & *diskMassDistribution%massEnclosedBySphere(fractionalRadius)
  end if
  return
end function Node_Component_Disk_Very_Simple_Size_Enclosed_Mass

double precision function Node_Component_Disk_Very_Simple_Size_Surface_Density(self,positionCylindrical,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the surface density at a given position for an standard disk.
  !!}
  use :: Coordinates                              , only : coordinateCylindrical
  use :: Galactic_Structure_Options               , only : componentTypeAll                          , componentTypeDisk         , massTypeAll                       , massTypeBaryonic            , &
          &                                                massTypeGalactic                          , massTypeGaseous           , massTypeStellar                   , weightByLuminosity          , &
          &                                                weightByMass
  use :: Mass_Distributions                       , only : massDistributionCylindrical
  use :: Node_Component_Disk_Very_Simple_Size_Data, only : Node_Component_Disk_Very_Simple_Size_Reset, diskMassDistribution      , lastUniqueID                      , radiusScaleDisk             , &
          &                                                radiusScaleDiskComputed                   , surfaceDensityCentralGas  , surfaceDensityCentralGasComputed  , surfaceDensityCentralStellar, &
          &                                                surfaceDensityCentralStellarComputed      , surfaceDensityCentralTotal, surfaceDensityCentralTotalComputed
  use :: Numerical_Constants_Math                 , only : Pi
  implicit none
  class           (nodeComponentDiskVerySimpleSize), intent(inout) :: self
  integer                                          , intent(in   ) :: componentType         , massType   , &
       &                                                              weightBy              , weightIndex
  double precision                                 , intent(in   ) :: positionCylindrical(3)
  type            (treeNode                       ), pointer       :: selfNode
  type            (stellarLuminosities            ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  type            (coordinateCylindrical          )                :: position

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Very_Simple_Size_Surface_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Check whether this is a new node.
  selfNode => self%host()
  if (selfNode%uniqueID() /= lastUniqueID) call Node_Component_Disk_Very_Simple_Size_Reset(selfNode%uniqueID())
  ! Determine disk radius.
  if (.not.radiusScaleDiskComputed) then
     radiusScaleDisk        =self%radius()
     radiusScaleDiskComputed=.true.
  end if
  ! Return zero if the disk has unphysical size.
  if (radiusScaleDisk <= 0.0d0) then
     Node_Component_Disk_Very_Simple_Size_Surface_Density=0.0d0
     return
  end if
    ! Determine mass type.
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        if (.not.surfaceDensityCentralTotalComputed  ) then
           surfaceDensityCentralTotal          =(self%massGas()+self%massStellar())/radiusScaleDisk**2
           surfaceDensityCentralTotalComputed  =.true.
        end if
        Node_Component_Disk_Very_Simple_Size_Surface_Density=surfaceDensityCentralTotal
     case (massTypeGaseous)
        if (.not.surfaceDensityCentralGasComputed    ) then
           surfaceDensityCentralGas            = self%massGas()                    /radiusScaleDisk**2
           surfaceDensityCentralGasComputed    =.true.
          end if
        Node_Component_Disk_Very_Simple_Size_Surface_Density=surfaceDensityCentralGas
     case (massTypeStellar)
        if (.not.surfaceDensityCentralStellarComputed) then
           surfaceDensityCentralStellar        =                self%massStellar() /radiusScaleDisk**2
           surfaceDensityCentralStellarComputed=.true.
        end if
        Node_Component_Disk_Very_Simple_Size_Surface_Density=surfaceDensityCentralStellar
     end select
  case (weightByLuminosity)
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Very_Simple_Size_Surface_Density=luminositiesDisk%luminosity(weightIndex)/2.0d0/Pi/radiusScaleDisk**2
     end select
  end select
  ! Return if no density.
  if (Node_Component_Disk_Very_Simple_Size_Surface_Density <= 0.0d0) return  
  ! Check that the mass distribution is cylindrical.
  select type (diskMassDistribution)
  class is (massDistributionCylindrical)
     ! Compute the surface density.
     call position%rSet(positionCylindrical(1)/self%radius())
     Node_Component_Disk_Very_Simple_Size_Surface_Density=         &
          & +Node_Component_Disk_Very_Simple_Size_Surface_Density  &
          & *diskMassDistribution%surfaceDensity(position)
  end select
  return
end function Node_Component_Disk_Very_Simple_Size_Surface_Density

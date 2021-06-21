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
Contains custom functions for the very simple disk component.
!!}

subroutine Node_Component_Disk_Very_Simple_Attach_Pipe(self)
  !!{
  Attach cooling pipes to the very simple disk component.
  !!}
  use :: Galacticus_Error, only : Galacticus_Error_Report
  implicit none
  class(nodeComponentDiskVerySimple), intent(in   ) :: self
  type (nodeComponentHotHalo       )                :: hotHalo
  !$GLC attributes unused :: self

  if     (                                                                                                 &
       &   hotHalo%      hotHaloCoolingMassRateIsAttached()                                                &
       &  .or.                                                                                             &
       &   hotHalo%hotHaloCoolingAbundancesRateIsAttached()                                                &
       & )                                                                                                 &
       & call Galacticus_Error_Report(                                                                     &
       &                              'expected to find unclaimed hot halo mass/abundances cooling pipe'// &
       &                              {introspection:location}                                             &
       &                             )
  call hotHalo%      hotHaloCoolingMassRateFunction(      DiskVerySimpleMassGasRateGeneric)
  call hotHalo%hotHaloCoolingAbundancesRateFunction(DiskVerySimpleAbundancesGasRateGeneric)
  return
end subroutine Node_Component_Disk_Very_Simple_Attach_Pipe

double precision function Node_Component_Disk_Very_Simple_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the mass within a given radius for an very simple disk.
  !!}
  use :: Galactic_Structure_Options, only : componentTypeAll       , componentTypeDisk, massTypeAll    , massTypeBaryonic, &
          &                                 massTypeGalactic       , massTypeGaseous  , massTypeStellar, radiusLarge     , &
          &                                 weightByMass
  use :: Galacticus_Error          , only : Galacticus_Error_Report
  implicit none
  class           (nodeComponentDiskVerySimple), intent(inout) :: self
  integer                                      , intent(in   ) :: componentType, massType, weightBy, weightIndex
  double precision                             , intent(in   ) :: radius
  !$GLC attributes unused :: weightIndex

  ! Return zero mass if the mass and component types do not match.
  Node_Component_Disk_Very_Simple_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Determine mass contributed.
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        Node_Component_Disk_Very_Simple_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous)
        Node_Component_Disk_Very_Simple_Enclosed_Mass=self%massGas()
     case (massTypeStellar)
        Node_Component_Disk_Very_Simple_Enclosed_Mass=               self%massStellar()
     end select
  case default
     call Galacticus_Error_Report('this component does not track luminosity'//{introspection:location})
  end select
  ! Return if no mass.
  if (Node_Component_Disk_Very_Simple_Enclosed_Mass <=       0.0d0) return
  ! Return if the total mass was requested.
  if (radius                                        >= radiusLarge) return
  ! Otherwise we have an error.
  call Galacticus_Error_Report('this component does not specify a mass profile'//{introspection:location})
  return
end function Node_Component_Disk_Very_Simple_Enclosed_Mass

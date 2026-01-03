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

!!{
Contains custom functions for the very simple disk component.
!!}

double precision function Node_Component_Disk_Very_Simple_Mass_Baryonic(self) result(massBaryonic)
  !!{
  Return the baryonic mass for the very simple disk component.
  !!}
  implicit none
  class(nodeComponentDiskVerySimple), intent(inout) :: self

  massBaryonic=+max(0.0d0,self%massStellar()) &
       &       +max(0.0d0,self%massGas    ())
  return
end function Node_Component_Disk_Very_Simple_Mass_Baryonic

double precision function Node_Component_Disk_Very_Simple_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the mass within a given radius for an very simple disk.
  !!}
  use :: Galactic_Structure_Options, only : componentTypeAll, componentTypeDisk           , massTypeAll            , massTypeBaryonic       , &
          &                                 massTypeGalactic, massTypeGaseous             , massTypeStellar        , radiusLarge            , &
          &                                 weightByMass    , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Error                     , only : Error_Report
  implicit none
  class           (nodeComponentDiskVerySimple ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: radius
  !$GLC attributes unused :: weightIndex

  ! Return zero mass if the mass and component types do not match.
  Node_Component_Disk_Very_Simple_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Determine mass contributed.
  select case (weightBy%ID)
  case (weightByMass      %ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
        Node_Component_Disk_Very_Simple_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous%ID)
        Node_Component_Disk_Very_Simple_Enclosed_Mass=self%massGas()
     case (massTypeStellar%ID)
        Node_Component_Disk_Very_Simple_Enclosed_Mass=               self%massStellar()
     end select
  case default
     call Error_Report('this component does not track luminosity'//{introspection:location})
  end select
  ! Return if no mass.
  if (Node_Component_Disk_Very_Simple_Enclosed_Mass <=       0.0d0) return
  ! Return if the total mass was requested.
  if (radius                                        >= radiusLarge) return
  ! Otherwise we have an error.
  call Error_Report('this component does not specify a mass profile'//{introspection:location})
  return
end function Node_Component_Disk_Very_Simple_Enclosed_Mass

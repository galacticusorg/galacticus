!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

function Node_Component_Disk_Very_Simple_Size_Mass_Distribution(self) result(massDistribution_)
  !!{
  Return the mass distribution for the very simple size disk component.
  !!}
   use :: Mass_Distributions                      , only : massDistributionClass      , massDistributionCylindricalScaler, massDistributionComposite, massDistributionList, &
        &                                                  massDistributionCylindrical
  use :: Node_Component_Disk_Very_Simple_Size_Data, only : massDistributionStellar_   , massDistributionGas_
  use :: Galactic_Structure_Options               , only : componentTypeDisk          , massTypeStellar                  , massTypeGaseous
  implicit none
  class          (massDistributionClass            ), pointer       :: massDistribution_
  class          (nodeComponentDiskVerySimpleSize  ), intent(inout) :: self
  type           (massDistributionCylindricalScaler), pointer       :: massDistributionStellar   , massDistributionGas
  type           (massDistributionComposite        ), pointer       :: massDistributionTotal
  type           (massDistributionList             ), pointer       :: massDistributionComponents
  double precision                                                  :: massStellar               , massGas            , &
       &                                                               radiusScale

  ! Get properties of the mass distribution and ensure they are physical.
  massStellar=max(0.0d0,self%massStellar())
  massGas    =max(0.0d0,self%massGas    ())
  radiusScale=          self%radius     ()
  if (radiusScale <= 0.0d0) then
     ! Disk has non-positive size. Set an arbitrary, positive size, but also set masses to zero. This results in an empty mass
     ! distribution.
     radiusScale=1.0d0
     massStellar=0.0d0
     massGas    =0.0d0
  end if
  ! Wrap the dimensionless mass distribution inside scaler classes to allow us to re-scale it to any disk system, and then composite those.
  allocate(massDistributionStellar        )
  allocate(massDistributionGas            )
  allocate(massDistributionTotal          )
  allocate(massDistributionComponents     )
  allocate(massDistributionComponents%next)
  massDistributionComponents     %massDistribution_ => massDistributionStellar
  massDistributionComponents%next%massDistribution_ => massDistributionGas
  select type (massDistributionStellar_)
  class is (massDistributionCylindrical)
     !![
     <referenceConstruct object="massDistributionStellar" constructor="massDistributionCylindricalScaler(factorScalingLength=radiusScale,factorScalingMass=massStellar,massDistribution_=massDistributionStellar_  )"/>
     !!]
  end select
  select type (massDistributionGas_    )
  class is (massDistributionCylindrical)
     !![
     <referenceConstruct object="massDistributionGas"     constructor="massDistributionCylindricalScaler(factorScalingLength=radiusScale,factorScalingMass=massGas    ,massDistribution_=massDistributionGas_      )"/>
     !!]
  end select
  !![
  <referenceConstruct    object="massDistributionTotal"   constructor="massDistributionComposite        (                                                                                massDistributionComponents)"/>
  <objectDestructor name="massDistributionStellar"/>
  <objectDestructor name="massDistributionGas"    />
  !!]
  nullify(massDistributionComponents)
  ! Return a pointer to the disk mass distribution.
  massDistribution_ => massDistributionTotal
  return
end function Node_Component_Disk_Very_Simple_Size_Mass_Distribution

double precision function Node_Component_Disk_Very_Simple_Size_Half_Mass_Radius(self) result(radiusHalfMass)
  !!{
  Return the half-mass radius of the very simple size disk.
  !!}
  use :: Error                                    , only : Error_Report
  use :: Mass_Distributions                       , only : massDistributionCylindrical
  use :: Node_Component_Disk_Very_Simple_Size_Data, only : massDistributionStellar_
  implicit none
  class(nodeComponentDiskVerySimpleSize), intent(inout) :: self

  select type (massDistributionStellar_)
  class is (massDistributionCylindrical)
     radiusHalfMass=+massDistributionStellar_%radiusHalfMass() &
          &         *self                    %radius        ()
  class default
     radiusHalfMass=0.0d0
     call Error_Report('disk mass distribution is not cylindrically-symmetric'//{introspection:location})
  end select
  return
end function Node_Component_Disk_Very_Simple_Size_Half_Mass_Radius

double precision function Node_Component_Disk_Very_Simple_Size_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex) result(massEnclosed)
  !!{
  Computes the mass within a given radius for a very simple size disk.
  !!}
  use :: Galactic_Structure_Options, only : weightByLuminosity   , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskVerySimpleSize), intent(inout) :: self
  type            (enumerationComponentTypeType   ), intent(in   ) :: componentType
  type            (enumerationMassTypeType        ), intent(in   ) :: massType
  type            (enumerationWeightByType        ), intent(in   ) :: weightBy
  integer                                          , intent(in   ) :: weightIndex
  double precision                                 , intent(in   ) :: radius
  class           (massDistributionClass          ), pointer       :: massDistribution_
  type            (stellarLuminosities            ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)

  massDistribution_ => self             %massDistribution    (                             )
  massEnclosed      =  massDistribution_%massEnclosedBySphere(radius,componentType,massType)
  if (massEnclosed > 0.0d0 .and. weightBy     == weightByLuminosity) then
     luminositiesDisk= self%luminositiesStellar()
     massEnclosed    =+massEnclosed                                         &
          &           *luminositiesDisk %luminosity(weightIndex           ) &
          &           /massDistribution_%massTotal (componentType,massType)
  end if
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Disk_Very_Simple_Size_Enclosed_Mass

double precision function Node_Component_Disk_Very_Simple_Size_Surface_Density(self,positionCylindrical,componentType,massType,weightBy,weightIndex) result(surfaceDensity)
  !!{
  Computes the surface density at a given position for a very simple size disk.
  !!}
  use :: Coordinates               , only : assignment(=)               , coordinateCylindrical
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType, weightByLuminosity 
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskVerySimpleSize), intent(inout) :: self
  type            (enumerationComponentTypeType   ), intent(in   ) :: componentType
  type            (enumerationMassTypeType        ), intent(in   ) :: massType
  type            (enumerationWeightByType        ), intent(in   ) :: weightBy
  integer                                          , intent(in   ) :: weightIndex
  double precision                                 , intent(in   ) :: positionCylindrical(3)
  class           (massDistributionClass          ), pointer       :: massDistribution_
  type            (stellarLuminosities            ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  type            (coordinateCylindrical          )                :: position

  position          =  positionCylindrical
  massDistribution_ => self               %massDistribution(                               )
  surfaceDensity    =  massDistribution_  %surfaceDensity  (position,componentType,massType)
  if (surfaceDensity > 0.0d0 .and. weightBy == weightByLuminosity) then
     luminositiesDisk= self%luminositiesStellar()
     surfaceDensity  =+surfaceDensity                                       &
          &           *luminositiesDisk %luminosity(weightIndex           ) &
          &           /massDistribution_%massTotal (componentType,massType)
  end if
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Disk_Very_Simple_Size_Surface_Density

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
Contains custom functions for the standard disk component.
!!}

subroutine Node_Component_Disk_Standard_Attach_Pipes(self)
  !!{
  Attach cooling pipes to the standard disk component.
  !!}
  use :: Error, only : Error_Report
  implicit none
  class(nodeComponentDiskStandard), intent(in   ) :: self
  type (nodeComponentHotHalo     )                :: hotHalo
  !$GLC attributes unused :: self

  if (hotHalo%hotHaloCoolingMassRateIsAttached           ()) &
       call Error_Report('expected to find unclaimed hot halo mass cooling pipe'            //{introspection:location})
  if (hotHalo%hotHaloCoolingAngularMomentumRateIsAttached()) &
       call Error_Report('expected to find unclaimed hot halo angular momentum cooling pipe'//{introspection:location})
  if (hotHalo%hotHaloCoolingAbundancesRateIsAttached     ()) &
       call Error_Report('expected to find unclaimed hot halo abundances cooling pipe'      //{introspection:location})
  call hotHalo%hotHaloCoolingMassRateFunction           (DiskStandardMassGasRateGeneric        )
  call hotHalo%hotHaloCoolingAngularMomentumRateFunction(DiskStandardAngularMomentumRateGeneric)
  call hotHalo%hotHaloCoolingAbundancesRateFunction     (DiskStandardAbundancesGasRateGeneric  )
  return
end subroutine Node_Component_Disk_Standard_Attach_Pipes

function Node_Component_Disk_Standard_Mass_Distribution(self) result(massDistribution_)
  !!{
  Return the mass distribution for the standard disk component.
  !!}
   use :: Mass_Distributions              , only : massDistributionClass      , massDistributionCylindricalScaler, massDistributionComposite, massDistributionList, &
        &                                          massDistributionCylindrical
  use :: Node_Component_Disk_Standard_Data, only : massDistributionStellar_   , massDistributionGas_
  use :: Galactic_Structure_Options       , only : componentTypeDisk          , massTypeStellar                  , massTypeGaseous
  implicit none
  class          (massDistributionClass            ), pointer       :: massDistribution_
  class          (nodeComponentDiskStandard        ), intent(inout) :: self
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
end function Node_Component_Disk_Standard_Mass_Distribution

double precision function Node_Component_Disk_Standard_Half_Mass_Radius(self) result(radiusHalfMass)
  !!{
  Return the half-mass radius of the standard disk.
  !!}
  use :: Error                            , only : Error_Report
  use :: Mass_Distributions               , only : massDistributionCylindrical
  use :: Node_Component_Disk_Standard_Data, only : massDistributionStellar_
  implicit none
  class(nodeComponentDiskStandard), intent(inout) :: self

  select type (massDistributionStellar_)
  class is (massDistributionCylindrical)
     radiusHalfMass=+massDistributionStellar_%radiusHalfMass() &
          &         *self                    %radius        ()
  class default
     radiusHalfMass=0.0d0
     call Error_Report('disk mass distribution is not cylindrically-symmetric'//{introspection:location})
  end select
  return
end function Node_Component_Disk_Standard_Half_Mass_Radius

double precision function Node_Component_Disk_Standard_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex) result(massEnclosed)
  !!{
  Computes the mass within a given radius for a standard disk.
  !!}
  use :: Galactic_Structure_Options, only : weightByLuminosity   , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: radius
  class           (massDistributionClass       ), pointer       :: massDistribution_
  type            (stellarLuminosities         ), save          :: luminositiesDisk
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
end function Node_Component_Disk_Standard_Enclosed_Mass

function Node_Component_Disk_Standard_Acceleration(self,positionCartesian,componentType,massType) result(acceleration)
  !!{
  Computes the gravitational acceleration at a given position for a standard disk.
  !!}
  use :: Coordinates               , only : assignment(=)               , coordinateCartesian
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  double precision                                             , dimension(3) :: acceleration
  class           (nodeComponentDiskStandard   ), intent(inout)               :: self
  double precision                              , intent(in   ), dimension(3) :: positionCartesian
  type            (enumerationComponentTypeType), intent(in   )               :: componentType
  type            (enumerationMassTypeType     ), intent(in   )               :: massType
  class           (massDistributionClass       ), pointer                     :: massDistribution_
  type            (coordinateCartesian         )                              :: position

  position          =  positionCartesian
  massDistribution_ => self             %massDistribution(                               )
  acceleration      =  massDistribution_%acceleration    (position,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Disk_Standard_Acceleration

function Node_Component_Disk_Standard_Tidal_Tensor(self,positionCartesian,componentType,massType) result(tidalTensor)
  !!{
  Computes the gravitational acceleration at a given position for a standard disk.
  !!}
  use :: Coordinates               , only : assignment(=)                 , coordinateCartesian
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType  , enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Tensors                   , only : tensorRank2Dimension3Symmetric
  implicit none
  type            (tensorRank2Dimension3Symmetric)                              :: tidalTensor
  class           (nodeComponentDiskStandard     ), intent(inout)               :: self
  double precision                                , intent(in   ), dimension(3) :: positionCartesian
  type            (enumerationComponentTypeType  ), intent(in   )               :: componentType
  type            (enumerationMassTypeType       ), intent(in   )               :: massType
  class           (massDistributionClass         ), pointer                     :: massDistribution_
  type            (coordinateCartesian           )                              :: position

  position          =  positionCartesian
  massDistribution_ => self             %massDistribution(                               )
  tidalTensor       =  massDistribution_%tidalTensor     (position,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Disk_Standard_Tidal_Tensor

double precision function Node_Component_Disk_Standard_Density(self,positionSpherical,componentType,massType,weightBy,weightIndex) result(density)
  !!{
  Computes the density at a given position for a standard disk.
  !!}
  use :: Coordinates               , only : assignment(=)                 , coordinateSpherical
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType  , enumerationMassTypeType, enumerationWeightByType, weightByLuminosity
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: positionSpherical(3)
  class           (massDistributionClass       ), pointer       :: massDistribution_
  type            (stellarLuminosities         ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  type            (coordinateSpherical         )                :: position

  position          =  positionSpherical
  massDistribution_ => self             %massDistribution(                               )
  density           =  massDistribution_%density         (position,componentType,massType)
  if (density > 0.0d0 .and. weightBy == weightByLuminosity) then
     luminositiesDisk= self%luminositiesStellar()
     density         =+density                                              &
          &           *luminositiesDisk %luminosity(weightIndex           ) &
          &           /massDistribution_%massTotal (componentType,massType)
  end if  
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Disk_Standard_Density

double precision function Node_Component_Disk_Standard_Density_Spherical_Average(self,radius,componentType,massType,weightBy,weightIndex) result(density)
  !!{
  Computes the density at a given position for a standard disk.
  !!}
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType, weightByLuminosity
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: radius
  type            (stellarLuminosities         ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  class           (massDistributionClass       ), pointer       :: massDistribution_

  massDistribution_ => self             %massDistribution       (                             )
  density           =  massDistribution_%densitySphericalAverage(radius,componentType,massType)
  if (density > 0.0d0 .and. weightBy == weightByLuminosity) then
     luminositiesDisk= self%luminositiesStellar()
     density         =+density                                              &
          &           *luminositiesDisk %luminosity(weightIndex           ) &
          &           /massDistribution_%massTotal (componentType,massType)
  end if
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
 end function Node_Component_Disk_Standard_Density_Spherical_Average

double precision function Node_Component_Disk_Standard_Potential(self,radius,componentType,massType,status) result(potential)
  !!{
  Compute the gravitational potential due to a standard disk.
  !!}
  use :: Coordinates               , only : assignment(=)               , coordinateCylindrical
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationStructureErrorCodeType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskStandard        ), intent(inout)           :: self
  type            (enumerationComponentTypeType     ), intent(in   )           :: componentType
  type            (enumerationMassTypeType          ), intent(in   )           :: massType
  double precision                                   , intent(in   )           :: radius
  type            (enumerationStructureErrorCodeType), intent(inout), optional :: status
  class           (massDistributionClass            ), pointer                 :: massDistribution_
  type            (coordinateCylindrical            )                          :: position
  !$GLC attributes unused :: status

  position          =  [radius,0.0d0,0.0d0]
  massDistribution_ => self             %massDistribution(                               )
  potential         =  massDistribution_%potential       (position,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Disk_Standard_Potential

double precision function Node_Component_Disk_Standard_Rotation_Curve(self,radius,componentType,massType) result(rotationCurve)
  !!{
  Computes the rotation curve at a given radius for a standard disk.
  !!}
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  double precision                              , intent(in   ) :: radius
  class           (massDistributionClass       ), pointer       :: massDistribution_

  massDistribution_ => self             %massDistribution(                             )
  rotationCurve     =  massDistribution_%rotationCurve   (radius,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Disk_Standard_Rotation_Curve

double precision function Node_Component_Disk_Standard_Rotation_Curve_Gradient(self,radius,componentType,massType) result(rotationCurveGradient)
  !!{
  Computes the rotation curve gradient for a standard disk.
  !!}
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  double precision                              , intent(in   ) :: radius
  class           (massDistributionClass       ), pointer       :: massDistribution_

  massDistribution_     => self             %massDistribution     (                             )
  rotationCurveGradient =  massDistribution_%rotationCurveGradient(radius,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Disk_Standard_Rotation_Curve_Gradient

double precision function Node_Component_Disk_Standard_Surface_Density(self,positionCylindrical,componentType,massType,weightBy,weightIndex) result(surfaceDensity)
  !!{
  Computes the surface density at a given position for a standard disk.
  !!}
  use :: Coordinates               , only : assignment(=)               , coordinateCylindrical
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType, weightByLuminosity 
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: positionCylindrical(3)
  class           (massDistributionClass       ), pointer       :: massDistribution_
  type            (stellarLuminosities         ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  type            (coordinateCylindrical       )                :: position

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
end function Node_Component_Disk_Standard_Surface_Density

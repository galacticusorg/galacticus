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
       call Error_Report('expected to find unclaimed hot halo mass cooling pipe'//{introspection:location})
  if (hotHalo%hotHaloCoolingAngularMomentumRateIsAttached()) &
       call Error_Report('expected to find unclaimed hot halo angular momentum cooling pipe'//{introspection:location})
  if (hotHalo%hotHaloCoolingAbundancesRateIsAttached     ()) &
       call Error_Report('expected to find unclaimed hot halo abundances cooling pipe'//{introspection:location})
  call hotHalo%hotHaloCoolingMassRateFunction           (DiskStandardMassGasRateGeneric        )
  call hotHalo%hotHaloCoolingAngularMomentumRateFunction(DiskStandardAngularMomentumRateGeneric)
  call hotHalo%hotHaloCoolingAbundancesRateFunction     (DiskStandardAbundancesGasRateGeneric  )
  return
end subroutine Node_Component_Disk_Standard_Attach_Pipes

double precision function Node_Component_Disk_Standard_Half_Mass_Radius(self)
  !!{
  Return the half-mass radius of the standard disk.
  !!}
  use :: Error                            , only : Error_Report
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  implicit none
  class(nodeComponentDiskStandard), intent(inout) :: self

  Node_Component_Disk_Standard_Half_Mass_Radius=self%radius()*diskMassDistribution%radiusHalfMass()
  return
end function Node_Component_Disk_Standard_Half_Mass_Radius

double precision function Node_Component_Disk_Standard_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the mass within a given radius for an standard disk.
  !!}
  use :: Galactic_Structure_Options       , only : componentTypeAll    , componentTypeDisk, massTypeAll                 , massTypeBaryonic       , &
          &                                        massTypeGalactic    , massTypeGaseous  , massTypeStellar             , radiusLarge            , &
          &                                        weightByLuminosity  , weightByMass     , enumerationComponentTypeType, enumerationMassTypeType, &
          &                                        enumerationWeightByType
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: radius
  double precision                                              :: radiusDisk      , fractionalRadius
  type            (stellarLuminosities         ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Standard_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Get the total mass.
  select case (weightBy%ID)
  case (weightByMass      %ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
        Node_Component_Disk_Standard_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous%ID)
        Node_Component_Disk_Standard_Enclosed_Mass=self%massGas()
     case (massTypeStellar%ID)
        Node_Component_Disk_Standard_Enclosed_Mass=               self%massStellar()
     end select
  case (weightByLuminosity%ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Standard_Enclosed_Mass   =luminositiesDisk%luminosity(weightIndex)
     end select
  end select
  ! Return if no mass.
  if (Node_Component_Disk_Standard_Enclosed_Mass <=       0.0d0) return
  ! Return if the total mass was requested.
  if (radius                                     >= radiusLarge) return
  ! Compute the actual mass.
  radiusDisk=self%radius()
  if (radiusDisk > 0.0d0) then
     fractionalRadius=radius/radiusDisk
     Node_Component_Disk_Standard_Enclosed_Mass=                         &
          & +Node_Component_Disk_Standard_Enclosed_Mass                  &
          & *diskMassDistribution%massEnclosedBySphere(fractionalRadius)
  end if
  return
end function Node_Component_Disk_Standard_Enclosed_Mass

function Node_Component_Disk_Standard_Acceleration(self,positionCartesian,componentType,massType)
  !!{
  Computes the gravitational acceleration at a given position for a standard disk.
  !!}
  use :: Coordinates                      , only : assignment(=)       , coordinateCartesian
  use :: Galactic_Structure_Options       , only : componentTypeAll    , componentTypeDisk  , massTypeAll                    , massTypeBaryonic            , &
       &                                           massTypeGalactic    , massTypeGaseous    , massTypeStellar                , enumerationComponentTypeType, &
       &                                           enumerationMassTypeType
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  use :: Numerical_Constants_Astronomical , only : gigaYear            , megaParsec         , gravitationalConstantGalacticus
  use :: Numerical_Constants_Prefixes     , only : kilo
  implicit none
  double precision                                             , dimension(3) :: Node_Component_Disk_Standard_Acceleration
  class           (nodeComponentDiskStandard   ), intent(inout)               :: self
  double precision                              , intent(in   ), dimension(3) :: positionCartesian
  type            (enumerationComponentTypeType), intent(in   )               :: componentType
  type            (enumerationMassTypeType     ), intent(in   )               :: massType
  double precision                                                            :: radius                                   , massTotal
  type            (coordinateCartesian         )                              :: positionScaleFree

  ! Return if the disk component is not selected.
  Node_Component_Disk_Standard_Acceleration=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Determine total mass.
  select case (massType%ID)
  case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
     massTotal=self%massGas()+self%massStellar()
  case (massTypeGaseous%ID)
     massTotal=self%massGas()
  case (massTypeStellar%ID)
     massTotal=               self%massStellar()
  case default
     massTotal=0.0d0
  end select
  ! Return if no mass.
  if (massTotal <= 0.0d0 .or. self%radius() <= 0.0d0) return
  ! Compute the acceleration.
  radius           =sqrt(sum(positionCartesian**2))
  positionScaleFree=positionCartesian/self%radius()
  Node_Component_Disk_Standard_Acceleration=+kilo                                                    &
         &                                  *gigaYear                                                &
         &                                  /megaParsec                                              &
         &                                  *gravitationalConstantGalacticus                         &
         &                                  *massTotal                                               &
         &                                  /self%radius()                                       **2 &
         &                                  *diskMassDistribution%acceleration(positionScaleFree)
  return
end function Node_Component_Disk_Standard_Acceleration

function Node_Component_Disk_Standard_Tidal_Tensor(self,positionCartesian,componentType,massType)
  !!{
  Computes the gravitational acceleration at a given position for a standard disk.
  !!}
  use :: Coordinates                      , only : assignment(=)                  , coordinateCartesian
  use :: Galactic_Structure_Options       , only : componentTypeAll               , componentTypeDisk  , massTypeAll    , massTypeBaryonic            , &
       &                                           massTypeGalactic               , massTypeGaseous    , massTypeStellar, enumerationComponentTypeType, &
       &                                           enumerationMassTypeType
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  use :: Numerical_Constants_Astronomical , only : gravitationalConstantGalacticus
  use :: Tensors                          , only : tensorRank2Dimension3Symmetric , tensorNullR2D3Sym  , operator(*)
  implicit none
  type            (tensorRank2Dimension3Symmetric)                              :: Node_Component_Disk_Standard_Tidal_Tensor
  class           (nodeComponentDiskStandard     ), intent(inout)               :: self
  double precision                                , intent(in   ), dimension(3) :: positionCartesian
  type            (enumerationComponentTypeType  ), intent(in   )               :: componentType
  type            (enumerationMassTypeType       ), intent(in   )               :: massType
  double precision                                                              :: radius                                   , massTotal
  type            (coordinateCartesian           )                              :: positionScaleFree

  ! Return if the disk component is not selected.q
  Node_Component_Disk_Standard_Tidal_Tensor=tensorNullR2D3Sym
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Determine total mass.
  select case (massType%ID)
  case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
     massTotal=self%massGas()+self%massStellar()
  case (massTypeGaseous%ID)
     massTotal=self%massGas()
  case (massTypeStellar%ID)
     massTotal=               self%massStellar()
  case default
     massTotal=0.0d0
  end select
  ! Return if no mass.
  if (massTotal <= 0.0d0 .or. self%radius() <= 0.0d0) return
  ! Compute the acceleration.
  radius           =sqrt(sum(positionCartesian**2))
  positionScaleFree=positionCartesian/self%radius()
  Node_Component_Disk_Standard_Tidal_Tensor=+gravitationalConstantGalacticus                        &
         &                                  *massTotal                                              &
         &                                  /self%radius()                                      **3 &
         &                                  *diskMassDistribution%tidalTensor(positionScaleFree)
  return
end function Node_Component_Disk_Standard_Tidal_Tensor

function Node_Component_Disk_Standard_Chandrasekhar_Integral(self,positionCartesian,velocityCartesian,radiusHalfMass,componentType,massType)
  !!{
  Computes the gravitational acceleration at a given position for a standard disk.
  !!}
  use :: Galactic_Structure_Options      , only : componentTypeAll                 , componentTypeDisk           , massTypeAll            , weightByMass         , &
       &                                          weightIndexNull                  , enumerationComponentTypeType, enumerationMassTypeType
  use :: Numerical_Constants_Math        , only : Pi
  use :: Coordinates                     , only : assignment(=)                    , coordinateSpherical         , coordinateCartesian    , coordinateCylindrical
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  use :: Mass_Distributions              , only : massDistributionGaussianEllipsoid
  use :: Linear_Algebra                  , only : vector                           , matrix             , assignment(=)
  implicit none
  double precision                                                  , dimension(3) :: Node_Component_Disk_Standard_Chandrasekhar_Integral
  class           (nodeComponentDiskStandard        ), intent(inout)               :: self
  double precision                                   , intent(in   ), dimension(3) :: positionCartesian                                          , velocityCartesian
  double precision                                   , intent(in   )               :: radiusHalfMass
  type            (enumerationComponentTypeType     ), intent(in   )               :: componentType
  type            (enumerationMassTypeType          ), intent(in   )               :: massType
  double precision                                   , parameter                   :: toomreQRadiusHalfMass                              =1.50d0  ! The Toomre Q-parameter at the disk half-mass radius (Benson et al.,
                                                                                                                                                  ! 2004 , https://ui.adsabs.harvard.edu/abs/2004MNRAS.351.1215B, Appendix A).
  double precision                                   , parameter                   :: toomreQFactor                                      =3.36d0  ! The factor appearing in the definition of the Toomre Q-parameter for
                                                                                                                                                  ! a stellar disk (Binney & Tremaine, eqn. 6.71).
  double precision                                                  , dimension(3) :: velocityDisk                                               , velocityRelative                , &
       &                                                                              positionSpherical                                          , positionSphericalMidplane       , &
       &                                                                              positionCartesianMidplane                                  , positionCylindricalMidplane     , &
       &                                                                              positionCylindricalHalfMass
  type            (massDistributionGaussianEllipsoid), save                        :: velocityDistribution
  logical                                            , save                        :: velocityDistributionInitialized                    =.false.
  !$omp threadprivate(velocityDistribution,velocityDistributionInitialized)
  type            (coordinateSpherical              )                              :: coordinatesSpherical
  type            (coordinateCartesian              )                              :: coordinatesCartesian
  type            (coordinateCylindrical            )                              :: coordinatesCylindrical
  double precision                                                                 :: velocityDispersionRadial                                   , velocityDispersionAzimuthal     , &
       &                                                                              velocityDispersionVertical                                 , velocityCircular                , &
       &                                                                              velocityCircularHalfMassRadius                             , velocityCircularSquaredGradient , &
       &                                                                              velocityCircularSquaredGradientHalfMassRadius              , density                         , &
       &                                                                              densityMidPlane                                            , densitySurface                  , &
       &                                                                              heightScale                                                , radiusMidplane                  , &
       &                                                                              frequencyCircular                                          , frequencyEpicyclic              , &
       &                                                                              frequencyCircularHalfMassRadius                            , frequencyEpicyclicHalfMassRadius, &
       &                                                                              densitySurfaceRadiusHalfMass                               , velocityDispersionRadialHalfMass, &
       &                                                                              velocityDispersionMaximum                                  , velocityRelativeMagnitude       , &
       &                                                                              factorSuppressionExtendedMass
  type            (matrix                           )                              :: rotation

  ! Return if the disk component is not selected.
  Node_Component_Disk_Standard_Chandrasekhar_Integral=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk) .or. self%radius() <= 0.0d0) return
  ! Construct the velocity vector of the disk rotation.
  positionCartesianMidplane                    =[positionCartesian(1),positionCartesian(2),0.0d0]
  coordinatesCartesian                         = positionCartesian
  coordinatesSpherical                         = coordinatesCartesian
  positionSpherical                            = coordinatesSpherical
  coordinatesCartesian                         = positionCartesianMidplane
  coordinatesSpherical                         = coordinatesCartesian
  coordinatesCylindrical                       = coordinatesCartesian 
  positionSphericalMidplane                    = coordinatesSpherical
  positionCylindricalMidplane                  = coordinatesCylindrical
  positionCylindricalHalfMass                  =[self%halfMassRadius(),0.0d0,0.0d0]
  radiusMidplane                               = coordinatesCylindrical%r()
  velocityCircular                             =Node_Component_Disk_Standard_Rotation_Curve         (self,     radiusMidplane  ,componentType,massType)
  velocityCircularSquaredGradient              =Node_Component_Disk_Standard_Rotation_Curve_Gradient(self,     radiusMidplane  ,componentType,massType)
  velocityCircularHalfMassRadius               =Node_Component_Disk_Standard_Rotation_Curve         (self,self%halfMassRadius(),componentType,massType)
  velocityCircularSquaredGradientHalfMassRadius=Node_Component_Disk_Standard_Rotation_Curve_Gradient(self,self%halfMassRadius(),componentType,massType)
  velocityDisk                                =+[positionCartesianMidplane(2),-positionCartesianMidplane(1),0.0d0] &
       &                                        /radiusMidplane                                                     &
       &                                        *velocityCircular
  ! Compute epicyclic frequency.
  frequencyCircular               =velocityCircular              /     radiusMidplane
  frequencyCircularHalfMassRadius =velocityCircularHalfMassRadius/self%halfMassRadius()
  frequencyEpicyclic              =sqrt(velocityCircularSquaredGradient              /     radiusMidplane  +2.0d0*frequencyCircular              **2)
  frequencyEpicyclicHalfMassRadius=sqrt(velocityCircularSquaredGradientHalfMassRadius/self%halfMassRadius()+2.0d0*frequencyCircularHalfMassRadius**2)
  ! Get disk structural properties.
  density                     =+Node_Component_Disk_Standard_Density        (self,positionSpherical          ,componentTypeDisk,massTypeAll,weightByMass,weightIndexNull)
  densityMidPlane             =+Node_Component_Disk_Standard_Density        (self,positionSphericalMidplane  ,componentTypeDisk,massTypeAll,weightByMass,weightIndexNull)
  densitySurface              =+Node_Component_Disk_Standard_Surface_Density(self,positionCylindricalMidplane,componentTypeDisk,massTypeAll,weightByMass,weightIndexNull)
  densitySurfaceRadiusHalfMass=+Node_Component_Disk_Standard_Surface_Density(self,positionCylindricalHalfMass,componentTypeDisk,massTypeAll,weightByMass,weightIndexNull)
  if (density <= 0.0d0) return
  heightScale                 =+0.5d0           &
       &                       *densitySurface  &
       &                       /densityMidPlane
  ! Compute normalization of the radial velocity dispersion.
  velocityDispersionRadialHalfMass=+toomreQFactor                    &
       &                           *gravitationalConstantGalacticus  &
       &                           *densitySurfaceRadiusHalfMass     &
       &                           *toomreQRadiusHalfMass            &
       &                           /frequencyEpicyclicHalfMassRadius
  ! Find the velocity dispersion components of the disk.
  velocityDispersionRadial   =+velocityDispersionRadialHalfMass                &
       &                      *exp(-     radiusMidPlane  /self%radius()/2.0d0) &
       &                      /exp(-self%halfMassRadius()/self%radius()/2.0d0)
  velocityDispersionAzimuthal=+velocityDispersionRadial*frequencyEpicyclic/2.0d0/frequencyCircular
  velocityDispersionVertical =+sqrt(Pi*gravitationalConstantGalacticus*densitySurface*heightScale)
  velocityDispersionMaximum  =+maxval([velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical])
  velocityDispersionRadial   =+velocityDispersionRadial   /velocityDispersionMaximum
  velocityDispersionAzimuthal=+velocityDispersionAzimuthal/velocityDispersionMaximum
  velocityDispersionVertical =+velocityDispersionVertical /velocityDispersionMaximum
  if (any([velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical] <= 0.0d0)) return
  ! Find the relative velocity of the perturber and the disk.
  velocityRelative=(velocityCartesian-velocityDisk)/velocityDispersionMaximum
  ! Handle limiting case of large relative velocity.
  velocityRelativeMagnitude=sqrt(sum(velocityRelative**2))
  ! Initialize the velocity distribution.
  rotation=reshape(                                                                               &
       &            [                                                                             &
       &             +positionCartesianMidplane(1),-positionCartesianMidplane(2),+0.0d0         , &
       &             +positionCartesianMidplane(2),+positionCartesianMidplane(1),+0.0d0         , &
       &             +0.0d0                       ,+0.0d0                       ,+radiusMidplane  &
       &            ]                                                                             &
       &           /radiusMidplane                                                              , &
       &           [3,3]                                                                          &
       &          )
  coordinatesCartesian=velocityRelative
  if (.not.velocityDistributionInitialized) then
     velocityDistribution           =massDistributionGaussianEllipsoid(scaleLength=[1.0d0,1.0d0,1.0d0],rotation=rotation,mass=1.0d0,dimensionless=.true.)
     velocityDistributionInitialized=.true.
  end if
  call velocityDistribution%initialize(scaleLength=[velocityDispersionRadial,velocityDispersionAzimuthal,velocityDispersionVertical],rotation=rotation)
  ! Compute suppression factor due to satellite being an extended mass distribution. This is largely untested - it is meant to
  ! simply avoid extremely large accelerations for subhalo close to the disk plane when that subhalo is much more extended than
  ! the disk.
  if (radiusHalfMass > heightScale) then
     factorSuppressionExtendedMass=+heightScale    &
          &                        /radiusHalfMass
  else
     factorSuppressionExtendedMass=+1.0d0
  end if
  ! Evaluate the integral.
  Node_Component_Disk_Standard_Chandrasekhar_Integral=+density                                                             &
       &                                              *velocityDistribution         %acceleration(coordinatesCartesian)    &
       &                                              /velocityDispersionMaximum                                       **2 &
       &                                              *factorSuppressionExtendedMass
  return
end function Node_Component_Disk_Standard_Chandrasekhar_Integral

double precision function Node_Component_Disk_Standard_Density(self,positionSpherical,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the density at a given position for an standard disk.
  !!}
  use :: Coordinates                      , only : assignment(=)       , coordinateSpherical
  use :: Galactic_Structure_Options       , only : componentTypeAll    , componentTypeDisk           , massTypeAll            , massTypeBaryonic       , &
          &                                        massTypeGalactic    , massTypeGaseous             , massTypeStellar        , weightByLuminosity     , &
          &                                        weightByMass        , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: positionSpherical(3)
  type            (stellarLuminosities         ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  type            (coordinateSpherical         )                :: position

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Standard_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Determine mass/luminosity type.
  select case (weightBy%ID)
  case (weightByMass      %ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
        Node_Component_Disk_Standard_Density=self%massGas()+self%massStellar()
     case (massTypeGaseous%ID)
        Node_Component_Disk_Standard_Density=self%massGas()
     case (massTypeStellar%ID)
        Node_Component_Disk_Standard_Density=               self%massStellar()
     end select
  case (weightByLuminosity%ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Standard_Density=luminositiesDisk%luminosity(weightIndex)
     end select
  end select
  ! Skip further calculation if mass or radius is zero.
  if (Node_Component_Disk_Standard_Density > 0.0d0 .and. self%radius() > 0.0d0) then
     ! Compute the actual density.
      position=[positionSpherical(1)/self%radius(),positionSpherical(2),positionSpherical(3)]
      Node_Component_Disk_Standard_Density=+Node_Component_Disk_Standard_Density      &
           &                               /self                %radius (        )**3 &
           &                               *diskMassDistribution%density(position)
   end if
   return
end function Node_Component_Disk_Standard_Density

double precision function Node_Component_Disk_Standard_Density_Spherical_Average(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the density at a given position for an standard disk.
  !!}
  use :: Galactic_Structure_Options       , only : componentTypeAll    , componentTypeDisk           , massTypeAll            , massTypeBaryonic       , &
          &                                        massTypeGalactic    , massTypeGaseous             , massTypeStellar        , weightByLuminosity     , &
          &                                        weightByMass        , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: radius
  type            (stellarLuminosities         ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Standard_Density_Spherical_Average=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Determine mass/luminosity type.
  select case (weightBy%ID)
  case (weightByMass      %ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
        Node_Component_Disk_Standard_Density_Spherical_Average=self%massGas()+self%massStellar()
     case (massTypeGaseous%ID)
        Node_Component_Disk_Standard_Density_Spherical_Average=self%massGas()
     case (massTypeStellar%ID)
        Node_Component_Disk_Standard_Density_Spherical_Average=               self%massStellar()
     end select
  case (weightByLuminosity%ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Standard_Density_Spherical_Average=luminositiesDisk%luminosity(weightIndex)
     end select
  end select
  ! Skip further calculation if mass or radius is zero.
  if (Node_Component_Disk_Standard_Density_Spherical_Average > 0.0d0 .and. self%radius() > 0.0d0) then
     ! Compute the actual density.
      Node_Component_Disk_Standard_Density_Spherical_Average=+Node_Component_Disk_Standard_Density_Spherical_Average                &
           &                                                 /self                %radius                 (                    )**3 &
           &                                                 *diskMassDistribution%densitySphericalAverage(radius/self%radius())
   end if
   return
 end function Node_Component_Disk_Standard_Density_Spherical_Average

double precision function Node_Component_Disk_Standard_Potential(self,radius,componentType,massType,status)
  !!{
  Compute the gravitational potential due to an standard disk.
  !!}
  use :: Coordinates                      , only : assignment(=)                  , coordinateCylindrical
  use :: Galactic_Structure_Options       , only : componentTypeAll               , componentTypeDisk           , radiusLarge            , weightByMass                     , &
          &                                        weightIndexNull                , enumerationComponentTypeType, enumerationMassTypeType, enumerationStructureErrorCodeType
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  use :: Numerical_Constants_Astronomical , only : gravitationalConstantGalacticus
  implicit none
  class           (nodeComponentDiskStandard        ), intent(inout)           :: self
  type            (enumerationComponentTypeType     ), intent(in   )           :: componentType
  type            (enumerationMassTypeType          ), intent(in   )           :: massType
  double precision                                   , intent(in   )           :: radius
  type            (enumerationStructureErrorCodeType), intent(inout), optional :: status
  double precision                                                             :: componentMass
  type            (coordinateCylindrical            )                          :: position
  !$GLC attributes unused :: status

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Standard_Potential=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Avoid an arithmetic exception at radius zero.
  if (radius <= 0.0d0) return
  ! Get the relevant mass of the disk.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull)
  if (componentMass <= 0.0d0) return
  ! Check for zero-sized disk.
  if (self%radius() <= 0.0d0) then
     ! Treat as a point mass.
     Node_Component_Disk_Standard_Potential=-gravitationalConstantGalacticus          &
          &                                 *componentMass                            &
          &                                 /radius
  else
     ! Compute the potential.
     position=[radius/self%radius(),0.0d0,0.0d0]
     Node_Component_Disk_Standard_Potential=+gravitationalConstantGalacticus          &
          &                                 *componentMass                            &
          &                                 /self%radius()                            &
          &                                 *diskMassDistribution%potential(position)
  end if
  return
end function Node_Component_Disk_Standard_Potential

double precision function Node_Component_Disk_Standard_Rotation_Curve(self,radius,componentType,massType)
  !!{
  Computes the rotation curve at a given radius for an standard disk.
  !!}
  use :: Galactic_Structure_Options       , only : radiusLarge                    , weightByMass, weightIndexNull, enumerationComponentTypeType, &
       &                                           enumerationMassTypeType
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  use :: Numerical_Constants_Astronomical , only : gravitationalConstantGalacticus
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  double precision                              , intent(in   ) :: radius
  double precision                                              :: componentMass, radiusDisk, fractionalRadius

  ! Set to zero by default.
  Node_Component_Disk_Standard_Rotation_Curve=0.0d0
  ! Get the mass of the disk.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull)
  if (componentMass <= 0.0d0) return
  ! Compute the actual velocity.
  radiusDisk=self%radius()
  if (radiusDisk > 0.0d0) then
     fractionalRadius=radius/radiusDisk
     Node_Component_Disk_Standard_Rotation_Curve=+sqrt(                                                &
          &                                            +gravitationalConstantGalacticus                &
          &                                            *componentMass                                  &
          &                                            /radiusDisk                                     &
          &                                           )                                                &
          &                                      *diskMassDistribution%rotationCurve(fractionalRadius)
  end if
  return
end function Node_Component_Disk_Standard_Rotation_Curve

double precision function Node_Component_Disk_Standard_Rotation_Curve_Gradient(self,radius,componentType,massType)
  !!{
  Computes the rotation curve gradient for an standard disk.
  !!}
  use :: Galactic_Structure_Options       , only : radiusLarge                    , weightByMass, weightIndexNull, enumerationComponentTypeType, &
       &                                           enumerationMassTypeType
  use :: Node_Component_Disk_Standard_Data, only : diskMassDistribution
  use :: Numerical_Constants_Astronomical , only : gravitationalConstantGalacticus
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  double precision                              , intent(in   ) :: radius
  double precision                                              :: radiusDisk   , fractionalRadius, &
       &                                                           componentMass

  ! Set to zero by default.
  Node_Component_Disk_Standard_Rotation_Curve_Gradient=0.0d0
  ! Return if radius is zero.
  if (radius <= 0.0d0) return
  ! Get the mass of the disk.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull)
  if (componentMass <= 0.0d0) return
  ! Compute the rotation curve gradient.
  radiusDisk=self%radius()
  if (radiusDisk > 0.0d0) then
     fractionalRadius=radius/radiusDisk
     Node_Component_Disk_Standard_Rotation_Curve_Gradient=+gravitationalConstantGalacticus                              &
          &                                               *componentMass                                                &
          &                                               /radiusDisk**2                                                &
          &                                               *diskMassDistribution%rotationCurveGradient(fractionalRadius)
  end if
  return
end function Node_Component_Disk_Standard_Rotation_Curve_Gradient

double precision function Node_Component_Disk_Standard_Surface_Density(self,positionCylindrical,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the surface density at a given position for an standard disk.
  !!}
  use :: Coordinates                      , only : coordinateCylindrical
  use :: Galactic_Structure_Options       , only : componentTypeAll                    , componentTypeDisk           , massTypeAll                       , massTypeBaryonic            , &
          &                                        massTypeGalactic                    , massTypeGaseous             , massTypeStellar                   , weightByLuminosity          , &
          &                                        weightByMass                        , enumerationComponentTypeType, enumerationMassTypeType           , enumerationWeightByType
  use :: Node_Component_Disk_Standard_Data, only : Node_Component_Disk_Standard_Reset  , diskMassDistribution        , lastUniqueID                      , radiusScaleDisk             , &
          &                                        radiusScaleDiskComputed             , surfaceDensityCentralGas    , surfaceDensityCentralGasComputed  , surfaceDensityCentralStellar, &
          &                                        surfaceDensityCentralStellarComputed, surfaceDensityCentralTotal  , surfaceDensityCentralTotalComputed
  use :: Numerical_Constants_Math         , only : Pi
  implicit none
  class           (nodeComponentDiskStandard   ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: positionCylindrical(3)
  type            (treeNode                    ), pointer       :: selfNode
  type            (stellarLuminosities         ), save          :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  type            (coordinateCylindrical       )                :: position

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Standard_Surface_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Check whether this is a new node.
  selfNode => self%host()
  if (selfNode%uniqueID() /= lastUniqueID) call Node_Component_Disk_Standard_Reset(selfNode%uniqueID())
  ! Determine disk radius.
  if (.not.radiusScaleDiskComputed) then
     radiusScaleDisk        =self%radius()
     radiusScaleDiskComputed=.true.
  end if
  ! Return zero if the disk has unphysical size.
  if (radiusScaleDisk <= 0.0d0) then
     Node_Component_Disk_Standard_Surface_Density=0.0d0
     return
  end if
  ! Determine mass type.
  select case (weightBy%ID)
  case (weightByMass      %ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
        if (.not.surfaceDensityCentralTotalComputed  ) then
           surfaceDensityCentralTotal          =(self%massGas()+self%massStellar())/radiusScaleDisk**2
           surfaceDensityCentralTotalComputed  =.true.
        end if
        Node_Component_Disk_Standard_Surface_Density=surfaceDensityCentralTotal
     case (massTypeGaseous%ID)
        if (.not.surfaceDensityCentralGasComputed    ) then
           surfaceDensityCentralGas            = self%massGas()                    /radiusScaleDisk**2
           surfaceDensityCentralGasComputed    =.true.
        end if
        Node_Component_Disk_Standard_Surface_Density=surfaceDensityCentralGas
     case (massTypeStellar%ID)
        if (.not.surfaceDensityCentralStellarComputed) then
           surfaceDensityCentralStellar        =                self%massStellar() /radiusScaleDisk**2
           surfaceDensityCentralStellarComputed=.true.
        end if
        Node_Component_Disk_Standard_Surface_Density=surfaceDensityCentralStellar
     end select
  case (weightByLuminosity%ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Standard_Surface_Density=luminositiesDisk%luminosity(weightIndex)/2.0d0/Pi/radiusScaleDisk**2
     end select
  end select
  ! Return if no surface density.
  if (Node_Component_Disk_Standard_Surface_Density <= 0.0d0) return
  ! Compute the surface density.
  call position%rSet(positionCylindrical(1)/radiusScaleDisk)
  Node_Component_Disk_Standard_Surface_Density=         &
       & +Node_Component_Disk_Standard_Surface_Density  &
       & *diskMassDistribution%surfaceDensity(position)
  return
end function Node_Component_Disk_Standard_Surface_Density

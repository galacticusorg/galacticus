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
Contains custom functions for the standard spheroid component.
!!}

function Node_Component_Spheroid_Standard_Mass_Distribution(self) result(massDistribution_)
  !!{
  Return the mass distribution for the standard spheroid component.
  !!}
  use :: Mass_Distributions                   , only : massDistributionClass      , massDistributionSphericalScaler, massDistributionComposite, massDistributionList, &
        &                                              massDistributionSpherical
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_   , massDistributionGas_
  use :: Galactic_Structure_Options           , only : componentTypeSpheroid      , massTypeStellar                , massTypeGaseous
  implicit none
  class           (massDistributionClass          ), pointer       :: massDistribution_
  class           (nodeComponentSpheroidStandard  ), intent(inout) :: self
  type            (massDistributionSphericalScaler), pointer       :: massDistributionStellar   , massDistributionGas
  type            (massDistributionComposite      ), pointer       :: massDistributionTotal
  type            (massDistributionList           ), pointer       :: massDistributionComponents
  double precision                                                 :: massStellar               , massGas            , &
       &                                                              radiusScale

  ! Get properties of the mass distribution and ensure they are physical.
  massStellar=max(0.0d0,self%massStellar())
  massGas    =max(0.0d0,self%massGas    ())
  radiusScale=          self%radius     ()
  if (radiusScale <= 0.0d0) then
     ! Spjeroid has non-positive size. Set an arbitrary, positive size, but also set masses to zero. This results in an empty mass
     ! distribution.
     radiusScale=1.0d0
     massStellar=0.0d0
     massGas    =0.0d0
  end if
  ! Wrap the dimensionless mass distribution inside scaler classes to allow us to re-scale it to any spheroid system, and then composite those.
  allocate(massDistributionStellar        )
  allocate(massDistributionGas            )
  allocate(massDistributionTotal          )
  allocate(massDistributionComponents     )
  allocate(massDistributionComponents%next)
  massDistributionComponents     %massDistribution_ => massDistributionStellar
  massDistributionComponents%next%massDistribution_ => massDistributionGas
  select type (massDistributionStellar_)
  class is (massDistributionSpherical)
     !![
     <referenceConstruct object="massDistributionStellar" constructor="massDistributionSphericalScaler(factorScalingLength=radiusScale,factorScalingMass=massStellar,massDistribution_=massDistributionStellar_  )"/>
     !!]
  end select
  select type (massDistributionGas_    )
  class is (massDistributionSpherical)
     !![
     <referenceConstruct object="massDistributionGas"     constructor="massDistributionSphericalScaler(factorScalingLength=radiusScale,factorScalingMass=massGas    ,massDistribution_=massDistributionGas_      )"/>
     !!]
  end select
  !![
  <referenceConstruct    object="massDistributionTotal"   constructor="massDistributionComposite      (                                                                                massDistributionComponents)"/>
  <objectDestructor name="massDistributionStellar"/>
  <objectDestructor name="massDistributionGas"    />
  !!]
  nullify(massDistributionComponents)
  ! Return a pointer to the spheroid mass distribution.
  massDistribution_ => massDistributionTotal
  return
end function Node_Component_Spheroid_Standard_Mass_Distribution

double precision function Node_Component_Spheroid_Standard_Half_Mass_Radius(self) result(radiusHalfMass)
  !!{
  Return the half-mass radius of the standard spheroid.
  !!}
  use :: Error                                , only : Error_Report
  use :: Mass_Distributions                   , only : massDistributionSpherical
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionStellar_
  implicit none
  class(nodeComponentSpheroidStandard), intent(inout) :: self


  select type (massDistributionStellar_)
  class is (massDistributionSpherical)
     radiusHalfMass=+massDistributionStellar_%radiusHalfMass() &
          &         *self                    %radius        ()
  class default
     radiusHalfMass=0.0d0
     call Error_Report('spheroid mass distribution is not spherically-symmetric'//{introspection:location})
  end select
  return
end function Node_Component_Spheroid_Standard_Half_Mass_Radius

double precision function Node_Component_Spheroid_Standard_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex) result(massEnclosed)
  !!{
  Computes the mass within a given radius for an standard spheroid.
  !!}
  use :: Galactic_Structure_Options, only : weightByLuminosity   , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  type            (enumerationWeightByType      ), intent(in   ) :: weightBy
  integer                                        , intent(in   ) :: weightIndex
  double precision                               , intent(in   ) :: radius
  class           (massDistributionClass        ), pointer       :: massDistribution_
  type            (stellarLuminosities          ), save          :: luminositiesSpheroid
  !$omp threadprivate(luminositiesSpheroid)
  
  massDistribution_ => self             %massDistribution    (                             )
  massEnclosed      =  massDistribution_%massEnclosedBySphere(radius,componentType,massType)
  if (massEnclosed > 0.0d0 .and. weightBy     == weightByLuminosity) then
     luminositiesSpheroid= self%luminositiesStellar()
     massEnclosed        =+massEnclosed                                            &
          &               *luminositiesSpheroid%luminosity(weightIndex           ) &
          &               /massDistribution_   %massTotal (componentType,massType)
  end if
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Spheroid_Standard_Enclosed_Mass

function Node_Component_Spheroid_Standard_Acceleration(self,positionCartesian,componentType,massType) result(acceleration)
  !!{
  Computes the gravitational acceleration at a given position for a standard spheroid.
  !!}
  use :: Coordinates               , only : assignment(=)               , coordinateCartesian
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  double precision                                              , dimension(3) :: acceleration
  class           (nodeComponentSpheroidStandard), intent(inout)               :: self
  double precision                               , intent(in   ), dimension(3) :: positionCartesian
  type            (enumerationComponentTypeType ), intent(in   )               :: componentType
  type            (enumerationMassTypeType      ), intent(in   )               :: massType
  class           (massDistributionClass        ), pointer                     :: massDistribution_
  type            (coordinateCartesian          )                              :: position

  position          =  positionCartesian
  massDistribution_ => self             %massDistribution(                               )
  acceleration      =  massDistribution_%acceleration    (position,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Spheroid_Standard_Acceleration

function Node_Component_Spheroid_Standard_Tidal_Tensor(self,positionCartesian,componentType,massType) result(tidalTensor)
  !!{
  Computes the gravitational acceleration at a given position for a standard spheroid.
  !!}
  use :: Coordinates               , only : assignment(=)                 , coordinateCartesian
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType  , enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Tensors                   , only : tensorRank2Dimension3Symmetric
  implicit none
  type            (tensorRank2Dimension3Symmetric)                              :: tidalTensor
  class           (nodeComponentSpheroidStandard ), intent(inout)               :: self
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
end function Node_Component_Spheroid_Standard_Tidal_Tensor

function Node_Component_Spheroid_Standard_Chandrasekhar_Integral(self,nodeSatellite,positionCartesian,velocityCartesian,componentType,massType) result(integralChandrasekhar)
  !!{
  Computes the Chandrasekhar integral at a given position for a standard spheroid.
  !!}
  use :: Coordinates               , only : coordinateCartesian  , assignment(=)
  use :: Galactic_Structure_Options, only : weightByMass         , weightIndexNull, enumerationComponentTypeType, enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  use :: Numerical_Constants_Math  , only : Pi
  implicit none
  double precision                                              , dimension(3) :: integralChandrasekhar
  class           (nodeComponentSpheroidStandard), intent(inout)               :: self
  type            (treeNode                     ), intent(inout)               :: nodeSatellite
  double precision                               , intent(in   ), dimension(3) :: positionCartesian           , velocityCartesian
  type            (enumerationComponentTypeType ), intent(in   )               :: componentType
  type            (enumerationMassTypeType      ), intent(in   )               :: massType
  double precision                                              , dimension(3) :: positionSpherical
  class           (massDistributionClass        ), pointer                     :: massDistribution_
  type            (coordinateCartesian          )                              :: position
  double precision                               , parameter                   :: XvMaximum            =10.0d0
  double precision                                                             :: radius                      , velocity          , &
       &                                                                          density                     , xV                , &
       &                                                                          velocityRotation            , velocityDispersion
  !$GLC attributes unused :: nodeSatellite
  
  position              =  positionCartesian
  massDistribution_     => self             %massDistribution()
  ! Compute the integral.
  integralChandrasekhar =  +0.0d0
  radius                =  +sqrt(sum(positionCartesian**2))
  velocity              =  +sqrt(sum(velocityCartesian**2))
  if (velocity <= 0.0d0) return
  positionSpherical     =  [radius,0.0d0,0.0d0]
  velocityRotation      =  massDistribution_%rotationCurve(radius  ,componentType,massType)
  density               =  massDistribution_%density      (position,componentType,massType)
  if (density  <= 0.0d0) return
  velocityDispersion    = +velocityRotation     &
       &                  /sqrt(2.0d0)
  xV                    = +velocity             &
       &                  /velocityDispersion   &
       &                  /sqrt(2.0d0)
  integralChandrasekhar = -density              &
       &                  *velocityCartesian    &
       &                  /velocity         **3
  if (Xv <= XvMaximum)                                &
       & integralChandrasekhar=+integralChandrasekhar &
       &                        *(                    &
       &                          +erf ( xV   )       &
       &                          -2.0d0              &
       &                          *      xV           &
       &                          *exp (-xV**2)       &
       &                          /sqrt( Pi   )       &
       &                        )
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Spheroid_Standard_Chandrasekhar_Integral

double precision function Node_Component_Spheroid_Standard_Rotation_Curve(self,radius,componentType,massType) result(rotationCurve)
  !!{
  Computes the rotation curve at a given radius for an standard spheroid.
  !!}
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  double precision                               , intent(in   ) :: radius
  class           (massDistributionClass        ), pointer       :: massDistribution_

  massDistribution_ => self             %massDistribution(                             )
  rotationCurve     =  massDistribution_%rotationCurve   (radius,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Spheroid_Standard_Rotation_Curve

double precision function Node_Component_Spheroid_Standard_Rotation_Curve_Gradient(self,radius,componentType,massType) result(rotationCurveGradient)
  !!{
  Computes the rotation curve gradient for an standard spheroid.
  !!}
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  double precision                               , intent(in   ) :: radius
  class           (massDistributionClass        ), pointer       :: massDistribution_

  massDistribution_     => self             %massDistribution     (                             )
  rotationCurveGradient =  massDistribution_%rotationCurveGradient(radius,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Spheroid_Standard_Rotation_Curve_Gradient

double precision function Node_Component_Spheroid_Standard_Density(self,positionSpherical,componentType,massType,weightBy,weightIndex) result(density)
  !!{
  Computes the density at a given position for an standard spheroid.
  !!}
  use :: Coordinates               , only : assignment(=)                 , coordinateSpherical
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType  , enumerationMassTypeType, enumerationWeightByType, weightByLuminosity
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  type            (enumerationWeightByType      ), intent(in   ) :: weightBy
  integer                                        , intent(in   ) :: weightIndex
  double precision                               , intent(in   ) :: positionSpherical(3)
  class           (massDistributionClass        ), pointer       :: massDistribution_
  type            (stellarLuminosities          ), save          :: luminositiesSpheroid
  !$omp threadprivate(luminositiesSpheroid)
  type            (coordinateSpherical          )                :: position


  position          =  positionSpherical
  massDistribution_ => self             %massDistribution(                               )
  density           =  massDistribution_%density         (position,componentType,massType)
  if (density > 0.0d0 .and. weightBy == weightByLuminosity) then
     luminositiesSpheroid= self                %luminositiesStellar(                      )
     density             =+density                                                          &
          &               *luminositiesSpheroid%luminosity         (weightIndex           ) &
          &               /massDistribution_   %massTotal          (componentType,massType)
  end if  
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Spheroid_Standard_Density

double precision function Node_Component_Spheroid_Standard_Density_Spherical_Average(self,radius,componentType,massType,weightBy,weightIndex) result(density)
  !!{
  Computes the density at a given position for an standard spheroid.
  !!}
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType, weightByLuminosity
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  type            (enumerationWeightByType      ), intent(in   ) :: weightBy
  integer                                        , intent(in   ) :: weightIndex
  double precision                               , intent(in   ) :: radius
  type            (stellarLuminosities          ), save          :: luminositiesSpheroid
  !$omp threadprivate(luminositiesSpheroid)
  class           (massDistributionClass        ), pointer       :: massDistribution_

  massDistribution_ => self             %massDistribution       (                             )
  density           =  massDistribution_%densitySphericalAverage(radius,componentType,massType)
  if (density > 0.0d0 .and. weightBy == weightByLuminosity) then
     luminositiesSpheroid= self                %luminositiesStellar(                      )
     density             =+density                                                          &
          &               *luminositiesSpheroid%luminosity         (weightIndex           ) &
          &               /massDistribution_   %massTotal          (componentType,massType)
  end if
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
 end function Node_Component_Spheroid_Standard_Density_Spherical_Average

double precision function Node_Component_Spheroid_Standard_Potential(self,radius,componentType,massType,status) result(potential)
  !!{
  Compute the gravitational potential due to an standard spheroid.
  !!}
  use :: Coordinates               , only : assignment(=)               , coordinateSpherical
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationStructureErrorCodeType
  use :: Mass_Distributions        , only : massDistributionClass
  implicit none
  class           (nodeComponentSpheroidStandard    ), intent(inout)           :: self
  type            (enumerationComponentTypeType     ), intent(in   )           :: componentType
  type            (enumerationMassTypeType          ), intent(in   )           :: massType
  double precision                                   , intent(in   )           :: radius
  type            (enumerationStructureErrorCodeType), intent(inout), optional :: status
  class           (massDistributionClass            ), pointer                 :: massDistribution_
  type            (coordinateSpherical              )                          :: position
  !$GLC attributes unused :: status

  position          =  [radius,0.0d0,0.0d0]
  massDistribution_ => self             %massDistribution(                               )
  potential         =  massDistribution_%potential       (position,componentType,massType)
  !![
  <objectDestructor name="massDistribution_"/>
  !!]
  return
end function Node_Component_Spheroid_Standard_Potential

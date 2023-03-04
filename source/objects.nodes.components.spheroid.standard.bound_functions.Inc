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

double precision function Node_Component_Spheroid_Standard_Half_Mass_Radius(self)
  !!{
  Return the half-mass radius of the standard spheroid.
  !!}
  use :: Error                                , only : Error_Report
  use :: Mass_Distributions                   , only : massDistributionSpherical
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionSpheroid
  implicit none
  class(nodeComponentSpheroidStandard), intent(inout) :: self

  select type (massDistributionSpheroid)
  class is (massDistributionSpherical)
     Node_Component_Spheroid_Standard_Half_Mass_Radius=self%radius()*massDistributionSpheroid%radiusHalfMass()
  class default
     Node_Component_Spheroid_Standard_Half_Mass_Radius=0.0d0
     call Error_Report('mass distribution should be spherical'//{introspection:location})
  end select
  return
end function Node_Component_Spheroid_Standard_Half_Mass_Radius

double precision function Node_Component_Spheroid_Standard_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the mass within a given radius for an standard spheroid.
  !!}
  use :: Galactic_Structure_Options           , only : componentTypeAll        , componentTypeSpheroid, massTypeAll                 , massTypeBaryonic       , &
          &                                            massTypeGalactic        , massTypeGaseous      , massTypeStellar             , radiusLarge            , &
          &                                            weightByLuminosity      , weightByMass         , enumerationComponentTypeType, enumerationMassTypeType, &
          &                                            enumerationWeightByType
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionSpheroid
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  type            (enumerationWeightByType      ), intent(in   ) :: weightBy
  integer                                        , intent(in   ) :: weightIndex
  double precision                               , intent(in   ) :: radius
  double precision                                               :: fractionalRadius    , radiusSpheroid
  type            (stellarLuminosities          ), save          :: luminositiesSpheroid
  !$omp threadprivate(luminositiesSpheroid)

  ! Return if the spheroid component is not selected.
  Node_Component_Spheroid_Standard_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
  ! Get the total mass.
  select case (weightBy%ID)
  case (weightByMass      %ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
        Node_Component_Spheroid_Standard_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous%ID)
        Node_Component_Spheroid_Standard_Enclosed_Mass=self%massGas()
     case (massTypeStellar%ID)
        Node_Component_Spheroid_Standard_Enclosed_Mass=               self%massStellar()
     end select
  case (weightByLuminosity%ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
        luminositiesSpheroid=self%luminositiesStellar()
        Node_Component_Spheroid_Standard_Enclosed_Mass       =luminositiesSpheroid%luminosity(weightIndex)
     end select
  end select
  ! Return if total mass was requested.
  if (radius                                         >= radiusLarge) return
  ! Return if mass is zero.
  if (Node_Component_Spheroid_Standard_Enclosed_Mass <=       0.0d0) return
  ! Compute actual mass.
  radiusSpheroid=self%radius()
  if (radiusSpheroid > 0.0d0) then
     fractionalRadius=radius/radiusSpheroid
     Node_Component_Spheroid_Standard_Enclosed_Mass   =Node_Component_Spheroid_Standard_Enclosed_Mass*massDistributionSpheroid%massEnclosedBySphere(fractionalRadius)
  end if
  return
end function Node_Component_Spheroid_Standard_Enclosed_Mass

function Node_Component_Spheroid_Standard_Acceleration(self,positionCartesian,componentType,massType)
  !!{
  Computes the gravitational acceleration at a given position for a standard spheroid.
  !!}
  use :: Galactic_Structure_Options      , only : componentTypeAll            , componentTypeSpheroid  , weightByMass                   , weightIndexNull, &
       &                                          enumerationComponentTypeType, enumerationMassTypeType
  use :: Numerical_Constants_Astronomical, only : gigaYear                    , megaParsec             , gravitationalConstantGalacticus
  use :: Numerical_Constants_Prefixes    , only : kilo
  implicit none
  double precision                                              , dimension(3) :: Node_Component_Spheroid_Standard_Acceleration
  class           (nodeComponentSpheroidStandard), intent(inout)               :: self
  double precision                               , intent(in   ), dimension(3) :: positionCartesian
  type            (enumerationComponentTypeType ), intent(in   )               :: componentType
  type            (enumerationMassTypeType      ), intent(in   )               :: massType
  double precision                                                             :: radius                                       , massEnclosed
  
  ! Return if the spheroid component is not selected.
  Node_Component_Spheroid_Standard_Acceleration=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
  ! Compute the acceleration.
  radius                                       =sqrt(sum(positionCartesian**2))
  massEnclosed                                 =Node_Component_Spheroid_Standard_Enclosed_Mass(self,radius,componentType,massType,weightByMass,weightIndexNull)
  Node_Component_Spheroid_Standard_Acceleration=-kilo                            &
         &                                      *gigaYear                        &
         &                                      /megaParsec                      &
         &                                      *gravitationalConstantGalacticus &
         &                                      *massEnclosed                    &
         &                                      *positionCartesian               &
         &                                      /radius**3
  return
end function Node_Component_Spheroid_Standard_Acceleration

function Node_Component_Spheroid_Standard_Tidal_Tensor(self,positionCartesian,componentType,massType)
  !!{
  Computes the gravitational tidal tensor at a given position for a standard spheroid.
  !!}
  use :: Coordinates                          , only : assignment(=)                  , coordinateCartesian
  use :: Galactic_Structure_Options           , only : componentTypeAll               , componentTypeSpheroid       , weightByMass           , weightIndexNull, &
       &                                               radiusLarge                    , enumerationComponentTypeType, enumerationMassTypeType
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionSpheroid
  use :: Numerical_Constants_Astronomical     , only : gravitationalConstantGalacticus
  use :: Tensors                              , only : tensorRank2Dimension3Symmetric , tensorNullR2D3Sym    , operator(*)
  implicit none
  type            (tensorRank2Dimension3Symmetric)                              :: Node_Component_Spheroid_Standard_Tidal_Tensor
  class           (nodeComponentSpheroidStandard ), intent(inout)               :: self
  double precision                                , intent(in   ), dimension(3) :: positionCartesian
  type            (enumerationComponentTypeType  ), intent(in   )               :: componentType
  type            (enumerationMassTypeType       ), intent(in   )               :: massType
  double precision                                                              :: massTotal
  type            (coordinateCartesian           )                              :: positionScaleFree

  ! Return if the spheroid component is not selected.
  Node_Component_Spheroid_Standard_Tidal_Tensor=tensorNullR2D3Sym
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
  ! Find the total mass.
  massTotal=Node_Component_Spheroid_Standard_Enclosed_Mass(self,radiusLarge,componentType,massType,weightByMass,weightIndexNull)
  if (massTotal <= 0.0d0 .or. self%radius() <= 0.0d0) return
  ! Compute the tidal tensor.
  positionScaleFree=positionCartesian/self%radius()
  Node_Component_Spheroid_Standard_Tidal_Tensor=+gravitationalConstantGalacticus                         &
       &                                        *massTotal                                               &
       &                                        *massDistributionSpheroid%tidalTensor(positionScaleFree) &
       &                                        /self%radius()**3
  return
end function Node_Component_Spheroid_Standard_Tidal_Tensor

function Node_Component_Spheroid_Standard_Chandrasekhar_Integral(self,nodeSatellite,positionCartesian,velocityCartesian,componentType,massType)
  !!{
  Computes the gravitational acceleration at a given position for a standard spheroid.
  !!}
  use :: Galactic_Structure_Options, only : componentTypeAll            , componentTypeSpheroid  , weightByMass, weightIndexNull, &
       &                                    enumerationComponentTypeType, enumerationMassTypeType
  use :: Numerical_Constants_Math  , only : Pi
  implicit none
  double precision                                              , dimension(3) :: Node_Component_Spheroid_Standard_Chandrasekhar_Integral
  class           (nodeComponentSpheroidStandard), intent(inout)               :: self
  type            (treeNode                     ), intent(inout)               :: nodeSatellite
  double precision                               , intent(in   ), dimension(3) :: positionCartesian                                             , velocityCartesian
  type            (enumerationComponentTypeType ), intent(in   )               :: componentType
  type            (enumerationMassTypeType      ), intent(in   )               :: massType
  double precision                                              , dimension(3) :: positionSpherical
  double precision                               , parameter                   :: XvMaximum                                              =10.0d0
  double precision                                                             :: radius                                                        , velocity          , &
       &                                                                          density                                                       , xV                , &
       &                                                                          velocityRotation                                              , velocityDispersion
  !$GLC attributes unused :: nodeSatellite
  
  ! Return if the spheroid component is not selected.
  Node_Component_Spheroid_Standard_Chandrasekhar_Integral=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
  ! Compute the integral.
  radius                                                  =  sqrt(sum(positionCartesian**2))
  velocity                                                =  sqrt(sum(velocityCartesian**2))
  if (velocity <= 0.0d0) return
  positionSpherical                                       =  [radius,0.0d0,0.0d0]
  velocityRotation                                        =  Node_Component_Spheroid_Standard_Rotation_Curve(self,radius           ,componentType,massType                             )
  density                                                 =  Node_Component_Spheroid_Standard_Density       (self,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
  if (density  <= 0.0d0) return
  velocityDispersion                                      = +velocityRotation   &
       &                                                    /sqrt(2.0d0)
  xV                                                      = +velocity           &
       &                                                    /velocityDispersion &
       &                                                    /sqrt(2.0d0)
  Node_Component_Spheroid_Standard_Chandrasekhar_Integral = -density              &
       &                                                    *velocityCartesian    &
       &                                                    /velocity         **3
  if (Xv <= XvMaximum)                                                                                                    &
       & Node_Component_Spheroid_Standard_Chandrasekhar_Integral=+Node_Component_Spheroid_Standard_Chandrasekhar_Integral &
       &                                                         *(                                                       &
       &                                                           +erf ( xV   )                                          &
       &                                                           -2.0d0                                                 &
       &                                                           *      xV                                              &
       &                                                           *exp (-xV**2)                                          &
       &                                                           /sqrt( Pi   )                                          &
       &                                                         )
  return
end function Node_Component_Spheroid_Standard_Chandrasekhar_Integral

double precision function Node_Component_Spheroid_Standard_Rotation_Curve(self,radius,componentType,massType)
  !!{
  Computes the rotation curve at a given radius for a standard spheroid.
  !!}
  use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull, enumerationComponentTypeType, enumerationMassTypeType
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  double precision                               , intent(in   ) :: radius
  double precision                                               :: componentMass

  ! Set to zero by default.
  Node_Component_Spheroid_Standard_Rotation_Curve=0.0d0

  ! Return immediately for non-positive radius.
  if (radius <= 0.0d0) return

  ! Compute if a spheroid is present.
  componentMass=self%enclosedMass(radius,componentType,massType,weightByMass,weightIndexNull)
  if (componentMass > 0.0d0) Node_Component_Spheroid_Standard_Rotation_Curve=sqrt(gravitationalConstantGalacticus&
       &*componentMass)/sqrt(radius)
  return
end function Node_Component_Spheroid_Standard_Rotation_Curve

double precision function Node_Component_Spheroid_Standard_Rotation_Curve_Gradient(self,radius,componentType,massType)
  !!{
  Computes the rotation curve gradient for the standard spheroid.
  !!}
  use :: Galactic_Structure_Options      , only : weightByMass                   , weightIndexNull, enumerationComponentTypeType, enumerationMassTypeType
  use :: Numerical_Constants_Math        , only : Pi
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  double precision                               , intent(in   ) :: radius
  double precision                                               :: componentDensity, componentMass, positionSpherical(3)

  ! Set to zero by default.
  Node_Component_Spheroid_Standard_Rotation_Curve_Gradient=0.0d0

  ! Return immediately for non-positive radius.
  if (radius <= 0.0d0) return

  ! Compute if a spheroid is present.
  positionSpherical= [radius,0.0d0,0.0d0]
  componentMass    =self%enclosedMass(radius           ,componentType,massType,weightByMass,weightIndexNull)
  componentDensity =self%density     (positionSpherical,componentType,massType,weightByMass,weightIndexNull)
  if (componentMass == 0.0d0 .or. componentDensity == 0.0d0) return
  Node_Component_Spheroid_Standard_Rotation_Curve_Gradient=  &
       &                  gravitationalConstantGalacticus    &
       &                 *(                                  &
       &                   -componentMass/radius**2          &
       &                   +4.0d0*Pi*radius*componentDensity &
       &                  )
  return
end function Node_Component_Spheroid_Standard_Rotation_Curve_Gradient

double precision function Node_Component_Spheroid_Standard_Density(self,positionSpherical,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the density at a given position for an standard spheroid.
  !!}
  use :: Coordinates                          , only : assignment(=)           , coordinateSpherical
  use :: Galactic_Structure_Options           , only : componentTypeAll        , componentTypeSpheroid       , massTypeAll            , massTypeBaryonic       , &
          &                                            massTypeGalactic        , massTypeGaseous             , massTypeStellar        , weightByLuminosity     , &
          &                                            weightByMass            , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionSpheroid
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  type            (enumerationWeightByType      ), intent(in   ) :: weightBy
  integer                                        , intent(in   ) :: weightIndex
  double precision                               , intent(in   ) :: positionSpherical(3)
  double precision                               , parameter     :: radiusHuge          =1.0d+100
  type            (coordinateSpherical          )                :: position
  type            (stellarLuminosities          ), save          :: luminositiesSpheroid
  !$omp threadprivate(luminositiesSpheroid)

  Node_Component_Spheroid_Standard_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return

  ! Get the spheroid component and check that it is of the standard class.
  select type (self)
     class is (nodeComponentSpheroidStandard)

     if (self%radius() <= 0.0d0 .or. self%radius() > radiusHuge) return
     select case (weightBy%ID)
     case (weightByMass      %ID)
        select case (massType%ID)
        case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
           Node_Component_Spheroid_Standard_Density=self%massGas()+self%massStellar()
        case (massTypeGaseous%ID)
           Node_Component_Spheroid_Standard_Density=self%massGas()
        case (massTypeStellar%ID)
           Node_Component_Spheroid_Standard_Density=               self%massStellar()
        end select
     case (weightByLuminosity%ID)
        select case (massType%ID)
        case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
           luminositiesSpheroid=self%luminositiesStellar()
           Node_Component_Spheroid_Standard_Density=luminositiesSpheroid%luminosity(weightIndex)
        end select
     end select
     ! Return if density is zero.
     if (Node_Component_Spheroid_Standard_Density <= 0.0d0) then
        Node_Component_Spheroid_Standard_Density=0.0d0
        return
     end if
     ! Compute actual density.
     position=[positionSpherical(1)/self%radius(),0.0d0,0.0d0]
     Node_Component_Spheroid_Standard_Density       =+Node_Component_Spheroid_Standard_Density                      &
          &                                          *massDistributionSpheroid                %density(position)
     if (Node_Component_Spheroid_Standard_Density > 0.0d0)                                                          &
          & Node_Component_Spheroid_Standard_Density=+Node_Component_Spheroid_Standard_Density                      &
          &                                          /self                                    %radius (        )**3
  end select
  return
end function Node_Component_Spheroid_Standard_Density

double precision function Node_Component_Spheroid_Standard_Density_Spherical_Average(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the spherically-averaged density at a given radius for an standard spheroid.
  !!}
  use :: Coordinates                          , only : assignment(=)           , coordinateSpherical
  use :: Galactic_Structure_Options           , only : componentTypeAll        , componentTypeSpheroid       , massTypeAll            , massTypeBaryonic       , &
          &                                            massTypeGalactic        , massTypeGaseous             , massTypeStellar        , weightByLuminosity     , &
          &                                            weightByMass            , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionSpheroid
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  type            (enumerationWeightByType      ), intent(in   ) :: weightBy
  integer                                        , intent(in   ) :: weightIndex
  double precision                               , intent(in   ) :: radius
  double precision                               , parameter     :: radiusHuge          =1.0d+100
  type            (coordinateSpherical          )                :: position
  type            (stellarLuminosities          ), save          :: luminositiesSpheroid
  !$omp threadprivate(luminositiesSpheroid)

  Node_Component_Spheroid_Standard_Density_Spherical_Average=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return

  ! Get the spheroid component and check that it is of the standard class.
  select type (self)
     class is (nodeComponentSpheroidStandard)

     if (self%radius() <= 0.0d0 .or. self%radius() > radiusHuge) return
     select case (weightBy%ID)
     case (weightByMass      %ID)
        select case (massType%ID)
        case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
           Node_Component_Spheroid_Standard_Density_Spherical_Average=self%massGas()+self%massStellar()
        case (massTypeGaseous%ID)
           Node_Component_Spheroid_Standard_Density_Spherical_Average=self%massGas()
        case (massTypeStellar%ID)
           Node_Component_Spheroid_Standard_Density_Spherical_Average=               self%massStellar()
        end select
     case (weightByLuminosity%ID)
        select case (massType%ID)
        case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
           luminositiesSpheroid=self%luminositiesStellar()
           Node_Component_Spheroid_Standard_Density_Spherical_Average=luminositiesSpheroid%luminosity(weightIndex)
        end select
     end select
     ! Return if density is zero.
     if (Node_Component_Spheroid_Standard_Density_Spherical_Average <= 0.0d0) then
        Node_Component_Spheroid_Standard_Density_Spherical_Average=0.0d0
        return
     end if
     ! Compute actual density.
     position=[radius/self%radius(),0.0d0,0.0d0]
     Node_Component_Spheroid_Standard_Density_Spherical_Average       =+Node_Component_Spheroid_Standard_Density_Spherical_Average     &
          &                                                            *massDistributionSpheroid                %density(position)
     if (Node_Component_Spheroid_Standard_Density_Spherical_Average > 0.0d0)                                                           &
          & Node_Component_Spheroid_Standard_Density_Spherical_Average=+Node_Component_Spheroid_Standard_Density_Spherical_Average     &
          &                                                             /self                                    %radius (        )**3
  end select
  return
end function Node_Component_Spheroid_Standard_Density_Spherical_Average

double precision function Node_Component_Spheroid_Standard_Potential(self,radius,componentType,massType,status)
  use :: Coordinates                          , only : assignment(=)                  , coordinateSpherical
  use :: Galactic_Structure_Options           , only : radiusLarge                    , weightByMass                     , weightIndexNull, enumerationComponentTypeType, &
       &                                               enumerationMassTypeType        , enumerationStructureErrorCodeType
  use :: Node_Component_Spheroid_Standard_Data, only : massDistributionSpheroid
  use :: Numerical_Constants_Astronomical     , only : gravitationalConstantGalacticus
  !!{
  Return the potential due to the standard spheroid.
  !!}
  class           (nodeComponentSpheroidStandard    ), intent(inout)           :: self
  type            (enumerationComponentTypeType     ), intent(in   )           :: componentType
  type            (enumerationMassTypeType          ), intent(in   )           :: massType
  double precision                                   , intent(in   )           :: radius
  type            (enumerationStructureErrorCodeType), intent(inout), optional :: status
  double precision                                                             :: componentMass
  type            (coordinateSpherical              )                          :: position
  !$GLC attributes unused :: status

  ! Set to zero by default.
  Node_Component_Spheroid_Standard_Potential=0.0d0

  ! Return immediately for non-positive radius.
  if (radius <= 0.0d0) return

  ! Get the spheroid component and check that it is of the standard class.
  select type (self)
  class is (nodeComponentSpheroidStandard)
     ! Compute if a spheroid is present.
     componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass&
          &,weightIndexNull)
     if (componentMass > 0.0d0 .and. self%radius() > 0.0d0) then
        position=[radius/self%radius(),0.0d0,0.0d0]
        Node_Component_Spheroid_Standard_Potential=(gravitationalConstantGalacticus*componentMass&
             &/self%radius())*massDistributionSpheroid%potential(position)
     end if
  end select
  return
end function Node_Component_Spheroid_Standard_Potential

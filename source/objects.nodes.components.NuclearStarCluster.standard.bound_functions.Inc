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
Contains custom functions for the standard nuclear star cluster component.
!!}

double precision function Node_Component_NSC_Standard_Velocity(self)
   !!{
   Return the circular velocity of the standard nuclear star cluster.
   !!}
   use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
   implicit none
   class(nodeComponentNSCStandard), intent(inout) :: self
   double precision                               :: radiusNSC, massNSC

   radiusNSC= self%radius()

   massNSC  = self%massGas()+self%massStellar()
   !Units corrected, G in galacticus is in units of Mpc \cdot M_\odot^{-1} \cdot (km/s)^2.
   !The velocity is returned in km/s units
   if (massNSC > 0.0d0 .and. radiusNSC >0.0d0) then
      Node_Component_NSC_Standard_Velocity = sqrt(gravitationalConstantGalacticus*massNSC/radiusNSC)
   else 
      Node_Component_NSC_Standard_Velocity = 0.0
   end if 
   return
end function Node_Component_NSC_Standard_Velocity

double precision function Node_Component_NSC_Standard_Angular_Momentum(self)
   !!{
   Return the specific angular momentum of the standard nuclear star cluster at a given radius.
   !!}
   implicit none
   class (nodeComponentNSCStandard), intent(inout) :: self
   double precision                                :: dynamicalMassNSC          , velocityNSC, &
         &                                            specificAngularMomentumNSC

   !Get the velocity and the dunamical mass of the nuclear star cluster.
   velocityNSC      = self%velocity()
   dynamicalMassNSC = self%massGas ()+self%massStellar()

   if (dynamicalMassNSC > 0.0d0) then
       Node_Component_NSC_Standard_Angular_Momentum = dynamicalMassNSC * velocityNSC
   else
       Node_Component_NSC_Standard_Angular_Momentum = 0.0d0
   end if
   return
end function Node_Component_NSC_Standard_Angular_Momentum

double precision function Node_Component_NSC_Standard_Radius(self)
   !!{
      Return the radius of the standard nuclear star cluster.
   !!}
   use :: Node_Component_NSC_Standard_Data, only: radiusNorm
   implicit none 
   class(nodeComponentNSCStandard), intent(inout) :: self 
   double precision               , save          :: M_0 = 1.0d6
   double precision                               :: dynamicalMassNSC

   dynamicalMassNSC = self%massGas()+self%massStellar()
   !The radius is retunerd in Mpc units.
   if (dynamicalMassNSC >0.0d0) then
      Node_Component_NSC_Standard_Radius = radiusNorm*sqrt(dynamicalMassNSC/M_0)
   else
      Node_Component_NSC_Standard_Radius = 0.0d0
   end if
   return
end function Node_Component_NSC_Standard_Radius

double precision function Node_Component_NSC_Standard_Half_Mass_Radius(self)
  !!{
  Return the half-mass radius of the standard nuclear star cluster.
  !!}
  use :: Error                           , only : Error_Report
  use :: Mass_Distributions              , only : massDistributionSpherical
  use :: Node_Component_NSC_Standard_Data, only : massDistributionNSC
  implicit none
  class(nodeComponentNSCStandard), intent(inout) :: self

 select type (massDistributionNSC)
  class is (massDistributionSpherical)
     Node_Component_NSC_Standard_Half_Mass_Radius=self%radius()*massDistributionNSC%radiusHalfMass()
  class default
     Node_Component_NSC_Standard_Half_Mass_Radius=0.0d0
     call Error_Report('mass distribution should be spherical'//{introspection:location})
  end select
  return
end function Node_Component_NSC_Standard_Half_Mass_Radius

double precision function Node_Component_NSC_Standard_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the mass within a given radius for an standard nuclear star cluster.
  !!}
  use :: Galactic_Structure_Options      , only : componentTypeAll       , componentTypeNSC, massTypeAll                 , massTypeBaryonic       , &
                            &                     massTypeGalactic       , massTypeGaseous , massTypeStellar             , radiusLarge            , &
                            &                     weightByLuminosity     , weightByMass    , enumerationComponentTypeType, enumerationMassTypeType, &
                            &                     enumerationWeightByType
  use :: Node_Component_NSC_Standard_Data, only : massDistributionNSC
  implicit none
  class           (nodeComponentNSCStandard    ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  type            (enumerationWeightByType     ), intent(in   ) :: weightBy
  integer                                       , intent(in   ) :: weightIndex
  double precision                              , intent(in   ) :: radius
  double precision                                              :: radiusNSC, fractionalRadius
  type            (stellarLuminosities         ), save          :: luminositiesNSC
  !$omp threadprivate(luminositiesNSC)

  ! Return immediately if nuclear star cluster component is not requested.
  Node_Component_NSC_Standard_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeNSC)) return
  ! Get the total mass.
  select case (weightBy%ID)
  case (weightByMass      %ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
        Node_Component_NSC_Standard_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous%ID)
        Node_Component_NSC_Standard_Enclosed_Mass=self%massGas()
     case (massTypeStellar%ID)
        Node_Component_NSC_Standard_Enclosed_Mass=               self%massStellar()
     end select
  case (weightByLuminosity%ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
        luminositiesNSC=self%luminositiesStellar()
        Node_Component_NSC_Standard_Enclosed_Mass   =luminositiesNSC%luminosity(weightIndex)
     end select
  end select
  ! Return if no mass.
  if (Node_Component_NSC_Standard_Enclosed_Mass <=       0.0d0) return
  ! Return if the total mass was requested.
  if (radius                                     >= radiusLarge) return
  ! Compute the actual mass.
  radiusNSC=self%radius()
  if (radiusNSC > 0.0d0) then
     fractionalRadius=radius/radiusNSC
     Node_Component_NSC_Standard_Enclosed_Mass=                         &
          & +Node_Component_NSC_Standard_Enclosed_Mass                  &
          & *massDistributionNSC%massEnclosedBySphere(fractionalRadius)
  end if
  return
end function Node_Component_NSC_Standard_Enclosed_Mass

function Node_Component_NSC_Standard_Acceleration(self,positionCartesian,componentType,massType)
  !!{
  Computes the gravitational acceleration at a given position for a standard nuclear star cluster.
  !!}
  use :: Coordinates                     , only : assignment(=)          , coordinateCartesian
  use :: Galactic_Structure_Options      , only : componentTypeAll       , componentTypeNSC   , massTypeAll     , massTypeBaryonic            , &
       &                                          massTypeGalactic       , massTypeGaseous    , massTypeStellar , enumerationComponentTypeType, &
       &                                          enumerationMassTypeType
  use :: Node_Component_NSC_Standard_Data, only : massDistributionNSC
  use :: Numerical_Constants_Astronomical, only : gigaYear               , megaParsec         , gravitationalConstantGalacticus
  use :: Numerical_Constants_Prefixes    , only : kilo
  implicit none
  double precision                                              , dimension(3)  :: Node_Component_NSC_Standard_Acceleration
  class           (nodeComponentNSCStandard     ), intent(inout)                :: self
  double precision                               , intent(in   ), dimension(3)  :: positionCartesian
  type            (enumerationComponentTypeType ), intent(in   )                :: componentType
  type            (enumerationMassTypeType      ), intent(in   )                :: massType
  double precision                                                              :: radius                                   , massTotal
  type            (coordinateCartesian          )                               :: positionScaleFree

  ! Return if the nuclear star cluster component is not selected.
  Node_Component_NSC_Standard_Acceleration=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeNSC)) return
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
  Node_Component_NSC_Standard_Acceleration=+kilo                                                     &
         &                                  *gigaYear                                                &
         &                                  /megaParsec                                              &
         &                                  *gravitationalConstantGalacticus                         &
         &                                  *massTotal                                               &
         &                                  /self%radius()                                       **2 &
         &                                  *massDistributionNSC%acceleration(positionScaleFree)
  return
end function Node_Component_NSC_Standard_Acceleration

function Node_Component_NSC_Standard_Tidal_Tensor(self,positionCartesian,componentType,massType)
  !!{
  Computes the gravitational acceleration at a given position for a standard nuclear star cluster.
  !!}
  use :: Coordinates                     , only : assignment(=)                  , coordinateCartesian
  use :: Galactic_Structure_Options      , only : componentTypeAll               , componentTypeNSC   , massTypeAll    , massTypeBaryonic            , &
                         &                        massTypeGalactic               , massTypeGaseous    , massTypeStellar, enumerationComponentTypeType, &
                         &                        enumerationMassTypeType
  use :: Node_Component_NSC_Standard_Data, only : massDistributionNSC
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  use :: Tensors                         , only : tensorRank2Dimension3Symmetric , tensorNullR2D3Sym  , operator(*)
  implicit none
  type            (tensorRank2Dimension3Symmetric)                              :: Node_Component_NSC_Standard_Tidal_Tensor
  class           (nodeComponentNSCStandard      ), intent(inout)               :: self
  double precision                                , intent(in   ), dimension(3) :: positionCartesian
  type            (enumerationComponentTypeType  ), intent(in   )               :: componentType
  type            (enumerationMassTypeType       ), intent(in   )               :: massType
  double precision                                                              :: radius                                  , massTotal
  type            (coordinateCartesian           )                              :: positionScaleFree

  ! Return if the nuclear star cluster component is not selected.
  Node_Component_NSC_Standard_Tidal_Tensor=tensorNullR2D3Sym
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeNSC)) return
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
  Node_Component_NSC_Standard_Tidal_Tensor=+gravitationalConstantGalacticus                        &
         &                                  *massTotal                                              &
         &                                  /self%radius()                                      **3 &
         &                                  *massDistributionNSC%tidalTensor(positionScaleFree)
  return
end function Node_Component_NSC_Standard_Tidal_Tensor

double precision function Node_Component_NSC_Standard_Density(self,positionSpherical,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the density at a given position for an standard nuclear star cluster.
  !!}
  use :: Coordinates                     , only : assignment(=)      , coordinateSpherical
  use :: Galactic_Structure_Options      , only : componentTypeAll   , componentTypeNSC            , massTypeAll            , massTypeBaryonic       , &
          &                                       massTypeGalactic   , massTypeGaseous             , massTypeStellar        , weightByLuminosity     , &
          &                                       weightByMass       , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Node_Component_NSC_Standard_Data, only : massDistributionNSC
  implicit none
  class           (nodeComponentNSCStandard     ), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  type            (enumerationWeightByType      ), intent(in   ) :: weightBy
  integer                                        , intent(in   ) :: weightIndex
  double precision                               , intent(in   ) :: positionSpherical(3)
  type            (stellarLuminosities          ), save          :: luminositiesNSC
  !$omp threadprivate(luminositiesNSC)
  type            (coordinateSpherical          )                :: position

  ! Return immediately if nuclear star cluster component is not requested.
  Node_Component_NSC_Standard_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeNSC)) return
  ! Determine mass/luminosity type.
  select case (weightBy%ID)
  case (weightByMass      %ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
        Node_Component_NSC_Standard_Density=self%massGas()+self%massStellar()
     case (massTypeGaseous%ID)
        Node_Component_NSC_Standard_Density=self%massGas()
     case (massTypeStellar%ID)
        Node_Component_NSC_Standard_Density=               self%massStellar()
     end select
  case (weightByLuminosity%ID)
     select case (massType%ID)
     case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
        luminositiesNSC=self%luminositiesStellar()
        Node_Component_NSC_Standard_Density=luminositiesNSC%luminosity(weightIndex)
     end select
  end select
  ! Skip further calculation if mass or radius is zero.
  if (Node_Component_NSC_Standard_Density > 0.0d0 .and. self%radius() > 0.0d0) then
     ! Compute the actual density.
      position=[positionSpherical(1)/self%radius(),positionSpherical(2),positionSpherical(3)]
      Node_Component_NSC_Standard_Density=+Node_Component_NSC_Standard_Density      &
           &                               /self                %radius (        )**3 &
           &                               *massDistributionNSC%density(position)
   end if
   return
end function Node_Component_NSC_Standard_Density


double precision function Node_Component_NSC_Standard_Density_Spherical_Average(self,radius,componentType,massType,weightBy,weightIndex)
  !!{
  Computes the spherically-averaged density at a given radius for an standard NSC.
  !!}
  use :: Coordinates                     , only : assignment(=)           , coordinateSpherical
  use :: Galactic_Structure_Options      , only : componentTypeAll        , componentTypeNSC            , massTypeAll            , massTypeBaryonic       , &
          &                                       massTypeGalactic        , massTypeGaseous             , massTypeStellar        , weightByLuminosity     , &
          &                                       weightByMass            , enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType
  use :: Node_Component_NSC_Standard_Data, only : massDistributionNSC
  implicit none
  class           (nodeComponentNSCStandard     ), intent(inout) :: self
  type            (enumerationComponentTypeType ), intent(in   ) :: componentType
  type            (enumerationMassTypeType      ), intent(in   ) :: massType
  type            (enumerationWeightByType      ), intent(in   ) :: weightBy
  integer                                        , intent(in   ) :: weightIndex
  double precision                               , intent(in   ) :: radius
  double precision                               , parameter     :: radiusHuge          =1.0d+100
  type            (coordinateSpherical          )                :: position
  type            (stellarLuminosities          ), save          :: luminositiesNSC
  !$omp threadprivate(luminositiesNSC)

  Node_Component_NSC_Standard_Density_Spherical_Average=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeNSC)) return

  ! Get the spheroid component and check that it is of the standard class.
  select type (self)
     class is (nodeComponentNSCStandard)

     if (self%radius() <= 0.0d0 .or. self%radius() > radiusHuge) return
     select case (weightBy%ID)
     case (weightByMass      %ID)
        select case (massType%ID)
        case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID)
           Node_Component_NSC_Standard_Density_Spherical_Average=self%massGas()+self%massStellar()
        case (massTypeGaseous%ID)
           Node_Component_NSC_Standard_Density_Spherical_Average=self%massGas()
        case (massTypeStellar%ID)
           Node_Component_NSC_Standard_Density_Spherical_Average=               self%massStellar()
        end select
     case (weightByLuminosity%ID)
        select case (massType%ID)
        case (massTypeAll%ID,massTypeBaryonic%ID,massTypeGalactic%ID,massTypeStellar%ID)
           luminositiesNSC=self%luminositiesStellar()
           Node_Component_NSC_Standard_Density_Spherical_Average=luminositiesNSC%luminosity(weightIndex)
        end select
     end select
     ! Return if density is zero.
     if (Node_Component_NSC_Standard_Density_Spherical_Average <= 0.0d0) then
         Node_Component_NSC_Standard_Density_Spherical_Average=0.0d0
        return
     end if
     ! Compute actual density.
     position=[radius/self%radius(),0.0d0,0.0d0]
     Node_Component_NSC_Standard_Density_Spherical_Average       =+Node_Component_NSC_Standard_Density_Spherical_Average     &
          &                                                            *massDistributionNSC            %density(position)
     if (Node_Component_NSC_Standard_Density_Spherical_Average > 0.0d0)                                                      &
          & Node_Component_NSC_Standard_Density_Spherical_Average=+Node_Component_NSC_Standard_Density_Spherical_Average     &
          &                                                             /self                          %radius (        )**3
  end select
  return
end function Node_Component_NSC_Standard_Density_Spherical_Average

double precision function Node_Component_NSC_Standard_Potential(self,radius,componentType,massType,status)
  !!{
  Compute the gravitational potential due to an standard nuclear star cluster.
  !!}
  use :: Coordinates                     , only : assignment(=)                  , coordinateCylindrical
  use :: Galactic_Structure_Options      , only : componentTypeAll               , componentTypeNSC            , radiusLarge            , weightByMass                     , &
          &                                       weightIndexNull                , enumerationComponentTypeType, enumerationMassTypeType, enumerationStructureErrorCodeType
  use :: Node_Component_NSC_Standard_Data, only : massDistributionNSC
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  implicit none
  class           (nodeComponentNSCStandard         ), intent(inout)           :: self
  type            (enumerationComponentTypeType     ), intent(in   )           :: componentType
  type            (enumerationMassTypeType          ), intent(in   )           :: massType
  double precision                                   , intent(in   )           :: radius
  type            (enumerationStructureErrorCodeType), intent(inout), optional :: status
  double precision                                                             :: componentMass
  type            (coordinateCylindrical            )                          :: position
  !$GLC attributes unused :: status

  ! Return immediately if nuclear star cluster component is not requested.
  Node_Component_NSC_Standard_Potential=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeNSC)) return
  ! Avoid an arithmetic exception at radius zero.
  if (radius <= 0.0d0) return
  ! Get the relevant mass of the nuclear star cluster.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull)
  if (componentMass <= 0.0d0) return
  ! Check for zero-sized nuclear star cluster.
  if (self%radius() <= 0.0d0) then
     ! Treat as a point mass.
     Node_Component_NSC_Standard_Potential=-gravitationalConstantGalacticus           &
          &                                 *componentMass                            &
          &                                 /radius
  else
     ! Compute the potential.
     position=[radius/self%radius(),0.0d0,0.0d0]
     Node_Component_NSC_Standard_Potential=+gravitationalConstantGalacticus           &
          &                                 *componentMass                            &
          &                                 /self%radius()                            &
          &                                 *massDistributionNSC%potential(position)
  end if
  return
end function Node_Component_NSC_Standard_Potential

double precision function Node_Component_NSC_Standard_Rotation_Curve(self,radius,componentType,massType)
  !!{
  Computes the rotation curve at a given radius for an standard nuclear star cluster.
  !!}
  use :: Galactic_Structure_Options      , only : radiusLarge            , weightByMass, weightIndexNull, enumerationComponentTypeType, &
       &                                          enumerationMassTypeType
  use :: Node_Component_NSC_Standard_Data, only : massDistributionNSC
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  implicit none
  class           (nodeComponentNSCStandard    ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  double precision                              , intent(in   ) :: radius
  double precision                                              :: componentMass, radiusNSC, fractionalRadius

  ! Set to zero by default.
  Node_Component_NSC_Standard_Rotation_Curve=0.0d0

  ! Return inmediately for non-positive radius.
  if (radius <= 0.0d0) return

  ! Gompute if a NSC is present.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull)
  if (componentMass > 0.0d0) Node_Component_NSC_Standard_Rotation_Curve=sqrt(gravitationalConstantGalacticus&
       &*componentMass)/sqrt(radius)
  return
end function Node_Component_NSC_Standard_Rotation_Curve

double precision function Node_Component_NSC_Standard_Rotation_Curve_Gradient(self,radius,componentType,massType)
  !!{
  Computes the rotation curve gradient for an standard NSC.
  !!}
  use :: Galactic_Structure_Options      , only : radiusLarge                    , weightByMass, weightIndexNull, enumerationComponentTypeType, &
       &                                          enumerationMassTypeType
  use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
  use :: Numerical_Constants_Math        , only : Pi

  implicit none
  class           (nodeComponentNSCStandard    ), intent(inout) :: self
  type            (enumerationComponentTypeType), intent(in   ) :: componentType
  type            (enumerationMassTypeType     ), intent(in   ) :: massType
  double precision                              , intent(in   ) :: radius
  double precision                                              :: componentDensity, componentMass, positionSpherical(3)

  ! Set to zero by default.
  Node_Component_NSC_Standard_Rotation_Curve_Gradient=0.0d0
  
  ! Return if radius is zero.
  if (radius <= 0.0d0) return
  
  ! Get the mass of the nuclear star cluster.
  positionSpherical= [radius,0.0d0,0.0d0]
  componentMass    =self%enclosedMass(radius           ,componentType,massType,weightByMass,weightIndexNull)
  componentDensity =self%density     (positionSpherical,componentType,massType,weightByMass,weightIndexNull)
  if (componentMass == 0.0d0 .or. componentDensity == 0.0d0) return
  Node_Component_NSC_Standard_Rotation_Curve_Gradient=       &
       &                  gravitationalConstantGalacticus    &
       &                 *(                                  &
       &                   -componentMass/radius**2          &
       &                   +4.0d0*Pi*radius*componentDensity &
       &                  )
  return
end function Node_Component_NSC_Standard_Rotation_Curve_Gradient

function Node_Component_NSC_Standard_Chandrasekhar_Integral(self,nodeSatellite,positionCartesian,velocityCartesian,componentType,massType)
  !!{
  Computes the gravitational acceleration at a given position for a standard NSC.
  !!}
  use :: Galactic_Structure_Options, only : componentTypeAll            , componentTypeNSC       , weightByMass, weightIndexNull, &
       &                                    enumerationComponentTypeType, enumerationMassTypeType
  use :: Numerical_Constants_Math  , only : Pi
  implicit none
  double precision                                             , dimension(3) :: Node_Component_NSC_Standard_Chandrasekhar_Integral
  class           (nodeComponentNSCStandard    ), intent(inout)               :: self
  type            (treeNode                    ), intent(inout)               :: nodeSatellite
  double precision                              , intent(in   ), dimension(3) :: positionCartesian                                 , velocityCartesian
  type            (enumerationComponentTypeType), intent(in   )               :: componentType
  type            (enumerationMassTypeType     ), intent(in   )               :: massType
  double precision                                             , dimension(3) :: positionSpherical
  double precision                              , parameter                   :: XvMaximum                                         =10.0d0
  double precision                                                            :: radius                                            , velocity          , &
       &                                                                         density                                           , xV                , &
       &                                                                         velocityRotation                                  , velocityDispersion
  !$GLC attributes unused :: nodeSatellite
  
  ! Return if the spheroid component is not selected.
  Node_Component_NSC_Standard_Chandrasekhar_Integral=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeNSC)) return
  ! Compute the integral.
  radius                                             =  sqrt(sum(positionCartesian**2))
  velocity                                           =  sqrt(sum(velocityCartesian**2))
  if (velocity <= 0.0d0) return
  positionSpherical                                  =  [radius,0.0d0,0.0d0]
  velocityRotation                                   =  Node_Component_NSC_Standard_Rotation_Curve(self,radius           ,componentType,massType                             )
  density                                            =  Node_Component_NSC_Standard_Density       (self,positionSpherical,componentType,massType,weightByMass,weightIndexNull)
  if (density  <= 0.0d0) return
  velocityDispersion                                 = +velocityRotation   &
       &                                               /sqrt(2.0d0)
  xV                                                 = +velocity           &
       &                                               /velocityDispersion &
       &                                               /sqrt(2.0d0)
  Node_Component_NSC_Standard_Chandrasekhar_Integral = -density            &
       &                                               *velocityCartesian  &
       &                                               /velocity         **3
  if (Xv <= XvMaximum)                                                                                          &
       & Node_Component_NSC_Standard_Chandrasekhar_Integral=+Node_Component_NSC_Standard_Chandrasekhar_Integral &
       &                                                         *(                                             &
       &                                                           +erf ( xV   )                                &
       &                                                           -2.0d0                                       &
       &                                                           *      xV                                    &
       &                                                           *exp (-xV**2)                                &
       &                                                           /sqrt( Pi   )                                &
       &                                                         )
  return
end function Node_Component_NSC_Standard_Chandrasekhar_Integral

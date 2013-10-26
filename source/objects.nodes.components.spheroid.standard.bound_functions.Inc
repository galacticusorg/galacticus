!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains custom functions for the standard spheroid component.

double precision function Node_Component_Spheroid_Standard_Half_Mass_Radius(self)
  !% Return the half-mass radius of the standard spheroid.
  use Node_Component_Spheroid_Standard_Data
  use Mass_Distributions
  implicit none
  class(nodeComponentSpheroidStandard), intent(inout) :: self

  select type (spheroidMassDistribution)
     class is (massDistributionSpherical)
     Node_Component_Spheroid_Standard_Half_Mass_Radius=self%radius()*spheroidMassDistribution%halfMassRadius()
  end select
  return
end function Node_Component_Spheroid_Standard_Half_Mass_Radius

double precision function Node_Component_Spheroid_Standard_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
  !% Computes the mass within a given radius for an standard spheroid.
  use Node_Component_Spheroid_Standard_Data
  use Galactic_Structure_Options
  use Node_Component_Spheroid_Standard_Data
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout)           :: self
  integer                                        , intent(in   )           :: componentType   , massType      , weightBy, weightIndex
  double precision                               , intent(in   )           :: radius
  logical                                        , intent(in   ), optional :: haloLoaded
  double precision                                                         :: fractionalRadius, spheroidRadius
  type            (stellarLuminosities          )                          :: luminositiesSpheroid

  ! Return if the spheroid component is not selected.
  Node_Component_Spheroid_Standard_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return
  ! Get the total mass.
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        Node_Component_Spheroid_Standard_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous)
        Node_Component_Spheroid_Standard_Enclosed_Mass=self%massGas()
     case (massTypeStellar)
        Node_Component_Spheroid_Standard_Enclosed_Mass=               self%massStellar()
     end select
  case (weightByLuminosity)
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
        luminositiesSpheroid=self%luminositiesStellar()
        Node_Component_Spheroid_Standard_Enclosed_Mass       =luminositiesSpheroid%luminosity(weightIndex)
     end select
  end select
  ! Return if total mass was requested.
  if (radius                                         >= radiusLarge) return
  ! Return if mass is zero.
  if (Node_Component_Spheroid_Standard_Enclosed_Mass <=       0.0d0) return
  ! Compute actual mass.
  spheroidRadius=self%radius()
  if (spheroidRadius > 0.0d0) then
     fractionalRadius=radius/spheroidRadius
     Node_Component_Spheroid_Standard_Enclosed_Mass   =Node_Component_Spheroid_Standard_Enclosed_Mass*spheroidMassDistribution%massEnclosedBySphere(fractionalRadius)
  end if
  return
end function Node_Component_Spheroid_Standard_Enclosed_Mass

double precision function Node_Component_Spheroid_Standard_Rotation_Curve(self,radius,componentType,massType,haloLoaded)
  !% Computes the rotation curve at a given radius for a standard spheroid.
  use Galactic_Structure_Options
  use Numerical_Constants_Physical
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout)           :: self
  integer                                        , intent(in   )           :: componentType, massType
  double precision                               , intent(in   )           :: radius
  logical                                        , intent(in   ), optional :: haloLoaded
  double precision                                                         :: componentMass

  ! Set to zero by default.
  Node_Component_Spheroid_Standard_Rotation_Curve=0.0d0

  ! Return immediately for non-positive radius.
  if (radius <= 0.0d0) return

  ! Compute if a spheroid is present.
  componentMass=self%enclosedMass(radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
  if (componentMass > 0.0d0) Node_Component_Spheroid_Standard_Rotation_Curve=sqrt(gravitationalConstantGalacticus&
       &*componentMass)/sqrt(radius)
  return
end function Node_Component_Spheroid_Standard_Rotation_Curve

double precision function Node_Component_Spheroid_Standard_Rotation_Curve_Gradient(self,radius,componentType,massType,haloLoaded)
  !% Computes the rotation curve gradient for the standard spheroid.
  use Galactic_Structure_Options
  use Numerical_Constants_Physical
  use Numerical_Constants_Prefixes
  use Numerical_Constants_Math
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout)           :: self
  integer                                        , intent(in   )           :: componentType   , massType
  double precision                               , intent(in   )           :: radius
  logical                                        , intent(in   ), optional :: haloLoaded
  double precision                                                         :: componentDensity, componentMass, positionSpherical(3)

  ! Set to zero by default.
  Node_Component_Spheroid_Standard_Rotation_Curve_Gradient=0.0d0

  ! Return immediatelt for non-positive radius.
  if (radius <= 0.0d0) return

  ! Compute if a spheroid is present.
  positionSpherical= [radius,0.0d0,0.0d0]
  componentMass    =self%enclosedMass(radius           ,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
  componentDensity =self%density     (positionSpherical,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
  if (componentMass == 0.0d0 .or. componentDensity == 0.0d0) return
  Node_Component_Spheroid_Standard_Rotation_Curve_Gradient=  &
       &                  gravitationalConstantGalacticus    &
       &                 *(                                  &
       &                   -componentMass/radius**2          &
       &                   +4.0d0*Pi*radius*componentDensity &
       &                  )
  return
end function Node_Component_Spheroid_Standard_Rotation_Curve_Gradient

double precision function Node_Component_Spheroid_Standard_Density(self,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
  !% Computes the density at a given position for an standard spheroid.
  use Galactic_Structure_Options
  use Numerical_Constants_Math
  use Coordinates
  use Node_Component_Spheroid_Standard_Data
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout)           :: self
  integer                                        , intent(in   )           :: componentType       , massType, weightBy, &
       &                                                                      weightIndex
  double precision                               , intent(in   )           :: positionSpherical(3)
  logical                                        , intent(in   ), optional :: haloLoaded
  type            (coordinateSpherical          )                          :: position
  type            (stellarLuminosities          )                          :: luminositiesSpheroid

  Node_Component_Spheroid_Standard_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeSpheroid)) return

  ! Get the spheroid component and check that it is of the standard class.
  select type (self)
     class is (nodeComponentSpheroidStandard)

     if (self%radius() <= 0.0d0) return
     select case (weightBy)
     case (weightByMass      )
        select case (massType)
        case (massTypeAll,massTypeBaryonic,massTypeGalactic)
           Node_Component_Spheroid_Standard_Density=self%massGas()+self%massStellar()
        case (massTypeGaseous)
           Node_Component_Spheroid_Standard_Density=self%massGas()
        case (massTypeStellar)
           Node_Component_Spheroid_Standard_Density=               self%massStellar()
        end select
     case (weightByLuminosity)
        select case (massType)
        case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
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
     Node_Component_Spheroid_Standard_Density=(Node_Component_Spheroid_Standard_Density/self%radius()**3)&
          &*spheroidMassDistribution%density(position)
  end select
  return
end function Node_Component_Spheroid_Standard_Density

double precision function Node_Component_Spheroid_Standard_Potential(self,radius,componentType,massType,haloLoaded)
  use Galactic_Structure_Options
  use Numerical_Constants_Physical
  use Coordinates
  use Node_Component_Spheroid_Standard_Data
  !% Return the potential due to the standard spheroid.
  class           (nodeComponentSpheroidStandard), intent(inout)           :: self
  integer                                        , intent(in   )           :: componentType, massType
  double precision                               , intent(in   )           :: radius
  logical                                        , intent(in   ), optional :: haloLoaded
  double precision                                                         :: componentMass
  type            (coordinateSpherical  )                                  :: position

  ! Set to zero by default.
  Node_Component_Spheroid_Standard_Potential=0.0d0

  ! Return immediately for non-positive radius.
  if (radius <= 0.0d0) return

  ! Get the spheroid component and check that it is of the standard class.
  select type (self)
     class is (nodeComponentSpheroidStandard)
        ! Compute if a spheroid is present.
     componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass&
          &,weightIndexNull,haloLoaded)
     if (componentMass > 0.0d0 .and. self%radius() > 0.0d0) then
        position=[radius/self%radius(),0.0d0,0.0d0]
        Node_Component_Spheroid_Standard_Potential=(gravitationalConstantGalacticus*componentMass&
             &/self%radius())*spheroidMassDistribution%potential(position)
     end if
  end select
  return
end function Node_Component_Spheroid_Standard_Potential

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
  implicit none
  class           (nodeComponentSpheroidStandard), intent(inout)           :: self
  integer                                        , intent(in   )           :: componentType   , massType      , weightBy, weightIndex
  double precision                               , intent(in   )           :: radius
  logical                                        , intent(in   ), optional :: haloLoaded
  double precision                                                         :: fractionalRadius, spheroidRadius

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
        Node_Component_Spheroid_Standard_Enclosed_Mass       =luminositiesSpheroid(weightIndex)
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

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains custom functions for the standard disk component.

subroutine Node_Component_Disk_Standard_Attach_Pipes(self)
  !% Attach cooling pipes to the standard disk component.
  use Galacticus_Error
  implicit none
  class(nodeComponentDiskStandard), intent(in   ) :: self
  type (nodeComponentHotHalo     )                :: hotHalo

  if (hotHalo%hotHaloCoolingMassRateIsAttached           ()) &
       call Galacticus_Error_Report('Node_Component_Disk_Standard_Attach_Pipes','expected to find unclaimed hot halo mass cooling pipe'            )
  if (hotHalo%hotHaloCoolingAngularMomentumRateIsAttached()) &
       call Galacticus_Error_Report('Node_Component_Disk_Standard_Attach_Pipes','expected to find unclaimed hot halo angular momentum cooling pipe')
  if (hotHalo%hotHaloCoolingAbundancesRateIsAttached     ()) &
       call Galacticus_Error_Report('Node_Component_Disk_Standard_Attach_Pipes','expected to find unclaimed hot halo abundances cooling pipe'      )
  call hotHalo%hotHaloCoolingMassRateFunction           (DiskStandardMassGasRateGeneric        )
  call hotHalo%hotHaloCoolingAngularMomentumRateFunction(DiskStandardAngularMomentumRateGeneric)
  call hotHalo%hotHaloCoolingAbundancesRateFunction     (DiskStandardAbundancesGasRateGeneric  )
  return
end subroutine Node_Component_Disk_Standard_Attach_Pipes

double precision function Node_Component_Disk_Standard_Half_Mass_Radius(self)
  !% Return the half-mass radius of the standard disk.
  use Node_Component_Disk_Standard_Data
  implicit none
  class(nodeComponentDiskStandard), intent(inout) :: self
double precision::a,b
  ! Check that the mass distribution is cylindrical.
  select type (diskMassDistribution)
  class is (massDistributionCylindrical)
     Node_Component_Disk_Standard_Half_Mass_Radius=self%radius()*diskMassDistribution%halfMassRadius()
  end select
  return
end function Node_Component_Disk_Standard_Half_Mass_Radius

double precision function Node_Component_Disk_Standard_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
  !% Computes the mass within a given radius for an standard disk.
  use Galactic_Structure_Options
  use Tables
  use Node_Component_Disk_Standard_Data
  implicit none
  class           (nodeComponentDiskStandard), intent(inout)           :: self
  integer                                    , intent(in   )           :: componentType   , massType        , weightBy, weightIndex
  double precision                           , intent(in   )           :: radius
  logical                                    , intent(in   ), optional :: haloLoaded
  double precision                                                     :: diskRadius      , fractionalRadius
  type            (stellarLuminosities      ), save                    :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Standard_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Get the total mass.
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        Node_Component_Disk_Standard_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous)
        Node_Component_Disk_Standard_Enclosed_Mass=self%massGas()
     case (massTypeStellar)
        Node_Component_Disk_Standard_Enclosed_Mass=               self%massStellar()
     end select
  case (weightByLuminosity)
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Standard_Enclosed_Mass   =luminositiesDisk%luminosity(weightIndex)
     end select
  end select
  ! Return if no mass.
  if (Node_Component_Disk_Standard_Enclosed_Mass <=       0.0d0) return
  ! Return if the total mass was requested.
  if (radius                                     >= radiusLarge) return
  ! Compute the actual mass.
  diskRadius=self%radius()
  if (diskRadius > 0.0d0) then
     fractionalRadius=radius/diskRadius
     Node_Component_Disk_Standard_Enclosed_Mass=                         &
          & +Node_Component_Disk_Standard_Enclosed_Mass                  &
          & *diskMassDistribution%massEnclosedBySphere(fractionalRadius)
  end if
  return
end function Node_Component_Disk_Standard_Enclosed_Mass

double precision function Node_Component_Disk_Standard_Density(self,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
  !% Computes the density at a given position for an standard disk.
  use Galactic_Structure_Options
  use Numerical_Constants_Math
  use Tables
  use Node_Component_Disk_Standard_Data
  use Coordinates
  implicit none
  class           (nodeComponentDiskStandard), intent(inout)           :: self
  integer                                    , intent(in   )           :: componentType       , massType   , &
       &                                                                  weightBy            , weightIndex
  double precision                           , intent(in   )           :: positionSpherical(3)
  logical                                    , intent(in   ), optional :: haloLoaded
  type            (stellarLuminosities      ), save                    :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  type            (coordinateSpherical      )                           :: position

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Standard_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Determine mass/luminosity type.
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        Node_Component_Disk_Standard_Density=self%massGas()+self%massStellar()
     case (massTypeGaseous)
        Node_Component_Disk_Standard_Density=self%massGas()
     case (massTypeStellar)
        Node_Component_Disk_Standard_Density=               self%massStellar()
     end select
  case (weightByLuminosity)
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
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

double precision function Node_Component_Disk_Standard_Potential(self,radius,componentType,massType,haloLoaded,status)
  !% Compute the gravitational potential due to an standard disk.
  use Numerical_Constants_Physical
  use Bessel_Functions
  use Galactic_Structure_Options
  use Tables
  use Node_Component_Disk_Standard_Data
  use Coordinates
  implicit none
  class           (nodeComponentDiskStandard), intent(inout)           :: self
  integer                                    , intent(in   )           :: componentType, massType
  double precision                           , intent(in   )           :: radius
  logical                                    , intent(in   ), optional :: haloLoaded
  integer                                    , intent(inout), optional :: status
  double precision                                                     :: componentMass
  type            (coordinateCylindrical    )                          :: position

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Standard_Potential=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Avoid an arithmetic exception at radius zero.
  if (radius <= 0.0d0) return
  ! Get the relevant mass of the disk.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
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

double precision function Node_Component_Disk_Standard_Rotation_Curve(self,radius,componentType,massType,haloLoaded)
  !% Computes the rotation curve at a given radius for an standard disk.
  use Tables
  use Node_Component_Disk_Standard_Data
  use Galactic_Structure_Options
  use Numerical_Constants_Physical
  implicit none
  class           (nodeComponentDiskStandard), intent(inout)           :: self
  integer                                    , intent(in   )           :: componentType         , massType
  double precision                           , intent(in   )           :: radius
  logical                                    , intent(in   ), optional :: haloLoaded
  double precision                                                     :: componentMass         , diskRadius, fractionalRadius

  ! Set to zero by default.
  Node_Component_Disk_Standard_Rotation_Curve=0.0d0
  ! Get the mass of the disk.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
  if (componentMass <= 0.0d0) return
  ! Check that the mass distribution is cylindrical.
  select type (diskMassDistribution)
  class is (massDistributionCylindrical)
     ! Compute the actual velocity.
     diskRadius=self%radius()
     if (diskRadius > 0.0d0) then        
        fractionalRadius=radius/diskRadius
        Node_Component_Disk_Standard_Rotation_Curve=+sqrt(                                                &
             &                                            +gravitationalConstantGalacticus                &
             &                                            *componentMass                                  &
             &                                            /diskRadius                                     &
             &                                           )                                                &
             &                                      *diskMassDistribution%rotationCurve(fractionalRadius)        
     end if
  end select
  return
end function Node_Component_Disk_Standard_Rotation_Curve

double precision function Node_Component_Disk_Standard_Rotation_Curve_Gradient(self,radius,componentType,massType,haloLoaded)
  !% Computes the rotation curve gradient for an standard disk.
  use Tables
  use Node_Component_Disk_Standard_Data
  use Galactic_Structure_Options
  use Numerical_Constants_Physical
  use Numerical_Constants_Prefixes
  implicit none
  class           (nodeComponentDiskStandard), intent(inout)           :: self
  integer                                    , intent(in   )           :: componentType, massType
  double precision                           , intent(in   )           :: radius
  logical                                    , intent(in   ), optional :: haloLoaded
  double precision                                                     :: diskRadius   , fractionalRadius, &
       &                                                                  componentMass

  ! Set to zero by default.
  Node_Component_Disk_Standard_Rotation_Curve_Gradient=0.0d0
  ! Return if radius is zero.
  if (radius <= 0.0d0) return
  ! Get the mass of the disk.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
  if (componentMass <= 0.0d0) return
  ! Check that the mass distribution is cylindrical.
  select type (diskMassDistribution)
  class is (massDistributionCylindrical)
     ! Compute the rotation curve gradient.
     diskRadius=self%radius()
     if (diskRadius > 0.0d0) then        
        fractionalRadius=radius/diskRadius
        Node_Component_Disk_Standard_Rotation_Curve_Gradient=+gravitationalConstantGalacticus                              &
             &                                               *componentMass                                                &
             &                                               /diskRadius**2                                                &
             &                                               *diskMassDistribution%rotationCurveGradient(fractionalRadius)  
     end if
  end select
  return
end function Node_Component_Disk_Standard_Rotation_Curve_Gradient

double precision function Node_Component_Disk_Standard_Surface_Density(self,positionCylindrical,componentType,massType&
     &,weightBy,weightIndex,haloLoaded)
  !% Computes the surface density at a given position for an standard disk.
  use Galactic_Structure_Options
  use Numerical_Constants_Math
  use Tables
  use Node_Component_Disk_Standard_Data
  use Coordinates
  implicit none
  class           (nodeComponentDiskStandard), intent(inout)           :: self
  integer                                    , intent(in   )           :: componentType         , massType   , &
       &                                                                  weightBy              , weightIndex
  double precision                           , intent(in   )           :: positionCylindrical(3)
  logical                                    , intent(in   ), optional :: haloLoaded
  type            (treeNode                 ), pointer                 :: selfNode
  type            (stellarLuminosities      ), save                    :: luminositiesDisk
  !$omp threadprivate(luminositiesDisk)
  type            (coordinateCylindrical    )                          :: position

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
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        if (.not.surfaceDensityCentralTotalComputed  ) then
           surfaceDensityCentralTotal          =(self%massGas()+self%massStellar())/radiusScaleDisk**2
           surfaceDensityCentralTotalComputed  =.true.
        end if
        Node_Component_Disk_Standard_Surface_Density=surfaceDensityCentralTotal
     case (massTypeGaseous)
        if (.not.surfaceDensityCentralGasComputed    ) then
           surfaceDensityCentralGas            = self%massGas()                    /radiusScaleDisk**2
           surfaceDensityCentralGasComputed    =.true.
        end if
        Node_Component_Disk_Standard_Surface_Density=surfaceDensityCentralGas
     case (massTypeStellar)
        if (.not.surfaceDensityCentralStellarComputed) then
           surfaceDensityCentralStellar        =                self%massStellar() /radiusScaleDisk**2
           surfaceDensityCentralStellarComputed=.true.
        end if
        Node_Component_Disk_Standard_Surface_Density=surfaceDensityCentralStellar
     end select
  case (weightByLuminosity)
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Standard_Surface_Density=luminositiesDisk%luminosity(weightIndex)/2.0d0/Pi/radiusScaleDisk**2
     end select
  end select
  ! Return if no density.
  if (Node_Component_Disk_Standard_Surface_Density <= 0.0d0) return
  ! Check that the mass distribution is cylindrical.
  select type (diskMassDistribution)
  class is (massDistributionCylindrical)
     ! Compute the surface density.
     position=[positionCylindrical(1)/self%radius(),0.0d0,0.0d0]
     Node_Component_Disk_Standard_Surface_Density=         &
          & +Node_Component_Disk_Standard_Surface_Density  &
          & *diskMassDistribution%surfaceDensity(position)  
  end select
  return
end function Node_Component_Disk_Standard_Surface_Density

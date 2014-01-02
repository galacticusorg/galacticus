!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains custom functions for the exponential disk component.

subroutine Node_Component_Disk_Exponential_Attach_Pipes(self)
  !% Attach cooling pipes to the exponential disk component.
  use Galacticus_Error
  implicit none
  class(nodeComponentDiskExponential), intent(in   ) :: self
  type (nodeComponentHotHalo        )                :: hotHalo

  if (hotHalo%hotHaloCoolingMassRateIsAttached           ()) &
       call Galacticus_Error_Report('Node_Component_Disk_Exponential_Attach_Pipes','expected to find unclaimed hot halo mass cooling pipe'            )
  if (hotHalo%hotHaloCoolingAngularMomentumRateIsAttached()) &
       call Galacticus_Error_Report('Node_Component_Disk_Exponential_Attach_Pipes','expected to find unclaimed hot halo angular momentum cooling pipe')
  if (hotHalo%hotHaloCoolingAbundancesRateIsAttached     ()) &
       call Galacticus_Error_Report('Node_Component_Disk_Exponential_Attach_Pipes','expected to find unclaimed hot halo abundances cooling pipe'      )
  call hotHalo%hotHaloCoolingMassRateFunction           (DiskExponentialMassGasRateGeneric        )
  call hotHalo%hotHaloCoolingAngularMomentumRateFunction(DiskExponentialAngularMomentumRateGeneric)
  call hotHalo%hotHaloCoolingAbundancesRateFunction     (DiskExponentialAbundancesGasRateGeneric  )
  return
end subroutine Node_Component_Disk_Exponential_Attach_Pipes

double precision function Node_Component_Disk_Exponential_Half_Mass_Radius(self)
  !% Return the half-mass radius of the exponential disk.
  implicit none
  class           (nodeComponentDiskExponential), intent(inout) :: self
  double precision                              , parameter     :: halfMassRadiusToScaleRadius=1.678346990d0

  Node_Component_Disk_Exponential_Half_Mass_Radius=self%radius()*halfMassRadiusToScaleRadius
  return
end function Node_Component_Disk_Exponential_Half_Mass_Radius

double precision function Node_Component_Disk_Exponential_Enclosed_Mass(self,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
  !% Computes the mass within a given radius for an exponential disk.
  use Galactic_Structure_Options
  use Tables
  use Node_Component_Disk_Exponential_Data
  implicit none
  class           (nodeComponentDiskExponential), intent(inout)           :: self
  integer                                       , intent(in   )           :: componentType   , massType        , weightBy, weightIndex
  double precision                              , intent(in   )           :: radius
  logical                                       , intent(in   ), optional :: haloLoaded
  double precision                                                        :: diskRadius      , fractionalRadius
  type            (stellarLuminosities         )                          :: luminositiesDisk

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Exponential_Enclosed_Mass=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Get the total mass.
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        Node_Component_Disk_Exponential_Enclosed_Mass=self%massGas()+self%massStellar()
     case (massTypeGaseous)
        Node_Component_Disk_Exponential_Enclosed_Mass=self%massGas()
     case (massTypeStellar)
        Node_Component_Disk_Exponential_Enclosed_Mass=               self%massStellar()
     end select
  case (weightByLuminosity)
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Exponential_Enclosed_Mass   =luminositiesDisk%luminosity(weightIndex)
     end select
  end select
  ! Return if no mass.
  if (Node_Component_Disk_Exponential_Enclosed_Mass <=       0.0d0) return
  ! Return if the total mass was requested.
  if (radius                                        >= radiusLarge) return
  ! Compute the actual mass.
  diskRadius=self%radius()
  if (diskRadius > 0.0d0) then
     fractionalRadius=radius/diskRadius
     Node_Component_Disk_Exponential_Enclosed_Mass=Node_Component_Disk_Exponential_Enclosed_Mass&
          &*Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless(fractionalRadius)
  end if
  return
end function Node_Component_Disk_Exponential_Enclosed_Mass

double precision function Node_Component_Disk_Exponential_Density(self,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
  !% Computes the density at a given position for an exponential disk.
  use Galactic_Structure_Options
  use Numerical_Constants_Math
  use Coordinate_Systems
  use Tables
  use Node_Component_Disk_Exponential_Data
  implicit none
  class           (nodeComponentDiskExponential)           , intent(inout)           :: self
  integer                                                  , intent(in   )           :: componentType                    , massType              , &
       &                                                                                weightBy                         , weightIndex
  double precision                                         , intent(in   )           :: positionSpherical      (3)
  logical                                                  , intent(in   ), optional :: haloLoaded
  double precision                              , parameter                          :: diskHeightToRadiusRatio   =0.1d0
  double precision                              , parameter                          :: coshArgumentMaximum       =50.0d0
  double precision                                                                   :: coshTerm                         , fractionalHeight      , &
       &                                                                                fractionalRadius                 , positionCylindrical(3)
  type            (stellarLuminosities         )                                     :: luminositiesDisk

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Exponential_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Determine mass/luminosity type.
  select case (weightBy)
  case (weightByMass      )
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic)
        Node_Component_Disk_Exponential_Density=self%massGas()+self%massStellar()
     case (massTypeGaseous)
        Node_Component_Disk_Exponential_Density=self%massGas()
     case (massTypeStellar)
        Node_Component_Disk_Exponential_Density=               self%massStellar()
     end select
  case (weightByLuminosity)
     select case (massType)
     case (massTypeAll,massTypeBaryonic,massTypeGalactic,massTypeStellar)
        luminositiesDisk=self%luminositiesStellar()
        Node_Component_Disk_Exponential_Density=luminositiesDisk%luminosity(weightIndex)
     end select
  end select
  ! Skip further calculation if mass or radius is zero.
  if (Node_Component_Disk_Exponential_Density > 0.0d0 .and. self%radius() > 0.0d0) then
     ! Compute the actual density.
     positionCylindrical=Coordinates_Spherical_To_Cylindrical(positionSpherical)
     fractionalRadius=positionCylindrical(1)/                         self%radius()
     fractionalHeight=positionCylindrical(3)/(diskHeightToRadiusRatio*self%radius())
     if (fractionalHeight > coshArgumentMaximum) then
        coshTerm=(2.0d0*exp(-0.5d0*fractionalHeight)/(1.0d0+exp(-fractionalHeight)))**2
     else
        coshTerm=1.0d0/cosh(0.5d0*fractionalHeight)**2
     end if
     Node_Component_Disk_Exponential_Density=Node_Component_Disk_Exponential_Density*exp(-fractionalRadius)*coshTerm/4.0d0&
          &/Pi/self%radius()**3/diskHeightToRadiusRatio
  end if
  return
end function Node_Component_Disk_Exponential_Density

double precision function Node_Component_Disk_Exponential_Potential(self,radius,componentType,massType,haloLoaded)
  !% Compute the gravitational potential due to an exponential disk.
  use Numerical_Constants_Physical
  use Bessel_Functions
  use Galactic_Structure_Options
  use Tables
  use Node_Component_Disk_Exponential_Data
  implicit none
  class           (nodeComponentDiskExponential), intent(inout)           :: self
  integer                                       , intent(in   )           :: componentType               , massType
  double precision                              , intent(in   )           :: radius
  logical                                       , intent(in   ), optional :: haloLoaded
  double precision                              , parameter               :: radiusLimitPointMass=500.0d0
  double precision                                                        :: componentMass               , correctionSmallRadius, &
       &                                                                     halfRadius

  ! Return immediately if disk component is not requested.
  Node_Component_Disk_Exponential_Potential=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Avoid an arithmetic exception at radius zero.
  if (radius <= 0.0d0) return
  ! Get the relevant mass of the disk.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
  if (componentMass <= 0.0d0) return
  ! If the radius is sufficiently large, treat the disk as a point mass.
  if (radius > radiusLimitPointMass*self%radius()) then
     Node_Component_Disk_Exponential_Potential=-gravitationalConstantGalacticus*componentMass/radius
  else
     ! Radius is sufficiently small to use the full calculation.
     ! Compute the potential. If the radius is lower than the height then approximate the disk mass as being spherically distributed.
     if (radius > heightToRadialScaleDisk*self%radius()) then
        halfRadius           =radius/2.0d0/self%radius()
        correctionSmallRadius=0.0d0
     else
        halfRadius           =heightToRadialScaleDisk/2.0d0
        correctionSmallRadius=-Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless(heightToRadialScaleDisk)&
             &/heightToRadialScaleDisk/self%radius()
        if (radius > 0.0d0) correctionSmallRadius=correctionSmallRadius&
             &+Node_Component_Disk_Exponential_Enclosed_Mass_Dimensionless(radius/self%radius())/radius
        correctionSmallRadius=-correctionSmallRadius*gravitationalConstantGalacticus*componentMass
     end if
     ! Compute the potential including the correction to small radii.
     Node_Component_Disk_Exponential_Potential=            &
          &             -0.5d0                             &
          &             *gravitationalConstantGalacticus   &
          &             *componentMass                     &
          &             /self%radius()                     &
          &             *(                                 &
          &                Bessel_Function_I0(halfRadius)  &
          &               *Bessel_Function_K1(halfRadius)  &
          &               -Bessel_Function_I1(halfRadius)  &
          &               *Bessel_Function_K0(halfRadius)  &
          &              )                                 &
          &             + correctionSmallRadius
  end if
  return
end function Node_Component_Disk_Exponential_Potential

double precision function Node_Component_Disk_Exponential_Rotation_Curve(self,radius,componentType,massType,haloLoaded)
  !% Computes the rotation curve at a given radius for an exponential disk.
  use Node_Component_Disk_Exponential_Data
  use Galactic_Structure_Options
  use Numerical_Constants_Physical
  implicit none
  class           (nodeComponentDiskExponential), intent(inout)           :: self
  integer                                       , intent(in   )           :: componentType         , massType
  double precision                              , intent(in   )           :: radius
  logical                                       , intent(in   ), optional :: haloLoaded
  double precision                                                        :: componentMass         , diskRadius, fractionalRadius, &
       &                                                                     fractionalRadiusFactor, halfRadius

  ! Set to zero by default.
  Node_Component_Disk_Exponential_Rotation_Curve=0.0d0
  ! Get the mass of the disk.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
  if (componentMass <= 0.0d0) return
  ! Compute the actual velocity.
  diskRadius=self%radius()
  if (diskRadius > 0.0d0) then
     fractionalRadius=radius/diskRadius
     if (fractionalRadius > fractionalRadiusMaximum) then
        ! Beyond some maximum radius, approximate the disk as a spherical distribution to avoid evaluating Bessel functions for
        ! very large arguments.
        Node_Component_Disk_Exponential_Rotation_Curve=sqrt(gravitationalConstantGalacticus*componentMass/radius)
     else
        ! We are often called at precisely one scale length. Use pre-computed factors in that case.
        if (fractionalRadius == 1.0d0) then
           !$omp critical (ExponentialDiskFactorCompute)
           if (.not.scaleLengthFactorSet) then
              halfRadius          =0.5d0
              scaleLengthFactor   =Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors(halfRadius)
              scaleLengthFactorSet=.true.
           end if
           !$omp end critical (ExponentialDiskFactorCompute)
           fractionalRadiusFactor=scaleLengthFactor
        else
           halfRadius=0.5d0*fractionalRadius
           fractionalRadiusFactor=Node_Component_Disk_Exponential_Rotation_Curve_Bessel_Factors(halfRadius)
        end if
        Node_Component_Disk_Exponential_Rotation_Curve=sqrt(2.0d0*(gravitationalConstantGalacticus*componentMass/diskRadius)&
             &*fractionalRadiusFactor)
     end if
  end if
  return
end function Node_Component_Disk_Exponential_Rotation_Curve

double precision function Node_Component_Disk_Exponential_Rotation_Curve_Gradient(self,radius,componentType,massType,haloLoaded)
  !% Computes the rotation curve gradient for an exponential disk.
  use Tables
  use Node_Component_Disk_Exponential_Data
  use Galactic_Structure_Options
  use Numerical_Constants_Physical
  use Numerical_Constants_Prefixes
  implicit none
  class           (nodeComponentDiskExponential), intent(inout)           :: self
  integer                                       , intent(in   )           :: componentType , massType
  double precision                              , intent(in   )           :: radius
  logical                                       , intent(in   ), optional :: haloLoaded
  double precision                                                        :: besselArgument, besselFactor, componentMass

  ! Set to zero by default.
  Node_Component_Disk_Exponential_Rotation_Curve_Gradient=0.0d0
  ! Return if radius is zero.
  if (radius <= 0.0d0) return
  ! Compute the enclosed mass.
  componentMass=self%enclosedMass(radiusLarge,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
  if (componentMass <= 0.0d0) return
  if (self%radius() <= 0.0d0) return
  besselArgument= radius        &
       &         /2.0d0         &
       &         /self%radius()
  if (2.0d0*besselArgument > fractionalRadiusMaximum) then
     ! Beyond some maximum radius, approximate the disk as a point mass to avoid evaluating Bessel functions for
     ! very large arguments.
     Node_Component_Disk_Exponential_Rotation_Curve_Gradient=-gravitationalConstantGalacticus*componentMass/radius**2
     return
  end if
  ! Checks for low radius and approximations.
  besselFactor=Node_Component_Disk_Exponential_Rttn_Crv_Grdnt_Bssl_Fctrs(besselArgument)
  Node_Component_Disk_Exponential_Rotation_Curve_Gradient=           &
       &                          gravitationalConstantGalacticus    &
       &                         *componentMass                      &
       &                         *besselFactor                       &
       &                         /self%radius()**2
  return
end function Node_Component_Disk_Exponential_Rotation_Curve_Gradient

double precision function Node_Component_Disk_Exponential_Surface_Density(self,positionCylindrical,componentType,massType&
     &,haloLoaded)
  !% Computes the surface density at a given position for an exponential disk.
  use Galactic_Structure_Options
  use Numerical_Constants_Math
  use Tables
  use Node_Component_Disk_Exponential_Data
  implicit none
  class           (nodeComponentDiskExponential), intent(inout)           :: self
  integer                                       , intent(in   )           :: componentType         , massType
  double precision                              , intent(in   )           :: positionCylindrical(3)
  logical                                       , intent(in   ), optional :: haloLoaded
  type            (treeNode                    ), pointer                 :: selfNode
  double precision                                                        :: fractionalRadius

    ! Return immediately if disk component is not requested.
  Node_Component_Disk_Exponential_Surface_Density=0.0d0
  if (.not.(componentType == componentTypeAll .or. componentType == componentTypeDisk)) return
  ! Check whether this is a new node.
  selfNode => self%host()
  if (selfNode%uniqueID() /= lastUniqueID) call Node_Component_Disk_Exponential_Reset(selfNode%uniqueID())
  ! Determine disk radius.
  if (.not.radiusScaleDiskComputed) then
     radiusScaleDisk        =self%radius()
     radiusScaleDiskComputed=.true.
  end if
  ! Return zero if the disk has unphysical size.
  if (radiusScaleDisk <= 0.0d0) then
     Node_Component_Disk_Exponential_Surface_Density=0.0d0
     return
  end if
  ! Determine mass type.
  select case (massType)
  case (massTypeAll,massTypeBaryonic,massTypeGalactic)
     if (.not.surfaceDensityCentralTotalComputed  ) then
        surfaceDensityCentralTotal          =(self%massGas()+self%massStellar())/2.0d0/Pi/radiusScaleDisk**2
        surfaceDensityCentralTotalComputed  =.true.
     end if
     Node_Component_Disk_Exponential_Surface_Density=surfaceDensityCentralTotal
  case (massTypeGaseous)
     if (.not.surfaceDensityCentralGasComputed    ) then
        surfaceDensityCentralGas            = self%massGas()                    /2.0d0/Pi/radiusScaleDisk**2
        surfaceDensityCentralGasComputed    =.true.
     end if
     Node_Component_Disk_Exponential_Surface_Density=surfaceDensityCentralGas
  case (massTypeStellar)
     if (.not.surfaceDensityCentralStellarComputed) then
        surfaceDensityCentralStellar        =                self%massStellar() /2.0d0/Pi/radiusScaleDisk**2
        surfaceDensityCentralStellarComputed=.true.
     end if
     Node_Component_Disk_Exponential_Surface_Density=surfaceDensityCentralStellar
  end select
  ! Return if no density.
  if (Node_Component_Disk_Exponential_Surface_Density <= 0.0d0) return
  ! Compute the actual density.
  fractionalRadius=positionCylindrical(1)/radiusScaleDisk
  Node_Component_Disk_Exponential_Surface_Density=Node_Component_Disk_Exponential_Surface_Density*exp(-fractionalRadius)
  return
end function Node_Component_Disk_Exponential_Surface_Density

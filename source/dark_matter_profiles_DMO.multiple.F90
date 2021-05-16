!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

  !% An implementation of multiple dark matter halo profiles which allow different profiles for the host and the satellite.

  !# <darkMatterProfileDMO name="darkMatterProfileDMOMultiple">
  !#  <description>
  !#   A dark matter profile DMO class in which the density profiles of the host halo and the satellite halo can be set separately
  !#   to any other {\normalfont \ttfamily darkMatterProfileDMO} available.
  !#  </description>
  !# </darkMatterProfileDMO>
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOMultiple
     !% A dark matter halo profile class implementing multiple dark matter halos which allow different profiles for the host and the satellite.
     private
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMOHost_ => null(), darkMatterProfileDMOSatellite_ => null()
   contains
     final     ::                                      multipleDestructor
     procedure :: density                           => multipleDensity
     procedure :: densityLogSlope                   => multipleDensityLogSlope
     procedure :: radiusEnclosingDensity            => multipleRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => multipleRadiusEnclosingMass
     procedure :: radialMoment                      => multipleRadialMoment
     procedure :: enclosedMass                      => multipleEnclosedMass
     procedure :: potential                         => multiplePotential
     procedure :: circularVelocity                  => multipleCircularVelocity
     procedure :: radiusCircularVelocityMaximum     => multipleRadiusCircularVelocityMaximum
     procedure :: circularVelocityMaximum           => multipleCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => multipleRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => multipleRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => multipleRotationNormalization
     procedure :: energy                            => multipleEnergy
     procedure :: energyGrowthRate                  => multipleEnergyGrowthRate
     procedure :: kSpace                            => multipleKSpace
     procedure :: freefallRadius                    => multipleFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => multipleFreefallRadiusIncreaseRate
  end type darkMatterProfileDMOMultiple

  interface darkMatterProfileDMOMultiple
     !% Constructors for the {\normalfont \ttfamily multiple} dark matter halo profile class.
     module procedure multipleConstructorParameters
     module procedure multipleConstructorInternal
  end interface darkMatterProfileDMOMultiple

contains

  function multipleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily multiple} dark matter halo profile class which takes a parameter set as input.
    implicit none
    type   (darkMatterProfileDMOMultiple)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    class  (darkMatterProfileDMOClass   ), pointer       :: darkMatterProfileDMOHost_, darkMatterProfileDMOSatellite_

    !# <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMOHost_"      parameterName="darkMatterProfileDMOHost"      source="parameters"/>
    !# <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMOSatellite_" parameterName="darkMatterProfileDMOSatellite" source="parameters"/>
    self=darkMatterProfileDMOMultiple(darkMatterProfileDMOHost_,darkMatterProfileDMOSatellite_)
    !# <objectDestructor name="darkMatterProfileDMOHost_"     />
    !# <objectDestructor name="darkMatterProfileDMOSatellite_"/>
    return
  end function multipleConstructorParameters

  function multipleConstructorInternal(darkMatterProfileDMOHost_,darkMatterProfileDMOSatellite_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily multiple} dark matter profile class.
    implicit none
    type            (darkMatterProfileDMOMultiple)                        :: self
    class           (darkMatterProfileDMOClass   ), intent(in   ), target :: darkMatterProfileDMOHost_, darkMatterProfileDMOSatellite_
    !# <constructorAssign variables="*darkMatterProfileDMOHost_, *darkMatterProfileDMOSatellite_"/>

    return
  end function multipleConstructorInternal

  subroutine multipleDestructor(self)
    !% Destructor for the {\normalfont \ttfamily multiple} dark matter halo profile class.
    implicit none
    type(darkMatterProfileDMOMultiple), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfileDMOHost_"     />
    !# <objectDestructor name="self%darkMatterProfileDMOSatellite_"/>
    return
  end subroutine multipleDestructor

  double precision function multipleDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    if (node%isSatellite()) then
       multipleDensity=self%darkMatterProfileDMOSatellite_%density(node,radius)
    else
       multipleDensity=self%darkMatterProfileDMOHost_     %density(node,radius)
    end if
    return
  end function multipleDensity

  double precision function multipleDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    if (node%isSatellite()) then
       multipleDensityLogSlope=self%darkMatterProfileDMOSatellite_%densityLogSlope(node,radius)
    else
       multipleDensityLogSlope=self%darkMatterProfileDMOHost_     %densityLogSlope(node,radius)
    end if
    return
  end function multipleDensityLogSlope

  double precision function multipleRadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: density

    if (node%isSatellite()) then
       multipleRadiusEnclosingDensity=self%darkMatterProfileDMOSatellite_%radiusEnclosingDensity(node,density)
    else
       multipleRadiusEnclosingDensity=self%darkMatterProfileDMOHost_     %radiusEnclosingDensity(node,density)
    end if
    return
  end function multipleRadiusEnclosingDensity

  double precision function multipleRadiusEnclosingMass(self,node,mass)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout), target  :: self
    type            (treeNode                    ), intent(inout), target  :: node
    double precision                              , intent(in   )          :: mass

    if (node%isSatellite()) then
       multipleRadiusEnclosingMass=self%darkMatterProfileDMOSatellite_%radiusEnclosingMass(node,mass)
    else
       multipleRadiusEnclosingMass=self%darkMatterProfileDMOHost_     %radiusEnclosingMass(node,mass)
    end if
    return
  end function multipleRadiusEnclosingMass

  double precision function multipleRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (node%isSatellite()) then
       multipleRadialMoment=self%darkMatterProfileDMOSatellite_%radialMoment(node,moment,radiusMinimum,radiusMaximum)
    else
       multipleRadialMoment=self%darkMatterProfileDMOHost_     %radialMoment(node,moment,radiusMinimum,radiusMaximum)
    end if
    return
  end function multipleRadialMoment

  double precision function multipleEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    if (node%isSatellite()) then
       multipleEnclosedMass=self%darkMatterProfileDMOSatellite_%enclosedMass(node,radius)
    else
       multipleEnclosedMass=self%darkMatterProfileDMOHost_     %enclosedMass(node,radius)
    end if
    return
  end function multipleEnclosedMass

  double precision function multiplePotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout)           :: self
    type            (treeNode                    ), intent(inout), target   :: node
    double precision                              , intent(in   )           :: radius
    integer                                       , intent(  out), optional :: status

    if (node%isSatellite()) then
       multiplePotential=self%darkMatterProfileDMOSatellite_%potential(node,radius,status)
    else
       multiplePotential=self%darkMatterProfileDMOHost_     %potential(node,radius,status)
    end if
    return
  end function multiplePotential

  double precision function multipleCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    if (node%isSatellite()) then
       multipleCircularVelocity=self%darkMatterProfileDMOSatellite_%circularVelocity(node,radius)
    else
       multipleCircularVelocity=self%darkMatterProfileDMOHost_     %circularVelocity(node,radius)
    end if
    return
  end function multipleCircularVelocity

  double precision function multipleRadiusCircularVelocityMaximum(self,node)
    !% Returns the radius (in Mpc) at which the maximum circular velocity is achieved in the dark matter profile of {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileDMOMultiple), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    if (node%isSatellite()) then
       multipleRadiusCircularVelocityMaximum=self%darkMatterProfileDMOSatellite_%radiusCircularVelocityMaximum(node)
    else
       multipleRadiusCircularVelocityMaximum=self%darkMatterProfileDMOHost_     %radiusCircularVelocityMaximum(node)
    end if
    return
  end function multipleRadiusCircularVelocityMaximum

  double precision function multipleCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileDMOMultiple), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    if (node%isSatellite()) then
       multipleCircularVelocityMaximum=self%darkMatterProfileDMOSatellite_%circularVelocityMaximum(node)
    else
       multipleCircularVelocityMaximum=self%darkMatterProfileDMOHost_     %circularVelocityMaximum(node)
    end if
    return
  end function multipleCircularVelocityMaximum

  double precision function multipleRadialVelocityDispersion(self,node,radius)
    !% Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    if (node%isSatellite()) then
       multipleRadialVelocityDispersion=self%darkMatterProfileDMOSatellite_%radialVelocityDispersion(node,radius)
    else
       multipleRadialVelocityDispersion=self%darkMatterProfileDMOHost_     %radialVelocityDispersion(node,radius)
    end if
    return
  end function multipleRadialVelocityDispersion

  double precision function multipleRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: specificAngularMomentum

    if (node%isSatellite()) then
       multipleRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMOSatellite_%radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
    else
       multipleRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMOHost_     %radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
    end if
    return
  end function multipleRadiusFromSpecificAngularMomentum

  double precision function multipleRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    implicit none
    class(darkMatterProfileDMOMultiple), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    if (node%isSatellite()) then
       multipleRotationNormalization=self%darkMatterProfileDMOSatellite_%rotationNormalization(node)
    else
       multipleRotationNormalization=self%darkMatterProfileDMOHost_     %rotationNormalization(node)
    end if
    return
  end function multipleRotationNormalization

  double precision function multipleEnergy(self,node)
    !% Return the energy of a multiple halo density profile.
    implicit none
    class(darkMatterProfileDMOMultiple), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    if (node%isSatellite()) then
       multipleEnergy=self%darkMatterProfileDMOSatellite_%energy(node)
    else
       multipleEnergy=self%darkMatterProfileDMOHost_     %energy(node)
    end if
    return
  end function multipleEnergy

  double precision function multipleEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of a multiple halo density profile.
    implicit none
    class(darkMatterProfileDMOMultiple), intent(inout)         :: self
    type (treeNode                    ), intent(inout), target :: node

    if (node%isSatellite()) then
       multipleEnergyGrowthRate=self%darkMatterProfileDMOSatellite_%energyGrowthRate(node)
    else
       multipleEnergyGrowthRate=self%darkMatterProfileDMOHost_     %energyGrowthRate(node)
    end if
    return
  end function multipleEnergyGrowthRate

  double precision function multipleKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the multiple density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: waveNumber

    if (node%isSatellite()) then
       multipleKSpace=self%darkMatterProfileDMOSatellite_%kSpace(node,waveNumber)
    else
       multipleKSpace=self%darkMatterProfileDMOHost_     %kSpace(node,waveNumber)
    end if
    return
  end function multipleKSpace

  double precision function multipleFreefallRadius(self,node,time)
    !% Returns the freefall radius in the multiple density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: time

    if (node%isSatellite()) then
       multipleFreefallRadius=self%darkMatterProfileDMOSatellite_%freefallRadius(node,time)
    else
       multipleFreefallRadius=self%darkMatterProfileDMOHost_     %freefallRadius(node,time)
    end if
    return
  end function multipleFreefallRadius

  double precision function multipleFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the multiple density profile at the specified {\normalfont
    !% \ttfamily time} (given in Gyr).
    implicit none
    class           (darkMatterProfileDMOMultiple), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: time

    if (node%isSatellite()) then
       multipleFreefallRadiusIncreaseRate=self%darkMatterProfileDMOSatellite_%freefallRadiusIncreaseRate(node,time)
    else
       multipleFreefallRadiusIncreaseRate=self%darkMatterProfileDMOHost_     %freefallRadiusIncreaseRate(node,time)
    end if
    return
  end function multipleFreefallRadiusIncreaseRate

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% An implementation of non-dark-matter-only dark matter halo profiles which are unchanged from their dark-matter-only counterpart.

  use :: Cosmology_Parameters    , only : cosmologyParameters , cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass

  !# <darkMatterProfile name="darkMatterProfileDarkMatterOnly">
  !#  <description>An implementation of non-dark-matter-only dark matter halo profiles which are unchanged from their dark-matter-only counterpart.</description>
  !# </darkMatterProfile>
  type, extends(darkMatterProfileClass) :: darkMatterProfileDarkMatterOnly
     !% A dark matter halo profile class implementing non-dark-matter-only dark matter halo profiles which are unchanged from their dark-matter-only counterpart.
     private
     class           (cosmologyParametersClass ), pointer :: cosmologyParameters_  => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     double precision                                     :: darkMatterFraction
   contains
     final     ::                                      darkMatterOnlyDestructor
     procedure :: density                           => darkMatterOnlyDensity
     procedure :: densityLogSlope                   => darkMatterOnlyDensityLogSlope
     procedure :: radiusEnclosingDensity            => darkMatterOnlyRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => darkMatterOnlyRadiusEnclosingMass
     procedure :: radialMoment                      => darkMatterOnlyRadialMoment
     procedure :: enclosedMass                      => darkMatterOnlyEnclosedMass
     procedure :: potential                         => darkMatterOnlyPotential
     procedure :: circularVelocity                  => darkMatterOnlyCircularVelocity
     procedure :: circularVelocityMaximum           => darkMatterOnlyCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => darkMatterOnlyRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => darkMatterOnlyRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => darkMatterOnlyRotationNormalization
     procedure :: energy                            => darkMatterOnlyEnergy
     procedure :: energyGrowthRate                  => darkMatterOnlyEnergyGrowthRate
     procedure :: kSpace                            => darkMatterOnlyKSpace
     procedure :: freefallRadius                    => darkMatterOnlyFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => darkMatterOnlyFreefallRadiusIncreaseRate
  end type darkMatterProfileDarkMatterOnly

  interface darkMatterProfileDarkMatterOnly
     !% Constructors for the {\normalfont \ttfamily darkMatterOnly} non-dark-matter-only dark matter halo profile class.
     module procedure darkMatterOnlyConstructorParameters
     module procedure darkMatterOnlyConstructorInternal
  end interface darkMatterProfileDarkMatterOnly

contains

  function darkMatterOnlyConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily darkMatterOnly} non-dark-matter-only dark matter halo profile class which takes
    !% a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDarkMatterOnly)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class(darkMatterProfileDMOClass      ), pointer       :: darkMatterProfileDMO_
    class(darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_

    !# <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !# <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    self=darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_" />
    !# <objectDestructor name="darkMatterHaloScale_" />
    !# <objectDestructor name="darkMatterProfileDMO_"/>
    return
  end function darkMatterOnlyConstructorParameters

  function darkMatterOnlyConstructorInternal(cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily darkMatterOnly} dark matter profile class.
    implicit none
    type (darkMatterProfileDarkMatterOnly)                        :: self
    class(cosmologyParametersClass       ), intent(in   ), target :: cosmologyParameters_
    class(darkMatterProfileDMOClass      ), intent(in   ), target :: darkMatterProfileDMO_
    class(darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    !# <constructorAssign variables="*cosmologyParameters_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>

    ! Evaluate the dark matter fraction.
    self%darkMatterFraction=+1.0d0                                   &
         &                  -self%cosmologyParameters_%OmegaBaryon() &
         &                  /self%cosmologyParameters_%OmegaMatter()
    return
  end function darkMatterOnlyConstructorInternal

  subroutine darkMatterOnlyDestructor(self)
    !% Destructor for the {\normalfont \ttfamily darkMatterOnly} dark matter halo profile class.
    implicit none
    type(darkMatterProfileDarkMatterOnly), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_" />
    !# <objectDestructor name="self%darkMatterHaloScale_" />
    !# <objectDestructor name="self%darkMatterProfileDMO_"/>
    return
  end subroutine darkMatterOnlyDestructor

  double precision function darkMatterOnlyDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    darkMatterOnlyDensity=+self%darkMatterFraction                         &
         &                *self%darkMatterProfileDMO_%density(node,radius)
    return
  end function darkMatterOnlyDensity

  double precision function darkMatterOnlyDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    darkMatterOnlyDensityLogSlope=self%darkMatterProfileDMO_%densityLogSlope(node,radius)
    return
  end function darkMatterOnlyDensityLogSlope

  double precision function darkMatterOnlyRadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout), target :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: density

    darkMatterOnlyRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity(node,density/self%darkMatterFraction)
    return
  end function darkMatterOnlyRadiusEnclosingDensity

  double precision function darkMatterOnlyRadiusEnclosingMass(self,node,mass)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout), target :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: mass

    darkMatterOnlyRadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass(node,mass/self%darkMatterFraction)
    return
  end function darkMatterOnlyRadiusEnclosingMass

  double precision function darkMatterOnlyRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout)           :: self
    type            (treeNode                       ), intent(inout)           :: node
    double precision                                 , intent(in   )           :: moment
    double precision                                 , intent(in   ), optional :: radiusMinimum, radiusMaximum

    darkMatterOnlyRadialMoment=+self%darkMatterFraction                                                          &
         &                     *self%darkMatterProfileDMO_%radialMoment(node,moment,radiusMinimum,radiusMaximum)
    return
  end function darkMatterOnlyRadialMoment

  double precision function darkMatterOnlyEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    darkMatterOnlyEnclosedMass=+self%darkMatterFraction                              &
         &                     *self%darkMatterProfileDMO_%enclosedMass(node,radius)
    return
  end function darkMatterOnlyEnclosedMass

  double precision function darkMatterOnlyPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout), target   :: self
    type            (treeNode                       ), intent(inout), target   :: node
    double precision                                 , intent(in   )           :: radius
    integer                                          , intent(  out), optional :: status

    darkMatterOnlyPotential=+self%darkMatterFraction                                  &
         &                  *self%darkMatterProfileDMO_%potential(node,radius,status)
    return
  end function darkMatterOnlyPotential

  double precision function darkMatterOnlyCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    darkMatterOnlyCircularVelocity=+sqrt(self%darkMatterFraction                                 ) &
         &                         *     self%darkMatterProfileDMO_%circularVelocity(node,radius)
    return
  end function darkMatterOnlyCircularVelocity

  double precision function darkMatterOnlyCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    darkMatterOnlyCircularVelocityMaximum=+sqrt(self%darkMatterFraction                                 ) &
         &                                *     self%darkMatterProfileDMO_%circularVelocityMaximum(node)
    return
  end function darkMatterOnlyCircularVelocityMaximum

  double precision function darkMatterOnlyRadialVelocityDispersion(self,node,radius)
    !% Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    darkMatterOnlyRadialVelocityDispersion=+sqrt(self%darkMatterFraction                                         ) &
         &                                 *     self%darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    return
  end function darkMatterOnlyRadialVelocityDispersion

  double precision function darkMatterOnlyRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout), target :: self
    type            (treeNode                       ), intent(inout)         :: node
    double precision                                 , intent(in   )         :: specificAngularMomentum

    darkMatterOnlyRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum(node,specificAngularMomentum/sqrt(self%darkMatterFraction))
    return
  end function darkMatterOnlyRadiusFromSpecificAngularMomentum

  double precision function darkMatterOnlyRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    implicit none
    class(darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    darkMatterOnlyRotationNormalization=self%darkMatterProfileDMO_%rotationNormalization(node)
    return
  end function darkMatterOnlyRotationNormalization

  double precision function darkMatterOnlyEnergy(self,node)
    !% Return the energy of the dark matter halo density profile.
    implicit none
    class(darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    darkMatterOnlyEnergy=+self%darkMatterFraction                **2 &
         &               *self%darkMatterProfileDMO_%energy(node)
    return
  end function darkMatterOnlyEnergy

  double precision function darkMatterOnlyEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of the dark matter halo density profile.
    implicit none
    class(darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    darkMatterOnlyEnergyGrowthRate=+self%darkMatterFraction                          **2 &
         &                         *self%darkMatterProfileDMO_%energyGrowthRate(node)
    return
  end function darkMatterOnlyEnergyGrowthRate

  double precision function darkMatterOnlyKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the dark matter halo density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout)         :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: waveNumber

    ! This is normalized by mass, so no need to include the dark matter fraction.
    darkMatterOnlyKSpace=+self%darkMatterProfileDMO_%kSpace(node,waveNumber)
    return
  end function darkMatterOnlyKSpace

  double precision function darkMatterOnlyFreefallRadius(self,node,time)
    !% Returns the freefall radius in the dark matter halo density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: time

    darkMatterOnlyFreefallRadius=self%darkMatterProfileDMO_%freefallRadius(node,time*sqrt(self%darkMatterFraction))
    return
  end function darkMatterOnlyFreefallRadius

  double precision function darkMatterOnlyFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the freefall radius in the dark matter halo density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    implicit none
    class           (darkMatterProfileDarkMatterOnly), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: time

    darkMatterOnlyFreefallRadiusIncreaseRate=+                                                                sqrt(self%darkMatterFraction)  &
         &                                   *self%darkMatterProfileDMO_%freefallRadiusIncreaseRate(node,time*sqrt(self%darkMatterFraction))
    return
  end function darkMatterOnlyFreefallRadiusIncreaseRate

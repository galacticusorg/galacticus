!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% An implementation of truncated dark matter halo profiles.

  use Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass, darkMatterHaloScale

  !# <darkMatterProfile name="darkMatterProfileTruncated">
  !#  <description>truncated dark matter halo profiles.</description>
  !# </darkMatterProfile>
  type, extends(darkMatterProfileClass) :: darkMatterProfileTruncated
     !% A dark matter halo profile class implementing truncated dark matter halos.
     private
     class           (darkMatterProfileClass  ), pointer :: darkMatterProfile_
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
     double precision                                    :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum
     logical                                             :: unimplementedIsFatal
   contains
     final                                             truncatedDestructor
     procedure :: density                           => truncatedDensity
     procedure :: densityLogSlope                   => truncatedDensityLogSlope
     procedure :: radiusEnclosingDensity            => truncatedRadiusEnclosingDensity
     procedure :: radialMoment                      => truncatedRadialMoment
     procedure :: enclosedMass                      => truncatedEnclosedMass
     procedure :: potential                         => truncatedPotential
     procedure :: circularVelocity                  => truncatedCircularVelocity
     procedure :: circularVelocityMaximum           => truncatedCircularVelocityMaximum
     procedure :: radiusFromSpecificAngularMomentum => truncatedRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => truncatedRotationNormalization
     procedure :: energy                            => truncatedEnergy
     procedure :: energyGrowthRate                  => truncatedEnergyGrowthRate
     procedure :: kSpace                            => truncatedKSpace
     procedure :: freefallRadius                    => truncatedFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => truncatedFreefallRadiusIncreaseRate     
  end type darkMatterProfileTruncated

  interface darkMatterProfileTruncated
     !% Constructors for the {\normalfont \ttfamily truncated} dark matter halo profile class.
     module procedure truncatedConstructorParameters
     module procedure truncatedConstructorInternal
  end interface darkMatterProfileTruncated

contains

  function truncatedConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily truncated} dark matter halo profile class which takes a parameter set as input.
    use Galacticus_Error
    use Input_Parameters
    implicit none
    type            (darkMatterProfileTruncated)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (darkMatterProfileClass    ), pointer       :: darkMatterProfile_
    class           (darkMatterHaloScaleClass  ), pointer       :: darkMatterHaloScale_
    logical                                                     :: unimplementedIsFatal
    double precision                                            :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum

    !# <inputParameter>
    !#   <name>unimplementedIsFatal</name>
    !#   <defaultValue>.true.</defaultValue>
    !#   <source>parameters</source>
    !#   <description>If {\normalfont \ttfamily true}, unimplemented features of the truncated dark matter profile cause fatal errors. Otherwise, the solution for an untruncated profile is returned.</description>
    !#   <type>boolean</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radiusFractionalTruncateMinimum</name>
    !#   <defaultValue>2.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>The minimum radius (in units of the virial radius) to begin truncating the density profile.</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radiusFractionalTruncateMaximum</name>
    !#   <defaultValue>4.0d0</defaultValue>
    !#   <source>parameters</source>
    !#   <description>The maximum radius (in units of the virial radius) to finish truncating the density profile.</description>
    !#   <type>float</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterProfile"   name="darkMatterProfile_"   source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=darkMatterProfileTruncated(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,unimplementedIsFatal,darkMatterProfile_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function truncatedConstructorParameters

  function truncatedConstructorInternal(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,unimplementedIsFatal,darkMatterProfile_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily truncated} dark matter profile class.
    implicit none
    type            (darkMatterProfileTruncated)                        :: self
    class           (darkMatterProfileClass    ), intent(in   ), target :: darkMatterProfile_
    class           (darkMatterHaloScaleClass  ), intent(in   ), target :: darkMatterHaloScale_
    double precision                            , intent(in   )         :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum
    logical                                     , intent(in   )         :: unimplementedIsFatal
    !# <constructorAssign variables="radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,unimplementedIsFatal,*darkMatterProfile_,*darkMatterHaloScale_"/>

    return
  end function truncatedConstructorInternal

  subroutine truncatedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily truncated} dark matter halo profile class.
    implicit none
    type(darkMatterProfileTruncated), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfile_"   />
    !# <objectDestructor name="self%darkMatterHaloScale_" />
    return
  end subroutine truncatedDestructor
  
  double precision function truncatedDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius
    double precision                                            :: radiusVirial, x

    radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
    if      (radius <= radiusVirial*self%radiusFractionalTruncateMinimum) then
       truncatedDensity=+self%darkMatterProfile_%density(node,radius)
    else if (radius >= radiusVirial*self%radiusFractionalTruncateMaximum) then
       truncatedDensity=+0.0d0
    else
       x               =+(     radius                         /radiusVirial-self%radiusFractionalTruncateMinimum) &
            &           /(self%radiusFractionalTruncateMaximum             -self%radiusFractionalTruncateMinimum)
       truncatedDensity=+self%darkMatterProfile_%density(node,radius) &
            &           *(                                            &
            &             +1.0d0                                      &
            &             -3.0d0*x**2                                 &
            &             +2.0d0*x**3                                 &
            &           )
    end if
    return
  end function truncatedDensity

  double precision function truncatedDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius

    if (self%unimplementedIsFatal) then
       truncatedDensityLogSlope=0.0d0
       call Galacticus_Error_Report('density logarithmic slope in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedDensityLogSlope=self%darkMatterProfile_%densityLogSlope(node,radius)
    end if
    return
  end function truncatedDensityLogSlope

  double precision function truncatedRadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout), target :: self
    type            (treeNode                  ), intent(inout), target :: node
    double precision                            , intent(in   )         :: density
    
    if (self%unimplementedIsFatal) then
       truncatedRadiusEnclosingDensity=0.0d0
       call Galacticus_Error_Report('radius enclosing density in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedRadiusEnclosingDensity=self%darkMatterProfile_%densityLogSlope(node,density)
    end if
    return
  end function truncatedRadiusEnclosingDensity

  double precision function truncatedRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout)           :: self
    type            (treeNode                  ), intent(inout)           :: node
    double precision                            , intent(in   )           :: moment
    double precision                            , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (self%unimplementedIsFatal) then
       truncatedRadialMoment=0.0d0
       call Galacticus_Error_Report('radial moment in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedRadialMoment=self%darkMatterProfile_%radialMoment(node,moment,radiusMinimum,radiusMaximum)
    end if
    return 
  end function truncatedRadialMoment

  double precision function truncatedEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use Numerical_Integration
    implicit none
    class           (darkMatterProfileTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius
    double precision                                            :: radiusVirial        , radiusMinimum, &
         &                                                         radiusMaximum
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
    if (radius <= radiusVirial*self%radiusFractionalTruncateMinimum) then
       truncatedEnclosedMass=+self%darkMatterProfile_%enclosedMass(node,radius)
    else
       radiusMinimum        =+           radiusVirial*self%radiusFractionalTruncateMinimum
       radiusMaximum        =+min(radius,radiusVirial*self%radiusFractionalTruncateMaximum)
       truncatedEnclosedMass=+Integrate(                                           &
            &                                             +radiusMinimum         , &
            &                                             +radiusMaximum         , &
            &                                              truncatedMassIntegrand, &
            &                                              integrandFunction     , &
            &                                              integrationWorkspace  , &
            &                           toleranceAbsolute=+0.0d0                 , &
            &                           toleranceRelative=+1.0d-9                  &
            &                          )                                           &
            &                 +self%darkMatterProfile_%enclosedMass(node,radiusMinimum)
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end if
    return

  contains

    double precision function truncatedMassIntegrand(radius)
      !% Integrand for mass in truncated dark matter profiles.
      use Numerical_Constants_Math
      implicit none
      double precision, intent(in   ) :: radius

      truncatedMassIntegrand=4.0d0*Pi*radius**2*self%density(node,radius)
      return
    end function truncatedMassIntegrand

  end function truncatedEnclosedMass
  
  double precision function truncatedPotential(self,node,radius,status)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    !% \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout)           :: self
    type            (treeNode                  ), intent(inout), pointer  :: node
    double precision                            , intent(in   )           :: radius
    integer                                     , intent(  out), optional :: status    
 
     if (self%unimplementedIsFatal) then
       truncatedPotential=0.0d0
       call Galacticus_Error_Report('potential in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedPotential=self%darkMatterProfile_%potential(node,radius,status)
    end if
   return
  end function truncatedPotential

  double precision function truncatedCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius

    if (self%unimplementedIsFatal) then
       truncatedCircularVelocity=0.0d0
       call Galacticus_Error_Report('circular velocity in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedCircularVelocity=self%darkMatterProfile_%circularVelocity(node,radius)
    end if
    return
  end function truncatedCircularVelocity

  double precision function truncatedCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTruncated), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedCircularVelocityMaximum=0.0d0
       call Galacticus_Error_Report('circular velocity maximum in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedCircularVelocityMaximum=self%darkMatterProfile_%circularVelocityMaximum(node)
    end if
    return
  end function truncatedCircularVelocityMaximum

  double precision function truncatedRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: node
    double precision                            , intent(in   )          :: specificAngularMomentum

    if (self%unimplementedIsFatal) then
       truncatedRadiusFromSpecificAngularMomentum=0.0d0
       call Galacticus_Error_Report('radius from specific angular momentum in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedRadiusFromSpecificAngularMomentum=self%darkMatterProfile_%radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
    end if
    return
  end function truncatedRadiusFromSpecificAngularMomentum

  double precision function truncatedRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTruncated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedRotationNormalization=0.0d0
       call Galacticus_Error_Report('rotation normalization in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedRotationNormalization=self%darkMatterProfile_%rotationNormalization(node)
    end if
    return
  end function truncatedRotationNormalization

  double precision function truncatedEnergy(self,node)
    !% Return the energy of a truncated halo density profile.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTruncated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedEnergy=0.0d0
       call Galacticus_Error_Report('energy in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedEnergy=self%darkMatterProfile_%energy(node)
    end if
    return
  end function truncatedEnergy

  double precision function truncatedEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of a truncated halo density profile.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileTruncated), intent(inout)         :: self
    type (treeNode                  ), intent(inout), target :: node

    if (self%unimplementedIsFatal) then
       truncatedEnergyGrowthRate=0.0d0
       call Galacticus_Error_Report('energy growth rate in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedEnergyGrowthRate=self%darkMatterProfile_%energyGrowthRate(node)
    end if
    return
  end function truncatedEnergyGrowthRate
  
  double precision function truncatedKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the truncated density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: node
    double precision                            , intent(in   )          :: waveNumber

    if (self%unimplementedIsFatal) then
       truncatedKSpace=0.0d0
       call Galacticus_Error_Report('k-space in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedKSpace=self%darkMatterProfile_%kSpace(node,waveNumber)
    end if
    return
  end function truncatedKSpace

  double precision function truncatedFreefallRadius(self,node,time)
    !% Returns the freefall radius in the truncated density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: time

    if (self%unimplementedIsFatal) then
       truncatedFreefallRadius=0.0d0
       call Galacticus_Error_Report('freefall radius in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedFreefallRadius=self%darkMatterProfile_%freefallRadius(node,time)
    end if
    return
  end function truncatedFreefallRadius

  double precision function truncatedFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the truncated density profile at the specified {\normalfont
    !% \ttfamily time} (given in Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: time

    if (self%unimplementedIsFatal) then
       truncatedFreefallRadiusIncreaseRate=0.0d0
       call Galacticus_Error_Report('freefall radius increase rate in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedFreefallRadiusIncreaseRate=self%darkMatterProfile_%freefallRadiusIncreaseRate(node,time)
    end if
    return
  end function truncatedFreefallRadiusIncreaseRate

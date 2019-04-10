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

  !# <darkMatterProfileDMO name="darkMatterProfileDMOTruncated">
  !#  <description>truncated dark matter halo profiles.</description>
  !# </darkMatterProfileDMO>
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOTruncated
     !% A dark matter halo profile class implementing truncated dark matter halos.
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_              => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_               => null()
     double precision                                     :: radiusFractionalTruncateMinimum             , radiusFractionalTruncateMaximum
     logical                                              :: unimplementedIsFatal
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8          )           :: lastUniqueID
     ! Stored values of computed quantities.
     double precision                                     :: enclosedMassTruncateMinimumPrevious         , enclosedMassTruncateMaximumPrevious, &
          &                                                  enclosingMassRadiusPrevious
   contains
     final                                             truncatedDestructor
     procedure :: calculationReset                  => truncatedCalculationReset
     procedure :: density                           => truncatedDensity
     procedure :: densityLogSlope                   => truncatedDensityLogSlope
     procedure :: radiusEnclosingDensity            => truncatedRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => truncatedRadiusEnclosingMass
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
  end type darkMatterProfileDMOTruncated

  interface darkMatterProfileDMOTruncated
     !% Constructors for the {\normalfont \ttfamily truncated} dark matter halo profile class.
     module procedure truncatedConstructorParameters
     module procedure truncatedConstructorInternal
  end interface darkMatterProfileDMOTruncated

contains

  function truncatedConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily truncated} dark matter halo profile class which takes a parameter set as input.
    use Galacticus_Error
    use Input_Parameters
    implicit none
    type            (darkMatterProfileDMOTruncated)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass    ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    logical                                                        :: unimplementedIsFatal
    double precision                                               :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum

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
    !# <objectBuilder class="darkMatterProfileDMO"   name="darkMatterProfileDMO_"   source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=darkMatterProfileDMOTruncated(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,unimplementedIsFatal,darkMatterProfileDMO_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterProfileDMO_"  />
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function truncatedConstructorParameters

  function truncatedConstructorInternal(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,unimplementedIsFatal,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily truncated} dark matter profile class.
    implicit none
    type            (darkMatterProfileDMOTruncated)                        :: self
    class           (darkMatterProfileDMOClass    ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    double precision                               , intent(in   )         :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum
    logical                                        , intent(in   )         :: unimplementedIsFatal
    !# <constructorAssign variables="radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,unimplementedIsFatal,*darkMatterProfileDMO_,*darkMatterHaloScale_"/>

    self%lastUniqueID=-1_kind_int8
    return
  end function truncatedConstructorInternal

  subroutine truncatedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily truncated} dark matter halo profile class.
    implicit none
    type(darkMatterProfileDMOTruncated), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterProfileDMO_"   />
    !# <objectDestructor name="self%darkMatterHaloScale_" />
    return
  end subroutine truncatedDestructor

  subroutine truncatedCalculationReset(self,node)
    !% Reset the dark matter profile calculation.
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    self%lastUniqueID                       =node%uniqueID()
    self%enclosingMassRadiusPrevious        =-1.0d0
    self%enclosedMassTruncateMinimumPrevious=-1.0d0
    self%enclosedMassTruncateMaximumPrevious=-1.0d0
    call self%darkMatterHaloScale_%calculationReset(node)
    return
  end subroutine truncatedCalculationReset
  
  double precision function truncatedDensity(self,node,radius)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    double precision                                               :: radiusVirial, x

    radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
    if      (radius <= radiusVirial*self%radiusFractionalTruncateMinimum) then
       truncatedDensity=+self%darkMatterProfileDMO_%density(node,radius)
    else if (radius >= radiusVirial*self%radiusFractionalTruncateMaximum) then
       truncatedDensity=+0.0d0
    else
       x               =+(     radius                         /radiusVirial-self%radiusFractionalTruncateMinimum) &
            &           /(self%radiusFractionalTruncateMaximum             -self%radiusFractionalTruncateMinimum)
       truncatedDensity=+self%darkMatterProfileDMO_%density(node,radius) &
            &           *(                                               &
            &             +1.0d0                                         &
            &             -3.0d0*x**2                                    &
            &             +2.0d0*x**3                                    &
            &           )
    end if
    return
  end function truncatedDensity

  double precision function truncatedDensityLogSlope(self,node,radius)
    !% Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius

    if (self%unimplementedIsFatal) then
       truncatedDensityLogSlope=0.0d0
       call Galacticus_Error_Report('density logarithmic slope in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedDensityLogSlope=self%darkMatterProfileDMO_%densityLogSlope(node,radius)
    end if
    return
  end function truncatedDensityLogSlope

  double precision function truncatedRadiusEnclosingDensity(self,node,density)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout), target :: self
    type            (treeNode                     ), intent(inout), target :: node
    double precision                               , intent(in   )         :: density
    
    if (self%unimplementedIsFatal) then
       truncatedRadiusEnclosingDensity=0.0d0
       call Galacticus_Error_Report('radius enclosing density in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity(node,density)
    end if
    return
  end function truncatedRadiusEnclosingDensity

  double precision function truncatedRadiusEnclosingMass(self,node,mass)
    !% Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    !% {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    use Root_Finder
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout), target :: self
    type            (treeNode                     ), intent(inout), target :: node
    double precision                               , intent(in   )         :: mass
    type            (rootFinder                   ), save                  :: finder
    !$omp threadprivate(finder)
    double precision                                                       :: radiusTruncateMinimum, radiusTruncateMaximum, &
         &                                                                    radiusVirial

    if (mass <= 0.0d0) then
       truncatedRadiusEnclosingMass=0.0d0
       return
    end if
    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Get the virial radius.
    radiusVirial         = self%darkMatterHaloScale_%virialRadius(node)
    ! Compute the radii where the truncation starts and ends.
    radiusTruncateMinimum= radiusVirial*self%radiusFractionalTruncateMinimum
    radiusTruncateMaximum= radiusVirial*self%radiusFractionalTruncateMaximum
    ! Compute the enclosed mass within the radii where the truncation starts and ends if required.
    if (self%enclosedMassTruncateMinimumPrevious < 0.0d0) then
       self%enclosedMassTruncateMinimumPrevious=self%enclosedMass(node,radiusTruncateMinimum)
    end if
    if (self%enclosedMassTruncateMaximumPrevious < 0.0d0) then
       self%enclosedMassTruncateMaximumPrevious=self%enclosedMass(node,radiusTruncateMaximum)
    end if
    ! If the given mass is smaller than the enclosed mass within the radius where the truncation starts,
    ! compute the radius from untruncated profile. If the given mass is larger than the enclosed mass
    ! within the radius where the truncation ends, return the maximum truncation radius. Otherwise, solve
    ! the radius numerically.
    if      (mass <= self%enclosedMassTruncateMinimumPrevious) then
       truncatedRadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass(node,mass)
    else if (mass >= self%enclosedMassTruncateMaximumPrevious) then
       truncatedRadiusEnclosingMass=radiusTruncateMaximum
    else
       ! Initialize the root finder.
       if (.not.finder%isInitialized()) then
          call finder%rangeExpand (                                                                 &
               &                       rangeExpandDownward          =0.5d0                        , &
               &                       rangeExpandUpward            =2.0d0                        , &
               &                       rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                       rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                       rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
          call finder%rootFunction(truncatedEnclosedMassRoot                         )
          call finder%tolerance   (toleranceAbsolute=0.0d0  ,toleranceRelative=1.0d-6)
       end if
       if (self%enclosingMassRadiusPrevious < 0.0d0) then
          self%enclosingMassRadiusPrevious = radiusVirial
       end if
       self%enclosingMassRadiusPrevious=finder%find(rootGuess=self%enclosingMassRadiusPrevious)
       truncatedRadiusEnclosingMass    =self%enclosingMassRadiusPrevious
    end if
    return
    contains
      double precision function truncatedEnclosedMassRoot(radius)
        !% Root function used in solving for the radius that encloses a given mass.
        implicit none
        double precision, intent(in) :: radius

        truncatedEnclosedMassRoot=self%enclosedMass(node,radius)-mass
        return
      end function truncatedEnclosedMassRoot
  end function truncatedRadiusEnclosingMass

  double precision function truncatedRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !% Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout)           :: self
    type            (treeNode                  ), intent(inout)           :: node
    double precision                            , intent(in   )           :: moment
    double precision                            , intent(in   ), optional :: radiusMinimum, radiusMaximum

    if (self%unimplementedIsFatal) then
       truncatedRadialMoment=0.0d0
       call Galacticus_Error_Report('radial moment in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedRadialMoment=self%darkMatterProfileDMO_%radialMoment(node,moment,radiusMinimum,radiusMaximum)
    end if
    return 
  end function truncatedRadialMoment

  double precision function truncatedEnclosedMass(self,node,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    !% units of Mpc).
    use FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use Numerical_Integration
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius
    double precision                                            :: radiusVirial        , radiusMinimum, &
         &                                                         radiusMaximum
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace

    radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
    if (radius <= radiusVirial*self%radiusFractionalTruncateMinimum) then
       truncatedEnclosedMass=+self%darkMatterProfileDMO_%enclosedMass(node,radius)
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
            &                 +self%darkMatterProfileDMO_%enclosedMass(node,radiusMinimum)
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
    class           (darkMatterProfileDMOTruncated), intent(inout)           :: self
    type            (treeNode                  ), intent(inout), pointer  :: node
    double precision                            , intent(in   )           :: radius
    integer                                     , intent(  out), optional :: status    
 
     if (self%unimplementedIsFatal) then
       truncatedPotential=0.0d0
       call Galacticus_Error_Report('potential in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedPotential=self%darkMatterProfileDMO_%potential(node,radius,status)
    end if
   return
  end function truncatedPotential

  double precision function truncatedCircularVelocity(self,node,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    !% {\normalfont \ttfamily radius} (given in units of Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: radius

    if (self%unimplementedIsFatal) then
       truncatedCircularVelocity=0.0d0
       call Galacticus_Error_Report('circular velocity in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedCircularVelocity=self%darkMatterProfileDMO_%circularVelocity(node,radius)
    end if
    return
  end function truncatedCircularVelocity

  double precision function truncatedCircularVelocityMaximum(self,node)
    !% Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedCircularVelocityMaximum=0.0d0
       call Galacticus_Error_Report('circular velocity maximum in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum(node)
    end if
    return
  end function truncatedCircularVelocityMaximum

  double precision function truncatedRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: node
    double precision                            , intent(in   )          :: specificAngularMomentum

    if (self%unimplementedIsFatal) then
       truncatedRadiusFromSpecificAngularMomentum=0.0d0
       call Galacticus_Error_Report('radius from specific angular momentum in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
    end if
    return
  end function truncatedRadiusFromSpecificAngularMomentum

  double precision function truncatedRotationNormalization(self,node)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedRotationNormalization=0.0d0
       call Galacticus_Error_Report('rotation normalization in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedRotationNormalization=self%darkMatterProfileDMO_%rotationNormalization(node)
    end if
    return
  end function truncatedRotationNormalization

  double precision function truncatedEnergy(self,node)
    !% Return the energy of a truncated halo density profile.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    if (self%unimplementedIsFatal) then
       truncatedEnergy=0.0d0
       call Galacticus_Error_Report('energy in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedEnergy=self%darkMatterProfileDMO_%energy(node)
    end if
    return
  end function truncatedEnergy

  double precision function truncatedEnergyGrowthRate(self,node)
    !% Return the rate of change of the energy of a truncated halo density profile.
    use Galacticus_Error
    implicit none
    class(darkMatterProfileDMOTruncated), intent(inout)         :: self
    type (treeNode                  ), intent(inout), target :: node

    if (self%unimplementedIsFatal) then
       truncatedEnergyGrowthRate=0.0d0
       call Galacticus_Error_Report('energy growth rate in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedEnergyGrowthRate=self%darkMatterProfileDMO_%energyGrowthRate(node)
    end if
    return
  end function truncatedEnergyGrowthRate
  
  double precision function truncatedKSpace(self,node,waveNumber)
    !% Returns the Fourier transform of the truncated density profile at the specified {\normalfont \ttfamily waveNumber}
    !% (given in Mpc$^{-1}$).
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: node
    double precision                            , intent(in   )          :: waveNumber

    if (self%unimplementedIsFatal) then
       truncatedKSpace=0.0d0
       call Galacticus_Error_Report('k-space in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedKSpace=self%darkMatterProfileDMO_%kSpace(node,waveNumber)
    end if
    return
  end function truncatedKSpace

  double precision function truncatedFreefallRadius(self,node,time)
    !% Returns the freefall radius in the truncated density profile at the specified {\normalfont \ttfamily time} (given in
    !% Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: time

    if (self%unimplementedIsFatal) then
       truncatedFreefallRadius=0.0d0
       call Galacticus_Error_Report('freefall radius in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedFreefallRadius=self%darkMatterProfileDMO_%freefallRadius(node,time)
    end if
    return
  end function truncatedFreefallRadius

  double precision function truncatedFreefallRadiusIncreaseRate(self,node,time)
    !% Returns the rate of increase of the freefall radius in the truncated density profile at the specified {\normalfont
    !% \ttfamily time} (given in Gyr).
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    class           (darkMatterProfileDMOTruncated), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    double precision                            , intent(in   ) :: time

    if (self%unimplementedIsFatal) then
       truncatedFreefallRadiusIncreaseRate=0.0d0
       call Galacticus_Error_Report('freefall radius increase rate in truncated dark matter profiles is not supported'//{introspection:location})
    else
       truncatedFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate(node,time)
    end if
    return
  end function truncatedFreefallRadiusIncreaseRate

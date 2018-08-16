!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% An implementation of dark matter halo scales based on virial density contrast.

  use Kind_Numbers
  use Tables
  use Cosmology_Parameters
  use Cosmology_Functions
  use Virial_Density_Contrast

  !# <darkMatterHaloScale name="darkMatterHaloScaleVirialDensityContrastDefinition">
  !#  <description>Dark matter halo scales derived from virial density contrasts.</description>
  !# </darkMatterHaloScale>
  type, extends(darkMatterHaloScaleClass) :: darkMatterHaloScaleVirialDensityContrastDefinition
     !% A dark matter halo scale contrast class using virial density contrasts.
     private
     class           (cosmologyParametersClass  ), pointer   :: cosmologyParameters_       => null()
     class           (cosmologyFunctionsClass   ), pointer   :: cosmologyFunctions_        => null()
     class           (virialDensityContrastClass), pointer   :: virialDensityContrast_     => null()
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8            )            :: lastUniqueID
     ! Record of whether or not halo scales have already been computed for this node.
     logical                                                 :: dynamicalTimescaleComputed          , virialRadiusComputed            , &
          &                                                     virialTemperatureComputed           , virialVelocityComputed
     ! Stored values of halo scales.
     double precision                                        :: dynamicalTimescaleStored            , virialRadiusStored              , &
          &                                                     virialTemperatureStored             , virialVelocityStored            , &
          &                                                     timePrevious                        , densityGrowthRatePrevious       , &
          &                                                     massPrevious
     ! Table for fast lookup of the mean density of halos.
     double precision                                        :: meanDensityTimeMaximum              , meanDensityTimeMinimum   =-1.0d0
     type            (table1DLogarithmicLinear  )            :: meanDensityTable
   contains
     final     ::                                        virialDensityContrastDefinitionDestructor
     procedure :: dynamicalTimescale                  => virialDensityContrastDefinitionDynamicalTimescale
     procedure :: virialVelocity                      => virialDensityContrastDefinitionVirialVelocity
     procedure :: virialVelocityGrowthRate            => virialDensityContrastDefinitionVirialVelocityGrowthRate
     procedure :: virialTemperature                   => virialDensityContrastDefinitionVirialTemperature
     procedure :: virialRadius                        => virialDensityContrastDefinitionVirialRadius
     procedure :: virialRadiusGradientLogarithmicMass => virialDensityContrastDefinitionVirialRadiusGradientLogMass
     procedure :: virialRadiusGrowthRate              => virialDensityContrastDefinitionVirialRadiusGrowthRate
     procedure :: meanDensity                         => virialDensityContrastDefinitionMeanDensity
     procedure :: meanDensityGrowthRate               => virialDensityContrastDefinitionMeanDensityGrowthRate
     procedure :: calculationReset                    => virialDensityContrastDefinitionCalculationReset
  end type darkMatterHaloScaleVirialDensityContrastDefinition

  interface darkMatterHaloScaleVirialDensityContrastDefinition
     !% Constructors for the {\normalfont \ttfamily virialDensityContrastDefinition} dark matter halo scales class.
     module procedure virialDensityContrastDefinitionParameters
     module procedure virialDensityContrastDefinitionInternal
  end interface darkMatterHaloScaleVirialDensityContrastDefinition

  integer, parameter :: virialDensityContrastDefinitionMeanDensityTablePointsPerDecade=100

contains

  function virialDensityContrastDefinitionParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily virialDensityContrastDefinition} dark matter halo scales class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type (darkMatterHaloScaleVirialDensityContrastDefinition), target        :: self
    type (inputParameters                                   ), intent(inout) :: parameters
    class(cosmologyParametersClass                          ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass                           ), pointer       :: cosmologyFunctions_
    class(virialDensityContrastClass                        ), pointer       :: virialDensityContrast_

    !# <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    !# <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    self=darkMatterHaloScaleVirialDensityContrastDefinition(cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function virialDensityContrastDefinitionParameters

  function virialDensityContrastDefinitionInternal(cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_) result(self)
    !% Default constructor for the {\normalfont \ttfamily virialDensityContrastDefinition} dark matter halo scales class.
    implicit none
    type (darkMatterHaloScaleVirialDensityContrastDefinition)                        :: self
    class(cosmologyParametersClass                          ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass                           ), intent(in   ), target :: cosmologyFunctions_
    class(virialDensityContrastClass                        ), intent(in   ), target :: virialDensityContrast_
    !# <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *virialDensityContrast_"/>

    self%lastUniqueID              =-1_kind_int8
    self%dynamicalTimescaleComputed=.false.
    self%virialRadiusComputed      =.false.
    self%virialTemperatureComputed =.false.
    self%virialVelocityComputed    =.false.
    self%meanDensityTimeMaximum    =-1.0d0
    self%meanDensityTimeMinimum    =-1.0d0
    self%timePrevious              =-1.0d0
    self%massPrevious              =-1.0d0
    return
  end function virialDensityContrastDefinitionInternal

  subroutine virialDensityContrastDefinitionDestructor(self)
    !% Destructir for the {\normalfont \ttfamily virialDensityContrastDefinition} dark matter halo scales class.
    implicit none
    type (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"  />
    !# <objectDestructor name="self%cosmologyFunctions_"   />
    !# <objectDestructor name="self%virialDensityContrast_"/>
    if (self%meanDensityTimeMinimum >= 0.0d0) call self%meanDensityTable%destroy()
    return
  end subroutine virialDensityContrastDefinitionDestructor

  subroutine virialDensityContrastDefinitionCalculationReset(self,node)
    !% Reset the halo scales calculation.
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node

    self%virialRadiusComputed      =.false.
    self%virialTemperatureComputed =.false.
    self%virialVelocityComputed    =.false.
    self%dynamicalTimescaleComputed=.false.
    self%lastUniqueID              =node%uniqueID()
    return
  end subroutine virialDensityContrastDefinitionCalculationReset

  double precision function virialDensityContrastDefinitionDynamicalTimescale(self,node)
    !% Returns the dynamical timescale for {\normalfont \ttfamily node}.
    use Numerical_Constants_Astronomical
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if halo dynamical timescale is already computed. Compute and store if not.
    if (.not.self%dynamicalTimescaleComputed) then
       self%dynamicalTimescaleComputed= .true.
       self%dynamicalTimescaleStored  = self%virialRadius  (node) &
            &                          /self%virialVelocity(node) &
            &                          *(megaParsec/kilo/gigaYear)
     end if
    ! Return the stored timescale.
    virialDensityContrastDefinitionDynamicalTimescale=self%dynamicalTimescaleStored
    return
  end function virialDensityContrastDefinitionDynamicalTimescale

  double precision function virialDensityContrastDefinitionVirialVelocity(self,node)
    !% Returns the virial velocity scale for {\normalfont \ttfamily node}.
    use Numerical_Constants_Physical
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic                                ), pointer       :: basic

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if virial velocity is already computed. Compute and store if not.
    if (.not.self%virialVelocityComputed) then
       ! Get the basic component.
       basic => node%basic()
       ! Compute the virial velocity.
       self%virialVelocityStored=sqrt(gravitationalConstantGalacticus*basic%mass() &
            &/self%virialRadius(node))
       ! Record that virial velocity has now been computed.
       self%virialVelocityComputed=.true.
    end if
    ! Return the stored virial velocity.
    virialDensityContrastDefinitionVirialVelocity=self%virialVelocityStored
    return
  end function virialDensityContrastDefinitionVirialVelocity

  double precision function virialDensityContrastDefinitionVirialVelocityGrowthRate(self,node)
    !% Returns the growth rate of the virial velocity scale for {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: node
    class(nodeComponentBasic                                )               , pointer :: basic

    ! Get the basic component.
    basic => node%basic()
    virialDensityContrastDefinitionVirialVelocityGrowthRate= &
         & +0.5d0                                            &
         & *  self%virialVelocity        (node)              &
         & *(                                                &
         &    basic%accretionRate        (    )              &
         &   /basic%mass                 (    )              &
         &   -self%virialRadiusGrowthRate(node)              &
         &   /self%virialRadius          (node)              &
         &  )
    return
  end function virialDensityContrastDefinitionVirialVelocityGrowthRate

  double precision function virialDensityContrastDefinitionVirialTemperature(self,node)
    !% Returns the virial temperature (in Kelvin) for {\normalfont \ttfamily node}.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if virial temperature is already computed. Compute and store if not.
    if (.not.self%virialTemperatureComputed) then
       self%virialTemperatureComputed=.true.
       self%virialTemperatureStored=0.5d0*atomicMassUnit*meanAtomicMassPrimordial*((kilo&
            &*self%virialVelocity(node))**2)/boltzmannsConstant
    end if
    ! Return the stored temperature.
    virialDensityContrastDefinitionVirialTemperature=self%virialTemperatureStored
    return
  end function virialDensityContrastDefinitionVirialTemperature

  double precision function virialDensityContrastDefinitionVirialRadius(self,node)
    !% Returns the virial radius scale for {\normalfont \ttfamily node}.
    use Numerical_Constants_Math
    use Math_Exponentiation
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic                                ), pointer       :: basic

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
    ! Check if virial radius is already computed. Compute and store if not.
    if (.not.self%virialRadiusComputed) then
       ! Get the basic component.
       basic => node%basic()
       ! Compute the virial radius.
       self%virialRadiusStored=cubeRoot(3.0d0*basic%mass()/4.0d0/Pi/self%meanDensity(node))
       ! Record that the virial radius has been computed.
       self%virialRadiusComputed=.true.
    end if
    ! Return the stored value.
    virialDensityContrastDefinitionVirialRadius=self%virialRadiusStored
    return
  end function virialDensityContrastDefinitionVirialRadius

  double precision function virialDensityContrastDefinitionVirialRadiusGradientLogMass(self,node)
    !% Returns the logarithmic gradient of virial radius with halo mass at fixed epoch for {\normalfont \ttfamily node}.
    use Numerical_Constants_Math
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: node
    !GCC$ attributes unused :: self, node
    
    ! Halos at given epoch have fixed density, so radius always grows as the cube-root of mass.
    virialDensityContrastDefinitionVirialRadiusGradientLogMass=1.0d0/3.0d0
    return
  end function virialDensityContrastDefinitionVirialRadiusGradientLogMass

  double precision function virialDensityContrastDefinitionVirialRadiusGrowthRate(self,node)
    !% Returns the growth rate of the virial radius scale for {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: node
    class(nodeComponentBasic)                                               , pointer :: basic

    ! Get the basic component.
    basic => node%basic()
    virialDensityContrastDefinitionVirialRadiusGrowthRate=(1.0d0/3.0d0)*self%virialRadius(node)&
         &*(basic%accretionRate()/basic%mass()-self%meanDensityGrowthRate(node)&
         &/virialDensityContrastDefinitionMeanDensity(self,node))
    return
  end function virialDensityContrastDefinitionVirialRadiusGrowthRate

  double precision function virialDensityContrastDefinitionMeanDensity(self,node)
    !% Returns the mean density for {\normalfont \ttfamily node}.
    implicit none
    class           (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type            (treeNode                                          ), intent(inout) :: node
    class           (nodeComponentBasic                                ), pointer       :: basic
    integer                                                                             :: i                   , meanDensityTablePoints
    double precision                                                                    :: time

    ! Get the basic component.
    basic => node%basic()
    ! Get the time at which this halo was last an isolated halo.
    time=basic%timeLastIsolated()
    if (time <= 0.0d0) time=basic%time()
    ! For mass-dependent virial density contrasts we must always recompute the result.
    if (self%virialDensityContrast_%isMassDependent()) then
       ! Get default objects.
       virialDensityContrastDefinitionMeanDensity =  +self%virialDensityContrast_%densityContrast(basic%mass(),time)    &
            &                                        *self%cosmologyParameters_  %OmegaMatter    (                 )    &
            &                                        *self%cosmologyParameters_  %densityCritical(                 )    &
            &                                        /self%cosmologyFunctions_   %expansionFactor(             time)**3
    else
       ! For non-mass-dependent virial density contrasts we can tabulate as a function of time.
       ! Retabulate the mean density vs. time if necessary.
       if (time < self%meanDensityTimeMinimum .or. time > self%meanDensityTimeMaximum) then
          if (self%meanDensityTimeMinimum <= 0.0d0) then
             self%meanDensityTimeMinimum=                                time/10.0d0
             self%meanDensityTimeMaximum=                                time* 2.0d0
          else
             self%meanDensityTimeMinimum=min(self%meanDensityTimeMinimum,time/10.0d0)
             self%meanDensityTimeMaximum=max(self%meanDensityTimeMaximum,time* 2.0d0)
          end if
          meanDensityTablePoints=int(log10(self%meanDensityTimeMaximum/self%meanDensityTimeMinimum)*dble(virialDensityContrastDefinitionMeanDensityTablePointsPerDecade))+1
          call self%meanDensityTable%destroy()
          call self%meanDensityTable%create(self%meanDensityTimeMinimum,self%meanDensityTimeMaximum,meanDensityTablePoints)
          do i=1,meanDensityTablePoints
             call self%meanDensityTable%populate                                                               &
                  & (                                                                                          &
                  &  +self%virialDensityContrast_%densityContrast(basic%mass(),self%meanDensityTable%x(i))     &
                  &  *self%cosmologyParameters_  %OmegaMatter    (                                       )     &
                  &  *self%cosmologyParameters_  %densityCritical(                                       )     &
                  &  /self%cosmologyFunctions_   %expansionFactor(             self%meanDensityTable%x(i))**3, &
                  &  i                                                                                         &
                  & )
          end do
       end if
       ! Return the stored value.
       virialDensityContrastDefinitionMeanDensity=self%meanDensityTable%interpolate(time)
    end if
    return
  end function virialDensityContrastDefinitionMeanDensity

  double precision function virialDensityContrastDefinitionMeanDensityGrowthRate(self,node)
    !% Returns the growth rate of the mean density for {\normalfont \ttfamily node}.
    implicit none
    class           (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type            (treeNode                                          ), intent(inout), pointer :: node
    class           (nodeComponentBasic                                )               , pointer :: basic
    double precision                                                                             :: expansionFactor    , time

    if (node%isSatellite()) then
       ! Satellite halo is not growing, return zero rate.
       virialDensityContrastDefinitionMeanDensityGrowthRate=0.0d0
    else
       ! Get the basic component.
       basic => node%basic()
       ! Get the time at which this halo was last an isolated halo.
       time=basic%timeLastIsolated()
       ! Check if the time is different from that one previously used.
       if (time /= self%timePrevious .or. basic%mass() /= self%massPrevious) then
          ! It is not, so recompute the density growth rate.
          self%timePrevious=time
          self%massPrevious=basic%mass()
          ! Get the expansion factor at this time.
          expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
          ! Compute growth rate of its mean density based on mean cosmological density and overdensity of a collapsing halo.
          self%densityGrowthRatePrevious=                                                                 &
               & self%meanDensity(node)                                                                   &
               & *(                                                                                       &
               &   +self%virialDensityContrast_%densityContrastRateOfChange(basic%mass(),time           ) &
               &   /self%virialDensityContrast_%densityContrast            (basic%mass(),time           ) &
               &   -3.0d0                                                                                 &
               &   *self%cosmologyFunctions_   %expansionRate              (             expansionFactor) &
               &  )
       end if
       ! Return the stored value.
       virialDensityContrastDefinitionMeanDensityGrowthRate=self%densityGrowthRatePrevious
    end if
    return
  end function virialDensityContrastDefinitionMeanDensityGrowthRate

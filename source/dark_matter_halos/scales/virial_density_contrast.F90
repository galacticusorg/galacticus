!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{RST
  An implementation of dark matter halo scales based on virial density contrast.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Kind_Numbers           , only : kind_int8
  use :: Tables                 , only : table1DLogarithmicLinear
  use :: Virial_Density_Contrast, only : virialDensityContrastClass
  
  !![
  <darkMatterHaloScale name="darkMatterHaloScaleVirialDensityContrastDefinition" recursive="yes" docformat="rst">
   <description>
   Dark matter halo scales derived from virial density contrasts.
   </description>
  </darkMatterHaloScale>
  !!]
  type, extends(darkMatterHaloScaleClass) :: darkMatterHaloScaleVirialDensityContrastDefinition
     !!{RST
     A dark matter halo scale contrast class using virial density contrasts.
     !!}
     private
     class           (cosmologyParametersClass  ), pointer :: cosmologyParameters_       => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_        => null()
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_     => null()
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8            )          :: lastUniqueID
     ! Record of whether or not halo scales have already been computed for this node.
     logical                                               :: timescaleDynamicalComputed          , radiusVirialComputed            , &
          &                                                   temperatureVirialComputed           , velocityVirialComputed
     ! Stored values of halo scales.
     double precision                                      :: timescaleDynamicalStored            , radiusVirialStored              , &
          &                                                   temperatureVirialStored             , velocityVirialStored            , &
          &                                                   timePrevious                        , densityGrowthRatePrevious       , &
          &                                                   massPrevious
     ! Table for fast lookup of the mean density of halos.
     double precision                                      :: densityMeanTimeMaximum              , densityMeanTimeMinimum   =-1.0d0
     type            (table1DLogarithmicLinear  )          :: densityMeanTable
   contains
     !![
     <methods docformat="rst">
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                                        virialDensityContrastDefinitionDestructor
     procedure :: autoHook                            => virialDensityContrastDefinitionAutoHook
     procedure :: calculationReset                    => virialDensityContrastDefinitionCalculationReset
     procedure :: timescaleDynamical                  => virialDensityContrastDefinitionDynamicalTimescale
     procedure :: velocityVirial                      => virialDensityContrastDefinitionVirialVelocity
     procedure :: velocityVirialGrowthRate            => virialDensityContrastDefinitionVirialVelocityGrowthRate
     procedure :: temperatureVirial                   => virialDensityContrastDefinitionVirialTemperature
     procedure :: radiusVirial                        => virialDensityContrastDefinitionVirialRadius
     procedure :: radiusVirialGradientLogarithmicMass => virialDensityContrastDefinitionVirialRadiusGradientLogMass
     procedure :: radiusVirialGrowthRate              => virialDensityContrastDefinitionVirialRadiusGrowthRate
     procedure :: densityMean                         => virialDensityContrastDefinitionMeanDensity
     procedure :: densityMeanGrowthRate               => virialDensityContrastDefinitionMeanDensityGrowthRate
  end type darkMatterHaloScaleVirialDensityContrastDefinition

  interface darkMatterHaloScaleVirialDensityContrastDefinition
     !!{RST
     Constructors for the :galacticus-class:`darkMatterHaloScaleVirialDensityContrastDefinition` dark matter halo scales class.
     !!}
     module procedure virialDensityContrastDefinitionConstructorParameters
     module procedure virialDensityContrastDefinitionConstructorInternal
  end interface darkMatterHaloScaleVirialDensityContrastDefinition

  integer, parameter :: meanDensityTablePointsPerDecade=100

contains

  recursive function virialDensityContrastDefinitionConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterHaloScaleVirialDensityContrastDefinition` dark matter halo scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (darkMatterHaloScaleVirialDensityContrastDefinition), target        :: self
    type   (inputParameters                                   ), intent(inout) :: parameters
    class  (cosmologyParametersClass                          ), pointer       :: cosmologyParameters_
    class  (cosmologyFunctionsClass                           ), pointer       :: cosmologyFunctions_
    class  (virialDensityContrastClass                        ), pointer       :: virialDensityContrast_

    !![
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=darkMatterHaloScaleVirialDensityContrastDefinition(cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function virialDensityContrastDefinitionConstructorParameters

  recursive function virialDensityContrastDefinitionConstructorInternal(cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_) result(self)
    !!{RST
    Default constructor for the ``virialDensityContrastDefinition`` dark matter halo scales class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type   (darkMatterHaloScaleVirialDensityContrastDefinition)                        :: self
    class  (cosmologyParametersClass                          ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass                           ), intent(in   ), target :: cosmologyFunctions_
    class  (virialDensityContrastClass                        ), intent(in   ), target :: virialDensityContrast_
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *virialDensityContrast_"/>
    !!]

    self%lastUniqueID              =-1_kind_int8
    self%timescaleDynamicalComputed=.false.
    self%radiusVirialComputed      =.false.
    self%temperatureVirialComputed =.false.
    self%velocityVirialComputed    =.false.
    self%densityMeanTimeMaximum    =-1.0d0
    self%densityMeanTimeMinimum    =-1.0d0
    self%timePrevious              =-1.0d0
    self%massPrevious              =-1.0d0
    return
  end function virialDensityContrastDefinitionConstructorInternal

  subroutine virialDensityContrastDefinitionAutoHook(self)
    !!{RST
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self

    call calculationResetEvent%attach(self,virialDensityContrastDefinitionCalculationReset,openMPThreadBindingAllLevels,label='darkMatterHaloScaleVirialDensityContrastDefinition')
    return
  end subroutine virialDensityContrastDefinitionAutoHook

  subroutine virialDensityContrastDefinitionDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`darkMatterHaloScaleVirialDensityContrastDefinition` dark matter halo scales class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    if (self%densityMeanTimeMinimum >= 0.0d0) call self%densityMeanTable%destroy()
    if (calculationResetEvent%isAttached(self,virialDensityContrastDefinitionCalculationReset)) call calculationResetEvent%detach(self,virialDensityContrastDefinitionCalculationReset)
    return
  end subroutine virialDensityContrastDefinitionDestructor

  subroutine virialDensityContrastDefinitionCalculationReset(self,node,uniqueID)
    !!{RST
    Reset the halo scales calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type   (treeNode                                          ), intent(inout) :: node
    integer(kind_int8                                         ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%radiusVirialComputed      =.false.
    self%temperatureVirialComputed =.false.
    self%velocityVirialComputed    =.false.
    self%timescaleDynamicalComputed=.false.
    self%lastUniqueID              =uniqueID
    return
  end subroutine virialDensityContrastDefinitionCalculationReset

  double precision function virialDensityContrastDefinitionDynamicalTimescale(self,node)
    !!{RST
    Returns the dynamical timescale for ``node``.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if halo dynamical timescale is already computed. Compute and store if not.
    if (.not.self%timescaleDynamicalComputed) then
       self%timescaleDynamicalComputed= .true.
       self%timescaleDynamicalStored  = self%radiusVirial  (node) &
            &                          /self%velocityVirial(node) &
            &                          *(megaParsec/kilo/gigaYear)
     end if
    ! Return the stored timescale.
    virialDensityContrastDefinitionDynamicalTimescale=self%timescaleDynamicalStored
    return
  end function virialDensityContrastDefinitionDynamicalTimescale

  double precision function virialDensityContrastDefinitionVirialVelocity(self,node)
    !!{RST
    Returns the virial velocity scale for ``node``.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic                                ), pointer       :: basic

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if virial velocity is already computed. Compute and store if not.
    if (.not.self%velocityVirialComputed) then
       ! Get the basic component.
       basic => node%basic()
       ! Compute the virial velocity.
       self%velocityVirialStored=sqrt(gravitationalConstant_internal*basic%mass() &
            &/self%radiusVirial(node))
       ! Record that virial velocity has now been computed.
       self%velocityVirialComputed=.true.
    end if
    ! Return the stored virial velocity.
    virialDensityContrastDefinitionVirialVelocity=self%velocityVirialStored
    return
  end function virialDensityContrastDefinitionVirialVelocity

  double precision function virialDensityContrastDefinitionVirialVelocityGrowthRate(self,node)
    !!{RST
    Returns the growth rate of the virial velocity scale for ``node``.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic                                ), pointer       :: basic

    ! Get the basic component.
    basic => node%basic()
    virialDensityContrastDefinitionVirialVelocityGrowthRate= &
         & +0.5d0                                            &
         & *  self%velocityVirial        (node)              &
         & *(                                                &
         &    basic%accretionRate        (    )              &
         &   /basic%mass                 (    )              &
         &   -self%radiusVirialGrowthRate(node)              &
         &   /self%radiusVirial          (node)              &
         &  )
    return
  end function virialDensityContrastDefinitionVirialVelocityGrowthRate

  double precision function virialDensityContrastDefinitionVirialTemperature(self,node)
    !!{RST
    Returns the virial temperature (in Kelvin) for ``node``.
    !!}
    use :: Numerical_Constants_Astronomical, only : meanAtomicMassPrimordial
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if virial temperature is already computed. Compute and store if not.
    if (.not.self%temperatureVirialComputed) then
       self%temperatureVirialComputed=.true.
       self%temperatureVirialStored=0.5d0*atomicMassUnit*meanAtomicMassPrimordial*((kilo&
            &*self%velocityVirial(node))**2)/boltzmannsConstant
    end if
    ! Return the stored temperature.
    virialDensityContrastDefinitionVirialTemperature=self%temperatureVirialStored
    return
  end function virialDensityContrastDefinitionVirialTemperature

  double precision function virialDensityContrastDefinitionVirialRadius(self,node)
    !!{RST
    Returns the virial radius scale for ``node``.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Math_Exponentiation     , only : cubeRoot
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic                                ), pointer       :: basic

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if virial radius is already computed. Compute and store if not.
    if (.not.self%radiusVirialComputed) then
       ! Get the basic component.
       basic => node%basic()
       ! Compute the virial radius.
       self%radiusVirialStored=cubeRoot(3.0d0*basic%mass()/4.0d0/Pi/self%densityMean(node))
       ! Record that the virial radius has been computed.
       self%radiusVirialComputed=.true.
    end if
    ! Return the stored value.
    virialDensityContrastDefinitionVirialRadius=self%radiusVirialStored
    return
  end function virialDensityContrastDefinitionVirialRadius

  double precision function virialDensityContrastDefinitionVirialRadiusGradientLogMass(self,node)
    !!{RST
    Returns the logarithmic gradient of virial radius with halo mass at fixed epoch for ``node``.
    !!}
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Halos at given epoch have fixed density, so radius always grows as the cube-root of mass.
    virialDensityContrastDefinitionVirialRadiusGradientLogMass=1.0d0/3.0d0
    return
  end function virialDensityContrastDefinitionVirialRadiusGradientLogMass

  double precision function virialDensityContrastDefinitionVirialRadiusGrowthRate(self,node)
    !!{RST
    Returns the growth rate of the virial radius scale for ``node``.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic)                                , pointer       :: basic

    ! Get the basic component.
    basic => node%basic()
    virialDensityContrastDefinitionVirialRadiusGrowthRate=(1.0d0/3.0d0)*self%radiusVirial(node)&
         &*(basic%accretionRate()/basic%mass()-self%densityMeanGrowthRate(node)&
         &/virialDensityContrastDefinitionMeanDensity(self,node))
    return
  end function virialDensityContrastDefinitionVirialRadiusGrowthRate

  double precision function virialDensityContrastDefinitionMeanDensity(self,node)
    !!{RST
    Returns the mean density for ``node``.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type            (treeNode                                          ), intent(inout) :: node
    class           (nodeComponentBasic                                ), pointer       :: basic
    integer                                                                             :: i    , densityMeanTablePoints
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
       if (time < self%densityMeanTimeMinimum .or. time > self%densityMeanTimeMaximum) then
          if (self%densityMeanTimeMinimum <= 0.0d0) then
             self%densityMeanTimeMinimum=                                time/2.0d0
             self%densityMeanTimeMaximum=                                time*2.0d0
          else
             self%densityMeanTimeMinimum=min(self%densityMeanTimeMinimum,time/2.0d0)
             self%densityMeanTimeMaximum=max(self%densityMeanTimeMaximum,time*2.0d0)
          end if
          densityMeanTablePoints=int(log10(self%densityMeanTimeMaximum/self%densityMeanTimeMinimum)*dble(meanDensityTablePointsPerDecade))+1
          call self%densityMeanTable%destroy()
          call self%densityMeanTable%create(self%densityMeanTimeMinimum,self%densityMeanTimeMaximum,densityMeanTablePoints)
          do i=1,densityMeanTablePoints
             call self%densityMeanTable%populate                                                               &
                  & (                                                                                          &
                  &  +self%virialDensityContrast_%densityContrast(basic%mass(),self%densityMeanTable%x(i))     &
                  &  *self%cosmologyParameters_  %OmegaMatter    (                                       )     &
                  &  *self%cosmologyParameters_  %densityCritical(                                       )     &
                  &  /self%cosmologyFunctions_   %expansionFactor(             self%densityMeanTable%x(i))**3, &
                  &  i                                                                                         &
                  & )
          end do
       end if
       ! Return the stored value.
       virialDensityContrastDefinitionMeanDensity=self%densityMeanTable%interpolate(time)
    end if
    return
  end function virialDensityContrastDefinitionMeanDensity

  double precision function virialDensityContrastDefinitionMeanDensityGrowthRate(self,node)
    !!{RST
    Returns the growth rate of the mean density for ``node``.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type            (treeNode                                          ), intent(inout) :: node
    class           (nodeComponentBasic                                ), pointer       :: basic
    double precision                                                                    :: expansionFactor, time

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
               & self%densityMean(node)                                                                   &
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

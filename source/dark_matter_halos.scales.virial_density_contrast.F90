!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  An implementation of dark matter halo scales based on virial density contrast.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Kind_Numbers           , only : kind_int8
  use :: Tables                 , only : table1DLogarithmicLinear
  use :: Virial_Density_Contrast, only : virialDensityContrastClass
  
  !![
  <darkMatterHaloScale name="darkMatterHaloScaleVirialDensityContrastDefinition" recursive="yes">
   <description>Dark matter halo scales derived from virial density contrasts.</description>
   <deepCopy>
    <ignore variables="recursiveSelf"/>
   </deepCopy>
  </darkMatterHaloScale>
  !!]
  type, extends(darkMatterHaloScaleClass) :: darkMatterHaloScaleVirialDensityContrastDefinition
     !!{
     A dark matter halo scale contrast class using virial density contrasts.
     !!}
     private
     logical                                                                       :: isRecursive                         , parentDeferred
     class           (darkMatterHaloScaleVirialDensityContrastDefinition), pointer :: recursiveSelf              => null()
     class           (cosmologyParametersClass                          ), pointer :: cosmologyParameters_       => null()
     class           (cosmologyFunctionsClass                           ), pointer :: cosmologyFunctions_        => null()
     class           (virialDensityContrastClass                        ), pointer :: virialDensityContrast_     => null()
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8                                    )          :: lastUniqueID
     ! Record of whether or not halo scales have already been computed for this node.
     logical                                                                       :: timescaleDynamicalComputed          , radiusVirialComputed            , &
          &                                                                           temperatureVirialComputed           , velocityVirialComputed
     ! Stored values of halo scales.
     double precision                                                              :: timescaleDynamicalStored            , radiusVirialStored              , &
          &                                                                           temperatureVirialStored             , velocityVirialStored            , &
          &                                                                           timePrevious                        , densityGrowthRatePrevious       , &
          &                                                                           massPrevious
     ! Table for fast lookup of the mean density of halos.
     double precision                                                              :: densityMeanTimeMaximum              , densityMeanTimeMinimum   =-1.0d0
     type            (table1DLogarithmicLinear                          )          :: densityMeanTable
   contains
     !![
     <methods>
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
     procedure :: deepCopy                            => virialDensityContrastDefinitionDeepCopy
     procedure :: deepCopyReset                       => virialDensityContrastDefinitionDeepCopyReset
     procedure :: deepCopyFinalize                    => virialDensityContrastDefinitionDeepCopyFinalize
  end type darkMatterHaloScaleVirialDensityContrastDefinition

  interface darkMatterHaloScaleVirialDensityContrastDefinition
     !!{
     Constructors for the {\normalfont \ttfamily virialDensityContrastDefinition} dark matter halo scales class.
     !!}
     module procedure virialDensityContrastDefinitionParameters
     module procedure virialDensityContrastDefinitionInternal
  end interface darkMatterHaloScaleVirialDensityContrastDefinition

  integer, parameter :: meanDensityTablePointsPerDecade=100

contains

  recursive function virialDensityContrastDefinitionParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily virialDensityContrastDefinition} dark matter halo scales class which takes a parameter set as input.
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
  end function virialDensityContrastDefinitionParameters

  recursive function virialDensityContrastDefinitionInternal(cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily virialDensityContrastDefinition} dark matter halo scales class.
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
    self%isRecursive               =.false.
    self%parentDeferred            =.false.
    return
  end function virialDensityContrastDefinitionInternal

  subroutine virialDensityContrastDefinitionAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self

    call calculationResetEvent%attach(self,virialDensityContrastDefinitionCalculationReset,openMPThreadBindingAllLevels,label='darkMatterHaloScaleVirialDensityContrastDefinition')
    return
  end subroutine virialDensityContrastDefinitionAutoHook

  subroutine virialDensityContrastDefinitionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily virialDensityContrastDefinition} dark matter halo scales class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self

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
    !!{
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
    !!{
    Returns the dynamical timescale for {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionDynamicalTimescale=self%recursiveSelf%timescaleDynamical(node)
       return
    end if
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
    !!{
    Returns the virial velocity scale for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic                                ), pointer       :: basic

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionVirialVelocity=self%recursiveSelf%velocityVirial(node)
       return
    end if
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
    !!{
    Returns the growth rate of the virial velocity scale for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic                                ), pointer       :: basic

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionVirialVelocityGrowthRate=self%recursiveSelf%velocityVirialGrowthRate(node)
       return
    end if
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
    !!{
    Returns the virial temperature (in Kelvin) for {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Astronomical, only : meanAtomicMassPrimordial
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionVirialTemperature=self%recursiveSelf%temperatureVirial(node)
       return
    end if
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
    !!{
    Returns the virial radius scale for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Math_Exponentiation     , only : cubeRoot
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic                                ), pointer       :: basic

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionVirialRadius=self%recursiveSelf%radiusVirial(node)
       return
    end if
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
    !!{
    Returns the logarithmic gradient of virial radius with halo mass at fixed epoch for {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionVirialRadiusGradientLogMass=self%recursiveSelf%radiusVirialGradientLogarithmicMass(node)
       return
    end if
    ! Halos at given epoch have fixed density, so radius always grows as the cube-root of mass.
    virialDensityContrastDefinitionVirialRadiusGradientLogMass=1.0d0/3.0d0
    return
  end function virialDensityContrastDefinitionVirialRadiusGradientLogMass

  double precision function virialDensityContrastDefinitionVirialRadiusGrowthRate(self,node)
    !!{
    Returns the growth rate of the virial radius scale for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type (treeNode                                          ), intent(inout) :: node
    class(nodeComponentBasic)                                , pointer       :: basic

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionVirialRadiusGrowthRate=self%recursiveSelf%radiusVirialGrowthRate(node)
       return
    end if
    ! Get the basic component.
    basic => node%basic()
    virialDensityContrastDefinitionVirialRadiusGrowthRate=(1.0d0/3.0d0)*self%radiusVirial(node)&
         &*(basic%accretionRate()/basic%mass()-self%densityMeanGrowthRate(node)&
         &/virialDensityContrastDefinitionMeanDensity(self,node))
    return
  end function virialDensityContrastDefinitionVirialRadiusGrowthRate

  double precision function virialDensityContrastDefinitionMeanDensity(self,node)
    !!{
    Returns the mean density for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type            (treeNode                                          ), intent(inout) :: node
    class           (nodeComponentBasic                                ), pointer       :: basic
    integer                                                                             :: i    , densityMeanTablePoints
    double precision                                                                    :: time

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionMeanDensity=self%recursiveSelf%densityMean(node)
       return
    end if
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
    !!{
    Returns the growth rate of the mean density for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    type            (treeNode                                          ), intent(inout) :: node
    class           (nodeComponentBasic                                ), pointer       :: basic
    double precision                                                                    :: expansionFactor, time

    ! Use recursive self if necessary.
    if (self%isRecursive) then
       virialDensityContrastDefinitionMeanDensityGrowthRate=self%recursiveSelf%densityMeanGrowthRate(node)
       return
    end if
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

  subroutine virialDensityContrastDefinitionDeepCopyReset(self)
    !!{
    Perform a deep copy reset of the object.
    !!}
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    
    self                           %   copiedSelf => null()
    if (.not.self%isRecursive) self%recursiveSelf => null()
    if (associated(self%cosmologyParameters_  )) call self%cosmologyParameters_  %deepCopyReset()
    if (associated(self%cosmologyFunctions_   )) call self%cosmologyFunctions_   %deepCopyReset()
    if (associated(self%virialDensityContrast_)) call self%virialDensityContrast_%deepCopyReset()
    return
  end subroutine virialDensityContrastDefinitionDeepCopyReset
  
  subroutine virialDensityContrastDefinitionDeepCopyFinalize(self)
    !!{
    Finalize a deep reset of the object.
    !!}
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self

    if (self%isRecursive) call virialDensityContrastFindParent(self)
    if (associated(self%cosmologyParameters_  )) call self%cosmologyParameters_  %deepCopyFinalize()
    if (associated(self%cosmologyFunctions_   )) call self%cosmologyFunctions_   %deepCopyFinalize()
    if (associated(self%virialDensityContrast_)) call self%virialDensityContrast_%deepCopyFinalize()
    return
  end subroutine virialDensityContrastDefinitionDeepCopyFinalize
  
  subroutine virialDensityContrastDefinitionDeepCopy(self,destination)
    !!{
    Perform a deep copy of the object.
    !!}
    use :: Error             , only : Error_Report
#ifdef OBJECTDEBUG
    use :: Display           , only : displayMessage            , verbosityLevelSilent
    use :: MPI_Utilities     , only : mpiSelf
    use :: Function_Classes  , only : debugReporting
    use :: ISO_Varying_String, only : operator(//)              , var_str
    use :: String_Handling   , only : operator(//)
#endif
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout), target :: self
    class(darkMatterHaloScaleClass                          ), intent(inout)         :: destination

    call self%darkMatterHaloScaleClass%deepCopy(destination)
    select type (destination)
    type is (darkMatterHaloScaleVirialDensityContrastDefinition)
       destination%isRecursive               =self%isRecursive
       destination%lastUniqueID              =self%lastUniqueID
       destination%timescaleDynamicalComputed=self%timescaleDynamicalComputed
       destination%radiusVirialComputed      =self%radiusVirialComputed
       destination%temperatureVirialComputed =self%temperatureVirialComputed
       destination%velocityVirialComputed    =self%velocityVirialComputed
       destination%timescaleDynamicalStored  =self%timescaleDynamicalStored
       destination%radiusVirialStored        =self%radiusVirialStored
       destination%temperatureVirialStored   =self%temperatureVirialStored
       destination%velocityVirialStored      =self%velocityVirialStored
       destination%timePrevious              =self%timePrevious
       destination%densityGrowthRatePrevious =self%densityGrowthRatePrevious
       destination%massPrevious              =self%massPrevious
       destination%densityMeanTimeMaximum    =self%densityMeanTimeMaximum
       destination%densityMeanTimeMinimum    =self%densityMeanTimeMinimum
       destination%densityMeanTable          =self%densityMeanTable
       destination%parentDeferred            =.false.
       if (self%isRecursive) then
          if (associated(self%recursiveSelf%recursiveSelf)) then
             ! If the parent self's recursiveSelf pointer is set, it indicates that it was deep-copied, and the pointer points to
             ! that copy. In that case we set the parent self of our destination to that copy.
             destination%recursiveSelf  => self%recursiveSelf%recursiveSelf
          else
             ! The parent self does not appear to have been deep-copied yet. Retain the same parent self pointer in our copy, but
             ! indicate that we need to look for the new parent later.
             destination%recursiveSelf  => self%recursiveSelf
             destination%parentDeferred =  .true.
          end if
       else
          ! This is a parent of a recursively-constructed object. Record the location of our copy so that it can be used to set
          ! the parent in deep copies of the child object.
          call virialDensityContrastDefinitionDeepCopyAssign(self,destination)
          destination%recursiveSelf     => null()
       end if       
       if (associated(self%cosmologyParameters_)) then
          if (associated(self%cosmologyParameters_%copiedSelf)) then
             select type(s => self%cosmologyParameters_%copiedSelf)
                class is (cosmologyParametersClass)
                destination%cosmologyParameters_ => s
                class default
                call Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%cosmologyParameters_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%cosmologyParameters_,mold=self%cosmologyParameters_)
             call self%cosmologyParameters_%deepCopy(destination%cosmologyParameters_)
             self%cosmologyParameters_%copiedSelf => destination%cosmologyParameters_
             call destination%cosmologyParameters_%autoHook()
          end if
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): cosmologyparameters : [destination] : ')//loc(destination)//' : '//loc(destination%cosmologyParameters_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
       if (associated(self%cosmologyFunctions_)) then
          if (associated(self%cosmologyFunctions_%copiedSelf)) then
             select type(s => self%cosmologyFunctions_%copiedSelf)
                class is (cosmologyFunctionsClass)
                destination%cosmologyFunctions_ => s
                class default
                call Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%cosmologyFunctions_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%cosmologyFunctions_,mold=self%cosmologyFunctions_)
             call self%cosmologyFunctions_%deepCopy(destination%cosmologyFunctions_)
             self%cosmologyFunctions_%copiedSelf => destination%cosmologyFunctions_
             call destination%cosmologyFunctions_%autoHook()
          end if
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): cosmologyfunctions : [destination] : ')//loc(destination)//' : '//loc(destination%cosmologyFunctions_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
       if (associated(self%virialDensityContrast_)) then
          if (associated(self%virialDensityContrast_%copiedSelf)) then
             select type(s => self%virialDensityContrast_%copiedSelf)
                class is (virialDensityContrastClass)
                destination%virialDensityContrast_ => s
                class default
                call Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%virialDensityContrast_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%virialDensityContrast_,mold=self%virialDensityContrast_)
             call self%virialDensityContrast_%deepCopy(destination%virialDensityContrast_)
             self%virialDensityContrast_%copiedSelf => destination%virialDensityContrast_
             call destination%virialDensityContrast_%autoHook()
          end if
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): virialdensitycontrast : [destination] : ')//loc(destination)//' : '//loc(destination%virialDensityContrast_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
       call destination%densityMeanTable%deepCopyActions()
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine virialDensityContrastDefinitionDeepCopy

  subroutine virialDensityContrastDefinitionDeepCopyAssign(self,destination)
    !!{
    Perform pointer assignment during a deep copy of the object.
    !!}
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)         :: self
    class(darkMatterHaloScaleClass                          ), intent(inout), target :: destination

    select type (destination)
    type is (darkMatterHaloScaleVirialDensityContrastDefinition)
       self%recursiveSelf => destination
    end select
    return
  end subroutine virialDensityContrastDefinitionDeepCopyAssign

  subroutine virialDensityContrastFindParent(self)
    !!{
    Find the deep-copied parent of a recursive child.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self

    if (self%parentDeferred) then
       if (associated(self%recursiveSelf%recursiveSelf)) then
          self%recursiveSelf => self%recursiveSelf%recursiveSelf
       else
         call Error_Report("recursive child's parent was not copied"//{introspection:location})
       end if
       self%parentDeferred=.false.
    end if
    return
  end subroutine virialDensityContrastFindParent

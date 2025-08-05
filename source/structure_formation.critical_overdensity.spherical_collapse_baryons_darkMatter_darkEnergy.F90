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
  An implementation of critical overdensity for collapse based on spherical collapse accounting for non-clustering of baryons.
  !!}

  use :: Cosmology_Parameters                 , only : cosmologyParameters                               , cosmologyParametersClass
  use :: Dark_Matter_Particles                , only : darkMatterParticle                                , darkMatterParticleClass
  use :: Intergalactic_Medium_Filtering_Masses, only : intergalacticMediumFilteringMass                  , intergalacticMediumFilteringMassClass
  use :: Spherical_Collapse_Solvers           , only : sphericalCollapseSolverBaryonsDarkMatterDarkEnergy, enumerationCllsnlssMttrDarkEnergyFixedAtType
  use :: Tables                               , only : table1D

  !![
  <criticalOverdensity name="criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy">
   <description>Critical overdensity for collapse based on the spherical collapse accounting for non-clustering of baryons.</description>
   <deepCopy>
    <functionClass variables="sphericalCollapseSolverClustered_, sphericalCollapseSolverUnclustered_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="sphericalCollapseSolverClustered_, sphericalCollapseSolverUnclustered_"/>
   </stateStorable>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy
     !!{
     A critical overdensity class based on spherical collapse accounting for non-clustering of baryons.
     !!}
     private
     logical                                                                           :: tableInitialized                  =  .false.
     double precision                                                                  :: tableClusteredTimeMinimum                   , tableClusteredTimeMaximum                    , &
          &                                                                               tableUnclusteredTimeMinimum                 , tableUnclusteredTimeMaximum
     double precision                                                                  :: normalization
     integer                                                                           :: tablePointsPerOctave
     logical                                                                           :: tableStore
     type            (enumerationCllsnlssMttrDarkEnergyFixedAtType      )              :: energyFixedAt
     class           (table1D                                           ), allocatable :: overdensityCriticalClustered                , overdensityCriticalUnclustered
     class           (darkMatterParticleClass                           ), pointer     :: darkMatterParticle_               => null()
     class           (cosmologyParametersClass                          ), pointer     :: cosmologyParameters_              => null()
     class           (intergalacticMediumFilteringMassClass             ), pointer     :: intergalacticMediumFilteringMass_ => null()
     type            (sphericalCollapseSolverBaryonsDarkMatterDarkEnergy), pointer     :: sphericalCollapseSolverClustered_ => null() , sphericalCollapseSolverUnclustered_ => null()
   contains
     !![
     <methods>
       <method description="Tabulate spherical collapse critical overdensity." method="retabulate" />
     </methods>
     !!]
     final     ::                    sphericalCollapseBrynsDrkMttrDrkEnrgyDestructor
     procedure :: value           => sphericalCollapseBrynsDrkMttrDrkEnrgyValue
     procedure :: gradientTime    => sphericalCollapseBrynsDrkMttrDrkEnrgyGradientTime
     procedure :: gradientMass    => sphericalCollapseBrynsDrkMttrDrkEnrgyGradientMass
     procedure :: retabulate      => sphericalCollapseBrynsDrkMttrDrkEnrgyRetabulate
     procedure :: isMassDependent => sphericalCollapseBrynsDrkMttrDrkEnrgyIsMassDependent
     procedure :: isNodeDependent => sphericalCollapseBrynsDrkMttrDrkEnrgyIsNodeDependent
     procedure :: isTreeDependent => sphericalCollapseBrynsDrkMttrDrkEnrgyIsTreeDependent
  end type criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy

  interface criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy
     !!{
     Constructors for the \refClass{criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy} critical overdensity for collapse class.
     !!}
     module procedure sphericalCollapseBrynsDrkMttrDrkEnrgyConstructorParameters
     module procedure sphericalCollapseBrynsDrkMttrDrkEnrgyConstructorInternal
  end interface criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy

contains

  function sphericalCollapseBrynsDrkMttrDrkEnrgyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy} critical overdensity class
    which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                                , inputParameters
    use :: Spherical_Collapse_Solvers, only : enumerationCllsnlssMttrDarkEnergyFixedAtEncode
    implicit none
    type            (criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy)                :: self
    type            (inputParameters                                         ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                                 ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass                                ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass                           ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                                 ), pointer       :: darkMatterParticle_
    class           (intergalacticMediumFilteringMassClass                   ), pointer       :: intergalacticMediumFilteringMass_
    double precision                                                                          :: normalization
    logical                                                                                   :: tableStore
    type            (varying_string                                          )                :: energyFixedAt
    integer                                                                                   :: tablePointsPerOctave

    !![
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>A normalizing factor to be applied to the critical overdensity.</description>
    </inputParameter>
    <inputParameter>
      <name>tableStore</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, store/restore the tabulated solution to/from file when possible.</description>
    </inputParameter>
    <inputParameter>
      <name>energyFixedAt</name>
      <defaultValue>var_str('turnaround')</defaultValue>
      <description>Selects the epoch at which the energy of a spherical top hat perturbation in a dark energy cosmology should be
        ``fixed'' for the purposes of computing virial density contrasts. (See the discussion in
        \citealt{percival_cosmological_2005}; \S8.)</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>tablePointsPerOctave</name>
      <source>parameters</source>
      <defaultValue>300</defaultValue>
      <description>The number of points per octave of time at which to tabulate solutions.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"               name="cosmologyFunctions_"               source="parameters"/>
    <objectBuilder class="cosmologyParameters"              name="cosmologyParameters_"              source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"         name="cosmologicalMassVariance_"         source="parameters"/>
    <objectBuilder class="darkMatterParticle"               name="darkMatterParticle_"               source="parameters"/>
    <objectBuilder class="intergalacticMediumFilteringMass" name="intergalacticMediumFilteringMass_" source="parameters"/>
    !!]
    self=criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy(cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,intergalacticMediumFilteringMass_,tableStore,tablePointsPerOctave,enumerationCllsnlssMttrDarkEnergyFixedAtEncode(char(energyFixedAt),includesPrefix=.false.),normalization)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"             />
    <objectDestructor name="cosmologyFunctions_"              />
    <objectDestructor name="cosmologicalMassVariance_"        />
    <objectDestructor name="darkMatterParticle_"              />
    <objectDestructor name="intergalacticMediumFilteringMass_"/>
    !!]
    return
  end function sphericalCollapseBrynsDrkMttrDrkEnrgyConstructorParameters

  function sphericalCollapseBrynsDrkMttrDrkEnrgyConstructorInternal(cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,intergalacticMediumFilteringMass_,tableStore,tablePointsPerOctave,energyFixedAt,normalization) result(self)
    !!{
    Internal constructor for the \refClass{criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy} critical overdensity class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleCDM
    use :: Error                , only : Error_Report
    implicit none
    type            (criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy)                          :: self
    class           (cosmologyFunctionsClass                                 ), target  , intent(in   ) :: cosmologyFunctions_
    class           (cosmologyParametersClass                                ), target  , intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass                           ), target  , intent(in   ) :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                                 ), target  , intent(in   ) :: darkMatterParticle_
    class           (intergalacticMediumFilteringMassClass                   ), target  , intent(in   ) :: intergalacticMediumFilteringMass_
    logical                                                                             , intent(in   ) :: tableStore
    integer                                                                             , intent(in   ) :: tablePointsPerOctave
    type            (enumerationCllsnlssMttrDarkEnergyFixedAtType            )          , intent(in   ) :: energyFixedAt
    double precision                                                          , optional, intent(in   ) :: normalization
    !![
    <optionalArgument name="normalization" defaultsTo="1.0d0" />
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_, *darkMatterParticle_, *intergalacticMediumFilteringMass_, tableStore, tablePointsPerOctave, energyFixedAt, normalization"/>
    !!]

    self%tableInitialized=.false.
    allocate(self%sphericalCollapseSolverClustered_  )
    allocate(self%sphericalCollapseSolverUnclustered_)
    !![
    <referenceConstruct isResult="yes" owner="self" object="sphericalCollapseSolverClustered_"   constructor="sphericalCollapseSolverBaryonsDarkMatterDarkEnergy(.true. ,self%tablePointsPerOctave,self%energyFixedAt,self%cosmologyParameters_,self%cosmologyFunctions_)"/>
    <referenceConstruct isResult="yes" owner="self" object="sphericalCollapseSolverUnclustered_" constructor="sphericalCollapseSolverBaryonsDarkMatterDarkEnergy(.false.,self%tablePointsPerOctave,self%energyFixedAt,self%cosmologyParameters_,self%cosmologyFunctions_)"/>
    !!]
    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Error_Report('critical overdensity expects a cold dark matter particle'//{introspection:location})
    end select
    return
  end function sphericalCollapseBrynsDrkMttrDrkEnrgyConstructorInternal

  subroutine sphericalCollapseBrynsDrkMttrDrkEnrgyDestructor(self)
    !!{
    Destructor for the \refClass{criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy} critical overdensity for collapse class.
    !!}
    implicit none
    type(criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"               />
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%darkMatterParticle_"                />
    <objectDestructor name="self%cosmologicalMassVariance_"          />
    <objectDestructor name="self%intergalacticMediumFilteringMass_"  />
    <objectDestructor name="self%sphericalCollapseSolverClustered_"  />
    <objectDestructor name="self%sphericalCollapseSolverUnclustered_"/>
    !!]
    if (self%tableInitialized) then
       call self%overdensityCriticalClustered  %destroy()
       call self%overdensityCriticalUnclustered%destroy()
       deallocate(self%overdensityCriticalClustered  )
       deallocate(self%overdensityCriticalUnclustered)
    end if
    return
  end subroutine sphericalCollapseBrynsDrkMttrDrkEnrgyDestructor

  subroutine sphericalCollapseBrynsDrkMttrDrkEnrgyRetabulate(self,time)
    !!{
    Recompute the look-up tables for critical overdensity for collapse.
    !!}
    implicit none
    class           (criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy), intent(inout) :: self
    double precision                                                          , intent(in   ) :: time
    logical                                                                                   :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableClusteredTimeMinimum .or. time > self%tableClusteredTimeMaximum .or. time < self%tableUnclusteredTimeMinimum .or. time > self%tableUnclusteredTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call self%sphericalCollapseSolverUnclustered_%criticalOverdensity(time,self%tableStore,self%overdensityCriticalUnclustered)
       call self%sphericalCollapseSolverClustered_  %criticalOverdensity(time,self%tableStore,self%overdensityCriticalClustered  )
       self%tableInitialized           =.true.
       self%tableClusteredTimeMinimum  =self%overdensityCriticalClustered  %x(+1)
       self%tableClusteredTimeMaximum  =self%overdensityCriticalClustered  %x(-1)
       self%tableUnclusteredTimeMinimum=self%overdensityCriticalUnclustered%x(+1)
       self%tableUnclusteredTimeMaximum=self%overdensityCriticalUnclustered%x(-1)
    end if
    return
  end subroutine sphericalCollapseBrynsDrkMttrDrkEnrgyRetabulate

  double precision function sphericalCollapseBrynsDrkMttrDrkEnrgyValue(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the critical overdensity at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    !!}
    implicit none
    class           (criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy), intent(inout)           :: self
    double precision                                                          , intent(in   ), optional :: time      , expansionFactor, &
         &                                                                                                 mass
    logical                                                                   , intent(in   ), optional :: collapsing
    type            (treeNode                                                ), intent(inout), optional :: node
    double precision                                                                                    :: time_     , interpolator
    !$GLC attributes unused :: node

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Construct an interpolation between the cases where baryons are clustered and unclustered. We use the same factor that
    ! appears in the suppression of baryonic accretion as a reasonable measure.
    interpolator=self%intergalacticMediumFilteringMass_%fractionBaryons(mass,time_)
    ! Interpolate the critical overdensity between the clustered and unclustered baryons case.
    sphericalCollapseBrynsDrkMttrDrkEnrgyValue=+(                                                                             &
         &                                       +self%overdensityCriticalClustered  %interpolate(time_)*       interpolator  &
         &                                       +self%overdensityCriticalUnclustered%interpolate(time_)*(1.0d0-interpolator) &
         &                                       )                                                                            &
         &                                     *  self%normalization

    return
  end function sphericalCollapseBrynsDrkMttrDrkEnrgyValue

  double precision function sphericalCollapseBrynsDrkMttrDrkEnrgyGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the time derivative of the critical overdensity at the given epoch, based spherical collapse in a matter plus
    cosmological constant universe.
    !!}
    implicit none
    class           (criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy), intent(inout)           :: self
    double precision                                                          , intent(in   ), optional :: time                    , expansionFactor, &
         &                                                                                                 mass
    logical                                                                   , intent(in   ), optional :: collapsing
    type            (treeNode                                                ), intent(inout), optional :: node
    double precision                                                                                    :: time_                   , interpolator   , &
         &                                                                                                 interpolatorRateOfChange
    !$GLC attributes unused :: node

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Construct an interpolation between the cases where baryons are clustered and unclustered. We use the same factor that
    ! appears in the suppression of baryonic accretion as a reasonable measure.
    interpolator            =self%intergalacticMediumFilteringMass_%fractionBaryons            (mass,time_)
    interpolatorRateOfChange=self%intergalacticMediumFilteringMass_%fractionBaryonsRateOfChange(mass,time_)
    ! Interpolate to get the expansion factor.
    sphericalCollapseBrynsDrkMttrDrkEnrgyGradientTime=+(                                                                                                   &
         &                                              +  self%overdensityCriticalClustered  %interpolateGradient(time_)*       interpolator              &
         &                                              +  self%overdensityCriticalUnclustered%interpolateGradient(time_)*(1.0d0-interpolator            ) &
         &                                              +(                                                                                                 &
         &                                                +self%overdensityCriticalClustered  %interpolate        (time_)                                  &
         &                                                -self%overdensityCriticalUnclustered%interpolate        (time_)                                  &
         &                                               )                                                                                                 &
         &                                              *                                                                        interpolatorRateOfChange  &
         &                                              )                                                                                                  &
         &                                              *  self%normalization
    return
  end function sphericalCollapseBrynsDrkMttrDrkEnrgyGradientTime

  double precision function sphericalCollapseBrynsDrkMttrDrkEnrgyGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to mass of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy), intent(inout)           :: self
    double precision                                                          , intent(in   ), optional :: time      , expansionFactor
    logical                                                                   , intent(in   ), optional :: collapsing
    double precision                                                          , intent(in   ), optional :: mass
    type            (treeNode                                                ), intent(inout), optional :: node
    double precision                                                                                    :: time_
    !$GLC attributes unused :: node

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    sphericalCollapseBrynsDrkMttrDrkEnrgyGradientMass=+(                                                                                &
         &                                              +self%overdensityCriticalClustered     %interpolate                (     time_) &
         &                                              +self%overdensityCriticalUnclustered   %interpolate                (     time_) &
         &                                             )                                                                                &
         &                                            *  self%intergalacticMediumFilteringMass_%fractionBaryonsGradientMass(mass,time_) &
         &                                            *  self%normalization
    return
  end function sphericalCollapseBrynsDrkMttrDrkEnrgyGradientMass

  logical function sphericalCollapseBrynsDrkMttrDrkEnrgyIsMassDependent(self)
    !!{
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy), intent(inout) :: self
    !$GLC attributes unused :: self

    sphericalCollapseBrynsDrkMttrDrkEnrgyIsMassDependent=.true.
    return
  end function sphericalCollapseBrynsDrkMttrDrkEnrgyIsMassDependent

  logical function sphericalCollapseBrynsDrkMttrDrkEnrgyIsNodeDependent(self)
    !!{
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy), intent(inout) :: self
    !$GLC attributes unused :: self

    sphericalCollapseBrynsDrkMttrDrkEnrgyIsNodeDependent=.false.
    return
  end function sphericalCollapseBrynsDrkMttrDrkEnrgyIsNodeDependent

  logical function sphericalCollapseBrynsDrkMttrDrkEnrgyIsTreeDependent(self)
    !!{
    Return whether the critical overdensity is tree dependent.
    !!}
    implicit none
    class(criticalOverdensitySphericalCollapseBrynsDrkMttrDrkEnrgy), intent(inout) :: self
    !$GLC attributes unused :: self

    sphericalCollapseBrynsDrkMttrDrkEnrgyIsTreeDependent=.false.
    return
  end function sphericalCollapseBrynsDrkMttrDrkEnrgyIsTreeDependent

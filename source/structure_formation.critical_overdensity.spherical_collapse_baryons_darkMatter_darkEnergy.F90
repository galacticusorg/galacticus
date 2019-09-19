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

  !% An implementation of critical overdensity for collapse based on spherical collapse accounting for non-clustering of baryons.

  use Tables                               , only : table1D
  use Cosmology_Parameters                 , only : cosmologyParameters                     , cosmologyParametersClass
  use Dark_Matter_Particles                , only : darkMatterParticle                      , darkMatterParticleClass
  use Intergalactic_Medium_Filtering_Masses, only : intergalacticMediumFilteringMass        , intergalacticMediumFilteringMassClass
  use Spherical_Collapse_Solvers           , only : sphericalCollapseSolverBaryonsDarkMatter

  !# <criticalOverdensity name="criticalOverdensitySphericalCollapseBaryonsDM">
  !#  <description>Critical overdensity for collapse based on the spherical collapse accounting for non-clustering of baryons.</description>
  !#  <deepCopy>
  !#   <functionClass variables="sphericalCollapseSolverClustered_, sphericalCollapseSolverUnclustered_"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="sphericalCollapseSolverClustered_, sphericalCollapseSolverUnclustered_"/>
  !#  </stateStorable>
  !# </criticalOverdensity>
  type, extends(criticalOverdensityClass) :: criticalOverdensitySphericalCollapseBaryonsDM
     !% A critical overdensity class based on spherical collapse accounting for non-clustering of baryons.
     private
     logical                                                                 :: tableInitialized
     double precision                                                        :: tableClusteredTimeMinimum                  , tableClusteredTimeMaximum                    , &
          &                                                                     tableUnclusteredTimeMinimum                , tableUnclusteredTimeMaximum
     double precision                                                        :: normalization
     logical                                                                 :: tableStore
     class           (table1D                                 ), allocatable :: overdensityCriticalClustered               , overdensityCriticalUnclustered
     class           (darkMatterParticleClass                 ), pointer     :: darkMatterParticle_               => null()
     class           (cosmologyParametersClass                ), pointer     :: cosmologyParameters_              => null()
     class           (intergalacticMediumFilteringMassClass   ), pointer     :: intergalacticMediumFilteringMass_ => null()
     type            (sphericalCollapseSolverBaryonsDarkMatter), pointer     :: sphericalCollapseSolverClustered_ => null(), sphericalCollapseSolverUnclustered_ => null()
   contains
     !@ <objectMethods>
     !@   <object>criticalOverdensitySphericalCollapseBaryonsDM</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate spherical collapse critical overdensity.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                    sphericalCollapseBaryonsDMDestructor
     procedure :: value           => sphericalCollapseBaryonsDMValue
     procedure :: gradientTime    => sphericalCollapseBaryonsDMGradientTime
     procedure :: gradientMass    => sphericalCollapseBaryonsDMGradientMass
     procedure :: retabulate      => sphericalCollapseBaryonsDMRetabulate
     procedure :: isMassDependent => sphericalCollapseBaryonsDMIsMassDependent
     procedure :: isNodeDependent => sphericalCollapseBaryonsDMIsNodeDependent
  end type criticalOverdensitySphericalCollapseBaryonsDM

  interface criticalOverdensitySphericalCollapseBaryonsDM
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseBaryonsDM} critical overdensity for collapse class.
     module procedure sphericalCollapseBaryonsDMConstructorParameters
     module procedure sphericalCollapseBaryonsDMConstructorInternal
  end interface criticalOverdensitySphericalCollapseBaryonsDM

contains

  function sphericalCollapseBaryonsDMConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseBaryonsDM} critical overdensity class
    !% which takes a parameter set as input.
    use Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (criticalOverdensitySphericalCollapseBaryonsDM)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                      ), pointer       :: cosmologyFunctions_    
    class           (cosmologyParametersClass                     ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass                ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                      ), pointer       :: darkMatterParticle_
    class           (intergalacticMediumFilteringMassClass        ), pointer       :: intergalacticMediumFilteringMass_
    double precision                                                               :: normalization
    logical                                                                        :: tableStore
    
    !# <inputParameter>
    !#   <name>normalization</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>A normalizing factor to be applied to the critical overdensity.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>tableStore</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If true, store/restore the tabulated solution to/from file when possible.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"               name="cosmologyFunctions_"               source="parameters"/>
    !# <objectBuilder class="cosmologyParameters"              name="cosmologyParameters_"              source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance"         name="cosmologicalMassVariance_"         source="parameters"/>
    !# <objectBuilder class="darkMatterParticle"               name="darkMatterParticle_"               source="parameters"/>
    !# <objectBuilder class="intergalacticMediumFilteringMass" name="intergalacticMediumFilteringMass_" source="parameters"/>
    self=criticalOverdensitySphericalCollapseBaryonsDM(cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,intergalacticMediumFilteringMass_,tableStore,normalization)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"             />
    !# <objectDestructor name="cosmologyFunctions_"              />
    !# <objectDestructor name="cosmologicalMassVariance_"        />
    !# <objectDestructor name="darkMatterParticle_"              />
    !# <objectDestructor name="intergalacticMediumFilteringMass_"/>
    return
  end function sphericalCollapseBaryonsDMConstructorParameters

  function sphericalCollapseBaryonsDMConstructorInternal(cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,intergalacticMediumFilteringMass_,tableStore,normalization) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseBaryonsDM} critical overdensity class.
    use Dark_Matter_Particles     , only : darkMatterParticleCDM
    use Galacticus_Error          , only : Galacticus_Error_Report
    use Spherical_Collapse_Solvers, only : matterDarkEnergyFixedAtUndefined
    implicit none
    type            (criticalOverdensitySphericalCollapseBaryonsDM)                          :: self
    class           (cosmologyFunctionsClass                      ), target  , intent(in   ) :: cosmologyFunctions_    
    class           (cosmologyParametersClass                     ), target  , intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass                ), target  , intent(in   ) :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                      ), target  , intent(in   ) :: darkMatterParticle_
    class           (intergalacticMediumFilteringMassClass        ), target  , intent(in   ) :: intergalacticMediumFilteringMass_
    logical                                                                  , intent(in   ) :: tableStore
    double precision                                               , optional, intent(in   ) :: normalization
    !# <optionalArgument name="normalization" defaultsTo="1.0d0" />
    !# <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_, *darkMatterParticle_, *intergalacticMediumFilteringMass_, tableStore, normalization"/>

    self%tableInitialized=.false.
    allocate(self%sphericalCollapseSolverClustered_  )
    allocate(self%sphericalCollapseSolverUnclustered_)
    !# <referenceConstruct isResult="yes" owner="self" object="sphericalCollapseSolverClustered_"   constructor="sphericalCollapseSolverBaryonsDarkMatter(.true. ,matterDarkEnergyFixedAtUndefined,self%cosmologyParameters_,self%cosmologyFunctions_)"/>
    !# <referenceConstruct isResult="yes" owner="self" object="sphericalCollapseSolverUnclustered_" constructor="sphericalCollapseSolverBaryonsDarkMatter(.false.,matterDarkEnergyFixedAtUndefined,self%cosmologyParameters_,self%cosmologyFunctions_)"/>
    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Galacticus_Error_Report('critical overdensity expects a cold dark matter particle'//{introspection:location})
    end select
    return
  end function sphericalCollapseBaryonsDMConstructorInternal
  
  subroutine sphericalCollapseBaryonsDMDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseBaryonsDM} critical overdensity for collapse class.
    implicit none
    type(criticalOverdensitySphericalCollapseBaryonsDM), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"               />
    !# <objectDestructor name="self%cosmologyFunctions_"                />
    !# <objectDestructor name="self%darkMatterParticle_"                />
    !# <objectDestructor name="self%cosmologicalMassVariance_"          />
    !# <objectDestructor name="self%intergalacticMediumFilteringMass_"  />
    !# <objectDestructor name="self%sphericalCollapseSolverClustered_"  />
    !# <objectDestructor name="self%sphericalCollapseSolverUnclustered_"/>
    if (self%tableInitialized) then
       call self%overdensityCriticalClustered  %destroy()
       call self%overdensityCriticalUnclustered%destroy()
       deallocate(self%overdensityCriticalClustered  )
       deallocate(self%overdensityCriticalUnclustered)
    end if
    return
  end subroutine sphericalCollapseBaryonsDMDestructor

  subroutine sphericalCollapseBaryonsDMRetabulate(self,time)
    !% Recompute the look-up tables for critical overdensity for collapse.
    implicit none
    class           (criticalOverdensitySphericalCollapseBaryonsDM), intent(inout) :: self
    double precision                                               , intent(in   ) :: time
    logical                                                                        :: remakeTable

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
  end subroutine sphericalCollapseBaryonsDMRetabulate

  double precision function sphericalCollapseBaryonsDMValue(self,time,expansionFactor,collapsing,mass,node)
    !% Return the critical overdensity at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    implicit none
    class           (criticalOverdensitySphericalCollapseBaryonsDM), intent(inout)           :: self
    double precision                                               , intent(in   ), optional :: time      , expansionFactor, &
         &                                                                                      mass
    logical                                                        , intent(in   ), optional :: collapsing
    type            (treeNode                                     ), intent(inout), optional :: node
    double precision                                                                         :: time_     , interpolator
    !GCC$ attributes unused :: node
    
    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Construct an interpolation between the cases where baryons are clustered and unclustered. We use the same factor that
    ! appears in the suppression of baryonic accretion as a reasonable measure.
    interpolator=self%intergalacticMediumFilteringMass_%fractionBaryons(mass,time_)
    ! Interpolate the critical overdensity between the clustered and unclustered baryons case.
    sphericalCollapseBaryonsDMValue=+(                                                                             &
         &                            +self%overdensityCriticalClustered  %interpolate(time_)*       interpolator  &
         &                            +self%overdensityCriticalUnclustered%interpolate(time_)*(1.0d0-interpolator) &
         &                            )                                                                            &
         &                          *  self%normalization

    return
  end function sphericalCollapseBaryonsDMValue

  double precision function sphericalCollapseBaryonsDMGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !% Return the time derivative of the critical overdensity at the given epoch, based spherical collapse in a matter plus
    !% cosmological constant universe.
    implicit none
    class           (criticalOverdensitySphericalCollapseBaryonsDM), intent(inout)           :: self
    double precision                                               , intent(in   ), optional :: time                    , expansionFactor, &
         &                                                                                      mass
    logical                                                        , intent(in   ), optional :: collapsing
    type            (treeNode                                     ), intent(inout), optional :: node
    double precision                                                                         :: time_                   , interpolator   , &
         &                                                                                      interpolatorRateOfChange
    !GCC$ attributes unused :: node

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Construct an interpolation between the cases where baryons are clustered and unclustered. We use the same factor that
    ! appears in the suppression of baryonic accretion as a reasonable measure.
    interpolator            =self%intergalacticMediumFilteringMass_%fractionBaryons            (mass,time_)
    interpolatorRateOfChange=self%intergalacticMediumFilteringMass_%fractionBaryonsRateOfChange(mass,time_)
    ! Interpolate to get the expansion factor.
    sphericalCollapseBaryonsDMGradientTime=+(                                                                                                   &
         &                                   +  self%overdensityCriticalClustered  %interpolateGradient(time_)*       interpolator              &
         &                                   +  self%overdensityCriticalUnclustered%interpolateGradient(time_)*(1.0d0-interpolator            ) &
         &                                   +(                                                                                                 &
         &                                     +self%overdensityCriticalClustered  %interpolate        (time_)                                  &
         &                                     -self%overdensityCriticalUnclustered%interpolate        (time_)                                  &
         &                                    )                                                                                                 &
         &                                   *                                                                        interpolatorRateOfChange  &
         &                                   )                                                                                                  &
         &                                   *  self%normalization
    return
  end function sphericalCollapseBaryonsDMGradientTime

  double precision function sphericalCollapseBaryonsDMGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !% Return the gradient with respect to mass of critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensitySphericalCollapseBaryonsDM), intent(inout)           :: self
    double precision                                               , intent(in   ), optional :: time      , expansionFactor
    logical                                                        , intent(in   ), optional :: collapsing
    double precision                                               , intent(in   ), optional :: mass
    type            (treeNode                                     ), intent(inout), optional :: node
    double precision                                                                         :: time_
    !GCC$ attributes unused :: node
    
    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    sphericalCollapseBaryonsDMGradientMass=+(                                                                                &
         &                                   +self%overdensityCriticalClustered     %interpolate                (     time_) &
         &                                   +self%overdensityCriticalUnclustered   %interpolate                (     time_) &
         &                                  )                                                                                &
         &                                 *  self%intergalacticMediumFilteringMass_%fractionBaryonsGradientMass(mass,time_) &
         &                                 *  self%normalization
    return
  end function sphericalCollapseBaryonsDMGradientMass

  logical function sphericalCollapseBaryonsDMIsMassDependent(self)
    !% Return whether the critical overdensity is mass dependent.
    implicit none
    class(criticalOverdensitySphericalCollapseBaryonsDM), intent(inout) :: self
    !GCC$ attributes unused :: self

    sphericalCollapseBaryonsDMIsMassDependent=.true.
    return
  end function sphericalCollapseBaryonsDMIsMassDependent

  logical function sphericalCollapseBaryonsDMIsNodeDependent(self)
    !% Return whether the critical overdensity is node dependent.
    implicit none
    class(criticalOverdensitySphericalCollapseBaryonsDM), intent(inout) :: self
    !GCC$ attributes unused :: self

    sphericalCollapseBaryonsDMIsNodeDependent=.false.
    return
  end function sphericalCollapseBaryonsDMIsNodeDependent

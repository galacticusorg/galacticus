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

  !% An implementation of critical overdensity for collapse based on spherical collapse in a
  !% matter plus cosmological constant universe.

  use Tables
  use Linear_Growth
  use Cosmology_Functions
  use Dark_Matter_Particles

  !# <criticalOverdensity name="criticalOverdensitySphericalCollapseMatterLambda">
  !#  <description>Critical overdensity for collapse based on the spherical collapse in a matter plus cosmological constant universe (see, for example, \citealt{percival_cosmological_2005}).</description>
  !# </criticalOverdensity>
  type, extends(criticalOverdensityClass) :: criticalOverdensitySphericalCollapseMatterLambda
     !% A critical overdensity class based on spherical collapse in a matter plus cosmological constant universe.
     private
     logical                                                :: tableInitialized    = .false.
     double precision                                       :: tableTimeMinimum             , tableTimeMaximum
     double precision                                       :: normalization
     logical                                                :: tableStore
     class           (table1D                ), allocatable :: overdensityCritical
     class           (linearGrowthClass      ), pointer     :: linearGrowth_       => null()
     class           (darkMatterParticleClass), pointer     :: darkMatterParticle_ => null()
   contains
     !@ <objectMethods>
     !@   <object>criticalOverdensitySphericalCollapseMatterLambda</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate spherical collapse critical overdensity.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                    sphericalCollapseMatterLambdaDestructor
     procedure :: value           => sphericalCollapseMatterLambdaValue
     procedure :: gradientTime    => sphericalCollapseMatterLambdaGradientTime
     procedure :: gradientMass    => sphericalCollapseMatterLambdaGradientMass
     procedure :: retabulate      => sphericalCollapseMatterLambdaRetabulate
     procedure :: isMassDependent => sphericalCollapseMatterLambdaIsMassDependent
  end type criticalOverdensitySphericalCollapseMatterLambda

  interface criticalOverdensitySphericalCollapseMatterLambda
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseMatterLambda} critical overdensity for collapse class.
     module procedure sphericalCollapseMatterLambdaConstructorParameters
     module procedure sphericalCollapseMatterLambdaConstructorInternal
  end interface criticalOverdensitySphericalCollapseMatterLambda

contains

  function sphericalCollapseMatterLambdaConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} critical overdensity class
    !% which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (criticalOverdensitySphericalCollapseMatterLambda)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                         ), pointer       :: cosmologyFunctions_    
    class           (linearGrowthClass                               ), pointer       :: linearGrowth_    
    class           (cosmologicalMassVarianceClass                   ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                         ), pointer       :: darkMatterParticle_
    double precision                                                                  :: normalization
    logical                                                                           :: tableStore
    
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
    !# <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !# <objectBuilder class="darkMatterParticle"       name="darkMatterParticle_"       source="parameters"/>
    self=criticalOverdensitySphericalCollapseMatterLambda(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,tableStore,normalization)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="linearGrowth_"            />
    !# <objectDestructor name="cosmologyFunctions_"      />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    !# <objectDestructor name="darkMatterParticle_"      />
    return
  end function sphericalCollapseMatterLambdaConstructorParameters

  function sphericalCollapseMatterLambdaConstructorInternal(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,tableStore,normalization) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} critical overdensity class.
    use Dark_Matter_Particles
    use Galacticus_Error
    implicit none
    type            (criticalOverdensitySphericalCollapseMatterLambda)                          :: self
    class           (cosmologyFunctionsClass                         ), target  , intent(in   ) :: cosmologyFunctions_    
    class           (linearGrowthClass                               ), target  , intent(in   ) :: linearGrowth_    
    class           (cosmologicalMassVarianceClass                   ), target  , intent(in   ) :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                         ), target  , intent(in   ) :: darkMatterParticle_
    logical                                                                     , intent(in   ) :: tableStore
    double precision                                                  , optional, intent(in   ) :: normalization
    !# <optionalArgument name="normalization" defaultsTo="1.0d0" />
    !# <constructorAssign variables="*linearGrowth_, *cosmologyFunctions_, *cosmologicalMassVariance_, *darkMatterParticle_, tableStore, normalization"/>
    
    self%tableInitialized=.false.
    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Galacticus_Error_Report('critical overdensity expects a cold dark matter particle'//{introspection:location})
    end select
    return
  end function sphericalCollapseMatterLambdaConstructorInternal
  
  subroutine sphericalCollapseMatterLambdaDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseMatterLambda} critical overdensity for collapse class.
    implicit none
    type(criticalOverdensitySphericalCollapseMatterLambda), intent(inout) :: self

    !# <objectDestructor name="self%linearGrowth_"            />
    !# <objectDestructor name="self%cosmologyFunctions_"      />
    !# <objectDestructor name="self%darkMatterParticle_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    if (self%tableInitialized) then
       call self%overdensityCritical%destroy()
       deallocate(self%overdensityCritical)
    end if
    return
  end subroutine sphericalCollapseMatterLambdaDestructor

  subroutine sphericalCollapseMatterLambdaRetabulate(self,time)
    !% Recompute the look-up tables for critical overdensity for collapse.
    use Spherical_Collapse_Matter_Lambda
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterLambda), intent(inout) :: self
    double precision                                                  , intent(in   ) :: time
    logical                                                                           :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableTimeMinimum .or. time > self%tableTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_Matter_Lambda_Critical_Overdensity_Tabulate(time,self%tableStore,self%overdensityCritical,self%cosmologyFunctions_,self%linearGrowth_)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%overdensityCritical%x(+1)
       self%tableTimeMaximum=self%overdensityCritical%x(-1)
    end if
    return
  end subroutine sphericalCollapseMatterLambdaRetabulate

  double precision function sphericalCollapseMatterLambdaValue(self,time,expansionFactor,collapsing,mass,node)
    !% Return the critical overdensity at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use Galacticus_Error
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                  , intent(in   ), optional :: time               , expansionFactor, &
         &                                                                                         mass
    logical                                                           , intent(in   ), optional :: collapsing
    type            (treeNode                                        ), intent(inout), optional :: node
    double precision                                                                            :: time_
    !GCC$ attributes unused :: mass, node
    
    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    sphericalCollapseMatterLambdaValue=+self%overdensityCritical%interpolate(time_) &
         &                             *self%normalization
    return
  end function sphericalCollapseMatterLambdaValue

  double precision function sphericalCollapseMatterLambdaGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !% Return the time derivative of the critical overdensity at the given epoch, based spherical collapse in a matter plus
    !% cosmological constant universe.
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                  , intent(in   ), optional :: time               , expansionFactor, &
         &                                                                                         mass
    logical                                                           , intent(in   ), optional :: collapsing
    type            (treeNode                                        ), intent(inout), optional :: node
    double precision                                                                            :: time_
    !GCC$ attributes unused :: mass, node

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    sphericalCollapseMatterLambdaGradientTime=+self%overdensityCritical%interpolateGradient(time_) &
         &                                    *self%normalization
    return
  end function sphericalCollapseMatterLambdaGradientTime

  double precision function sphericalCollapseMatterLambdaGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !% Return the gradient with respect to mass of critical overdensity at the given time and mass.
    use Linear_Growth
    use Cosmology_Functions
    implicit none
    class           (criticalOverdensitySphericalCollapseMatterLambda), intent(inout)           :: self
    double precision                                                  , intent(in   ), optional :: time      , expansionFactor
    logical                                                           , intent(in   ), optional :: collapsing
    double precision                                                  , intent(in   ), optional :: mass
    type            (treeNode                                        ), intent(inout), optional :: node
    !GCC$ attributes unused :: self, time, expansionFactor, collapsing, mass, node
    
    sphericalCollapseMatterLambdaGradientMass=0.0d0
    return
  end function sphericalCollapseMatterLambdaGradientMass

  logical function sphericalCollapseMatterLambdaIsMassDependent(self)
    !% Return whether the critical overdensity is mass dependent.
    implicit none
    class(criticalOverdensitySphericalCollapseMatterLambda), intent(inout) :: self
    !GCC$ attributes unused :: self

    sphericalCollapseMatterLambdaIsMassDependent=.false.
    return
  end function sphericalCollapseMatterLambdaIsMassDependent

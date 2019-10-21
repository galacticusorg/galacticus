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

  use :: Dark_Matter_Particles     , only : darkMatterParticleClass
  use :: Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt
  use :: Tables                    , only : table1D

  !# <criticalOverdensity name="criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt">
  !#  <description>Critical overdensity for collapse based on the spherical collapse in a matter plus cosmological constant universe (see, for example, \citealt{percival_cosmological_2005}).</description>
  !#  <deepCopy>
  !#   <functionClass variables="sphericalCollapseSolver_"/>
  !#  </deepCopy>
  !#  <stateStorable>
  !#   <functionClass variables="sphericalCollapseSolver_"/>
  !#  </stateStorable>
  !# </criticalOverdensity>
  type, extends(criticalOverdensityClass) :: criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt
     !% A critical overdensity class based on spherical collapse in a matter plus cosmological constant universe.
     private
     logical                                                                         :: tableInitialized         = .false.
     double precision                                                                :: tableTimeMinimum                  , tableTimeMaximum
     double precision                                                                :: normalization
     logical                                                                         :: tableStore
     class           (table1D                                         ), allocatable :: overdensityCritical
     class           (darkMatterParticleClass                         ), pointer     :: darkMatterParticle_      => null()
     class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt), pointer     :: sphericalCollapseSolver_ => null()
   contains
     !@ <objectMethods>
     !@   <object>criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate spherical collapse critical overdensity.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                    sphericalCollapseCllsnlssMttrCsmlgclCnstntDestructor
     procedure :: value           => sphericalCollapseCllsnlssMttrCsmlgclCnstntValue
     procedure :: gradientTime    => sphericalCollapseCllsnlssMttrCsmlgclCnstntGradientTime
     procedure :: gradientMass    => sphericalCollapseCllsnlssMttrCsmlgclCnstntGradientMass
     procedure :: retabulate      => sphericalCollapseCllsnlssMttrCsmlgclCnstntRetabulate
     procedure :: isMassDependent => sphericalCollapseCllsnlssMttrCsmlgclCnstntIsMassDependent
     procedure :: isNodeDependent => sphericalCollapseCllsnlssMttrCsmlgclCnstntIsNodeDependent
  end type criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt

  interface criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrCsmlgclCnstnt} critical overdensity for collapse class.
     module procedure sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorParameters
     module procedure sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorInternal
  end interface criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt

contains

  function sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrCsmlgclCnstnt} critical overdensity class
    !% which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt)                :: self
    type            (inputParameters                                              ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                                      ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass                                            ), pointer       :: linearGrowth_
    class           (cosmologicalMassVarianceClass                                ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                                      ), pointer       :: darkMatterParticle_
    double precision                                                                               :: normalization
    logical                                                                                        :: tableStore

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
    self=criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,tableStore,normalization)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="linearGrowth_"            />
    !# <objectDestructor name="cosmologyFunctions_"      />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    !# <objectDestructor name="darkMatterParticle_"      />
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorParameters

  function sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorInternal(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,tableStore,normalization) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrCsmlgclCnstnt} critical overdensity class.
    use :: Dark_Matter_Particles, only : darkMatterParticleCDM  , darkMatterParticleClass
    use :: Galacticus_Error     , only : Galacticus_Error_Report
    implicit none
    type            (criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt)                          :: self
    class           (cosmologyFunctionsClass                                      ), target  , intent(in   ) :: cosmologyFunctions_
    class           (linearGrowthClass                                            ), target  , intent(in   ) :: linearGrowth_
    class           (cosmologicalMassVarianceClass                                ), target  , intent(in   ) :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                                      ), target  , intent(in   ) :: darkMatterParticle_
    logical                                                                                  , intent(in   ) :: tableStore
    double precision                                                               , optional, intent(in   ) :: normalization
    !# <optionalArgument name="normalization" defaultsTo="1.0d0" />
    !# <constructorAssign variables="*linearGrowth_, *cosmologyFunctions_, *cosmologicalMassVariance_, *darkMatterParticle_, tableStore"/>

    self%normalization   =normalization_
    self%tableInitialized=.false.
    allocate(sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt :: self%sphericalCollapseSolver_)
    select type (sphericalCollapseSolver_ => self%sphericalCollapseSolver_)
    type is (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)
       !# <referenceConstruct isResult="yes" object="sphericalCollapseSolver_" constructor="sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt(self%cosmologyFunctions_,self%linearGrowth_)"/>
    end select
    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Galacticus_Error_Report('critical overdensity expects a cold dark matter particle'//{introspection:location})
    end select
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntConstructorInternal

  subroutine sphericalCollapseCllsnlssMttrCsmlgclCnstntDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseCllsnlssMttrCsmlgclCnstnt} critical overdensity for collapse class.
    implicit none
    type(criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout) :: self

    !# <objectDestructor name="self%linearGrowth_"            />
    !# <objectDestructor name="self%cosmologyFunctions_"      />
    !# <objectDestructor name="self%darkMatterParticle_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    !# <objectDestructor name="self%sphericalCollapseSolver_" />
    if (self%tableInitialized) then
       call self%overdensityCritical%destroy()
       deallocate(self%overdensityCritical)
    end if
    return
  end subroutine sphericalCollapseCllsnlssMttrCsmlgclCnstntDestructor

  subroutine sphericalCollapseCllsnlssMttrCsmlgclCnstntRetabulate(self,time)
    !% Recompute the look-up tables for critical overdensity for collapse.
    implicit none
    class           (criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout) :: self
    double precision                                                               , intent(in   ) :: time
    logical                                                                                        :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableTimeMinimum .or. time > self%tableTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call self%sphericalCollapseSolver_%criticalOverdensity(time,self%tableStore,self%overdensityCritical)
       self%tableInitialized=.true.
       self%tableTimeMinimum=self%overdensityCritical%x(+1)
       self%tableTimeMaximum=self%overdensityCritical%x(-1)
    end if
    return
  end subroutine sphericalCollapseCllsnlssMttrCsmlgclCnstntRetabulate

  double precision function sphericalCollapseCllsnlssMttrCsmlgclCnstntValue(self,time,expansionFactor,collapsing,mass,node)
    !% Return the critical overdensity at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                               , intent(in   ), optional :: time      , expansionFactor, &
         &                                                                                                      mass
    logical                                                                        , intent(in   ), optional :: collapsing
    type            (treeNode                                                     ), intent(inout), optional :: node
    double precision                                                                                         :: time_
    !GCC$ attributes unused :: mass, node

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the critical overdensity.
    sphericalCollapseCllsnlssMttrCsmlgclCnstntValue=+self%overdensityCritical%interpolate(time_) &
         &                                          *self%normalization
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntValue

  double precision function sphericalCollapseCllsnlssMttrCsmlgclCnstntGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !% Return the time derivative of the critical overdensity at the given epoch, based spherical collapse in a matter plus
    !% cosmological constant universe.
    implicit none
    class           (criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                               , intent(in   ), optional :: time      , expansionFactor, &
         &                                                                                                      mass
    logical                                                                        , intent(in   ), optional :: collapsing
    type            (treeNode                                                     ), intent(inout), optional :: node
    double precision                                                                                         :: time_
    !GCC$ attributes unused :: mass, node

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    sphericalCollapseCllsnlssMttrCsmlgclCnstntGradientTime=+self%overdensityCritical%interpolateGradient(time_) &
         &                                                 *self%normalization
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntGradientTime

  double precision function sphericalCollapseCllsnlssMttrCsmlgclCnstntGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !% Return the gradient with respect to mass of critical overdensity at the given time and mass.
    implicit none
    class           (criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout)           :: self
    double precision                                                               , intent(in   ), optional :: time      , expansionFactor
    logical                                                                        , intent(in   ), optional :: collapsing
    double precision                                                               , intent(in   ), optional :: mass
    type            (treeNode                                                     ), intent(inout), optional :: node
    !GCC$ attributes unused :: self, time, expansionFactor, collapsing, mass, node

    sphericalCollapseCllsnlssMttrCsmlgclCnstntGradientMass=0.0d0
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntGradientMass

  logical function sphericalCollapseCllsnlssMttrCsmlgclCnstntIsMassDependent(self)
    !% Return whether the critical overdensity is mass dependent.
    implicit none
    class(criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout) :: self
    !GCC$ attributes unused :: self

    sphericalCollapseCllsnlssMttrCsmlgclCnstntIsMassDependent=.false.
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntIsMassDependent

  logical function sphericalCollapseCllsnlssMttrCsmlgclCnstntIsNodeDependent(self)
    !% Return whether the critical overdensity is node dependent.
    implicit none
    class(criticalOverdensitySphericalCollapseCllsnlssMttrCsmlgclCnstnt), intent(inout) :: self
    !GCC$ attributes unused :: self

    sphericalCollapseCllsnlssMttrCsmlgclCnstntIsNodeDependent=.false.
    return
  end function sphericalCollapseCllsnlssMttrCsmlgclCnstntIsNodeDependent

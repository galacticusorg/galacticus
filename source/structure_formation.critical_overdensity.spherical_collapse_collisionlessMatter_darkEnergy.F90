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

  !!{
  An implementation of critical overdensity for collapse based on spherical collapse in a
  matter plus dark energy universe.
  !!}

  !![
  <criticalOverdensity name="criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy">
   <description>Critical overdensity for collapse based on the spherical collapse in a matter plus dark energy universe.</description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensitySphericalCollapseClsnlssMttrCsmlgclCnstnt) :: criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy
     !!{
     A dark matter halo virial density contrast class based on spherical collapse in a matter plus dark energy universe.
     !!}
     private
  end type criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy

  interface criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy
     !!{
     Constructors for the \refClass{criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy} critical overdensity for collapse class.
     !!}
     module procedure sphericalCollapseClsnlssMttrDrkEnrgyConstructorParameters
     module procedure sphericalCollapseClsnlssMttrDrkEnrgyConstructorInternal
  end interface criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy

contains

  function sphericalCollapseClsnlssMttrDrkEnrgyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy} critical overdensity class
    which takes a parameter set as input.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticle, darkMatterParticleClass
    use :: Error                , only : Error_Report
    use :: Input_Parameters     , only : inputParameter    , inputParameters
    implicit none
    type            (criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy)                :: self
    type            (inputParameters                                        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                                ), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass                                      ), pointer       :: linearGrowth_
    class           (cosmologicalMassVarianceClass                          ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                                ), pointer       :: darkMatterParticle_
    double precision                                                                         :: normalization            , countTimeCollapsePerUnit
    logical                                                                                  :: tableStore

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
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"       name="darkMatterParticle_"       source="parameters"/>
    !!]
    if (parameters%isPresent('countTimeCollapsePerUnit')) then
       !![
       <inputParameter>
	 <name>countTimeCollapsePerUnit</name>
	 <source>parameters</source>
	 <description>The number of points per unit $w(t)=\delta_\mathrm{c}(t)/D(t)$ at which to tabulate the time of collapse.</description>
       </inputParameter>
       !!]
    end if
    !![
    <conditionalCall>
      <call>self=criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,tableStore,normalization{conditions})</call>
      <argument name="countTimeCollapsePerUnit" value="countTimeCollapsePerUnit" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="darkMatterParticle_"      />
    !!]
    return
  end function sphericalCollapseClsnlssMttrDrkEnrgyConstructorParameters

  function sphericalCollapseClsnlssMttrDrkEnrgyConstructorInternal(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,tableStore,normalization,countTimeCollapsePerUnit) result(self)
    !!{
    Internal constructor for the \refClass{criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy} critical overdensity class.
    !!}
    use :: Dark_Matter_Particles     , only : darkMatterParticleCDM                 , darkMatterParticleClass
    use :: Error                     , only : Error_Report
    use :: Spherical_Collapse_Solvers, only : cllsnlssMttrDarkEnergyFixedAtUndefined, sphericalCollapseSolverCllsnlssMttrDarkEnergy
    implicit none
    type            (criticalOverdensitySphericalCollapseClsnlssMttrDrkEnrgy)                          :: self
    class           (cosmologyFunctionsClass                                ), target  , intent(in   ) :: cosmologyFunctions_
    class           (linearGrowthClass                                      ), target  , intent(in   ) :: linearGrowth_
    class           (cosmologicalMassVarianceClass                          ), target  , intent(in   ) :: cosmologicalMassVariance_
    class           (darkMatterParticleClass                                ), target  , intent(in   ) :: darkMatterParticle_
    logical                                                                            , intent(in   ) :: tableStore
    double precision                                                         , optional, intent(in   ) :: normalization            , countTimeCollapsePerUnit
    !![
    <optionalArgument name="normalization" defaultsTo="1.0d0" />
    <constructorAssign variables="*linearGrowth_, *cosmologyFunctions_, *cosmologicalMassVariance_, *darkMatterParticle_, normalization, tableStore, countTimeCollapsePerUnit"/>
    !!]

    self%tableInitialized=.false.
    allocate(sphericalCollapseSolverCllsnlssMttrDarkEnergy :: self%sphericalCollapseSolver_)
    select type (sphericalCollapseSolver_ => self%sphericalCollapseSolver_)
    type is (sphericalCollapseSolverCllsnlssMttrDarkEnergy)
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="sphericalCollapseSolver_" object="sphericalCollapseSolver_" constructor="sphericalCollapseSolverCllsnlssMttrDarkEnergy(cllsnlssMttrDarkEnergyFixedAtUndefined,self%cosmologyFunctions_,self%linearGrowth_)"/>
       !!]
    end select
    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Error_Report('critical overdensity expects a cold dark matter particle'//{introspection:location})
    end select
    return
  end function sphericalCollapseClsnlssMttrDrkEnrgyConstructorInternal

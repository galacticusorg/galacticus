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
Implements a critical overdensity for collapse using the revised spherical collapse model for decaying
dark matter (DDM) of :cite:t:`montandon_decaying_2026`.
!!}
  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <criticalOverdensity name="criticalOverdensityDecayingDarkMatter" docformat="rst">
   <description>
   A critical overdensity for collapse class implementing the revised spherical collapse model for
   decaying dark matter (:term:`DDM`) of :cite:t:`montandon_decaying_2026`. The decay of dark matter
   into a massive daughter (with velocity kick :math:`v_k`) plus dark radiation induces a continuous
   mass loss during collapse, which softens the potential and delays collapse, raising the linear
   critical overdensity for collapse in a mass-dependent way. This class wraps another (:term:`CDM`)
   critical overdensity class and multiplies it by the ratio
   :math:`\delta_\mathrm{c}^\mathrm{fit}(M_0)/\delta_\mathrm{c}^\mathrm{EdS}`, where
   :math:`\delta_\mathrm{c}^\mathrm{fit}(M_0)` is the DDM critical overdensity fitting function of
   :cite:t:`montandon_decaying_2026` (their eqs. 40, 42, 43) and
   :math:`\delta_\mathrm{c}^\mathrm{EdS}\approx1.686` (their eq. 32). Because the correction is a ratio
   normalized to the Einstein--de Sitter value, the :term:`CDM` limit (:math:`v_k\rightarrow0` or
   infinite lifetime) recovers the wrapped critical overdensity exactly. The decay lifetime and
   velocity kick are taken from a :cite:t:`montandon_decaying_2026` decaying dark matter particle
   (:galacticus-class:`darkMatterParticleDecayingDarkMatter`).

   .. warning::

      In the model of :cite:t:`montandon_decaying_2026` all DDM physics is encoded through the critical
      overdensity (and, for the halo mass function, a mass remapping); the variance
      :math:`\sigma(M)` is computed from the *unmodified* :math:`\Lambda`CDM linear power
      spectrum. This class must therefore be used with a standard :math:`\Lambda`CDM
      ``cosmologicalMassVariance``/``transferFunction`` --- combining it with a suppressed
      (DDM/:term:`WDM`) transfer function would double-count the small-scale suppression.
   </description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensityDecayingDarkMatter
     !!{RST
     A critical overdensity for collapse class implementing the decaying dark matter model of
     :cite:t:`montandon_decaying_2026`.
     !!}
     private
     ! Note: cosmologyFunctions_, cosmologicalMassVariance_, and linearGrowth_ are members of the base
     ! criticalOverdensityClass (used by its collapsingMass/timeOfCollapse methods), and are set here via
     ! constructorAssign rather than declared as our own members.
     class           (criticalOverdensityClass), pointer :: criticalOverdensity_ => null()
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
     double precision                                    :: lifetime                            , velocityKick
     ! Cache of epoch-dependent (but mass-independent) quantities.
     double precision                                    :: timeCached           =  -huge(0.0d0)
     double precision                                    :: deltaCLargeCached                   , deltaCSmallCached, &
          &                                                 massScale1Cached                    , deltaCEdSCached
   contains
     final     ::                    decayingDarkMatterDestructor
     procedure :: value           => decayingDarkMatterValue
     procedure :: gradientTime    => decayingDarkMatterGradientTime
     procedure :: gradientMass    => decayingDarkMatterGradientMass
     procedure :: isMassDependent => decayingDarkMatterIsMassDependent
     procedure :: isNodeDependent => decayingDarkMatterIsNodeDependent
     procedure :: isTreeDependent => decayingDarkMatterIsTreeDependent
  end type criticalOverdensityDecayingDarkMatter

  interface criticalOverdensityDecayingDarkMatter
     !!{RST
     Constructors for the :galacticus-class:`criticalOverdensityDecayingDarkMatter` critical overdensity for collapse class.
     !!}
     module procedure decayingDarkMatterConstructorParameters
     module procedure decayingDarkMatterConstructorInternal
  end interface criticalOverdensityDecayingDarkMatter

contains

  function decayingDarkMatterConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`criticalOverdensityDecayingDarkMatter` critical overdensity for collapse class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (criticalOverdensityDecayingDarkMatter)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(criticalOverdensityClass             ), pointer       :: criticalOverdensity_
    class(cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass        ), pointer       :: cosmologicalMassVariance_
    class(linearGrowthClass                    ), pointer       :: linearGrowth_
    class(darkMatterParticleClass              ), pointer       :: darkMatterParticle_

    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="darkMatterParticle"       name="darkMatterParticle_"       source="parameters"/>
    !!]
    self=criticalOverdensityDecayingDarkMatter(criticalOverdensity_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="darkMatterParticle_"      />
    !!]
    return
  end function decayingDarkMatterConstructorParameters

  function decayingDarkMatterConstructorInternal(criticalOverdensity_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,darkMatterParticle_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`criticalOverdensityDecayingDarkMatter` critical overdensity for collapse class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleDecayingDarkMatter
    use :: Error                , only : Error_Report
    implicit none
    type (criticalOverdensityDecayingDarkMatter)                        :: self
    class(criticalOverdensityClass             ), target, intent(in   ) :: criticalOverdensity_
    class(cosmologyFunctionsClass              ), target, intent(in   ) :: cosmologyFunctions_
    class(cosmologicalMassVarianceClass        ), target, intent(in   ) :: cosmologicalMassVariance_
    class(linearGrowthClass                    ), target, intent(in   ) :: linearGrowth_
    class(darkMatterParticleClass              ), target, intent(in   ) :: darkMatterParticle_
    !![
    <constructorAssign variables="*criticalOverdensity_, *cosmologyFunctions_, *cosmologicalMassVariance_, *linearGrowth_, *darkMatterParticle_"/>
    !!]

    ! Extract the decay lifetime and velocity kick from the (decaying dark matter) particle. Note that we
    ! must select on our own (assigned) pointer to the particle, as the accessor methods require an
    ! `intent(inout)` object.
    select type (particle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime    =particle_%lifetime    ()
       self%velocityKick=particle_%velocityKick()
    class default
       call Error_Report('a decaying dark matter particle ([darkMatterParticleDecayingDarkMatter]) is required'//{introspection:location})
    end select
    return
  end function decayingDarkMatterConstructorInternal

  subroutine decayingDarkMatterDestructor(self)
    !!{RST
    Destructor for the decayingDarkMatter critical overdensity for collapse class.
    !!}
    implicit none
    type(criticalOverdensityDecayingDarkMatter), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%darkMatterParticle_"      />
    !!]
    return
  end subroutine decayingDarkMatterDestructor

  double precision function decayingDarkMatterCorrection(self,timeValue,mass) result(correction)
    !!{RST
    Return the DDM correction factor :math:`\delta_\mathrm{c}^\mathrm{fit}(M_0)/\delta_\mathrm{c}^\mathrm{EdS}`
    at the given cosmic time and (Lagrangian) mass. Epoch-dependent (mass-independent) quantities are
    cached, since the halo mass function evaluates this at many masses for a single epoch.
    !!}
    use :: Decaying_Dark_Matter_Spherical_Collapse, only : decayingDarkMatterEpsilon    , decayingDarkMatterGammaTilde            , &
         &                                                 decayingDarkMatterJIntegral  , decayingDarkMatterDeltaCLarge           , &
         &                                                 decayingDarkMatterDeltaCSmall, decayingDarkMatterMassScale1            , &
         &                                                 decayingDarkMatterDeltaCFit  , decayingDarkMatterCriticalOverdensityEdS
    implicit none
    class           (criticalOverdensityDecayingDarkMatter), intent(inout) :: self
    double precision                                       , intent(in   ) :: timeValue , mass
    double precision                                                       :: gammaTilde, epsilon, jIntegral

    ! Update the cache of epoch-dependent quantities if the epoch has changed.
    if (timeValue /= self%timeCached) then
       gammaTilde            =decayingDarkMatterGammaTilde            (     timeValue   ,self%lifetime              )
       epsilon               =decayingDarkMatterEpsilon               (                  self%velocityKick          )
       jIntegral             =decayingDarkMatterJIntegral             (     gammaTilde                              )
       self%deltaCLargeCached=decayingDarkMatterDeltaCLarge           (     gammaTilde  ,     epsilon     ,jIntegral)
       self%deltaCSmallCached=decayingDarkMatterDeltaCSmall           (     gammaTilde                              )
       self%massScale1Cached =decayingDarkMatterMassScale1            (self%velocityKick,     gammaTilde  ,timeValue)
       self%deltaCEdSCached  =decayingDarkMatterCriticalOverdensityEdS(                                             )
       self%timeCached       =timeValue
    end if
    correction=+decayingDarkMatterDeltaCFit(mass,self%deltaCLargeCached,self%deltaCSmallCached,self%massScale1Cached) &
         &     /self%deltaCEdSCached
    return
  end function decayingDarkMatterCorrection

  double precision function decayingDarkMatterEpoch(self,time,expansionFactor) result(timeValue)
    !!{RST
    Resolve the cosmic time (in Gyr) from either an explicit ``time`` or an ``expansionFactor``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (criticalOverdensityDecayingDarkMatter), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: time, expansionFactor

    if      (present(time           )) then
       timeValue=time
    else if (present(expansionFactor)) then
       timeValue=self%cosmologyFunctions_%cosmicTime(expansionFactor)
    else
       timeValue=0.0d0
       call Error_Report('either time or expansionFactor must be provided'//{introspection:location})
    end if
    return
  end function decayingDarkMatterEpoch

  double precision function decayingDarkMatterValue(self,time,expansionFactor,collapsing,mass,node)
    !!{RST
    Return the DDM critical overdensity for collapse at the given time and mass.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (criticalOverdensityDecayingDarkMatter), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: time      , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    double precision                                       , intent(in   ), optional :: mass
    type            (treeNode                             ), intent(inout), optional :: node
    double precision                                                                 :: timeValue
    !$GLC attributes unused :: node

    if (.not.present(mass)) call Error_Report('mass is required for this critical overdensity class'//{introspection:location})
    timeValue              =+decayingDarkMatterEpoch(self,time,expansionFactor)
    decayingDarkMatterValue=+self%criticalOverdensity_%value                       (     time     ,expansionFactor,collapsing,mass) &
         &                  *                          decayingDarkMatterCorrection(self,timeValue                           ,mass)
    return
  end function decayingDarkMatterValue

  double precision function decayingDarkMatterGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{RST
    Return the gradient with respect to time of the DDM critical overdensity at the given time and mass.
    The derivative of the (weakly time-dependent) DDM correction factor is evaluated by central finite
    differences, and combined with the wrapped critical overdensity's own time gradient via the product
    rule.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (criticalOverdensityDecayingDarkMatter), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: time             , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    double precision                                       , intent(in   ), optional :: mass
    type            (treeNode                             ), intent(inout), optional :: node
    double precision                                       , parameter               :: fractionalStep=1.0d-4
    double precision                                                                 :: timeValue            , timeStep          , &
         &                                                                              correction           , gradientCorrection
    !$GLC attributes unused :: node

    if (.not.present(mass)) call Error_Report('mass is required for this critical overdensity class'//{introspection:location})
    timeValue         =+decayingDarkMatterEpoch(self,time,expansionFactor)
    timeStep          =+fractionalStep &
         &             *timeValue
    correction        =+  decayingDarkMatterCorrection(self,timeValue         ,mass)
    gradientCorrection=+(                                                             &
         &               +decayingDarkMatterCorrection(self,timeValue+timeStep,mass)  &
         &               -decayingDarkMatterCorrection(self,timeValue-timeStep,mass)  &
         &              )                                                             &
         &             /(2.0d0*timeStep)
    decayingDarkMatterGradientTime=+self%criticalOverdensity_%gradientTime(time,expansionFactor,collapsing,mass) &
         &                         *correction                                                                   &
         &                         +self%criticalOverdensity_%value       (time,expansionFactor,collapsing,mass) &
         &                         *gradientCorrection
    return
  end function decayingDarkMatterGradientTime

  double precision function decayingDarkMatterGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{RST
    Return the gradient with respect to mass of the DDM critical overdensity at the given time and mass.
    The derivative of the mass-dependent DDM correction factor is evaluated by central finite
    differences, and combined with the wrapped critical overdensity's own mass gradient via the product
    rule.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (criticalOverdensityDecayingDarkMatter), intent(inout)           :: self
    double precision                                       , intent(in   ), optional :: time                 , expansionFactor
    logical                                                , intent(in   ), optional :: collapsing
    double precision                                       , intent(in   ), optional :: mass
    type            (treeNode                             ), intent(inout), optional :: node
    double precision                                       , parameter               :: fractionalStep=1.0d-4
    double precision                                                                 :: timeValue            , massStep          , &
         &                                                                              correction           , gradientCorrection
    !$GLC attributes unused :: node

    if (.not.present(mass)) call Error_Report('mass is required for this critical overdensity class'//{introspection:location})
    timeValue         =+decayingDarkMatterEpoch(self,time,expansionFactor)
    massStep          =+fractionalStep &
         &             *mass
    correction        =+  decayingDarkMatterCorrection(self,timeValue,mass         )
    gradientCorrection=+(                                                            &
         &               +decayingDarkMatterCorrection(self,timeValue,mass+massStep) &
         &               -decayingDarkMatterCorrection(self,timeValue,mass-massStep) &
         &              )                                                            &
         &             /(2.0d0*massStep)
    decayingDarkMatterGradientMass=+self%criticalOverdensity_%gradientMass(time,expansionFactor,collapsing,mass) &
         &                         *correction                                                                   &
         &                         +self%criticalOverdensity_%value       (time,expansionFactor,collapsing,mass) &
         &                         *gradientCorrection
    return
  end function decayingDarkMatterGradientMass

  logical function decayingDarkMatterIsMassDependent(self)
    !!{RST
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensityDecayingDarkMatter), intent(inout) :: self
    !$GLC attributes unused :: self

    decayingDarkMatterIsMassDependent=.true.
    return
  end function decayingDarkMatterIsMassDependent

  logical function decayingDarkMatterIsNodeDependent(self)
    !!{RST
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensityDecayingDarkMatter), intent(inout) :: self

    decayingDarkMatterIsNodeDependent=self%criticalOverdensity_%isNodeDependent()
    return
  end function decayingDarkMatterIsNodeDependent

  logical function decayingDarkMatterIsTreeDependent(self)
    !!{RST
    Return whether the critical overdensity is tree dependent.
    !!}
    implicit none
    class(criticalOverdensityDecayingDarkMatter), intent(inout) :: self

    decayingDarkMatterIsTreeDependent=self%criticalOverdensity_%isTreeDependent()
    return
  end function decayingDarkMatterIsTreeDependent

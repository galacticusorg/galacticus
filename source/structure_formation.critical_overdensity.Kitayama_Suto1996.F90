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
  Implements a critical overdensity class based on the fitting functions of
  \cite{kitayama_semianalytic_1996}.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <criticalOverdensity name="criticalOverdensityKitayamaSuto1996">
   <description>
    A critical overdensity class based on the fitting functions of \cite{kitayama_semianalytic_1996}, which is therefore valid
    only for flat cosmological models (an error will be reported in non-flat models). Specifically,
    \begin{equation}
     \delta_\mathrm{crit} (t) = {3 (12 \pi)^{2/3} \over 20} [1 + 0.0123 \log_{10}\left\{\Omega_\mathrm{matter}(t)\right\}]/D(t).
    \end{equation}
   </description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensityKitayamaSuto1996
     !!{
     A critical overdensity class based on the fitting functions of \cite{kitayama_semianalytic_1996}.
     !!}
     private
     double precision                                   :: timePrevious                 , valuePrevious
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
    contains
     final     ::                    kitayamaSuto1996Destructor
     procedure :: value           => kitayamaSuto1996Value
     procedure :: gradientTime    => kitayamaSuto1996GradientTime
     procedure :: gradientMass    => kitayamaSuto1996GradientMass
     procedure :: isMassDependent => kitayamaSuto1996IsMassDependent
     procedure :: isNodeDependent => kitayamaSuto1996IsNodeDependent
  end type criticalOverdensityKitayamaSuto1996

  interface criticalOverdensityKitayamaSuto1996
     !!{
     Constructors for the \refClass{criticalOverdensityKitayamaSuto1996} critical overdensity class.
     !!}
     module procedure kitayamaSuto1996ConstructorParameters
     module procedure kitayamaSuto1996ConstructorInternal
  end interface criticalOverdensityKitayamaSuto1996

contains

  function kitayamaSuto1996ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{criticalOverdensityKitayamaSuto1996} critical overdensity class
    which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (criticalOverdensityKitayamaSuto1996)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class(linearGrowthClass                  ), pointer       :: linearGrowth_
    class(cosmologicalMassVarianceClass      ), pointer       :: cosmologicalMassVariance_
    class(darkMatterParticleClass            ), pointer       :: darkMatterParticle_

    !![
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"       name="darkMatterParticle_"       source="parameters"/>
    !!]
    self=criticalOverdensityKitayamaSuto1996(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="darkMatterParticle_"      />
    !!]
    return
  end function kitayamaSuto1996ConstructorParameters

  function kitayamaSuto1996ConstructorInternal(linearGrowth_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the \refClass{criticalOverdensityKitayamaSuto1996} critical overdensity class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleCDM, darkMatterParticleClass
    use :: Error                , only : Error_Report
    implicit none
    type (criticalOverdensityKitayamaSuto1996)                        :: self
    class(cosmologyFunctionsClass            ), target, intent(in   ) :: cosmologyFunctions_
    class(linearGrowthClass                  ), target, intent(in   ) :: linearGrowth_
    class(cosmologicalMassVarianceClass      ), target, intent(in   ) :: cosmologicalMassVariance_
    class(darkMatterParticleClass            ), target, intent(in   ) :: darkMatterParticle_
    !![
    <constructorAssign variables="*linearGrowth_, *cosmologyFunctions_, *cosmologicalMassVariance_, *darkMatterParticle_"/>
    !!]

    self%timePrevious=-1.0d0
    ! Require that the dark matter be cold dark matter.
    select type (darkMatterParticle_)
    class is (darkMatterParticleCDM)
       ! Cold dark matter particle - this is as expected.
    class default
       call Error_Report('critical overdensity expects a cold dark matter particle'//{introspection:location})
    end select
   return
  end function kitayamaSuto1996ConstructorInternal

  subroutine kitayamaSuto1996Destructor(self)
    !!{
    Destructor for the \refClass{criticalOverdensityKitayamaSuto1996} critical overdensity class.
    !!}
    implicit none
    type(criticalOverdensityKitayamaSuto1996), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%linearGrowth_"            />
    <objectDestructor name="self%darkMatterParticle_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine kitayamaSuto1996Destructor

  double precision function kitayamaSuto1996Value(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the critical overdensity at the given time and mass.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (criticalOverdensityKitayamaSuto1996), intent(inout)           :: self
    double precision                                     , intent(in   ), optional :: time               , expansionFactor
    logical                                              , intent(in   ), optional :: collapsing
    double precision                                     , intent(in   ), optional :: mass
    type            (treeNode                           ), intent(inout), optional :: node
    double precision                                                               :: time_
    !$GLC attributes unused :: mass, node

    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    if (time_ /= self%timePrevious)                                                                       &
         & self%valuePrevious=+(3.0d0*(12.0d0*Pi)**(2.0d0/3.0d0)/20.0d0)                                  &
         &                    *(1.0d0+0.0123d0*log10(self%cosmologyFunctions_%omegaMatterEpochal(time_)))
    kitayamaSuto1996Value=self%valuePrevious
    return
  end function kitayamaSuto1996Value

  double precision function kitayamaSuto1996GradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to time of critical overdensity at the given time and mass.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (criticalOverdensityKitayamaSuto1996), intent(inout)           :: self
    double precision                                     , intent(in   ), optional :: time      , expansionFactor
    logical                                              , intent(in   ), optional :: collapsing
    double precision                                     , intent(in   ), optional :: mass
    type            (treeNode                           ), intent(inout), optional :: node
    double precision                                                               :: time_     , expansionFactor_
    !$GLC attributes unused :: mass, node

    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_,expansionFactorOut=expansionFactor_)
    kitayamaSuto1996GradientTime=+(3.0d0*(12.0d0*Pi)**(2.0d0/3.0d0)/20.0d0)                          &
         &                       *(                                                                  &
         &                         +0.0123d0*self%cosmologyFunctions_%omegaMatterRateOfChange(time_) &
         &                         /         self%cosmologyFunctions_%omegaMatterEpochal     (time_) &
         &                         /log(10.0d0)  &
         &                        )
    return
  end function kitayamaSuto1996GradientTime

  double precision function kitayamaSuto1996GradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to mass of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityKitayamaSuto1996), intent(inout)           :: self
    double precision                                     , intent(in   ), optional :: time      , expansionFactor
    logical                                              , intent(in   ), optional :: collapsing
    double precision                                     , intent(in   ), optional :: mass
    type            (treeNode                           ), intent(inout), optional :: node
    !$GLC attributes unused :: self, time, expansionFactor, collapsing, mass, node

    kitayamaSuto1996GradientMass=0.0d0
    return
  end function kitayamaSuto1996GradientMass

  logical function kitayamaSuto1996IsMassDependent(self)
    !!{
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensityKitayamaSuto1996), intent(inout) :: self
    !$GLC attributes unused :: self

    kitayamaSuto1996IsMassDependent=.false.
    return
  end function kitayamaSuto1996IsMassDependent

  logical function kitayamaSuto1996IsNodeDependent(self)
    !!{
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensityKitayamaSuto1996), intent(inout) :: self
    !$GLC attributes unused :: self

    kitayamaSuto1996IsNodeDependent=.false.
    return
  end function kitayamaSuto1996IsNodeDependent

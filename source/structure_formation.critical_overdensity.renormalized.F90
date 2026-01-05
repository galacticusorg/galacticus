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
Implements a critical overdensity class which renormalizes another class based on the ratio of two mass variance classes. This is
intended to allow different window functions to be used for $\sigma(M)$ while retaining the same ratio
$\delta_\mathrm{c}/\sigma(M)$ (and, therefore, the same halo mass function) on a mass scale $M_\mathrm{match}$.
!!}

  !![
  <criticalOverdensity name="criticalOverdensityRenormalize">
    <description>
      A critical overdensity class which renormalizes another class based on the ratio of two mass variance classes. This is
      intended to allow different window functions to be used for $\sigma(M)$ while retaining the same ratio
      $\delta_\mathrm{c}/\sigma(M)$ (and, therefore, the same halo mass function) on a mass scale $M_\mathrm{match}$.

      Specifically, the renormalized critical overdensity, $\delta^\prime_\mathrm{c}$, is given by $\delta^\prime_\mathrm{c} = n
      \delta_\mathrm{c}$, where $\delta_\mathrm{c}$ is the original critical overdensity, and      
      \begin{equation}
      n=\sigma(M_\mathrm{match})/\sigma_\mathrm{r}(M_\mathrm{match}),
      \end{equation}      
      with $\sigma(M)$ being the cosmological mass variance (computed with whatever window function is required),
      $\sigma_\mathrm{r}(M)$ is the reference cosmological mass variance (typically computed using a top-hat window function), and
      $M_\mathrm{match}$ is the mass at which to match the mass variance. The mass variances are evaluated at the present epoch.

      The matching scale is given by $M_\mathrm{match}=${\normalfont \ttfamily [massMatch]} if this parameter is
      present. Otherwise $M_\mathrm{match}=M_*$ is used, where $\sigma(M_*)=\delta_\mathrm{crit}$, computed at the present epoch.
   </description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensityRenormalize
     !!{
     A critical overdensity class which renormalizes another class based on the ratio of two mass variance classes. This is
     intended to allow different window functions to be used for $\sigma(M)$ while retaining the same ratio
     $\delta_\mathrm{c}/\sigma(M)$ (and, therefore, the same halo mass function) on a mass scale $M_\mathrm{match}$.
     !!}
     private
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_               => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVarianceReference_ => null()
     double precision                                         :: massMatch                                   , normalization
    contains
     final     ::                    renormalizeDestructor
     procedure :: value           => renormalizeValue
     procedure :: gradientTime    => renormalizeGradientTime
     procedure :: gradientMass    => renormalizeGradientMass
     procedure :: isMassDependent => renormalizeIsMassDependent
     procedure :: isNodeDependent => renormalizeIsNodeDependent
     procedure :: isTreeDependent => renormalizeIsTreeDependent
  end type criticalOverdensityRenormalize

  interface criticalOverdensityRenormalize
     !!{
     Constructors for the \refClass{criticalOverdensityRenormalize} critical overdensity class.
     !!}
     module procedure renormalizeConstructorParameters
     module procedure renormalizeConstructorInternal
  end interface criticalOverdensityRenormalize

contains

  function renormalizeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{criticalOverdensityRenormalize} critical overdensity class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (criticalOverdensityRenormalize)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (criticalOverdensityClass      ), pointer       :: criticalOverdensity_
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass ), pointer       :: cosmologicalMassVariance_, cosmologicalMassVarianceReference_
    class           (linearGrowthClass             ), pointer       :: linearGrowth_
    double precision                                                :: massMatch
    
    ! Check and read parameters.
    if (parameters%isPresent('massMatch')) then
       !![
       <inputParameter>
	 <name>massMatch</name>
	 <source>parameters</source>
	 <description>The mass scale at which to renormalize.</description>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"                                                                 source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"                                                                  source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_"                                                            source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVarianceReference_" parameterName="cosmologicalMassVarianceReference" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"                                                                        source="parameters"/>
    <conditionalCall>
      <call>self=criticalOverdensityRenormalize(criticalOverdensity_,cosmologyFunctions_,cosmologicalMassVariance_,cosmologicalMassVarianceReference_,linearGrowth_{conditions})</call>
      <argument name="massMatch" value="massMatch" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"              />
    <objectDestructor name="cosmologyFunctions_"               />
    <objectDestructor name="cosmologicalMassVariance_"         />
    <objectDestructor name="cosmologicalMassVarianceReference_"/>
    <objectDestructor name="linearGrowth_"                     />
    !!]
    return
  end function renormalizeConstructorParameters

  function renormalizeConstructorInternal(criticalOverdensity_,cosmologyFunctions_,cosmologicalMassVariance_,cosmologicalMassVarianceReference_,linearGrowth_,massMatch) result(self)
    !!{
    Internal constructor for the \refClass{criticalOverdensityRenormalize} critical overdensity class.
    !!}
    implicit none
    type            (criticalOverdensityRenormalize)                          :: self
    class           (criticalOverdensityClass      ), target  , intent(in   ) :: criticalOverdensity_
    class           (cosmologyFunctionsClass       ), target  , intent(inout) :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass ), target  , intent(inout) :: cosmologicalMassVariance_, cosmologicalMassVarianceReference_
    class           (linearGrowthClass             ), target  , intent(in   ) :: linearGrowth_
    double precision                                , optional, intent(in   ) :: massMatch
    double precision                                                          :: massMatch_
    !![
    <constructorAssign variables="massMatch, *criticalOverdensity_, *cosmologyFunctions_, *cosmologicalMassVariance_, *cosmologicalMassVarianceReference_, *linearGrowth_"/>
    !!]

    ! Determine the mass scale at which to renormalize.
    if (present(massMatch)) then
       massMatch_=massMatch
    else
       massMatch_=self%criticalOverdensity_%collapsingMass(expansionFactor=1.0d0)
    end if
    ! Compute the normalization factor.
    self%normalization=+cosmologicalMassVariance_         %rootVariance(massMatch_,cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)) &
         &             /cosmologicalMassVarianceReference_%rootVariance(massMatch_,cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0))
    return
  end function renormalizeConstructorInternal

  subroutine renormalizeDestructor(self)
    !!{
    Destructor for the \refClass{criticalOverdensityRenormalize} critical overdensity class.
    !!}
    implicit none
    type(criticalOverdensityRenormalize), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"              />
    <objectDestructor name="self%cosmologyFunctions_"               />
    <objectDestructor name="self%cosmologicalMassVariance_"         />
    <objectDestructor name="self%cosmologicalMassVarianceReference_"/>
    <objectDestructor name="self%linearGrowth_"                     />
    !!]
    return
  end subroutine renormalizeDestructor

  double precision function renormalizeValue(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the critical overdensity for collapse at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityRenormalize), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: time      , expansionFactor
    logical                                         , intent(in   ), optional :: collapsing
    double precision                                , intent(in   ), optional :: mass
    type            (treeNode                      ), intent(inout), optional :: node

    renormalizeValue=+self%criticalOverdensity_%value(time,expansionFactor,collapsing,mass,node) &
         &           *self%normalization
    return
  end function renormalizeValue

  double precision function renormalizeGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to time of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityRenormalize), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: time      , expansionFactor
    logical                                         , intent(in   ), optional :: collapsing
    double precision                                , intent(in   ), optional :: mass
    type            (treeNode                      ), intent(inout), optional :: node

    renormalizeGradientTime=+self%criticalOverdensity_%gradientTime(time,expansionFactor,collapsing,mass,node) &
         &                  *self%normalization
    return
  end function renormalizeGradientTime

  double precision function renormalizeGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to mass of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityRenormalize), intent(inout)           :: self
    double precision                                , intent(in   ), optional :: time      , expansionFactor
    logical                                         , intent(in   ), optional :: collapsing
    double precision                                , intent(in   ), optional :: mass
    type            (treeNode                      ), intent(inout), optional :: node

    renormalizeGradientMass=+self%criticalOverdensity_%gradientMass(time,expansionFactor,collapsing,mass,node) &
         &                  *self%normalization
    return
  end function renormalizeGradientMass

  logical function renormalizeIsMassDependent(self)
    !!{
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensityRenormalize), intent(inout) :: self

    renormalizeIsMassDependent=self%criticalOverdensity_%isMassDependent()
    return
  end function renormalizeIsMassDependent

  logical function renormalizeIsNodeDependent(self)
    !!{
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensityRenormalize), intent(inout) :: self

    renormalizeIsNodeDependent=self%criticalOverdensity_%isNodeDependent()
    return
  end function renormalizeIsNodeDependent

  logical function renormalizeIsTreeDependent(self)
    !!{
    Return whether the critical overdensity is tree dependent.
    !!}
    implicit none
    class(criticalOverdensityRenormalize), intent(inout) :: self

    renormalizeIsTreeDependent=self%criticalOverdensity_%isTreeDependent()
    return
  end function renormalizeIsTreeDependent

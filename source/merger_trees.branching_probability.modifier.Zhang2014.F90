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
  Implements a merger tree branching probability rate modifier which uses the model of :cite:t:`zhang_dark-matter_2014` to account for environmental-dependence.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass, haloEnvironmentClass
  use :: Kind_Numbers              , only : kind_int8
  
  !![
  <mergerTreeBranchingProbabilityModifier name="mergerTreeBranchingProbabilityModifierZhang2014" docformat="rst">
   <description>
   Provides a merger tree branching probability rate modifier which uses the model of :cite:t:`zhang_dark-matter_2014` to account for envionmental-dependence of the merger rate. Specifically, the modifier is:

   .. math::

      1 + \kappa \beta \alpha - \sqrt{2\pi} \kappa \nu (1-\alpha)^{3/2} + \pi \kappa \left(\nu^2\frac{\delta_\mathrm{d}}{\delta_\mathrm{c}}-1\right) (1-\alpha)^{3/2} \exp\left[\frac{\nu^2}{2}\right] \hbox{erfc}\left[\frac{\nu}{\sqrt{2}}\right]

   where :math:`\kappa` characterizes the non-Markovian behavior in the excursion set random walks, and whose value depends on the shape of the window function :cite:p:`zhang_dark-matter_2014`, :math:`\alpha = S_\mathrm{d}/S_\mathrm{p}` with :math:`S_\mathrm{d}` being the variance on the scale of the descendant (i.e., parent) halo, and :math:`S_\mathrm{p}` being the variance on the scale of the progenitor (i.e., child) halo,

   .. math::

      \beta = -2 + \frac{(1-\alpha)^{3/2}}{2\alpha}\ln\left(\frac{1+\sqrt{1-\alpha}}{1-\sqrt{1-\alpha}}\right)+\frac{1}{\alpha}+2\alpha,

   :math:`\nu = \delta_\mathrm{c}/\sqrt{S_\mathrm{d}}`, :math:`\delta_\mathrm{e}` is the overdensity on the scale of the environment (evaluated at :math:`z=0`), and :math:`\delta_\mathrm{c}` is the critical overdensity at the time of the descendant halo.

   Note that this implementation uses the solutions for the limit of large environment scale (:math:`S_\mathrm{e} \rightarrow 0`; equation 30 of :cite:author:`zhang_dark-matter_2014` :cite:year:`zhang_dark-matter_2014`). A warning is emitted if :math:`S_\mathrm{e}` is sufficiently large to break this assumption.
   </description>
  </mergerTreeBranchingProbabilityModifier>
  !!]
  type, extends(mergerTreeBranchingProbabilityModifierClass) :: mergerTreeBranchingProbabilityModifierZhang2014
     !!{RST
     A merger tree branching probability rate modifier which uses the model of :cite:t:`zhang_dark-matter_2014` to account for environmental-dependence.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (haloEnvironmentClass         ), pointer :: haloEnvironment_          => null()
     double precision                                         :: kappa                              , overdensityEnvironment, &
          &                                                      overdensityCriticalParent          , timeParent
     integer         (kind_int8                    )          :: uniqueID
   contains
     final     ::                 zhang2014Destructor
     procedure :: rateModifier => zhang2014RateModifier
  end type mergerTreeBranchingProbabilityModifierZhang2014

  interface mergerTreeBranchingProbabilityModifierZhang2014
     !!{RST
     Constructors for the :galacticus-class:`mergerTreeBranchingProbabilityModifierZhang2014` merger tree branching probability rate class.
     !!}
     module procedure zhang2014ConstructorParameters
     module procedure zhang2014ConstructorInternal
  end interface mergerTreeBranchingProbabilityModifierZhang2014

contains

  function zhang2014ConstructorParameters(parameters) result(self)
    !!{RST
    A constructor for the ``zhang2014`` merger tree branching probability rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBranchingProbabilityModifierZhang2014)                :: self
    type            (inputParameters                                ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                        ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass                  ), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass                       ), pointer       :: criticalOverdensity_
    class           (haloEnvironmentClass                           ), pointer       :: haloEnvironment_
    double precision                                                                 :: kappa

    !![
    <inputParameter docformat="rst">
      <name>kappa</name>
      <defaultValue>0.44d0</defaultValue>
      <defaultSource>
      :cite:p:`zhang_dark-matter_2014`
      </defaultSource>
      <description>
      The parameter characterizing non-Markovian behavior in excursion set random walks :cite:p:`zhang_dark-matter_2014`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="haloEnvironment"          name="haloEnvironment_"          source="parameters"/>
    !!]
    self=mergerTreeBranchingProbabilityModifierZhang2014(kappa,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,haloEnvironment_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="haloEnvironment_"         />
    !!]
    return
  end function zhang2014ConstructorParameters

  function zhang2014ConstructorInternal(kappa,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,haloEnvironment_) result(self)
    !!{RST
    Default constructor for the ``zhang2014`` merger tree branching probability rate class.
    !!}
    use :: Error, only : Warn
    implicit none
    type            (mergerTreeBranchingProbabilityModifierZhang2014)                        :: self
    class           (cosmologyFunctionsClass                        ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass                  ), intent(in   ), target :: cosmologicalMassVariance_
    class           (criticalOverdensityClass                       ), intent(in   ), target :: criticalOverdensity_
    class           (haloEnvironmentClass                           ), intent(in   ), target :: haloEnvironment_
    double precision                                                 , intent(in   )         :: kappa
    double precision                                                 , parameter             :: varianceLarge            =0.1d0
    double precision                                                                         :: massEnvironment                , timePresent, &
         &                                                                                      varianceEnvironment
    character       (len= 5                                         )                        :: labelVariance
    character       (len=12                                         )                        :: labelMass
    !![
    <constructorAssign variables="kappa, *cosmologyFunctions_, *cosmologicalMassVariance_, *criticalOverdensity_, *haloEnvironment_"/>
    !!]

    ! Initialization memoization state.
    self%uniqueID=-1_kind_int8
    ! Validate the environment scale.
    massEnvironment    =self%haloEnvironment_         %environmentMass(                                                )
    timePresent        =self%cosmologyFunctions_      %cosmicTime     (                     expansionFactor=1.0d0      )
    varianceEnvironment=self%cosmologicalMassVariance_%rootVariance   (mass=massEnvironment,time           =timePresent)**2
    if (varianceEnvironment > varianceLarge) then
       write (labelVariance,'(f05.2)') varianceEnvironment
       write (labelMass    ,'(e12.6)') massEnvironment
       call Warn('this class assume the environment variance Sₑ ≪ 1, but found Sₑ = '//labelVariance//' [Mₑ = '//labelMass//' M☉]')
    end if
    return
  end function zhang2014ConstructorInternal

  subroutine zhang2014Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`mergerTreeBranchingProbabilityModifierZhang2014` class.
    !!}
    implicit none
    type(mergerTreeBranchingProbabilityModifierZhang2014), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%haloEnvironment_"         />
    !!]
    return
  end subroutine zhang2014Destructor

  double precision function zhang2014RateModifier(self,nodeParent,massParent,sigmaParent,sigmaChild,timeParent) result(modifier)
    !!{RST
    Returns a modifier for merger tree branching rates using the model of :cite:t:`zhang_dark-matter_2014` to account for environmental-dependence.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Galacticus_Nodes        , only : nodeComponentBasic
    implicit none
    class           (mergerTreeBranchingProbabilityModifierZhang2014), intent(inout) :: self
    type            (treeNode                                       ), intent(inout) :: nodeParent
    double precision                                                 , intent(in   ) :: sigmaChild , timeParent, &
         &                                                                              sigmaParent, massParent
    class           (nodeComponentBasic                             ), pointer       :: basicParent
    double precision                                                                 :: peakHeight , alpha     , &
         &                                                                              beta       , wParent
    
    ! Evaluate the environmental and critical overdensities.
    if (nodeParent%uniqueID() /= self%uniqueID .or. timeParent /= self%timeParent) then
       basicParent => nodeParent%basic()
       ! Store the w "time" parameter used in merger tree building.
       wParent=basicParent%time()
       ! Replace with the actual time.
       call basicParent%timeSet(timeParent)
       ! Compute environmental and critical overdensities.
       self%overdensityEnvironment   =self%haloEnvironment_    %overdensityLinear(                                     nodeParent)
       self%overdensityCriticalParent=self%criticalOverdensity_%value            (time=timeParent,mass=massParent,node=nodeParent)
       ! Restore the w time parameter.
       call basicParent%timeSet(wParent)
       ! Record memoization state.
       self%uniqueID  =nodeParent%uniqueID()
       self%timeParent=timeParent
    end if
    ! Evaluate other parameters needed for the rate modifier.
    peakHeight=+self%overdensityCriticalParent  &
         &     /sqrt(sigmaParent              )    
    alpha     =(             &
         &      +sigmaParent &
         &      /sigmaChild  &
         &     )**2    
    beta      =- 2.0d0                         &
         &     +(1.0d0-alpha)**1.5d0           &
         &     / 2.0d0/alpha                   &
         &     *log(                           &
         &          +(1.0d0+sqrt(1.0d0-alpha)) &
         &          /(1.0d0-sqrt(1.0d0-alpha)) &
         &         )                           &
         &     +1.0d0/alpha                    &
         &     +2.0d0*alpha
    ! Evaluate the rate modifier function.
    modifier=+1.0d0                            &
         &   +self%kappa                       &
         &   *beta                             &
         &   *alpha                            &
         &   -sqrt(                            &
         &         +2.0d0                      &
         &         *Pi                         &
         &        )                            &
         &   *self%kappa                       &
         &   *peakHeight                       &
         &   *(1.0d0-alpha)**1.5d0             &
         &   +Pi                               &
         &   *self%kappa                       &
         &   *(                                &
         &     +     peakHeight**2             &
         &     *self%overdensityEnvironment    &
         &     /self%overdensityCriticalParent &
         &     -1.0d0                          &
         &    )                                &
         &   *(1.0d0-alpha)**1.5d0             &
         &   *exp (peakHeight**2/     2.0d0 )  & 
         &   *erfc(peakHeight   /sqrt(2.0d0))
    return
  end function zhang2014RateModifier

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
  Implementation of a cooling rate class in which the cooling rate scales with the peak circular velocity in the halo.
  !!}

  use :: Cosmology_Functions     , only : cosmologyFunctions       , cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Math_Exponentiation     , only : fastExponentiator

  !![
  <coolingRate name="coolingRateVelocityMaximumScaling">
   <description>
    A cooling rate class in which the cooling rate scales with the peak circular velocity in the halo. Specifically, the
    cooling rate is given by
    \begin{equation}
    \dot{M}_\mathrm{cool} = M_\mathrm{hot}/\tau_\mathrm{cool}(V_\mathrm{max,halo},z) ,
    \end{equation}
    where 
    \begin{equation}
    \tau_\mathrm{cool}=\hbox{max}\left[ \tau_\mathrm{infall} \left({V_\mathrm{max} \over 200
    \hbox{km/s}}\right)^{-\gamma_\mathrm{infall}} (1+z)^{\alpha_\mathrm{infall}} \left( 1 + \exp\left[
    {\log_{10}(V_\mathrm{max}/(1+z)^{\delta_\mathrm{infall}}\mathcal{V}_\mathrm{infall})] \over \Delta \log_{10}
    \mathcal{V}_\mathrm{infall}}\right]\right)^{\beta_\mathrm{infall}}, \tau_\mathrm{infall,min} \right],
    \end{equation}
    with $\tau_\mathrm{infall}=${\normalfont \ttfamily [timescale]}, $\tau_\mathrm{infall,min}=${\normalfont \ttfamily
    [timescaleMinimum]}, $\alpha_\mathrm{infall}=${\normalfont \ttfamily [exponentRedshift]},
    $\beta_\mathrm{infall}=${\normalfont \ttfamily [exponentCutOff]}, $\gamma_\mathrm{infall}=${\normalfont \ttfamily
    [exponentVelocity]}, $\delta_\mathrm{infall}=${\normalfont \ttfamily [velocityCutOffExponentRedshift},
    $\mathcal{V}_\mathrm{infall}=${\normalfont \ttfamily [velocityCutOff]}, and $\Delta \log_{10}
    \mathcal{V}_\mathrm{infall}=${\normalfont \ttfamily [widthCutOff]}.
   </description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateVelocityMaximumScaling
     !!{
     Implementation of cooling rate class in which the cooling rate scales with the peak circular velocity in the halo.
     !!}
     private
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_   => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     ! Parameters controlling the cooling rate.
     double precision                                   :: timescale                         , exponentRedshift              , &
          &                                                widthCutOff                       , velocityCutOff                , &
          &                                                exponentCutOff                    , exponentVelocity              , &
          &                                                timescaleMinimum                  , velocityCutOffExponentRedshift, &
          &                                                normalization
     ! Fast exponentiation tables for rapid computation of the cooling rate.
     type            (fastExponentiator      )          :: exponentiatorVelocity             , exponentiatorExpansionFactor  , &
          &                                                exponentiatorExpansionFactorCutOff
     ! Stored values for rapid re-use.
     double precision                                   :: expansionFactorPrevious           , velocityMaximumPrevious       , &
          &                                                coolingRateStored
   contains
     final     ::         velocityMaximumScalingDestructor
     procedure :: rate => velocityMaximumScalingRate
  end type coolingRateVelocityMaximumScaling

  interface coolingRateVelocityMaximumScaling
     !!{
     Constructors for the velocity maximum scaling cooling rate class.
     !!}
     module procedure velocityMaximumScalingConstructorParameters
     module procedure velocityMaximumScalingConstructorInternal
  end interface coolingRateVelocityMaximumScaling

contains

  function velocityMaximumScalingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the velocity maximum scaling cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingRateVelocityMaximumScaling)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass        ), pointer       :: darkMatterProfileDMO_
    double precision                                                   :: timeScale            , timescaleMinimum              , &
         &                                                                exponentRedshift     , exponentVelocity              , &
         &                                                                velocityCutOff       , velocityCutOffExponentRedshift, &
         &                                                                widthCutOff          , exponentCutOff

    !![
    <inputParameter>
      <name>timescale</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The timescale (in Gyr) for cooling in low mass halos at $z=0$ in the velocity maximum scaling scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>timescaleMinimum</name>
      <source>parameters</source>
      <defaultValue>0.001d0</defaultValue>
      <description>The minimum timescale (in Gyr) for cooling the velocity maximum scaling scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentRedshift</name>
      <source>parameters</source>
      <defaultValue>-1.5d0</defaultValue>
      <description>The exponent of $(1+z)$ in the cooling timescale for low mass halos in the velocity maximum scaling scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentVelocity</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of velocity in the cooling timescale for low mass halos in the velocity maximum scaling scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityCutOff</name>
      <source>parameters</source>
      <defaultValue>200.0d0</defaultValue>
      <description>The halo maximum velocity scale appearing in the exponential term for cooling timescale in the velocity maximum scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityCutOffExponentRedshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of $(1+z)$ in the velocity scale appearing in the exponential term for cooling timescale in the velocity maximum scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>widthCutOff</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The width appearing in the exponential term for cooling timescale in the velocity maximum scaling scaling cooling rate model.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentCutOff</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The exponent appearing in the exponential term for cooling timescale in the velocity maximum scaling scaling cooling rate model.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=coolingRateVelocityMaximumScaling(timeScale,timescaleMinimum,exponentRedshift,exponentVelocity,velocityCutOff,velocityCutOffExponentRedshift,widthCutOff,exponentCutOff,cosmologyFunctions_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function velocityMaximumScalingConstructorParameters

  function velocityMaximumScalingConstructorInternal(timeScale,timescaleMinimum,exponentRedshift,exponentVelocity,velocityCutOff,velocityCutOffExponentRedshift,widthCutOff,exponentCutOff,cosmologyFunctions_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the velocity maximum scaling cooling rate class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List          , Error_Report
    use :: Galacticus_Nodes, only : defaultBasicComponent   , defaultHotHaloComponent
    implicit none
    type            (coolingRateVelocityMaximumScaling)                        :: self
    double precision                                   , intent(in   )         :: timeScale                    , timescaleMinimum              , &
         &                                                                        exponentRedshift             , exponentVelocity              , &
         &                                                                        velocityCutOff               , velocityCutOffExponentRedshift, &
         &                                                                        widthCutOff                  , exponentCutOff
    class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                   , parameter             :: velocityNormalization=200.0d0
    !![
    <constructorAssign variables="timeScale, timescaleMinimum, exponentRedshift, exponentVelocity, velocityCutOff, velocityCutOffExponentRedshift, widthCutOff, exponentCutOff, *cosmologyFunctions_, *darkMatterProfileDMO_"/>
    !!]

    ! Check that the properties we need are gettable.
    if (.not.defaultHotHaloComponent%massIsGettable())                                                            &
         & call Error_Report(                                                                                     &
         &                   'Hot halo component must have gettable mass.'                                     // &
         &                   Component_List(                                                                      &
         &                                  'hotHalo'                                                          ,  &
         &                                   defaultHotHaloComponent%massAttributeMatch(requireGettable=.true.)   &
         &                                 )                                                                   // &
         &                   {introspection:location}                                                             &
         &                  )
    if     (                                                                                                      &
         &  .not.(                                                                                                &
         &         defaultBasicComponent%massIsGettable()                                                         &
         &        .and.                                                                                           &
         &         defaultBasicComponent%timeIsGettable()                                                         &
         &       )                                                                                                &
         & ) call Error_Report(                                                                                   &
         &                     'Basic component must have gettable mass and time.'//                              &
         &                     Component_List(                                                                    &
         &                                    'basic'                                                          ,  &
         &                                     defaultBasicComponent%massAttributeMatch(requireGettable=.true.)   &
         &                                    .intersection.                                                      &
         &                                     defaultBasicComponent%timeAttributeMatch(requireGettable=.true.)   &
         &                                   )                                                                 // &
         &                     {introspection:location}                                                           &
         &                    )
    ! Compute normalization.
    self%normalization=+1.0d0                                        &
         &             /self%timescale                               &
         &             /velocityNormalization**self%exponentVelocity
    ! Initialize exponentiators.
    self%exponentiatorVelocity             =fastExponentiator(1.0d+0,1.0d+3,self%exponentVelocity              ,1.0d+1,abortOutsideRange=.false.)
    self%exponentiatorExpansionFactor      =fastExponentiator(1.0d-3,1.0d+0,self%exponentRedshift              ,1.0d+3,abortOutsideRange=.false.)
    self%exponentiatorExpansionFactorCutOff=fastExponentiator(1.0d-3,1.0d+0,self%velocityCutOffExponentRedshift,1.0d+3,abortOutsideRange=.false.)
    ! Initialize stored solutions.
    self%expansionFactorPrevious=-1.0d0
    self%velocityMaximumPrevious=-1.0d0
    self%coolingRateStored      =-1.0d0
    return
  end function velocityMaximumScalingConstructorInternal

  subroutine velocityMaximumScalingDestructor(self)
    !!{
    Destructor for the velocity maximum scaling cooling rate class.
    !!}
    implicit none
    type(coolingRateVelocityMaximumScaling), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine velocityMaximumScalingDestructor

  double precision function velocityMaximumScalingRate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for a model in which this rate scales with the maximum circular velocity of the halo.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic   , nodeComponentHotHalo, treeNode
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (coolingRateVelocityMaximumScaling), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: node
    double precision                                   , parameter     :: expArgumentMaximum=100.0d0
    class           (nodeComponentBasic               ), pointer       :: basic
    class           (nodeComponentHotHalo             ), pointer       :: hotHalo
    class           (massDistributionClass            ), pointer       :: massDistribution_
    double precision                                                   :: expFactor                 , expansionFactor, &
         &                                                                expArgument               , velocityMaximum

    ! Compute expansion factor and maximum velocity.
    massDistribution_ => self                    %darkMatterProfileDMO_%get                         (node        )
    basic             => node                                          %basic                       (            )
    hotHalo           => node                                          %hotHalo                     (            )
    expansionFactor   =  self%cosmologyFunctions_                      %expansionFactor             (basic%time())
    velocityMaximum   =  massDistribution_                             %velocityRotationCurveMaximum(            )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (expansionFactor /= self%expansionFactorPrevious .or. velocityMaximum /= self%velocityMaximumPrevious) then
       expArgument=log10(                                                                       &
            &            +velocityMaximum                                                       &
            &            /self%velocityCutOff                                                   &
            &            *self%exponentiatorExpansionFactorCutOff%exponentiate(expansionFactor) &
            &           )                                                                       &
            &      /self%widthCutOff
       if (expArgument < expArgumentMaximum) then
          expFactor=(1.0d0+exp(+expArgument))**(-self%exponentCutOff)
       else
          expFactor=       exp(+expArgument   *(-self%exponentCutOff))
       end if
       self%coolingRateStored=+self%normalization                                              &
            &                 *self%exponentiatorExpansionFactor%exponentiate(expansionFactor) &
            &                 *self%exponentiatorVelocity       %exponentiate(velocityMaximum) &
            &                 *expFactor
       self%coolingRateStored=min(                                                        &
            &                     self%coolingRateStored                                , &
            &                     1.0d0/self%timescaleMinimum                             &
            &                    )
       self%expansionFactorPrevious=expansionFactor
       self%velocityMaximumPrevious=velocityMaximum
    end if
    velocityMaximumScalingRate=hotHalo%mass()*self%coolingRateStored
    return
  end function velocityMaximumScalingRate


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
Implements a critical overdensity for collapse the \gls{wdm} modifier of \cite{barkana_constraints_2001}.
!!}
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Dark_Matter_Particles  , only : darkMatterParticleClass
  use :: Numerical_Interpolation, only : interpolator

  !![
  <criticalOverdensity name="criticalOverdensityBarkana2001WDM" defaultThreaedPrivate="yes">
   <description>
  !!]

  !![
    A critical overdensity for collapse class based on the \gls{wdm} modifier of \cite{barkana_constraints_2001} applied to
    some other, \gls{cdm} critical overdensity class. Specifically, the critical overdensity is multiplied by a factor
    \begin{equation}
     \exp\left[ \left({M_\mathrm{J} \over 8 M}\right)^{1.40}+\left({M_\mathrm{J} \over 8 M}\right)^{0.45} \right],
    \end{equation}
    where $M$ is the mass in question, $M_\mathrm{J}$ is the effective Jeans mass of the warm dark matter as defined by
    \citeauthor{barkana_constraints_2001}~[\citeyear{barkana_constraints_2001}; their eqn.~10]:
    \begin{equation}
    M_\mathrm{J} = 3.06 \times 10^8 \left( {1+z_\mathrm{eq} \over 3000}\right)^{1.5} \left({\Omega_\mathrm{M} h_0^2 \over
    0.15}\right)^{1/2} \left({g_\mathrm{X} \over 1.5} \right)^{-1} \left({m_\mathrm{X}/1.0\hbox{keV}}\right)^{-4},
    \end{equation}
    the redshift of matter-radiation equality is given by
    \begin{equation}
    z_\mathrm{eq} = 3600 \left({\Omega_\mathrm{M} h_0^2 \over 0.15}\right)-1,
    \end{equation}
    and $g_\mathrm{X}$ and $m_\mathrm{X}$ are the effective number of degrees of freedom and the mass of the warm dark matter
    particle respectively. This fitting function has been found the fit the numerical results of
    \cite{barkana_constraints_2001} well.
   </description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensityBarkana2001WDM
     !!{
     A critical overdensity for collapse class which modifies another transfer function using the \gls{wdm} modifier of \cite{barkana_constraints_2001}.
     !!}
     private
     class           (criticalOverdensityClass), pointer                   :: criticalOverdensityCDM => null()
     class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_   => null()
     class           (darkMatterParticleClass ), pointer                   :: darkMatterParticle_    => null()
     logical                                                               :: useFittingFunction
     double precision                                                      :: jeansMass
     integer                                                               :: deltaTableCount
     double precision                          , allocatable, dimension(:) :: deltaTableDelta                 , deltaTableMass
     type            (interpolator            )                            :: interpolator_
   contains
     final     ::                    barkana2001WDMDestructor
     procedure :: value           => barkana2001WDMValue
     procedure :: gradientTime    => barkana2001WDMGradientTime
     procedure :: gradientMass    => barkana2001WDMGradientMass
     procedure :: isMassDependent => barkana2001WDMIsMassDependent
     procedure :: isNodeDependent => barkana2001WDMIsNodeDependent
     procedure :: isTreeDependent => barkana2001WDMIsTreeDependent
  end type criticalOverdensityBarkana2001WDM

  interface criticalOverdensityBarkana2001WDM
     !!{
     Constructors for the \refClass{criticalOverdensityBarkana2001WDM} critical overdensity for collapse class.
     !!}
     module procedure barkana2001WDMConstructorParameters
     module procedure barkana2001WDMConstructorInternal
  end interface criticalOverdensityBarkana2001WDM

  ! Parameters of fitting functions.
  double precision, parameter :: fitParameterA            =+2.40000d0
  double precision, parameter :: fitParameterB            =+0.10000d0
  double precision, parameter :: fitParameterC            =+0.04000d0
  double precision, parameter :: fitParameterD            =+2.30000d0
  double precision, parameter :: fitParameterE            =+0.31687d0
  double precision, parameter :: fitParameterF            =+0.80900d0
  double precision, parameter :: smallMassLogarithmicSlope=-1.93400d0

contains

  function barkana2001WDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{criticalOverdensityBarkana2001WDM} critical overdensity for collapse class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (criticalOverdensityBarkana2001WDM)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (criticalOverdensityClass         ), pointer       :: criticalOverdensityCDM
    class           (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass    ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterParticleClass          ), pointer       :: darkMatterParticle_
    class           (linearGrowthClass                ), pointer       :: linearGrowth_
    logical                                                            :: useFittingFunction

    !![
    <inputParameter>
      <name>useFittingFunction</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether the warm dark matter critical overdensity mass scaling should be computed from a fitting function or from tabulated data.</description>
    </inputParameter>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensityCDM"    source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"       name="darkMatterParticle_"       source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !!]
    ! Call the internal constructor
    self=barkana2001WDMConstructorInternal(criticalOverdensityCDM,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,linearGrowth_,useFittingFunction)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensityCDM"   />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="darkMatterParticle_"      />
    <objectDestructor name="linearGrowth_"            />
    !!]
    return
  end function barkana2001WDMConstructorParameters

  function barkana2001WDMConstructorInternal(criticalOverdensityCDM,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,linearGrowth_,useFittingFunction) result(self)
    !!{
    Internal constructor for the \refClass{criticalOverdensityBarkana2001WDM} critical overdensity for collapse class.
    !!}
    use :: Cosmology_Parameters   , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles  , only : darkMatterParticleWDMThermal
    use :: FoX_DOM                , only : destroy                     , node
    use :: Error                  , only : Error_Report
    use :: Input_Paths            , only : inputPath                   , pathTypeDataStatic
    use :: IO_XML                 , only : XML_Array_Read              , XML_Get_First_Element_By_Tag_Name, XML_Parse
    use :: Numerical_Interpolation, only : GSL_Interp_CSpline
    use :: Table_Labels           , only : extrapolationTypeFix
    implicit none
    type            (criticalOverdensityBarkana2001WDM)                        :: self
    class           (criticalOverdensityClass         ), target, intent(in   ) :: criticalOverdensityCDM
    class           (cosmologyParametersClass         ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass          ), target, intent(in   ) :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass    ), target, intent(in   ) :: cosmologicalMassVariance_
    class           (darkMatterParticleClass          ), target, intent(in   ) :: darkMatterParticle_
    class           (linearGrowthClass                ), target, intent(in   ) :: linearGrowth_
    logical                                                    , intent(in   ) :: useFittingFunction
    type            (node                             ), pointer               :: doc                              , element
    double precision                                                           :: matterRadiationEqualityRedshift
    integer                                                                    :: ioStatus
    !![
    <constructorAssign variables="*criticalOverdensityCDM, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_, *darkMatterParticle_, *linearGrowth_, useFittingFunction"/>
    !!]

    ! Compute corresponding Jeans mass.
    matterRadiationEqualityRedshift=+3600.0d0                                                          &
         &                          *(                                                                 &
         &                            +self%cosmologyParameters_%OmegaMatter   (                  )    &
         &                            *self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                            /0.15d0                                                          &
         &                          )                                                                  &
         &                          -1.0d0
    select type (darkMatterParticleWDMThermal_ => self%darkMatterParticle_)
    class is (darkMatterParticleWDMThermal)
       self%jeansMass=+3.06d8                                                                    &
            &         *((1.0d0+matterRadiationEqualityRedshift)/3000.0d0)                **1.5d0 &
            &         *sqrt(                                                                     &
            &               +self%cosmologyParameters_%OmegaMatter   (                  )        &
            &               *self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2     &
            &               /0.15d0                                                              &
            &              )                                                                     &
            &         /(darkMatterParticleWDMThermal_%degreesOfFreedomEffective()/1.5d0)         &
            &         /(darkMatterParticleWDMThermal_%mass                     ()/1.0d0) **4
    class default
       call Error_Report('critical overdensity expects a thermal warm dark matter particle'//{introspection:location})
    end select
    ! Read in the tabulated critical overdensity scaling.
    !$omp critical (FoX_DOM_Access)
    doc => XML_Parse(char(inputPath(pathTypeDataStatic))//"darkMatter/criticalOverdensityWarmDarkMatterBarkana.xml",iostat=ioStatus)
    if (ioStatus /= 0) call Error_Report('unable to find or parse the tabulated data'//{introspection:location})
    ! Extract the datum lists.
    element    => XML_Get_First_Element_By_Tag_Name(doc,"mass" )
    call XML_Array_Read(element,"datum",self%deltaTableMass )
    element    => XML_Get_First_Element_By_Tag_Name(doc,"delta")
    call XML_Array_Read(element,"datum",self%deltaTableDelta)
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    self%deltaTableCount=size(self%deltaTableMass)
    ! Convert tabulations to logarithmic versions.
    self%deltaTableMass =log(self%deltaTableMass )
    self%deltaTableDelta=log(self%deltaTableDelta)
    ! Build interpolator.
    self%interpolator_=interpolator(self%deltaTableMass,self%deltaTableDelta,interpolationType=GSL_Interp_CSpline,extrapolationType=extrapolationTypeFix)
    return
  end function barkana2001WDMConstructorInternal

  subroutine barkana2001WDMDestructor(self)
    !!{
    Destructor for the barkana2001WDM critical overdensity for collapse class.
    !!}
    implicit none
    type(criticalOverdensityBarkana2001WDM), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%criticalOverdensityCDM"   />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%darkMatterParticle_"      />
    <objectDestructor name="self%linearGrowth_"            />
    !!]
    return
  end subroutine barkana2001WDMDestructor

  double precision function barkana2001WDMValue(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Returns a mass scaling for critical overdensities based on the results of \cite{barkana_constraints_2001}. This method
    assumes that their results for the original collapse barrier (i.e. the critical overdensity, and which they call $B_0$)
    scale with the effective Jeans mass of the warm dark matter particle as computed using their eqn.~(10).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (criticalOverdensityBarkana2001WDM), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: time                         , expansionFactor
    logical                                            , intent(in   ), optional :: collapsing
    double precision                                   , intent(in   ), optional :: mass
    type            (treeNode                         ), intent(inout), optional :: node
    double precision                                   , parameter               :: massScaleFreeMinimum=- 10.0d0
    double precision                                   , parameter               :: massScaleFreeLarge  =+100.0d0
    double precision                                                             :: exponentialFit               , massScaleFree   , &
         &                                                                          powerLawFit                  , smoothTransition
    !$GLC attributes unused :: node

    ! Validate.
    if (.not.present(mass)) call Error_Report('mass is required for this critical overdensity class'//{introspection:location})
    ! Determine the scale-free mass.
    massScaleFree=log(mass/self%jeansMass)
    ! Compute the mass scaling via a fitting function or interpolation in tabulated results.
    if (self%useFittingFunction) then
       ! Impose a minimum to avoid divergence in fit.
       massScaleFree   =max(massScaleFree,massScaleFreeMinimum)
       ! Check for very large masses, and simply set critical overdensity multiplier to unity in such cases.
       if (massScaleFree > massScaleFreeLarge*abs(log(fitParameterE))/fitParameterF) then
          barkana2001WDMValue=1.0d0
       else
          smoothTransition=+1.0d0                    &
               &           /(                        &
               &             +1.0d0                  &
               &             +exp(                   &
               &                  +(                 &
               &                    +massScaleFree   &
               &                    +fitParameterA   &
               &                   )                 &
               &                  /  fitParameterB   &
               &                 )                   &
               &            )
          powerLawFit     =+     fitParameterC       &
               &           /exp(                     &
               &                +fitParameterD       &
               &                *massScaleFree       &
               &               )
          if (smoothTransition < 1.0d0) then
             exponentialFit=+exp(                    &
                  &              +     fitParameterE &
                  &              /exp(               &
                  &                   +fitParameterF &
                  &                   *massScaleFree &
                  &                  )               &
                  &             )
          else
             exponentialFit=+0.0d0
          end if
          barkana2001WDMValue=+       smoothTransition *powerLawFit    &
               &              +(1.0d0-smoothTransition)*exponentialFit
       end if
    else
       if      (massScaleFree > self%deltaTableMass(self%deltaTableCount)) then
          barkana2001WDMValue=1.0d0
       else if (massScaleFree < self%deltaTableMass(                   1)) then
          barkana2001WDMValue=exp(                           &
               &                  +  self%deltaTableDelta(1) &
               &                  +(                         &
               &                    +massScaleFree           &
               &                    -self%deltaTableMass (1) &
               &                   )                         &
               &                  *smallMassLogarithmicSlope &
               &                 )
       else
          barkana2001WDMValue=exp(self%interpolator_%interpolate(massScaleFree))
       end if
    end if
    barkana2001WDMValue=barkana2001WDMValue*self%criticalOverdensityCDM%value(time,expansionFactor,collapsing,mass)
    return
  end function barkana2001WDMValue

  double precision function barkana2001WDMGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to time of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityBarkana2001WDM), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: time      , expansionFactor
    logical                                            , intent(in   ), optional :: collapsing
    double precision                                   , intent(in   ), optional :: mass
    type            (treeNode                         ), intent(inout), optional :: node
    !$GLC attributes unused :: node

    barkana2001WDMGradientTime=+self                       %value       (time,expansionFactor,collapsing,mass) &
         &                     *self%criticalOverdensityCDM%gradientTime(time,expansionFactor,collapsing,mass) &
         &                     /self%criticalOverdensityCDM%value       (time,expansionFactor,collapsing,mass)
    return
  end function barkana2001WDMGradientTime

  double precision function barkana2001WDMGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to mass of critical overdensity at the given time and mass.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (criticalOverdensityBarkana2001WDM), intent(inout)           :: self
    double precision                                   , intent(in   ), optional :: time                    , expansionFactor
    logical                                            , intent(in   ), optional :: collapsing
    double precision                                   , intent(in   ), optional :: mass
    type            (treeNode                         ), intent(inout), optional :: node
    double precision                                                             :: exponentialFit          , exponentialFitGradient, &
         &                                                                          massScaleFree           , powerLawFit           , &
         &                                                                          powerLawFitGradient     , smoothTransition      , &
         &                                                                          smoothTransitionGradient
    !$GLC attributes unused :: node

    ! Validate.
    if (.not.present(mass)) call Error_Report('mass is required for this critical overdensity class'//{introspection:location})
    ! Determine the scale-free mass.
    massScaleFree=log(mass/self%jeansMass)
    ! Compute the mass scaling via a fitting function or interpolation in tabulated results.
    if (self%useFittingFunction) then
       smoothTransition          =+1.0d0                     &
            &                     /(                         &
            &                       +1.0d0                   &
            &                       +exp(                    &
            &                            +(                  &
            &                              +massScaleFree    &
            &                              +fitParameterA    &
            &                             )                  &
            &                            /  fitParameterB    &
            &                           )                    &
            &                      )
       powerLawFit               =+     fitParameterC        &
            &                     /exp(                      &
            &                          +fitParameterD        &
            &                          *massScaleFree        &
            &                         )
       exponentialFit            =+exp(                      &
            &                          +     fitParameterE   &
            &                          /exp(                 &
            &                               +fitParameterF   &
            &                               *massScaleFree   &
            &                              )                 &
            &                         )
       powerLawFitGradient       =-fitParameterD             &
            &                     *powerLawFit               &
            &                     /exp(                      &
            &                          +massScaleFree        &
            &                         )
       exponentialFitGradient    =-exponentialFit            &
            &                     *fitParameterF             &
            &                     *fitParameterE             &
            &                     /exp(                      &
            &                          +(                    &
            &                            +1.0d0              &
            &                            +fitParameterF      &
            &                           )                    &
            &                          *massScaleFree        &
            &                         )
       smoothTransitionGradient  =-exp(                      &
            &                          +(                    &
            &                            +     massScaleFree &
            &                            +     fitParameterA &
            &                           )                    &
            &                          /       fitParameterB &
            &                         )                      &
            &                        /         fitParameterB &
            &                        /(                      &
            &                          +1.0d0                &
            &                          +exp(                 &
            &                               +(               &
            &                                 +massScaleFree &
            &                                 +fitParameterA &
            &                                )               &
            &                               /  fitParameterB &
            &                              )                 &
            &                         )**2                   &
            &                        /exp(                   &
            &                             +massScaleFree     &
            &                            )
       barkana2001WDMGradientMass=+(                                                                                         &
            &                       +       smoothTransition *powerLawFitGradient   +smoothTransitionGradient*powerLawFit    &
            &                       +(1.0d0-smoothTransition)*exponentialFitGradient-smoothTransitionGradient*exponentialFit &
            &                      )                                                                                         &
            &                     /self%jeansMass
    else
       if      (massScaleFree > self%deltaTableMass(self%deltaTableCount)) then
          barkana2001WDMGradientMass=+0.0d0
       else if (massScaleFree < self%deltaTableMass(                   1)) then
          barkana2001WDMGradientMass=+smallMassLogarithmicSlope                                                     &
               &                     *self%value(time,expansionFactor,collapsing,mass)                              &
               &                     /mass
       else
          barkana2001WDMGradientMass=+self%interpolator_%derivative(                                massScaleFree) &
               &                     *self               %value    (time,expansionFactor,collapsing,mass         ) &
               &                     /                                                              mass
       end if
    end if
    ! Include gradient from CDM critical overdensity.
    barkana2001WDMGradientMass=+barkana2001WDMGradientMass                                                     &
         &                     *self%criticalOverdensityCDM%value       (time,expansionFactor,collapsing,mass) &
         &                     +self                       %value       (time,expansionFactor,collapsing,mass) &
         &                     /self%criticalOverdensityCDM%value       (time,expansionFactor,collapsing,mass) &
         &                     *self%criticalOverdensityCDM%gradientMass(time,expansionFactor,collapsing,mass)
    return
  end function barkana2001WDMGradientMass

  logical function barkana2001WDMIsMassDependent(self)
    !!{
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensityBarkana2001WDM), intent(inout) :: self
    !$GLC attributes unused :: self

    barkana2001WDMIsMassDependent=.true.
    return
  end function barkana2001WDMIsMassDependent

  logical function barkana2001WDMIsNodeDependent(self)
    !!{
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensityBarkana2001WDM), intent(inout) :: self
    !$GLC attributes unused :: self

    barkana2001WDMIsNodeDependent=.false.
    return
  end function barkana2001WDMIsNodeDependent

  logical function barkana2001WDMIsTreeDependent(self)
    !!{
    Return whether the critical overdensity is tree dependent.
    !!}
    implicit none
    class(criticalOverdensityBarkana2001WDM), intent(inout) :: self
    !$GLC attributes unused :: self

    barkana2001WDMIsTreeDependent=.false.
    return
  end function barkana2001WDMIsTreeDependent

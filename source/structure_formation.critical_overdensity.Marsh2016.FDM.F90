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
  Implements a critical overdensity for collapse using the \gls{fdm} modifier of \cite{marsh_warmAndFuzzy_2016}.
  !!}
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Dark_Matter_Particles  , only : darkMatterParticleClass
  use :: Numerical_Interpolation, only : interpolator

  !![
  <criticalOverdensity name="criticalOverdensityMarsh2016FDM">
   <description>
    A critical overdensity for collapse class based on the \gls{fdm} modifier of \cite{marsh_warmAndFuzzy_2016} applied to
    some other, \gls{cdm} critical overdensity class. Specifically, the critical overdensity is multiplied by a factor
    \begin{equation}
     h_\mathrm{F}(x) \exp [a_3 x^{-a_4}]+[1-h_\mathrm{F}(x)]\exp [a_5 x^{-a_6}],
    \end{equation}
    where $x=M/M_J$ with $M$ the mass in question, $M_\mathrm{J}$ the effective Jeans mass of the fuzzy dark matter as defined by
    \citeauthor{marsh_warmAndFuzzy_2016}~[\citeyear{marsh_warmAndFuzzy_2016}; their eqn.~18]:
    \begin{equation}
     M_\mathrm{J} = a_1\times 10^8 \left(\frac{m_a}{10^{-22}\text{ eV}}\right)^{-3/2}\left(\frac{\Omega_m h^2}{0.14} \right)^{1/4} h^{-1}M_{\odot},
    \end{equation}
    the function $h_\mathrm{F}$ is given by
    \begin{equation}
     h_\mathrm{F}(x) =(1/2)\{1-\tanh [M_J(x-a_2)]\}.
    \end{equation}
    The best-fit parameters $a_i$ are
    $\{a_1,a_2,a_3,a_4,a_5,a_6\}=\{3.4,1.0,1.8,0.5,1.7,0.9\}$.
   </description>
  </criticalOverdensity>
  !!]
  type, extends(criticalOverdensityClass) :: criticalOverdensityMarsh2016FDM
     !!{
     A critical overdensity for collapse class which modifies another transfer function using the \gls{fdm} modifier of \cite{marsh_warmAndFuzzy_2016}.
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
     final     ::                    marsh2016FDMDestructor
     procedure :: value           => marsh2016FDMValue
     procedure :: gradientTime    => marsh2016FDMGradientTime
     procedure :: gradientMass    => marsh2016FDMGradientMass
     procedure :: isMassDependent => marsh2016FDMIsMassDependent
     procedure :: isNodeDependent => marsh2016FDMIsNodeDependent
  end type criticalOverdensityMarsh2016FDM

  interface criticalOverdensityMarsh2016FDM
     !!{
     Constructors for the {\normalfont \ttfamily marsh2016FDM} critical overdensity for collapse class.
     !!}
     module procedure marsh2016FDMConstructorParameters
     module procedure marsh2016FDMConstructorInternal
  end interface criticalOverdensityMarsh2016FDM

  ! Parameters of fitting functions.
  double precision, parameter :: fitParameterA1=+3.4d0
  double precision, parameter :: fitParameterA2=+1.0d0
  double precision, parameter :: fitParameterA3=+1.8d0
  double precision, parameter :: fitParameterA4=+0.5d0
  double precision, parameter :: fitParameterA5=+1.7d0
  double precision, parameter :: fitParameterA6=+0.9d0

contains

  function marsh2016FDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily marsh2016FDM} critical overdensity for collapse class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (criticalOverdensityMarsh2016FDM)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (criticalOverdensityClass       ), pointer       :: criticalOverdensityCDM
    class           (cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterParticleClass        ), pointer       :: darkMatterParticle_
    class           (linearGrowthClass              ), pointer       :: linearGrowth_
    logical                                                          :: useFittingFunction

    !![
    <inputParameter>
      <name>useFittingFunction</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether the fuzzy dark matter critical overdensity mass scaling should be computed from a fitting function or from tabulated data.</description>
    </inputParameter>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensityCDM"    source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"       name="darkMatterParticle_"       source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !!]
    ! Call the internal constructor
    self=marsh2016FDMConstructorInternal(criticalOverdensityCDM,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,linearGrowth_,useFittingFunction)
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
  end function marsh2016FDMConstructorParameters

  function marsh2016FDMConstructorInternal(criticalOverdensityCDM,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,darkMatterParticle_,linearGrowth_,useFittingFunction) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily marsh2016FDM} critical overdensity for collapse class.
    !!}
    use :: Cosmology_Parameters        , only : hubbleUnitsLittleH
    use :: Dark_Matter_Particles       , only : darkMatterParticleFuzzyDarkMatter
    use :: FoX_DOM                     , only : destroy                          , node                             , parseFile
    use :: Error                       , only : Error_Report
    use :: Input_Paths                 , only : inputPath                        , pathTypeDataStatic
    use :: IO_XML                      , only : XML_Array_Read                   , XML_Get_First_Element_By_Tag_Name
    use :: Numerical_Interpolation     , only : GSL_Interp_CSpline
    use :: Table_Labels                , only : extrapolationTypeFix
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    type            (criticalOverdensityMarsh2016FDM)                        :: self
    class           (criticalOverdensityClass       ), target, intent(in   ) :: criticalOverdensityCDM
    class           (cosmologyParametersClass       ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass        ), target, intent(in   ) :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass  ), target, intent(in   ) :: cosmologicalMassVariance_
    class           (darkMatterParticleClass        ), target, intent(in   ) :: darkMatterParticle_
    class           (linearGrowthClass              ), target, intent(in   ) :: linearGrowth_
    double precision                                                         :: m22
    logical                                                  , intent(in   ) :: useFittingFunction
    type            (node                           ), pointer               :: doc                      , element
    integer                                                                  :: ioStatus
    !![
    <constructorAssign variables="*criticalOverdensityCDM, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_, *darkMatterParticle_, *linearGrowth_, useFittingFunction"/>
    !!]

    ! Compute corresponding Jeans mass.
    select type (darkMatterParticleFDM_ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       m22           =self%darkMatterParticle_%mass()*kilo/1.0d-22
       self%jeansMass=+fitParameterA1                                        &
            &         *1.0d8                                                             &
            &         *m22**(-3.0d0/2.0d0)                                               &
            &         *(                                                                 &
            &           +self%cosmologyParameters_%OmegaMatter   (                  )    &
            &           *self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
            &           /0.14d0                                                          &
            &          )**(1.0d0/4.0d0)                                                  &
            &         /  self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    class default
       call Error_Report('critical overdensity expects a fuzzy dark matter particle'//{introspection:location})
    end select
    if (.not.self%useFittingFunction) then
       ! Read in the tabulated critical overdensity scaling.
       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(inputPath(pathTypeDataStatic))//"darkMatter/criticalOverdensityFuzzyDarkMatterMarsh.xml",iostat=ioStatus)
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
    end if
    return
  end function marsh2016FDMConstructorInternal

  subroutine marsh2016FDMDestructor(self)
    !!{
    Destructor for the marsh2016FDM critical overdensity for collapse class.
    !!}
    implicit none
    type(criticalOverdensityMarsh2016FDM), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%criticalOverdensityCDM"   />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%darkMatterParticle_"      />
    <objectDestructor name="self%linearGrowth_"            />
    !!]
    return
  end subroutine marsh2016FDMDestructor

  double precision function marsh2016FDMValue(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Returns a mass scaling for critical overdensities based on the results of \cite{marsh_warmAndFuzzy_2016}. This method
    assumes that their results for the critical overdensity scale with the Jeans mass of the fuzzy dark matter particle
    as computed using their eqn.~(18).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (criticalOverdensityMarsh2016FDM), intent(inout)           :: self
    double precision                                 , intent(in   ), optional :: time                       , expansionFactor
    logical                                          , intent(in   ), optional :: collapsing
    double precision                                 , intent(in   ), optional :: mass
    type            (treeNode                       ), intent(inout), optional :: node
    double precision                                 , parameter               :: massScaleFreeMinimum=2.0d-3
    double precision                                 , parameter               :: massScaleFreeLarge  =1.0d50
    double precision                                                           :: hF                         , massScaleFree  , &
         &                                                                        exponentialFit1            , exponentialFit2
    !$GLC attributes unused :: node

    ! Validate.
    if (.not.present(mass)) call Error_Report('mass is required for this critical overdensity class'//{introspection:location})
    ! Determine the scale-free mass.
    massScaleFree=mass/self%jeansMass
    ! Compute the mass scaling via a fitting function or interpolation in tabulated results.
    if (self%useFittingFunction) then
       ! Impose a minimum to avoid divergence in fit.
       massScaleFree=max(massScaleFree,massScaleFreeMinimum)
       ! Check for very large masses, and simply set critical overdensity multiplier to unity in such cases.
       if (massScaleFree > massScaleFreeLarge) then
          marsh2016FDMValue=1.0d0
       else
          hF               =+0.5d0                                &
               &            *(                                    &
               &              +1.0d0                              &
               &              -tanh(                              &
               &                    +self%jeansMass               &
               &                    *(                            &
               &                      +massScaleFree              &
               &                      -fitParameterA2 &
               &                     )                            &
               &                   )                              &
               &             )
          exponentialFit1  =exp(                                              &
               &                +fitParameterA3                   &
               &                *massScaleFree**(-fitParameterA4) &
               &               )
          exponentialFit2  =exp(                                              &
               &                +fitParameterA5                   &
               &                *massScaleFree**(-fitParameterA6) &
               &               )
          marsh2016FDMValue=        hF *exponentialFit1 &
               &            +(1.0d0-hF)*exponentialFit2
       end if
    else
       if      (massScaleFree > self%deltaTableMass(self%deltaTableCount)) then
          marsh2016FDMValue=1.0d0
       else if (massScaleFree < self%deltaTableMass(                   1)) then
          marsh2016FDMValue=self%deltaTableMass(1)
       else
          marsh2016FDMValue=self%interpolator_%interpolate(massScaleFree)
       end if
    end if
    marsh2016FDMValue=marsh2016FDMValue*self%criticalOverdensityCDM%value(time,expansionFactor,collapsing,mass)
    return
  end function marsh2016FDMValue

  double precision function marsh2016FDMGradientTime(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to time of critical overdensity at the given time and mass.
    !!}
    implicit none
    class           (criticalOverdensityMarsh2016FDM), intent(inout)           :: self
    double precision                                 , intent(in   ), optional :: time      , expansionFactor
    logical                                          , intent(in   ), optional :: collapsing
    double precision                                 , intent(in   ), optional :: mass
    type            (treeNode                       ), intent(inout), optional :: node
    !$GLC attributes unused :: node

    marsh2016FDMGradientTime=+self                       %value       (time,expansionFactor,collapsing,mass) &
         &                   *self%criticalOverdensityCDM%gradientTime(time,expansionFactor,collapsing,mass) &
         &                   /self%criticalOverdensityCDM%value       (time,expansionFactor,collapsing,mass)
    return
  end function marsh2016FDMGradientTime

  double precision function marsh2016FDMGradientMass(self,time,expansionFactor,collapsing,mass,node)
    !!{
    Return the gradient with respect to mass of critical overdensity at the given time and mass.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (criticalOverdensityMarsh2016FDM), intent(inout)           :: self
    double precision                                 , intent(in   ), optional :: time                       , expansionFactor
    logical                                          , intent(in   ), optional :: collapsing
    double precision                                 , intent(in   ), optional :: mass
    type            (treeNode                       ), intent(inout), optional :: node
    double precision                                 , parameter               :: massScaleFreeMinimum=2.0d-3
    double precision                                 , parameter               :: massScaleFreeLarge  =1.0d50
    double precision                                                           :: hF                         , dhF_dM                 , &
         &                                                                        exponentialFit1            , exponentialFit2        , &
         &                                                                        exponentialFitGradient1    , exponentialFitGradient2, &
         &                                                                        massScaleFree
    !$GLC attributes unused :: node

    ! Validate.
    if (.not.present(mass)) call Error_Report('mass is required for this critical overdensity class'//{introspection:location})
    ! Determine the scale-free mass.
    massScaleFree=mass/self%jeansMass
    ! Compute the mass scaling via a fitting function or interpolation in tabulated results.
    if (self%useFittingFunction) then
       ! Impose a minimum to avoid divergence in fit.
       massScaleFree=max(massScaleFree,massScaleFreeMinimum)
       ! Check for very large masses, and simply set critical overdensity multiplier to unity in such cases.
       if (massScaleFree > massScaleFreeLarge) then
          marsh2016FDMGradientMass=0.0d0
       else
          hF                      =+0.5d0                                  &
               &                   *(                                      &
               &                     +1.0d0                                &
               &                     -tanh(                                &
               &                           +self%jeansMass                 &
               &                           *(                              &
               &                             +massScaleFree                &
               &                             -fitParameterA2               &
               &                            )                              &
               &                          )                                &
               &                    )
          dhF_dM                  =-0.5d0                                  &
               &                   *self%jeansMass                         &
               &                   /(                                      &
               &                     cosh(                                 &
               &                          +self%jeansMass                  &
               &                          *(                               &
               &                            +massScaleFree                 &
               &                            -fitParameterA2                &
               &                           )                               &
               &                         )                                 &
               &                    )**2.0d0
          exponentialFit1         =exp(                                    &
               &                       +fitParameterA3                     &
               &                       *massScaleFree**(-fitParameterA4)   &
               &                      )
          exponentialFit2         =exp(                                    &
               &                       +fitParameterA5                     &
               &                       *massScaleFree**(-fitParameterA6)   &
               &                      )
          exponentialFitGradient1 =-fitParameterA3                         &
               &                   *fitParameterA4                         &
               &                   *massScaleFree**(-fitParameterA4-1.0d0) &
               &                   *exponentialFit1
          exponentialFitGradient2 =-fitParameterA5                         &
               &                   *fitParameterA6                         &
               &                   *massScaleFree**(-fitParameterA6-1.0d0) &
               &                   *exponentialFit2
          marsh2016FDMGradientMass=+(                                      &
               &                     +dhF_dM                               &
               &                     *(exponentialFit1-exponentialFit2)    &
               &                     +       hF                            &
               &                     *exponentialFitGradient1              &
               &                     +(1.0d0-hF)                           &
               &                     *exponentialFitGradient2              &
               &                    )                                      &
               &                   /self%jeansMass
       end if
    else
       if      (massScaleFree > self%deltaTableMass(self%deltaTableCount)) then
          marsh2016FDMGradientMass=+0.0d0
       else if (massScaleFree < self%deltaTableMass(                   1)) then
          marsh2016FDMGradientMass=+self%interpolator_%derivative(self%deltaTableMass(1)) &
               &                   /self%jeansMass
       else
          marsh2016FDMGradientMass=+self%interpolator_%derivative(massScaleFree         ) &
               &                   /self%jeansMass
       end if
    end if
    ! Include gradient from CDM critical overdensity.
    marsh2016FDMGradientMass=+marsh2016FDMGradientMass                                                       &
         &                   *self%criticalOverdensityCDM%value       (time,expansionFactor,collapsing,mass) &
         &                   +self                       %value       (time,expansionFactor,collapsing,mass) &
         &                   /self%criticalOverdensityCDM%value       (time,expansionFactor,collapsing,mass) &
         &                   *self%criticalOverdensityCDM%gradientMass(time,expansionFactor,collapsing,mass)
    return
  end function marsh2016FDMGradientMass

  logical function marsh2016FDMIsMassDependent(self)
    !!{
    Return whether the critical overdensity is mass dependent.
    !!}
    implicit none
    class(criticalOverdensityMarsh2016FDM), intent(inout) :: self
    !$GLC attributes unused :: self

    marsh2016FDMIsMassDependent=.true.
    return
  end function marsh2016FDMIsMassDependent

  logical function marsh2016FDMIsNodeDependent(self)
    !!{
    Return whether the critical overdensity is node dependent.
    !!}
    implicit none
    class(criticalOverdensityMarsh2016FDM), intent(inout) :: self
    !$GLC attributes unused :: self

    marsh2016FDMIsNodeDependent=.false.
    return
  end function marsh2016FDMIsNodeDependent

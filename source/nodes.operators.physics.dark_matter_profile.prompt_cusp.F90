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
  Implements a node operator class that evaluates the properties of prompt cusps following the model of
  \cite{delos_cusp-halo_2025}.
  !!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Power_Spectra           , only : powerSpectrumClass
  use :: Linear_Growth           , only : linearGrowthClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass

  !![
  <nodeOperator name="nodeOperatorDarkMatterProfilePromptCusps">
   <description>
    A node operator class that evaluates the properties of prompt cusps following the model of \cite{delos_cusp-halo_2025}, with a
    log-normal scatter of $\mu \exp(-1/\sigma_0)$~dex added to the cusp amplitude, where $\mu=${\normalfont \ttfamily
    [coefficientScatter]}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfilePromptCusps
     !!{     
     A node operator class that evaluates the properties of prompt cusps following the model of \cite{delos_cusp-halo_2025}.
     !!}
     private
     class           (powerSpectrumClass        ), pointer                   :: powerSpectrum_              => null()
     class           (linearGrowthClass         ), pointer                   :: linearGrowth_               => null()
     class           (cosmologyParametersClass  ), pointer                   :: cosmologyParameters_        => null()
     class           (cosmologyFunctionsClass   ), pointer                   :: cosmologyFunctions_         => null()
     class           (darkMatterHaloScaleClass  ), pointer                   :: darkMatterHaloScale_        => null()
     class           (darkMatterProfileDMOClass ), pointer                   :: darkMatterProfileDMO_       => null()
     class           (virialDensityContrastClass), pointer                   :: virialDensityContrast_      => null()
     double precision                            , allocatable, dimension(:) :: sigma_
     logical                                                                 :: growthIsWavenumberDependent          , nonConvergenceIsFatal
     double precision                                                        :: alpha                                , beta                 , &
          &                                                                     kappa                                , C                    , &
          &                                                                     p                                    , coefficientScatter
     integer                                                                 :: promptCuspMassID                     , promptCuspAmplitudeID, &
          &                                                                     promptCuspNFWYID                     , promptCuspNFWScaleID , &
          &                                                                     promptCuspNFWGrowthRateID
   contains
     !![
     <methods>
       <method method="sigma" description="Evaluate $\sigma_j^2 = \int_0^\infty \frac{\mathrm{d}k}{k} \mathcal{P}(k,t) k^{2j}$ where $\mathcal{P}(k) = k^3 P(k) / 2 \pi^2$ is the dimensionless form of the power spectrum."/>
     </methods>
     !!]
     final     ::                                        darkMatterProfilePromptCuspsDestructor
     procedure :: nodeTreeInitialize                  => darkMatterProfilePromptCuspsNodeTreeInitialize
     procedure :: nodeInitialize                      => darkMatterProfilePromptCuspsNodeInitialize
     procedure :: nodePromote                         => darkMatterProfilePromptCuspsNodePromote
     procedure :: differentialEvolutionSolveAnalytics => darkMatterProfilePromptCuspsSolveAnalytics
     procedure :: sigma                               => darkMatterProfilePromptCuspsNodeSigma
  end type nodeOperatorDarkMatterProfilePromptCusps
  
  interface nodeOperatorDarkMatterProfilePromptCusps
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfilePromptCusps} node operator class.
     !!}
     module procedure darkMatterProfilePromptCuspsConstructorParameters
     module procedure darkMatterProfilePromptCuspsConstructorInternal
  end interface nodeOperatorDarkMatterProfilePromptCusps

  ! Submodule-scope variables used in root-finding.
  class           (nodeOperatorDarkMatterProfilePromptCusps), pointer   :: self_
  double precision                                                      :: sigma0Collapse          , time_                    , &
       &                                                                   expansionFactor_        , concentrationFactorTarget, &
       &                                                                   massHalo_               , amplitudeCusp_           , &
       &                                                                   concentration_          , radiusScale_
  integer                                                               :: j_ 
  !$omp threadprivate(self_,sigma0Collapse,time_,expansionFactor_,concentrationFactorTarget,massHalo_,concentration_,radiusScale_,amplitudeCusp_,j_)

  ! Maximum allowed value of the y-parameter in the cusp-NFW profile. Values of 1 or greater are not valid. We limit here to a
  ! value close to 1.
  double precision                                          , parameter :: yMaximum        =0.999d0
  
contains
  
  function darkMatterProfilePromptCuspsConstructorParameters(parameters) result(self)
    !!{    
    Constructor for the \refClass{nodeOperatorDarkMatterProfilePromptCusps} node operator class which takes a parameter set as
    input.    
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfilePromptCusps)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (linearGrowthClass                       ), pointer       :: linearGrowth_
    class           (powerSpectrumClass                      ), pointer       :: powerSpectrum_
    class           (cosmologyParametersClass                ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass               ), pointer       :: darkMatterProfileDMO_
    class           (virialDensityContrastClass              ), pointer       :: virialDensityContrast_
    logical                                                                   :: nonConvergenceIsFatal
    double precision                                                          :: alpha                 , beta              , &
         &                                                                       kappa                 , C                 , &
         &                                                                       p                     , coefficientScatter

    !![
    <inputParameter>
      <name>nonConvergenceIsFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, failure to converge on a solution for the scale radius, $r_\mathrm{s}$, will result in a fatal error. Otherwise, only warnings are issued.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <defaultValue>24.0d0</defaultValue>
      <defaultSource>\citep[][Table 3]{delos_cusp-halo_2025}</defaultSource>
      <description>The coefficient, $\alpha$ of the cusp amplitude, $A$, in the peak-cusp connection of the \cite{delos_cusp-halo_2025} prompt cusp model.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>7.3d0</defaultValue>
      <defaultSource>\citep[][Table 3]{delos_cusp-halo_2025}</defaultSource>
      <description>The coefficient, $\beta$, of the cusp mass, $m$, in the peak-cusp connection of the \cite{delos_cusp-halo_2025} prompt cusp model.</description>
    </inputParameter>
    <inputParameter>
      <name>C</name>
      <source>parameters</source>
      <defaultValue>0.8d0</defaultValue>
      <defaultSource>\citep[][Table 3]{delos_cusp-halo_2025}</defaultSource>
      <description>The coefficient, $C$, of the cusp $A$--$m$ relation in the \cite{delos_cusp-halo_2025} prompt cusp model.</description>
    </inputParameter>
    <inputParameter>
      <name>p</name>
      <source>parameters</source>
      <defaultValue>1.9d0</defaultValue>
      <defaultSource>\citep[][Table 3]{delos_cusp-halo_2025}</defaultSource>
      <description>The exponent, $p$, of the cusp $A$--$m$ relation in the \cite{delos_cusp-halo_2025} prompt cusp model.</description>
    </inputParameter>
    <inputParameter>
      <name>kappa</name>
      <source>parameters</source>
      <defaultValue>4.5d0</defaultValue>
      <defaultSource>\citep[][Table 3]{delos_cusp-halo_2025}</defaultSource>
      <description>The parameter, $\kappa$, of the mass growth factor in the \cite{delos_cusp-halo_2025} prompt cusp model.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientScatter</name>
      <source>parameters</source>
      <defaultValue>0.195d0</defaultValue>
      <defaultSource>Delos (private communication)</defaultSource>
      <description>The parameter, $\mu$, in the expression for the scatter in cusp amplitude.</description>
    </inputParameter>
    <objectBuilder class="linearGrowth"          name="linearGrowth_"          source="parameters"/>
    <objectBuilder class="powerSpectrum"         name="powerSpectrum_"         source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfilePromptCusps(nonConvergenceIsFatal,alpha,beta,kappa,C,p,coefficientScatter,linearGrowth_,powerSpectrum_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="linearGrowth_"         />
    <objectDestructor name="powerSpectrum_"        />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function darkMatterProfilePromptCuspsConstructorParameters

  function darkMatterProfilePromptCuspsConstructorInternal(nonConvergenceIsFatal,alpha,beta,kappa,C,p,coefficientScatter,linearGrowth_,powerSpectrum_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfilePromptCusps} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfilePromptCusps)                        :: self
    class           (linearGrowthClass                       ), intent(in   ), target :: linearGrowth_
    class           (powerSpectrumClass                      ), intent(in   ), target :: powerSpectrum_
    class           (cosmologyParametersClass                ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass               ), intent(in   ), target :: darkMatterProfileDMO_
    class           (virialDensityContrastClass              ), intent(in   ), target :: virialDensityContrast_
    logical                                                   , intent(in   )         :: nonConvergenceIsFatal
    double precision                                          , intent(in   )         :: alpha                 , beta              , &
         &                                                                               kappa                 , C                 , &
         &                                                                               p                     , coefficientScatter
    !![
    <constructorAssign variables="nonConvergenceIsFatal, alpha, beta, kappa, C, p, coefficientScatter, *linearGrowth_, *powerSpectrum_, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *virialDensityContrast_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="promptCuspAmplitude"       id="self%promptCuspAmplitudeID"     isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspMass"            id="self%promptCuspMassID"          isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWY"            id="self%promptCuspNFWYID"          isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWScale"        id="self%promptCuspNFWScaleID"      isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWGrowthRateID" id="self%promptCuspNFWGrowthRateID" isEvolvable="no" isCreator="yes"/>
    !!]
    self%growthIsWavenumberDependent=self%linearGrowth_%isWavenumberDependent()
    return
  end function darkMatterProfilePromptCuspsConstructorInternal

  subroutine darkMatterProfilePromptCuspsDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorDarkMatterProfilePromptCusps} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfilePromptCusps), intent(inout) :: self

    !![
    <objectDestructor name="self%linearGrowth_"         />
    <objectDestructor name="self%powerSpectrum_"        />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%darkMatterHaloScale_"  />
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine darkMatterProfilePromptCuspsDestructor

  subroutine darkMatterProfilePromptCuspsNodeTreeInitialize(self,node)
    !!{
    Initialize dark matter profile prompt cusp properties.
    !!}
    use :: Calculations_Resets                 , only : Calculations_Reset
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Error                               , only : Error_Report                       , Warn
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , nodeComponentDarkMatterProfile
    use :: Lambert_Ws                          , only : Lambert_W0
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Comparison                , only : Values_Agree
    use :: Root_Finder                         , only : rootFinder                         , rangeExpandMultiplicative     , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    type            (treeNode                                )               , pointer :: nodeChild
    class           (nodeComponentBasic                      )               , pointer :: basic
    class           (nodeComponentDarkMatterProfile          )               , pointer :: darkMatterProfile
    integer                                                   , parameter              :: iterationCountMaximum    =100
    double precision                                          , parameter              :: toleranceRelative        =1.0d-3
    type            (rootFinder                              )               , save    :: finderCollapse                   , finderRadius
    logical                                                                  , save    :: finderCollapseInitialized=.false., finderRadiusInitialized=.false.
    !$omp threadprivate(finderCollapse,finderRadius,finderCollapseInitialized,finderRadiusInitialized)
    double precision                                                         , save    :: errorFractionalMaximum=  0.0d+0
    logical                                                                            :: computeCusp
    integer                                                                            :: iterationCount
    double precision                                                                   :: sigma0                           , sigma2                         , &
         &                                                                                densityMean                      , densityMeanCollapse            , &
         &                                                                                sigma2Collapse                   , timeCollapse                   , &
         &                                                                                amplitude                        , mass                           , &
         &                                                                                y                                , radiusScale                    , &
         &                                                                                concentration                    , densityScale                   , &
         &                                                                                radiusScalePrevious              , mass200Critical                , &
         &                                                                                gamma                            , zeta                           , &
         &                                                                                radiusMinus2                     , densityContrast                , &
         &                                                                                scatterRandom                    , errorFractional
    
    ! Compute cusp properties for leaf nodes.
    computeCusp=.not.associated(node%firstChild)
    if (.not.computeCusp) then
       ! For non-leaf nodes, use the cusp properties from the leaf node.
       nodeChild => node
       do while (associated(nodeChild%firstChild))
          nodeChild => nodeChild%firstChild
       end do
       darkMatterProfile =>  nodeChild        %darkMatterProfile        (                          )
       amplitude         =   darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspAmplitudeID)
       mass              =   darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspMassID     )
       y                 =   darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspNFWYID     )
    else
       ! For leaf nodes, assume no cusp on the initial iteration.
       y                 =   0.0d0
    end if
    ! Initialize the root finder.
    if (.not.finderCollapseInitialized) then
       finderCollapse=rootFinder(                                                             &
            &                    rootFunction                 =timeCollapseRoot             , &
            &                    toleranceRelative            =1.0d-6                       , &
            &                    rangeExpandUpward            =2.0d+0                       , &
            &                    rangeExpandDownward          =0.5d+0                       , &
            &                    rangeExpandType              =rangeExpandMultiplicative    , &
            &                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                    rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
            &                   )
      finderCollapseInitialized=.true.
    end if
    ! Evaluate the required integrals over the power spectrum at the time of this node.
    self_             =>  self
    basic             =>  node                                  %basic               (                   )
    darkMatterProfile =>  node                                  %darkMatterProfile   (                   )
    sigma0            =   self                                  %sigma               (0,     basic%time())
    sigma2            =   self                                  %sigma               (2,     basic%time())
    densityMean       =   self             %cosmologyFunctions_ %matterDensityEpochal(  time=basic%time())
    ! We assume (following Delos 2025) that the "scale radius" that has been set is actually r₋₂. We must iteratively solve for
    ! the actual scale radius. We use rₛ=r₋₂ as our initial guess.
    radiusMinus2       =darkMatterProfile%scale()
    radiusScale        =radiusMinus2
    radiusScalePrevious=huge(0.0d0)
    iterationCount     =0
    scatterRandom      =0.0d0
    ! Begin iteration.
    do while (.not.Values_Agree(radiusScale,radiusScalePrevious,relTol=toleranceRelative) .and. iterationCount < iterationCountMaximum)
       iterationCount     =iterationCount+1
       radiusScalePrevious=radiusScale
       if (computeCusp) then
          ! Compute the mass following the 200 times critical density definition as used by Delos (2025; note that Delos actually
          ! states 200 times mean density, but these are Einstein-de Sitter cosmologies where mean and critical are equivalent -
          ! when comparing with other simulations Sten Delos uses 200 times critical [private communication]).
          call Calculations_Reset(node)
          densityContrast=+200.0d0                                                                    & ! 200 ρ_crit/ρ_mean since the density contrast is
               &          /self%cosmologyFunctions_%OmegaMatterEpochal(time=basic%timeLastIsolated())   ! relative to mean density.
          mass200Critical=Dark_Matter_Profile_Mass_Definition(                                                    &
               &                                              node                  =node                       , &
               &                                              densityContrast       =densityContrast            , &
               &                                              cosmologyParameters_  =self%cosmologyParameters_  , &
               &                                              cosmologyFunctions_   =self%cosmologyFunctions_   , &
               &                                              virialDensityContrast_=self%virialDensityContrast_, &
               &                                              darkMatterProfileDMO_ =self%darkMatterProfileDMO_ , &
               &                                              useLastIsolatedTime   =.true.                       &
               &                                             )
          ! Find the collapse time of this node.
          sigma0Collapse     =+            (2.0d0*self%p-1.0d0)*self%kappa/3.0d0              &
               &              /Lambert_W0(                                                    &
               &                          +(2.0d0*self%p-1.0d0)*self%kappa/3.0d0              &
               &                          *(self%beta*self%C**2/self%alpha**2)**(1.0d0/3.0d0) &
               &                          *(                                                  &
               &                            +exp(self%kappa/sigma0)                           &
               &                            *mass200Critical                                  &
               &                            /densityMean                                      &
               &                            /(sigma0/sigma2)**1.5d0                           &
               &                           )**((2.0d0*self%p-1.0d0)/3.0d0)                    &
               &                         )
          timeCollapse       =finderCollapse            %find(                  rootGuess=basic%time        ())
          sigma2Collapse     =self                      %sigma               (2,time     =      timeCollapse  )
          densityMeanCollapse=self  %cosmologyFunctions_%matterDensityEpochal(  time     =      timeCollapse  )
          ! Evaluate the prompt cusp parameters.
          amplitude=self%alpha*(self%alpha/self%C/self%beta**self%p)**(1.0d0/(2.0d0*self%p-1.0d0))*densityMeanCollapse*sigma0Collapse**((9.0d0-6.0d0*self%p)/(4.0d0-8.0d0*self%p))/sigma2Collapse**0.75d0
          mass     =self%beta *(self%alpha/self%C/self%beta**self%p)**(2.0d0/(2.0d0*self%p-1.0d0))*densityMeanCollapse*sigma0Collapse**((9.0d0-6.0d0*self%p)/(2.0d0-4.0d0*self%p))/sigma2Collapse**1.50d0
          ! Add scatter to the cusp amplitude.
          if (self%coefficientScatter > 0.0d0) then
             if (iterationCount == 1) scatterRandom=node%hostTree%randomNumberGenerator_%standardNormalSample()
             amplitude=+amplitude                         &
                  &    *10.0d0**(                         &
                  &              +self%coefficientScatter &
                  &              *exp(                    &
                  &                   -1.0d0              &
                  &                   /sigma0             &
                  &                  )                    &
                  &              *scatterRandom           &
                  &             )
          end if
       end if
       ! Adjust the initial concentration (and r₋₂ radius) if necessary. The equation in footnote 9 of Delos (2025;
       ! https://ui.adsabs.harvard.edu/abs/2025arXiv250603240D) gives a condition that must be satisfied by the concentration
       ! c=rᵥ/r₋₂. If that is not satisfied, adjust r₋₂ to satisfy that relation.
       concentration            =+self %darkMatterHaloScale_%radiusVirial(node)    &
            &                    /                           radiusMinus2
       concentrationFactorTarget=+self %darkMatterHaloScale_%densityMean (node)    &
            &                    *basic                     %mass        (    )    &
            &                    /                           amplitude         **2
       if (concentrationTargetRoot(concentration) < 0.0d0) then
          if (.not.finderRadiusInitialized) then
             finderRadius=rootFinder(                                                             &
                  &                  rootFunction                 =concentrationTargetRoot      , &
                  &                  toleranceRelative            =1.0d-6                       , &
                  &                  rangeExpandUpward            =2.0d+0                       , &
                  &                  rangeExpandDownward          =0.5d+0                       , &
                  &                  rangeExpandType              =rangeExpandMultiplicative    , &
                  &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
                  &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
                  &                 )
             finderRadiusInitialized=.true.
          end if
          ! Find a new concentration and use it to update r₋₂ (including some small tolerance to avoid being exactly at the
          ! limit).
          concentration= finderRadius                     %find         (concentration)
          radiusMinus2 =+self        %darkMatterHaloScale_%radiusVirial (node         ) &
               &        /                                  concentration                &
               &        *                                  yMaximum     **(2.0d0/3.0d0)
          if (iterationCount == 1) then
             ! If this adjustment in r₋₂ is made on the first iteration, reset our initial guess for the scale radius to ensure
             ! that we start from a stable point when seeking an iterative solution.
             radiusScale        =radiusMinus2
             radiusScalePrevious=radiusScale
          end if
       end if
       ! Compute the normalization of the cusp-NFW profile. Handle the case of y=0 here (which is assumed on the first iteration).
       concentration=+self%darkMatterHaloScale_%radiusVirial(node) &
            &        /                          radiusScale
       if (y > 0.0d0) then
          ! Cusp-NFW (y>0) case.
          densityScale=+basic%mass()                                                                                 &
               &       /radiusScale**3                                                                               &
               &       /4.0d0                                                                                        &
               &       /Pi                                                                                           &
               &       /(                                                                                            &
               &         + 2.0d0                       *asinh(sqrt(concentration            )/               y     ) &
               &         -(2.0d0-y**2)/sqrt(1.0d0-y**2)*atanh(sqrt(concentration*(1.0d0-y**2)/(concentration+y**2))) &
               &         -sqrt(concentration*(concentration+y**2))/(1.0d0+concentration)                             &
               &        )
       else
          ! NFW (y=0) case.
          densityScale=+basic%mass()                             &
               &       /4.0d0                                    &
               &       /Pi                                       &                                 
               &       /radiusScale**3                           &
               &       /(                                        &
               &         +              log(1.0d0+concentration) &
               &         -concentration/   (1.0d0+concentration) &
               &        )
       end if
       ! Compute a new scale radius to obtain the target r₋₂, using equations (23) from Delos (2025).
       gamma      =+amplitude           &
            &      /densityScale        &
            &      /radiusMinus2**1.5d0
       zeta       =+(                         &
            &        +        1.00d0          &
            &        +       18.00d0*gamma**2 &
            &        +        0.75d0*gamma    &
            &        *sqrt(                   &
            &              + 72.00d0          &
            &              +564.00d0*gamma**2 &
            &              +  6.00d0*gamma**4 &
            &             )                   &
            &       )**(1.0d0/3.0d0)
       radiusScale=+(                                     &
            &        +(+1.0d0/3.0d0-gamma**2/2.0d0)/zeta  &
            &        +(+1.0d0      +zeta          )/3.0d0 &
            &       )                                     &
            &      *radiusMinus2
       ! Enforce a maximum scale radius (minimum concentration) as required to ensure that y<1 in the cusp-NFW profile.
       radiusScale=max(radiusScale,(amplitude/densityScale/yMaximum)**(2.0d0/3.0d0))
       ! Evaluate the prompt cusp y-parameter.
       y          =+amplitude           &
            &      /densityScale        &
            &      /radiusScale **1.5d0
       ! Store prompt cusp parameters.
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspAmplitudeID,amplitude  )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspMassID     ,mass       )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID     ,y          )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWScaleID ,radiusScale)
    end do
    ! Check for convergence.
    if (.not.Values_Agree(radiusScale,radiusScalePrevious,relTol=toleranceRelative)) then
       if (self%nonConvergenceIsFatal) then
          call Error_Report('failed to converge when seeking solution for rₛ'//{introspection:location})
       else
          block
            character(len=12        ) :: label
            type     (varying_string) :: message
            errorFractional=+abs(+radiusScale-radiusScalePrevious) &
                 &          /   (+radiusScale+radiusScalePrevious) &
                 &          /0.5d0
            if (errorFractional > errorFractionalMaximum) then
               write (label,'(e12.6)') errorFractional
               message='failed to converge when seeking solution for rₛ (relative error = '//label//')'
               call Warn(message)
               errorFractionalMaximum=errorFractional
            end if
          end block
       end if
    end if    
    return
  end subroutine darkMatterProfilePromptCuspsNodeTreeInitialize

  double precision function timeCollapseRoot(timeCollapse)
    !!{
    Root function used to find the time of collapse.
    !!}
    implicit none
    double precision, intent(in   ) :: timeCollapse

    timeCollapseRoot=+self_%sigma         (0,timeCollapse) &
         &           -      sigma0Collapse
    return
  end function timeCollapseRoot

  double precision function concentrationTargetRoot(concentration)
    !!{
    Implements the equation in footnote~9 of \cite{delos_cusp-halo_2025}. Used in solving for the minimum allowed concentration in
    a cusp-NFW density profile.    
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: concentration

    concentrationTargetRoot=concentrationFactorTarget-384.0d0*Pi*(asinh(sqrt(concentration/2.0d0))-sqrt(concentration/(2.0d0+concentration)))**2/concentration**3
    return
  end function concentrationTargetRoot
  
  double precision function darkMatterProfilePromptCuspsNodeSigma(self,j,time) result(sigma)
    !!{
    Evaluate the integral
    \begin{equation}
      \sigma_j^2(t) = \int_0^\infty \frac{\mathrm{d}k}{k} \mathcal{P}(k,t) k^{2j},
    \end{equation}
    where $\mathcal{P}(k) = k^3 P(k) / 2 \pi^2$ is the dimensionless form of the power spectrum.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Comparison , only : Values_Agree
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout), target       :: self
    integer                                                   , intent(in   )               :: j
    double precision                                          , intent(in   )               :: time
    double precision                                          , parameter                   :: wavenumberPhysicalMinimum           =1.0d-2 , wavenumberPhysicalMaximum           =1.0d+2, &
         &                                                                                     toleranceRelative                   =1.0d-3
    type            (integrator                              ), save                        :: integrator_
    logical                                                   , save                        :: integratorInitialized               =.false.
    !$omp threadprivate(integrator_,integratorInitialized)
    double precision                                          , allocatable  , dimension(:) :: sigmaTmp
    double precision                                                                        :: wavenumberPhysicalLogarithmicMinimum        , wavenumberPhysicalLogarithmicMaximum, &
         &                                                                                     sigmaPrevious

    ! Initialize an integrator if necessary.
    if (.not.integratorInitialized) then       
       integrator_          =integrator(integrand,toleranceRelative=toleranceRelative)
       integratorInitialized=.true.
    end if
    ! Is growth is not wavenumber dependent, check if we have computed the value of sigma for this j.
    if (.not.self%growthIsWavenumberDependent) then
       if (.not.allocated(self%sigma_)) then
          allocate(self%sigma_(0:j))
          self%sigma_                    =-huge(0.0d0)
       else if (ubound(self%sigma_,dim=1) < j) then
          call move_alloc(self%sigma_,sigmaTmp)
          allocate(self%sigma_(0:j))
          self%sigma_                    =-huge(0.0d0)
          self%sigma_(0:size(sigmaTmp)-1)= sigmaTmp
       end if
    end if
    ! Determine if we need to perform the integral.
    if (self%growthIsWavenumberDependent .or. self%sigma_(j) < 0.0d0) then
       ! Make an initial guess for the range of wavenumbers over which to integrate the power spectrum.
       if (self%growthIsWavenumberDependent) then
          time_                             =  time
          expansionFactor_                  =  self%cosmologyFunctions_%expansionFactor(time            )
       else
          expansionFactor_                  =  1.0d0
          time_                             =  self%cosmologyFunctions_%cosmicTime     (expansionFactor_)
       end if
       self_                                => self
       j_                                   =  j
       wavenumberPhysicalLogarithmicMinimum =  log(wavenumberPhysicalMinimum)
       wavenumberPhysicalLogarithmicMaximum =  log(wavenumberPhysicalMaximum)
       sigmaPrevious                        =  -huge(0.0d0)
       sigma                                =  +     0.0d0
       ! Expand the range over wavenumbers integrated over until the integral is sufficiently well converged.
       do while (.not.Values_Agree(sigma,sigmaPrevious,relTol=toleranceRelative))
          sigmaPrevious                       =+sigma
          sigma                               =+sqrt(                                                            &
               &                                     integrator_%integrate(                                      &
               &                                                           wavenumberPhysicalLogarithmicMinimum, &
               &                                                           wavenumberPhysicalLogarithmicMaximum  &
               &                                                          )                                      &
               &                                    )
          wavenumberPhysicalLogarithmicMinimum=+wavenumberPhysicalLogarithmicMinimum-1.0d0
          wavenumberPhysicalLogarithmicMaximum=+wavenumberPhysicalLogarithmicMaximum+1.0d0
       end do
       if (.not.self%growthIsWavenumberDependent) self%sigma_(j)=sigma
    end if
    ! If growth is wavenumber independent simply scale sigma to the current time (accounting for both linear growth and the fact
    ! that these σ are defined in physical, not comoving coordinates).
    if (.not.self%growthIsWavenumberDependent)                           &
         & sigma=+self                    %sigma_         (     j   )    &
         &       *self%linearGrowth_      %value          (time=time)    &
         &       /self%cosmologyFunctions_%expansionFactor(time=time)**j
    return
  end function darkMatterProfilePromptCuspsNodeSigma
  
  double precision function integrand(wavenumberPhysicalLogarithmic)
    !!{
    Integrand used to compute the quantity $\sigma^2_j(t)$ in the prompt cusps model of \cite{delos_cusp-halo_2025}.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumberPhysicalLogarithmic
    double precision                :: wavenumberComoving           , wavenumberPhysical
    
    wavenumberPhysical=+exp(wavenumberPhysicalLogarithmic)
    wavenumberComoving=+expansionFactor_   &
         &             *wavenumberPhysical
    integrand         =+self_%powerSpectrum_%powerDimensionless(wavenumberComoving        ,time_) &
         &             *                                        wavenumberPhysical**(2*j_)
    return
  end function integrand

  subroutine darkMatterProfilePromptCuspsNodeInitialize(self,node)
    !!{
    Compute the rate of growth of dark matter profile scale radius assuming a constant growth rate.
    !!}
    use :: Display         , only : displayBlue       , displayGreen                  , displayYellow, displayBold, &
         &                          displayReset
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    class           (nodeComponentBasic                      )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile          )               , pointer :: darkMatterProfile, darkMatterProfileParent
    double precision                                                                   :: timeInterval

    ! Set the growth rate for the scale radius.
    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       call Error_Report(                                                                                                                                                                                            &
            &            displayBold()//'darkMatterProfile'//displayReset()//' component must be created prior to initialization cusp-NFW scale radius interpolation'                                   //char(10)// &
            &            '   For example, by using the following nodeOperator'                                                                                                                          //char(10)// &
            &            '    <'//displayBlue()//'nodeOperator'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"darkMatterProfileScaleSet"'//displayReset()//'/>'          // &
            &            {introspection:location}                                                                                                                                                                    &
            &           )
    class default
       if (node%isPrimaryProgenitor()) then
          ! Node is the primary progenitor, so compute the scale radius growth rate.
          basic        =>  node              %basic()
          basicParent  =>  node       %parent%basic()
          timeInterval =  +basicParent       %time () &
               &          -basic             %time ()
          if (timeInterval > 0.0d0) then
             darkMatterProfileParent => node%parent%darkMatterProfile()
             call darkMatterProfile%floatRank0MetaPropertySet(                                                                                 &
                  &                                              self                   %promptCuspNFWGrowthRateID                           , &
                  &                                           +(                                                                               &
                  &                                             +darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWScaleID)  &
                  &                                             -darkMatterProfile      %floatRank0MetaPropertyGet(self%promptCuspNFWScaleID)  &
                  &                                            )                                                                               &
                  &                                           /                          timeInterval                                          &
                  &                                          )
          else
             ! Time interval is non-positive - assume zero growth rate.
             call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWGrowthRateID,0.0d0)
          end if
       else
          ! Node is a non-primary progenitor - assume zero growth rate.
          call    darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWGrowthRateID,0.0d0)
       end if
    end select
    return
  end subroutine darkMatterProfilePromptCuspsNodeInitialize

  subroutine darkMatterProfilePromptCuspsNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the scale radius
    growth rate of {\normalfont \ttfamily node} to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfilePromptCusps), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile, darkMatterProfileParent

    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID         ,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWYID         ))
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWGrowthRateID,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWGrowthRateID))
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWScaleID     ,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWScaleID     ))
    return
  end subroutine darkMatterProfilePromptCuspsNodePromote

  subroutine darkMatterProfilePromptCuspsSolveAnalytics(self,node,time)
    !!{
    Compute the value of the $y$-parameter in the prompt cusp.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rootFinder        , rangeExpandMultiplicative     , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: time
    class           (nodeComponentBasic                      ), pointer       :: basic                    , basicParent
    class           (nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile        , darkMatterProfileParent
    type            (rootFinder                              ), save          :: finder
    logical                                                   , save          :: finderInitialized=.false.
    !$omp threadprivate(finder,finderInitialized)
    double precision                                                          :: densityScale             , y

    ! Primary progenitors do not evolve.
    if (.not.node%isPrimaryProgenitor()) return
    ! Initialize the root finder.
    if (.not.finderInitialized) then
       finder=rootFinder(                                                             &
            &            rootFunction                 =densityNormalizationRoot     , &
            &            toleranceRelative            =1.0d-6                       , &
            &            rangeExpandUpward            =2.0d+0                       , &
            &            rangeExpandDownward          =0.5d+0                       , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
            &           )
      finderInitialized=.true.
    end if
    ! Solve for the density normalization.
    basic                   =>  node                                        %basic                    (                              )
    basicParent             =>  node                   %parent              %basic                    (                              )
    darkMatterProfile       =>  node                                        %darkMatterProfile        (                              )
    darkMatterProfileParent =>  node                   %parent              %darkMatterProfile        (                              )
    massHalo_               =   basic                                       %mass                     (                              )
    amplitudeCusp_          =  +darkMatterProfile                           %floatRank0MetaPropertyGet(self%promptCuspAmplitudeID    )
    y                       =  +darkMatterProfile                           %floatRank0MetaPropertyGet(self%promptCuspNFWYID         )
    radiusScale_            =  +darkMatterProfileParent                     %floatRank0MetaPropertyGet(self%promptCuspNFWScaleID     ) &
         &                     +(                                                                                                      &
         &                       +                                           time                                                      &
         &                       -basicParent                               %time                     (                              ) &
         &                      )                                                                                                      &
         &                     *darkMatterProfile                           %floatRank0MetaPropertyGet(self%promptCuspNFWGrowthRateID)
    concentration_          =  +self                   %darkMatterHaloScale_%radiusVirial             (node                          ) &
         &                     /                                             radiusScale_
    densityScale            =  +massHalo_                                                                                      &
         &                     /radiusScale_**3                                                                                &
         &                     /4.0d0                                                                                          &
         &                     /Pi                                                                                             &
         &                     /(                                                                                              &
         &                       + 2.0d0                       *asinh(sqrt(concentration_            )/                y     ) &
         &                       -(2.0d0-y**2)/sqrt(1.0d0-y**2)*atanh(sqrt(concentration_*(1.0d0-y**2)/(concentration_+y**2))) &
         &                       -sqrt(concentration_*(concentration_+y**2))/(1.0d0+concentration_)                            &
         &                      )
    densityScale            =   finder%find(rootGuess=densityScale)
    ! Compute the cusp y parameter.
    y                     =min(                      &
         &                     +amplitudeCusp_       &
         &                     /densityScale         &
         &                     /radiusScale_**1.5d0, &
         &                     +yMaximum             &
         &                    )
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID    ,y           )
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWScaleID,radiusScale_)
    return
  end subroutine darkMatterProfilePromptCuspsSolveAnalytics

  double precision function densityNormalizationRoot(densityScale)
    !!{
    Root function used in finding the density normalization for cusp-NFW density profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: densityScale
    double precision                :: y

    y                 =min(                      &
         &                 +amplitudeCusp_       &
         &                 /densityScale         &
         &                 /radiusScale_**1.5d0, &
         &                 +yMaximum             &
         &                )
    densityNormalizationRoot=+massHalo_                                                                                      &
         &                   /radiusScale_**3                                                                                &
         &                   /4.0d0                                                                                          &
         &                   /Pi                                                                                             &
         &                   /(                                                                                              &
         &                     + 2.0d0                       *asinh(sqrt(concentration_            )/                y     ) &
         &                     -(2.0d0-y**2)/sqrt(1.0d0-y**2)*atanh(sqrt(concentration_*(1.0d0-y**2)/(concentration_+y**2))) &
         &                     -sqrt(concentration_*(concentration_+y**2))/(1.0d0+concentration_)                            &
         &                    )                                                                                              &
         &                   -densityScale
    return
  end function densityNormalizationRoot

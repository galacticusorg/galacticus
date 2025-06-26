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

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Power_Spectra          , only : powerSpectrumClass
  
  !![
  <nodeOperator name="nodeOperatorDarkMatterProfilePromptCusps">
   <description>
    A node operator class that evaluates the properties of prompt cusps following the model of \cite{delos_cusp-halo_2025}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfilePromptCusps
     !!{     
     A node operator class that evaluates the properties of prompt cusps following the model of \cite{delos_cusp-halo_2025}.
     !!}
     private
     class           (powerSpectrumClass      ), pointer :: powerSpectrum_       => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: alpha                         , beta                 , &
          &                                                 kappa                         , C                    , &
          &                                                 p
     integer                                             :: promptCuspMassID              , promptCuspAmplitudeID, &
          &                                                 promptCuspNFWYID
   contains
     !![
     <methods>
       <method method="sigma" description="Evaluate $\sigma_j^2 = \int_0^\infty \frac{\mathrm{d}k}{k} \mathcal{P}(k,t) k^{2j}$ where $\mathcal{P}(k) = k^3 P(k) / 2 \pi^2$ is the dimensionless form of the power spectrum."/>
     </methods>
     !!]
     final     ::                       darkMatterProfilePromptCuspsConstructorDestructor
     procedure :: nodeTreeInitialize => darkMatterProfilePromptCuspsNodeTreeInitialize
     procedure :: sigma              => darkMatterProfilePromptCuspsNodeSigma
  end type nodeOperatorDarkMatterProfilePromptCusps
  
  interface nodeOperatorDarkMatterProfilePromptCusps
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfilePromptCusps} node operator class.
     !!}
     module procedure darkMatterProfilePromptCuspsConstructorParameters
     module procedure darkMatterProfilePromptCuspsConstructorInternal
  end interface nodeOperatorDarkMatterProfilePromptCusps

  ! Submodule-scope variables used in root-finding.
  class           (nodeOperatorDarkMatterProfilePromptCusps), pointer :: self_
  double precision                                                    :: sigma0Collapse  , time_ ,&
       &                                                                 expansionFactor_
  integer                                                             :: j_
  !$omp threadprivate(self_,sigma0Collapse,time_,expansionFactor_,j_)
  
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
    class           (powerSpectrumClass                      ), pointer       :: powerSpectrum_
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    double precision                                                          :: alpha               , beta, &
         &                                                                       kappa               , C   , &
         &                                                                       p

    !![
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
    <objectBuilder class="powerSpectrum"      name="powerSpectrum_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfilePromptCusps(alpha,beta,kappa,C,p,powerSpectrum_,cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="powerSpectrum_"      />
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function darkMatterProfilePromptCuspsConstructorParameters

  function darkMatterProfilePromptCuspsConstructorInternal(alpha,beta,kappa,C,p,powerSpectrum_,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfilePromptCusps} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfilePromptCusps)                        :: self
    class           (powerSpectrumClass                      ), intent(in   ), target :: powerSpectrum_
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                          , intent(in   )         :: alpha               , beta, &
         &                                                                               kappa               , C   , &
         &                                                                               p
    !![
    <constructorAssign variables="alpha, beta, kappa, C, p, *powerSpectrum_, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="promptCuspAmplitude" id="self%promptCuspAmplitudeID" isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspMass"      id="self%promptCuspMassID"      isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWY"      id="self%promptCuspNFWYID"      isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function darkMatterProfilePromptCuspsConstructorInternal

  subroutine darkMatterProfilePromptCuspsConstructorDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorDarkMatterProfilePromptCusps} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfilePromptCusps), intent(inout) :: self

    !![
    <objectDestructor name="self%powerSpectrum_"      />
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine darkMatterProfilePromptCuspsConstructorDestructor

  subroutine darkMatterProfilePromptCuspsNodeTreeInitialize(self,node)
    !!{
    Initialize dark matter profile prompt cusp properties.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    use :: Lambert_Ws              , only : Lambert_W0
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rootFinder        , rangeExpandMultiplicative     , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    class           (nodeComponentBasic                      )               , pointer :: basic
    class           (nodeComponentDarkMatterProfile          )               , pointer :: darkMatterProfile
    type            (rootFinder                              )               , save    :: finder
    logical                                                                  , save    :: finderInitialized=.false.
    !$omp threadprivate(finder,finderInitialized)
    double precision                                                                   :: sigma0                   , sigma2             , &
         &                                                                                densityMean              , densityMeanCollapse, &
         &                                                                                sigma2Collapse           , timeCollapse       , &
         &                                                                                amplitude                , mass               , &
         &                                                                                y                        , radiusScale        , &
         &                                                                                concentration            , densityScale
    
    ! Only initialize prompt cusp properties in leaf nodes.
    if (associated(node%firstChild)) return
    ! Initialize the root finder.
    if (.not.finderInitialized) then
       finder=rootFinder(&
            &            rootFunction                 =timeCollapseRoot             , &
            &            toleranceRelative            =1.0d-6                       , &
            &            rangeExpandUpward            =2.0d+0                       , &
            &            rangeExpandDownward          =0.5d+0                       , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
            &           )
       finderInitialized=.true.
    end if
    ! Evaluate the required integrals over the power spectrum at the time of this node.
    self_             =>  self
    basic             =>  node                                  %basic               (                   )
    darkMatterProfile =>  node                                  %darkMatterProfile   (                   )
    sigma0            =   self                                  %sigma               (0,     basic%time())
    sigma2            =   self                                  %sigma               (2,     basic%time())
    densityMean       =   self             %cosmologyFunctions_ %matterDensityEpochal(  time=basic%time())
    radiusScale       =   darkMatterProfile%scale()
    concentration     =  +self             %darkMatterHaloScale_%radiusVirial        (       node        ) &
         &               /                                       radiusScale
    densityScale      =  +basic                                 %mass                (                   ) &
         &               /4.0d0                                                                            &
         &               /Pi                                                                               &                                 
         &               /radiusScale**3                                                                   &
         &               /(                                                                                &
         &                 +              log(1.0d0+concentration)                                         &
         &                 -concentration/   (1.0d0+concentration)                                         &
         &                )
    ! Find the collapse time of this node.
    sigma0Collapse     =+            (2.0d0*self%p-1.0d0)*self%kappa/3.0d0              &
         &              /Lambert_W0(                                                    &
         &                          +(2.0d0*self%p-1.0d0)*self%kappa/3.0d0              &
         &                          *(self%beta*self%C**2/self%alpha**2)**(1.0d0/3.0d0) &
         &                          *(                                                  &
         &                            +exp(self%kappa/sigma0)                           &
         &                            *basic%mass()                                     &
         &                            /densityMean                                      &
         &                            /(sigma0/sigma2)**1.5d0                           &
         &                           )**((2.0d0*self%p-1.0d0)/3.0d0)                    &
         &                         )
    timeCollapse       =finder                    %find(                  rootGuess=basic%time        ())
    sigma2Collapse     =self                      %sigma               (2,time     =      timeCollapse  )
    densityMeanCollapse=self  %cosmologyFunctions_%matterDensityEpochal(  time     =      timeCollapse  )
    ! Evaluate the prompt cusp parameters.
    amplitude=self%alpha*(self%alpha/self%C/self%beta**self%p)**(1.0d0/(2.0d0*self%p-1.0d0))*densityMeanCollapse*sigma0Collapse**((9.0d0-6.0d0*self%p)/(4.0d0-8.0d0*self%p))/sigma2Collapse**0.75d0
    mass     =self%beta *(self%alpha/self%C/self%beta**self%p)**(2.0d0/(2.0d0*self%p-1.0d0))*densityMeanCollapse*sigma0Collapse**((9.0d0-6.0d0*self%p)/(2.0d0-4.0d0*self%p))/sigma2Collapse**1.50d0
    y        =+amplitude                        &
         &    /densityScale                     &
         &    /darkMatterProfile%scale()**1.5d0
    ! Store prompt cusp parameters.
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspAmplitudeID,amplitude)
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspMassID     ,mass     )
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID     ,y        )
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
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout), target :: self
    integer                                                   , intent(in   )         :: j
    double precision                                          , intent(in   )         :: time
    double precision                                          , parameter             :: wavenumberPhysicalMinimum           =1.0d-2 , wavenumberPhysicalMaximum           =1.0d+2, &
         &                                                                               toleranceRelative                   =1.0d-3
    type            (integrator                              ), save                  :: integrator_
    logical                                                   , save                  :: integratorInitialized               =.false.
    !$omp threadprivate(integrator_,integratorInitialized)
    double precision                                                                  :: wavenumberPhysicalLogarithmicMinimum        , wavenumberPhysicalLogarithmicMaximum, &
         &                                                                               sigmaPrevious

    ! Initialize an integrator if necessary.
    if (.not.integratorInitialized) then       
       integrator_          =integrator(integrand,toleranceRelative=toleranceRelative)
       integratorInitialized=.true.
    end if
    ! Make an initial guess for the range of wavenumbers over which to integrate the power spectrum.
    self_                                => self
    j_                                   =  j
    time_                                =  time
    expansionFactor_                     =  self%cosmologyFunctions_%expansionFactor(time)
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
    integrand         =+self_%powerSpectrum_%powerDimensionless(wavenumberComoving       ,time_) &
         &             *                                        wavenumberPhysical**(2*j_)
    return
  end function integrand

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

  !+    Contributions to this file made by: Andrew Benson.

  !!{
  An implementation of cosmological velocity field computed using a filtered power spectrum.
  !!}
  use :: Cosmology_Parameters           , only : cosmologyParametersClass
  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Correlation_Functions_Two_Point, only : correlationFunctionTwoPointClass
  use :: Linear_Growth                  , only : linearGrowthClass
  use :: Power_Spectra                  , only : powerSpectrumClass
  use :: Power_Spectrum_Window_Functions, only : powerSpectrumWindowFunctionClass

  !![
  <cosmologicalVelocityField name="cosmologicalVelocityFieldFilteredPower">
   <description>
    Cosmological velocity field computed by filtering the linear theory power spectrum. The growth factor for velocities is
    $D_\mathrm{v}(t) = a(t) D(t) H(t) f(t)$, where $D(t)$ is the usual growth factor for density, and $f(t) = \mathrm{d}\log D
    / \mathrm{d} \log a$. Note that the factor of $D(t)$ does not explicitly appear in expressions for the velocity dispersion
    since it is included in the linear theory power spectrum appearing in those expressions.
   </description>
  </cosmologicalVelocityField>
  !!]
  type, extends(cosmologicalVelocityFieldClass) :: cosmologicalVelocityFieldFilteredPower
     !!{
     A cosmological mass variance class computing variance from a filtered power spectrum.
     !!}
     private
     class           (cosmologyParametersClass        ), pointer :: cosmologyParameters_         => null()
     class           (cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_          => null()
     class           (powerSpectrumClass              ), pointer :: powerSpectrum_               => null()
     class           (linearGrowthClass               ), pointer :: linearGrowth_                => null()
     class           (powerSpectrumWindowFunctionClass), pointer :: powerSpectrumWindowFunction_ => null()
     class           (correlationFunctionTwoPointClass), pointer :: correlationFunctionTwoPoint_ => null()
     double precision                                            :: wavenumberMaximum
   contains
     !![
     <methods>
      <method description="Compute the function $\sigma_j^2(m) = {1 \over 2 \pi^2} \int_0^\infty \mathrm{d}k k^{2+2j} P(k) W^2[kR(m)]$, e.g. \cite[][unnumbered equation following eqn.~8]{sheth_peculiar_2001}."                    method="sigmaJ"        />
      <method description="Compute the peak correction term for the velocity dispersion of halos of given {\normalfont \ttfamily mass}, e.g. \cite[][eqn.~8]{sheth_peculiar_2001}, and \cite[][eqn. 4.26]{bardeen_statistics_1986}." method="peakCorrection"/>
     </methods>
     !!]
     final     ::                                     filteredPowerDestructor
     procedure :: velocityRadialMeanPairwise       => filteredPowerVelocityRadialMeanPairwise
     procedure :: velocityDispersion1D             => filteredPowerVelocityDispersion1D
     procedure :: velocityDispersion1DHaloPairwise => filteredPowerVelocityDispersion1DHaloPairwise
     procedure :: sigmaJ                           => filteredPowerSigmaJ
     procedure :: peakCorrection                   => filteredPowerPeakCorrection
  end type cosmologicalVelocityFieldFilteredPower

  interface cosmologicalVelocityFieldFilteredPower
     !!{
     Constructors for the \refClass{cosmologicalVelocityFieldFilteredPower} cosmological mass variance class.
     !!}
     module procedure filteredPowerConstructorParameters
     module procedure filteredPowerConstructorInternal
  end interface cosmologicalVelocityFieldFilteredPower

contains

  function filteredPowerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{cosmologicalVelocityFieldFilteredPower} cosmological mass variance class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (cosmologicalVelocityFieldFilteredPower)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyParametersClass              ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class           (powerSpectrumClass                    ), pointer       :: powerSpectrum_
    class           (powerSpectrumWindowFunctionClass      ), pointer       :: powerSpectrumWindowFunction_
    class           (linearGrowthClass                     ), pointer       :: linearGrowth_
    class           (correlationFunctionTwoPointClass      ), pointer       :: correlationFunctionTwoPoint_
    double precision                                                        :: wavenumberMaximum

    !![
    <inputParameter>
      <name>wavenumberMaximum</name>
      <source>parameters</source>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>The maximum wavenumber to which to integrate the power spectrum. By default this in infinite. It can be useful to set a smaller value for power spectra with small-scale cut offs to avoid convergence issues in the integrals.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"         name="cosmologyParameters_"         source="parameters"/>
    <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters"/>
    <objectBuilder class="powerSpectrum"               name="powerSpectrum_"               source="parameters"/>
    <objectBuilder class="powerSpectrumWindowFunction" name="powerSpectrumWindowFunction_" source="parameters"/>
    <objectBuilder class="linearGrowth"                name="linearGrowth_"                source="parameters"/>
    <objectBuilder class="correlationFunctionTwoPoint" name="correlationFunctionTwoPoint_" source="parameters"/>
    !!]
    self=filteredPowerConstructorInternal(wavenumberMaximum,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,powerSpectrum_,powerSpectrumWindowFunction_,correlationFunctionTwoPoint_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"        />
    <objectDestructor name="cosmologyFunctions_"         />
    <objectDestructor name="linearGrowth_"               />
    <objectDestructor name="powerSpectrum_"              />
    <objectDestructor name="powerSpectrumWindowFunction_"/>
    <objectDestructor name="correlationFunctionTwoPoint_"/>
    !!]
    return
  end function filteredPowerConstructorParameters

  function filteredPowerConstructorInternal(wavenumberMaximum,cosmologyParameters_,cosmologyFunctions_,linearGrowth_,powerSpectrum_,powerSpectrumWindowFunction_,correlationFunctionTwoPoint_) result(self)
    !!{
    Internal constructor for the \refClass{cosmologicalVelocityFieldFilteredPower} linear growth class.
    !!}
    implicit none
    type            (cosmologicalVelocityFieldFilteredPower)                        :: self
    class           (cosmologyParametersClass              ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class           (powerSpectrumClass                    ), intent(in   ), target :: powerSpectrum_
    class           (powerSpectrumWindowFunctionClass      ), intent(in   ), target :: powerSpectrumWindowFunction_
    class           (linearGrowthClass                     ), intent(in   ), target :: linearGrowth_
    class           (correlationFunctionTwoPointClass      ), intent(in   ), target :: correlationFunctionTwoPoint_
    double precision                                        , intent(in   )         :: wavenumberMaximum
    !![
    <constructorAssign variables="wavenumberMaximum, *cosmologyParameters_, *cosmologyFunctions_, *linearGrowth_, *powerSpectrum_, *powerSpectrumWindowFunction_, *correlationFunctionTwoPoint_"/>
    !!]

    return
  end function filteredPowerConstructorInternal

  subroutine filteredPowerDestructor(self)
    !!{
    Destructor for the \refClass{cosmologicalVelocityFieldFilteredPower} linear growth class.
    !!}
    implicit none
    type   (cosmologicalVelocityFieldFilteredPower), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"        />
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%linearGrowth_"               />
    <objectDestructor name="self%powerSpectrum_"              />
    <objectDestructor name="self%powerSpectrumWindowFunction_"/>
    <objectDestructor name="self%correlationFunctionTwoPoint_"/>
    !!]
    return
  end subroutine filteredPowerDestructor

  double precision function filteredPowerVelocityRadialMeanPairwise(self,separation,time,includeHubbleFlow)
    !!{
    Return the mean radial velocity (averaged over all positions) at a given {\normalfont \ttfamily separation} and
    {\normalfont \ttfamily time}. If {\normalfont \ttfamily includeHubbleFlow} is {\normalfont \ttfamily true} then the Hubble
    flow is included, otherwise only the peculiar component of the mean radial velocity is computed.
    !!}
    implicit none
    class           (cosmologicalVelocityFieldFilteredPower), intent(inout) :: self
    double precision                                        , intent(in   ) :: separation        , time
    logical                                                 , intent(in   ) :: includeHubbleFlow
    double precision                                                        :: separationComoving, hubbleFlowTerm

    separationComoving=+separation                                     &
         &             /self%cosmologyFunctions_%expansionFactor(time)
    if (includeHubbleFlow) then
       hubbleFlowTerm=1.0d0
    else
       hubbleFlowTerm=0.0d0
    end if
    filteredPowerVelocityRadialMeanPairwise=+(                                                                                                   &
         &                                    +  hubbleFlowTerm                                                                                  &
         &                                    -  2.0d0                                                                                           &
         &                                    /  3.0d0                                                                                           &
         &                                    *  self%linearGrowth_               %logarithmicDerivativeExpansionFactor(                   time) &
         &                                    *  self%correlationFunctionTwoPoint_%correlationVolumeAveraged           (separationComoving,time) &
         &                                    /(                                                                                                 &
         &                                      +1.0d0                                                                                           &
         &                                      +self%correlationFunctionTwoPoint_%correlation                         (separationComoving,time) &
         &                                     )                                                                                                 &
         &                                   )                                                                                                   &
         &                                  *    self%cosmologyFunctions_         %HubbleParameterEpochal              (                   time) &
         &                                  *                                                                            separation
    return
  end function filteredPowerVelocityRadialMeanPairwise

  double precision function filteredPowerVelocityDispersion1D(self,mass,time)
    !!{
    Return the 1-D dispersion of the velocity field smoothed over in a spherical region containing the given {\normalfont
    \ttfamily mass} on average.
    !!}
    implicit none
    class           (cosmologicalVelocityFieldFilteredPower), intent(inout) :: self
    double precision                                        , intent(in   ) :: mass, time
    
    filteredPowerVelocityDispersion1D=+self%cosmologyFunctions_%hubbleParameterEpochal              (     time   ) &
         &                            *self%linearGrowth_      %logarithmicDerivativeExpansionFactor(     time   ) &
         &                            *self%cosmologyFunctions_%expansionFactor                     (     time   ) &
         &                            *self                    %sigmaJ                              (mass,time,-1) &
         &                            /sqrt(3.0d0)
    return
  end function filteredPowerVelocityDispersion1D

  double precision function filteredPowerSigmaJ(self,mass,time,j)
    !!{
    Compute the quantity:
    \begin{equation}
     \sigma_j^2(m) = {1 \over 2 \pi^2} \int_0^\infty \mathrm{d}k k^{2+2j} P(k) W^2[kR(m)],
    \end{equation}
    e.g. \cite[][unnumbered equation following eqn.~8]{sheth_peculiar_2001}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator, GSL_Integ_Gauss15
    implicit none
    implicit none
    class           (cosmologicalVelocityFieldFilteredPower), intent(inout) :: self
    double precision                                        , intent(in   ) :: mass             , time
    integer                                                 , intent(in   ) :: j
    type            (integrator                            )                :: integrator_
    double precision                                                        :: wavenumberMinimum, wavenumberMaximum, &
         &                                                                     sigmaJSquared    , radiusTopHat

    integrator_        = integrator(powerIntegrand,toleranceRelative=1.0d-5,integrationRule=GSL_Integ_Gauss15)
    radiusTopHat       =+(                                        &
         &                +(                                      &
         &                  +3.0d0                                &
         &                  /4.0d0                                &
         &                  /Pi                                   &
         &                 )                                      &
         &                *mass                                   &
         &                /self%cosmologyParameters_%OmegaMatter    () &
         &                /self%cosmologyParameters_%densityCritical() &
         &               )**(1.0d0/3.0d0)
    wavenumberMinimum  =+    0.0d0
    wavenumberMaximum  = min(1.0d3/radiusTopHat,self%powerSpectrumWindowFunction_%wavenumberMaximum(mass))
    sigmaJSquared      =+integrator_%integrate(wavenumberMinimum,wavenumberMaximum) &
         &              /2.0d0                                                      &
         &              /Pi**2
    filteredPowerSigmaJ=+sqrt(sigmaJSquared)
    return

  contains
    
    double precision function powerIntegrand(wavenumber)
      !!{
      Integrand appearing in the definition of $\sigma_j(m)$, e.g. 
      \cite[][unnumbered equation following eqn.~8]{sheth_peculiar_2001}.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      powerIntegrand=+  self%powerSpectrum_              %power(wavenumber,time)        &
           &         *(                                                                 &
           &           +self%powerSpectrumWindowFunction_%value(wavenumber,mass)        &
           &           *                                        wavenumber      **(1+j) &
           &          )**2    
      return
    end function powerIntegrand
    
  end function filteredPowerSigmaJ
  
  double precision function filteredPowerVelocityDispersion1DHaloPairwise(self,mass1,mass2,separation,time)
    !!{
    Compute the linear theory velocity dispersion between two patches along line of sight connecting them (this is therefore relevant to the 1D
    dispersion), including the effect of correlated velocities between the patches. The integral evaluated here is the sum of the
    integral for the linear velocity dispersion of each halo \citep[][eqn.~8]{sheth_peculiar_2001} and the correlation integral
    \citep[][eqns.~28 \& 29]{sheth_linear_2001}. We evaluate the integrals together (instead of separately) as at small radii
    they cancel to high order---evaluating them separately leads to inaccurate estimates of the velocity dispersion.
    !!}
    use :: Error                   , only : Warn      , errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator, GSL_Integ_Gauss15
    implicit none
    class           (cosmologicalVelocityFieldFilteredPower), intent(inout) :: self
    double precision                                        , intent(in   ) :: mass1                    , mass2             , &
         &                                                                     separation               , time
    double precision                                        , parameter     :: radiusMinimum    =1.0d-3
    double precision                                        , parameter     :: radiusMaximum    =1.0d+3
    logical                                                 , save          :: warningIssued    =.false.
    type            (integrator                            )                :: integrator_
    double precision                                                        :: wavenumberMinimum        , wavenumberMaximum , &
         &                                                                     peakCorrection1          , peakCorrection2   , &
         &                                                                     radiusTopHat             , separationComoving
    integer                                                                 :: status

    integrator_                                  = integrator(velocityDispersionIntegrand,toleranceRelative=1.0d-2,integrationRule=GSL_Integ_Gauss15)
    radiusTopHat                                 =+(                                         &
         &                                          +(                                       &
         &                                            +3.0d0                                 &
         &                                            /4.0d0                                 &
         &                                            /Pi                                    &
         &                                           )                                       &
         &                                          *min(mass1,mass2)                        &
         &                                          /self%cosmologyParameters_%OmegaMatter    ()  &
         &                                          /self%cosmologyParameters_%densityCritical()  &
         &                                         )**(1.0d0/3.0d0)
    separationComoving                           =+separation                                     &
         &                                        /self%cosmologyFunctions_%expansionFactor(time)
    wavenumberMinimum                            =+                               radiusMinimum/max(radiusTopHat,separationComoving)
    wavenumberMaximum                            =+min(self%wavenumberMaximum,min(radiusMaximum/min(radiusTopHat,separationComoving),self%powerSpectrumWindowFunction_%wavenumberMaximum(min(mass1,mass2))))
    peakCorrection1                              =self%peakCorrection(mass1,time)
    peakCorrection2                              =self%peakCorrection(mass2,time)
    filteredPowerVelocityDispersion1DHaloPairwise=+sqrt(                                                                   &
         &                                              +integrator_%integrate(wavenumberMinimum,wavenumberMaximum,status) &
         &                                              /2.0d0                                                             &
         &                                              /Pi**2                                                             &
         &                                             )                                                                   &
         &                                        *self%cosmologyFunctions_%hubbleParameterEpochal              (time)     &
         &                                        *self%linearGrowth_      %logarithmicDerivativeExpansionFactor(time)     &
         &                                        *self%cosmologyFunctions_%expansionFactor                     (time)
    if (status /= errorStatusSuccess .and. .not.warningIssued) then
       call Warn('integration failure in halo-halo velocity dispersion calculation'//{introspection:location})
       warningIssued=.true.
    end if
    return

  contains
    
    double precision function velocityDispersionIntegrand(wavenumber)
      !!{
      Integrand for the velocity dispersion.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber
      
      velocityDispersionIntegrand=+    self%powerSpectrum_              %power(wavenumber                   ,time)                        &
           &                      *(                                                                                                      &
           &                        +(                                                                                                    &
           &                          +self%powerSpectrumWindowFunction_%value(wavenumber                   ,mass1)**2*peakCorrection1**2 &
           &                          +self%powerSpectrumWindowFunction_%value(wavenumber                   ,mass2)**2*peakCorrection2**2 &
           &                         )                                                                                                    &
           &                        /3.0d0                                                                                                & ! Convert from 3D to 1D dispersion.
           &                        -2.0d0                                                                                                &
           &                        *  self%powerSpectrumWindowFunction_%value(wavenumber                   ,mass1)   *peakCorrection1    &
           &                        *  self%powerSpectrumWindowFunction_%value(wavenumber                   ,mass2)   *peakCorrection2    &
           &                        *  K                                      (wavenumber*separationComoving      )                       &
           &                       )
      return
    end function velocityDispersionIntegrand
    
  end function filteredPowerVelocityDispersion1DHaloPairwise

  double precision function K(x)
    !!{
    Correlation window function appearing in the definition of the correlation in velocity between two patches along line of sight connecting
    them, $\psi(m_1,m_2|r)$, e.g. \cite[][eqn.~28]{sheth_linear_2001}.
    !!}
    implicit none
    double precision, intent(in   ) :: x
    
    K      =+  sin(x)    &
         &  /      x     &
         &  -2.0d0       &
         &  /      x **3 &
         &  *(           &
         &    +sin(x)    &
         &    -    x     &
         &    *cos(x)    &
         &  )
    return
  end function K

  double precision function filteredPowerPeakCorrection(self,mass,time)
    !!{
    Compute the peak correction term for the velocity dispersion of halos of given {\normalfont \ttfamily mass},
    e.g. \cite[][eqn.~8]{sheth_peculiar_2001}, and \cite[][eqn. 4.26]{bardeen_statistics_1986}.
    !!}
    implicit none
    class           (cosmologicalVelocityFieldFilteredPower), intent(inout)    :: self
    double precision                                        , intent(in   )    :: mass  , time

    filteredPowerPeakCorrection=+sqrt(                              &
         &                            +1.0d0                        &
         &                            -self%sigmaJ(mass,time, 0)**4 &
         &                            /self%sigmaJ(mass,time,+1)**2 &
         &                            /self%sigmaJ(mass,time,-1)**2 &
         &                           )
    return
  end function filteredPowerPeakCorrection

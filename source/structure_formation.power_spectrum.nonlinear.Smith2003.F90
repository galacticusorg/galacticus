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
Implements a nonlinear power spectrum class in which the nonlinear power spectrum is computed using the
algorithm of \cite{smith_stable_2003}.
!!}

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Power_Spectra       , only : powerSpectrumClass

  !![
  <powerSpectrumNonlinear name="powerSpectrumNonlinearSmith2003">
   <description>Provides a nonlinear power spectrum class in which the power spectrum is computed using the algorithm of \cite{smith_stable_2003}.</description>
  </powerSpectrumNonlinear>
  !!]
  type, extends(powerSpectrumNonlinearClass) :: powerSpectrumNonlinearSmith2003
     !!{
     A nonlinear power spectrum class in which the power spectrum is computed using the algorithm of \cite{smith_stable_2003}.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_      => null()
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_     => null()
     class           (powerSpectrumClass      ), pointer :: powerSpectrum_           => null()
     logical                                             :: includePeacockCorrection          , includeQuasiLinearPower, &
          &                                                 includeHaloPower
     double precision                                    :: timePrevious                      , wavenumberNonLinear    , &
          &                                                 effectiveIndex                    , spectralCurvature      , &
          &                                                 alpha                             , beta                   , &
          &                                                 a                                 , b                      , &
          &                                                 c                                 , gamma                  , &
          &                                                 mu                                , nu                     , &
          &                                                 f1                                , f2                     , &
          &                                                 f3
   contains
     !![
     <methods>
       <method description="Compute the fitting function coefficients at the given time." method="coefficients" />
     </methods>
     !!]
     final     ::                 smith2003Destructor
     procedure :: value        => smith2003Value
     procedure :: coefficients => smith2003Coefficients
  end type powerSpectrumNonlinearSmith2003

  interface powerSpectrumNonlinearSmith2003
     !!{
     Constructors for the \refClass{powerSpectrumNonlinearSmith2003} nonlinear power spectrum class.
     !!}
     module procedure smith2003ConstructorParameters
     module procedure smith2003ConstructorInternal
  end interface powerSpectrumNonlinearSmith2003

  ! Submodule-scope variables used in root finding and integrations.
  class           (powerSpectrumNonlinearSmith2003), pointer :: self_
  double precision                                           :: time_, radiusGaussian
  !$omp threadprivate(self_,time_,radiusGaussian)
  
contains

  function smith2003ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the smith2003 nonlinear power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (powerSpectrumNonlinearSmith2003)                        :: self
    type   (inputParameters                ), target, intent(inout) :: parameters
    class  (cosmologyFunctionsClass        ), pointer               :: cosmologyFunctions_
    class  (cosmologyParametersClass       ), pointer               :: cosmologyParameters_
    class  (powerSpectrumClass             ), pointer               :: powerSpectrum_
    logical                                                         :: includePeacockCorrection, includeQuasiLinearPower, &
          &                                                            includeHaloPower

    !![
    <inputParameter>
      <name>includePeacockCorrection</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <defaultSource>\href{https://www.roe.ac.uk/~jap/haloes/}{{\normalfont \ttfamily https://www.roe.ac.uk/\~jap/haloes/}}</defaultSource>
      <description>If true, include the correction proposed on John Peacock's \href{https://www.roe.ac.uk/~jap/haloes/}{web page}.</description>
    </inputParameter>
    <inputParameter>
      <name>includeQuasiLinearPower</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, include quasi-linear contribution to the power spectrum.</description>
    </inputParameter>
    <inputParameter>
      <name>includeHaloPower</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, include halo contribution to the power spectrum.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="powerSpectrum"       name="powerSpectrum_"       source="parameters"/>
    !!]
    ! Call the internal constructor.
    self=smith2003ConstructorInternal(includePeacockCorrection,includeQuasiLinearPower,includeHaloPower,cosmologyParameters_,cosmologyFunctions_,powerSpectrum_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="powerSpectrum_"      />
    !!]
    return
  end function smith2003ConstructorParameters

  function smith2003ConstructorInternal(includePeacockCorrection,includeQuasiLinearPower,includeHaloPower,cosmologyParameters_,cosmologyFunctions_,powerSpectrum_) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumNonlinearSmith2003} nonlinear power spectrum class.
    !!}
    implicit none
    type   (powerSpectrumNonlinearSmith2003)                        :: self
    logical                                 , intent(in   )         :: includePeacockCorrection, includeQuasiLinearPower, &
          &                                                            includeHaloPower
    class  (cosmologyParametersClass       ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class  (powerSpectrumClass             ), intent(in   ), target :: powerSpectrum_
    !![
    <constructorAssign variables="includePeacockCorrection, ,includeQuasiLinearPower, includeHaloPower, *cosmologyParameters_, *cosmologyFunctions_, *powerSpectrum_"/>
    !!]

    self%timePrevious=-1.0d0
    return
  end function smith2003ConstructorInternal

  subroutine smith2003Destructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumNonlinearSmith2003} nonlinear power spectrum class.
    !!}
    implicit none
    type(powerSpectrumNonlinearSmith2003), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%powerSpectrum_"      />
    !!]
    return
  end subroutine smith2003Destructor

  double precision function smith2003Value(self,wavenumber,time)
    !!{
    Return a nonlinear power spectrum equal using the algorithm of \cite{smith_stable_2003}.
    !!}
    use :: Cosmology_Parameters    , only : hubbleUnitsLittleH
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(powerSpectrumNonlinearSmith2003), intent(inout) :: self
    double precision                      , intent(in   ) :: time              , wavenumber
    double precision                                      :: deltaLinearSquared, deltaHaloPrimedSquared     , &
         &                                                   deltaHaloSquared  , powerSpectrumQuasiLinear   , &
         &                                                   powerSpectrumHalo , y                          , &
         &                                                   f                 , wavenumberPeacockCorrection

    ! Compute fitting function coefficients.
    call self%coefficients(time)
    ! Evaluate scaled wavenumber.
    y      =+     wavenumber          &
         &  /self%wavenumberNonlinear
    ! Evaluate the quasi-linear power spectrum.
    if (self%includeQuasiLinearPower) then
       deltaLinearSquared      =+self%powerSpectrum_%powerDimensionless(wavenumber,time)
       f                       =+y   /4.0d0                                       &
            &                   +y**2/8.0d0
       powerSpectrumQuasiLinear=+self%powerSpectrum_%power(wavenumber,time)       & ! Smith et al. (2003; eqn. 48).
            &                   *(1.0d0+           deltaLinearSquared)**self%beta &
            &                   /(1.0d0+self%alpha*deltaLinearSquared)            &
            &                   *exp(-f)
    else
       powerSpectrumQuasiLinear=+0.0d0
    end if
    ! Evaluate the halo power spectrum.
    if (self%includeHaloPower) then
       deltaHaloPrimedSquared  =+   self%a        *y **(3.0d0*self%f1)    & ! Smith et al. (2003; eqn. 50).
            &                   /(                                        &
            &                     +1.0d0                                  &
            &                     + self%b        *y **       self%f2     &
            &                     +(self%c*self%f3*y)**(3.0d0-self%gamma) &
            &                    )
       deltaHaloSquared        =+deltaHaloPrimedSquared & ! Smith et al. (2003; eqn. 51).
            &                   /(                      &
            &                     +1.0d0                &
            &                     +self%mu/y            &
            &                     +self%nu/y**2         &
            &                    )
       powerSpectrumHalo       =+deltaHaloSquared       &
            &                   *(2.0d0*Pi)**3          &
            &                   / 4.0d0/Pi              &
            &                   /wavenumber**3
       ! Include the Peacock correction if requested.
       if (self%includePeacockCorrection) then
          wavenumberPeacockCorrection=+10.0d0                                                       &
               &                      *self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
          powerSpectrumHalo          =+powerSpectrumHalo                                            &
               &                      *(1.0d0+2.0d0*(wavenumber/wavenumberPeacockCorrection)**2)    &
               &                      /(1.0d0+      (wavenumber/wavenumberPeacockCorrection)**2)
       end if
    else
       powerSpectrumHalo       =+0.0d0
    end if
    ! Sum the two contributions to the power spectrum.
    smith2003Value          =+powerSpectrumQuasiLinear &
         &                   +powerSpectrumHalo
    return
  end function smith2003Value
  
  subroutine smith2003Coefficients(self,time)
    !!{
    Evaluate the fitting function coefficients at the given time.
    !!}
    use :: Numerical_Comparison , only : Values_Agree
    use :: Numerical_Integration, only : integrator
    use :: Root_Finder          , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (powerSpectrumNonlinearSmith2003), intent(inout), target :: self
    double precision                                 , intent(in   )         :: time
    double precision                                 , parameter             :: wavenumberMultiplier=1.0d2
    type            (rootFinder                     )                        :: finder
    type            (integrator                     )                        :: integratorIndex           , integratorCurvature

    if (time /= self%timePrevious) then
       ! Record the time for which coefficients were computed.
       self%timePrevious=time
       ! Set a sub-module scope pointers.
       self_ => self
       time_ =  time
       ! Solve for the nonlinear wavenumber.
       finder=rootFinder()
       call finder%rootFunction(                                                             &
            &                                                 wavenumberNonLinearRoot        &
            &                  )
       call finder%rangeExpand (                                                             &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandType              =rangeExpandMultiplicative    , &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
            &                  )
       self%wavenumberNonlinear=finder%find(rootGuess=0.3d0)
       ! Compute effective spectral index (Smith et al. 2003; eqn. 59) and spectral curvature (Smith et al. 2003; eqn. 60).
       radiusGaussian=1.0d0/self%wavenumberNonlinear
       integratorIndex       = integrator                   (integrandEffectiveIndex   ,toleranceRelative=1.0d-6)
       integratorCurvature   = integrator                   (integrandSpectralCurvature,toleranceRelative=1.0d-6)
       self%effectiveIndex   =+integratorIndex    %integrate(0.0d0,self%wavenumberNonLinear*wavenumberMultiplier) &
            &                 *2.0d0                                                                              &
            &                 -3.0d0
       self%spectralCurvature=+integratorCurvature%integrate(0.0d0,self%wavenumberNonLinear*wavenumberMultiplier) &
            &                 *4.0d0                                                                              &
            &                 +(                                                                                  &
            &                   +3.0d0                                                                            &
            &                   +self%effectiveIndex                                                              &
            &                 )**2

       ! Compute cosmology-dependent parameters.
       if      (                                                                                  &
            &    Values_Agree(self%cosmologyParameters_%OmegaCurvature() ,   0.0d0,absTol=1.0d-3) &
            &  ) then
          ! For flat universe.
          self%f1=self%cosmologyFunctions_%omegaMatterEpochal(time=time)**(-0.0307d0)
          self%f2=self%cosmologyFunctions_%omegaMatterEpochal(time=time)**(-0.0585d0)
          self%f3=self%cosmologyFunctions_%omegaMatterEpochal(time=time)**(+0.0743d0)
       else if (                                                                                  &
            &    Values_Agree(self%cosmologyParameters_%OmegaDarkEnergy() ,  0.0d0,absTol=1.0d-3) &
            &   .and.                                                                             &
            &                 self%cosmologyParameters_%OmegaMatter    () <= 1.0d0                &
            &  ) then
          ! Zero cosmological constant.
          self%f1=self%cosmologyFunctions_%omegaMatterEpochal(time=time)**(-0.0732d0)
          self%f2=self%cosmologyFunctions_%omegaMatterEpochal(time=time)**(-0.1432d0)
          self%f3=self%cosmologyFunctions_%omegaMatterEpochal(time=time)**(+0.0725d0)
       else
          call Error_Report('no fitting functions available for this cosmology'//{introspection:location})
       end if
       ! Evaluate fitting function coefficients - using the fits from Smith et al. (2003; eqns C9 through C16).
       self%alpha=         +1.3884d0+0.3700d0*self%effectiveIndex-0.1452d0*self%effectiveIndex**2
       self%beta =         +0.8291d0+0.9854d0*self%effectiveIndex+0.3401d0*self%effectiveIndex**2
       self%a    =10.0d0**(+1.4861d0+1.8369d0*self%effectiveIndex+1.6762d0*self%effectiveIndex**2+0.7940d0*self%effectiveIndex**3+0.1670d0*self%effectiveIndex**4-0.6206d0*self%spectralCurvature)
       self%b    =10.0d0**(+0.9463d0+0.9466d0*self%effectiveIndex+0.3084d0*self%effectiveIndex**2                                                                -0.9400d0*self%spectralCurvature)
       self%c    =10.0d0**(-0.2807d0+0.6669d0*self%effectiveIndex+0.3214d0*self%effectiveIndex**2                                                                -0.0793d0*self%spectralCurvature)
       self%gamma=         +0.8649d0+0.2989d0*self%effectiveIndex+0.1631d0*self%spectralCurvature
       self%mu   =10.0d0**(-3.5442d0+0.1908d0*self%effectiveIndex                                                                                                                                )
       self%nu   =10.0d0**(+0.9589d0+1.2857d0*self%effectiveIndex                                                                                                                                )
    end if
    return
  end subroutine smith2003Coefficients

  double precision function wavenumberNonLinearRoot(wavenumberNonLinear)
    !!{
    Root function used to find the non-linear wavenumber, using the definition of \cite[][eqns.~53 \& 54]{smith_stable_2003}.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: wavenumberNonLinear
    double precision            , parameter     :: wavenumberMultiplier=1.0d2
    type            (integrator)                :: integrator_

    integrator_            = integrator           (integrandSigmaSquared,toleranceRelative=1.0d-6)
    radiusGaussian         =+1.0d0                                                                 &
         &                  /wavenumberNonLinear
    wavenumberNonLinearRoot=+integrator_%integrate(0.0d0,wavenumberNonLinear*wavenumberMultiplier) &
         &                  -1.0d0
    return
  end function wavenumberNonLinearRoot
  
  double precision function integrandSigmaSquared(wavenumber)
    !!{
    Integrand used to compute $\sigma^2(R)$ \citep[][eqn.~54]{smith_stable_2003}.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber

    if (wavenumber > 0.0d0) then
       integrandSigmaSquared=+self_%powerSpectrum_%powerDimensionless(wavenumber,time_) &
            &                *exp(                                                      &
            &                     -(                                                    &
            &                       +wavenumber                                         &
            &                       *radiusGaussian                                     &
            &                      )**2                                                 &
            &                    )                                                      &
            &                /wavenumber
    else
       integrandSigmaSquared=+0.0d0
    end if
    return
  end function integrandSigmaSquared

  double precision function integrandEffectiveIndex(wavenumber)
    !!{
    Integrand used to compute the effective spectral index\citep[][eqn.~59]{smith_stable_2003}.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber

    if (wavenumber > 0.0d0) then
       integrandEffectiveIndex=+self_%powerSpectrum_%powerDimensionless(wavenumber,time_) &
            &                  *     (                                                    &
            &                         +wavenumber                                         &
            &                         *radiusGaussian                                     &
            &                        )**2                                                 &
            &                  *exp(                                                      &
            &                       -(                                                    &
            &                         +wavenumber                                         &
            &                         *radiusGaussian                                     &
            &                        )**2                                                 &
            &                      )                                                      &
            &                  /wavenumber
    else
       integrandEffectiveIndex=+0.0d0
    end if
    return
  end function integrandEffectiveIndex

  double precision function integrandSpectralCurvature(wavenumber)
    !!{
    Integrand used to compute the spectral curvature \citep[][eqn.~60]{smith_stable_2003}.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber

    if (wavenumber > 0.0d0) then
       integrandSpectralCurvature=+self_%powerSpectrum_%powerDimensionless(wavenumber,time_) &
            &                     *(                                                         &
            &                       +(                                                       &
            &                         +wavenumber                                            &
            &                         *radiusGaussian                                        &
            &                        )**2                                                    &
            &                       -(                                                       &
            &                         +wavenumber                                            &
            &                         *radiusGaussian                                        &
            &                        )**4                                                    &
            &                      )                                                         &
            &                  *exp(                                                         &
            &                       -(                                                       &
            &                         +wavenumber                                            &
            &                         *radiusGaussian                                        &
            &                        )**2                                                    &
            &                      )                                                         &
            &                  /wavenumber
    else
       integrandSpectralCurvature=+0.0d0
    end if
    return
  end function integrandSpectralCurvature

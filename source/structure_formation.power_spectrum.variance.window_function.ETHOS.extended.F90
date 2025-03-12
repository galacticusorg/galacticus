!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implements a generalization of the ETHOS power spectrum window function class from \cite{bohr_halo_2021}.
  !!}
  
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredClass
  use :: Numerical_Interpolation             , only : interpolator

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionETHOSExtended">
   <description>
     A generalization of the ETHOS window function for filtering of power spectra from \cite{bohr_halo_2021}. The window function
     has the same functional form
     \begin{equation}
      W(kR) = \left\{ \begin{array}{ll} 1 &amp; \hbox{if } \frac{kR}{c_\mathrm{W}} &lt; x_\mathrm{min} \\ \frac{1}{1+ \left(\frac{kR}{c_\mathrm{W}} - x_\mathrm{min}\right)^\beta} &amp; \hbox{otherwise,} \end{array} \right.
     \end{equation}
     but the parameters $c_\mathrm{W}$ and $\beta$ are now scale dependent following     
     \begin{equation}
     x = x_0 x_1^{n-n_0}
     \end{equation}     
     where $x$ refers to either $c_\mathrm{W}$ or $\beta$, $n = \mathrm{d}\log P / \mathrm{d} \log k$ is the logarithmic
     derivative of the linear theory power spectrum, and $n_0 = -2.6$ is a convenient zero-point.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionETHOS) :: powerSpectrumWindowFunctionETHOSExtended
     !!{
     A generalization of the ETHOS power spectrum window function class.
     !!}
     private
     class           (powerSpectrumPrimordialTransferredClass), pointer                   :: powerSpectrumPrimordialTransferred_ => null()
     double precision                                                                     :: cW0                                          , beta0                      , &
          &                                                                                  cW1                                          , beta1                      , &
          &                                                                                  wavenumberScaledMinimum_                     , powerSpectrumSmoothingWidth, &
          &                                                                                  wavenumberMinimum_                           , wavenumberMaximum_         , &
          &                                                                                  timePrevious
     double precision                                         , allocatable, dimension(:) :: logPowerSpectra
     type            (interpolator                           ), allocatable               :: powerSpectrumSmoothed
   contains
     !![
     <methods>
       <method method="powerSpectrumSlopeSmoothed" description="Compute the slope of the smoothed power spectrum."/>
     </methods>
     !!]
     final     ::                               ETHOSExtendedDestructor
     procedure :: cW                         => ETHOSExtendedCW
     procedure :: beta                       => ETHOSExtendedBeta
     procedure :: wavenumberScaledMinimum    => ETHOSExtendedWavenumberScaledMinimum
     procedure :: powerSpectrumSlopeSmoothed => ETHOSExtendedPowerSpectrumSlopeSmoothed
  end type powerSpectrumWindowFunctionETHOSExtended

  interface powerSpectrumWindowFunctionETHOSExtended
     !!{
     Constructors for the ETHOS power spectrum window function class.
     !!}
     module procedure ETHOSExtendedConstructorParameters
     module procedure ETHOSExtendedConstructorInternal
  end interface powerSpectrumWindowFunctionETHOSExtended

  ! Zero-point in logarithmic slope of the power spectrum. Arbitrary, but convenient.
  double precision, parameter :: logarithmicDerivativeReference=-2.6d0

  ! Maximum allowed value for parameters.
  double precision, parameter :: parameterValueMaximum         =huge(0.0d0)/1.0d30
  double precision, parameter :: exponentPowerMaximum          =1.0d2

  ! Submodule-scope variables used in integration.
  class(powerSpectrumWindowFunctionETHOSExtended), pointer :: self_
  double precision :: time_, logWavenumbers_
  !$omp threadprivate(self_,time_,logWavenumbers_)
  
contains

  function ETHOSExtendedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ETHOS  power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionETHOSExtended)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (powerSpectrumPrimordialTransferredClass ), pointer       :: powerSpectrumPrimordialTransferred_
    class           (cosmologyParametersClass                ), pointer       :: cosmologyParameters_
    double precision                                                          :: cW0                                , beta0                      , &
         &                                                                       cW1                                , beta1                      , &
         &                                                                       wavenumberScaledMinimum            , powerSpectrumSmoothingWidth

    !![
    <inputParameter>
      <name>cW0</name>
      <source>parameters</source>
      <defaultValue>3.78062835d0</defaultValue>
      <description>The parameter $c_\mathrm{W,0}$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>beta0</name>
      <source>parameters</source>
      <defaultValue>3.4638743d0</defaultValue>
      <description>The parameter $\beta_0$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>cW1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $c_\mathrm{W,1}$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>beta1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $\beta_1$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>wavenumberScaledMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $x_\mathrm{min}$ in the generalized ETHOS power spectrum window function.</description>
    </inputParameter>
    <inputParameter>
      <name>powerSpectrumSmoothingWidth</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The width (in natural logarithm of wavenumber) over which to smooth the power spectrum when estimating the power spectrum slope.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"/>
    <objectBuilder class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferred_" source="parameters"/>
    !!]
    self=powerSpectrumWindowFunctionETHOSExtended(cW0,cW1,beta0,beta1,wavenumberScaledMinimum,powerSpectrumSmoothingWidth,cosmologyParameters_,powerSpectrumPrimordialTransferred_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"               />
    <objectDestructor name="powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end function ETHOSExtendedConstructorParameters

  function ETHOSExtendedConstructorInternal(cW0,cW1,beta0,beta1,wavenumberScaledMinimum_,powerSpectrumSmoothingWidth,cosmologyParameters_,powerSpectrumPrimordialTransferred_) result(self)
    !!{
    Internal constructor for the ETHOS power spectrum window function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (powerSpectrumWindowFunctionETHOSExtended)                        :: self
    double precision                                          , intent(in   )         :: cW0                                , beta0                      , &
         &                                                                               cW1                                , beta1                      , &
         &                                                                               wavenumberScaledMinimum_           , powerSpectrumSmoothingWidth
    class           (cosmologyParametersClass                ), intent(in   ), target :: cosmologyParameters_
    class           (powerSpectrumPrimordialTransferredClass ), intent(in   ), target :: powerSpectrumPrimordialTransferred_
    !![
    <constructorAssign variables="cW0, cW1, beta0, beta1, wavenumberScaledMinimum_, powerSpectrumSmoothingWidth, *cosmologyParameters_, *powerSpectrumPrimordialTransferred_"/>
    !!]

    self%timePrevious      =-huge(0.0d0)
    self%wavenumberMinimum_=+huge(0.0d0)
    self%wavenumberMaximum_=-huge(0.0d0)
    return
  end function ETHOSExtendedConstructorInternal

  subroutine ETHOSExtendedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily ETHOS} window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self

    !![
    <objectDestructor name="self%powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end subroutine ETHOSExtendedDestructor
  
  double precision function ETHOSExtendedCW(self,wavenumber,time) result(cW)
    !!{
    Compute the $c_\mathrm{W}$ parameter for the extended ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber   , time
    double precision                                                          :: exponentPower

    exponentPower=min(                                                                                      &
         &            max(                                                                                  &
         &                +self%powerSpectrumSlopeSmoothed(wavenumber,time)-logarithmicDerivativeReference, &
         &                -exponentPowerMaximum                                                             &
         &               )                                                                                , &
         &                +exponentPowerMaximum                                                             &
         &           )
    if (exponent(self%cW0)+exponent(self%cW1)*exponentPower < exponent(parameterValueMaximum)) then
       cW     =+self%cW0                &
            &  *self%cW1**exponentPower
    else
       cW     =+parameterValueMaximum
    end if
    return
  end function ETHOSExtendedCW
  
  double precision function ETHOSExtendedBeta(self,wavenumber,time) result(beta)
    !!{
    Compute the $\beta$ parameter for the extended ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber   , time
    double precision                                                          :: exponentPower

    exponentPower=min(                                                                                      &
         &            max(                                                                                  &
         &                +self%powerSpectrumSlopeSmoothed(wavenumber,time)-logarithmicDerivativeReference, &
         &                -exponentPowerMaximum                                                             &
         &               )                                                                                , &
         &                +exponentPowerMaximum                                                             &
         &           )
    if (exponent(self%beta0)+exponent(self%beta1)*exponentPower < maxExponent(parameterValueMaximum)) then
       beta   =+self%beta0                &
            &  *self%beta1**exponentPower
    else
       beta   =+parameterValueMaximum
    end if
    return
  end function ETHOSExtendedBeta
  
  double precision function ETHOSExtendedWavenumberScaledMinimum(self,wavenumber,time) result(wavenumberScaledMinimum)
    !!{
    Compute the $\beta$ parameter for the extended ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber, time
    !$GLC attributes unused :: wavenumber, time

    wavenumberScaledMinimum=self%wavenumberScaledMinimum_
    return
  end function ETHOSExtendedWavenumberScaledMinimum
  
  double precision function ETHOSExtendedPowerSpectrumSlopeSmoothed(self,wavenumber,time) result(slope)
    !!{
    Compute the logarithmic derivate of the power spectrum after smoothing.
    !!}
    use :: Numerical_Ranges     , only : Make_Range, rangeTypeLinear
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (powerSpectrumWindowFunctionETHOSExtended), intent(inout), target      :: self
    double precision                                          , intent(in   )              :: wavenumber                    , time
    double precision                                          , parameter                  :: countPerInterval       =10.0d0
    double precision                                          , dimension(:) , allocatable :: logWavenumbers                , logPowerSpectra
    integer                                                                                :: countPoints                   , i
    double precision                                                                       :: logWavenumberLimitLower       , logWavenumberLimitUpper
    logical                                                                                :: remakeTable

    if (.not.allocated(self%powerSpectrumSmoothed) .or. time /= self%timePrevious) then
       remakeTable=.true.
       if (allocated(self%powerSpectrumSmoothed)) deallocate(self%powerSpectrumSmoothed)
    else
       remakeTable= wavenumber < self%wavenumberMinimum_ &
            &      .or.                                  &
            &       wavenumber > self%wavenumberMaximum_
    end if
    if (remakeTable) then
       block
         type            (integrator) :: integrator_
         double precision             :: wavenumberMinimum_, wavenumberMaximum_
         integer                      :: countNewLower     , countNewUpper

         if (.not.allocated(self%powerSpectrumSmoothed)) then
            self%wavenumberMinimum_=min(self%wavenumberMinimum_,wavenumber/2.0d0)
            self%wavenumberMaximum_=max(self%wavenumberMaximum_,wavenumber*2.0d0)
            ! Allocate storage.
            countPoints  =int(log(self%wavenumberMaximum_/self%wavenumberMinimum_)/self%powerSpectrumSmoothingWidth*countPerInterval)+1
            countNewLower=0
            countNewUpper=0
            allocate(logPowerSpectra(countPoints))
            logPowerSpectra=-huge(0.0d0)
         else
            wavenumberMinimum_=min(self%wavenumberMinimum_,wavenumber/2.0d0)
            wavenumberMaximum_=max(self%wavenumberMaximum_,wavenumber*2.0d0)
            ! Determine how many points the table must be extended by in each direction to span the new required range.
            countNewLower=0
            countNewUpper=0
            if (self%wavenumberMinimum_ > wavenumberMinimum_) countNewLower=int(+log(self%wavenumberMinimum_/wavenumberMinimum_)/self%powerSpectrumSmoothingWidth*countPerInterval+1.0d0)
            if (self%wavenumberMaximum_ < wavenumberMaximum_) countNewUpper=int(-log(self%wavenumberMaximum_/wavenumberMaximum_)/self%powerSpectrumSmoothingWidth*countPerInterval+1.0d0)
            ! Adjust the limits of the table by an integer number of steps.
            self%wavenumberMinimum_=self%wavenumberMinimum_/exp(dble(countNewLower)/countPerInterval*self%powerSpectrumSmoothingWidth)
            self%wavenumberMaximum_=self%wavenumberMaximum_*exp(dble(countNewUpper)/countPerInterval*self%powerSpectrumSmoothingWidth)
            ! Allocate storage.
            countPoints=size(self%logPowerSpectra)+countNewLower+countNewUpper
            allocate(logPowerSpectra(countPoints))
            logPowerSpectra=-huge(0.0d0)
            ! Populate the table with pre-existing results.
            logPowerSpectra(countNewLower+1:countNewLower+size(self%logPowerSpectra))=self%logPowerSpectra
         end if
         allocate(logWavenumbers(countPoints))
         logWavenumbers =  Make_Range(log(self%wavenumberMinimum_),log(self%wavenumberMaximum_),countPoints,rangeTypeLinear)
         integrator_    =  integrator(integrandSmoothing,toleranceRelative=1.0d-3)
         self_          => self
         time_          =  time
         do i=1,countPoints
            if (logPowerSpectra(i) >= 0.0d0) cycle
            logWavenumbers_           =logWavenumbers(i)
            logWavenumberLimitLower   =logWavenumbers(i)-10.0d0*self%powerSpectrumSmoothingWidth
            logWavenumberLimitUpper   =logWavenumbers(i)+10.0d0*self%powerSpectrumSmoothingWidth
            logPowerSpectra        (i)=integrator_%integrate(logWavenumberLimitLower,logWavenumberLimitUpper)
            if (logPowerSpectra(i) > 0.0d0) then
               logPowerSpectra(i)=log(logPowerSpectra(i  ))
            else
               logPowerSpectra(i)=    logPowerSpectra(i-1)
            end if
         end do
         if (allocated(self%logPowerSpectra      )) deallocate(self%logPowerSpectra      )
         if (allocated(self%powerSpectrumSmoothed)) deallocate(self%powerSpectrumSmoothed)
         allocate(self%logPowerSpectra      (countPoints))
         allocate(self%powerSpectrumSmoothed             )
         self%logPowerSpectra      =                            logPowerSpectra
         self%powerSpectrumSmoothed=interpolator(logWavenumbers,logPowerSpectra)
         self%timePrevious         =time
       end block
    end if
    slope=self%powerSpectrumSmoothed%derivative(log(wavenumber))
    return
  end function ETHOSExtendedPowerSpectrumSlopeSmoothed
  
  double precision function integrandSmoothing(logWavenumber)
    !!{
    Integrand function used in smoothing the power spectrum.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: logWavenumber
    
    integrandSmoothing=+self_%powerSpectrumPrimordialTransferred_%power(exp(logWavenumber),time_) &
         &             *exp(                                                                      &
         &                  -0.5d0                                                                &
         &                  *(                                                                    &
         &                    +(                                                                  &
         &                      +logWavenumber                                                    &
         &                      -logWavenumbers_                                                  &
         &                     )                                                                  &
         &                    /self_%powerSpectrumSmoothingWidth                                  &
         &                   )**2                                                                 &
         &                 )                                                                      &
         &             /sqrt(2.0d0*Pi)                                                            &
         &             /       self_%powerSpectrumSmoothingWidth
    return
  end function integrandSmoothing

  

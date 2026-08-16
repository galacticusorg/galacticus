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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a generalization of the ETHOS power spectrum window function class from :cite:t:`bohr_halo_2021`.
  !!}
  
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredClass
  use :: Numerical_Interpolation             , only : interpolator

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionETHOSExtended" docformat="rst">
   <description>
   A generalization of the ETHOS window function for filtering of power spectra from :cite:t:`bohr_halo_2021`. The window function has the same functional form

   .. math::

      W(kR) = \left\{ \begin{array}{ll} 1 &amp; \hbox{if } \frac{kR}{c_\mathrm{W}} &lt; x_\mathrm{min} \\ \left[ 1+ \left(\frac{kR}{c_\mathrm{W}} - x_\mathrm{min}\right)^\beta \right]^{-1} &amp; \hbox{otherwise,} \end{array} \right.

   with :math:`x_\mathrm{min}=` ``[wavenumberScaledMinimum]``, but the parameters :math:`c_\mathrm{W}` and :math:`\beta` are now scale dependent following

   .. math::

      x = x_0 x_1^{n-n_0}

   where :math:`x` refers to either :math:`c_\mathrm{W}` or :math:`\beta` (controlled by parameters ``[cW0]``, ``[cW1]``, ``[beta0]``, and ``[beta1]``,, :math:`n = \mathrm{d}\log \tilde{P} / \mathrm{d} \log k` is the logarithmic derivative of the smoothed linear theory power spectrum, :math:`\tilde{P}`, and :math:`n_0 = -2.6` is a convenient zero-point. The smoothed power spectrum is defined as

   .. math::

      \tilde{P}(k) = \int_{-\infty}^\infty \frac{1}{\sqrt{2 \pi} \sigma} \exp\left[-\frac{1}{2} \left(\frac{\log k - \log k^\prime}{\sigma}\right)^2\right] P(k^\prime) \mathrm{d} \log k^\prime,

   with the :math:`\log k` terms using natural logarithms, and :math:`\sigma=` ``[powerSpectrumSmoothingWidth]``.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionETHOS) :: powerSpectrumWindowFunctionETHOSExtended
     !!{RST
     A generalization of the ETHOS power spectrum window function class.
     !!}
     private
     class           (powerSpectrumPrimordialTransferredClass), pointer                   :: powerSpectrumPrimordialTransferred_ => null()
     double precision                                                                     :: cW0                                          , beta0                      , &
          &                                                                                  cW1                                          , beta1                      , &
          &                                                                                  wavenumberScaledMinimum_                     , powerSpectrumSmoothingWidth, &
          &                                                                                  wavenumberMinimum_                           , wavenumberMaximum_         , &
          &                                                                                  timePrevious
     integer                                                                              :: indexMinimum_                                , indexMaximum_
     double precision                                         , allocatable, dimension(:) :: slopes
     type            (interpolator                           ), allocatable               :: slopeSmoothed
   contains
     !![
     <methods docformat="rst">
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
     !!{RST
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
    !!{RST
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
    <inputParameter docformat="rst">
      <name>cW0</name>
      <source>parameters</source>
      <defaultValue>3.78062835d0</defaultValue>
      <description>
      The parameter :math:`c_\mathrm{W,0}` in the generalized ETHOS power spectrum window function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>beta0</name>
      <source>parameters</source>
      <defaultValue>3.4638743d0</defaultValue>
      <description>
      The parameter :math:`\beta_0` in the generalized ETHOS power spectrum window function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>cW1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The parameter :math:`c_\mathrm{W,1}` in the generalized ETHOS power spectrum window function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>beta1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The parameter :math:`\beta_1` in the generalized ETHOS power spectrum window function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavenumberScaledMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The parameter :math:`x_\mathrm{min}` in the generalized ETHOS power spectrum window function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>powerSpectrumSmoothingWidth</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The width (in natural logarithm of wavenumber) over which to smooth the power spectrum when estimating the power spectrum slope.
      </description>
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
    !!{RST
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
    self%indexMinimum_     =+huge(0    )
    self%indexMaximum_     =-huge(0    )
    return
  end function ETHOSExtendedConstructorInternal

  subroutine ETHOSExtendedDestructor(self)
    !!{RST
    Destructor for the ``ETHOS`` window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self

    !![
    <objectDestructor name="self%powerSpectrumPrimordialTransferred_"/>
    !!]
    return
  end subroutine ETHOSExtendedDestructor
  
  double precision function ETHOSExtendedCW(self,wavenumber,time) result(cW)
    !!{RST
    Compute the :math:`c_\mathrm{W}` parameter for the extended ETHOS window function.
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
    !!{RST
    Compute the :math:`\beta` parameter for the extended ETHOS window function.
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
    !!{RST
    Compute the :math:`\beta` parameter for the extended ETHOS window function.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionETHOSExtended), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber, time
    !$GLC attributes unused :: wavenumber, time

    wavenumberScaledMinimum=self%wavenumberScaledMinimum_
    return
  end function ETHOSExtendedWavenumberScaledMinimum
  
  double precision function ETHOSExtendedPowerSpectrumSlopeSmoothed(self,wavenumber,time) result(slope)
    !!{RST
    Compute the logarithmic derivative of the power spectrum after smoothing.

    The derivative is evaluated analytically, by differentiating the smoothing kernel under the integral sign, so that

    .. math::

       \frac{\mathrm{d} \log \tilde{P}}{\mathrm{d} \log k} = \frac{1}{\tilde{P}(k)} \int_{-\infty}^\infty \left( \frac{\log k^\prime - \log k}{\sigma^2} \right) \frac{1}{\sqrt{2 \pi} \sigma} \exp\left[-\frac{1}{2} \left(\frac{\log k - \log k^\prime}{\sigma}\right)^2\right] P(k^\prime) \mathrm{d} \log k^\prime.

    The slope, rather than the smoothed power spectrum itself, is then tabulated and interpolated. Tabulating
    :math:`\log \tilde{P}` and differentiating the interpolant instead would amplify the numerical error in each tabulated
    point by the reciprocal of the tabulation interval, making the slope (and so the window function that depends on it)
    noisy on the scale of that interval.
    !!}
    use :: Error                , only : Error_Report
    use :: Numerical_Integration, only : integrator
    implicit none
    class           (powerSpectrumWindowFunctionETHOSExtended), intent(inout), target      :: self
    double precision                                          , intent(in   )              :: wavenumber                   , time
    double precision                                          , parameter                  :: countPerInterval       =1.0d2
    double precision                                          , parameter                  :: toleranceSlope         =1.0d-6
    double precision                                          , parameter                  :: widthsSmoothing        =1.0d1
    double precision                                          , dimension(:) , allocatable :: logWavenumbers               , slopes
    logical                                                   , dimension(:) , allocatable :: computed
    integer                                                                                :: countPoints                  , i
    double precision                                                                       :: logWavenumberLimitLower      , logWavenumberLimitUpper
    logical                                                                                :: remakeTable                  , retainTable

    if (.not.allocated(self%slopeSmoothed) .or. time /= self%timePrevious) then
       ! No table exists, or the table that does exist was built for a different time - in either case any tabulated slopes
       ! must be recomputed.
       remakeTable=.true.
       retainTable=.false.
       if (allocated(self%slopeSmoothed)) deallocate(self%slopeSmoothed)
    else
       ! A table for this time exists - it need be extended only if the requested wavenumber lies outside of it.
       remakeTable= wavenumber < self%wavenumberMinimum_ &
            &      .or.                                  &
            &       wavenumber > self%wavenumberMaximum_
       retainTable=.true.
    end if
    if (remakeTable) then
       block
         type            (integrator) :: integratorPower       , integratorPowerGradient
         double precision             :: separationLogarithmic , powerSmoothed          , &
              &                          powerSmoothedGradient
         integer                      :: indexMinimum          , indexMaximum           , &
              &                          offset

         ! Tabulate on an absolute lattice in log(wavenumber), of spacing set by the number of points per smoothing width, and
         ! anchored at a wavenumber of unity. Since the lattice does not depend on the range spanned by the table, points
         ! already computed remain at precisely the wavenumbers at which they were computed when the table is extended.
         separationLogarithmic=+self%powerSpectrumSmoothingWidth &
              &                /     countPerInterval
         indexMinimum         =min(self%indexMinimum_,int(floor  (log(wavenumber/2.0d0)/separationLogarithmic)))
         indexMaximum         =max(self%indexMaximum_,int(ceiling(log(wavenumber*2.0d0)/separationLogarithmic)))
         countPoints          =+indexMaximum &
              &                -indexMinimum &
              &                +1
         allocate(logWavenumbers(countPoints))
         allocate(slopes        (countPoints))
         allocate(computed      (countPoints))
         slopes  =0.0d0
         computed=.false.
         do i=1,countPoints
            logWavenumbers(i)=dble(indexMinimum+i-1)*separationLogarithmic
         end do
         ! Carry over any pre-existing results.
         if (retainTable.and.allocated(self%slopes)) then
            offset=+self%indexMinimum_ &
                 & -     indexMinimum
            slopes  (offset+1:offset+size(self%slopes))=self%slopes
            computed(offset+1:offset+size(self%slopes))=.true.
         end if
         ! Construct integrators for the smoothed power spectrum and its gradient. An absolute tolerance is required for the
         ! gradient integral as it passes through zero wherever the smoothed power spectrum peaks. It is set, for each point in
         ! turn, such that the resulting tolerance in the slope is `toleranceSlope`.
         integratorPower         =  integrator(integrandSmoothing        ,toleranceRelative=toleranceSlope)
         integratorPowerGradient =  integrator(integrandSmoothingGradient,toleranceRelative=toleranceSlope)
         self_                   => self
         time_                   =  time
         do i=1,countPoints
            if (computed(i)) cycle
            logWavenumbers_        =logWavenumbers(i)
            logWavenumberLimitLower=logWavenumbers(i)-widthsSmoothing*self%powerSpectrumSmoothingWidth
            logWavenumberLimitUpper=logWavenumbers(i)+widthsSmoothing*self%powerSpectrumSmoothingWidth
            powerSmoothed          =integratorPower%integrate(logWavenumberLimitLower,logWavenumberLimitUpper)
            ! If the smoothed power underflowed to zero the slope is indeterminate - such points are filled in below.
            if (powerSmoothed <= 0.0d0) cycle
            call integratorPowerGradient%toleranceSet(toleranceAbsolute=toleranceSlope*powerSmoothed,toleranceRelative=toleranceSlope)
            powerSmoothedGradient  =integratorPowerGradient%integrate(logWavenumberLimitLower,logWavenumberLimitUpper)
            slopes         (i)     =+powerSmoothedGradient &
                 &                  /powerSmoothed
            computed       (i)     =.true.
         end do
         ! Fill any points at which the smoothed power underflowed with the slope of the nearest point at which it did not.
         if (.not.any(computed)) call Error_Report('smoothed power spectrum underflowed at all tabulated wavenumbers'//{introspection:location})
         do i=2,countPoints
            if (.not.computed(i).and.computed(i-1)) then
               slopes  (i)=slopes(i-1)
               computed(i)=.true.
            end if
         end do
         do i=countPoints-1,1,-1
            if (.not.computed(i).and.computed(i+1)) then
               slopes  (i)=slopes(i+1)
               computed(i)=.true.
            end if
         end do
         if (allocated(self%slopes       )) deallocate(self%slopes       )
         if (allocated(self%slopeSmoothed)) deallocate(self%slopeSmoothed)
         allocate(self%slopes       (countPoints))
         allocate(self%slopeSmoothed             )
         self%slopes            =                            slopes
         self%slopeSmoothed     =interpolator(logWavenumbers,slopes)
         self%indexMinimum_     =indexMinimum
         self%indexMaximum_     =indexMaximum
         self%wavenumberMinimum_=exp(logWavenumbers(          1))
         self%wavenumberMaximum_=exp(logWavenumbers(countPoints))
         self%timePrevious      =time
       end block
    end if
    slope=self%slopeSmoothed%interpolate(log(wavenumber))
    return
  end function ETHOSExtendedPowerSpectrumSlopeSmoothed

  double precision function integrandSmoothing(logWavenumber)
    !!{RST
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

  double precision function integrandSmoothingGradient(logWavenumber)
    !!{RST
    Integrand function used in computing the gradient, with respect to :math:`\log k`, of the smoothed power spectrum. This is
    simply the integrand used in smoothing the power spectrum multiplied by the logarithmic derivative of the smoothing kernel.
    !!}
    implicit none
    double precision, intent(in   ) :: logWavenumber

    integrandSmoothingGradient=+(                                          &
         &                       +logWavenumber                            &
         &                       -logWavenumbers_                          &
         &                      )                                          &
         &                     /self_%powerSpectrumSmoothingWidth**2       &
         &                     *integrandSmoothing(logWavenumber)
    return
  end function integrandSmoothingGradient

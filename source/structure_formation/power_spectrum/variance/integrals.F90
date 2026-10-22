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

!+    Contributions to this file made by: Andrew Benson, Christoph Behrens, Xiaolong Du, Claude.

!!{RST
Contains a module which integrates over a filtered power spectrum to find the variance of the density field.
!!}

module Power_Spectrum_Variance_Integrals
  !!{RST
  Integrates over a filtered power spectrum to find the variance of the density field. The integration is shared by the
  :galacticus-class:`cosmologicalMassVarianceFilteredPower` class (to tabulate :math:`\sigma(M)`) and the
  :galacticus-class:`powerSpectrumNormalizerSigma8` class (to normalize the power spectrum to :math:`\sigma_8`), so that the two
  evaluate :math:`\sigma` identically.
  !!}
  private
  public :: Power_Spectrum_Variance_Integrate, Power_Spectrum_Variance_Wavenumber_Half_Mode, powerSpectrumVarianceIntegrand

  abstract interface
     double precision function powerSpectrumVarianceIntegrand(x)
       !!{RST
       Interface for integrands of the variance integral. For the linear integrand ``x`` is the wavenumber, :math:`k`, and the
       integrand is :math:`P(k) W^2(k) k^2`. For the logarithmic integrand ``x`` is :math:`\log k`, and the integrand is
       :math:`P(k) W^2(k) k^3`.
       !!}
       double precision, intent(in   ) :: x
     end function powerSpectrumVarianceIntegrand
  end interface

contains

  double precision function Power_Spectrum_Variance_Integrate(integrand,integrandLogarithmic,wavenumberMinimum,wavenumberMaximum,radiusTopHat,wavenumberHalfMode,wavenumbersLocalMinima,toleranceRelative,integrationFailureIsFatal,className,classValue,toleranceName) result(rootVariance)
    !!{RST
    Return the root-variance, :math:`\sigma = [(2\pi^2)^{-1} \int P(k) W^2(k) k^2 \mathrm{d}k]^{1/2}`, of a filtered power
    spectrum, integrating from ``wavenumberMinimum`` to ``wavenumberMaximum``.

    The integral is split into subsets of the wavenumber range:

    1. at a wavenumber corresponding to the smallest scale at which baryon acoustic oscillation features are significant
       (unless the upper limit of the integral is already below that wavenumber). This allows the oscillatory part of the
       integral to be computed more accurately, without affecting the non-oscillatory part at larger wavenumbers, and leads to
       an overall more accurate and robust determination of :math:`\sigma`;
    2. at a modest multiple of the inverse of ``radiusTopHat``. We would expect the window function to cut off above this
       scale, so it is useful to have a split just above this cut off;
    3. at a multiple of ``wavenumberHalfMode`` (if positive). For non-cold dark matter models with power suppressed on small
       scales this resolves :math:`\sigma` better near its saturation value, which is particularly useful if the transfer
       function oscillates on small scales;
    4. at each of ``wavenumbersLocalMinima`` (typically the local minima of the transfer function), which allows integration
       to accurately capture each oscillation of an oscillating transfer function.

    Where integration over :math:`k` fails numerically, or gives a zero result, over an interval with a non-zero lower limit
    the integral over that interval is instead evaluated over :math:`\log k`, which places more points at small wavenumber.

    The names ``className``, ``classValue``, and ``toleranceName`` identify the parameter controlling ``toleranceRelative``, and
    are used only in error messages.
    !!}
    use, intrinsic :: ISO_C_Binding           , only : c_size_t
    use            :: Display                 , only : displayIndent     , displayUnindent              , displayMessage, verbosityLevelStandard, &
         &                                             verbosityLevelWarn, enumerationVerbosityLevelType, displayReset  , displayGreen          , &
         &                                             displayBlue
    use            :: Error                   , only : Error_Report      , Warn
    use            :: Interface_GSL           , only : GSL_EBadTol       , GSL_ETol                     , GSL_ERound    , GSL_Success           , &
         &                                             GSL_EMaxIter      , gslErrorDecode
    use            :: ISO_Varying_String      , only : varying_string    , var_str                      , operator(//)  , assignment(=)
    use            :: Numerical_Constants_Math, only : Pi
    use            :: Numerical_Integration   , only : GSL_Integ_Gauss15 , integrator
    use            :: Sorting                 , only : sort
    use            :: String_Handling         , only : operator(//)      , stringXMLFormat
    implicit none
    procedure       (powerSpectrumVarianceIntegrand)                              :: integrand                               , integrandLogarithmic
    double precision                                , intent(in   )               :: wavenumberMinimum                       , wavenumberMaximum     , &
         &                                                                           radiusTopHat                            , wavenumberHalfMode    , &
         &                                                                           toleranceRelative
    double precision                                , intent(in   ), dimension(:) :: wavenumbersLocalMinima
    logical                                         , intent(in   )               :: integrationFailureIsFatal
    character       (len=*                         ), intent(in   )               :: className                               , classValue            , &
         &                                                                           toleranceName
    double precision                                , parameter                   :: wavenumberBAO            =5.0d0 ! The wavenumber above which baryon acoustic oscillations are small - used to split the integral allowing the oscillating part to be handled robustly.
    integer         (c_size_t                      ), parameter                   :: intervalsMaximum         =10000_c_size_t
    double precision                                , allocatable  , dimension(:) :: wavenumbers                             , wavenumbersRestricted
    double precision                                                              :: integral                                , integralInterval      , &
         &                                                                           wavenumberLower                         , wavenumberUpper
    integer                                                                       :: i                                       , status
    logical                                                                       :: computeLogarithmically
    type            (integrator                    )                              :: integrator_                             , integratorLogarithmic_

    integrator_           =integrator(integrand           ,toleranceRelative=toleranceRelative,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=intervalsMaximum)
    integratorLogarithmic_=integrator(integrandLogarithmic,toleranceRelative=toleranceRelative,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=intervalsMaximum)
    ! Assemble the points at which to split the integral.
    allocate(wavenumbers(4+size(wavenumbersLocalMinima)))
    wavenumbers(1:4                )=[                           &
         &                                   wavenumberMinimum , &
         &                            +3.0d0*wavenumberHalfMode, &
         &                            +3.0d0/radiusTopHat      , &
         &                                   wavenumberBAO       &
         &                           ]
    wavenumbers(5:size(wavenumbers))=wavenumbersLocalMinima
    ! Restrict the intervals to those between the minimum and maximum wavenumbers, and then sort them.
    wavenumbersRestricted=pack(wavenumbers,wavenumbers >= wavenumberMinimum .and. wavenumbers < wavenumberMaximum)
    call sort(wavenumbersRestricted)
    ! Iterate over intervals, accumulating the integral.
    integral=0.0d0
    do i=1,size(wavenumbersRestricted)
       wavenumberLower=wavenumbersRestricted(i)
       if (i < size(wavenumbersRestricted)) then
          wavenumberUpper=wavenumbersRestricted(i+1)
       else
          wavenumberUpper=wavenumberMaximum
       end if
       integralInterval=integrator_%integrate(wavenumberLower,wavenumberUpper,status=status)
       ! Decide if we need to try a logarithmic integral instead.
       if     (                         &
            &   (                       &
            &     status == GSL_EBadTol &
            &    .or.                   &
            &     status == GSL_ETol    &
            &    .or.                   &
            &     status == GSL_ERound  &
            &   )                       &
            &  .and.                    &
            &   wavenumberLower > 0.0d0 &
            & ) then
          ! The integration failed numerically, and the lower limit is non-zero, so try a logarithmic integral instead.
          computeLogarithmically=.true.
       else if (status /= GSL_Success) then
          ! Integration failed for some other reason, report an error.
          computeLogarithmically=.false.
          if (integrationFailureIsFatal) then
             block
               type     (varying_string) :: message
               character(len=12        ) :: label
               message=var_str('integration over interval failed [error number: ')//status//'; "'//gslErrorDecode(status)//'"]'//char(10)// &
                    &  displayGreen()//'HELP:'//displayReset()//' try increasing the value of <'//displayBlue()//toleranceName
               write (label,'(e12.6)') toleranceRelative
               message=message//stringXMLFormat('<'//className//' value="'//classValue//'">**B<'//toleranceName//' value="'//trim(adjustl(label))//'"/>**C</'//className//'>',indentInitial=3)//char(10)
               call Error_Report(message//{introspection:location})
             end block
          end if
       else if (integralInterval <= 0.0d0 .and. wavenumberLower > 0.0d0) then
          ! Integration gave a zero result, and the lower limit is non-zero. This may be because the upper limit is large and
          ! the power is confined to small wavenumbers near the lower limit. This can happen, for example, if attempting to
          ! compute σ(M) for mass scales far below the cut off for power spectra with a small-scale cut off. In such cases
          ! attempt to evaluate the integral again, but integrating over log(wavenumber) such that more points in the integrand
          ! are placed at small wavenumber.
          computeLogarithmically=.true.
       else
          ! The integration was successful, no need to try a logarithmic integral.
          computeLogarithmically=.false.
       end if
       ! Recompute the integral using a logarithmic integral if necessary.
       if (computeLogarithmically) then
          integralInterval=integratorLogarithmic_%integrate(log(wavenumberLower),log(wavenumberUpper),status=status)
          if ((status == GSL_EMaxIter .or. status == GSL_ERound) .and. integralInterval < toleranceRelative*integral) then
             ! Maximum iterations were exceeded, but the integrand is tiny - so proceed.
          else if (status /= GSL_Success) then
             if (integrationFailureIsFatal)                                                        &
                  & call Error_Report(var_str       ('variance integral failed [error number: ')// &
                  &                                  status                                     // &
                  &                           '; "'                                             // &
                  &                   gslErrorDecode(status                                    )// &
                  &                           '"]'                                              // &
                  &                   {introspection:location}                                     &
                  &                  )
          else if (integralInterval <= 0.0d0) then
             ! The integrated power is still non-positive - check the integrand at the lower limit. If it is positive, we have a
             ! problem. Otherwise, the power in this interval seems to be actually zero, so proceed.
             if (integrand(wavenumberLower) > 0.0d0) then
                if (integrand(wavenumberUpper) > 0.0d0) then
                   ! The integrand is non-zero at the upper limit also. A zero integral should therefore not be possible. First,
                   ! check if the integrand is so small that, when multiplied by the typical step size, it underflows to
                   ! zero. If it is, this is acceptable, and no error needs to be emitted. If such underflow did not occur, then
                   ! the integral should, definitely be positive. As it is not, this is an error situation.
                   if     (                                                   &
                        &    +max(                                            &
                        &         +integrandLogarithmic(log(wavenumberUpper)),&
                        &         +integrandLogarithmic(log(wavenumberLower)) &
                        &         )                                           &
                        &    *    (                                           &
                        &          +                    log(wavenumberUpper)  &
                        &          -                    log(wavenumberLower)  &
                        &         )                                           &
                        &    /10.0d0                                          &
                        &   >                                                 &
                        &    + 0.0d0                                          &
                        &  .and.                                              &
                        &   integrationFailureIsFatal                         &
                        & ) call Error_Report('no power in interval integrand - unexpected'//{introspection:location})
                else
                   ! The integrand is zero at the upper limit - try reducing the upper limit until we get a non-zero result.
                   do while (wavenumberUpper > wavenumberLower .and. integrand(wavenumberUpper) <= 0.0d0)
                      integralInterval=integratorLogarithmic_%integrate(log(wavenumberLower),log(wavenumberUpper),status=status)
                      if (status == GSL_Success .and. integralInterval > 0.0d0) exit
                      wavenumberUpper=sqrt(wavenumberUpper*wavenumberLower)
                   end do
                   ! If we still failed to get a non-zero power (or if the integral failed) - there's nothing more we can do.
                   if (integralInterval <= 0.0d0 .or. status /= GSL_Success) then
                      block
                        type     (varying_string               ) :: message
                        character(len=13                       ) :: label
                        integer                                  :: j
                        type     (enumerationVerbosityLevelType) :: verbosityLevel
                        if (integrationFailureIsFatal) then
                           verbosityLevel=verbosityLevelStandard
                        else
                           verbosityLevel=verbosityLevelWarn
                        end if
                        call displayIndent ('Integration report:'  ,verbosityLevel)
                        call displayIndent ('Wavenumber intervals:',verbosityLevel)
                        call displayMessage('i    kₘᵢₙ,ᵢ  kₘₐₓ,ᵢ'    ,verbosityLevel)
                        do j=1,size(wavenumbersRestricted)
                           wavenumberLower=wavenumbersRestricted(j)
                           if (j < size(wavenumbersRestricted)) then
                              wavenumberUpper=wavenumbersRestricted(j+1)
                           else
                              wavenumberUpper=wavenumberMaximum
                           end if
                           write (label,'(i4)') j
                           message=trim(label)
                           write (label,'(e13.6)') wavenumberLower
                           message=message//" "//label
                           write (label,'(e13.6)') wavenumberUpper
                           message=message//" "//label
                           if (j == i) message=message//" ← failed"
                           call displayMessage(message,verbosityLevel)
                        end do
                        call displayUnindent('done' ,verbosityLevel)
                        write (label,'(e13.6)') radiusTopHat
                        message="Top-hat radius = "//label
                        call displayMessage (message,verbosityLevel)
                        write (label,'(e13.6)') integral
                        message="Current integral = "//label
                        call displayMessage (message,verbosityLevel)
                        call displayUnindent('done' ,verbosityLevel)
                      end block
                      if (integrationFailureIsFatal) then
                         call Error_Report('variance integration failed - unexpected'//{introspection:location})
                      else
                         call Warn        ('variance integration failed - unexpected'                          )
                      end if
                   end if
                end if
             end if
          end if
       end if
       ! Accumulate the integral from this region.
       integral=+integral         &
            &   +integralInterval
    end do
    ! Apply factors of 2 and π to compute the root variance.
    rootVariance=+sqrt(           &
         &             +integral  &
         &             /2.0d0     &
         &             /Pi**2     &
         &            )
    return
  end function Power_Spectrum_Variance_Integrate

  double precision function Power_Spectrum_Variance_Wavenumber_Half_Mode(transferFunction_,cosmologyParameters_) result(wavenumberHalfMode)
    !!{RST
    Return the wavenumber corresponding to the half-mode mass of the given transfer function, :math:`k_\mathrm{hm} = \pi [4 \pi
    \Omega_\mathrm{m} \rho_\mathrm{crit} / 3 M_\mathrm{hm}]^{1/3}`, for use as a point at which to split the variance integral
    (see ``Power_Spectrum_Variance_Integrate``). A negative value is returned if the transfer function has no half-mode mass.
    !!}
    use :: Cosmology_Parameters    , only : cosmologyParametersClass
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    use :: Transfer_Functions      , only : transferFunctionClass
    implicit none
    class           (transferFunctionClass   ), intent(inout) :: transferFunction_
    class           (cosmologyParametersClass), intent(inout) :: cosmologyParameters_
    double precision                                          :: halfModeMass
    integer                                                   :: status

    halfModeMass=transferFunction_%halfModeMass(status)
    if (status == errorStatusSuccess) then
       wavenumberHalfMode=+(                                        &
            &               +4.0d0                                  &
            &               *Pi                                     &
            &               /3.0d0                                  &
            &               *cosmologyParameters_%OmegaMatter    () &
            &               *cosmologyParameters_%densityCritical() &
            &               /halfModeMass                           &
            &              )**(1.0d0/3.0d0)                         &
            &             *Pi
    else
       wavenumberHalfMode=-1.0d0
    end if
    return
  end function Power_Spectrum_Variance_Wavenumber_Half_Mode

end module Power_Spectrum_Variance_Integrals

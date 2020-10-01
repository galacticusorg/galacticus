!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements a two-point correlation function class in which the correlation function is found by Fourier
!% transforming a power spectrum.
  
  use :: Numerical_Interpolation , only : interpolator
  use :: Power_Spectra           , only : powerSpectrumClass
  use :: Power_Spectra_Nonlinear , only : powerSpectrumNonlinearClass

  !# <correlationFunctionTwoPoint name="correlationFunctionTwoPointPowerSpectrumTransform">
  !#  <description>Provides a two-point correlation function class in which the correlation function is found by Fourier transforming a power spectrum.</description>
  !# </correlationFunctionTwoPoint>
  type, extends(correlationFunctionTwoPointClass) :: correlationFunctionTwoPointPowerSpectrumTransform
     !% A two-point correlation function class in which the correlation function is found by Fourier transforming a power spectrum.
     private
     class           (powerSpectrumClass         ), pointer     :: powerSpectrum_                  => null()
     class           (powerSpectrumNonlinearClass), pointer     :: powerSpectrumNonlinear_         => null()
     type            (interpolator               ), allocatable :: interpolator_                            , interpolatorVolumeAveraged_ 
     double precision                                           :: separationMinimum                        , separationMaximum              , &
          &                                                        separationMinimumVolumeAveraged          , separationMaximumVolumeAveraged, &
          &                                                        time                                     , timeVolumeAveraged
  contains
     final     ::                              powerSpectrumTransformDestructor
     procedure :: correlation               => powerSpectrumTransformCorrelation
     procedure :: correlationVolumeAveraged => powerSpectrumTransformCorrelationVolumeAveraged
  end type correlationFunctionTwoPointPowerSpectrumTransform

  interface correlationFunctionTwoPointPowerSpectrumTransform
     !% Constructors for the {\normalfont \ttfamily powerSpectrumTransform} two-point correlation function class.
     module procedure powerSpectrumTransformConstructorParameters
     module procedure powerSpectrumTransformConstructorInternal
  end interface correlationFunctionTwoPointPowerSpectrumTransform

contains

  function powerSpectrumTransformConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily powerSpectrumTransform} two-point correlation function class which takes a
    !% parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (correlationFunctionTwoPointPowerSpectrumTransform)                :: self
    type (inputParameters                                  ), intent(inout) :: parameters
    class(powerSpectrumNonlinearClass                      ), pointer       :: powerSpectrumNonlinear_
    class(powerSpectrumClass                               ), pointer       :: powerSpectrum_

    powerSpectrumNonlinear_ => null()
    powerSpectrum_          => null()
    if (parameters%isPresent('powerSpectrumNonlinearMethod')) then
       !# <objectBuilder class="powerSpectrumNonlinear" name="powerSpectrumNonlinear_" source="parameters"/>
    end if
    if (parameters%isPresent('powerSpectrumMethod'         )) then
       !# <objectBuilder class="powerSpectrum"          name="powerSpectrum_"          source="parameters"/>
    end if
    !# <conditionalCall>
    !#  <call>self=correlationFunctionTwoPointPowerSpectrumTransform({conditions})</call>
    !#  <argument name="powerSpectrumNonlinear_" value="powerSpectrumNonlinear_" parameterPresent="parameters" parameterName="powerSpectrumNonlinearMethod"/>
    !#  <argument name="powerSpectrum_"          value="powerSpectrum_"          parameterPresent="parameters" parameterName="powerSpectrumMethod"         />
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="powerSpectrumNonlinear_"/>
    !# <objectDestructor name="powerSpectrum_"         />
    return
  end function powerSpectrumTransformConstructorParameters

  function powerSpectrumTransformConstructorInternal(powerSpectrumNonlinear_,powerSpectrum_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily powerSpectrumTransform} two-point correlation function class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type (correlationFunctionTwoPointPowerSpectrumTransform)                                  :: self
    class(powerSpectrumNonlinearClass                      ), intent(in   ), target, optional :: powerSpectrumNonlinear_
    class(powerSpectrumClass                               ), intent(in   ), target, optional :: powerSpectrum_
    !# <constructorAssign variables="*powerSpectrumNonlinear_, *powerSpectrum_"/>
    
    if      (      present(powerSpectrumNonlinear_).and.present(powerSpectrum_) ) then
       call Galacticus_Error_Report('provide either a linear or non-linear power spectrum, not both'//{introspection:location})
    else if (.not.(present(powerSpectrumNonlinear_).or. present(powerSpectrum_))) then
       call Galacticus_Error_Report('provide either a linear or non-linear power spectrum'          //{introspection:location})
    end if
    ! Initialize correlation range.
    self%time                           =-huge(0.0d0)
    self%timeVolumeAveraged             =-huge(0.0d0)
    self%separationMinimum              =+huge(0.0d0)
    self%separationMaximum              =-huge(0.0d0)
    self%separationMinimumVolumeAveraged=+huge(0.0d0)
    self%separationMaximumVolumeAveraged=-huge(0.0d0)
    return
  end function powerSpectrumTransformConstructorInternal

  subroutine powerSpectrumTransformDestructor(self)
    !% Destructor for the {\normalfont \ttfamily powerSpectrumTransform} two-point correlation function class.
    implicit none
    type(correlationFunctionTwoPointPowerSpectrumTransform), intent(inout) :: self

    !# <objectDestructor name="self%powerSpectrumNonlinear_"/>
    !# <objectDestructor name="self%powerSpectrum_"         />
    return
  end subroutine powerSpectrumTransformDestructor

  double precision function powerSpectrumTransformCorrelation(self,separation,time)
    !% Return a two-point correlation function by Fourier transforming a power spectrum.
    use :: FFTLogs                 , only : FFTLogSineTransform, fftLogForward 
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Make_Range         , rangeTypeLogarithmic
    implicit none
    class           (correlationFunctionTwoPointPowerSpectrumTransform), intent(inout)             :: self
    double precision                                                   , intent(in   )             :: time                      , separation
    double precision                                                   , allocatable, dimension(:) :: wavenumbers               , powerSpectrum    , &
         &                                                                                            correlation               , separations
    integer                                                            , parameter                 :: wavenumbersPerDecade=125
    double precision                                                   , parameter                 :: wavenumbersRange    =1.0d4
    double precision                                                                               :: wavenumberMinimum         , wavenumberMaximum
    integer         (c_size_t                                         )                            :: countWavenumbers          , i

    if (time /= self%time .or. separation < self%separationMinimum .or. separation > self%separationMaximum) then
       if (self%time < 0.0d0) then
          self%separationMinimum=separation
          self%separationMaximum=separation
       end if
       wavenumberMinimum=1.0d0/max(separation*wavenumbersRange,self%separationMaximum)
       wavenumberMaximum=1.0d0/min(separation/wavenumbersRange,self%separationMinimum)
       countWavenumbers =int(log10(wavenumberMaximum/wavenumberMinimum)*dble(wavenumbersPerDecade),c_size_t)
       allocate(wavenumbers  (countWavenumbers))
       allocate(powerSpectrum(countWavenumbers))
       allocate(correlation  (countWavenumbers))
       allocate(separations  (countWavenumbers))
       wavenumbers=Make_Range(wavenumberMinimum,wavenumberMaximum,int(countWavenumbers),rangeTypeLogarithmic)
       do i=1,countWavenumbers
          if (associated(self%powerSpectrum_)) then
             powerSpectrum(i)=self%powerSpectrum_         %power(wavenumbers(i),time)
          else
             powerSpectrum(i)=self%powerSpectrumNonLinear_%value(wavenumbers(i),time)
          end if
       end do
       call FFTLogSineTransform(                &
            &                    wavenumbers  , &
            &                    separations  , &
            &                   +powerSpectrum  &
            &                   *wavenumbers    &
            &                   * 4.0d0*Pi      &
            &                   /(2.0d0*Pi)**3, &
            &                    correlation  , &
            &                    fftLogForward  &
            &                  )
       correlation=correlation/separations
       if (allocated(self%interpolator_)) deallocate(self%interpolator_)
       allocate(self%interpolator_)
       self%interpolator_    =interpolator(separations,correlation)
       self%time             =time
       self%separationMinimum=separations(               1)
       self%separationMaximum=separations(countWavenumbers)
    end if
    powerSpectrumTransformCorrelation=self%interpolator_%interpolate(separation)
    return
  end function powerSpectrumTransformCorrelation
  
  double precision function powerSpectrumTransformCorrelationVolumeAveraged(self,separation,time)
    !% Return a volume-averaged two-point correlation function by Fourier transforming a power spectrum. The volume-averaged
    !% two-point correlation function is defined as:
    !% \begin{equation}
    !%  \bar{\xi}(r) = \int_0^2 \mathrm{d} r 4 \pi r^2 \xi(r) /  \int_0^2 \mathrm{d} r 4 \pi r^2.
    !% \end{equation}
    !% Since
    !% \begin{equation}
    !%  \xi(r) = \int \mathrm{d}k {P(k) \over (2 \pi)^3} 4 \pi {k^2 \over k r} \sin (k r),
    !% \end{equation}
    !% then
    !% \begin{equation}
    !%  \bar{\xi}(r) = 3 \int \mathrm{d}k {P(k) \over (2 \pi)^3} 4 \pi {k^2 \over (k r)^2} \left[ {\sin (k r) \over k r} - \cos(k r) \right].
    !% \end{equation}
    use :: FFTLogs                 , only : FFTLog     , fftLogForward
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Interpolation , only : interpolator
    use :: Numerical_Ranges        , only : Make_Range , rangeTypeLogarithmic
    implicit none
    class           (correlationFunctionTwoPointPowerSpectrumTransform), intent(inout)             :: self
    double precision                                                   , intent(in   )             :: time                      , separation
    double precision                                                   , allocatable, dimension(:) :: wavenumbers               , powerSpectrum    , &
         &                                                                                            correlation               , separations
    integer                                                            , parameter                 :: wavenumbersPerDecade=125
    double precision                                                   , parameter                 :: wavenumbersRange    =1.0d4
    double precision                                                                               :: wavenumberMinimum         , wavenumberMaximum
    integer         (c_size_t                                         )                            :: countWavenumbers          , i

    if (time /= self%timeVolumeAveraged .or. separation < self%separationMinimumVolumeAveraged .or. separation > self%separationMaximumVolumeAveraged) then
       if (self%timeVolumeAveraged < 0.0d0) then
          self%separationMinimumVolumeAveraged=separation
          self%separationMaximumVolumeAveraged=separation
       end if
       wavenumberMinimum=1.0d0/max(separation*wavenumbersRange,self%separationMaximumVolumeAveraged)
       wavenumberMaximum=1.0d0/min(separation/wavenumbersRange,self%separationMinimumVolumeAveraged)
       countWavenumbers =int(log10(wavenumberMaximum/wavenumberMinimum)*dble(wavenumbersPerDecade),c_size_t)
       allocate(wavenumbers  (countWavenumbers))
       allocate(powerSpectrum(countWavenumbers))
       allocate(correlation  (countWavenumbers))
       allocate(separations  (countWavenumbers))
       wavenumbers=Make_Range(wavenumberMinimum,wavenumberMaximum,int(countWavenumbers),rangeTypeLogarithmic)
       do i=1,countWavenumbers
          if (associated(self%powerSpectrum_)) then
             powerSpectrum(i)=self%powerSpectrum_         %power(wavenumbers(i),time)
          else
             powerSpectrum(i)=self%powerSpectrumNonLinear_%value(wavenumbers(i),time)
          end if
       end do
       ! In FFTLog(), we use μ=3/2 since:
       !   J_{3/2}(x) = √(2/π) [ sin(x)/x - cos(x) ] / √x
       ! So we must multiply by √(kr) √(π/2) so that we have:
       !   √(kr) √(π/2) J_{3/2}(x) = sin(x)/x - cos(x).    
       call FFTLog(                     &
            &       wavenumbers       , &
            &       separations       , &
            &      +powerSpectrum       &
            &      *sqrt(wavenumbers)   &
            &      * 3.0d0              &
            &      * 4.0d0*Pi           &
            &      /(2.0d0*Pi)**3     , &
            &       correlation       , &
            &        1.5d0            , &
            &       fftLogForward       &
            &     )
       correlation=+      correlation    &
            &      /      separations**2 &
            &      *sqrt(                &
            &            +separations    &
            &            *Pi             &
            &            /2.0d0          &
            &           )
       if (allocated(self%interpolatorVolumeAveraged_)) deallocate(self%interpolatorVolumeAveraged_)
       allocate(self%interpolatorVolumeAveraged_)
       self%interpolatorVolumeAveraged_    =interpolator(separations,correlation)
       self%timeVolumeAveraged             =time
       self%separationMinimumVolumeAveraged=separations(               1)
       self%separationMaximumVolumeAveraged=separations(countWavenumbers)
    end if
    powerSpectrumTransformCorrelationVolumeAveraged=self%interpolatorVolumeAveraged_%interpolate(separation)
    return
  end function powerSpectrumTransformCorrelationVolumeAveraged

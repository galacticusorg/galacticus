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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  A power spectrum normalizer which normalizes to a specified value of :math:`\sigma_8`.
  !!}

  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters           , only : cosmologyParametersClass
  use :: Power_Spectrum_Window_Functions, only : powerSpectrumWindowFunctionClass
  use :: Transfer_Functions             , only : transferFunctionClass

  !![
  <powerSpectrumNormalizer name="powerSpectrumNormalizerSigma8" docformat="rst">
   <description>
   A power spectrum normalizer which normalizes the power spectrum such that the linear theory root-variance of the density field
   at the present day in spheres of radius :math:`8h^{-1}` Mpc equals ``[sigma_8]``:

   .. math::

      \sigma_8^2 = {1 \over 2 \pi^2} \int_0^\infty P(k) W^2(k) k^2 \mathrm{d}k,

   where :math:`W(k)` is the Fourier transform of a real-space top-hat window function of radius :math:`8h^{-1}` Mpc (as per the
   definition of :math:`\sigma_8`), unless a different window function is explicitly defined via the
   ``[powerSpectrumWindowFunctionTopHat]`` parameter. The relative tolerance of the integration can be set via the
   ``[tolerance]`` parameter.
   </description>
  </powerSpectrumNormalizer>
  !!]
  type, extends(powerSpectrumNormalizerClass) :: powerSpectrumNormalizerSigma8
     !!{RST
     A power spectrum normalizer which normalizes to a specified value of :math:`\sigma_8`.
     !!}
     private
     class           (cosmologyParametersClass        ), pointer :: cosmologyParameters_         => null()
     class           (cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_          => null()
     class           (transferFunctionClass           ), pointer :: transferFunction_            => null()
     class           (powerSpectrumWindowFunctionClass), pointer :: powerSpectrumWindowFunction_ => null()
     double precision                                            :: sigma8                                , tolerance, &
          &                                                         wavenumberHalfMode
     logical                                                     :: integrationFailureIsFatal
   contains
     !![
     <methods docformat="rst">
       <method description="Add the parameters and objects common to this class and its extensions to a descriptor." method="descriptorCommon"/>
     </methods>
     !!]
     final     ::                     sigma8Destructor
     procedure :: normalization    => sigma8Normalization
     procedure :: descriptor       => sigma8Descriptor
     procedure :: descriptorCommon => sigma8DescriptorCommon
  end type powerSpectrumNormalizerSigma8

  interface powerSpectrumNormalizerSigma8
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumNormalizerSigma8` power spectrum normalizer class.
     !!}
     module procedure sigma8ConstructorParameters
     module procedure sigma8ConstructorInternal
  end interface powerSpectrumNormalizerSigma8

  ! Radius for σ₈ normalization in Mpc/h.
  double precision, parameter :: radiusNormalization=8.0d0

contains

  function sigma8ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumNormalizerSigma8` power spectrum normalizer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumNormalizerSigma8   )                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (transferFunctionClass           ), pointer       :: transferFunction_
    class           (powerSpectrumWindowFunctionClass), pointer       :: powerSpectrumWindowFunction_
    double precision                                                  :: sigma8                   , tolerance
    logical                                                           :: integrationFailureIsFatal

    !![
    <inputParameter docformat="rst">
      <name>sigma_8</name>
      <source>parameters</source>
      <variable>sigma8</variable>
      <defaultValue>0.8111d0</defaultValue>
      <defaultSource>
      (:cite:author:`planck_collaboration_planck_2018` :cite:year:`planck_collaboration_planck_2018`; TT,TE,EE\ :math:`+`\ lowE\ :math:`+`\ lensing)
      </defaultSource>
      <description>
      The fractional mass fluctuation in the linear density field at the present day in spheres of radius 8 Mpc/h.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>tolerance</name>
      <source>parameters</source>
      <defaultValue>1.0d-6</defaultValue>
      <description>
      The relative tolerance to use in integrating over the linear power spectrum to compute :math:`\sigma_8`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>integrationFailureIsFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true any failure in the integral to compute :math:`\sigma_8` is treated as a fatal error. Otherwise only a warning is issued.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunction_"    source="parameters"/>
    !!]
    if (parameters%isPresent('powerSpectrumWindowFunctionTopHat')) then
       !![
       <objectBuilder class="powerSpectrumWindowFunction" name="powerSpectrumWindowFunction_" source="parameters" parameterName="powerSpectrumWindowFunctionTopHat"/>
       !!]
       self=powerSpectrumNormalizerSigma8(sigma8,tolerance,integrationFailureIsFatal,cosmologyParameters_,cosmologyFunctions_,transferFunction_,powerSpectrumWindowFunction_)
       !![
       <objectDestructor name="powerSpectrumWindowFunction_"/>
       !!]
    else
       self=powerSpectrumNormalizerSigma8(sigma8,tolerance,integrationFailureIsFatal,cosmologyParameters_,cosmologyFunctions_,transferFunction_                             )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunction_"   />
    !!]
    return
  end function sigma8ConstructorParameters

  function sigma8ConstructorInternal(sigma8,tolerance,integrationFailureIsFatal,cosmologyParameters_,cosmologyFunctions_,transferFunction_,powerSpectrumWindowFunction_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumNormalizerSigma8` power spectrum normalizer class. If no
    ``powerSpectrumWindowFunction_`` is given a real-space top-hat window function is used.
    !!}
    use :: Power_Spectrum_Variance_Integrals, only : Power_Spectrum_Variance_Wavenumber_Half_Mode
    use :: Power_Spectrum_Window_Functions  , only : powerSpectrumWindowFunctionTopHat
    implicit none
    type            (powerSpectrumNormalizerSigma8   )                                  :: self
    double precision                                  , intent(in   )                   :: sigma8                      , tolerance
    logical                                           , intent(in   )                   :: integrationFailureIsFatal
    class           (cosmologyParametersClass        ), intent(in   ), target           :: cosmologyParameters_
    class           (cosmologyFunctionsClass         ), intent(in   ), target           :: cosmologyFunctions_
    class           (transferFunctionClass           ), intent(in   ), target           :: transferFunction_
    class           (powerSpectrumWindowFunctionClass), intent(in   ), target, optional :: powerSpectrumWindowFunction_
    !![
    <constructorAssign variables="sigma8, tolerance, integrationFailureIsFatal, *cosmologyParameters_, *cosmologyFunctions_, *transferFunction_, *powerSpectrumWindowFunction_"/>
    !!]

    if (.not.present(powerSpectrumWindowFunction_)) then
       allocate(powerSpectrumWindowFunctionTopHat :: self%powerSpectrumWindowFunction_)
       select type (powerSpectrumWindowFunction__ => self%powerSpectrumWindowFunction_)
       type is (powerSpectrumWindowFunctionTopHat)
          !![
          <referenceConstruct isResult="yes" object="powerSpectrumWindowFunction__" constructor="powerSpectrumWindowFunctionTopHat(cosmologyParameters_)"/>
          !!]
       end select
    end if
    self%wavenumberHalfMode=Power_Spectrum_Variance_Wavenumber_Half_Mode(self%transferFunction_,self%cosmologyParameters_)
    return
  end function sigma8ConstructorInternal

  subroutine sigma8Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`powerSpectrumNormalizerSigma8` power spectrum normalizer class.
    !!}
    implicit none
    type(powerSpectrumNormalizerSigma8), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"        />
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%transferFunction_"           />
    <objectDestructor name="self%powerSpectrumWindowFunction_"/>
    !!]
    return
  end subroutine sigma8Destructor

  double precision function sigma8Normalization(self,powerUnnormalized) result(normalization)
    !!{RST
    Return the factor by which the given power spectrum must be multiplied to normalize it to the specified :math:`\sigma_8`.
    !!}
    use :: Cosmology_Parameters             , only : hubbleUnitsLittleH
    use :: Numerical_Constants_Math         , only : Pi
    use :: Power_Spectrum_Variance_Integrals, only : Power_Spectrum_Variance_Integrate
    implicit none
    class           (powerSpectrumNormalizerSigma8           ), intent(inout), target       :: self
    procedure       (powerSpectrumNormalizerPowerUnnormalized)                              :: powerUnnormalized
    double precision                                          , allocatable  , dimension(:) :: wavenumbersLocalMinima
    double precision                                                                        :: massNormalization     , radiusNormalization_, &
         &                                                                                     timePresent           , wavenumberMaximum   , &
         &                                                                                     sigma8Unnormalized

    ! Find the mass, and corresponding top-hat radius, at which σ₈ is defined. The radius is computed from the mass (rather than
    ! taken directly as 8 Mpc/h) to match exactly the mass-to-radius conversion used in computing σ(M) (see
    ! cosmologicalMassVarianceFilteredPower).
    massNormalization   =+(                                                               &
         &                 +4.0d0                                                         &
         &                 /3.0d0                                                         &
         &                 *Pi                                                            &
         &                )                                                               &
         &               *  self%cosmologyParameters_%OmegaMatter    (                  ) &
         &               *  self%cosmologyParameters_%densityCritical(                  ) &
         &               *(                                                               &
         &                 +radiusNormalization                                           &
         &                 /self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH) &
         &                )**3
    radiusNormalization_=+(                                             &
         &                 +(                                           &
         &                   +3.0d0                                     &
         &                   /4.0d0                                     &
         &                   /Pi                                        &
         &                  )                                           &
         &                 *massNormalization                           &
         &                 /self%cosmologyParameters_%OmegaMatter    () &
         &                 /self%cosmologyParameters_%densityCritical() &
         &                )**(1.0d0/3.0d0)
    timePresent         =+self%cosmologyFunctions_%cosmicTime(1.0d0)
    wavenumberMaximum   =min(1.0d3/radiusNormalization_,self%powerSpectrumWindowFunction_%wavenumberMaximum(massNormalization))
    call self%transferFunction_%wavenumbersLocalMinima(wavenumbersLocalMinima)
    sigma8Unnormalized  =Power_Spectrum_Variance_Integrate(                                                    &
         &                                                 integrand                =integrand               , &
         &                                                 integrandLogarithmic     =integrandLogarithmic    , &
         &                                                 wavenumberMinimum        =0.0d0                   , &
         &                                                 wavenumberMaximum        =wavenumberMaximum       , &
         &                                                 radiusTopHat             =radiusNormalization_    , &
         &                                                 wavenumberHalfMode       =self%wavenumberHalfMode , &
         &                                                 wavenumbersLocalMinima   =wavenumbersLocalMinima  , &
         &                                                 toleranceRelative        =self%tolerance          , &
         &                                                 integrationFailureIsFatal=self%integrationFailureIsFatal, &
         &                                                 className                ='powerSpectrumNormalizer', &
         &                                                 classValue               ='sigma8'                , &
         &                                                 toleranceName            ='tolerance'               &
         &                                                )
    normalization       =+(                    &
         &                 +self%sigma8        &
         &                 /sigma8Unnormalized &
         &                )**2
    return

  contains

    double precision function integrand(wavenumber)
      !!{RST
      Integrand used in computing the variance in the window function from the unnormalized power spectrum. Factors of 2 and
      :math:`\pi` are included elsewhere.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      integrand=+  powerUnnormalized                       (wavenumber                  ,timePresent) &
           &    *(                                                                                    &
           &      +self%powerSpectrumWindowFunction_%value(wavenumber,massNormalization,timePresent) &
           &      *                                        wavenumber                                 &
           &     )**2
      return
    end function integrand

    double precision function integrandLogarithmic(wavenumberLogarithmic)
      !!{RST
      Integrand used in computing the variance in the window function from the unnormalized power spectrum. This version
      integrates with respect to :math:`\log(k)`. Factors of 2 and :math:`\pi` are included elsewhere.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumberLogarithmic
      double precision                :: wavenumber

      wavenumber          =+exp(wavenumberLogarithmic)
      integrandLogarithmic=+  powerUnnormalized                       (wavenumber                  ,timePresent) &
           &               *(                                                                                    &
           &                 +self%powerSpectrumWindowFunction_%value(wavenumber,massNormalization,timePresent) &
           &                 *                                        wavenumber                                 &
           &                )**2                                                                                 &
           &               *                                          wavenumber
      return
    end function integrandLogarithmic

  end function sigma8Normalization

  subroutine sigma8Descriptor(self,descriptor,includeClass,includeFileModificationTimes,parameterName)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters, descriptorParameterName
    implicit none
    class    (powerSpectrumNormalizerSigma8), intent(inout)           :: self
    type     (inputParameters              ), intent(inout)           :: descriptor
    logical                                 , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=*                        ), intent(in   ), optional :: parameterName
    character(len=18                       )                          :: parameterLabel
    type     (inputParameters              )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter(descriptorParameterName('powerSpectrumNormalizer',parameterName),'sigma8')
    parameters=descriptor%subparameters(descriptorParameterName('powerSpectrumNormalizer',parameterName))
    write (parameterLabel,'(e17.10)') self%sigma8
    call parameters%addParameter('sigma_8',trim(adjustl(parameterLabel)))
    call self%descriptorCommon(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine sigma8Descriptor

  subroutine sigma8DescriptorCommon(self,parameters,includeClass,includeFileModificationTimes)
    !!{RST
    Add the parameters and objects common to the :galacticus-class:`powerSpectrumNormalizerSigma8` class and its extensions to
    the given descriptor ``parameters``. The window function is added under its explicit parameter name,
    ``[powerSpectrumWindowFunctionTopHat]``, so that the descriptor recreates this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (powerSpectrumNormalizerSigma8), intent(inout)           :: self
    type     (inputParameters              ), intent(inout)           :: parameters
    logical                                 , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                       )                          :: parameterLabel

    write (parameterLabel,'(e17.10)') self%tolerance
    call parameters%addParameter('tolerance'                ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(l1)'    ) self%integrationFailureIsFatal
    call parameters%addParameter('integrationFailureIsFatal',trim(adjustl(parameterLabel)))
    call self%cosmologyParameters_        %descriptor(parameters,includeClass,includeFileModificationTimes                                                )
    call self%cosmologyFunctions_         %descriptor(parameters,includeClass,includeFileModificationTimes                                                )
    call self%transferFunction_           %descriptor(parameters,includeClass,includeFileModificationTimes                                                )
    call self%powerSpectrumWindowFunction_%descriptor(parameters,includeClass,includeFileModificationTimes,parameterName='powerSpectrumWindowFunctionTopHat')
    return
  end subroutine sigma8DescriptorCommon

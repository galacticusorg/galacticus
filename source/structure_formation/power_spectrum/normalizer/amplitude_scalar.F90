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
  A power spectrum normalizer which normalizes to a specified amplitude of the primordial scalar power spectrum.
  !!}

  !![
  <powerSpectrumNormalizer name="powerSpectrumNormalizerAmplitudeScalar" docformat="rst">
   <description>
   A power spectrum normalizer which normalizes the power spectrum to a specified amplitude, :math:`A_\mathrm{s}=`\
   ``[amplitudeScalar]``, of the primordial scalar power spectrum. The amplitude is converted to the corresponding value of
   :math:`\sigma_8` by running the CLASS Boltzmann code for the cosmological parameters of the model, and the power
   spectrum is then normalized to that :math:`\sigma_8` exactly as by the :galacticus-class:`powerSpectrumNormalizerSigma8`
   class (whose ``[tolerance]``, ``[integrationFailureIsFatal]``, and ``[powerSpectrumWindowFunctionTopHat]`` parameters apply
   here also).
   </description>
  </powerSpectrumNormalizer>
  !!]
  type, extends(powerSpectrumNormalizerSigma8) :: powerSpectrumNormalizerAmplitudeScalar
     !!{RST
     A power spectrum normalizer which normalizes to a specified amplitude of the primordial scalar power spectrum.
     !!}
     private
     double precision :: amplitudeScalar
   contains
     procedure :: descriptor => amplitudeScalarDescriptor
  end type powerSpectrumNormalizerAmplitudeScalar

  interface powerSpectrumNormalizerAmplitudeScalar
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumNormalizerAmplitudeScalar` power spectrum normalizer class.
     !!}
     module procedure amplitudeScalarConstructorParameters
     module procedure amplitudeScalarConstructorInternal
  end interface powerSpectrumNormalizerAmplitudeScalar

contains

  function amplitudeScalarConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumNormalizerAmplitudeScalar` power spectrum normalizer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumNormalizerAmplitudeScalar)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyParametersClass              ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class           (transferFunctionClass                 ), pointer       :: transferFunction_
    class           (powerSpectrumWindowFunctionClass      ), pointer       :: powerSpectrumWindowFunction_
    double precision                                                        :: amplitudeScalar          , tolerance
    logical                                                                 :: integrationFailureIsFatal

    !![
    <inputParameter docformat="rst">
      <name>amplitudeScalar</name>
      <source>parameters</source>
      <description>
      The amplitude of the primordial scalar power spectrum, :math:`A_\mathrm{s}`, such that :math:`P_\chi(k) = A_\mathrm{s} (k/k_\mathrm{s0})^{n_\mathrm{s}-1}` with :math:`k_\mathrm{s0}=0.05` Mpc\ :math:`^{-1}`.
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
       self=powerSpectrumNormalizerAmplitudeScalar(amplitudeScalar,tolerance,integrationFailureIsFatal,cosmologyParameters_,cosmologyFunctions_,transferFunction_,powerSpectrumWindowFunction_)
       !![
       <objectDestructor name="powerSpectrumWindowFunction_"/>
       !!]
    else
       self=powerSpectrumNormalizerAmplitudeScalar(amplitudeScalar,tolerance,integrationFailureIsFatal,cosmologyParameters_,cosmologyFunctions_,transferFunction_                             )
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunction_"   />
    !!]
    return
  end function amplitudeScalarConstructorParameters

  function amplitudeScalarConstructorInternal(amplitudeScalar,tolerance,integrationFailureIsFatal,cosmologyParameters_,cosmologyFunctions_,transferFunction_,powerSpectrumWindowFunction_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumNormalizerAmplitudeScalar` power spectrum normalizer class. If
    no ``powerSpectrumWindowFunction_`` is given a real-space top-hat window function is used.
    !!}
    use :: Interfaces_CLASS, only : Interface_CLASS_Normalization
    implicit none
    type            (powerSpectrumNormalizerAmplitudeScalar)                                  :: self
    double precision                                        , intent(in   )                   :: amplitudeScalar          , tolerance
    logical                                                 , intent(in   )                   :: integrationFailureIsFatal
    class           (cosmologyParametersClass              ), intent(in   ), target           :: cosmologyParameters_
    class           (cosmologyFunctionsClass               ), intent(in   ), target           :: cosmologyFunctions_
    class           (transferFunctionClass                 ), intent(in   ), target           :: transferFunction_
    class           (powerSpectrumWindowFunctionClass      ), intent(in   ), target, optional :: powerSpectrumWindowFunction_

    self%powerSpectrumNormalizerSigma8=powerSpectrumNormalizerSigma8(-1.0d0,tolerance,integrationFailureIsFatal,cosmologyParameters_,cosmologyFunctions_,transferFunction_,powerSpectrumWindowFunction_)
    self%amplitudeScalar              =amplitudeScalar
    ! Convert the scalar amplitude to σ₈. Interface_CLASS_Normalization() returns σ₈²/Aₛ.
    self%sigma8                       =sqrt(amplitudeScalar*Interface_CLASS_Normalization(self%cosmologyParameters_))
    return
  end function amplitudeScalarConstructorInternal

  subroutine amplitudeScalarDescriptor(self,descriptor,includeClass,includeFileModificationTimes,parameterName)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters, descriptorParameterName
    implicit none
    class    (powerSpectrumNormalizerAmplitudeScalar), intent(inout)           :: self
    type     (inputParameters                       ), intent(inout)           :: descriptor
    logical                                          , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=*                                 ), intent(in   ), optional :: parameterName
    character(len=18                                )                          :: parameterLabel
    type     (inputParameters                       )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter(descriptorParameterName('powerSpectrumNormalizer',parameterName),'amplitudeScalar')
    parameters=descriptor%subparameters(descriptorParameterName('powerSpectrumNormalizer',parameterName))
    write (parameterLabel,'(e17.10)') self%amplitudeScalar
    call parameters%addParameter('amplitudeScalar',trim(adjustl(parameterLabel)))
    call self%descriptorCommon(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine amplitudeScalarDescriptor

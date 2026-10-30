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
  A power spectrum normalizer which matches the power spectrum to a reference power spectrum at a specified wavenumber.
  !!}

  use :: Cosmology_Functions                 , only : cosmologyFunctionsClass
  use :: Power_Spectra_Primordial_Transferred, only : powerSpectrumPrimordialTransferredClass

  !![
  <powerSpectrumNormalizer name="powerSpectrumNormalizerMatchedWavenumber" docformat="rst">
   <description>
   A power spectrum normalizer which normalizes the power spectrum such that, at the present day, it equals a reference power
   spectrum at the wavenumber ``[wavenumberReference]``. The reference power spectrum is the
   :galacticus-class:`powerSpectrumPrimordialTransferred` defined in the ``[reference]`` parameter section, which must also
   define the normalizer of the reference power spectrum (a :galacticus-class:`powerSpectrumNormalizer`, either directly within
   the ``[reference]`` section or within its :galacticus-class:`powerSpectrumPrimordialTransferred`). For example:

   .. code-block:: xml

      &lt;powerSpectrumNormalizer value="matchedWavenumber"&gt;
        &lt;wavenumberReference value="100.0"/&gt;
        &lt;reference&gt;
          &lt;powerSpectrumPrimordialTransferred value="simple"/&gt;
          &lt;powerSpectrumNormalizer value="sigma8"&gt;
            &lt;sigma_8 value="0.8228"/&gt;
          &lt;/powerSpectrumNormalizer&gt;
        &lt;/reference&gt;
      &lt;/powerSpectrumNormalizer&gt;

   This is useful, for example, where the power spectrum of a model excludes modes (such as those which do not fit within a
   simulation volume) such that normalizing it to :math:`\sigma_8` would result in too much power on small scales.
   </description>
  </powerSpectrumNormalizer>
  !!]
  type, extends(powerSpectrumNormalizerClass) :: powerSpectrumNormalizerMatchedWavenumber
     !!{RST
     A power spectrum normalizer which matches the power spectrum to a reference power spectrum at a specified wavenumber.
     !!}
     private
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_                          => null()
     class           (powerSpectrumPrimordialTransferredClass), pointer :: powerSpectrumPrimordialTransferredReference_ => null()
     double precision                                                   :: wavenumberReference
   contains
     final     ::                  matchedWavenumberDestructor
     procedure :: normalization => matchedWavenumberNormalization
     procedure :: descriptor    => matchedWavenumberDescriptor
  end type powerSpectrumNormalizerMatchedWavenumber

  interface powerSpectrumNormalizerMatchedWavenumber
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumNormalizerMatchedWavenumber` power spectrum normalizer class.
     !!}
     module procedure matchedWavenumberConstructorParameters
     module procedure matchedWavenumberConstructorInternal
  end interface powerSpectrumNormalizerMatchedWavenumber

contains

  function matchedWavenumberConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumNormalizerMatchedWavenumber` power spectrum normalizer class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumNormalizerMatchedWavenumber)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    type            (inputParameters                         )                :: referenceParameters                         , referenceTransferredParameters
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (powerSpectrumPrimordialTransferredClass ), pointer       :: powerSpectrumPrimordialTransferredReference_
    double precision                                                          :: wavenumberReference

    if (.not.parameters%isPresent('reference',requireValue=.false.)) call Error_Report('parameters must contain a "reference" section'//{introspection:location})
    referenceParameters=parameters%subParameters('reference',requireValue=.false.)
    if (.not.referenceParameters%isPresent('powerSpectrumPrimordialTransferred')) call Error_Report('"reference" section must explicitly define a "powerSpectrumPrimordialTransferred"'//{introspection:location})
    ! The reference power spectrum must have its own normalizer. Otherwise, the search for it would proceed up the parameter tree
    ! and find this normalizer, which is currently being built.
    referenceTransferredParameters=referenceParameters%subParameters('powerSpectrumPrimordialTransferred')
    if     (                                                                         &
         &  .not.referenceParameters           %isPresent('powerSpectrumNormalizer') &
         &  .and.                                                                    &
         &  .not.referenceTransferredParameters%isPresent('powerSpectrumNormalizer') &
         & ) call Error_Report('"reference" section must explicitly define a "powerSpectrumNormalizer" for the reference power spectrum'//{introspection:location})
    !![
    <inputParameter docformat="rst">
      <name>wavenumberReference</name>
      <source>parameters</source>
      <description>
      The reference wavenumber (in units of Mpc\ :math:`^{-1}`) at which the amplitude of the power spectrum is matched to that of the reference power spectrum.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                          source="parameters"         />
    <objectBuilder class="powerSpectrumPrimordialTransferred" name="powerSpectrumPrimordialTransferredReference_" source="referenceParameters"/>
    !!]
    self=powerSpectrumNormalizerMatchedWavenumber(wavenumberReference,cosmologyFunctions_,powerSpectrumPrimordialTransferredReference_)
    !![
    <inputParametersValidate source="parameters" extraAllowedNames="reference"/>
    <objectDestructor name="cosmologyFunctions_"                         />
    <objectDestructor name="powerSpectrumPrimordialTransferredReference_"/>
    !!]
    return
  end function matchedWavenumberConstructorParameters

  function matchedWavenumberConstructorInternal(wavenumberReference,cosmologyFunctions_,powerSpectrumPrimordialTransferredReference_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumNormalizerMatchedWavenumber` power spectrum normalizer class.
    !!}
    implicit none
    type            (powerSpectrumNormalizerMatchedWavenumber)                        :: self
    double precision                                          , intent(in   )         :: wavenumberReference
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (powerSpectrumPrimordialTransferredClass ), intent(in   ), target :: powerSpectrumPrimordialTransferredReference_
    !![
    <constructorAssign variables="wavenumberReference, *cosmologyFunctions_, *powerSpectrumPrimordialTransferredReference_"/>
    !!]

    return
  end function matchedWavenumberConstructorInternal

  subroutine matchedWavenumberDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`powerSpectrumNormalizerMatchedWavenumber` power spectrum normalizer class.
    !!}
    implicit none
    type(powerSpectrumNormalizerMatchedWavenumber), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"                         />
    <objectDestructor name="self%powerSpectrumPrimordialTransferredReference_"/>
    !!]
    return
  end subroutine matchedWavenumberDestructor

  double precision function matchedWavenumberNormalization(self,powerUnnormalized) result(normalization)
    !!{RST
    Return the factor by which the given power spectrum must be multiplied to match the reference power spectrum at the reference
    wavenumber at the present day.
    !!}
    implicit none
    class           (powerSpectrumNormalizerMatchedWavenumber), intent(inout), target :: self
    procedure       (powerSpectrumNormalizerPowerUnnormalized)                        :: powerUnnormalized
    double precision                                                                  :: timePresent

    timePresent  =+self%cosmologyFunctions_                         %cosmicTime       (1.0d0                               )
    normalization=+self%powerSpectrumPrimordialTransferredReference_%power            (self%wavenumberReference,timePresent) &
         &        /                                                  powerUnnormalized(self%wavenumberReference,timePresent)
    return
  end function matchedWavenumberNormalization

  subroutine matchedWavenumberDescriptor(self,descriptor,includeClass,includeFileModificationTimes,parameterName)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters, descriptorParameterName
    implicit none
    class    (powerSpectrumNormalizerMatchedWavenumber), intent(inout)           :: self
    type     (inputParameters                         ), intent(inout)           :: descriptor
    logical                                            , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=*                                   ), intent(in   ), optional :: parameterName
    character(len=18                                  )                          :: parameterLabel
    type     (inputParameters                         )                          :: parameters    , referenceParameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter(descriptorParameterName('powerSpectrumNormalizer',parameterName),'matchedWavenumber')
    parameters=descriptor%subparameters(descriptorParameterName('powerSpectrumNormalizer',parameterName))
    write (parameterLabel,'(e17.10)') self%wavenumberReference
    call parameters%addParameter('wavenumberReference',trim(adjustl(parameterLabel)))
    call self%cosmologyFunctions_%descriptor(parameters,includeClass,includeFileModificationTimes)
    call parameters%addParameter('reference','')
    referenceParameters=parameters%subparameters('reference')
    call self%powerSpectrumPrimordialTransferredReference_%descriptor(referenceParameters,includeClass,includeFileModificationTimes)
    return
  end subroutine matchedWavenumberDescriptor

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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module of global functions which operate on ``powerSpectrumNormalizer`` objects.
!!}

module Power_Spectrum_Normalizer_Utilities
  !!{RST
  Provides global functions which operate on :galacticus-class:`powerSpectrumNormalizer` objects.

  The :galacticus-class:`powerSpectrumPrimordialTransferred` classes hold their normalizer as an unlimited polymorphic
  (``class(*)``) object, and operate on it only through these functions (accessed via ``Functions_Global``), other than for
  reference counting (which is done through the ``functionClass`` base type, so needs no knowledge of the normalizer class). This avoids a
  circular module dependency: some normalizers depend on modules (such as ``Power_Spectrum_Window_Functions``, or
  ``Power_Spectra_Primordial_Transferred`` itself) which in turn depend on ``Power_Spectra_Primordial_Transferred``.
  !!}
  use :: Power_Spectrum_Normalizers, only : powerSpectrumNormalizer, powerSpectrumNormalizerClass
  private
  public :: Power_Spectrum_Normalizer_Construct         , Power_Spectrum_Normalizer_Release  , Power_Spectrum_Normalizer_Normalization  , &
       &    Power_Spectrum_Normalizer_Descriptor        , Power_Spectrum_Normalizer_Deep_Copy, Power_Spectrum_Normalizer_Deep_Copy_Reset, &
       &    Power_Spectrum_Normalizer_Deep_Copy_Finalize, Power_Spectrum_Normalizer_Allowed_Parameters

contains

  !![
  <functionGlobal>
   <unitName>Power_Spectrum_Normalizer_Construct</unitName>
   <type>class(*), pointer</type>
   <module>Input_Parameters, only : inputParameters</module>
   <arguments>type(inputParameters), intent(inout), target :: parameters</arguments>
  </functionGlobal>
  !!]
  function Power_Spectrum_Normalizer_Construct(parameters) result(self)
    !!{RST
    Build a :galacticus-class:`powerSpectrumNormalizer` object from the given ``parameters``, returning it as an unlimited
    polymorphic pointer. As for any object built from parameters, the normalizer is found in ``parameters``, or in any of its
    parents, and a default normalizer is built if none is found.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class(*                           ), pointer               :: self
    type (inputParameters             ), intent(inout), target :: parameters
    class(powerSpectrumNormalizerClass), pointer               :: powerSpectrumNormalizer_

    !![
    <objectBuilder class="powerSpectrumNormalizer" name="powerSpectrumNormalizer_" source="parameters"/>
    !!]
    self => powerSpectrumNormalizer_
    nullify(powerSpectrumNormalizer_)
    return
  end function Power_Spectrum_Normalizer_Construct

  !![
  <functionGlobal>
   <unitName>Power_Spectrum_Normalizer_Release</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout), pointer :: self</arguments>
  </functionGlobal>
  !!]
  subroutine Power_Spectrum_Normalizer_Release(self)
    !!{RST
    Release a reference to a :galacticus-class:`powerSpectrumNormalizer` object (i.e. decrement its reference count, destroying
    it if no references remain), and nullify the pointer.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*                           ), intent(inout), pointer :: self
    class(powerSpectrumNormalizerClass)               , pointer :: powerSpectrumNormalizer_

    if (.not.associated(self)) return
    select type (self)
    class is (powerSpectrumNormalizerClass)
       powerSpectrumNormalizer_ => self
    class default
       call Error_Report('object must be of `powerSpectrumNormalizerClass` class'//{introspection:location})
    end select
    !![
    <objectDestructor name="powerSpectrumNormalizer_"/>
    !!]
    nullify(self)
    return
  end subroutine Power_Spectrum_Normalizer_Release

  !![
  <functionGlobal>
   <unitName>Power_Spectrum_Normalizer_Normalization</unitName>
   <type>double precision</type>
   <module>Power_Spectrum_Normalizer_Interfaces, only : powerSpectrumNormalizerPowerUnnormalized</module>
   <arguments>class    (*                                       ), intent(inout) :: self</arguments>
   <arguments>procedure(powerSpectrumNormalizerPowerUnnormalized)                :: powerUnnormalized</arguments>
  </functionGlobal>
  !!]
  double precision function Power_Spectrum_Normalizer_Normalization(self,powerUnnormalized) result(normalization)
    !!{RST
    Return the factor by which the power spectrum given by ``powerUnnormalized`` must be multiplied to normalize it, as computed
    by the given :galacticus-class:`powerSpectrumNormalizer` object.
    !!}
    use :: Error                               , only : Error_Report
    use :: Power_Spectrum_Normalizer_Interfaces, only : powerSpectrumNormalizerPowerUnnormalized
    implicit none
    class    (*                                       ), intent(inout) :: self
    procedure(powerSpectrumNormalizerPowerUnnormalized)                :: powerUnnormalized

    select type (self)
    class is (powerSpectrumNormalizerClass)
       normalization=self%normalization(powerUnnormalized)
    class default
       normalization=0.0d0
       call Error_Report('object must be of `powerSpectrumNormalizerClass` class'//{introspection:location})
    end select
    return
  end function Power_Spectrum_Normalizer_Normalization

  !![
  <functionGlobal>
   <unitName>Power_Spectrum_Normalizer_Descriptor</unitName>
   <type>void</type>
   <module>Input_Parameters, only : inputParameters</module>
   <arguments>class  (*              ), intent(inout)           :: self                                      </arguments>
   <arguments>type   (inputParameters), intent(inout)           :: descriptor                                </arguments>
   <arguments>logical                 , intent(in   ), optional :: includeClass, includeFileModificationTimes</arguments>
  </functionGlobal>
  !!]
  subroutine Power_Spectrum_Normalizer_Descriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{RST
    Add the descriptor of the given :galacticus-class:`powerSpectrumNormalizer` object to ``descriptor``.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    class  (*              ), intent(inout)           :: self
    type   (inputParameters), intent(inout)           :: descriptor
    logical                 , intent(in   ), optional :: includeClass, includeFileModificationTimes

    select type (self)
    class is (powerSpectrumNormalizerClass)
       call self%descriptor(descriptor,includeClass,includeFileModificationTimes)
    class default
       call Error_Report('object must be of `powerSpectrumNormalizerClass` class'//{introspection:location})
    end select
    return
  end subroutine Power_Spectrum_Normalizer_Descriptor

  !![
  <functionGlobal>
   <unitName>Power_Spectrum_Normalizer_Deep_Copy</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout), target  :: self       </arguments>
   <arguments>class(*), intent(  out), pointer :: destination</arguments>
  </functionGlobal>
  !!]
  subroutine Power_Spectrum_Normalizer_Deep_Copy(self,destination)
    !!{RST
    Return in ``destination`` a deep copy of the given :galacticus-class:`powerSpectrumNormalizer` object. If the object has
    already been copied during the current deep copy operation, ``destination`` points to (and acquires a reference to) that copy,
    so that objects shared in the original remain shared in the copy.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*                           ), intent(inout), target  :: self
    class(*                           ), intent(  out), pointer :: destination
    class(powerSpectrumNormalizerClass)               , pointer :: copy

    select type (self)
    class is (powerSpectrumNormalizerClass)
       if (associated(self%copiedSelf)) then
          select type (copiedSelf => self%copiedSelf)
          class is (powerSpectrumNormalizerClass)
             copy => copiedSelf
          class default
             copy => null()
             call Error_Report('copiedSelf has incorrect type'//{introspection:location})
          end select
          call copy%referenceCountIncrement()
       else
          allocate(copy,mold=self)
          call self%deepCopy(copy)
          self%copiedSelf => copy
          call copy%autoHook()
       end if
       destination => copy
    class default
       destination => null()
       call Error_Report('object must be of `powerSpectrumNormalizerClass` class'//{introspection:location})
    end select
    return
  end subroutine Power_Spectrum_Normalizer_Deep_Copy

  !![
  <functionGlobal>
   <unitName>Power_Spectrum_Normalizer_Deep_Copy_Reset</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine Power_Spectrum_Normalizer_Deep_Copy_Reset(self)
    !!{RST
    Reset the deep copy state of the given :galacticus-class:`powerSpectrumNormalizer` object.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*), intent(inout) :: self

    select type (self)
    class is (powerSpectrumNormalizerClass)
       call self%deepCopyReset()
    class default
       call Error_Report('object must be of `powerSpectrumNormalizerClass` class'//{introspection:location})
    end select
    return
  end subroutine Power_Spectrum_Normalizer_Deep_Copy_Reset

  !![
  <functionGlobal>
   <unitName>Power_Spectrum_Normalizer_Deep_Copy_Finalize</unitName>
   <type>void</type>
   <arguments>class(*), intent(inout) :: self</arguments>
  </functionGlobal>
  !!]
  subroutine Power_Spectrum_Normalizer_Deep_Copy_Finalize(self)
    !!{RST
    Finalize the deep copy of the given :galacticus-class:`powerSpectrumNormalizer` object.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*), intent(inout) :: self

    select type (self)
    class is (powerSpectrumNormalizerClass)
       call self%deepCopyFinalize()
    class default
       call Error_Report('object must be of `powerSpectrumNormalizerClass` class'//{introspection:location})
    end select
    return
  end subroutine Power_Spectrum_Normalizer_Deep_Copy_Finalize

  !![
  <functionGlobal>
   <unitName>Power_Spectrum_Normalizer_Allowed_Parameters</unitName>
   <type>void</type>
   <module>ISO_Varying_String, only : varying_string</module>
   <arguments>class    (*             ), intent(inout)                            :: self             </arguments>
   <arguments>type     (varying_string), intent(inout), dimension(:), allocatable :: allowedParameters</arguments>
   <arguments>character(len=*         ), intent(in   )                            :: sourceName       </arguments>
   <arguments>logical                  , intent(in   )                            :: objectsOnly      </arguments>
  </functionGlobal>
  !!]
  subroutine Power_Spectrum_Normalizer_Allowed_Parameters(self,allowedParameters,sourceName,objectsOnly)
    !!{RST
    Add the names of parameters allowed for the given :galacticus-class:`powerSpectrumNormalizer` object to
    ``allowedParameters``.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : varying_string
    implicit none
    class    (*             ), intent(inout)                            :: self
    type     (varying_string), intent(inout), dimension(:), allocatable :: allowedParameters
    character(len=*         ), intent(in   )                            :: sourceName
    logical                  , intent(in   )                            :: objectsOnly

    select type (self)
    class is (powerSpectrumNormalizerClass)
       call self%allowedParameters(allowedParameters,sourceName,objectsOnly)
    class default
       call Error_Report('object must be of `powerSpectrumNormalizerClass` class'//{introspection:location})
    end select
    return
  end subroutine Power_Spectrum_Normalizer_Allowed_Parameters

end module Power_Spectrum_Normalizer_Utilities

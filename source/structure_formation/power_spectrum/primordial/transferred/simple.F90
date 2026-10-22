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
  A simple transferred primordial power spectrum class.
  !!}

  use :: Linear_Growth           , only : linearGrowthClass
  use :: Power_Spectra_Primordial, only : powerSpectrumPrimordialClass
  use :: Transfer_Functions      , only : transferFunctionClass

  !![
  <powerSpectrumPrimordialTransferred name="powerSpectrumPrimordialTransferredSimple" docformat="rst">
   <description>
   Implements a simple transferred primordial power spectrum, :math:`P(k,t) = A T^2(k) D^2(k,t) P_\mathrm{prim}(k)`, where
   :math:`T(k)` is the transfer function (see :galacticus-class:`transferFunctionClass`), :math:`D(k,t)` is the linear growth
   factor (see :galacticus-class:`linearGrowthClass`) relative to that at the epoch at which the transfer function is defined,
   :math:`P_\mathrm{prim}(k)` is the primordial power spectrum (see :galacticus-class:`powerSpectrumPrimordialClass`), and the
   normalization, :math:`A`, is determined by a :galacticus-class:`powerSpectrumNormalizerClass` object (for example, to give
   a specified value of :math:`\sigma_8`).

   In a universe with non-zero curvature, :math:`K`, the transfer function excludes the factor :math:`1-3K/k^2` which arises from the curvature term in the Poisson equation (see :galacticus-class:`transferFunctionClass`), and this factor is not reapplied here. The resulting power spectrum is therefore the curvature-free spectrum, valid for :math:`k \gg |K|^{1/2}`. On larger scales the power spectrum is gauge-dependent, and integrals over it (such as the variance, :math:`\sigma^2(M)`) are computed as flat-space Fourier integrals, which do not account for the curvature of space. Reapplying the factor would cause such integrals to diverge as :math:`k \rightarrow 0`, since the power spectrum would then scale as :math:`k^{n_\mathrm{s}-4}` for :math:`k \ll |K|^{1/2}`. On the scales relevant to galaxy formation the difference is negligible (about :math:`4\times10^{-6}` in :math:`\sigma_8` for :math:`\Omega_k=0.1`).
   </description>
   <deepCopy>
     <ignore variables="powerSpectrumNormalizer_"/>
   </deepCopy>
   <assignment>
     <functionClass variables="powerSpectrumNormalizer_"/>
   </assignment>
  </powerSpectrumPrimordialTransferred>
  !!]
  type, extends(powerSpectrumPrimordialTransferredClass) :: powerSpectrumPrimordialTransferredSimple
     !!{RST
     A simple transferred primordial power spectrum class.
     !!}
     private
     class           (transferFunctionClass       ), pointer :: transferFunction_                => null()
     class           (linearGrowthClass           ), pointer :: linearGrowth_                    => null()
     class           (powerSpectrumPrimordialClass), pointer :: powerSpectrumPrimordial_         => null()
     ! The normalizer is held as an unlimited polymorphic object, and operated on through global functions, to avoid a circular
     ! module dependency (see Power_Spectrum_Normalizer_Utilities).
     class           (*                           ), pointer :: powerSpectrumNormalizer_         => null()
     double precision                                        :: timeTransferFunction                      , linearGrowthTransferFunction, &
          &                                                     timePrevious                              , linearGrowthPrevious        , &
          &                                                     normalization_
     logical                                                 :: linearGrothIsWavenumberDependent
   contains
     !![
     <methods docformat="rst">
       <method description="Return the unnormalized transferred primordial power spectrum at the given ``wavenumber`` and ``time``." method="powerUnnormalized"/>
       <method description="Return the normalization of the power spectrum, computing it if necessary."                               method="normalization"    />
     </methods>
     !!]
     final     ::                                simpleDestructor
     procedure :: power                       => simplePower
     procedure :: powerUnnormalized           => simplePowerUnnormalized
     procedure :: normalization               => simpleNormalization
     procedure :: logarithmicDerivative       => simpleLogarithmicDerivative
     procedure :: growthIsWavenumberDependent => simpleGrowthIsWavenumberDependent
     procedure :: descriptor                  => simpleDescriptor
     procedure :: deepCopy                    => simpleDeepCopy
     procedure :: deepCopyReset               => simpleDeepCopyReset
     procedure :: deepCopyFinalize            => simpleDeepCopyFinalize
  end type powerSpectrumPrimordialTransferredSimple

  interface powerSpectrumPrimordialTransferredSimple
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumPrimordialTransferredSimple` transferred primordial power spectrum class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface powerSpectrumPrimordialTransferredSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumPrimordialTransferredSimple` transferred primordial power spectrum class which takes a parameter set as input.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Construct_, Power_Spectrum_Normalizer_Release_
    use :: Input_Parameters, only : inputParameter                      , inputParameters
    implicit none
    type (powerSpectrumPrimordialTransferredSimple)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(transferFunctionClass                   ), pointer       :: transferFunction_
    class(powerSpectrumPrimordialClass            ), pointer       :: powerSpectrumPrimordial_
    class(linearGrowthClass                       ), pointer       :: linearGrowth_
    class(*                                       ), pointer       :: powerSpectrumNormalizer_

    !![
    <objectBuilder class="powerSpectrumPrimordial" name="powerSpectrumPrimordial_" source="parameters"/>
    <objectBuilder class="transferFunction"        name="transferFunction_"        source="parameters"/>
    <objectBuilder class="linearGrowth"            name="linearGrowth_"            source="parameters"/>
    !!]
    powerSpectrumNormalizer_ => Power_Spectrum_Normalizer_Construct_(parameters)
    self=powerSpectrumPrimordialTransferredSimple(powerSpectrumPrimordial_,transferFunction_,linearGrowth_,powerSpectrumNormalizer_)
    !![
    <inputParametersValidate source="parameters" extraAllowedNames="powerSpectrumNormalizer"/>
    <objectDestructor name="powerSpectrumPrimordial_"/>
    <objectDestructor name="transferFunction_"       />
    <objectDestructor name="linearGrowth_"           />
    !!]
    call Power_Spectrum_Normalizer_Release_(powerSpectrumNormalizer_)
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(powerSpectrumPrimordial_,transferFunction_,linearGrowth_,powerSpectrumNormalizer_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumPrimordialTransferredSimple` transferred primordial power
    spectrum class. The ``powerSpectrumNormalizer_`` argument must be of :galacticus-class:`powerSpectrumNormalizerClass` class.
    !!}
    implicit none
    type (powerSpectrumPrimordialTransferredSimple)                        :: self
    class(powerSpectrumPrimordialClass            ), intent(in   ), target :: powerSpectrumPrimordial_
    class(transferFunctionClass                   ), intent(in   ), target :: transferFunction_
    class(linearGrowthClass                       ), intent(in   ), target :: linearGrowth_
    class(*                                       ), intent(in   ), target :: powerSpectrumNormalizer_
    !![
    <constructorAssign variables="*powerSpectrumPrimordial_, *transferFunction_, *linearGrowth_, *powerSpectrumNormalizer_"/>
    !!]

    self       %timeTransferFunction            =self%transferFunction_%epochTime            (                              )
    self       %linearGrothIsWavenumberDependent=self%linearGrowth_    %isWaveNumberDependent(                              )
    if (.not.self%linearGrothIsWavenumberDependent) &
         & self%linearGrowthTransferFunction    =self%linearGrowth_    %value                (time=self%timeTransferFunction)
    self       %timePrevious                    =-huge(0.0d0)
    self       %linearGrowthPrevious            =-huge(0.0d0)
    self       %normalization_                  =-1.0d0
    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`powerSpectrumPrimordialTransferredSimple` transferred primordial power spectrum class.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Release_
    implicit none
    type(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%transferFunction_"       />
    <objectDestructor name="self%powerSpectrumPrimordial_"/>
    <objectDestructor name="self%linearGrowth_"           />
    !!]
    if (associated(self%powerSpectrumNormalizer_)) call Power_Spectrum_Normalizer_Release_(self%powerSpectrumNormalizer_)
    return
  end subroutine simpleDestructor

  double precision function simplePower(self,wavenumber,time)
    !!{RST
    Return the transferred primordial power spectrum at the given ``wavenumber``.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredSimple), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber, time

    simplePower=+self%normalization    (                ) &
         &      *self%powerUnnormalized(wavenumber,time)
    return
  end function simplePower

  double precision function simplePowerUnnormalized(self,wavenumber,time)
    !!{RST
    Return the unnormalized transferred primordial power spectrum at the given ``wavenumber``.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredSimple), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber  , time
    double precision                                                          :: linearGrowth

    if (self%linearGrothIsWavenumberDependent) then
       linearGrowth=+self%linearGrowth_%value(wavenumber=wavenumber,time=     time                ) &
            &       /self%linearGrowth_%value(wavenumber=wavenumber,time=self%timeTransferFunction)
    else
       if (time /= self%timePrevious) then
          self%timePrevious        =                                                      time
          self%linearGrowthPrevious=+self%linearGrowth_%value                       (time=time) &
               &                    /self              %linearGrowthTransferFunction
       end if
       linearGrowth=self%linearGrowthPrevious
    end if
    simplePowerUnnormalized=+(                                                                   &
         &                    +self%transferFunction_       %value       (wavenumber=wavenumber) &
         &                    *                              linearGrowth                        &
         &                   )**2                                                                &
         &                  *  self%powerSpectrumPrimordial_%power       (wavenumber=wavenumber)
    return
  end function simplePowerUnnormalized

  double precision function simpleNormalization(self)
    !!{RST
    Return the normalization of the power spectrum, computing it (using the normalizer) on the first call.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Normalization_
    implicit none
    class(powerSpectrumPrimordialTransferredSimple), intent(inout), target :: self

    if (self%normalization_ < 0.0d0) self%normalization_=Power_Spectrum_Normalizer_Normalization_(self%powerSpectrumNormalizer_,powerUnnormalized)
    simpleNormalization=self%normalization_
    return

  contains

    double precision function powerUnnormalized(wavenumber,time)
      !!{RST
      Return the unnormalized transferred primordial power spectrum, for use by the normalizer.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber, time

      powerUnnormalized=self%powerUnnormalized(wavenumber,time)
      return
    end function powerUnnormalized

  end function simpleNormalization

  double precision function simpleLogarithmicDerivative(self,wavenumber,time)
    !!{RST
    Return the logarithmic derivative of the transferred primordial power spectrum at the given ``wavenumber``.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredSimple), intent(inout) :: self
    double precision                                          , intent(in   ) :: wavenumber, time

    simpleLogarithmicDerivative=+2.0d0*self%transferFunction_       %logarithmicDerivative          (wavenumber=wavenumber                               ) &
         &                      +2.0d0*self%linearGrowth_           %logarithmicDerivativeWavenumber(wavenumber=wavenumber,time=     time                ) &
         &                      -2.0d0*self%linearGrowth_           %logarithmicDerivativeWavenumber(wavenumber=wavenumber,time=self%timeTransferFunction) &
         &                      +      self%powerSpectrumPrimordial_%logarithmicDerivative          (wavenumber=wavenumber                               )
    return
  end function simpleLogarithmicDerivative

  logical function simpleGrowthIsWavenumberDependent(self)
    !!{RST
    Return true if the growth of the power spectrum is wavenumber-dependent.
    !!}
    implicit none
    class(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self

    simpleGrowthIsWavenumberDependent=self%linearGrowth_%isWavenumberDependent()
    return
  end function simpleGrowthIsWavenumberDependent

  subroutine simpleDescriptor(self,descriptor,includeClass,includeFileModificationTimes,parameterName)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Descriptor_
    use :: Input_Parameters, only : inputParameters                      , descriptorParameterName
    implicit none
    class    (powerSpectrumPrimordialTransferredSimple), intent(inout)           :: self
    type     (inputParameters                         ), intent(inout)           :: descriptor
    logical                                            , intent(in   ), optional :: includeClass, includeFileModificationTimes
    character(len=*                                   ), intent(in   ), optional :: parameterName
    type     (inputParameters                         )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter(descriptorParameterName('powerSpectrumPrimordialTransferred',parameterName),'simple')
    parameters=descriptor%subparameters(descriptorParameterName('powerSpectrumPrimordialTransferred',parameterName))
    call self%powerSpectrumPrimordial_%descriptor(parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    call self%transferFunction_       %descriptor(parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    call self%linearGrowth_           %descriptor(parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    call Power_Spectrum_Normalizer_Descriptor_(self%powerSpectrumNormalizer_,parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    return
  end subroutine simpleDescriptor

  subroutine simpleDeepCopyReset(self)
    !!{RST
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Deep_Copy_Reset_
    implicit none
    class(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self

    self%copiedSelf => null()
    if (associated(self%transferFunction_       )) call self%transferFunction_       %deepCopyReset()
    if (associated(self%linearGrowth_           )) call self%linearGrowth_           %deepCopyReset()
    if (associated(self%powerSpectrumPrimordial_)) call self%powerSpectrumPrimordial_%deepCopyReset()
    if (associated(self%powerSpectrumNormalizer_)) call Power_Spectrum_Normalizer_Deep_Copy_Reset_(self%powerSpectrumNormalizer_)
    return
  end subroutine simpleDeepCopyReset

  subroutine simpleDeepCopyFinalize(self)
    !!{RST
    Finalize a deep copy of the object.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Deep_Copy_Finalize_
    implicit none
    class(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self

    if (associated(self%transferFunction_       )) call self%transferFunction_       %deepCopyFinalize()
    if (associated(self%linearGrowth_           )) call self%linearGrowth_           %deepCopyFinalize()
    if (associated(self%powerSpectrumPrimordial_)) call self%powerSpectrumPrimordial_%deepCopyFinalize()
    if (associated(self%powerSpectrumNormalizer_)) call Power_Spectrum_Normalizer_Deep_Copy_Finalize_(self%powerSpectrumNormalizer_)
    return
  end subroutine simpleDeepCopyFinalize

  subroutine simpleDeepCopy(self,destination)
    !!{RST
    Perform a deep copy of the object.
    !!}
    use :: Error           , only : Error_Report
    use :: Function_Classes, only : functionClass
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Deep_Copy_
    implicit none
    class  (powerSpectrumPrimordialTransferredSimple), intent(inout), target :: self
    class  (powerSpectrumPrimordialTransferredClass ), intent(inout)         :: destination
    integer                                                                  :: referenceCount

    call self%deepCopy_(destination)
    select type (destination)
    type is (powerSpectrumPrimordialTransferredSimple)
       ! The (shallow) assignment of self to destination in deepCopy_() acquired a reference to our normalizer. Release that
       ! reference, as the destination is given its own (deep) copy of the normalizer.
       if (associated(destination%powerSpectrumNormalizer_)) then
          select type (powerSpectrumNormalizer_ => destination%powerSpectrumNormalizer_)
          class is (functionClass)
             referenceCount=powerSpectrumNormalizer_%referenceCountDecrement()
          end select
       end if
       nullify(destination%powerSpectrumNormalizer_)
       if (associated(self%powerSpectrumNormalizer_)) call Power_Spectrum_Normalizer_Deep_Copy_(self%powerSpectrumNormalizer_,destination%powerSpectrumNormalizer_)
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine simpleDeepCopy

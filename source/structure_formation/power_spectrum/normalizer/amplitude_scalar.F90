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

  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Transfer_Functions  , only : transferFunctionClass

  !![
  <powerSpectrumNormalizer name="powerSpectrumNormalizerAmplitudeScalar" docformat="rst">
   <description>
   A power spectrum normalizer which normalizes the power spectrum to a specified amplitude, :math:`A_\mathrm{s}=`\
   ``[amplitudeScalar]``, of the dimensionless power spectrum of primordial comoving curvature perturbations at the pivot
   wavenumber :math:`k_*=`\ ``[wavenumberPivot]``, :math:`\mathcal{P}_\mathcal{R}(k_*)=A_\mathrm{s}`. The linear matter power
   spectrum at the present day is

   .. math::

      P(k) = 2 \pi^2 \mathcal{T}_0^2 k \mathcal{P}_\mathcal{R}(k) T^2(k),

   where :math:`T(k)` is the transfer function (normalized to unity as :math:`k \rightarrow 0`), and
   :math:`\mathcal{T}_0=\lim_{k\rightarrow 0}[-\delta_\mathrm{m}(k)/k^2]` is the large-scale amplitude of the matter density
   perturbation at the present day per unit primordial curvature perturbation. :math:`\mathcal{T}_0` is computed by the CLASS
   Boltzmann code for the cosmological parameters of the model, so includes the growth of perturbations through the radiation
   era as computed by CLASS. It does not depend on the primordial power spectrum, nor on the transfer function used by the
   model - the shape of the primordial power spectrum (including its spectral index and any running) is that of the model, and
   :math:`T(k)` is that of the model's transfer function.

   Since the unnormalized power spectrum is :math:`\propto k \mathcal{P}_\mathcal{R}(k) T^2(k)`, the normalization factor is
   found by matching the above expression to it at :math:`k_*`. Only the relative growth of perturbations,
   :math:`D(t)/D(t_0)`, is taken from the model's linear growth factor.
   </description>
  </powerSpectrumNormalizer>
  !!]
  type, extends(powerSpectrumNormalizerClass) :: powerSpectrumNormalizerAmplitudeScalar
     !!{RST
     A power spectrum normalizer which normalizes to a specified amplitude of the primordial scalar power spectrum.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class           (transferFunctionClass   ), pointer :: transferFunction_    => null()
     double precision                                    :: amplitudeScalar               , wavenumberPivot
   contains
     final     ::                  amplitudeScalarDestructor
     procedure :: normalization => amplitudeScalarNormalization
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
    double precision                                                        :: amplitudeScalar     , wavenumberPivot

    !![
    <inputParameter docformat="rst">
      <name>amplitudeScalar</name>
      <source>parameters</source>
      <description>
      The amplitude of the dimensionless primordial curvature power spectrum at the pivot wavenumber, :math:`A_\mathrm{s}=\mathcal{P}_\mathcal{R}(k_*)`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>wavenumberPivot</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <description>
      The pivot wavenumber, :math:`k_*` (in units of Mpc\ :math:`^{-1}`), at which the amplitude of the primordial curvature power spectrum is specified.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunction_"    source="parameters"/>
    !!]
    self=powerSpectrumNormalizerAmplitudeScalar(amplitudeScalar,wavenumberPivot,cosmologyParameters_,cosmologyFunctions_,transferFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunction_"   />
    !!]
    return
  end function amplitudeScalarConstructorParameters

  function amplitudeScalarConstructorInternal(amplitudeScalar,wavenumberPivot,cosmologyParameters_,cosmologyFunctions_,transferFunction_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumNormalizerAmplitudeScalar` power spectrum normalizer class.
    !!}
    implicit none
    type            (powerSpectrumNormalizerAmplitudeScalar)                        :: self
    double precision                                        , intent(in   )         :: amplitudeScalar     , wavenumberPivot
    class           (cosmologyParametersClass              ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class           (transferFunctionClass                 ), intent(in   ), target :: transferFunction_
    !![
    <constructorAssign variables="amplitudeScalar, wavenumberPivot, *cosmologyParameters_, *cosmologyFunctions_, *transferFunction_"/>
    !!]

    return
  end function amplitudeScalarConstructorInternal

  subroutine amplitudeScalarDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`powerSpectrumNormalizerAmplitudeScalar` power spectrum normalizer class.
    !!}
    implicit none
    type(powerSpectrumNormalizerAmplitudeScalar), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%transferFunction_"   />
    !!]
    return
  end subroutine amplitudeScalarDestructor

  double precision function amplitudeScalarNormalization(self,powerUnnormalized) result(normalization)
    !!{RST
    Return the factor by which the given power spectrum must be multiplied to normalize it to the specified amplitude of the
    primordial scalar power spectrum.
    !!}
    use :: Interfaces_CLASS        , only : Interface_CLASS_Transfer_Function_Amplitude
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumNormalizerAmplitudeScalar  ), intent(inout), target :: self
    procedure       (powerSpectrumNormalizerPowerUnnormalized)                        :: powerUnnormalized
    double precision                                                                  :: amplitudeTransfer, timePresent

    amplitudeTransfer=Interface_CLASS_Transfer_Function_Amplitude(self%cosmologyParameters_)
    timePresent      =self%cosmologyFunctions_%cosmicTime(1.0d0)
    normalization    =+2.0d0                                           &
         &            *Pi                                         **2  &
         &            *amplitudeTransfer                          **2  &
         &            *self%wavenumberPivot                            &
         &            *self%amplitudeScalar                            &
         &            *self%transferFunction_%value(self%wavenumberPivot)**2 &
         &            /powerUnnormalized(self%wavenumberPivot,timePresent)
    return
  end function amplitudeScalarNormalization

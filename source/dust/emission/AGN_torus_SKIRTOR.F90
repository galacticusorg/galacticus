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
  Implements a spectrum of thermal emission from the dusty torus around an AGN using the SKIRTOR library of models.
  !!}

  use :: ISO_Varying_String     , only : varying_string
  use :: Numerical_Interpolation, only : interpolator

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumAGNTorusSKIRTOR" docformat="rst">
   <description>
   Thermal emission from the dusty torus around an AGN, taken from the SKIRTOR library of radiative transfer models of a
   two-phase (clumpy plus smooth) torus :cite:p:`stalevski_3d_2012,stalevski_dust_2016`. The model is selected by the
   edge-on optical depth at :math:`9.7\,\mu\hbox{m}` (``[depthOpticalEdgeOn]``), the exponents of the radial
   (``[exponentDensityRadial]``) and polar (``[exponentDensityPolar]``) gradients of dust density, the ratio of outer to
   inner radius (``[ratioRadii]``), and the angle between the equatorial plane and the edge of the torus
   (``[angleTorus]``), each of which must be one of the values tabulated.

   The emission of the model, averaged over all directions, sets the shape of the spectrum, which is normalized to the
   luminosity absorbed. The shape is therefore fixed: it does not follow the angle of the torus chosen for each galaxy by
   :galacticus-class:`dustAttenuationAGNTorusSKIRTOR`, nor depend on the inclination at which the torus is seen. Between
   the tabulated wavelengths the spectrum is interpolated as a power law, and it is normalized by integrating that
   interpolant analytically, so that it integrates to exactly the absorbed luminosity. It is zero outside the range of
   wavelength over which the model emits.

   The temperature of the dust in the torus is set by its distance from the accretion disk, not by the luminosity it
   absorbs, so the mass of dust is not used. Give the torus phase a zero ``fractionsMassPhase`` in
   :galacticus-class:`nodePropertyExtractorSEDDustEmission`.

   The library is read from ``[fileName]``, in the ``dust/agnTorus`` directory of the dynamic datasets path; if it is not
   present it is downloaded from ``[url]``.
   </description>
  </dustEmissionSpectrum>
  !!]
  type, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumAGNTorusSKIRTOR
     !!{RST
     A spectrum of thermal emission from the dusty torus around an AGN using the SKIRTOR library of models.
     !!}
     private
     type            (varying_string)                            :: fileName              , url
     double precision                                            :: depthOpticalEdgeOn    , exponentDensityRadial, &
          &                                                         exponentDensityPolar  , ratioRadii           , &
          &                                                         angleTorus
     ! The tabulated wavelengths (in Å) over which the model emits, and the logarithm of the spectrum at each, normalized
     ! to unit luminosity (so in Hz⁻¹).
     double precision                , allocatable, dimension(:) :: wavelength            , spectrumLogarithmic
     type            (interpolator  )                            :: interpolatorWavelength
   contains
     procedure :: luminosity => agnTorusSKIRTORLuminosity
  end type dustEmissionSpectrumAGNTorusSKIRTOR

  interface dustEmissionSpectrumAGNTorusSKIRTOR
     !!{RST
     Constructors for the :galacticus-class:`dustEmissionSpectrumAGNTorusSKIRTOR` dust emission spectrum class.
     !!}
     module procedure agnTorusSKIRTORConstructorParameters
     module procedure agnTorusSKIRTORConstructorInternal
  end interface dustEmissionSpectrumAGNTorusSKIRTOR

contains

  function agnTorusSKIRTORConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustEmissionSpectrumAGNTorusSKIRTOR` dust emission spectrum class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: ISO_Varying_String, only : var_str
    implicit none
    type            (dustEmissionSpectrumAGNTorusSKIRTOR)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    type            (varying_string                     )                :: fileName            , url
    double precision                                                     :: depthOpticalEdgeOn  , exponentDensityRadial, &
         &                                                                  exponentDensityPolar, ratioRadii           , &
         &                                                                  angleTorus

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <defaultValue>var_str('SKIRTOR_2016-7-18.hdf5')</defaultValue>
      <description>
      The name of the file holding the SKIRTOR tabulation, in the ``dust/agnTorus`` directory of the dynamic datasets
      path.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>url</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>
      The URL from which to download ``fileName`` if it is not already present. If ``none``, the file must be supplied
      by other means.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>depthOpticalEdgeOn</name>
      <defaultValue>7.0d0</defaultValue>
      <defaultSource>The default of the SKIRTOR module of CIGALE :cite:p:`boquien_cigale_2019`</defaultSource>
      <description>
      The optical depth at :math:`9.7\,\mu\hbox{m}` of the torus seen edge-on. One of 3, 5, 7, 9 or 11.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentDensityRadial</name>
      <defaultValue>1.0d0</defaultValue>
      <defaultSource>The default of the SKIRTOR module of CIGALE :cite:p:`boquien_cigale_2019`</defaultSource>
      <description>
      The exponent, :math:`p`, of the radial gradient of dust density, :math:`\rho \propto r^{-p}`. One of 0, 0.5, 1 or
      1.5.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentDensityPolar</name>
      <defaultValue>1.0d0</defaultValue>
      <defaultSource>The default of the SKIRTOR module of CIGALE :cite:p:`boquien_cigale_2019`</defaultSource>
      <description>
      The index, :math:`q`, of the polar gradient of dust density, :math:`\rho \propto \exp(-q|\cos\theta|)`. One of 0,
      0.5, 1 or 1.5.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>ratioRadii</name>
      <defaultValue>20.0d0</defaultValue>
      <defaultSource>The default of the SKIRTOR module of CIGALE :cite:p:`boquien_cigale_2019`</defaultSource>
      <description>
      The ratio of the outer to the inner radius of the torus. One of 10, 20 or 30.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>angleTorus</name>
      <defaultValue>40.0d0</defaultValue>
      <defaultSource>The default of the SKIRTOR module of CIGALE :cite:p:`boquien_cigale_2019`</defaultSource>
      <description>
      The angle, in degrees, between the equatorial plane and the edge of the torus. One of 10, 20, 30, 40, 50, 60, 70
      or 80.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustEmissionSpectrumAGNTorusSKIRTOR(fileName,url,depthOpticalEdgeOn,exponentDensityRadial,exponentDensityPolar,ratioRadii,angleTorus)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function agnTorusSKIRTORConstructorParameters

  function agnTorusSKIRTORConstructorInternal(fileName,url,depthOpticalEdgeOn,exponentDensityRadial,exponentDensityPolar,ratioRadii,angleTorus) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustEmissionSpectrumAGNTorusSKIRTOR` dust emission spectrum class.
    The model is read, restricted to the range of wavelength over which it emits, and normalized to unit luminosity.
    !!}
    use :: Dust_AGN_Torus_SKIRTOR          , only : skirtorRead
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : degreesToRadians
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Units       , only : metersToAngstroms
    use :: Table_Labels                    , only : extrapolationTypeFix
    implicit none
    type            (dustEmissionSpectrumAGNTorusSKIRTOR)                                :: self
    type            (varying_string                     ), intent(in   )                 :: fileName                   , url
    double precision                                     , intent(in   )                 :: depthOpticalEdgeOn         , exponentDensityRadial, &
         &                                                                                  exponentDensityPolar       , ratioRadii           , &
         &                                                                                  angleTorus
    ! Tolerance on the match between the requested angle of the torus and a tabulated one, in radians.
    double precision                                     , parameter                     :: toleranceAngle      =1.0d-6
    double precision                                     , allocatable  , dimension(:  ) :: angleTorusModel            , wavelength           , &
         &                                                                                  inclination                , luminosityLogarithmic
    double precision                                     , allocatable  , dimension(:,:) :: emissionDustMean
    logical                                              , allocatable  , dimension(:  ) :: isEmitting
    double precision                                                                     :: normalization
    integer                                                                              :: indexModel
    !![
    <constructorAssign variables="fileName, url, depthOpticalEdgeOn, exponentDensityRadial, exponentDensityPolar, ratioRadii, angleTorus"/>
    !!]

    call skirtorRead(fileName,url,depthOpticalEdgeOn,exponentDensityRadial,exponentDensityPolar,ratioRadii,angleTorusModel,wavelength,inclination,emissionDustMean=emissionDustMean)
    if (all(abs(angleTorusModel-angleTorus*degreesToRadians) > toleranceAngle)) &
         & call Error_Report('`angleTorus` is not one of the tabulated values'//{introspection:location})
    indexModel=minloc(abs(angleTorusModel-angleTorus*degreesToRadians),dim=1)
    ! Restrict to the range over which the model emits. Emission is non-zero over a single contiguous range.
    isEmitting=emissionDustMean(:,indexModel) > 0.0d0
    if (count(isEmitting) < 2) call Error_Report('the model emits at fewer than two tabulated wavelengths'//{introspection:location})
    self%wavelength      =pack(wavelength                                          ,isEmitting)
    luminosityLogarithmic=pack(log(max(emissionDustMean(:,indexModel),tiny(0.0d0))),isEmitting)
    ! Convert from νL_ν to L_ν, which is proportional to λ νL_ν, and normalize to unit luminosity. The spectrum is
    ! interpolated as a power law between the tabulated wavelengths, and normalized by integrating that interpolant.
    normalization           =dustEmissionIntegralPowerLaw(self%wavelength,luminosityLogarithmic)
    self%spectrumLogarithmic=+luminosityLogarithmic                             &
         &                   +log(self%wavelength/speedLight/metersToAngstroms) &
         &                   -log(normalization)
    self%interpolatorWavelength=interpolator(log(self%wavelength),extrapolationType=extrapolationTypeFix)
    return
  end function agnTorusSKIRTORConstructorInternal

  function agnTorusSKIRTORLuminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Return the luminosity per unit frequency emitted by the torus, normalized to the luminosity it absorbs.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (dustEmissionSpectrumAGNTorusSKIRTOR), intent(inout)                               :: self
    double precision                                     , intent(in   ), dimension(:                ) :: wavelengths
    double precision                                     , intent(in   )                               :: luminosityAbsorbed, massDust, &
         &                                                                                                time
    double precision                                                    , dimension(size(wavelengths)) :: luminosity
    double precision                                                    , dimension(0:1              ) :: weights
    integer         (c_size_t                           )                                              :: j
    integer                                                                                            :: i
    !$GLC attributes unused :: massDust, time

    luminosity=0.0d0
    if (luminosityAbsorbed <= 0.0d0) return
    do i=1,size(wavelengths)
       if (wavelengths(i) < self%wavelength(1) .or. wavelengths(i) > self%wavelength(size(self%wavelength))) cycle
       call self%interpolatorWavelength%linearFactors(log(wavelengths(i)),j,weights)
       luminosity(i)=+luminosityAbsorbed                            &
            &        *exp(                                          &
            &             +weights(0)*self%spectrumLogarithmic(j  ) &
            &             +weights(1)*self%spectrumLogarithmic(j+1) &
            &            )
    end do
    return
  end function agnTorusSKIRTORLuminosity

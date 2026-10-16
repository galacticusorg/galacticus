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
  Implements the template spectra of thermal emission from dust of :cite:t:`dale_infrared_2002`.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumDaleHelou2002" docformat="rst">
   <description>
   The template infrared spectra of star-forming galaxies of :cite:t:`dale_infrared_2002`. Each is built from the
   emission of dust heated by radiation fields of intensity :math:`U`, relative to the local interstellar radiation
   field, combined over a power-law distribution of dust mass,

   .. math::

      \mathrm{d}M_\mathrm{dust} \propto U^{-\alpha}\,\mathrm{d}U, \qquad 0.3 \le U \le 10^5,

   following :cite:t:`dale_infrared_2001`. The single parameter ``alpha`` sets the shape: small values describe galaxies
   heated intensely, with a warm infrared peak, and large values quiescent galaxies, with a cool one. Normal
   star-forming galaxies span :math:`1 \lesssim \alpha \lesssim 2.5`. The 64 tabulated templates, spanning
   :math:`0.0625 \le \alpha \le 4`, are interpolated linearly in :math:`\alpha`, in the logarithm of the luminosity.

   The templates carry no absolute scale that relates the luminosity to the mass of dust, so the mass of dust is not
   used: the shape is fixed by ``alpha``, and the spectrum is normalized to the absorbed luminosity.

   Only the part of each template which is emission from dust is used, from :math:`3` to :math:`1100\,\mu\hbox{m}`,
   and the spectrum is zero outside that range. The distributed templates extend further. Shortward of
   :math:`3\,\mu\hbox{m}` they include a component whose shape does not depend on :math:`\alpha`---presumably
   starlight---which carries from 5 to 38% of the tabulated power; including it would count starlight twice. Beyond
   about :math:`1\,\hbox{mm}` they are extended to radio wavelengths using the far-infrared--radio correlation, which
   describes synchrotron and free-free emission, not emission from dust.

   Between the tabulated wavelengths the spectrum is interpolated as a power law. It is normalized by integrating that
   interpolant analytically, so that it integrates to exactly the absorbed luminosity whatever wavelengths it is
   evaluated at.
   </description>
  </dustEmissionSpectrum>
  !!]
  type, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumDaleHelou2002
     !!{RST
     The template spectra of thermal emission from dust of :cite:t:`dale_infrared_2002`.
     !!}
     private
     double precision                                          :: alpha
     ! The tabulated wavelengths (in Å) within the dust emission range, and the logarithm of the spectrum at each,
     ! normalized to unit luminosity (so in Hz⁻¹).
     double precision              , allocatable, dimension(:) :: wavelength            , spectrumLogarithmic
     type            (interpolator)                            :: interpolatorWavelength
   contains
     procedure :: luminosity => daleHelou2002Luminosity
  end type dustEmissionSpectrumDaleHelou2002

  interface dustEmissionSpectrumDaleHelou2002
     !!{RST
     Constructors for the :galacticus-class:`dustEmissionSpectrumDaleHelou2002` dust emission spectrum class.
     !!}
     module procedure daleHelou2002ConstructorParameters
     module procedure daleHelou2002ConstructorInternal
  end interface dustEmissionSpectrumDaleHelou2002

contains

  function daleHelou2002ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustEmissionSpectrumDaleHelou2002` dust emission spectrum class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustEmissionSpectrumDaleHelou2002)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: alpha

    !![
    <inputParameter docformat="rst">
      <name>alpha</name>
      <defaultValue>2.0d0</defaultValue>
      <description>
      The exponent, :math:`\alpha`, of the distribution of dust mass over the intensity of the heating radiation field,
      :math:`\mathrm{d}M_\mathrm{dust} \propto U^{-\alpha}\,\mathrm{d}U`, which selects the template. Must lie between
      0.0625 and 4.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustEmissionSpectrumDaleHelou2002(alpha)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function daleHelou2002ConstructorParameters

  function daleHelou2002ConstructorInternal(alpha) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustEmissionSpectrumDaleHelou2002` dust emission spectrum class. The
    templates are read, the one selected by ``alpha`` is formed, restricted to the range of dust emission, and
    normalized to unit luminosity.
    !!}
    use :: Error                       , only : Error_Report
    use :: HDF5_Access                 , only : hdf5Access
    use :: Input_Paths                 , only : inputPath           , pathTypeDataStatic
    use :: IO_HDF5                     , only : hdf5File
    use :: ISO_Varying_String          , only : char
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms
    use :: Table_Labels                , only : extrapolationTypeFix
    implicit none
    type            (dustEmissionSpectrumDaleHelou2002)                                :: self
    double precision                                   , intent(in   )                 :: alpha
    double precision                                   , allocatable  , dimension(:  ) :: wavelength           , alphas               , &
         &                                                                                luminosityLogarithmic
    double precision                                   , allocatable  , dimension(:,:) :: luminosity
    logical                                            , allocatable  , dimension(:  ) :: isDust
    double precision                                                                   :: wavelengthDustMinimum, wavelengthDustMaximum, &
         &                                                                                weight               , normalization
    integer                                                                            :: i
    !![
    <constructorAssign variables="alpha"/>
    !!]

    !$ call hdf5Access%set()
    hdf5ReadScope : block
      type(hdf5File) :: file
      file=hdf5File(char(inputPath(pathTypeDataStatic)//'dust/emission/daleHelou2002.hdf5'),readOnly=.true.)
      call file%readDataset  ('wavelength'           ,wavelength           )
      call file%readDataset  ('alpha'                ,alphas               )
      call file%readDataset  ('luminosity'           ,luminosity           )
      call file%readAttribute('wavelengthDustMinimum',wavelengthDustMinimum)
      call file%readAttribute('wavelengthDustMaximum',wavelengthDustMaximum)
    end block hdf5ReadScope
    !$ call hdf5Access%unset()
    ! The dataset is written (alpha, wavelength) as a row-major reader sees it, so the Fortran dimensions are the reverse.
    if (any(shape(luminosity) /= [size(wavelength),size(alphas)])) &
         & call Error_Report('`luminosity` does not have the shape implied by the axes'//{introspection:location})
    if (alpha < alphas(1) .or. alpha > alphas(size(alphas)))       &
         & call Error_Report('`alpha` lies outside the tabulated range of 0.0625 to 4' //{introspection:location})
    ! Interpolate linearly in alpha, in the logarithm of the luminosity.
    i                    =max(1,min(size(alphas)-1,count(alphas <= alpha)))
    weight               =+(alpha      -alphas(i)) &
         &                /(alphas(i+1)-alphas(i))
    luminosityLogarithmic=+(1.0d0-weight)*log(luminosity(:,i  )) &
         &                +       weight *log(luminosity(:,i+1))
    ! Restrict to the range of dust emission.
    isDust                  =wavelength >= wavelengthDustMinimum .and. wavelength <= wavelengthDustMaximum
    self%wavelength         =pack(wavelength           ,isDust)
    luminosityLogarithmic   =pack(luminosityLogarithmic,isDust)
    ! Convert from νL_ν to L_ν, which is proportional to λ νL_ν, and normalize to unit luminosity. The spectrum is
    ! interpolated as a power law between the tabulated wavelengths, and normalized by integrating that interpolant.
    normalization           =dustEmissionIntegralPowerLaw(self%wavelength,luminosityLogarithmic)
    self%spectrumLogarithmic=+luminosityLogarithmic                             &
         &                   +log(self%wavelength/speedLight/metersToAngstroms) &
         &                   -log(normalization)
    self%interpolatorWavelength=interpolator(log(self%wavelength),extrapolationType=extrapolationTypeFix)
    return
  end function daleHelou2002ConstructorInternal

  function daleHelou2002Luminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Return the luminosity per unit frequency emitted by the dust, normalized to the luminosity it absorbs.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (dustEmissionSpectrumDaleHelou2002), intent(inout)                               :: self
    double precision                                   , intent(in   ), dimension(:                ) :: wavelengths
    double precision                                   , intent(in   )                               :: luminosityAbsorbed, massDust, &
         &                                                                                              time
    double precision                                                  , dimension(size(wavelengths)) :: luminosity
    double precision                                                  , dimension(0:1              ) :: weights
    integer         (c_size_t                         )                                              :: j
    integer                                                                                          :: i
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
  end function daleHelou2002Luminosity

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
  Implements the spectra of thermal emission from dust of :cite:t:`draine_infrared_2007`.
  !!}

  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator

  !![
  <enumeration docformat="rst">
   <name>draineLi2007GrainModel</name>
   <description>
   Enumerates the grain models of :cite:t:`draine_infrared_2007`, named for the galaxy whose extinction they reproduce
   and the fraction of dust mass in PAHs, :math:`q_\mathrm{PAH}`.
   </description>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="milkyWay00"            />
   <entry label="milkyWay10"            />
   <entry label="milkyWay20"            />
   <entry label="milkyWay30"            />
   <entry label="milkyWay40"            />
   <entry label="milkyWay50"            />
   <entry label="milkyWay60"            />
   <entry label="largeMagellanicCloud00"/>
   <entry label="largeMagellanicCloud05"/>
   <entry label="largeMagellanicCloud10"/>
   <entry label="smallMagellanicCloud"  />
  </enumeration>
  !!]

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumDraineLi2007" docformat="rst">
   <description>
   Thermal emission from the silicate-graphite-PAH dust models of :cite:t:`draine_infrared_2007`, which follow the
   stochastic heating of small grains and PAHs as well as the equilibrium heating of large grains. The dust is heated
   by starlight with the spectrum of the local interstellar radiation field, scaled by an intensity :math:`U`, with
   the mass of dust distributed over intensity as

   .. math::

      \frac{\mathrm{d}M_\mathrm{dust}}{\mathrm{d}U} = (1-\gamma) M_\mathrm{dust} \delta(U-U_\mathrm{min}) + \gamma M_\mathrm{dust} \frac{U^{-2}}{U_\mathrm{min}^{-1}-U_\mathrm{max}^{-1}}, \qquad U_\mathrm{min} \le U \le U_\mathrm{max}:

   a fraction :math:`1-\gamma` of the dust in the diffuse interstellar medium, heated by :math:`U_\mathrm{min}`, and the
   rest in regions of more intense heating, such as photodissociation regions.

   The ``grainModel`` sets the grain composition and size distribution: seven Milky Way models, with PAH mass fractions
   :math:`q_\mathrm{PAH}` from 0.47 to 4.58% (``milkyWay00`` to ``milkyWay60``), three for the Large Magellanic
   Cloud, and one for the Small Magellanic Cloud. The fraction ``fractionMassPowerLaw`` is :math:`\gamma`, and
   ``intensityMaximum`` is :math:`U_\mathrm{max}`, one of the tabulated :math:`10^3`, :math:`10^4`, :math:`10^5`, or
   :math:`10^6`.

   The dust radiates the power it absorbs, so the luminosity it emits per unit mass is proportional to the mean
   intensity heating it,

   .. math::

      \frac{L_\mathrm{abs}}{M_\mathrm{dust}} = P_0 \langle U \rangle, \qquad \langle U \rangle = (1-\gamma) U_\mathrm{min} + \gamma \frac{\ln(U_\mathrm{max}/U_\mathrm{min})}{U_\mathrm{min}^{-1}-U_\mathrm{max}^{-1}},

   where :math:`P_0 \approx 136\,L_\odot\,M_\odot^{-1}` for ``milkyWay60``. By default :math:`U_\mathrm{min}` is found
   from this energy balance, using the absorbed luminosity and the mass of dust given, and the tabulated power of each
   model and dust-to-hydrogen mass ratio (:cite:t:`draine_infrared_2007`, Table 3). Setting ``intensityMinimum``
   instead fixes :math:`U_\mathrm{min}`, and the mass of dust is then not used.

   The models are tabulated for :math:`0.1 \le U_\mathrm{min} \le 25`, and spectra are interpolated between them
   linearly in :math:`\ln U_\mathrm{min}`. Dust heated more weakly or more strongly than that range allows takes the
   shape of the spectrum at its nearest end. The emitted luminosity is still exactly that absorbed, but the spectrum
   peaks at longer wavelengths than it should: for :math:`\gamma = 0.01` and :math:`U_\mathrm{max} = 10^6` the upper
   limit is :math:`L_\mathrm{abs}/M_\mathrm{dust} \approx 3.7 \times 10^3\,L_\odot\,M_\odot^{-1}`, which intensely
   star-forming galaxies can exceed. The cosmic microwave background is not included as a source of heating.

   The spectra span :math:`1\,\mu\hbox{m}` to :math:`1\,\hbox{cm}`, and are zero outside that range. Between the
   tabulated wavelengths each is interpolated as a power law, and normalized by integrating that interpolant
   analytically, so that the spectrum integrates to exactly the absorbed luminosity whatever wavelengths it is
   evaluated at.
   </description>
  </dustEmissionSpectrum>
  !!]
  type, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumDraineLi2007
     !!{RST
     The spectra of thermal emission from dust of :cite:t:`draine_infrared_2007`.
     !!}
     private
     type            (enumerationDraineLi2007GrainModelType)                              :: grainModel
     double precision                                                                     :: fractionMassPowerLaw          , intensityMinimum_     , &
          &                                                                                  intensityMaximum
     ! The tabulated wavelengths (in Å) and minimum intensities, the logarithm of the luminosity emitted per unit dust
     ! mass (in L☉/M☉) at each minimum intensity, and the logarithm of the spectrum at each wavelength and minimum
     ! intensity, normalized to unit luminosity (so in Hz⁻¹).
     double precision                                       , allocatable, dimension(:  ) :: wavelength                    , intensitiesMinimum    , &
          &                                                                                  luminositySpecificLogarithmic
     double precision                                       , allocatable, dimension(:,:) :: spectrumLogarithmic
     type            (interpolator                         )                              :: interpolatorWavelength        , interpolatorLuminosity
     ! For a fixed minimum intensity, the bracketing tabulated intensity and the interpolating weights.
     integer         (c_size_t                             )                              :: indexIntensityFixed
     double precision                                                    , dimension(0:1) :: weightsIntensityFixed
   contains
     !![
     <methods docformat="rst">
       <method method="intensityMinimum" description="Return the minimum intensity, :math:`U_\mathrm{min}`, of the radiation heating dust of the given mass absorbing the given luminosity, limited to the tabulated range."/>
     </methods>
     !!]
     procedure :: luminosity       => draineLi2007Luminosity
     procedure :: intensityMinimum => draineLi2007IntensityMinimum
     procedure :: weights          => draineLi2007Weights
  end type dustEmissionSpectrumDraineLi2007

  interface dustEmissionSpectrumDraineLi2007
     !!{RST
     Constructors for the :galacticus-class:`dustEmissionSpectrumDraineLi2007` dust emission spectrum class.
     !!}
     module procedure draineLi2007ConstructorParameters
     module procedure draineLi2007ConstructorInternal
  end interface dustEmissionSpectrumDraineLi2007

contains

  function draineLi2007ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustEmissionSpectrumDraineLi2007` dust emission spectrum class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: ISO_Varying_String, only : char          , var_str        , varying_string
    implicit none
    type            (dustEmissionSpectrumDraineLi2007)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    type            (varying_string                  )                :: grainModel
    double precision                                                  :: fractionMassPowerLaw, intensityMinimum_, &
         &                                                               intensityMaximum

    !![
    <inputParameter docformat="rst">
      <name>grainModel</name>
      <defaultValue>var_str('milkyWay60')</defaultValue>
      <description>
      The grain model: one of ``milkyWay00``, ``milkyWay10``, ``milkyWay20``, ``milkyWay30``, ``milkyWay40``,
      ``milkyWay50``, ``milkyWay60``, ``largeMagellanicCloud00``, ``largeMagellanicCloud05``,
      ``largeMagellanicCloud10``, or ``smallMagellanicCloud``.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>fractionMassPowerLaw</name>
      <defaultValue>0.01d0</defaultValue>
      <description>
      The fraction, :math:`\gamma`, of the dust mass heated by a power-law distribution of intensities, rather than by
      the minimum intensity alone.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>intensityMinimum</name>
      <variable>intensityMinimum_</variable>
      <defaultValue>-1.0d0</defaultValue>
      <description>
      If positive, the minimum intensity, :math:`U_\mathrm{min}`, of the radiation heating the dust, overriding the
      value which energy balance would give. The mass of dust is then not used. Must lie between 0.1 and 25.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>intensityMaximum</name>
      <defaultValue>1.0d6</defaultValue>
      <description>
      The maximum intensity, :math:`U_\mathrm{max}`, of the power-law distribution of intensities: one of
      :math:`10^3`, :math:`10^4`, :math:`10^5`, or :math:`10^6`.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustEmissionSpectrumDraineLi2007(enumerationDraineLi2007GrainModelEncode(char(grainModel),includesPrefix=.false.),fractionMassPowerLaw,intensityMinimum_,intensityMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function draineLi2007ConstructorParameters

  function draineLi2007ConstructorInternal(grainModel,fractionMassPowerLaw,intensityMinimum_,intensityMaximum) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustEmissionSpectrumDraineLi2007` dust emission spectrum class. The
    spectra of the selected grain model are read, combined for the given :math:`\gamma` and :math:`U_\mathrm{max}` at
    each tabulated :math:`U_\mathrm{min}`, and normalized to unit luminosity.
    !!}
    use :: Error                           , only : Error_Report
    use :: HDF5_Access                     , only : hdf5Access
    use :: Input_Paths                     , only : inputPath           , pathTypeDataStatic
    use :: IO_HDF5                         , only : hdf5File            , hdf5Group
    use :: ISO_Varying_String              , only : char
    use :: Numerical_Constants_Astronomical, only : luminositySolar     , massSolar
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Units       , only : ergs                , metersToAngstroms
    use :: Table_Labels                    , only : extrapolationTypeFix
    implicit none
    type            (dustEmissionSpectrumDraineLi2007     )                                :: self
    type            (enumerationDraineLi2007GrainModelType), intent(in   )                 :: grainModel
    double precision                                       , intent(in   )                 :: fractionMassPowerLaw   , intensityMinimum_, &
         &                                                                                    intensityMaximum
    type            (interpolator                         )                                :: interpolatorIntensity
    double precision                                       , allocatable, dimension(:    ) :: intensitiesMaximum     , powerSingle      , &
         &                                                                                    emissivityLogarithmic
    double precision                                       , allocatable, dimension(:,:  ) :: emissivitySingle       , powerPowerLaw
    double precision                                       , allocatable, dimension(:,:,:) :: emissivityPowerLaw
    double precision                                                                       :: massDustPerHydrogen
    integer                                                                                :: i                      , indexMaximum
    !![
    <constructorAssign variables="grainModel, fractionMassPowerLaw, intensityMinimum_, intensityMaximum"/>
    !!]

    if (fractionMassPowerLaw < 0.0d0 .or. fractionMassPowerLaw > 1.0d0) &
         & call Error_Report('`fractionMassPowerLaw` must lie between 0 and 1'//{introspection:location})
    !$ call hdf5Access%set()
    hdf5ReadScope: block
      type(hdf5File ) :: file
      type(hdf5Group) :: group
      file =hdf5File(char(inputPath(pathTypeDataStatic)//'dust/emission/draineLi2007.hdf5'),readOnly=.true.)
      call file %readDataset  ('wavelength'         ,self%wavelength         )
      call file %readDataset  ('intensityMinimum'   ,self%intensitiesMinimum )
      call file %readDataset  ('intensityMaximum'   ,     intensitiesMaximum )
      group=file%openGroup(char(enumerationDraineLi2007GrainModelDecode(grainModel,includePrefix=.false.)))
      call group%readDataset  ('emissivitySingle'   ,     emissivitySingle   )
      call group%readDataset  ('emissivityPowerLaw' ,     emissivityPowerLaw )
      call group%readDataset  ('powerSingle'        ,     powerSingle        )
      call group%readDataset  ('powerPowerLaw'      ,     powerPowerLaw      )
      call group%readAttribute('massDustPerHydrogen',     massDustPerHydrogen)
    end block hdf5ReadScope
    !$ call hdf5Access%unset()
    ! Datasets are written with wavelength last, as a row-major reader sees them, so the Fortran dimensions are the
    ! reverse of that.
    if (any(shape(emissivitySingle  ) /= [size(self%wavelength),size(self%intensitiesMinimum)                         ])) &
         & call Error_Report('`emissivitySingle` does not have the shape implied by the axes'  //{introspection:location})
    if (any(shape(emissivityPowerLaw) /= [size(self%wavelength),size(self%intensitiesMinimum),size(intensitiesMaximum)])) &
         & call Error_Report('`emissivityPowerLaw` does not have the shape implied by the axes'//{introspection:location})
    ! Find the tabulated maximum intensity.
    indexMaximum=0
    do i=1,size(intensitiesMaximum)
       if (abs(log(intensitiesMaximum(i)/intensityMaximum)) < 1.0d-6) indexMaximum=i
    end do
    if (indexMaximum == 0) call Error_Report('`intensityMaximum` must be one of the tabulated values 10³, 10⁴, 10⁵, or 10⁶'//{introspection:location})
    ! Combine the single-intensity and power-law spectra, which are each per H nucleon of the same dust, weighted by the
    ! fractions of dust mass heated in each way. Normalize each to unit luminosity, converting from ν dP/dν to L_ν, which
    ! is proportional to λ ν dP/dν.
    allocate(self%spectrumLogarithmic          (size(self%wavelength),size(self%intensitiesMinimum)))
    allocate(self%luminositySpecificLogarithmic(                      size(self%intensitiesMinimum)))
    do i=1,size(self%intensitiesMinimum)
       emissivityLogarithmic=log(                                                                   &
            &                    +(1.0d0-fractionMassPowerLaw)*emissivitySingle  (:,i             ) &
            &                    +       fractionMassPowerLaw *emissivityPowerLaw(:,i,indexMaximum) &
            &                   )
       self%spectrumLogarithmic(:,i)=+emissivityLogarithmic                                                    &
            &                        +log(self%wavelength/speedLight/metersToAngstroms)                        &
            &                        -log(dustEmissionIntegralPowerLaw(self%wavelength,emissivityLogarithmic))
       ! The power radiated per unit dust mass, in L☉/M☉, from the tabulated power per H nucleon.
       self%luminositySpecificLogarithmic(i)=log(                                                                 &
            &                                    +(                                                               &
            &                                      +(1.0d0-fractionMassPowerLaw)*powerSingle  (i             )    &
            &                                      +       fractionMassPowerLaw *powerPowerLaw(i,indexMaximum)    &
            &                                     )                                                               &
            &                                    *ergs                                                            &
            &                                    /massDustPerHydrogen                                             &
            &                                    /massHydrogenAtom                                                &
            &                                    *massSolar                                                       &
            &                                    /luminositySolar                                                 &
            &                                   )
    end do
    if (any(self%luminositySpecificLogarithmic(2:) <= self%luminositySpecificLogarithmic(:size(self%intensitiesMinimum)-1))) &
         & call Error_Report('luminosity per unit dust mass does not increase with minimum intensity'//{introspection:location})
    self%interpolatorWavelength=interpolator(log(self%wavelength                   ),extrapolationType=extrapolationTypeFix)
    self%interpolatorLuminosity=interpolator(    self%luminositySpecificLogarithmic ,extrapolationType=extrapolationTypeFix)
    ! For a fixed minimum intensity, find the interpolating weights once.
    self%indexIntensityFixed  =0_c_size_t
    self%weightsIntensityFixed=0.0d0
    if (intensityMinimum_ > 0.0d0) then
       if (intensityMinimum_ < self%intensitiesMinimum(1) .or. intensityMinimum_ > self%intensitiesMinimum(size(self%intensitiesMinimum))) &
            & call Error_Report('`intensityMinimum` lies outside the tabulated range of 0.1 to 25'//{introspection:location})
       interpolatorIntensity=interpolator(log(self%intensitiesMinimum))
       call interpolatorIntensity%linearFactors(log(intensityMinimum_),self%indexIntensityFixed,self%weightsIntensityFixed)
    end if
    return
  end function draineLi2007ConstructorInternal

  subroutine draineLi2007Weights(self,luminosityAbsorbed,massDust,index,weights)
    !!{RST
    Find the tabulated minimum intensities bracketing that of dust of mass ``massDust`` absorbing a luminosity
    ``luminosityAbsorbed``, and the weights, linear in :math:`\ln U_\mathrm{min}`, with which to interpolate between them.

    The luminosity per unit dust mass is tabulated at each :math:`U_\mathrm{min}`, and interpolating
    :math:`\ln U_\mathrm{min}` linearly in the logarithm of that luminosity gives weights which are the same in both.
    Beyond the tabulation the nearest end is used.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (dustEmissionSpectrumDraineLi2007), intent(inout)                 :: self
    double precision                                  , intent(in   )                 :: luminosityAbsorbed, massDust
    integer         (c_size_t                        ), intent(  out)                 :: index
    double precision                                  , intent(  out), dimension(0:1) :: weights
    double precision                                                                  :: luminositySpecific
    integer         (c_size_t                        )                                :: countIntensities

    if (self%intensityMinimum_ > 0.0d0) then
       index  =self%indexIntensityFixed
       weights=self%weightsIntensityFixed
       return
    end if
    if (massDust <= 0.0d0) &
         & call Error_Report('dust absorbs luminosity but has no mass, so no heating intensity can be found'//{introspection:location})
    countIntensities  =size(self%intensitiesMinimum,kind=c_size_t)
    luminositySpecific=log(luminosityAbsorbed/massDust)
    if      (luminositySpecific <= self%luminositySpecificLogarithmic(1               )) then
       index  =1_c_size_t
       weights=[1.0d0,0.0d0]
    else if (luminositySpecific >= self%luminositySpecificLogarithmic(countIntensities)) then
       index  =countIntensities-1_c_size_t
       weights=[0.0d0,1.0d0]
    else
       call self%interpolatorLuminosity%linearFactors(luminositySpecific,index,weights)
    end if
    return
  end subroutine draineLi2007Weights

  double precision function draineLi2007IntensityMinimum(self,luminosityAbsorbed,massDust) result(intensityMinimum)
    !!{RST
    Return the minimum intensity, :math:`U_\mathrm{min}`, of the radiation heating dust of mass ``massDust`` (in
    :math:`M_\odot`) absorbing a luminosity ``luminosityAbsorbed`` (in :math:`L_\odot`), limited to the tabulated range.
    !!}
    implicit none
    class           (dustEmissionSpectrumDraineLi2007), intent(inout)  :: self
    double precision                                  , intent(in   )  :: luminosityAbsorbed, massDust
    double precision                                  , dimension(0:1) :: weights
    integer         (c_size_t                        )                 :: index

    call self%weights(luminosityAbsorbed,massDust,index,weights)
    intensityMinimum=exp(                                                  &
         &               +weights(0)*log(self%intensitiesMinimum(index  )) &
         &               +weights(1)*log(self%intensitiesMinimum(index+1)) &
         &              )
    return
  end function draineLi2007IntensityMinimum

  function draineLi2007Luminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Return the luminosity per unit frequency emitted by the dust, normalized to the luminosity it absorbs.

    The normalized spectra at the two bracketing minimum intensities are each interpolated in wavelength and then
    combined with weights summing to one, so the result integrates to exactly the absorbed luminosity.
    !!}
    implicit none
    class           (dustEmissionSpectrumDraineLi2007), intent(inout)                               :: self
    double precision                                  , intent(in   ), dimension(:                ) :: wavelengths
    double precision                                  , intent(in   )                               :: luminosityAbsorbed, massDust        , &
         &                                                                                             time
    double precision                                                 , dimension(size(wavelengths)) :: luminosity
    double precision                                                 , dimension(0:1              ) :: weightsIntensity  , weightsWavelength
    integer         (c_size_t                        )                                              :: indexIntensity    , indexWavelength
    integer                                                                                         :: i
    !$GLC attributes unused :: time

    luminosity=0.0d0
    if (luminosityAbsorbed <= 0.0d0) return
    call self%weights(luminosityAbsorbed,massDust,indexIntensity,weightsIntensity)
    do i=1,size(wavelengths)
       if (wavelengths(i) < self%wavelength(1) .or. wavelengths(i) > self%wavelength(size(self%wavelength))) cycle
       call self%interpolatorWavelength%linearFactors(log(wavelengths(i)),indexWavelength,weightsWavelength)
       luminosity(i)=+luminosityAbsorbed                                                                       &
            &        *(                                                                                        &
            &          +weightsIntensity(0)                                                                    &
            &          *exp(                                                                                   &
            &               +weightsWavelength(0)*self%spectrumLogarithmic(indexWavelength  ,indexIntensity  ) &
            &               +weightsWavelength(1)*self%spectrumLogarithmic(indexWavelength+1,indexIntensity  ) &
            &              )                                                                                   &
            &          +weightsIntensity(1)                                                                    &
            &          *exp(                                                                                   &
            &               +weightsWavelength(0)*self%spectrumLogarithmic(indexWavelength  ,indexIntensity+1) &
            &               +weightsWavelength(1)*self%spectrumLogarithmic(indexWavelength+1,indexIntensity+1) &
            &              )                                                                                   &
            &         )
    end do
    return
  end function draineLi2007Luminosity

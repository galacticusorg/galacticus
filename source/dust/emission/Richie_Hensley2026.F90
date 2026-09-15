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
  Implements emission from PAHs in the single-photon approximation of :cite:t:`richie_pah_2026`.
  !!}

  !![
  <enumeration docformat="rst">
   <name>richieHensley2026SizeDistribution</name>
   <description>
   Enumerates the PAH size distributions of :cite:t:`draine_excitation_2021` for which emission is tabulated.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="small"   />
   <entry label="standard"/>
   <entry label="large"   />
  </enumeration>
  <enumeration docformat="rst">
   <name>richieHensley2026IonizationFunction</name>
   <description>
   Enumerates the PAH ionization functions of :cite:t:`draine_excitation_2021` for which emission is tabulated.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="low"     />
   <entry label="standard"/>
   <entry label="high"    />
  </enumeration>
  !!]

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumRichieHensley2026" docformat="rst">
   <description>
   Emission from polycyclic aromatic hydrocarbons (PAHs) heated by the light the dust absorbs, in the single-photon
   approximation of :cite:t:`richie_pah_2026`, together with emission from the larger grains which absorb the rest of
   that light.

   A PAH is small enough that it cools completely between the absorption of one photon and the next, so its emission is
   a sum of the spectra which follow single absorptions, and depends on the spectrum of the heating light, not only on
   its intensity. :cite:t:`richie_pah_2026` tabulate those spectra for neutral and ionized PAHs of sizes from 3.5 to
   100 Å, absorbing photons in bins 1% wide in wavelength from 912 Å to 10.1 :math:`\mu\hbox{m}`. Here they are used
   integrated over the size distribution and ionization function of PAHs selected by ``sizeDistribution`` and
   ``ionizationFunction`` (:cite:t:`draine_excitation_2021`), per unit energy absorbed in each bin.

   The absorbed light is shared between PAHs and other grains in proportion to how strongly each absorbs it. In each bin
   PAHs absorb a fraction

   .. math::

      f_\mathrm{PAH}(\lambda) = \frac{\kappa_\mathrm{abs,PAH}(\lambda)}{\kappa_\mathrm{abs,PAH}(\lambda)+\kappa_\mathrm{abs,astrodust}(\lambda)},

   where :math:`\kappa_\mathrm{abs,PAH}` is the absorption opacity of the tabulated PAH population, the same one whose
   emission is tabulated, and :math:`\kappa_\mathrm{abs,astrodust}` is that of the astrodust grains of
   :cite:t:`hensley_astrodust_2023`, extrapolated as a power law below 0.1 :math:`\mu\hbox{m}`. For the interstellar
   radiation field PAHs absorb about 30% of the light, and about 40% of that from a population 10 Myr old. The rest of
   the absorbed light, including any outside the tabulated range, heats the larger grains, whose emission is given by
   the ``dustEmissionSpectrum`` of this class (a :galacticus-class:`dustEmissionSpectrumBlackBodyModified`, for
   example), which is given the whole mass of dust. Since the models of :cite:t:`draine_infrared_2007` already include
   PAHs, a :galacticus-class:`dustEmissionSpectrumDraineLi2007` should not be used for the larger grains.

   The emission of PAHs integrates to exactly the energy they absorb. It is tabulated in bins of wavelength of
   resolution 500 from 1 to 20 :math:`\mu\hbox{m}`, and 100 elsewhere from 0.1 :math:`\mu\hbox{m}` to 1 cm, and is
   integrated exactly over any interval of wavelength requested.

   The single-photon approximation agrees with calculations including multi-photon heating to within 5% below
   20 :math:`\mu\hbox{m}` for radiation fields of the intensity of the local interstellar radiation field, within 10%
   below 10 :math:`\mu\hbox{m}` for intensities up to 100 times that, and within 10% below 6 :math:`\mu\hbox{m}` up to
   :math:`10^4` times that (:cite:t:`richie_pah_2026`). It is poor beyond 20 :math:`\mu\hbox{m}`, where larger PAHs
   dominate, and in intense radiation fields; neither is checked here.

   The emission depends on the spectrum of the absorbed light, so this class requires it: it can be used through
   ``luminosityIntegrated``, as by :galacticus-class:`nodePropertyExtractorSEDDustEmission`, but not through
   ``luminosity``, which is given only the total absorbed luminosity.
   </description>
  </dustEmissionSpectrum>
  !!]
  type, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumRichieHensley2026
     !!{RST
     Emission from PAHs in the single-photon approximation of :cite:t:`richie_pah_2026`, together with emission from
     larger grains.
     !!}
     private
     class           (dustEmissionSpectrumClass                         ), pointer                     :: dustEmissionSpectrum_ => null()
     type            (enumerationRichieHensley2026SizeDistributionType  )                              :: sizeDistribution
     type            (enumerationRichieHensley2026IonizationFunctionType)                              :: ionizationFunction
     ! Edges of the bins of absorbed and emitted wavelength (in Å), the fraction of the absorbed light taken by PAHs in
     ! each absorbed bin, and the emission, the mean of λ p_λ per unit ln λ in each emitted bin per unit energy absorbed
     ! in each absorbed bin (the first index runs over emitted bins, the second over absorbed bins).
     double precision                                                    , allocatable, dimension(:  ) :: wavelengthsAbsorbed            , wavelengthsEmitted, &
          &                                                                                               fractionPAH
     double precision                                                    , allocatable, dimension(:,:) :: emission
   contains
     final     ::                         richieHensley2026Destructor
     procedure :: luminosity           => richieHensley2026Luminosity
     procedure :: luminosityIntegrated => richieHensley2026LuminosityIntegrated
  end type dustEmissionSpectrumRichieHensley2026

  interface dustEmissionSpectrumRichieHensley2026
     !!{RST
     Constructors for the :galacticus-class:`dustEmissionSpectrumRichieHensley2026` dust emission spectrum class.
     !!}
     module procedure richieHensley2026ConstructorParameters
     module procedure richieHensley2026ConstructorInternal
  end interface dustEmissionSpectrumRichieHensley2026

contains

  function richieHensley2026ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustEmissionSpectrumRichieHensley2026` dust emission spectrum class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: ISO_Varying_String, only : char          , var_str        , varying_string
    implicit none
    type (dustEmissionSpectrumRichieHensley2026)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(dustEmissionSpectrumClass            ), pointer       :: dustEmissionSpectrum_
    type (varying_string                       )                :: sizeDistribution     , ionizationFunction

    !![
    <inputParameter docformat="rst">
      <name>sizeDistribution</name>
      <defaultValue>var_str('standard')</defaultValue>
      <description>
      The PAH size distribution of :cite:t:`draine_excitation_2021`: one of ``small``, ``standard``, or ``large``.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>ionizationFunction</name>
      <defaultValue>var_str('standard')</defaultValue>
      <description>
      The PAH ionization function of :cite:t:`draine_excitation_2021`: one of ``low``, ``standard``, or ``high``.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustEmissionSpectrum" name="dustEmissionSpectrum_" source="parameters"/>
    !!]
    self=dustEmissionSpectrumRichieHensley2026(                                                                                                       &
         &                                     enumerationRichieHensley2026SizeDistributionEncode  (char(sizeDistribution  ),includesPrefix=.false.), &
         &                                     enumerationRichieHensley2026IonizationFunctionEncode(char(ionizationFunction),includesPrefix=.false.), &
         &                                     dustEmissionSpectrum_                                                                                  &
         &                                    )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="dustEmissionSpectrum_"/>
    !!]
    return
  end function richieHensley2026ConstructorParameters

  function richieHensley2026ConstructorInternal(sizeDistribution,ionizationFunction,dustEmissionSpectrum_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustEmissionSpectrumRichieHensley2026` dust emission spectrum class.
    The tabulated emission and absorption opacities of the selected PAH population are read, and the fraction of
    absorbed light taken by PAHs found.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: Input_Paths       , only : inputPath   , pathTypeDataStatic
    use :: IO_HDF5           , only : hdf5File    , hdf5Group
    use :: ISO_Varying_String, only : char        , operator(//)      , var_str, varying_string
    implicit none
    type            (dustEmissionSpectrumRichieHensley2026             )                                :: self
    type            (enumerationRichieHensley2026SizeDistributionType  ), intent(in   )                 :: sizeDistribution
    type            (enumerationRichieHensley2026IonizationFunctionType), intent(in   )                 :: ionizationFunction
    class           (dustEmissionSpectrumClass                         ), intent(in   ), target         :: dustEmissionSpectrum_
    double precision                                                    , allocatable  , dimension(:  ) :: opacityPAH           , opacityAstrodust
    type            (varying_string                                    )                                :: label
    !![
    <constructorAssign variables="sizeDistribution, ionizationFunction, *dustEmissionSpectrum_"/>
    !!]

    ! Construct the name of the group holding the selected population.
    select case (sizeDistribution%ID)
    case (richieHensley2026SizeDistributionSmall     %ID)
       label=var_str('sizeSmall'   )
    case (richieHensley2026SizeDistributionStandard  %ID)
       label=var_str('sizeStandard')
    case (richieHensley2026SizeDistributionLarge     %ID)
       label=var_str('sizeLarge'   )
    case default
       call Error_Report('unrecognized size distribution'//{introspection:location})
    end select
    select case (ionizationFunction%ID)
    case (richieHensley2026IonizationFunctionLow     %ID)
       label=label//'IonizationLow'
    case (richieHensley2026IonizationFunctionStandard%ID)
       label=label//'IonizationStandard'
    case (richieHensley2026IonizationFunctionHigh    %ID)
       label=label//'IonizationHigh'
    case default
       call Error_Report('unrecognized ionization function'//{introspection:location})
    end select
    ! Read the tabulation. The HDF5 objects are scoped within a block so that they are finalized before the lock is
    ! released.
    !$ call hdf5Access%set()
    hdf5ReadScope: block
      type(hdf5File ) :: file
      type(hdf5Group) :: group
      file =hdf5File(char(inputPath(pathTypeDataStatic)//'dust/emission/richieHensley2026.hdf5'),readOnly=.true.)
      call file %readDataset('wavelengthAbsorbedEdges'   ,self%wavelengthsAbsorbed)
      call file %readDataset('wavelengthEmittedEdges'    ,self%wavelengthsEmitted )
      call file %readDataset('opacityAbsorptionAstrodust',     opacityAstrodust   )
      group=file%openGroup(char(label))
      call group%readDataset('opacityAbsorptionPAH'      ,     opacityPAH         )
      call group%readDataset('emission'                  ,self%emission           )
    end block hdf5ReadScope
    !$ call hdf5Access%unset()
    ! The emission dataset is written (absorbed, emitted) as a row-major reader sees it, so the Fortran dimensions are the
    ! reverse of that.
    if     (                                                                                                                 &
         &   size(opacityPAH      ) /= size(self%wavelengthsAbsorbed)-1                                                      &
         &  .or.                                                                                                             &
         &   size(opacityAstrodust) /= size(self%wavelengthsAbsorbed)-1                                                      &
         &  .or.                                                                                                             &
         &   any(shape(self%emission) /= [size(self%wavelengthsEmitted)-1,size(self%wavelengthsAbsorbed)-1])                 &
         & ) call Error_Report('tabulated PAH emission does not have the shape implied by its axes'//{introspection:location})
    self%fractionPAH=opacityPAH/(opacityPAH+opacityAstrodust)
    return
  end function richieHensley2026ConstructorInternal

  subroutine richieHensley2026Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustEmissionSpectrumRichieHensley2026` dust emission spectrum class.
    !!}
    implicit none
    type(dustEmissionSpectrumRichieHensley2026), intent(inout) :: self

    !![
    <objectDestructor name="self%dustEmissionSpectrum_"/>
    !!]
    return
  end subroutine richieHensley2026Destructor

  function richieHensley2026Luminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Report that the emission of PAHs can not be found from the total absorbed luminosity alone.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (dustEmissionSpectrumRichieHensley2026), intent(inout)                  :: self
    double precision                                       , intent(in   ), dimension(:   ) :: wavelengths
    double precision                                       , intent(in   )                  :: luminosityAbsorbed, massDust, &
         &                                                                                     time
    double precision                                       , dimension(size(wavelengths))   :: luminosity
    !$GLC attributes unused :: self, luminosityAbsorbed, massDust, time

    luminosity=0.0d0
    call Error_Report('emission from PAHs depends on the spectrum of the absorbed light, so must be found through `luminosityIntegrated`'//{introspection:location})
    return
  end function richieHensley2026Luminosity

  function richieHensley2026LuminosityIntegrated(self,wavelengthsMinimum,wavelengthsMaximum,wavelengthsHeating,luminositiesAbsorbed,massDust,time) result(integral)
    !!{RST
    Return :math:`\int L_\nu\,\mathrm{d}\ln\lambda` over each of the given intervals of wavelength, from PAHs and from the
    larger grains.

    The luminosity absorbed in each interval of heating wavelength is taken to be spread uniformly in :math:`\ln\lambda`
    across it, and so shared among the tabulated absorbed bins it overlaps. PAHs take the fraction
    :math:`f_\mathrm{PAH}` of each share, and the rest is returned to the interval from which it came to heat the larger
    grains. The PAH emission is constant in :math:`\lambda p_\lambda` per unit :math:`\ln\lambda` across each emitted bin,
    so :math:`L_\nu = \lambda p_\lambda/\nu` there, and its integral over :math:`\ln\lambda` across any part of the bin is
    that constant times the range of :math:`\lambda` spanned, divided by :math:`c`.
    !!}
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms
    implicit none
    class           (dustEmissionSpectrumRichieHensley2026), intent(inout)                                             :: self
    double precision                                       , intent(in   ), dimension(:                              ) :: wavelengthsMinimum   , wavelengthsMaximum  , &
         &                                                                                                                wavelengthsHeating   , luminositiesAbsorbed
    double precision                                       , intent(in   )                                             :: massDust             , time
    double precision                                                      , dimension(size(wavelengthsMinimum     )  ) :: integral
    double precision                                                      , dimension(size(luminositiesAbsorbed   )  ) :: luminositiesRemaining
    double precision                                                      , dimension(size(self%fractionPAH       )  ) :: luminositiesPAH
    double precision                                                      , dimension(size(self%wavelengthsEmitted)-1) :: emitted
    double precision                                                                                                   :: widthLogarithmic     , overlap             , &
         &                                                                                                                share
    integer                                                                                                            :: i                    , j                   , &
         &                                                                                                                countAbsorbed        , countEmitted

    countAbsorbed        =size(self%fractionPAH       )
    countEmitted         =size(self%wavelengthsEmitted)-1
    luminositiesRemaining=luminositiesAbsorbed
    luminositiesPAH      =0.0d0
    ! Share the absorbed light among the tabulated absorbed bins, and between PAHs and larger grains.
    do j=1,size(luminositiesAbsorbed)
       if (luminositiesAbsorbed(j) <= 0.0d0                                                                                              ) cycle
       if (wavelengthsHeating(j+1) <= self%wavelengthsAbsorbed(1) .or. wavelengthsHeating(j) >= self%wavelengthsAbsorbed(countAbsorbed+1)) cycle
       widthLogarithmic=log(wavelengthsHeating(j+1)/wavelengthsHeating(j))
       do i=richieHensley2026Bin(self%wavelengthsAbsorbed,wavelengthsHeating(j)),countAbsorbed
          if (self%wavelengthsAbsorbed(i) >= wavelengthsHeating(j+1)) exit
          overlap=log(min(wavelengthsHeating(j+1),self%wavelengthsAbsorbed(i+1))/max(wavelengthsHeating(j),self%wavelengthsAbsorbed(i)))
          if (overlap <= 0.0d0) cycle
          share                   =+luminositiesAbsorbed (j) &
               &                   *overlap                  &
               &                   /widthLogarithmic
          luminositiesPAH      (i)=+luminositiesPAH      (i) &
               &                   +self%fractionPAH     (i) &
               &                   *share
          luminositiesRemaining(j)=+luminositiesRemaining(j) &
               &                   -self%fractionPAH     (i) &
               &                   *share
       end do
    end do
    luminositiesRemaining=max(luminositiesRemaining,0.0d0)
    ! The emission of PAHs, λ p_λ per unit ln λ in each emitted bin, in L☉.
    emitted =matmul(self%emission,luminositiesPAH)
    ! Integrate over each requested interval.
    integral=0.0d0
    if (any(luminositiesPAH > 0.0d0)) then
       do j=1,size(wavelengthsMinimum)
          if (wavelengthsMaximum(j) <= self%wavelengthsEmitted(1) .or. wavelengthsMinimum(j) >= self%wavelengthsEmitted(countEmitted+1)) cycle
          do i=richieHensley2026Bin(self%wavelengthsEmitted,wavelengthsMinimum(j)),countEmitted
             if (self%wavelengthsEmitted(i) >= wavelengthsMaximum(j)) exit
             overlap=+min(wavelengthsMaximum(j),self%wavelengthsEmitted(i+1)) &
                  &  -max(wavelengthsMinimum(j),self%wavelengthsEmitted(i  ))
             if (overlap <= 0.0d0) cycle
             integral(j)=+integral(j)       &
                  &      +emitted (i)       &
                  &      *overlap           &
                  &      /speedLight        &
                  &      /metersToAngstroms
          end do
       end do
    end if
    ! Add the emission of the larger grains, heated by the rest of the absorbed light.
    integral=+integral                                                                                                                                      &
         &   +self%dustEmissionSpectrum_%luminosityIntegrated(wavelengthsMinimum,wavelengthsMaximum,wavelengthsHeating,luminositiesRemaining,massDust,time)
    return
  end function richieHensley2026LuminosityIntegrated

  integer function richieHensley2026Bin(edges,wavelength) result(index)
    !!{RST
    Return the index of the bin, bounded by the ascending ``edges``, which contains ``wavelength``, limited to the range
    of bins, by bisection.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: edges
    double precision, intent(in   )               :: wavelength
    integer                                       :: lower     , upper, &
         &                                           middle

    lower=1
    upper=size(edges)
    if (wavelength <= edges(lower)) then
       index=1
       return
    end if
    if (wavelength >= edges(upper)) then
       index=size(edges)-1
       return
    end if
    do while (upper-lower > 1)
       middle=(lower+upper)/2
       if (edges(middle) <= wavelength) then
          lower=middle
       else
          upper=middle
       end if
    end do
    index=lower
    return
  end function richieHensley2026Bin

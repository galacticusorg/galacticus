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
  Implements a dust attenuation class describing the dusty torus around an AGN using the SKIRTOR library of models.
  !!}

  use :: Black_Hole_Torus_Covering_Fractions, only : blackHoleTorusCoveringFractionClass
  use :: Galactic_Inclinations              , only : galacticInclinationClass
  use :: ISO_Varying_String                 , only : varying_string

  !![
  <dustAttenuation name="dustAttenuationAGNTorusSKIRTOR" docformat="rst">
   <description>
   The dusty torus around the accretion disk of an AGN, using the SKIRTOR library of radiative transfer models of a
   two-phase (clumpy plus smooth) torus :cite:p:`stalevski_3d_2012,stalevski_dust_2016`. The models are selected by the
   edge-on optical depth at :math:`9.7\,\mu\hbox{m}` (``[depthOpticalEdgeOn]``), the exponents of the radial
   (``[exponentDensityRadial]``) and polar (``[exponentDensityPolar]``) gradients of dust density, and the ratio of outer
   to inner radius (``[ratioRadii]``), each of which must be one of the values tabulated. The angle of the torus is set
   for each galaxy from its covering fraction (see below).

   **Transmission.** SKIRTOR tabulates, for each model and inclination :math:`i`, the flux of the accretion disk seen
   directly (:math:`F_\mathrm{direct}`), seen after scattering by the dust (:math:`F_\mathrm{scattered}`), and as it
   would be with no dust (:math:`F_\mathrm{transparent}`). The transmission is

   .. math::

      T(\lambda,i) = \frac{F_\mathrm{direct}(\lambda,i)}{F_\mathrm{transparent}(\lambda,i)} + \frac{F_\mathrm{scattered}(\lambda,i)}{\langle F_\mathrm{transparent}(\lambda) \rangle},

   where :math:`\langle\rangle` denotes an average over all directions. The first term is the transmission along the
   line of sight, which depends only on the dust, and so applies to an accretion disk which emits isotropically---as
   assumed here---although the accretion disk of SKIRTOR does not :cite:p:`netzer_quasar_1987`. The second is light
   scattered towards the observer: since that was emitted in other directions it is counted relative to the mean flux of
   the accretion disk. Light scattered by the torus is therefore not counted as absorbed. At wavelengths shorter than
   about :math:`0.01\,\mu\hbox{m}` Monte Carlo noise in the scattered flux of some models can make the transmission,
   averaged over orientation, exceed one; where it does, the scattered term is reduced to make that average exactly one,
   so that the absorbed luminosity is never negative.

   :math:`T` is interpolated linearly in inclination and in the logarithm of wavelength, and held fixed beyond the
   tabulated range of wavelength (:math:`0.001` to :math:`1000\,\mu\hbox{m}`).

   **Covering fraction.** The covering fraction of each model is taken to be the fraction of the sky obscured in the
   :math:`V` band,

   .. math::

      f_\mathrm{cov} = 1 - \int_0^{\pi/2} \frac{F_\mathrm{direct}(\lambda_\mathrm{V},i)}{F_\mathrm{transparent}(\lambda_\mathrm{V},i)} \sin i\,\mathrm{d}i,

   with the same interpolation in inclination. This exceeds :math:`\sin(\hbox{angle of the torus})`, the covering
   fraction of the tabulated geometry, because the tabulation, spaced by :math:`10^\circ` in inclination, can not resolve
   the edge of the torus, which interpolation spreads over that interval; defined this way, the covering fraction is
   that of the attenuation actually applied. For each galaxy the two models whose covering fractions bracket that given
   by the ``[blackHoleTorusCoveringFraction]`` object are combined linearly, which---since :math:`f_\mathrm{cov}` is
   linear in :math:`T`---gives exactly that covering fraction. A covering fraction outside the range spanned by the
   models (for the default parameters, :math:`0.26` to :math:`0.98`) is moved to its nearer end. The covering fraction
   must increase with the angle of the torus; for two sets of parameters (``[depthOpticalEdgeOn]=3``,
   ``[exponentDensityRadial]=0``, ``[exponentDensityPolar]=1.5`` and ``[ratioRadii]=10`` or :math:`30`) Monte Carlo
   noise makes it decrease slightly between the two largest angles but one, and those sets can not be used.

   **What is attenuated.** As for :galacticus-class:`dustAttenuationAGNTorus`, only the continuum emission of the
   accretion disk of the central black hole (``emissionSourceAccretionDisk`` emitted by ``componentTypeBlackHole``) is
   attenuated, and the axis of the torus is taken to be that of its host galaxy, seen at the inclination given by the
   ``[galacticInclination]`` object (or by the ``inclination`` passed in when averaging over orientation).

   **Tabulation.** The library is read from ``[fileName]``, in the ``dust/agnTorus`` directory of the dynamic datasets
   path; if it is not present it is downloaded from ``[url]``. The file is made from the SKIRTOR library, as distributed
   by its authors, by ``dataPreparation/convertSKIRTOR.py`` of the ``galacticusDevTools`` repository
   (https://github.com/galacticusorg/galacticusDevTools).
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationAGNTorusSKIRTOR
     !!{RST
     A dust attenuation class describing the dusty torus around an AGN using the SKIRTOR library of models.
     !!}
     private
     class           (galacticInclinationClass           ), pointer                       :: galacticInclination_            => null()
     class           (blackHoleTorusCoveringFractionClass), pointer                       :: blackHoleTorusCoveringFraction_ => null()
     type            (varying_string                     )                                :: fileName                                 , url
     double precision                                                                     :: depthOpticalEdgeOn                       , exponentDensityRadial, &
          &                                                                                  exponentDensityPolar                     , ratioRadii
     logical                                                                              :: inclinationAvailable
     ! The tabulation: wavelengths (in Å) and their logarithms, inclinations (in radians), the transmission indexed
     ! (wavelength, inclination, model), and the covering fraction of each model.
     double precision                                     , allocatable, dimension(:    ) :: wavelength                               , wavelengthLogarithmic, &
          &                                                                                  inclination                              , coveringFractions
     double precision                                     , allocatable, dimension(:,:,:) :: transmissionTabulated
   contains
     final     ::                           agnTorusSKIRTORDestructor
     procedure :: transmission           => agnTorusSKIRTORTransmission
     procedure :: isOrientationDependent => agnTorusSKIRTORIsOrientationDependent
  end type dustAttenuationAGNTorusSKIRTOR

  interface dustAttenuationAGNTorusSKIRTOR
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationAGNTorusSKIRTOR` dust attenuation class.
     !!}
     module procedure agnTorusSKIRTORConstructorParameters
     module procedure agnTorusSKIRTORConstructorInternal
  end interface dustAttenuationAGNTorusSKIRTOR

contains

  function agnTorusSKIRTORConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationAGNTorusSKIRTOR` dust attenuation class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: ISO_Varying_String, only : var_str
    implicit none
    type            (dustAttenuationAGNTorusSKIRTOR     )                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (galacticInclinationClass           ), pointer       :: galacticInclination_
    class           (blackHoleTorusCoveringFractionClass), pointer       :: blackHoleTorusCoveringFraction_
    type            (varying_string                     )                :: fileName                       , url
    double precision                                                     :: depthOpticalEdgeOn             , exponentDensityRadial, &
         &                                                                  exponentDensityPolar           , ratioRadii

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
    <objectBuilder class="galacticInclination"            name="galacticInclination_"            source="parameters"/>
    <objectBuilder class="blackHoleTorusCoveringFraction" name="blackHoleTorusCoveringFraction_" source="parameters"/>
    !!]
    self=dustAttenuationAGNTorusSKIRTOR(fileName,url,depthOpticalEdgeOn,exponentDensityRadial,exponentDensityPolar,ratioRadii,galacticInclination_,blackHoleTorusCoveringFraction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticInclination_"           />
    <objectDestructor name="blackHoleTorusCoveringFraction_"/>
    !!]
    return
  end function agnTorusSKIRTORConstructorParameters

  function agnTorusSKIRTORConstructorInternal(fileName,url,depthOpticalEdgeOn,exponentDensityRadial,exponentDensityPolar,ratioRadii,galacticInclination_,blackHoleTorusCoveringFraction_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationAGNTorusSKIRTOR` dust attenuation class. The models
    are read, their transmissions formed, and their covering fractions found.
    !!}
    use :: Dust_AGN_Torus_SKIRTOR, only : skirtorRead , skirtorWeightsInclination
    use :: Error                 , only : Error_Report
    implicit none
    type            (dustAttenuationAGNTorusSKIRTOR     )                                :: self
    type            (varying_string                     ), intent(in   )                 :: fileName                                , url
    double precision                                     , intent(in   )                 :: depthOpticalEdgeOn                      , exponentDensityRadial    , &
         &                                                                                  exponentDensityPolar                    , ratioRadii
    class           (galacticInclinationClass           ), intent(in   ), target         :: galacticInclination_
    class           (blackHoleTorusCoveringFractionClass), intent(in   ), target         :: blackHoleTorusCoveringFraction_
    ! Wavelength (in Å) of the V band, at which the covering fraction is measured.
    double precision                                     , parameter                     :: wavelengthV                    =5500.0d0
    double precision                                     , allocatable, dimension(:,:,:) :: transmissionDirect                      , transmissionScattered
    double precision                                     , allocatable, dimension(:    ) :: angleTorus                              , weights
    double precision                                                                     :: transmissionMeanDirect                  , transmissionMeanScattered, &
         &                                                                                  weightWavelength
    integer                                                                              :: i                                       , j                        , &
         &                                                                                  indexWavelength
    !![
    <constructorAssign variables="fileName, url, depthOpticalEdgeOn, exponentDensityRadial, exponentDensityPolar, ratioRadii, *galacticInclination_, *blackHoleTorusCoveringFraction_"/>
    !!]

    self%inclinationAvailable=self%galacticInclination_%isAvailable()
    call skirtorRead(fileName,url,depthOpticalEdgeOn,exponentDensityRadial,exponentDensityPolar,ratioRadii,angleTorus,self%wavelength,self%inclination,transmissionDirect=transmissionDirect,transmissionScattered=transmissionScattered)
    if (size(angleTorus) < 2) call Error_Report('at least two models, differing in the angle of the torus, are needed'//{introspection:location})
    self%wavelengthLogarithmic=log(self%wavelength)
    weights                   =skirtorWeightsInclination(self%inclination)
    ! Limit the scattered light so that the transmission averaged over orientation never exceeds one.
    do j=1,size(angleTorus)
       do i=1,size(self%wavelength)
          transmissionMeanDirect   =sum(weights*transmissionDirect   (i,:,j))
          transmissionMeanScattered=sum(weights*transmissionScattered(i,:,j))
          if (transmissionMeanDirect+transmissionMeanScattered > 1.0d0)                                          &
               & transmissionScattered(i,:,j)=+transmissionScattered(i,:,j)                                      &
               &                              *max(0.0d0,1.0d0-transmissionMeanDirect)/transmissionMeanScattered
       end do
    end do
    self%transmissionTabulated=transmissionDirect+transmissionScattered
    ! Find the covering fraction of each model from its direct transmission in the V band.
    indexWavelength =max(1,min(size(self%wavelength)-1,count(self%wavelength <= wavelengthV)))
    weightWavelength=+(     log(wavelengthV                        )-self%wavelengthLogarithmic(indexWavelength)) &
         &           /(self%wavelengthLogarithmic(indexWavelength+1)-self%wavelengthLogarithmic(indexWavelength))
    allocate(self%coveringFractions(size(angleTorus)))
    do j=1,size(angleTorus)
       self%coveringFractions(j)=+1.0d0                                                                     &
            &                    -sum(                                                                      &
            &                         weights                                                               &
            &                         *(                                                                    &
            &                           +(1.0d0-weightWavelength)*transmissionDirect(indexWavelength  ,:,j) &
            &                           +       weightWavelength *transmissionDirect(indexWavelength+1,:,j) &
            &                          )                                                                    &
            &                        )
    end do
    if (any(self%coveringFractions(2:size(angleTorus)) <= self%coveringFractions(1:size(angleTorus)-1))) &
         & call Error_Report('the covering fraction of these SKIRTOR models does not increase with the angle of the torus, so a model can not be chosen to match a covering fraction - choose different parameters'//{introspection:location})
    return
  end function agnTorusSKIRTORConstructorInternal

  subroutine agnTorusSKIRTORDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationAGNTorusSKIRTOR` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationAGNTorusSKIRTOR), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticInclination_"           />
    <objectDestructor name="self%blackHoleTorusCoveringFraction_"/>
    !!]
    return
  end subroutine agnTorusSKIRTORDestructor

  function agnTorusSKIRTORTransmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission through the SKIRTOR torus. Only accretion disk emission from the central black hole is
    attenuated; everything else is transmitted in full. The models, and inclination, to interpolate between are found
    once for all parcels.
    !!}
    use :: Dust_Attenuation_Descriptors, only : emissionSourceAccretionDisk
    use :: Error                       , only : Error_Report
    use :: Galactic_Structure_Options  , only : componentTypeBlackHole
    implicit none
    class           (dustAttenuationAGNTorusSKIRTOR), intent(inout)                               :: self
    type            (treeNode                      ), intent(inout), target                       :: node
    type            (emissionDescriptor            ), intent(in   ), dimension(:                ) :: descriptors
    double precision                                , intent(in   ), optional                     :: inclination
    double precision                                               , dimension(size(descriptors)) :: transmission
    double precision                                               , dimension(2,2              ) :: weightsModel
    double precision                                                                              :: inclination_    , coveringFraction     , &
         &                                                                                           weightModel     , weightInclination    , &
         &                                                                                           weightWavelength, wavelengthLogarithmic
    logical                                                                                       :: computed
    integer                                                                                       :: i               , indexModel           , &
         &                                                                                           indexInclination, indexWavelength      , &
         &                                                                                           j               , k

    transmission=1.0d0
    computed    =.false.
    do i=1,size(descriptors)
       if     (                                                             &
            &   descriptors(i)%componentType /= componentTypeBlackHole      &
            &  .or.                                                         &
            &   descriptors(i)%sourceType    /= emissionSourceAccretionDisk &
            & ) cycle
       if (.not.computed) then
          if (present(inclination)) then
             inclination_=inclination
          else if (self%inclinationAvailable) then
             inclination_=self%galacticInclination_%inclination(node)
          else
             inclination_=0.0d0
             call Error_Report('this attenuator depends on orientation, but no inclination is available: either set `galacticInclination` to a class which supplies one, or wrap this attenuator in `inclinationAveraged`'//{introspection:location})
          end if
          ! Find the models bracketing the required covering fraction, and the weight of each.
          coveringFraction =max(                                                                 &
               &                    self%coveringFractions(1                           )       , &
               &                min(                                                             &
               &                    self%coveringFractions(size(self%coveringFractions))       , &
               &                    self%blackHoleTorusCoveringFraction_%coveringFraction(node)  &
               &                   )                                                             &
               &               )
          indexModel       =max(1,min(size(self%coveringFractions)-1,count(self%coveringFractions <= coveringFraction)))
          weightModel      =+(     coveringFraction               -self%coveringFractions(indexModel)) &
               &            /(self%coveringFractions(indexModel+1)-self%coveringFractions(indexModel))
          ! Find the bracketing inclinations, and the weight of each.
          inclination_     =max(self%inclination(1),min(self%inclination(size(self%inclination)),inclination_))
          indexInclination =max(1,min(size(self%inclination)-1,count(self%inclination <= inclination_)))
          weightInclination=+(inclination_                        -self%inclination(indexInclination)) &
               &            /(self%inclination(indexInclination+1)-self%inclination(indexInclination))
          weightsModel(1,1)=(1.0d0-weightInclination)*(1.0d0-weightModel)
          weightsModel(2,1)=       weightInclination *(1.0d0-weightModel)
          weightsModel(1,2)=(1.0d0-weightInclination)*       weightModel
          weightsModel(2,2)=       weightInclination *       weightModel
          computed         =.true.
       end if
       ! Interpolate in the logarithm of wavelength, holding the transmission fixed beyond the tabulated range.
       wavelengthLogarithmic=max(self%wavelengthLogarithmic(1),min(self%wavelengthLogarithmic(size(self%wavelengthLogarithmic)),log(descriptors(i)%wavelength)))
       indexWavelength      =max(1,min(size(self%wavelengthLogarithmic)-1,count(self%wavelengthLogarithmic <= wavelengthLogarithmic)))
       weightWavelength     =+(     wavelengthLogarithmic                   -self%wavelengthLogarithmic(indexWavelength)) &
            &                /(self%wavelengthLogarithmic(indexWavelength+1)-self%wavelengthLogarithmic(indexWavelength))
       transmission(i)=0.0d0
       do k=1,2
          do j=1,2
             transmission(i)=+transmission(i)                                                                                              &
                  &          +weightsModel(j,k)                                                                                            &
                  &          *(                                                                                                            &
                  &            +(1.0d0-weightWavelength)*self%transmissionTabulated(indexWavelength  ,indexInclination+j-1,indexModel+k-1) &
                  &            +       weightWavelength *self%transmissionTabulated(indexWavelength+1,indexInclination+j-1,indexModel+k-1) &
                  &           )
          end do
       end do
    end do
    return
  end function agnTorusSKIRTORTransmission

  logical function agnTorusSKIRTORIsOrientationDependent(self) result(isOrientationDependent)
    !!{RST
    Return true: the transmission through the torus depends on the inclination at which it is seen.
    !!}
    implicit none
    class(dustAttenuationAGNTorusSKIRTOR), intent(inout) :: self
    !$GLC attributes unused :: self

    isOrientationDependent=.true.
    return
  end function agnTorusSKIRTORIsOrientationDependent

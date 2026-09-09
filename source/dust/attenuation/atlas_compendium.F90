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
  Implements a dust attenuation class using the dust compendium of :cite:t:`benson_compendium_2018`.
  !!}

  use :: Galactic_Inclinations         , only : galacticInclinationClass
  use :: ISO_Varying_String            , only : varying_string
  use :: Numerical_Interpolation       , only : interpolator
  use :: Numerical_Interpolation_MultiD, only : interpolatorMultiD

  !![
  <enumeration docformat="rst">
   <name>compendiumSpheroidProfile</name>
   <description>
   Enumerates the spheroid density profiles for which dust compendium tabulations have been computed.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="hernquist"/>
   <entry label="jaffe"    />
  </enumeration>
  !!]

  !![
  <dustAttenuation name="dustAttenuationAtlasCompendium" docformat="rst">
   <description>
   Dust attenuation from the dust compendium of :cite:t:`benson_compendium_2018`: radiative transfer solutions for simple
   galactic geometries, tabulating the fraction of light escaping a galaxy as a function of wavelength, inclination,
   optical depth, and---for the spheroid---the size of the spheroid relative to the disk.

   This is a *non-separable* attenuator, exactly as :galacticus-class:`dustAttenuationAtlasFerrara2000` is: the
   tabulated escape fraction is the result of a radiative transfer calculation through a specified geometry, and is
   not the exponential of a single optical depth times a wavelength-dependent curve, so it takes no
   ``dustExtinctionCurve``. It also means the transmission can slightly exceed unity where scattering redirects more
   light into the line of sight than the dust removes from it.

   The tabulations are published as a set of HDF5 files, one per combination of grain properties and geometry, and
   are described in the :doc:`dust compendium datasets &lt;/manuals/user-guide/data/dust-compendium-datasets&gt;` section
   of the user guide. ``fileName`` selects one. If that file is not already present it is downloaded from ``url``
   and cached under the dynamic datasets path, so a file is fetched once and reused thereafter.

   Two quantities are supplied per galaxy rather than tabulated:

   * The :math:`V`-band optical depth through the center of the galaxy, perpendicular to the plane. Unlike
     :galacticus-class:`dustAttenuationAtlasFerrara2000`, which delegates this to a
     :galacticus-class:`dustAttenuationScreen`, it is computed here from the tabulation's *own* opacity, read from
     the ``opacity`` attribute of the file:

     .. math::

        \tau_\mathrm{V} = \kappa_\mathrm{V} \, f_\mathrm{dust:metals} \, \Sigma_\mathrm{Z},

     with :math:`\kappa_\mathrm{V}` the :math:`V`-band opacity per unit mass of dust used in the radiative transfer
     calculation, :math:`f_\mathrm{dust:metals}` the dust-to-metals ratio (``dustToMetalsRatio``), and
     :math:`\Sigma_\mathrm{Z} = Z M_\mathrm{gas} / 2 \pi r_\mathrm{d}^2` the central surface density of gas-phase
     metals of the disk. Using the opacity which produced the table is what makes the optical depth mean the same
     thing to the model as it does to the tabulation; a screen calibrated to the Milky Way, as
     :galacticus-class:`dustAttenuationScreenSurfaceDensityMetals` is, need not agree with it.

     As in :galacticus-class:`dustAttenuationAtlasFerrara2000`, the optical depth is *always* that of the disk,
     whichever component is being attenuated: the dust lies in the disk, and a spheroid is reddened by the disk's
     dust.

   * The inclination, from a :galacticus-class:`galacticInclinationClass` object, or from the ``inclination``
     argument when one is imposed, as :galacticus-class:`dustAttenuationInclinationAveraged` does. One or the other
     must be available, and an error is reported if neither is.

   * The size of the spheroid, in the units the ``spheroidScaleRadial`` axis is tabulated against. Note that this
     axis does *not* mean the same thing as the corresponding axis of
     :galacticus-class:`dustAttenuationAtlasFerrara2000`: :cite:t:`ferrara_atlas_1999` tabulate against the spheroid
     *effective* radius, whereas the compendium tabulates against the spheroid *scale* radius. ``spheroidProfile``
     names the profile the tabulation was computed for---``hernquist`` for most of them, ``jaffe`` for those
     computed to match :cite:t:`ferrara_atlas_1999`---and the model's spheroid is matched to it on half-mass radius,
     which is the radius that means the same thing whatever profile either side assumes, then converted to that
     profile's scale radius. For a Hernquist profile the half-mass radius is :math:`(1+\sqrt{2})` times the scale
     radius, and for a Jaffe profile the two are equal.

   Beyond the largest tabulated optical depth the transmission is extrapolated as
   :math:`T = \exp(c_0 + c_1 \ln \tau_\mathrm{V})`, using coefficients tabulated alongside the attenuations. Setting
   ``extrapolateOpticalDepth`` to false instead holds the transmission at its value at the tabulation boundary.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationAtlasCompendium
     !!{RST
     A dust attenuation class using the dust compendium of :cite:t:`benson_compendium_2018`.
     !!}
     private
     class           (galacticInclinationClass                ), pointer                         :: galacticInclination_          => null()
     type            (varying_string                          )                                  :: fileName                                , url
     type            (enumerationCompendiumSpheroidProfileType)                                  :: spheroidProfile
     double precision                                                                            :: dustToMetalsRatio                       , opacity                         , &
          &                                                                                         radiusSpheroidHalfMassToScale
     logical                                                                                     :: extrapolateOpticalDepth                 , inclinationAvailable
     ! Grid axes. Wavelengths are in microns and inclinations in degrees, as tabulated; the spheroid axis is the
     ! spheroid scale radius in units of the disk scale length.
     double precision                                          , allocatable, dimension(:      ) :: wavelength                              , inclination                     , &
          &                                                                                         depthOptical                            , radiusSpheroid
     ! Tabulated transmission. The datasets are written (wavelength, inclination, opticalDepth[, radius]) as a
     ! row-major reader sees them, so the Fortran dimensions are the reverse of that.
     double precision                                          , allocatable, dimension(:,:,:  ) :: transmissionDisk
     double precision                                          , allocatable, dimension(:,:,:,:) :: transmissionSpheroid
     ! Coefficients of the high-optical-depth extrapolation, split into their constant and logarithmic terms at read
     ! time so that each is a contiguous array the interpolators can be handed directly.
     double precision                                          , allocatable, dimension(:,:    ) :: extrapolationDiskConstant               , extrapolationDiskLogarithmic
     double precision                                          , allocatable, dimension(:,:,:  ) :: extrapolationSpheroidConstant           , extrapolationSpheroidLogarithmic
     type            (interpolatorMultiD                      )                                  :: interpolatorDisk                        , interpolatorSpheroid            , &
          &                                                                                         interpolatorDiskExtrapolate             , interpolatorSpheroidExtrapolate
     ! The axis interpolators are retained as well as being handed to the multilinear interpolators above. Only the
     ! wavelength changes between the parcels of a galaxy, so bracketing the other axes once per galaxy and reusing
     ! the result saves that work on every parcel.
     type            (interpolator                            )                                  :: interpolatorWavelength                  , interpolatorInclination         , &
          &                                                                                         interpolatorDepthOptical                , interpolatorRadiusSpheroid
   contains
     final     ::                      atlasCompendiumDestructor
     procedure :: transmission      => atlasCompendiumTransmission
     procedure :: request           => atlasCompendiumRequest
     procedure :: supportsComponent => atlasCompendiumSupportsComponent
  end type dustAttenuationAtlasCompendium

  interface dustAttenuationAtlasCompendium
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationAtlasCompendium` dust attenuation class.
     !!}
     module procedure atlasCompendiumConstructorParameters
     module procedure atlasCompendiumConstructorInternal
  end interface dustAttenuationAtlasCompendium

contains

  function atlasCompendiumConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationAtlasCompendium` dust attenuation class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: ISO_Varying_String, only : char          , var_str        , varying_string
    implicit none
    type            (dustAttenuationAtlasCompendium)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (galacticInclinationClass      ), pointer       :: galacticInclination_
    type            (varying_string                )                :: fileName               , url, &
         &                                                             spheroidProfile
    double precision                                                :: dustToMetalsRatio
    logical                                                         :: extrapolateOpticalDepth

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <description>
      The name of the dust compendium tabulation file to use. If no file of this name is present under the
      ``dust/compendium`` directory of the dynamic datasets path it is downloaded from ``url``.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>url</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>
      The URL from which to download ``fileName`` if it is not already present. If ``none``, the file must be
      supplied by other means.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>dustToMetalsRatio</name>
      <defaultValue>0.44d0</defaultValue>
      <defaultSource>Approximately correct for the Milky Way (e.g. :cite:t:`popping_dust_2017`).</defaultSource>
      <description>
      The fraction of the mass of metals which is in dust, used with the opacity of the tabulation to convert a
      surface density of metals into an optical depth.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>extrapolateOpticalDepth</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, transmissions for optical depths beyond the largest tabulated value are extrapolated using the
      coefficients supplied with the tabulation. If false, they are held at the value at the tabulation boundary.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>spheroidProfile</name>
      <defaultValue>var_str('hernquist')</defaultValue>
      <description>
      The spheroid density profile for which the tabulation was computed: ``hernquist``, or ``jaffe`` for
      tabulations computed to match :cite:t:`ferrara_atlas_1999`. This sets how a model galaxy's spheroid is mapped
      onto the ``spheroidScaleRadial`` axis.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="galacticInclination" name="galacticInclination_" source="parameters"/>
    !!]
    self=dustAttenuationAtlasCompendium(fileName,url,dustToMetalsRatio,extrapolateOpticalDepth,enumerationCompendiumSpheroidProfileEncode(char(spheroidProfile),includesPrefix=.false.),galacticInclination_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticInclination_"/>
    !!]
    return
  end function atlasCompendiumConstructorParameters

  function atlasCompendiumConstructorInternal(fileName,url,dustToMetalsRatio,extrapolateOpticalDepth,spheroidProfile,galacticInclination_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationAtlasCompendium` dust attenuation class. The
    tabulation is located---downloading it if necessary---read once, here, and interpolators built over it.
    !!}
    use :: Error             , only : Error_Report
    use :: File_Utilities    , only : Directory_Make      , File_Exists       , File_Lock     , File_Unlock   , &
          &                           lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: Input_Paths       , only : inputPath           , pathTypeDataDynamic
    use :: IO_HDF5           , only : hdf5File
    use :: ISO_Varying_String, only : char                , operator(//)      , operator(==)  , varying_string
    use :: System_Download   , only : download
    use :: Table_Labels      , only : extrapolationTypeFix
    implicit none
    type            (dustAttenuationAtlasCompendium          )                                  :: self
    type            (varying_string                          ), intent(in   )                   :: fileName                    , url
    double precision                                          , intent(in   )                   :: dustToMetalsRatio
    logical                                                   , intent(in   )                   :: extrapolateOpticalDepth
    type            (enumerationCompendiumSpheroidProfileType), intent(in   )                   :: spheroidProfile
    class           (galacticInclinationClass                ), intent(in   ), target           :: galacticInclination_
    type            (hdf5File                                )                                  :: file
    type            (lockDescriptor                          )                                  :: lock
    type            (varying_string                          )                                  :: pathFile                    , pathDirectory
    integer                                                                                     :: status
    double precision                                          , allocatable, dimension(:,:,:  ) :: extrapolationDisk
    double precision                                          , allocatable, dimension(:,:,:,:) :: extrapolationSpheroid
    type            (interpolator                            )             , dimension(3      ) :: interpolatorsDisk           , interpolatorsSpheroidExtrapolate
    type            (interpolator                            )             , dimension(4      ) :: interpolatorsSpheroid
    type            (interpolator                            )             , dimension(2      ) :: interpolatorsDiskExtrapolate
    !![
    <constructorAssign variables="fileName, url, dustToMetalsRatio, extrapolateOpticalDepth, spheroidProfile, *galacticInclination_"/>
    !!]

    ! An inclination must be available, either per galaxy or imposed by an attenuator averaging over orientation.
    self%inclinationAvailable=self%galacticInclination_%isAvailable()
    ! The ratio of half-mass to scale radius for the profile the tabulation was computed with. A model galaxy's
    ! spheroid is matched to the tabulation on half-mass radius, then divided by this to reach the scale radius that
    ! the `spheroidScaleRadial` axis is measured in.
    select case (self%spheroidProfile%ID)
    case (compendiumSpheroidProfileHernquist%ID)
       self%radiusSpheroidHalfMassToScale=1.0d0+sqrt(2.0d0)
    case (compendiumSpheroidProfileJaffe    %ID)
       self%radiusSpheroidHalfMassToScale=1.0d0
    case default
       self%radiusSpheroidHalfMassToScale=0.0d0
       call Error_Report('unrecognized spheroid profile'//{introspection:location})
    end select
    ! Locate the tabulation, downloading it if we do not already have it. The lock makes concurrent processes wait
    ! for the first of them to finish the download rather than each starting one of their own.
    pathDirectory=inputPath(pathTypeDataDynamic)//'dust/compendium'
    pathFile     =pathDirectory//'/'//self%fileName
    if (.not.File_Exists(pathFile)) then
       if (self%url == 'none') call Error_Report('the compendium tabulation `'//self%fileName//'` is not present, and no `url` was given from which to download it'//{introspection:location})
       call Directory_Make(pathDirectory)
       call File_Lock     (char(pathFile),lock,lockIsShared=.false.)
       if (.not.File_Exists(pathFile)) then
          call download(char(self%url),char(pathFile),status=status)
          if (status /= 0 .or. .not.File_Exists(pathFile)) then
             call File_Unlock(lock)
             call Error_Report('unable to download the compendium tabulation from `'//self%url//'`'//{introspection:location})
          end if
       end if
       call File_Unlock(lock)
    end if
    !$ call hdf5Access%set()
    file=hdf5File(char(pathFile),readOnly=.true.)
    if (.not.file%hasAttribute('opacity')) then
       !$ call hdf5Access%unset()
       call Error_Report('`'//self%fileName//'` has no `opacity` attribute, so is not a dust compendium tabulation'//{introspection:location})
    end if
    call file%readAttribute('opacity'                          ,self%opacity              )
    call file%readDataset  ('wavelength'                       ,self%wavelength           )
    call file%readDataset  ('inclination'                      ,self%inclination          )
    call file%readDataset  ('opticalDepth'                     ,self%depthOptical         )
    call file%readDataset  ('spheroidScaleRadial'              ,self%radiusSpheroid       )
    call file%readDataset  ('attenuationDisk'                  ,self%transmissionDisk     )
    call file%readDataset  ('attenuationSpheroid'              ,self%transmissionSpheroid )
    call file%readDataset  ('extrapolationCoefficientsDisk'    ,     extrapolationDisk    )
    call file%readDataset  ('extrapolationCoefficientsSpheroid',     extrapolationSpheroid)
    !$ call hdf5Access%unset()
    ! Check that the tables have the shape the axes imply. A transposed read would otherwise show up much later as
    ! quietly wrong attenuation.
    if (any(shape(self%transmissionDisk    ) /= [size(self%depthOptical  ),size(self%inclination ),size(self%wavelength)                       ])) &
         & call Error_Report('`attenuationDisk` does not have the shape implied by the axes'                  //{introspection:location})
    if (any(shape(self%transmissionSpheroid) /= [size(self%radiusSpheroid),size(self%depthOptical),size(self%inclination),size(self%wavelength)])) &
         & call Error_Report('`attenuationSpheroid` does not have the shape implied by the axes'              //{introspection:location})
    if (any(shape(     extrapolationDisk   ) /= [size(self%inclination   ),size(self%wavelength  ),2                                           ])) &
         & call Error_Report('`extrapolationCoefficientsDisk` does not have the shape implied by the axes'    //{introspection:location})
    if (any(shape(     extrapolationSpheroid) /= [size(self%radiusSpheroid),size(self%inclination ),size(self%wavelength),2                    ])) &
         & call Error_Report('`extrapolationCoefficientsSpheroid` does not have the shape implied by the axes'//{introspection:location})
    ! Split the extrapolation coefficients into their constant and logarithmic terms, so that each is a contiguous
    ! array rather than a strided section of a larger one.
    self%extrapolationDiskConstant       =extrapolationDisk    (:,:  ,1)
    self%extrapolationDiskLogarithmic    =extrapolationDisk    (:,:  ,2)
    self%extrapolationSpheroidConstant   =extrapolationSpheroid(:,:,:,1)
    self%extrapolationSpheroidLogarithmic=extrapolationSpheroid(:,:,:,2)
    ! Build interpolators, ordered from the most rapidly varying dimension of the table outward. Optical depth and
    ! spheroid size are interpolated logarithmically, being tabulated on geometric grids; wavelength likewise. Values
    ! outside the tabulated ranges are held at the boundary.
    interpolatorsDisk                   (1)=interpolator(log(self%depthOptical  ),extrapolationType=extrapolationTypeFix)
    interpolatorsDisk                   (2)=interpolator(    self%inclination    ,extrapolationType=extrapolationTypeFix)
    interpolatorsDisk                   (3)=interpolator(log(self%wavelength    ),extrapolationType=extrapolationTypeFix)
    interpolatorsSpheroid               (1)=interpolator(log(self%radiusSpheroid),extrapolationType=extrapolationTypeFix)
    interpolatorsSpheroid               (2)=interpolatorsDisk    (1)
    interpolatorsSpheroid               (3)=interpolatorsDisk    (2)
    interpolatorsSpheroid               (4)=interpolatorsDisk    (3)
    interpolatorsDiskExtrapolate        (1)=interpolatorsDisk    (2)
    interpolatorsDiskExtrapolate        (2)=interpolatorsDisk    (3)
    interpolatorsSpheroidExtrapolate    (1)=interpolatorsSpheroid(1)
    interpolatorsSpheroidExtrapolate    (2)=interpolatorsDisk    (2)
    interpolatorsSpheroidExtrapolate    (3)=interpolatorsDisk    (3)
    self%interpolatorDisk                  =interpolatorMultiD(interpolatorsDisk               )
    self%interpolatorSpheroid              =interpolatorMultiD(interpolatorsSpheroid           )
    self%interpolatorDiskExtrapolate       =interpolatorMultiD(interpolatorsDiskExtrapolate    )
    self%interpolatorSpheroidExtrapolate   =interpolatorMultiD(interpolatorsSpheroidExtrapolate)
    self%interpolatorDepthOptical          =interpolatorsDisk    (1)
    self%interpolatorInclination           =interpolatorsDisk    (2)
    self%interpolatorWavelength            =interpolatorsDisk    (3)
    self%interpolatorRadiusSpheroid        =interpolatorsSpheroid(1)
    return
  end function atlasCompendiumConstructorInternal

  subroutine atlasCompendiumDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationAtlasCompendium` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationAtlasCompendium), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticInclination_"/>
    !!]
    return
  end subroutine atlasCompendiumDestructor

  double precision function atlasCompendiumDepthOpticalV(self,node) result(depthOpticalV)
    !!{RST
    Return the :math:`V`-band optical depth through the center of the disk, perpendicular to its plane, computed
    from the opacity of the tabulation itself.

    The tabulation records the opacity per unit mass of dust which was used in the radiative transfer calculation
    which produced it, so combining that with a dust-to-metals ratio and the surface density of gas-phase metals
    places a model galaxy on the tabulation using the same definition of optical depth that the tabulation was built
    with.
    !!}
    use :: Galactic_Structure_Options      , only : componentTypeDisk
    use :: Numerical_Constants_Astronomical, only : massSolar        , megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : hecto            , kilo
    implicit none
    class           (dustAttenuationAtlasCompendium), intent(inout)         :: self
    type            (treeNode                      ), intent(inout), target :: node
    double precision                                                        :: massGas             , radius, &
         &                                                                     metallicity         ,         &
         &                                                                     densitySurfaceMetals

    call componentGasProperties(node,componentTypeDisk,massGas,radius,metallicity)
    ! A disk with no gas, or no size, has no dust.
    if (massGas <= 0.0d0 .or. radius <= 0.0d0) then
       depthOpticalV=0.0d0
       return
    end if
    ! Central surface density of metals, in g/cm².
    densitySurfaceMetals=+metallicity      &
         &               *massGas          &
         &               *massSolar        &
         &               *kilo             &
         &               /2.0d0            &
         &               /Pi               &
         &               /(                &
         &                 +radius         &
         &                 *megaParsec     &
         &                 *hecto          &
         &                )**2
    depthOpticalV       =+self%opacity           &
         &               *self%dustToMetalsRatio &
         &               *densitySurfaceMetals
    return
  end function atlasCompendiumDepthOpticalV

  function atlasCompendiumTransmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission of each parcel, interpolated in the compendium tabulation.

    The optical depth and the inclination are properties of the galaxy rather than of a parcel, so both are obtained
    once here and reused across every parcel, as is the spheroid size if any parcel needs it. Whether the optical
    depth lies beyond the tabulation---and so whether the transmission is interpolated or extrapolated---is likewise
    a property of the galaxy, and is decided once.
    !!}
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Error                           , only : Error_Report
    use            :: Galactic_Structure_Options      , only : componentTypeDisk, componentTypeSpheroid
    use            :: Numerical_Constants_Astronomical, only : degreesToRadians
    implicit none
    class           (dustAttenuationAtlasCompendium), intent(inout)                               :: self
    type            (treeNode                      ), intent(inout), target                       :: node
    type            (emissionDescriptor            ), intent(in   ), dimension(:                ) :: descriptors
    double precision                                , intent(in   ), optional                     :: inclination
    double precision                                               , dimension(size(descriptors)) :: transmission
    ! The tabulation is in microns, while parcels report their wavelength in Angstroms.
    double precision                                , parameter                                   :: micronsPerAngstrom        =1.0d-4
    double precision                                                                              :: depthOptical                     , inclination_       , &
         &                                                                                           radiusSpheroid                   , logDepth           , &
         &                                                                                           inclinationDegrees               , coefficientConstant, &
         &                                                                                           coefficientLogarithmic
    logical                                                                                       :: radiusSpheroidComputed           , extrapolating
    integer                                                                                       :: i
    ! Bracketing indices and linear weights, per dimension, in the order the tables are laid out: for the disk
    ! (optical depth, inclination, wavelength), and for the spheroid (spheroid size, optical depth, inclination,
    ! wavelength). The extrapolation coefficients carry no optical depth axis, so drop it. Only the wavelength entry
    ! changes between parcels.
    integer         (c_size_t                      )                , dimension(    3)            :: indicesDisk
    double precision                                                , dimension(0:1,3)            :: weightsDisk
    integer         (c_size_t                      )                , dimension(    4)            :: indicesSpheroid
    double precision                                                , dimension(0:1,4)            :: weightsSpheroid
    integer         (c_size_t                      )                , dimension(    2)            :: indicesDiskExtrapolate
    double precision                                                , dimension(0:1,2)            :: weightsDiskExtrapolate
    integer         (c_size_t                      )                , dimension(    3)            :: indicesSpheroidExtrapolate
    double precision                                                , dimension(0:1,3)            :: weightsSpheroidExtrapolate
    integer         (c_size_t                      )                                              :: indexInclination                 , indexWavelength   , &
         &                                                                                           indexRadiusSpheroid              , indexDepthOptical
    double precision                                                , dimension(0:1  )            :: weightInclination                , weightWavelength  , &
         &                                                                                           weightRadiusSpheroid             , weightDepthOptical

    ! The dust lies in the disk in this model, and a spheroid is reddened by the disk's dust, so the optical depth is
    ! always that of the disk.
    depthOptical=atlasCompendiumDepthOpticalV(self,node)
    if (present(inclination)) then
       inclination_=inclination
    else if (self%inclinationAvailable) then
       inclination_=self%galacticInclination_%inclination(node)
    else
       inclination_=0.0d0
       call Error_Report('this attenuator depends on orientation, but no inclination is available: either set `galacticInclination` to a class which supplies one, or wrap this attenuator in `inclinationAveraged`'//{introspection:location})
    end if
    ! The tabulation is in degrees, and interpolated logarithmically in optical depth. A galaxy with no dust at all
    ! transmits everything, so guard the logarithm rather than letting it overflow.
    inclinationDegrees=+inclination_     &
         &             /degreesToRadians
    if (depthOptical <= 0.0d0) then
       transmission=1.0d0
       return
    end if
    logDepth              =log(depthOptical)
    extrapolating         =  self%extrapolateOpticalDepth                              &
         &                 .and.                                                       &
         &                   depthOptical > self%depthOptical(size(self%depthOptical))
    radiusSpheroidComputed=.false.
    radiusSpheroid        =0.0d0
    ! Bracket the axes which are properties of the galaxy rather than of a parcel, once. Only the set of factors
    ! belonging to the interpolator which will actually be used is maintained.
    call self%interpolatorInclination%linearFactors(inclinationDegrees,indexInclination,weightInclination)
    if (extrapolating) then
       indicesDiskExtrapolate    (  1  )=indexInclination
       weightsDiskExtrapolate    (:,1  )=weightInclination
       indicesSpheroidExtrapolate(  2  )=indexInclination
       weightsSpheroidExtrapolate(:,2  )=weightInclination
    else
       call self%interpolatorDepthOptical%linearFactors(logDepth,indexDepthOptical,weightDepthOptical)
       indicesDisk               (  1  )=indexDepthOptical
       weightsDisk               (:,1  )=weightDepthOptical
       indicesDisk               (  2  )=indexInclination
       weightsDisk               (:,2  )=weightInclination
       indicesSpheroid           (  2:3)=indicesDisk       (  1:2)
       weightsSpheroid           (:,2:3)=weightsDisk       (:,1:2)
    end if
    do i=1,size(descriptors)
       call self%interpolatorWavelength%linearFactors(log(descriptors(i)%wavelength*micronsPerAngstrom),indexWavelength,weightWavelength)
       if (extrapolating) then
          indicesDiskExtrapolate    (  2)=indexWavelength
          weightsDiskExtrapolate    (:,2)=weightWavelength
          indicesSpheroidExtrapolate(  3)=indexWavelength
          weightsSpheroidExtrapolate(:,3)=weightWavelength
       else
          indicesDisk               (  3)=indexWavelength
          weightsDisk               (:,3)=weightWavelength
          indicesSpheroid           (  4)=indexWavelength
          weightsSpheroid           (:,4)=weightWavelength
       end if
       if      (descriptors(i)%componentType == componentTypeDisk    ) then
          if (extrapolating) then
             coefficientConstant   =self%interpolatorDiskExtrapolate%interpolateFactors(self%extrapolationDiskConstant   ,indicesDiskExtrapolate,weightsDiskExtrapolate)
             coefficientLogarithmic=self%interpolatorDiskExtrapolate%interpolateFactors(self%extrapolationDiskLogarithmic,indicesDiskExtrapolate,weightsDiskExtrapolate)
             transmission(i)       =exp(                        &
                  &                     +coefficientConstant    &
                  &                     +coefficientLogarithmic &
                  &                     *logDepth               &
                  &                    )
          else
             transmission(i)       =self%interpolatorDisk%interpolateFactors(self%transmissionDisk,indicesDisk,weightsDisk)
          end if
       else if (descriptors(i)%componentType == componentTypeSpheroid) then
          if (.not.radiusSpheroidComputed) then
             ! Clamp into the tabulated range before taking a logarithm. A galaxy may have no spheroid, or no disk
             ! to measure one against, giving a ratio of zero whose logarithm would trap; and the interpolator holds
             ! values at the boundary in any case, so nothing is lost by clamping here rather than there.
             radiusSpheroid        =max(                                     &
                  &                     +radiusSpheroidRelative(node)        &
                  &                     /self%radiusSpheroidHalfMassToScale, &
                  &                     +minval(self%radiusSpheroid)         &
                  &                    )
             radiusSpheroidComputed=.true.
             call self%interpolatorRadiusSpheroid%linearFactors(log(radiusSpheroid),indexRadiusSpheroid,weightRadiusSpheroid)
             indicesSpheroid           (  1)=indexRadiusSpheroid
             weightsSpheroid           (:,1)=weightRadiusSpheroid
             indicesSpheroidExtrapolate(  1)=indexRadiusSpheroid
             weightsSpheroidExtrapolate(:,1)=weightRadiusSpheroid
          end if
          if (extrapolating) then
             coefficientConstant   =self%interpolatorSpheroidExtrapolate%interpolateFactors(self%extrapolationSpheroidConstant   ,indicesSpheroidExtrapolate,weightsSpheroidExtrapolate)
             coefficientLogarithmic=self%interpolatorSpheroidExtrapolate%interpolateFactors(self%extrapolationSpheroidLogarithmic,indicesSpheroidExtrapolate,weightsSpheroidExtrapolate)
             transmission(i)       =exp(                        &
                  &                     +coefficientConstant    &
                  &                     +coefficientLogarithmic &
                  &                     *logDepth               &
                  &                    )
          else
             transmission(i)       =self%interpolatorSpheroid%interpolateFactors(self%transmissionSpheroid,indicesSpheroid,weightsSpheroid)
          end if
       else
          transmission(i)=1.0d0
          call Error_Report('this tabulation covers only disk and spheroid components'//{introspection:location})
       end if
    end do
    return
  end function atlasCompendiumTransmission

  function atlasCompendiumRequest(self) result(request)
    !!{RST
    Return the decomposition required: the attenuation differs between disk and spheroid, but depends on nothing
    else which varies within a component.
    !!}
    implicit none
    type (decompositionRequest          )                :: request
    class(dustAttenuationAtlasCompendium), intent(inout) :: self
    !$GLC attributes unused :: self

    request%resolveComponents =.true.
    request%resolveMetallicity=.false.
    request%resolveRadius     =.false.
    return
  end function atlasCompendiumRequest

  logical function atlasCompendiumSupportsComponent(self,componentType) result(supportsComponent)
    !!{RST
    Return true only for the disk and spheroid, which are the components this tabulation covers.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid
    implicit none
    class(dustAttenuationAtlasCompendium), intent(inout) :: self
    type (enumerationComponentTypeType  ), intent(in   ) :: componentType
    !$GLC attributes unused :: self

    supportsComponent=  componentType == componentTypeDisk     &
         &            .or.                                     &
         &              componentType == componentTypeSpheroid
    return
  end function atlasCompendiumSupportsComponent

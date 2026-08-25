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
  Implements a dust attenuation class using the atlas of Ferrara et al. (1999).
  !!}

  use :: Galactic_Inclinations         , only : galacticInclinationClass
  use :: Numerical_Interpolation       , only : interpolator
  use :: Numerical_Interpolation_MultiD, only : interpolatorMultiD

  !![
  <enumeration docformat="rst">
   <name>ferrara2000Dust</name>
   <description>
   Enumerates the dust types tabulated in the atlas of :cite:t:`ferrara_atlas_1999`.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="milkyWay"            />
   <entry label="smallMagellanicCloud"/>
  </enumeration>
  !!]

  !![
  <dustAttenuation name="dustAttenuationAtlasFerrara2000" docformat="rst">
   <description>
   Dust attenuation from the radiative transfer atlas of :cite:t:`ferrara_atlas_1999`, which tabulates the fraction
   of light escaping a galaxy as a function of wavelength, inclination, optical depth, and---for the spheroid---the
   size of the spheroid relative to the disk.

   This is a *non-separable* attenuator. The tabulated escape fraction is the result of a radiative transfer
   calculation through a specified geometry, and is not the exponential of a single optical depth times a
   wavelength-dependent curve, so it takes no ``dustExtinctionCurve``---the wavelength dependence is in the table,
   which is the point of tabulating it. It also means the transmission can slightly exceed unity at low optical
   depth and low inclination, where scattering redirects more light into the line of sight than the dust removes
   from it.

   The tabulation is selected by two parameters: ``dust``, the extinction curve of the grains, either ``milkyWay``
   or ``smallMagellanicCloud``; and ``heightRatio``, the ratio of the dust scale height to the stellar scale height,
   one of 0.4, 1.0, or 2.5.

   Two quantities are supplied per galaxy rather than tabulated:

   * The optical depth, obtained from a :galacticus-class:`dustAttenuationScreen` object through its
     ``depthOpticalV`` method. Delegating this keeps the Milky Way calibration of
     :galacticus-class:`dustAttenuationScreenSurfaceDensityMetals` as the single home of that normalization, and
     lets a fixed depth be substituted for testing. It is *always* evaluated for the disk, whichever component is
     being attenuated: in this model the dust lies in the disk and a spheroid is reddened by the disk's dust, so
     asking for a spheroid's own optical depth would compute a surface density of a component holding no dust here.

     The atlas defines its optical depth as that through the center of the galaxy perpendicular to the plane, in the
     :math:`V` band. For an exponential disk :math:`\Sigma_0 = M/2\pi r_\mathrm{d}^2`, which is what
     ``screenSurfaceDensityMetals`` computes, so the definitions agree.

   * The inclination, from a :galacticus-class:`galacticInclinationClass` object, or from the ``inclination``
     argument when one is imposed, as :galacticus-class:`dustAttenuationInclinationAveraged` does. One or the other
     must be available, and an error is reported if neither is.
   </description>
  </dustAttenuation>
  !!]
  type, extends(dustAttenuationClass) :: dustAttenuationAtlasFerrara2000
     !!{RST
     A dust attenuation class using the atlas of :cite:t:`ferrara_atlas_1999`.
     !!}
     private
     class           (galacticInclinationClass       ), pointer                       :: galacticInclination_   => null()
     class           (dustAttenuationScreen          ), pointer                       :: dustAttenuation_       => null()
     type            (enumerationFerrara2000DustType )                                :: dust
     double precision                                                                 :: heightRatio
     logical                                                                          :: inclinationAvailable
     ! Grid axes. Wavelengths are in Angstroms and inclinations in degrees, as tabulated; the spheroid axis is the
     ! spheroid half-mass radius in units of the disk scale length.
     double precision                                 , allocatable, dimension(:      ) :: wavelength           , inclination        , &
          &                                                                                depthOptical         , radiusSpheroid
     ! Tabulated transmission. The datasets are written (wavelength, inclination, opticalDepth[, radius]) as a
     ! row-major reader sees them, so the Fortran dimensions are the reverse of that.
     double precision                                 , allocatable, dimension(:,:,:  ) :: transmissionDisk
     double precision                                 , allocatable, dimension(:,:,:,:) :: transmissionSpheroid
     type            (interpolatorMultiD             )                                  :: interpolatorDisk     , interpolatorSpheroid
   contains
     final     ::                      atlasFerrara2000Destructor
     procedure :: transmission      => atlasFerrara2000Transmission
     procedure :: request           => atlasFerrara2000Request
     procedure :: supportsComponent => atlasFerrara2000SupportsComponent
  end type dustAttenuationAtlasFerrara2000

  interface dustAttenuationAtlasFerrara2000
     !!{RST
     Constructors for the :galacticus-class:`dustAttenuationAtlasFerrara2000` dust attenuation class.
     !!}
     module procedure atlasFerrara2000ConstructorParameters
     module procedure atlasFerrara2000ConstructorInternal
  end interface dustAttenuationAtlasFerrara2000

contains

  function atlasFerrara2000ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustAttenuationAtlasFerrara2000` dust attenuation class which takes a
    parameter set as input.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: ISO_Varying_String, only : char          , var_str       , varying_string
    implicit none
    type            (dustAttenuationAtlasFerrara2000)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (galacticInclinationClass       ), pointer       :: galacticInclination_
    class           (dustAttenuationClass           ), pointer       :: dustAttenuation_
    class           (dustAttenuationScreen          ), pointer       :: dustAttenuationScreen_
    type            (varying_string                 )                :: dust
    double precision                                                 :: heightRatio

    !![
    <inputParameter docformat="rst">
      <name>dust</name>
      <defaultValue>var_str('milkyWay')</defaultValue>
      <description>
      The extinction curve of the dust grains: ``milkyWay`` or ``smallMagellanicCloud``.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>heightRatio</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The ratio of the dust scale height to the stellar scale height. Tabulations exist for 0.4, 1.0, and 2.5.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="galacticInclination" name="galacticInclination_" source="parameters"/>
    <objectBuilder class="dustAttenuation"     name="dustAttenuation_"     source="parameters"/>
    !!]
    ! The optical depth is taken from a screen attenuator, so that its calibration lives in one place.
    select type (dustAttenuation_)
    class is (dustAttenuationScreen)
       dustAttenuationScreen_ => dustAttenuation_
    class default
       call Error_Report('the `dustAttenuation` supplied must be a `screen` class: it is used only as a source of the V-band optical depth'//{introspection:location})
    end select
    self=dustAttenuationAtlasFerrara2000(enumerationFerrara2000DustEncode(char(dust),includesPrefix=.false.),heightRatio,galacticInclination_,dustAttenuationScreen_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticInclination_"/>
    <objectDestructor name="dustAttenuation_"    />
    !!]
    return
  end function atlasFerrara2000ConstructorParameters

  function atlasFerrara2000ConstructorInternal(dust,heightRatio,galacticInclination_,dustAttenuation_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustAttenuationAtlasFerrara2000` dust attenuation class. The
    tabulation is read once, here, and interpolators built over it.
    !!}
    use :: Error                   , only : Error_Report
    use :: HDF5_Access             , only : hdf5Access
    use :: Input_Paths             , only : inputPath                , pathTypeDataStatic
    use :: IO_HDF5                 , only : hdf5File
    use :: ISO_Varying_String      , only : char                     , operator(//)     , varying_string
    use :: Table_Labels            , only : extrapolationTypeFix
    implicit none
    type            (dustAttenuationAtlasFerrara2000)                        :: self
    type            (enumerationFerrara2000DustType ), intent(in   )         :: dust
    double precision                                 , intent(in   )         :: heightRatio
    class           (galacticInclinationClass       ), intent(in   ), target :: galacticInclination_
    class           (dustAttenuationScreen          ), intent(in   ), target :: dustAttenuation_
    type            (hdf5File                       )                        :: file
    type            (varying_string                 )                        :: fileName             , labelDust
    character       (len=3                          )                        :: labelRatio
    type            (interpolator                   ), dimension(3)          :: interpolatorsDisk
    type            (interpolator                   ), dimension(4)          :: interpolatorsSpheroid
    !![
    <constructorAssign variables="dust, heightRatio, *galacticInclination_, *dustAttenuation_"/>
    !!]

    ! An inclination must be available, either per galaxy or imposed by an attenuator averaging over orientation.
    self%inclinationAvailable=self%galacticInclination_%isAvailable()
    ! Select the tabulation.
    select case (self%dust%ID)
    case (ferrara2000DustMilkyWay            %ID)
       labelDust='MilkyWay'
    case (ferrara2000DustSmallMagellanicCloud%ID)
       labelDust='SmallMagellanicCloud'
    case default
       labelDust=''
       call Error_Report('unrecognized dust type'//{introspection:location})
    end select
    if      (abs(self%heightRatio-0.4d0) < 1.0d-3) then
       labelRatio='0.4'
    else if (abs(self%heightRatio-1.0d0) < 1.0d-3) then
       labelRatio='1.0'
    else if (abs(self%heightRatio-2.5d0) < 1.0d-3) then
       labelRatio='2.5'
    else
       labelRatio='   '
       call Error_Report('`heightRatio` must be one of the tabulated values 0.4, 1.0, or 2.5'//{introspection:location})
    end if
    fileName=inputPath(pathTypeDataStatic)//'dust/atlasFerrara2000/attenuations_'//labelDust//'_dustHeightRatio'//labelRatio//'.hdf5'
    !$ call hdf5Access%set()
    file=hdf5File(char(fileName),readOnly=.true.)
    call file%readDataset('wavelength'         ,self%wavelength          )
    call file%readDataset('inclination'        ,self%inclination         )
    call file%readDataset('opticalDepth'       ,self%depthOptical        )
    call file%readDataset('spheroidScaleRadial',self%radiusSpheroid      )
    call file%readDataset('attenuationDisk'    ,self%transmissionDisk    )
    call file%readDataset('attenuationSpheroid',self%transmissionSpheroid)
    !$ call hdf5Access%unset()
    ! Check that the tables have the shape the axes imply. A transposed read would otherwise show up much later as
    ! quietly wrong attenuation.
    if (any(shape(self%transmissionDisk    ) /= [size(self%depthOptical),size(self%inclination),size(self%wavelength)                      ])) &
         & call Error_Report('`attenuationDisk` does not have the shape implied by the axes'    //{introspection:location})
    if (any(shape(self%transmissionSpheroid) /= [size(self%radiusSpheroid),size(self%depthOptical),size(self%inclination),size(self%wavelength)])) &
         & call Error_Report('`attenuationSpheroid` does not have the shape implied by the axes'//{introspection:location})
    ! Build interpolators, ordered from the most rapidly varying dimension of the table outward. Optical depth and
    ! spheroid size are interpolated logarithmically, being tabulated on geometric grids; wavelength likewise. Values
    ! outside the tabulated ranges are held at the boundary.
    interpolatorsDisk    (1)=interpolator(log(self%depthOptical  ),extrapolationType=extrapolationTypeFix)
    interpolatorsDisk    (2)=interpolator(    self%inclination    ,extrapolationType=extrapolationTypeFix)
    interpolatorsDisk    (3)=interpolator(log(self%wavelength    ),extrapolationType=extrapolationTypeFix)
    interpolatorsSpheroid(1)=interpolator(log(self%radiusSpheroid),extrapolationType=extrapolationTypeFix)
    interpolatorsSpheroid(2)=interpolatorsDisk(1)
    interpolatorsSpheroid(3)=interpolatorsDisk(2)
    interpolatorsSpheroid(4)=interpolatorsDisk(3)
    self%interpolatorDisk    =interpolatorMultiD(interpolatorsDisk    )
    self%interpolatorSpheroid=interpolatorMultiD(interpolatorsSpheroid)
    return
  end function atlasFerrara2000ConstructorInternal

  subroutine atlasFerrara2000Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`dustAttenuationAtlasFerrara2000` dust attenuation class.
    !!}
    implicit none
    type(dustAttenuationAtlasFerrara2000), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticInclination_"  />
    <objectDestructor name="self%dustAttenuation_"      />
    !!]
    return
  end subroutine atlasFerrara2000Destructor

  function atlasFerrara2000Transmission(self,node,descriptors,inclination) result(transmission)
    !!{RST
    Return the transmission of each parcel, interpolated in the atlas.

    The optical depth and the inclination are properties of the galaxy rather than of a parcel, so both are obtained
    once here and reused across every parcel, as is the spheroid size if any parcel needs it.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galactic_Structure_Options      , only : componentTypeDisk    , componentTypeSpheroid
    use :: Numerical_Constants_Astronomical, only : degreesToRadians
    implicit none
    class           (dustAttenuationAtlasFerrara2000), intent(inout)                               :: self
    type            (treeNode                       ), intent(inout), target                       :: node
    type            (emissionDescriptor             ), intent(in   ), dimension(:                ) :: descriptors
    double precision                                 , intent(in   ), optional                     :: inclination
    double precision                                                , dimension(size(descriptors)) :: transmission
    double precision                                                                               :: depthOptical      , inclination_, &
         &                                                                                            radiusSpheroid    , logDepth    , &
         &                                                                                            inclinationDegrees
    logical                                                                                        :: radiusSpheroidComputed
    integer                                                                                        :: i

    ! The dust lies in the disk in this model, and a spheroid is reddened by the disk's dust, so the optical depth is
    ! always that of the disk.
    depthOptical=self%dustAttenuation_%depthOpticalV(node,componentTypeDisk)
    if (present(inclination)) then
       inclination_=inclination
    else if (self%inclinationAvailable) then
       inclination_=self%galacticInclination_%inclination(node)
    else
       inclination_=0.0d0
       call Error_Report('this attenuator depends on orientation, but no inclination is available: either set `galacticInclination` to a class which supplies one, or wrap this attenuator in `inclinationAveraged`'//{introspection:location})
    end if
    ! The atlas is tabulated in degrees, and interpolated logarithmically in optical depth. A galaxy with no dust at
    ! all transmits everything, so guard the logarithm rather than letting it overflow.
    inclinationDegrees=+inclination_     &
         &             /degreesToRadians
    if (depthOptical <= 0.0d0) then
       transmission=1.0d0
       return
    end if
    logDepth              =log(depthOptical)
    radiusSpheroidComputed=.false.
    radiusSpheroid        =0.0d0
    do i=1,size(descriptors)
       if      (descriptors(i)%componentType == componentTypeDisk    ) then
          transmission(i)=self%interpolatorDisk    %interpolate(self%transmissionDisk    ,[logDepth,inclinationDegrees,log(descriptors(i)%wavelength)])
       else if (descriptors(i)%componentType == componentTypeSpheroid) then
          if (.not.radiusSpheroidComputed) then
             ! Clamp into the tabulated range before taking a logarithm. A galaxy may have no spheroid, or no disk
             ! to measure one against, giving a ratio of zero whose logarithm would trap; and the interpolator holds
             ! values at the boundary in any case, so nothing is lost by clamping here rather than there.
             radiusSpheroid        =max(atlasFerrara2000RadiusSpheroid(node),minval(self%radiusSpheroid))
             radiusSpheroidComputed=.true.
          end if
          transmission(i)=self%interpolatorSpheroid%interpolate(self%transmissionSpheroid,[log(radiusSpheroid),logDepth,inclinationDegrees,log(descriptors(i)%wavelength)])
       else
          transmission(i)=1.0d0
          call Error_Report('this atlas tabulates only disk and spheroid components'//{introspection:location})
       end if
    end do
    return
  end function atlasFerrara2000Transmission

  double precision function atlasFerrara2000RadiusSpheroid(node) result(radiusSpheroid)
    !!{RST
    Return the size of the spheroid in the units the atlas is tabulated against: its half-mass radius, in units of
    the disk scale length.

    :cite:t:`ferrara_atlas_1999` model spheroids as Jaffe profiles and tabulate against the spheroid effective
    radius. For a Jaffe profile the enclosed mass is :math:`M(r)/M = (r/r_0)/(1+r/r_0)`, so the half-mass radius is
    exactly the scale radius :math:`r_0`; the tabulated axis is therefore a half-mass radius. A model galaxy's
    spheroid will in general follow some other profile, so it is matched on that same physical radius rather than on
    a scale radius, whose meaning differs between profiles.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid
    implicit none
    type            (treeNode             ), intent(inout), target :: node
    class           (nodeComponentDisk    )               , pointer :: disk
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    double precision                                                :: radiusDisk

    disk       => node    %disk  ()
    spheroid   => node    %spheroid()
    radiusDisk =  disk    %radius()
    if (radiusDisk <= 0.0d0) then
       ! With no disk there is no scale to measure the spheroid against, and no dust either, so the value is
       ! immaterial: return the smallest tabulated size.
       radiusSpheroid=0.0d0
       return
    end if
    radiusSpheroid=+spheroid%radius() &
         &         /         radiusDisk
    return
  end function atlasFerrara2000RadiusSpheroid

  function atlasFerrara2000Request(self) result(request)
    !!{RST
    Return the decomposition required: the attenuation differs between disk and spheroid, but depends on nothing
    else which varies within a component.
    !!}
    implicit none
    type (decompositionRequest              )                :: request
    class(dustAttenuationAtlasFerrara2000   ), intent(inout) :: self
    !$GLC attributes unused :: self

    request%resolveComponents =.true.
    request%resolveMetallicity=.false.
    request%resolveRadius     =.false.
    return
  end function atlasFerrara2000Request

  logical function atlasFerrara2000SupportsComponent(self,componentType) result(supportsComponent)
    !!{RST
    Return true only for the disk and spheroid, which are the components this atlas tabulates.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid
    implicit none
    class(dustAttenuationAtlasFerrara2000), intent(inout) :: self
    type (enumerationComponentTypeType   ), intent(in   ) :: componentType
    !$GLC attributes unused :: self

    supportsComponent=  componentType == componentTypeDisk     &
         &            .or.                                     &
         &              componentType == componentTypeSpheroid
    return
  end function atlasFerrara2000SupportsComponent

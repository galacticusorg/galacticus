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
  Implements a spectrum of emission from the dusty torus around an AGN using the CLUMPY library of models.
  !!}

  use :: ISO_Varying_String     , only : varying_string
  use :: Numerical_Interpolation, only : interpolator

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumAGNTorusCLUMPY" docformat="rst">
   <description>
   Emission from the dusty torus around an AGN, taken from the CLUMPY library of models of a torus of optically thick
   clouds of dust :cite:p:`nenkova_agn_2008-1,nenkova_agn_2008`. The clouds are distributed as in
   :galacticus-class:`dustAttenuationAGNTorus`, with a mean number :math:`N(\beta) = N_0 \exp(-\beta^2/\sigma^2)` met
   along a ray at an angle :math:`\beta` from the equatorial plane, where :math:`N_0` is ``[countCloudsEquatorial]`` and
   :math:`\sigma` is ``[widthAngular]``. Each cloud has an optical depth ``[depthOpticalCloudV]`` in the :math:`V` band,
   and the number of clouds per unit length falls as :math:`r^{-q}`, with :math:`q` given by ``[exponentRadial]``, out to
   ``[ratioRadii]`` times the radius at which dust sublimates. The model includes light scattered by the clouds as well
   as their thermal emission, matching :galacticus-class:`dustAttenuationAGNTorus`, which counts scattered light as
   absorbed. This class is intended to be used with that one, so that, with matching ``[countCloudsEquatorial]`` and
   ``[depthOpticalCloudV]``, the torus absorbing the light of the AGN and the torus emitting it are the same.

   The emission of the model, averaged over all directions, sets the shape of the spectrum, which is normalized to the
   luminosity absorbed. The shape is therefore fixed: it does not follow the angular width of the torus chosen for each
   galaxy by :galacticus-class:`dustAttenuationAGNTorus`, nor depend on the inclination at which the torus is seen. The
   parameters need not be tabulated values: the spectrum is interpolated linearly between the tabulated models, in
   :math:`N_0`, :math:`\sigma`, and :math:`q`, and in the logarithms of the optical depth of a cloud and of the ratio of
   radii. It is an error for a parameter to lie outside the range tabulated. Between the tabulated wavelengths the
   spectrum is interpolated as a power law, and it is normalized by integrating that interpolant analytically, so that
   it integrates to exactly the absorbed luminosity. Because it includes scattered light, the spectrum extends down to
   :math:`0.01\,\mu\hbox{m}`; a grid of wavelength for the emission of dust that begins at longer wavelengths misses that
   part of it (about 5% of the luminosity for the default parameters, with the grid beginning at
   :math:`0.1\,\mu\hbox{m}`).

   The temperature of the dust in the torus is set by its distance from the accretion disk, not by the luminosity it
   absorbs, so the mass of dust is not used. Give the torus phase a zero ``fractionsMassPhase`` in
   :galacticus-class:`nodePropertyExtractorSEDDustEmission`.

   The library is read from ``[fileName]``, in the ``dust/agnTorus`` directory of the dynamic datasets path; if it is not
   present it is downloaded from ``[url]``. It is made from the CLUMPY library by ``dataPreparation/convertCLUMPY.py`` of
   the ``galacticusDevTools`` repository (https://github.com/galacticusorg/galacticusDevTools). By default it is
   downloaded from the copy archived on Zenodo (https://doi.org/10.5281/zenodo.23420867), which is redistributed with the
   permission of the authors of CLUMPY.
   </description>
  </dustEmissionSpectrum>
  !!]
  type, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumAGNTorusCLUMPY
     !!{RST
     A spectrum of emission from the dusty torus around an AGN using the CLUMPY library of models.
     !!}
     private
     type            (varying_string)                            :: fileName              , url
     double precision                                            :: countCloudsEquatorial , widthAngular       , &
          &                                                         depthOpticalCloudV    , ratioRadii         , &
          &                                                         exponentRadial
     ! The tabulated wavelengths (in Å) over which the model emits, and the logarithm of the spectrum at each, normalized
     ! to unit luminosity (so in Hz⁻¹).
     double precision                , allocatable, dimension(:) :: wavelength            , spectrumLogarithmic
     type            (interpolator  )                            :: interpolatorWavelength
   contains
     procedure :: luminosity => agnTorusCLUMPYLuminosity
  end type dustEmissionSpectrumAGNTorusCLUMPY

  interface dustEmissionSpectrumAGNTorusCLUMPY
     !!{RST
     Constructors for the :galacticus-class:`dustEmissionSpectrumAGNTorusCLUMPY` dust emission spectrum class.
     !!}
     module procedure agnTorusCLUMPYConstructorParameters
     module procedure agnTorusCLUMPYConstructorInternal
  end interface dustEmissionSpectrumAGNTorusCLUMPY

contains

  function agnTorusCLUMPYConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustEmissionSpectrumAGNTorusCLUMPY` dust emission spectrum class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: ISO_Varying_String, only : var_str
    implicit none
    type            (dustEmissionSpectrumAGNTorusCLUMPY)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    type            (varying_string                    )                :: fileName             , url
    double precision                                                    :: countCloudsEquatorial, widthAngular, &
         &                                                                 depthOpticalCloudV   , ratioRadii  , &
         &                                                                 exponentRadial

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <defaultValue>var_str('CLUMPY_201410_tvavg.hdf5')</defaultValue>
      <description>
      The name of the file holding the CLUMPY tabulation, in the ``dust/agnTorus`` directory of the dynamic datasets
      path.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>url</name>
      <defaultValue>var_str('https://zenodo.org/records/23420867/files/CLUMPY_201410_tvavg.hdf5')</defaultValue>
      <description>
      The URL from which to download ``fileName`` if it is not already present. If ``none``, the file must be supplied
      by other means.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>countCloudsEquatorial</name>
      <defaultValue>5.0d0</defaultValue>
      <defaultSource>The default of :galacticus-class:`dustAttenuationAGNTorus`</defaultSource>
      <description>
      The mean number of clouds, :math:`N_0`, along a ray in the equatorial plane of the torus. Between 1 and 15.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>depthOpticalCloudV</name>
      <defaultValue>30.0d0</defaultValue>
      <defaultSource>The default of :galacticus-class:`dustAttenuationAGNTorus`</defaultSource>
      <description>
      The optical depth in the :math:`V` band of a single cloud. Between 10 and 300.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>widthAngular</name>
      <defaultValue>30.0d0</defaultValue>
      <description>
      The angular width, :math:`\sigma`, in degrees, of the distribution of clouds about the equatorial plane. Between 15
      and 70.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>ratioRadii</name>
      <defaultValue>20.0d0</defaultValue>
      <defaultSource>The default ratio of radii of :galacticus-class:`dustEmissionSpectrumAGNTorusSKIRTOR`</defaultSource>
      <description>
      The ratio, :math:`Y`, of the outer radius of the torus to the radius at which dust sublimates. Between 5 and 100.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentRadial</name>
      <defaultValue>1.0d0</defaultValue>
      <defaultSource>The default radial exponent of :galacticus-class:`dustEmissionSpectrumAGNTorusSKIRTOR`</defaultSource>
      <description>
      The exponent, :math:`q`, of the radial distribution of clouds, whose number per unit length is proportional to
      :math:`r^{-q}`. Between 0 and 3.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustEmissionSpectrumAGNTorusCLUMPY(fileName,url,countCloudsEquatorial,widthAngular,depthOpticalCloudV,ratioRadii,exponentRadial)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function agnTorusCLUMPYConstructorParameters

  function agnTorusCLUMPYConstructorInternal(fileName,url,countCloudsEquatorial,widthAngular,depthOpticalCloudV,ratioRadii,exponentRadial) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustEmissionSpectrumAGNTorusCLUMPY` dust emission spectrum class.
    The models bracketing the requested parameters are read and interpolated between, and the result restricted to the
    range of wavelength over which it emits and normalized to unit luminosity.
    !!}
    use            :: Error                       , only : Error_Report
    use            :: File_Utilities              , only : Directory_Make      , File_Exists        , File_Lock     , File_Remove, &
         &                                                 File_Unlock         , lockDescriptor
    use            :: HDF5                        , only : hsize_t
    use            :: HDF5_Access                 , only : hdf5Access
    use            :: Input_Paths                 , only : inputPath           , pathTypeDataDynamic
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: IO_HDF5                     , only : hdf5File
    use            :: ISO_Varying_String          , only : char                , operator(//)       , operator(==)  , var_str
    use            :: Numerical_Constants_Physical, only : speedLight
    use            :: Numerical_Constants_Units   , only : metersToAngstroms   , micronsToAngstroms
    use            :: System_Download             , only : download
    use            :: Table_Labels                , only : extrapolationTypeFix
    implicit none
    type            (dustEmissionSpectrumAGNTorusCLUMPY)                                          :: self
    type            (varying_string                    ), intent(in   )                           :: fileName                  , url
    double precision                                    , intent(in   )                           :: countCloudsEquatorial     , widthAngular         , &
         &                                                                                           depthOpticalCloudV        , ratioRadii           , &
         &                                                                                           exponentRadial
    double precision                                    , allocatable  , dimension(:          )   :: wavelength                , luminosityLogarithmic, &
         &                                                                                           countCloudsEquatorialModel, widthAngularModel    , &
         &                                                                                           depthOpticalCloudVModel   , ratioRadiiModel      , &
         &                                                                                           exponentRadialModel       , emission
    double precision                                    , allocatable  , dimension(:,:,:,:,:,:)   :: emissionModels
    logical                                             , allocatable  , dimension(:          )   :: isEmitting
    integer         (hsize_t                           )               , dimension(6          )   :: readBegin                 , readCount
    double precision                                                   , dimension(5          )   :: fraction
    type            (lockDescriptor                    )                                          :: lock
    type            (varying_string                    )                                          :: pathFile                  , pathDirectory
    integer         (c_size_t                          )                                          :: sizeFile
    double precision                                                                              :: normalization             , weight
    integer                                                                                       :: status                    , iCloud               , &
         &                                                                                           iWidth                    , iDepth               , &
         &                                                                                           iRatio                    , iExponent
    !![
    <constructorAssign variables="fileName, url, countCloudsEquatorial, widthAngular, depthOpticalCloudV, ratioRadii, exponentRadial"/>
    !!]

    ! Locate the tabulation, downloading it if we do not already have it. The lock makes concurrent processes wait for
    ! the first of them to finish the download rather than each starting one of their own.
    pathDirectory=inputPath(pathTypeDataDynamic)//'dust/agnTorus'
    pathFile     =pathDirectory//'/'//fileName
    if (.not.File_Exists(pathFile)) then
       if (url == 'none') call Error_Report(                                                                                                                     &
            &                               var_str('the CLUMPY tabulation `')//pathFile//'` is not present, and no `url` was given from which to download it'// &
            &                               ' - it can be made from the CLUMPY library with `dataPreparation/convertCLUMPY.py` of galacticusDevTools'         // &
            &                               {introspection:location}                                                                                             &
            &                              )
       call Directory_Make(pathDirectory                          )
       call File_Lock     (pathFile     ,lock,lockIsShared=.false.)
       if (.not.File_Exists(pathFile)) then
          call download(char(url),char(pathFile),retries=3,retryWait=10,status=status)
          ! Judge the download by what actually landed on disk, not by the reported status alone.
          sizeFile=0_c_size_t
          if (File_Exists(pathFile)) inquire(file=char(pathFile),size=sizeFile)
          if (status /= 0 .or. sizeFile <= 0_c_size_t) then
             if (File_Exists(pathFile)) call File_Remove(pathFile)
             call File_Unlock(lock)
             call Error_Report('unable to download the CLUMPY tabulation from `'//url//'`'//{introspection:location})
          end if
       end if
       call File_Unlock(lock)
    end if
    ! Read the axes of the tabulation, find the models which bracket the requested parameters, and read only those. The
    ! tabulation is indexed (wavelength, optical depth, width, exponent, ratio of radii, number of clouds).
    !$ call hdf5Access%set()
    hdf5ReadScope : block
      type(hdf5File) :: file
      file=hdf5File(char(pathFile),readOnly=.true.)
      call file%readDataset('wavelength'           ,wavelength                )
      call file%readDataset('countCloudsEquatorial',countCloudsEquatorialModel)
      call file%readDataset('widthAngular'         ,widthAngularModel         )
      call file%readDataset('depthOpticalCloudV'   ,depthOpticalCloudVModel   )
      call file%readDataset('ratioRadii'           ,ratioRadiiModel           )
      call file%readDataset('exponentRadial'       ,exponentRadialModel       )
      call bracket('countCloudsEquatorial',countCloudsEquatorialModel,countCloudsEquatorial,.false.,iCloud   ,fraction(5))
      call bracket('ratioRadii'           ,ratioRadiiModel           ,ratioRadii           ,.true. ,iRatio   ,fraction(4))
      call bracket('exponentRadial'       ,exponentRadialModel       ,exponentRadial       ,.false.,iExponent,fraction(3))
      call bracket('widthAngular'         ,widthAngularModel         ,widthAngular         ,.false.,iWidth   ,fraction(2))
      call bracket('depthOpticalCloudV'   ,depthOpticalCloudVModel   ,depthOpticalCloudV   ,.true. ,iDepth   ,fraction(1))
      readBegin=int([1               ,iDepth,iWidth,iExponent,iRatio,iCloud],kind=hsize_t)
      readCount=int([size(wavelength),2     ,2     ,2        ,2     ,2     ],kind=hsize_t)
      call file%readDataset('emissionDustMean',emissionModels,readBegin=readBegin,readCount=readCount)
    end block hdf5ReadScope
    !$ call hdf5Access%unset()
    wavelength=wavelength*micronsToAngstroms
    ! Interpolate multilinearly between the bracketing models. The emission is λL_λ per unit luminosity of the AGN, so this
    ! is a weighted mixture of the models.
    allocate(emission(size(wavelength)))
    emission=0.0d0
    do iCloud=0,1
       do iRatio=0,1
          do iExponent=0,1
             do iWidth=0,1
                do iDepth=0,1
                   weight  =+merge(fraction(5),1.0d0-fraction(5),iCloud    == 1) &
                        &   *merge(fraction(4),1.0d0-fraction(4),iRatio    == 1) &
                        &   *merge(fraction(3),1.0d0-fraction(3),iExponent == 1) &
                        &   *merge(fraction(2),1.0d0-fraction(2),iWidth    == 1) &
                        &   *merge(fraction(1),1.0d0-fraction(1),iDepth    == 1)
                   emission=+emission                                                                 &
                        &   +weight*emissionModels(:,iDepth+1,iWidth+1,iExponent+1,iRatio+1,iCloud+1)
                end do
             end do
          end do
       end do
    end do
    ! Restrict to the range over which the model emits.
    isEmitting=emission > 0.0d0
    if (count(isEmitting) < 2) call Error_Report('the model emits at fewer than two tabulated wavelengths'//{introspection:location})
    if (any(.not.isEmitting(findloc(isEmitting,.true.,dim=1):findloc(isEmitting,.true.,dim=1,back=.true.)))) &
         & call Error_Report('the model does not emit over a contiguous range of wavelength'//{introspection:location})
    self%wavelength      =pack(wavelength   ,isEmitting)
    luminosityLogarithmic=pack(log(emission),isEmitting)
    ! Convert from λL_λ to L_ν, which is proportional to λ λL_λ, and normalize to unit luminosity. The spectrum is
    ! interpolated as a power law between the tabulated wavelengths, and normalized by integrating that interpolant.
    normalization           =dustEmissionIntegralPowerLaw(self%wavelength,luminosityLogarithmic)
    self%spectrumLogarithmic=+luminosityLogarithmic                             &
         &                   +log(self%wavelength/speedLight/metersToAngstroms) &
         &                   -log(normalization)
    self%interpolatorWavelength=interpolator(log(self%wavelength),extrapolationType=extrapolationTypeFix)
    return

  contains

    subroutine bracket(name,values,value,logarithmic,index,fraction)
      !!{RST
      Find the ``index`` of the tabulated ``values`` (in increasing order) such that ``value`` lies between it and the
      next, and the ``fraction`` of the way from the one to the other at which it lies, measured in the logarithm if
      ``logarithmic`` is true.
      !!}
      implicit none
      character       (len=*), intent(in   )               :: name
      double precision       , intent(in   ), dimension(:) :: values
      double precision       , intent(in   )               :: value
      logical                , intent(in   )               :: logarithmic
      integer                , intent(  out)               :: index
      double precision       , intent(  out)               :: fraction
      ! Tolerance on a value lying at the end of the tabulated range.
      double precision       , parameter                   :: tolerance   =1.0d-6
      character       (len=16)                             :: labelMinimum       , labelMaximum

      if     (                                                                             &
           &   value < values(1           )-tolerance*max(1.0d0,abs(values(1           ))) &
           &  .or.                                                                         &
           &   value > values(size(values))+tolerance*max(1.0d0,abs(values(size(values)))) &
           & ) then
         write (labelMinimum,'(g0.4)') values(1           )
         write (labelMaximum,'(g0.4)') values(size(values))
         !$ call hdf5Access%unset()
         call Error_Report(                                                                                                        &
              &            var_str('`')//name//'` lies outside the range tabulated by CLUMPY, from '//trim(labelMinimum)//' to '// &
              &            trim(labelMaximum)                                                                                   // &
              &            {introspection:location}                                                                                &
              &           )
      end if
      index=1
      do while (index < size(values)-1 .and. value > values(index+1))
         index=index+1
      end do
      if (logarithmic) then
         fraction=log(value/values(index))/log(values(index+1)/values(index))
      else
         fraction=   (value-values(index))/   (values(index+1)-values(index))
      end if
      fraction=min(max(fraction,0.0d0),1.0d0)
      return
    end subroutine bracket

  end function agnTorusCLUMPYConstructorInternal

  function agnTorusCLUMPYLuminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Return the luminosity per unit frequency emitted by the torus, normalized to the luminosity it absorbs.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (dustEmissionSpectrumAGNTorusCLUMPY), intent(inout)                               :: self
    double precision                                    , intent(in   ), dimension(:                ) :: wavelengths
    double precision                                    , intent(in   )                               :: luminosityAbsorbed, massDust, &
         &                                                                                               time
    double precision                                                   , dimension(size(wavelengths)) :: luminosity
    double precision                                                   , dimension(0:1              ) :: weights
    integer         (c_size_t                          )                                              :: j
    integer                                                                                           :: i
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
  end function agnTorusCLUMPYLuminosity

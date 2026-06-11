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

  !!{RST
  Implements a class for intergalactic background light.
  !!}

  use            :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator

  !![
  <radiationField name="radiationFieldIntergalacticBackgroundFile" docformat="rst">
   <description>
   A radiation field class for intergalactic background light with properties read from file. The flux is determined by linearly interpolating to the required time and wavelength. The XML or HDF5 file to read is specified by ``[fileName]``. An example of the required file structure for XML files is:

   .. code-block:: none

      &lt;spectrum&gt;
        &lt;URL&gt;http://adsabs.harvard.edu/abs/1996ApJ...461...20H&lt;/URL&gt;
        &lt;description&gt;Cosmic background radiation spectrum from quasars alone.&lt;/description&gt;
        &lt;reference&gt;Haardt, F. &amp; Madau, P. 1996, ApJ, 461, 20&lt;/reference&gt;
        &lt;source&gt;Francesco Haardt on Aug 6 2005, via Cloudy 08.00&lt;/source&gt;
        &lt;wavelengths&gt;
          &lt;datum&gt;0.0002481&lt;/datum&gt;
          &lt;datum&gt;0.001489&lt;/datum&gt;
          .
          .
          .
          &lt;units&gt;Angstroms&lt;/units&gt;
        &lt;/wavelengths&gt;
        &lt;spectra&gt;
          &lt;datum&gt;7.039E-49&lt;/datum&gt;
          &lt;datum&gt;8.379E-48&lt;/datum&gt;
          &lt;datum&gt;1.875E-39&lt;/datum&gt;
          &lt;datum&gt;7.583E-38&lt;/datum&gt;
          .
          .
          .
          &lt;redshift&gt;0&lt;/redshift&gt;
          &lt;units&gt;erg cm^-2 s^-1 Hz^-1 sr^-1&lt;/units&gt;
        &lt;/spectra&gt;
      &lt;/spectrum&gt;

   The optional ``URL``, ``description``, ``reference`` and ``source`` elements can be used to give the provenance of the data. The ``wavelengths`` element should contain a set of ``datum`` elements each containing a wavelength (in increasing order) at which the spectrum will be tabulated. Wavelengths must be given in Angstroms. Multiple ``spectra`` elements can be given, each specifying the spectrum at a redshift as given in the ``redshift`` element. Each ``spectra`` element must contain an array of ``datum`` elements that gives the spectrum at each wavelength listed in the ``wavelength`` element. Spectra must be in units of erg cm\ :math:`^{-2}` s\ :math:`^{-1}` Hz\ :math:`^{-1}` sr\ :math:`^{-1}`.
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </radiationField>
  !!]
  type, extends(radiationFieldIntergalacticBackground) :: radiationFieldIntergalacticBackgroundFile
     !!{RST
     A radiation field class for intergalactic background light with properties read from file.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer                     :: cosmologyFunctions_     => null()
     type            (varying_string         )                              :: fileName
     double precision                                                       :: time_
     integer                                                                :: spectraTimesCount                , spectraWavelengthsCount
     double precision                         , allocatable, dimension(:  ) :: spectraTimes                     , spectraWavelengths
     double precision                         , allocatable, dimension(:,:) :: spectra
     type            (interpolator           )                              :: interpolatorWavelengths          , interpolatorTimes
     integer         (c_size_t               )                              :: iTime
     double precision                         , dimension(0:1)              :: hTime
   contains
     final     ::                      intergalacticBackgroundFileDestructor
     procedure :: flux              => intergalacticBackgroundFileFlux
     procedure :: time              => intergalacticBackgroundFileTime
     procedure :: timeSet           => intergalacticBackgroundFileTimeSet
     procedure :: timeDependentOnly => intergalacticBackgroundFileTimeDependentOnly
  end type radiationFieldIntergalacticBackgroundFile

  interface radiationFieldIntergalacticBackgroundFile
     !!{RST
     Constructors for the ``radiationFieldIntergalacticBackgroundFile`` radiation field class.
     !!}
     module procedure intergalacticBackgroundFileConstructorParameters
     module procedure intergalacticBackgroundFileConstructorInternal
  end interface radiationFieldIntergalacticBackgroundFile

  ! Current file format version for intergalactic background radiation files.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function intergalacticBackgroundFileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``radiationFieldIntergalacticBackgroundFile`` radiation field class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (radiationFieldIntergalacticBackgroundFile)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_
    type (varying_string                           )                :: fileName

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <description>
      The name of the file from which to read intergalactic background light properties.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=radiationFieldIntergalacticBackgroundFile(fileName,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function intergalacticBackgroundFileConstructorParameters

  function intergalacticBackgroundFileConstructorInternal(fileName,cosmologyFunctions_) result(self)
    !!{RST
    Internal constructor for the ``radiationFieldIntergalacticBackgroundFile`` radiation field class.
    !!}
    use :: Array_Utilities, only : Array_Is_Monotonic               , Array_Reverse        , directionIncreasing
    use :: FoX_DOM        , only : destroy                          , extractDataContent   , node
    use :: Error          , only : Error_Report
    use :: HDF5_Access    , only : hdf5Access
    use :: IO_HDF5        , only : hdf5Object
    use :: IO_XML         , only : XML_Array_Read                   , XML_Array_Read_Static, XML_Count_Elements_By_Tag_Name, XML_Get_Elements_By_Tag_Name, &
         &                         XML_Get_First_Element_By_Tag_Name, XML_Parse            , xmlNodeList
    use :: Table_Labels   , only : extrapolationTypeZero
    implicit none
    type            (radiationFieldIntergalacticBackgroundFile)                                :: self
    type            (varying_string                           ), intent(in   )                 :: fileName
    class           (cosmologyFunctionsClass                  ), intent(in   ), target         :: cosmologyFunctions_
    type            (node                                     ), pointer                       :: doc                , datum     , &
         &                                                                                        spectrum           , wavelength
    type            (xmlNodeList                              ), allocatable  , dimension(:  ) :: spectraList
    double precision                                           , allocatable  , dimension(:,:) :: spectraTmp
    integer                                                                                    :: fileFormatVersion  , iSpectrum , &
         &                                                                                        status             , jSpectrum
    logical                                                                                    :: timesIncreasing
    type            (hdf5Object                               )                                :: file
    !![
    <constructorAssign variables="fileName, *cosmologyFunctions_"/>
    !!]

    ! Determine if we are reading and XML or HDF5 file.
    if (extract(fileName,len(fileName)-3,len(fileName)) == ".xml") then
       ! XML file
       !$omp critical (FoX_DOM_Access)
       doc => XML_Parse(char(fileName),iostat=status)
       if (status /= 0) call Error_Report('Unable to find or parse data file'//{introspection:location})
       ! Check the file format version of the file.
       datum => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
       call extractDataContent(datum,fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Error_Report('file format version is out of date'//{introspection:location})
       ! Get a list of all spectra.
       call XML_Get_Elements_By_Tag_Name(doc,"spectra",spectraList)
       ! Get the wavelengths.
       wavelength => XML_Get_First_Element_By_Tag_Name(doc,"wavelengths")
       call XML_Array_Read(wavelength,"datum",self%spectraWavelengths)
       self%spectraWavelengthsCount=size(self%spectraWavelengths)
       ! Allocate array for spectra.
       self%spectraTimesCount=size(spectraList)
       allocate(self%spectra     (self%spectraWavelengthsCount,self%spectraTimesCount))
       allocate(self%spectraTimes(self%spectraTimesCount))
       ! Read times.
       do iSpectrum=1,self%spectraTimesCount
          ! Get the data.
          spectrum => spectraList(iSpectrum-1)%element
          ! Extract the redshift.
          call extractDataContent(XML_Get_First_Element_By_Tag_Name(spectrum,"redshift"),self%spectraTimes(iSpectrum))
          ! Convert redshift to a time.
          self%spectraTimes(iSpectrum)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%spectraTimes(iSpectrum)))
       end do
       ! Check if the times are monotonically ordered.
       if (.not.Array_Is_Monotonic(self%spectraTimes)) call Error_Report('spectra must be monotonically ordered in time'//{introspection:location})
       timesIncreasing=Array_Is_Monotonic(self%spectraTimes,direction=directionIncreasing)
       ! Reverse times if necessary.
       if (.not.timesIncreasing) self%spectraTimes=Array_Reverse(self%spectraTimes)
       ! Read spectra into arrays.
       do iSpectrum=1,self%spectraTimesCount
          ! Determine where to store this spectrum, depending on whether the times were stored in increasing or decreasing order.
          if (timesIncreasing) then
             jSpectrum=                        +iSpectrum
          else
             jSpectrum=self%spectraTimesCount+1-iSpectrum
          end if
          ! Get the data.
          spectrum => spectraList(iSpectrum-1)%element
          ! Check that we have the correct number of data.
          if (XML_Count_Elements_By_Tag_Name(spectrum,"datum") /= self%spectraWavelengthsCount) call Error_Report('all spectra must contain the same number of wavelengths'//{introspection:location})
          ! Extract the data.
          call XML_Array_Read_Static(spectrum,"datum",self%spectra(:,jSpectrum))
       end do
       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
    else if (extract(fileName,len(fileName)-4,len(fileName)) == ".hdf5") then
       ! HDF5 file.
       !$ call hdf5Access%set()
       file=hdf5Object(char(self%fileName),readOnly=.true.)
       ! Check the file format version of the file.
       call file%readAttribute('fileFormat',fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Error_Report('file format version is out of date'//{introspection:location})
       ! Get the wavelengths.
       call file%readDataset('wavelengths',self%spectraWavelengths)
       self%spectraWavelengthsCount=size(self%spectraWavelengths)
       ! Read redshifts.
       call file%readDataset('redshifts',self%spectraTimes)
       self%spectraTimesCount=size(self%spectraTimes)
       do iSpectrum=1,self%spectraTimesCount
          ! Convert redshift to a time.
          self%spectraTimes(iSpectrum)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%spectraTimes(iSpectrum)))
       end do
       ! Check if the times are monotonically ordered.
       if (.not.Array_Is_Monotonic(self%spectraTimes)) call Error_Report('spectra must be monotonically ordered in time'//{introspection:location})
       timesIncreasing=Array_Is_Monotonic(self%spectraTimes,direction=directionIncreasing)
       ! Reverse times if necessary.
       if (.not.timesIncreasing) self%spectraTimes=Array_Reverse(self%spectraTimes)
       ! Read spectra.
       call file%readDataset('spectra',self%spectra)
       if (.not.timesIncreasing) then
          spectraTmp=self%spectra
          do iSpectrum=1,self%spectraTimesCount
             self%spectra(:,iSpectrum)=spectraTmp(:,self%spectraTimesCount+1-iSpectrum)
          end do
       end if
       !$ call hdf5Access%unset()
    else
       call Error_Report('unrecognized file format'//{introspection:location})
    end if
    self%interpolatorWavelengths=interpolator(self%spectraWavelengths,extrapolationType=extrapolationTypeZero)
    self%interpolatorTimes      =interpolator(self%spectraTimes      ,extrapolationType=extrapolationTypeZero)
    return
  end function intergalacticBackgroundFileConstructorInternal

  subroutine intergalacticBackgroundFileDestructor(self)
    !!{RST
    Destructor for the ``radiationFieldIntergalacticBackgroundFile`` radiation field class.
    !!}
    implicit none
    type(radiationFieldIntergalacticBackgroundFile), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine intergalacticBackgroundFileDestructor

  double precision function intergalacticBackgroundFileFlux(self,wavelength,node)
    !!{RST
    Return the flux in the intergalactic background radiation field.
    !!}
    implicit none
    class           (radiationFieldIntergalacticBackgroundFile), intent(inout)  :: self
    double precision                                           , intent(in   )  :: wavelength
    type            (treeNode                                 ), intent(inout)  :: node
    double precision                                           , dimension(0:1) :: hWavelength
    integer         (c_size_t                                 )                 :: iWavelength, jTime, &
         &                                                                         jWavelength
    !$GLC attributes unused :: node

    intergalacticBackgroundFileFlux=0.0d0
    ! Return if out of range.
    if     (                                                                    &
         &   self%time_ < self%spectraTimes      (1                           ) &
         &  .or.                                                                &
         &   self%time_ > self%spectraTimes      (self%spectraTimesCount      ) &
         &  .or.                                                                &
         &   wavelength < self%spectraWavelengths(1                           ) &
         &  .or.                                                                &
         &   wavelength > self%spectraWavelengths(self%spectraWavelengthsCount) &
         & ) return
    ! Find interpolation in the array of wavelengths.
    call self%interpolatorWavelengths%linearFactors(wavelength,iWavelength,hWavelength)
    do jTime=0,1
       do jWavelength=0,1
          if (self%iTime+jTime > 0)                                                                          &
               & intergalacticBackgroundFileFlux=+intergalacticBackgroundFileFlux                            &
               &                                 +self%hTime      (                        jTime           ) &
               &                                 *     hWavelength(jWavelength                             ) &
               &                                 *self%spectra    (jWavelength+iWavelength,jTime+self%iTime)
       end do
    end do
    return
  end function intergalacticBackgroundFileFlux

  double precision function intergalacticBackgroundFileTime(self)
    !!{RST
    Set the time of the intergalactic background radiation field.
    !!}
    implicit none
    class(radiationFieldIntergalacticBackgroundFile), intent(inout) :: self

    intergalacticBackgroundFileTime=self%time_
    return
  end function intergalacticBackgroundFileTime

  subroutine intergalacticBackgroundFileTimeSet(self,time)
    !!{RST
    Set the time of the intergalactic background radiation field.
    !!}
    implicit none
    class           (radiationFieldIntergalacticBackgroundFile), intent(inout) :: self
    double precision                                           , intent(in   ) :: time

    self%time_=time
    ! Find interpolating factors in the array of times.
    call self%interpolatorTimes%linearFactors(time,self%iTime,self%hTime)
    return
  end subroutine intergalacticBackgroundFileTimeSet

  logical function intergalacticBackgroundFileTimeDependentOnly(self)
    !!{RST
    Return true as this radiation field depends on time only.
    !!}
    implicit none
    class(radiationFieldIntergalacticBackgroundFile), intent(inout) :: self
    !$GLC attributes unused :: self

    intergalacticBackgroundFileTimeDependentOnly=.true.
    return
  end function intergalacticBackgroundFileTimeDependentOnly

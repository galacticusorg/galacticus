!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !+    Contributions to this file made by: Alex Merson.

  !!{
  Implements a file-based stellar population spectra class.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <deepCopyActions class="spectralTable">
   <spectralTable>
    <methodCall method="interpolatorsDeepCopy"/>
   </spectralTable>
  </deepCopyActions>
  !!]

  type spectralTable
     !!{
     Structure to hold tabulated stellar population data.
     !!}
     ! The spectra tables.
     integer                                                       :: agesCount             , metallicityCount       , &
          &                                                           wavelengthsCount
     double precision              , allocatable, dimension(:    ) :: ages                  , metallicities          , &
          &                                                           wavelengths
     double precision              , allocatable, dimension(:,:,:) :: table
     type            (interpolator)                                :: interpolatorAge       , interpolatorMetallicity, &
          &                                                           interpolatorWavelength
   contains
     !![
     <methods>
       <method description="Perform deep copy actions on interpolators." method="interpolatorsDeepCopy" />
     </methods>
     !!]
     procedure :: interpolatorsDeepCopy => spectralTableInterpolatorsDeepCopy
  end type spectralTable

  !![
  <stellarPopulationSpectra name="stellarPopulationSpectraFile">
   <description>
    A stellar population spectra class which computes spectra via interpolation in a tabulation read from file. This should be
    an HDF5 file with the following structure:
    \begin{verbatim}
     ages                     Dataset {ageCount}
     metallicities            Dataset {metallicityCount}
     spectra                  Dataset {metallicityCount, ageCount, metallicityCount}
     wavelengths              Dataset {wavelengthCount}
    \end{verbatim}
    where the datasets contain the tabulated ages (in Gyr), metallicities (logarithmic, relative to Solar), wavelengths (in
    \AA) and spectra (in $L_\odot$ Hz$^{-1}$).
  
    Scripts to convert the data provided by \cite{maraston_evolutionary_2005} and \cite{bruzual_stellar_2003} into \glc's
    format are provided in the {\normalfont \ttfamily scripts/ssps} folder.
   </description>
   <stateStorable>
    <exclude variables="spectra, forceZeroMetallicity, fileName, fileRead"/>
   </stateStorable>
   <runTimeFileDependencies paths="fileName"/>
  </stellarPopulationSpectra>
  !!]
  type, extends(stellarPopulationSpectraClass) :: stellarPopulationSpectraFile
     !!{
     A stellar population spectra class which interpolates spectra given in a file.
     !!}
     private
     type   (spectralTable ) :: spectra
     logical                 :: forceZeroMetallicity, fileRead
     type   (varying_string) :: fileName
   contains
     !![
     <methods>
       <method description="Read the named stellar population spectra file." method="readFile" />
     </methods>
     !!]
     procedure :: readFile           => fileReadFile
     procedure :: luminosity         => fileLuminosity
     procedure :: tabulation         => fileTabulation
     procedure :: wavelengths        => fileWavelengths
     procedure :: wavelengthInterval => fileWavelengthInterval
  end type stellarPopulationSpectraFile

  interface stellarPopulationSpectraFile
     !!{
     Constructors for the file stellar spectra class.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface stellarPopulationSpectraFile

  ! The current file format version.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function fileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the file stellar spectra class which takes a parameter set as input.
    !!}
    implicit none
    type   (stellarPopulationSpectraFile)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    type   (varying_string              )                :: fileName
    logical                                              :: forceZeroMetallicity

    !![
    <inputParameter>
      <name>forceZeroMetallicity</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>Force the use of zero metallicity (or lowest metallicity available) for all stellar populations.</description>
    </inputParameter>
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file from which to read spectra.</description>
    </inputParameter>
    !!]
    self=stellarPopulationSpectraFile(forceZeroMetallicity,char(fileName))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fileConstructorParameters

  function fileConstructorInternal(forceZeroMetallicity,fileName) result(self)
    !!{
    Internal constructor for the file stellar spectra class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type     (stellarPopulationSpectraFile)                :: self
    logical                                , intent(in   ) :: forceZeroMetallicity
    character(len=*                       ), intent(in   ) :: fileName
    !![
    <constructorAssign variables="forceZeroMetallicity, fileName"/>
    !!]

    self%fileRead=.false.
    return
  end function fileConstructorInternal

  double precision function fileLuminosity(self,abundancesStellar,age,wavelength,status)
    !!{
    Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population with composition {\normalfont \ttfamily abundances}, of the
    given {\normalfont \ttfamily age} (in Gyr) and the specified {\normalfont \ttfamily wavelength} (in Angstroms). This is found by interpolating in tabulated
    spectra.
    !!}
    use            :: Abundances_Structure, only : Abundances_Get_Metallicity           , abundances            , logMetallicityZero, max, &
          &                                        metallicityTypeLogarithmicByMassSolar
    use            :: Error               , only : Error_Report                         , errorStatusInputDomain, errorStatusSuccess
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    implicit none
    class           (stellarPopulationSpectraFile), intent(inout)            :: self
    type            (abundances                  ), intent(in   )            :: abundancesStellar
    double precision                              , intent(in   )            :: age                        , wavelength
    integer                                       , intent(  out) , optional :: status
    double precision                              , dimension(0:1)           :: hAge                       , hMetallicity, &
         &                                                                      hWavelength
    double precision                              , parameter                :: metallicityTolerance=0.01d0
    integer         (c_size_t                    )                           :: iAge                       , iMetallicity, &
         &                                                                      iWavelength
    integer                                                                  :: jWavelength                , jMetallicity, &
         &                                                                      jAge
    double precision                                                         :: metallicity
    type            (varying_string              ), save                     :: message
    !$omp threadprivate(message)
    character       (len=12                      )                           :: metallicityLabel           , label

    ! Ensure that the file has been read.
    call self%readFile()
    ! Set status code if present.
    if (present(status)) status=errorStatusSuccess
    ! Check for out of range conditions.
    if (age > self%spectra%ages(self%spectra%agesCount)) then
       write (label,'(e12.4)') age
       message='age ['//trim(label)//'] exceeds the maximum tabulated ['
       write (label,'(e12.4)') self%spectra%ages(self%spectra%agesCount)
       message=message//trim(label)//']'
       call Error_Report(message//{introspection:location})
    end if
    if (self%forceZeroMetallicity) then
       metallicity=logMetallicityZero
    else
       metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    end if
    if (metallicity > self%spectra%metallicities(self%spectra%metallicityCount)+metallicityTolerance) then
       if (present(status)) then
          metallicity=self%spectra%metallicities(self%spectra%metallicityCount)
          status     =errorStatusInputDomain
       else
          write (metallicityLabel,'(f12.6)') metallicity
          message='metallicity ['//trim(adjustl(metallicityLabel))//'] exceeds the maximum tabulated ['
          write (metallicityLabel,'(f12.6)') self%spectra%metallicities(self%spectra%metallicityCount)
          message=message//trim(adjustl(metallicityLabel))//']'
          call Error_Report(message//{introspection:location})
       end if
    end if
    ! Assume zero flux outside of the tabulated wavelength range.
    if     (                                                                      &
         &   wavelength < self%spectra%wavelengths(                            1) &
         & .or.                                                                   &
         &   wavelength > self%spectra%wavelengths(self%spectra%wavelengthsCount) &
         & ) then
       fileLuminosity=0.0d0
       return
    end if
    ! Get the interpolating factors.    
    call self%spectra%interpolatorAge       %linearFactors(age       ,iAge       ,hAge       )
    call self%spectra%interpolatorWavelength%linearFactors(wavelength,iWavelength,hWavelength)
    if      (                                                                         &
         &    metallicity == logMetallicityZero                                       &
         &   .or.                                                                     &
         &    metallicity <  self%spectra%metallicities(                           1) &
         &  ) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (                                                                         &
         &   metallicity >  self%spectra%metallicities(self%spectra%metallicityCount) &
         &  ) then
       iMetallicity=self%spectra%metallicityCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       call self%spectra%interpolatorMetallicity%linearFactors(metallicity,iMetallicity,hMetallicity)
    end if
    ! Do the interpolation.
    fileLuminosity=0.0d0
    do jAge=0,1
       do jWavelength=0,1
          do jMetallicity=0,1
             fileLuminosity=+fileLuminosity                                &
                  &         +self                                          &
                  &              %spectra                                  &
                  &                      %table(                           &
                  &                             iWavelength +jWavelength , &
                  &                             iAge        +jAge        , &
                  &                             iMetallicity+jMetallicity  &
                  &                            )                           &
                  &         *hAge              (             jAge        ) &
                  &         *hWavelength       (             jWavelength ) &
                  &         *hMetallicity      (             jMetallicity)
           end do
        end do
     end do
     ! Prevent interpolation from returning negative fluxes.
     fileLuminosity=max(fileLuminosity,0.0d0)
    return
  end function fileLuminosity

  subroutine fileReadFile(self)
    !!{
    Read a file of simple stellar population spectra.
    !!}
    use :: Error       , only : Error_Report
    use :: HDF5_Access , only : hdf5Access
    use :: IO_HDF5     , only : hdf5Object
    use :: Table_Labels, only : extrapolationTypeExtrapolate, extrapolationTypeZero, extrapolationTypeFix
    implicit none
    class  (stellarPopulationSpectraFile), intent(inout) :: self
    integer                                              :: fileFormatVersion
    type   (hdf5Object                  )                :: spectraFile

    ! Read file if necessary.
    if (.not.self%fileRead) then
       !$ call hdf5Access%set()
       ! Open the HDF5 file.
       spectraFile=hdf5Object(self%fileName,readOnly=.true.)
       ! Check that this file has the correct format.
       call spectraFile%readAttribute('fileFormat',fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Error_Report('format of stellar tracks file is out of date'//{introspection:location})
       ! Read the wavelengths array.
       call spectraFile%readDataset('wavelengths'                 ,self%spectra%wavelengths  )
       self%spectra%wavelengthsCount=size(self%spectra%wavelengths  )
       ! Read the ages array.
       call spectraFile%readDataset('ages'                        ,self%spectra%ages         )
       self%spectra%agesCount       =size(self%spectra%ages         )
       ! Read the metallicity array.
       call spectraFile%readDataset('metallicities'               ,self%spectra%metallicities)
       self%spectra%metallicityCount=size(self%spectra%metallicities)
       ! Read the spectra.
       call spectraFile%readDataset('spectra'                     ,self%spectra%table        )
       !$ call hdf5Access%unset()
       self%fileRead=.true.
       ! Build interpolators.
       self%spectra%interpolatorAge        =interpolator(self%spectra%ages         ,extrapolationType=extrapolationTypeExtrapolate)
       self%spectra%interpolatorWavelength =interpolator(self%spectra%wavelengths  ,extrapolationType=extrapolationTypeZero       )
       self%spectra%interpolatorMetallicity=interpolator(self%spectra%metallicities,extrapolationType=extrapolationTypeFix        )
    end if
    return
  end subroutine fileReadFile

  subroutine fileTabulation(self,agesCount,metallicitiesCount,ages,metallicity)
    !!{
    Return a tabulation of ages and metallicities at which stellar spectra should be tabulated.
    !!}
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    class           (stellarPopulationSpectraFile)                           , intent(inout) :: self
    integer                                                                  , intent(  out) :: agesCount, metallicitiesCount
    double precision                              , allocatable, dimension(:), intent(  out) :: ages     , metallicity

    ! Ensure that the file has been read.
    call self%readFile()
    ! Return the relevant data.
    agesCount         =self%spectra%agesCount
    metallicitiesCount=self%spectra%metallicityCount
    allocate(ages       (agesCount         ))
    allocate(metallicity(metallicitiesCount))
    ages              =         self%spectra%ages
    metallicity       =(10.0d0**self%spectra%metallicities)*metallicitySolar
    return
  end subroutine fileTabulation

  subroutine fileWavelengths(self,wavelengthsCount,wavelengths)
    !!{
    Return a tabulation of wavelengths at which stellar spectra should be tabulated.
    !!}
    implicit none
    class           (stellarPopulationSpectraFile)                           , intent(inout) :: self
    integer                                                                  , intent(  out) :: wavelengthsCount
    double precision                              , allocatable, dimension(:), intent(  out) :: wavelengths

    ! Ensure that the file has been read.
    call self%readFile()
    ! Return the relevant data.
    wavelengthsCount=self%spectra  %wavelengthsCount
    allocate(wavelengths(wavelengthsCount))
    wavelengths     =self%spectra  %wavelengths
    return
  end subroutine fileWavelengths

  double precision function fileWavelengthInterval(self,wavelength)
    !!{
    Return a tabulation of wavelengths at which stellar spectra should be tabulated.
    !!}
    implicit none
    class           (stellarPopulationSpectraFile)                           , intent(inout) :: self
    double precision                                                         , intent(in   ) :: wavelength
    integer                                                                                  :: wavelengthsCount
    double precision                              , allocatable, dimension(:)                :: wavelengths     , wavelengthDifference

    ! Ensure that the file is read.
    call self%readFile()
    ! Return the relevant data.
    wavelengthsCount=self%spectra  %wavelengthsCount
    allocate(wavelengths(wavelengthsCount))
    wavelengths     =self%spectra  %wavelengths
    ! Check if wavelength inside range
    if     (                                            &
         &   wavelength < wavelengths(               1) &
         &  .or.                                        &
         &   wavelength > wavelengths(wavelengthsCount) &
         & ) then
       fileWavelengthInterval=-999.9d0
    else
       ! Compute difference in wavelength at position of interest
       allocate(wavelengthDifference(wavelengthsCount))
       wavelengthDifference  =+       wavelengths &
            &                 -       wavelength
       fileWavelengthInterval=+minval(wavelengths,dim=1,mask=wavelengthDifference >  0.0d0) &
            &                 -maxval(wavelengths,dim=1,mask=wavelengthDifference <= 0.0d0)
    end if
    deallocate(wavelengths)
    if(allocated(wavelengthDifference)) deallocate(wavelengthDifference)
    return
  end function fileWavelengthInterval

  subroutine spectralTableInterpolatorsDeepCopy(self)
    !!{
    Perform deep copy actions on interpolators.
    !!}
    implicit none
    class(spectralTable), intent(inout) :: self

    call self%interpolatorAge        %GSLReallocate()
    call self%interpolatorMetallicity%GSLReallocate()
    call self%interpolatorWavelength %GSLReallocate()
    return
  end subroutine spectralTableInterpolatorsDeepCopy

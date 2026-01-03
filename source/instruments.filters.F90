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

!+    Contributions to this file made by:  Alex Merson.

!!{
Contains a module which implements calculations of filter response curves.
!!}

module Instruments_Filters
  !!{
  Implements calculations of filter response curves.
  !!}
  use :: ISO_Varying_String , only : varying_string
  use :: Locks              , only : ompReadWriteLock
  implicit none
  private
  public :: Filter_Get_Index, Filter_Response            , Filter_Extent     , Filter_Vega_Offset      , &
       &    Filter_Name     , Filter_Wavelength_Effective, Filters_Initialize, Filter_Response_Function

  type filterType
     !!{
     A structure which holds filter response curves.
     !!}
     integer                                                     :: nPoints
     double precision                , allocatable, dimension(:) :: response           , wavelength
     type            (varying_string)                            :: name
     logical                                                     :: vegaOffsetAvailable, wavelengthEffectiveAvailable
     double precision                                            :: vegaOffset         , wavelengthEffective
  end type filterType

  ! Array to hold filter data.
  type   (filterType      ), allocatable, dimension(:) :: filterResponses
  integer                                              :: countFilterResponses    =0

  ! Read/write lock used to control access to the filter responses.
  type   (ompReadWriteLock)                            :: lock
  logical                                              :: lockInitialized         =.false.

  ! Options controlling output.
  logical                                              :: filtersConstructedOutput
  
contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Filters_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Filters_Initialize(parameters)
    !!{
    Initialize the module.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    
    !![
    <inputParameter>
      <name>filtersConstructedOutput</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, any filters constructed internally will be written to file.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    return
  end subroutine Filters_Initialize
  
  integer function Filter_Get_Index(filterName)
    !!{
    Return the index for the specified filter, loading that filter if necessary.
    !!}
    use :: ISO_Varying_String, only : operator(==)
    implicit none
    type   (varying_string), intent(in   ) :: filterName
    integer                                :: iFilter
    logical                                :: isLocked

    ! Initialize lock if necessary.
    if (.not.lockInitialized) then
       !$omp critical (instrumentsFiltersLock)
       if (.not.lockInitialized) then
          lock           =ompReadWriteLock()
          lockInitialized=.true.
       end if
       !$omp end critical (instrumentsFiltersLock)
    end if
    ! See if we already have this filter loaded. If not, load it.
    isLocked=lock%owned()
    if (.not.isLocked) call lock%setRead()
    if (countFilterResponses == 0) then
       call Filter_Response_Load(filterName)
       Filter_Get_Index=countFilterResponses
    else
       Filter_Get_Index=0
       do iFilter=1,countFilterResponses
         if (filterResponses(iFilter)%name == filterName) then
             Filter_Get_Index=iFilter
             exit
          end if
       end do
       if (Filter_Get_Index == 0) then
          call Filter_Response_Load(filterName)
          Filter_Get_Index=countFilterResponses
       end if
    end if
    if (.not.isLocked) call lock%unsetRead()
    return
  end function Filter_Get_Index

  function Filter_Name(filterIndex)
    !!{
    Return the name of the specified filter.
    !!}
    implicit none
    type   (varying_string)                :: Filter_Name
    integer                , intent(in   ) :: filterIndex
    logical                                :: isLocked

    isLocked=lock%owned()
    if (.not.isLocked) call lock%setRead()
    Filter_Name=filterResponses(filterIndex)%name
    if (.not.isLocked) call lock%unsetRead()
    return
  end function Filter_Name

  function Filter_Extent(filterIndex)
    !!{
    Return an array containing the minimum and maximum wavelengths tabulated for this specified filter.
    !!}
    implicit none
    double precision, dimension(2)  :: Filter_Extent
    integer         , intent(in   ) :: filterIndex
    logical                         :: isLocked

    isLocked=lock%owned()
    if (.not.isLocked) call lock%setRead()
    Filter_Extent(1)=filterResponses(filterIndex)%wavelength(1)
    Filter_Extent(2)=filterResponses(filterIndex)%wavelength(filterResponses(filterIndex)%nPoints)
    if (.not.isLocked) call lock%unsetRead()
    return
  end function Filter_Extent

  subroutine Filter_Response_Load(filterName)
    !!{
    Load a filter response curve.
    !!}
    use :: File_Utilities           , only : File_Exists            , Directory_Make
    use :: FoX_DOM                  , only : DOMException           , Node               , destroy           , getExceptionCode, &
         &                                   inException
    use :: FoX_WXML                 , only : xml_AddCharacters      , xml_Close          , xml_EndElement    , xml_NewElement  , &
          &                                  xml_OpenFile           , xmlf_t
    use :: Error                    , only : Error_Report
    use :: Input_Paths              , only : inputPath              , pathTypeDataDynamic, pathTypeDataStatic
    use :: HII_Region_Emission_Lines, only : emissionLineWavelength
    use :: IO_XML                   , only : XML_Array_Read         , XML_Parse
    use :: ISO_Varying_String       , only : assignment(=)          , char               , operator(//)      , operator(==)    , &
         &                                   extract
    use :: String_Handling          , only : String_Split_Words     , operator(//)
    implicit none
    type            (varying_string), intent(in   )               :: filterName
    type            (Node          ), pointer                     :: doc
    type            (filterType    ), allocatable  , dimension(:) :: filterResponsesTemporary
    type            (varying_string)               , dimension(4) :: specialFilterWords
    double precision                , parameter                   :: cutOffResolution        =1.0d4
    integer                                                       :: filterIndex                   , ioErr                     , &
         &                                                           i
    type            (varying_string)                              :: filterFileName                , filterDescription         , &
         &                                                           errorMessage
    type            (DOMException  )                              :: exception
    type            (xmlf_t        )                              :: filterDoc
    logical                                                       :: parseSuccess                  , filterConstructed
    character       (len=64        )                              :: word                          , label
    double precision                                              :: centralWavelength             , resolution                , &
         &                                                           filterWidth

    call lock%setWrite(haveReadLock=.true.)
    ! Allocate space for this filter.
    if (countFilterResponses == 0) then
       allocate(filterResponses(1))
    else if (countFilterResponses == size(filterResponses)) then
       call Move_Alloc(filterResponses,filterResponsesTemporary)
       allocate(filterResponses(2*size(filterResponsesTemporary)))
       filterResponses(1:size(filterResponsesTemporary))=filterResponsesTemporary
       deallocate(filterResponsesTemporary)       
    end if
    ! Index in array to load into.
    countFilterResponses=countFilterResponses+1
    filterIndex         =countFilterResponses
    ! Store the name of the filter.
    filterResponses(filterIndex)%name=filterName
    ! Check for special filters.
    if (extract(filterName,1,22)=="fixedResolutionTopHat_") then
       ! Construct a top-hat filter. Extract central wavelength and resolution.
       call String_Split_Words(specialFilterWords,char(filterName),separator="_")
       word=char(specialFilterWords(2))
       read (word,*) centralWavelength
       word=char(specialFilterWords(3))
       read (word,*) resolution
       filterResponses(filterIndex)%nPoints            =4
       allocate(filterResponses(filterIndex)%wavelength(4))
       allocate(filterResponses(filterIndex)%response  (4))
       filterResponses(filterIndex)%wavelength         =                                                                  &
            & [                                                                                                           &
            &  centralWavelength*(sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution/(1.0d0+1.0d0/cutOffResolution), &
            &  centralWavelength*(sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution                               , &
            &  centralWavelength*(sqrt(4.0d0*resolution**2+1.0d0)+1.0d0)/2.0d0/resolution                               , &
            &  centralWavelength*(sqrt(4.0d0*resolution**2+1.0d0)+1.0d0)/2.0d0/resolution*(1.0d0+1.0d0/cutOffResolution)  &
            & ]
       filterResponses(filterIndex)%response           =                                                                  &
            & [                                                                                                           &
            &  0.0d0                                                                                                    , &
            &  1.0d0                                                                                                    , &
            &  1.0d0                                                                                                    , &
            &  0.0d0                                                                                                      &
            & ]
       filterConstructed=.true.
       write (label,'(f16.8,":",f16.12)') centralWavelength,resolution
       filterDescription="Top-hat filter; wavelength:resolution = "//trim(label)
    else if (extract(filterName,1,25)=="adaptiveResolutionTopHat_") then
       ! Construct an SED top-hat filter. Extract central wavelength and top hat width.
       call String_Split_Words(specialFilterWords,char(filterName),separator="_")
       word=char(specialFilterWords(2))
       read (word,*) centralWavelength
       word=char(specialFilterWords(3))
       read (word,*) filterWidth
       filterResponses(filterIndex)%nPoints            =4
       allocate(filterResponses(filterIndex)%wavelength(4))
       allocate(filterResponses(filterIndex)%response  (4))
       filterResponses(filterIndex)%wavelength         =                                                              &
            & [                                                                                                       &
            &  centralWavelength-filterWidth/2.0d0-filterWidth/100.0d0                                              , &
            &  centralWavelength-filterWidth/2.0d0                                                                  , &
            &  centralWavelength+filterWidth/2.0d0                                                                  , &
            &  centralWavelength+filterWidth/2.0d0+filterWidth/100.0d0                                                &
            & ]
       filterResponses(filterIndex)%response           =                                                              &
            & [                                                                                                       &
            &  0.0d0                                                                                                , &
            &  1.0d0                                                                                                , &
            &  1.0d0                                                                                                , &
            &  0.0d0                                                                                                  &
            & ]
       filterConstructed=.true.
       write (label,'(f16.8,":",f16.8)') centralWavelength,filterWidth
       filterDescription="Top-hat filter; wavelength:width = "//trim(label)
    else if (extract(filterName,1,21)=="emissionLineContinuum") then
       ! Construct a top-hat filter for calculating equivalent width of specified emission line.
       ! From filter name extract line name (to determine central wavelength) and resolution.
       call String_Split_Words(specialFilterWords,char(filterName),separator="_")
       ! Check whether filter is centered on emission line or is offset and extract wavelength
       ! and resolution accordingly
       word=char(specialFilterWords(1))
       if (word=="emissionLineContinuumCentral") then
          word=char(specialFilterWords(2))
          centralWavelength=emissionLineWavelength(word)
          word=char(specialFilterWords(3))
          read (word,*) resolution
       else if (word=="emissionLineContinuumBracketed") then
          word=char(specialFilterWords(3))
          read (word,*) centralWavelength
          word=char(specialFilterWords(4))
          read (word,*) resolution
       end if
       filterResponses(filterIndex)%nPoints            =4
       allocate(filterResponses(filterIndex)%wavelength(4))
       allocate(filterResponses(filterIndex)%response  (4))
       filterResponses(filterIndex)%wavelength         =                                                                &
            & [                                                                                                         &
            &  centralWavelength*(sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution/(1.0+1.0d0/cutOffResolution), &
            &  centralWavelength*(sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution                             , &
            &  centralWavelength*(sqrt(4.0d0*resolution**2+1.0d0)+1.0d0)/2.0d0/resolution                             , &
            &  centralWavelength*(sqrt(4.0d0*resolution**2+1.0d0)+1.0d0)/2.0d0/resolution*(1.0+1.0d0/cutOffResolution)  &
            & ]
       filterResponses(filterIndex)%response           =                                                                &
            & [                                                                                                         &
            &  0.0d0                                                                                                  , &
            &  1.0d0                                                                                                  , &
            &  1.0d0                                                                                                  , &
            &  0.0d0                                                                                                    &
            & ]
       filterConstructed=.true.
       write (label,'(f16.8,":",f16.12)') centralWavelength,resolution
       filterDescription="Emission line continuum filter; wavelength:resolution = "//trim(label)
    else
       ! Construct a file name for the filter.
       filterFileName=char(inputPath(pathTypeDataStatic))//'filters/'//filterName//'.xml'
       if (.not.File_Exists(filterFileName)) then
          filterFileName=char(inputPath(pathTypeDataDynamic))//'filters/'//filterName//'.xml'
          if (.not.File_Exists(filterFileName)) call Error_Report('filter file for filter "'//filterName//'" can not be found'//{introspection:location})
       end if
       ! Parse the XML file.
       !$omp critical (FoX_DOM_Access)
       doc => XML_Parse(char(filterFileName),iostat=ioErr,ex=exception)
       parseSuccess=.true.
       errorMessage=''
       if (inException(exception)) then
          parseSuccess=.false.
          errorMessage=errorMessage//char(10)//'exception raised [code='//getExceptionCode(exception)//']'
       end if
       if (ioErr /= 0) then
          parseSuccess=.false.
          errorMessage=errorMessage//char(10)//'I/O error [code='//ioErr//']'
       end if
       if (.not.parseSuccess)                                                      &
            & call Error_Report(                                                   &
            &                   'unable to read or parse filter response file: '// &
            &                   char(filterFileName)                            // &
            &                        errorMessage                               // &
            &                   {introspection:location}                           &
            &                  )
       ! Extract wavelengths and filter response.
       call XML_Array_Read(doc,"datum",filterResponses(filterIndex)%wavelength,filterResponses(filterIndex)%response)
       filterResponses(filterIndex)%nPoints=size(filterResponses(filterIndex)%wavelength)
       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
       filterConstructed=.false.
    end if
    ! Mark effective wavelength and Vega offset as unavailable.
    filterResponses(filterIndex)%wavelengthEffectiveAvailable=.false.
    filterResponses(filterIndex)%         vegaOffsetAvailable=.false.
    call lock%unsetWrite(haveReadLock=.true.)
    ! Write out the filter file if necessary.
    call Directory_Make(inputPath(pathTypeDataDynamic)//'filters')
    filterFileName=inputPath(pathTypeDataDynamic)//'filters/'//filterName//'.xml'
    if (filtersConstructedOutput .and. filterConstructed .and. .not.File_Exists(filterFileName)) then
       call xml_OpenFile     (char(filterFileName),filterDoc)
       call xml_NewElement   (filterDoc,"filter"     )
       call xml_NewElement   (filterDoc,"description")
       call xml_AddCharacters(filterDoc,char(filterDescription        ))
       call xml_EndElement   (filterDoc,"description")
       call xml_NewElement   (filterDoc,"name"       )
       call xml_AddCharacters(filterDoc,char(filterName               ))
       call xml_EndElement   (filterDoc,"name"       )
       call xml_NewElement   (filterDoc,"origin"     )
       call xml_AddCharacters(filterDoc,     "Automatically generated" )
       call xml_EndElement   (filterDoc,"origin"     )
       call xml_NewElement   (filterDoc,"response"   )
       do i=1,size(filterResponses(filterIndex)%response)
          call xml_NewElement   (filterDoc,"datum")
          write (label,'(f16.8,1x,f16.12)') filterResponses(filterIndex)%wavelength(i),filterResponses(filterIndex)%response(i)
          call xml_AddCharacters(filterDoc,trim(adjustl(label)))
          call xml_EndElement   (filterDoc,"datum")
       end do
       call xml_EndElement   (filterDoc,"response"           )
       call xml_NewElement   (filterDoc,"effectiveWavelength")
       write (label,'(f16.8)') Filter_Wavelength_Effective(filterIndex)
       call xml_AddCharacters(filterDoc,trim(adjustl(label)))
       call xml_EndElement   (filterDoc,"effectiveWavelength")
       call xml_NewElement   (filterDoc,"vegaOffset"         )
       write (label,'(f16.12)') Filter_Vega_Offset        (filterIndex)
       call xml_AddCharacters(filterDoc,trim(adjustl(label)))
       call xml_EndElement   (filterDoc,"vegaOffset"         )
       call xml_EndElement   (filterDoc,"filter"             )
       call xml_Close        (filterDoc                      )
    end if
    return
  end subroutine Filter_Response_Load

  !![
  <outputFileClose function="Filters_Output"/>
  !!]
  subroutine Filters_Output()
    !!{
    Output accumulated filter data to file.
    !!}
    use :: Output_HDF5              , only : outputFile
    use :: HDF5_Access              , only : hdf5Access
    use :: IO_HDF5                  , only : hdf5Object
    use :: Numerical_Constants_Units, only : metersToAngstroms
    implicit none
    type            (hdf5Object) :: filtersGroup       , dataset
    integer                      :: i
    double precision             :: wavelengthEffective
    
    if (countFilterResponses == 0) return
    ! Ensure effective wavelengths are set.
    do i=1,countFilterResponses
       wavelengthEffective=Filter_Wavelength_Effective(i)
    end do
    !$ call hdf5Access%set()
    filtersGroup=outputFile%openGroup('Filters','Properties of filters used.')
    call filtersGroup%writeDataset(filterResponses(1:countFilterResponses)%name               ,'name'               ,'Filter name.'                                               )
    call filtersGroup%writeDataset(filterResponses(1:countFilterResponses)%wavelengthEffective,'wavelengthEffective','Effective wavelength of filter [Å].',datasetReturned=dataset)
    call dataset%writeAttribute("Angstroms [Å]"        ,"units"    )
    call dataset%writeAttribute(1.0d0/metersToAngstroms,"unitsInSI")
    call dataset     %close()
    call filtersGroup%close()
    !$ call hdf5Access%unset()
    return
  end subroutine Filters_Output
  
  function Filter_Response_Function(filterIndex) result(interpolator_)
    !!{
    Return the filter response function (as an interpolator) as a function of wavelength (specified in Angstroms). Note that we follow the
    convention of \cite{hogg_k_2002} and assume that the filter response gives the fraction of incident photons received by the
    detector at a given wavelength, multiplied by the relative photon response (which will be 1 for a photon-counting detector
    such as a CCD, or proportional to the photon energy for a bolometer/calorimeter type detector.
    !!}
    use :: Numerical_Interpolation, only : interpolator
    implicit none
    type   (interpolator), pointer       :: interpolator_
    integer              , intent(in   ) :: filterIndex
    logical                              :: isLocked

    allocate(interpolator_)
    isLocked=lock%owned()
    if (.not.isLocked) call lock%setRead()
    interpolator_=interpolator(filterResponses(filterIndex)%wavelength,filterResponses(filterIndex)%response)
    if (.not.isLocked) call lock%unsetRead()
    return
  end function Filter_Response_Function

  double precision function Filter_Response(filterIndex,wavelength)
    !!{
    Return the filter response function at the given {\normalfont \ttfamily wavelength} (specified in Angstroms).
    !!}
    use :: Numerical_Interpolation, only : interpolator
    implicit none
    integer                       , intent(in   ) :: filterIndex
    double precision              , intent(in   ) :: wavelength
    type            (interpolator), pointer       :: interpolator_

    interpolator_   => Filter_Response_Function(filterIndex)
    Filter_Response =  interpolator_%interpolate(wavelength)
    deallocate(interpolator_)
    return
  end function Filter_Response

  double precision function Filter_Wavelength_Effective(filterIndex)
    !!{
    Return the effective wavelength for the specified filter.
    !!}
    implicit none
    integer, intent(in   ) :: filterIndex
    logical                :: isLocked

    isLocked=lock%owned()
    if (.not.isLocked) call lock%setRead()
    if (.not.filterResponses(filterIndex)%wavelengthEffectiveAvailable) then
       ! No effective wavelength was supplied - compute it directly.
       filterResponses(filterIndex)%wavelengthEffectiveAvailable=.true.
       filterResponses(filterIndex)%wavelengthEffective         =+sum(filterResponses(filterIndex)%wavelength*filterResponses(filterIndex)%response) &
            &                                                    /sum(                                        filterResponses(filterIndex)%response)
    end if
    Filter_Wavelength_Effective=filterResponses(filterIndex)%wavelengthEffective
    if (.not.isLocked) call lock%unsetRead()
    return
  end function Filter_Wavelength_Effective

  double precision function Filter_Vega_Offset(indexFilter)
    !!{
    Compute the Vega-AB offset for the given filter.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: FoX_DOM                , only : node                 , destroy
    use            :: Error                  , only : Error_Report
    use            :: Input_Paths            , only : inputPath            , pathTypeDataStatic
    use            :: IO_XML                 , only : XML_Array_Read       , XML_Parse
    use            :: ISO_Varying_String     , only : var_str              , char              , operator(//), assignment(=), &
         &                                            varying_string
    use            :: Numerical_Integration  , only : GSL_Integ_Gauss15    , integrator
    use            :: Numerical_Interpolation, only : interpolator
    use            :: Table_Labels           , only : extrapolationTypeZero
    implicit none
    integer                         , intent(in   )              :: indexFilter
    type            (integrator    )               , allocatable :: integratorVegaBuserV_        , integratorABBuserV_  , &
         &                                                          integratorVegaFilter_        , integratorABFilter_
    type            (interpolator  ), pointer                    :: interpolatorBuserV_          , interpolatorFilter_
    type            (interpolator  ), save                       :: interpolatorVega_
    type            (node          ), pointer                    :: doc
    logical                         , save                       :: vegaLoaded           =.false.
    double precision                , dimension(2)               :: wavelengthRangeBuserV        , wavelengthRangeFilter
    double precision                , dimension(:) , allocatable :: wavelengthVega               , spectrumVega
    integer                                                      :: indexBuserV                  , ioErr
    logical                                                      :: isLocked
    type            (varying_string)                             :: fileName

    isLocked=lock%owned()
    if (.not.isLocked) call lock%setRead()
    if (.not.filterResponses(indexFilter)%vegaOffsetAvailable) then
       if (.not.isLocked) call lock%unsetRead()
       ! Get the Buser_V reference filter.
       indexBuserV         =  Filter_Get_Index        (var_str('Buser_V'   ))
       interpolatorBuserV_ => Filter_Response_Function(         indexBuserV )
       ! Get the interpolator for the given filter.
       interpolatorFilter_ => Filter_Response_Function(         indexFilter )
       ! Load the reference Vega spectrum.
       if (.not.vegaLoaded) then
          !$omp critical (loadVegaSpectrum)
          if (.not.vegaLoaded) then
             fileName=inputPath(pathTypeDataStatic)//'stellarAstrophysics/vega/A0V_Castelli.xml'
             !$omp critical (FoX_DOM_Access)
             doc => XML_Parse(char(fileName),iostat=ioErr)
             if (ioErr /= 0) call Error_Report('failed to read Vega spectrum'//{introspection:location})
             call XML_Array_Read(doc,"datum",wavelengthVega,spectrumVega)
             call destroy(doc)
             !$omp end critical (FoX_DOM_Access)
             interpolatorVega_=interpolator(wavelengthVega,spectrumVega,extrapolationType=extrapolationTypeZero)
             vegaLoaded       =.true.
          end if
          !$omp end critical (loadVegaSpectrum)
       end if
       ! Get wavelength ranges for integration.
       wavelengthRangeBuserV=Filter_Extent(indexBuserV)
       wavelengthRangeFilter=Filter_Extent(indexFilter)
       ! Build integrators.
       if (.not.isLocked) call lock%setRead()
       allocate(integratorVegaBuserV_)
       allocate(integratorVegaFilter_)
       allocate(integratorABBuserV_  )
       allocate(integratorABFilter_  )
       integratorVegaBuserV_=integrator(integrandVegaBuserV,toleranceRelative=1.0d-1,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=10000_c_size_t)
       integratorVegaFilter_=integrator(integrandVegaFilter,toleranceRelative=1.0d-1,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=10000_c_size_t)
       integratorABBuserV_  =integrator(integrandABBuserV  ,toleranceRelative=1.0d-1,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=10000_c_size_t)
       integratorABFilter_  =integrator(integrandABFilter  ,toleranceRelative=1.0d-1,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=10000_c_size_t)
       ! Compute the Vega offset.
       filterResponses(indexFilter)%vegaOffset         =+2.5d0                                                                                     &
            &                                           *log10(                                                                                    &
            &                                                  +integratorVegaFilter_%integrate(wavelengthRangeFilter(1),wavelengthRangeFilter(2)) &
            &                                                  *integratorABBuserV_  %integrate(wavelengthRangeBuserV(1),wavelengthRangeBuserV(2)) &
            &                                                  /integratorVegaBuserV_%integrate(wavelengthRangeBuserV(1),wavelengthRangeBuserV(2)) &
            &                                                  /integratorABFilter_  %integrate(wavelengthRangeFilter(1),wavelengthRangeFilter(2)) &
            &                                                 )
       filterResponses(indexFilter)%vegaOffsetAvailable=.true.
       ! Clean up.
       deallocate(interpolatorBuserV_  )
       deallocate(interpolatorFilter_  )
       deallocate(integratorVegaBuserV_)
       deallocate(integratorVegaFilter_)
       deallocate(integratorABBuserV_  )
       deallocate(integratorABFilter_  )
    end if
    ! Return the offset.
    Filter_Vega_Offset=filterResponses(indexFilter)%vegaOffset
    if (.not.isLocked) call lock%unsetRead()
    return

  contains

    double precision function integrandVegaBuserV(wavelength) result(flux)
      !!{
      Integrand used in calculation of Vega-AB offsets.
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength

      flux   =+interpolatorVega_  %interpolate(wavelength) &
           &  *interpolatorBuserV_%interpolate(wavelength)
      return
    end function integrandVegaBuserV
    
    double precision function integrandABBuserV(wavelength) result(flux)
      !!{
      Integrand used in calculation of Vega-AB offsets.
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength

      flux   =+interpolatorBuserV_%interpolate(wavelength)    &
           &  /                                wavelength **2
      return
    end function integrandABBuserV
    
    double precision function integrandVegaFilter(wavelength) result(flux)
      !!{
      Integrand used in calculation of Vega-AB offsets.
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength

      flux   =+interpolatorVega_  %interpolate(wavelength) &
           &  *interpolatorFilter_%interpolate(wavelength)
      return
    end function integrandVegaFilter
    
    double precision function integrandABFilter(wavelength) result(flux)
      !!{
      Integrand used in calculation of Vega-AB offsets.
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength

      flux   =+interpolatorFilter_%interpolate(wavelength)    &
           &  /                                wavelength **2
      return
    end function integrandABFilter
    
  end function Filter_Vega_Offset
  
end module Instruments_Filters

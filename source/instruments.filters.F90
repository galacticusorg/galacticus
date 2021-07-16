!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  public :: Filter_Get_Index, Filter_Response            , Filter_Extent , Filter_Vega_Offset      , &
       &    Filter_Name     , Filter_Wavelength_Effective, Filters_Output, Filter_Response_Function

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
  integer                                              :: countFilterResponses=0

  ! Read/write lock used to control access to the filter responeses.
  type   (ompReadWriteLock)                            :: lock
  logical                                              :: lockInitialized=.false.
  
contains

  integer function Filter_Get_Index(filterName)
    !!{
    Return the index for the specified filter, loading that filter if necessary.
    !!}
    use :: ISO_Varying_String, only : operator(==)
    implicit none
    type   (varying_string), intent(in   ) :: filterName
    integer                                :: iFilter

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
    call lock%setRead()
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
    call lock%unsetRead()
    return
  end function Filter_Get_Index

  function Filter_Name(filterIndex)
    !!{
    Return the name of the specified filter.
    !!}
    implicit none
    type   (varying_string)                :: Filter_Name
    integer                , intent(in   ) :: filterIndex

    call lock%setRead()
    Filter_Name=filterResponses(filterIndex)%name
    call lock%unsetRead()
    return
  end function Filter_Name

  function Filter_Extent(filterIndex)
    !!{
    Return an array containing the minimum and maximum wavelengths tabulated for this specified filter.
    !!}
    implicit none
    double precision, dimension(2)  :: Filter_Extent
    integer         , intent(in   ) :: filterIndex

    call lock%setRead()
    Filter_Extent(1)=filterResponses(filterIndex)%wavelength(1)
    Filter_Extent(2)=filterResponses(filterIndex)%wavelength(filterResponses(filterIndex)%nPoints)
    call lock%unsetRead()
    return
  end function Filter_Extent

  subroutine Filter_Response_Load(filterName)
    !!{
    Load a filter response curve.
    !!}
    use :: File_Utilities           , only : File_Exists
    use :: FoX_dom                  , only : DOMException           , Node                             , destroy           , getExceptionCode                          , &
         &                                   inException
    use :: Galacticus_Error         , only : Galacticus_Error_Report
    use :: Galacticus_Paths         , only : galacticusPath         , pathTypeDataDynamic              , pathTypeDataStatic
    use :: HII_Region_Emission_Lines, only : emissionLineWavelength
    use :: IO_XML                   , only : XML_Array_Read         , XML_Get_First_Element_By_Tag_Name, XML_Path_Exists   , extractDataContent => extractDataContentTS, &
         &                                   XML_Parse
    use :: ISO_Varying_String       , only : assignment(=)          , char                             , operator(//)      , operator(==)                              , &
         &                                   extract
    use :: Memory_Management        , only : Memory_Usage_Record    , allocateArray
    use :: String_Handling          , only : String_Split_Words     , operator(//)
    implicit none
    type            (varying_string), intent(in   )               :: filterName
    type            (Node          ), pointer                     :: doc
    type            (filterType    ), allocatable  , dimension(:) :: filterResponsesTemporary
    type            (varying_string)               , dimension(4) :: specialFilterWords
    type            (node          ), pointer                     :: vegaElement                   , wavelengthEffectiveElement
    double precision                , parameter                   :: cutOffResolution        =1.0d4
    integer                                                       :: filterIndex                   , ioErr
    type            (varying_string)                              :: filterFileName                , errorMessage
    type            (DOMException  )                              :: exception
    logical                                                       :: parseSuccess
    character       (len=64        )                              :: word
    double precision                                              :: centralWavelength             , resolution                , &
         &                                                           filterWidth

    call lock%setWrite(haveReadLock=.true.)
    ! Allocate space for this filter.
    if (countFilterResponses == 0) then
       allocate(filterResponses(1))
       call Memory_Usage_Record(sizeof(filterResponses))
    else if (countFilterResponses == size(filterResponses)) then
       call Move_Alloc(filterResponses,filterResponsesTemporary)
       allocate(filterResponses(2*size(filterResponsesTemporary)))
       filterResponses(1:size(filterResponsesTemporary))=filterResponsesTemporary
       deallocate(filterResponsesTemporary)       
       call Memory_Usage_Record(sizeof(filterResponses(1)),blockCount=0)
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
       filterResponses(filterIndex)%vegaOffsetAvailable=.false.
       filterResponses(filterIndex)%vegaOffset         =0.0d0
       filterResponses(filterIndex)%nPoints            =4
       call allocateArray(filterResponses(filterIndex)%wavelength,[4])
       call allocateArray(filterResponses(filterIndex)%response  ,[4])
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
    else if (extract(filterName,1,25)=="adaptiveResolutionTopHat_") then
       ! Construct an SED top-hat filter. Extract central wavelength and top hat width.
       call String_Split_Words(specialFilterWords,char(filterName),separator="_")
       word=char(specialFilterWords(2))
       read (word,*) centralWavelength
       word=char(specialFilterWords(3))
       read (word,*) filterWidth
       filterResponses(filterIndex)%vegaOffsetAvailable=.false.
       filterResponses(filterIndex)%vegaOffset         =0.0d0
       filterResponses(filterIndex)%nPoints            =4
       call allocateArray(filterResponses(filterIndex)%wavelength,[4])
       call allocateArray(filterResponses(filterIndex)%response  ,[4])
       filterResponses(filterIndex)%wavelength         =                                                                  &
            & [                                                                                                           &
            &  centralWavelength - filterWidth/2.0d0 - filterWidth/100.0d0                                              , &
            &  centralWavelength - filterWidth/2.0d0                                                                    , &
            &  centralWavelength + filterWidth/2.0d0                                                                    , &
            &  centralWavelength + filterWidth/2.0d0 + filterWidth/100.0d0                                                &
            & ]
       filterResponses(filterIndex)%response           =                                                                  &
            & [                                                                                                           &
            &  0.0d0                                                                                                    , &
            &  1.0d0                                                                                                    , &
            &  1.0d0                                                                                                    , &
            &  0.0d0                                                                                                      &
            & ]
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
       filterResponses(filterIndex)%vegaOffsetAvailable=.false.
       filterResponses(filterIndex)%vegaOffset         =0.0d0
       filterResponses(filterIndex)%nPoints            =4
       call allocateArray(filterResponses(filterIndex)%wavelength,[4])
       call allocateArray(filterResponses(filterIndex)%response  ,[4])
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
    else
       ! Construct a file name for the filter.
       filterFileName=char(galacticusPath(pathTypeDataStatic))//'filters/'//filterName//'.xml'
       if (.not.File_Exists(filterFileName)) then
          filterFileName=char(galacticusPath(pathTypeDataDynamic))//'filters/'//filterName//'.xml'
          if (.not.File_Exists(filterFileName)) call Galacticus_Error_Report('filter file for filter "'//filterName//'" can not be found'//{introspection:location})
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
       if (.not.parseSuccess)                                                                 &
            & call Galacticus_Error_Report(                                                   &
            &                              'unable to read or parse filter response file: '// &
            &                              char(filterFileName)                            // &
            &                                   errorMessage                               // &
            &                              {introspection:location}                           &
            &                             )
       ! Extract wavelengths and filter response.
       call XML_Array_Read(doc,"datum",filterResponses(filterIndex)%wavelength,filterResponses(filterIndex)%response)
       filterResponses(filterIndex)%nPoints=size(filterResponses(filterIndex)%wavelength)
       ! Extract the Vega offset.
       filterResponses(filterIndex)%vegaOffsetAvailable=XML_Path_Exists(doc,"vegaOffset")
       if (filterResponses(filterIndex)%vegaOffsetAvailable) then
          vegaElement => XML_Get_First_Element_By_Tag_Name(doc,"vegaOffset")
          call extractDataContent(vegaElement,filterResponses(filterIndex)%vegaOffset)
       else
          filterResponses(filterIndex)%vegaOffset=0.0d0
       end if
       ! Extract the effective wavelength.
       filterResponses(filterIndex)%wavelengthEffectiveAvailable=XML_Path_Exists(doc,"effectiveWavelength")
       if (filterResponses(filterIndex)%wavelengthEffectiveAvailable) then
          wavelengthEffectiveElement => XML_Get_First_Element_By_Tag_Name(doc,"effectiveWavelength")
          call extractDataContent(wavelengthEffectiveElement,filterResponses(filterIndex)%wavelengthEffective)
       else
          filterResponses(filterIndex)%wavelengthEffective=0.0d0
       end if
       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
    end if
    ! No effective wavelength was supplied - compute it directly.
    filterResponses(filterIndex)%wavelengthEffectiveAvailable=.true.
    filterResponses(filterIndex)%wavelengthEffective         =+sum(filterResponses(filterIndex)%wavelength*filterResponses(filterIndex)%response) &
         &                                                    /sum(                                        filterResponses(filterIndex)%response)
    call lock%unsetWrite(haveReadLock=.true.)
    return
  end subroutine Filter_Response_Load

  !![
  <hdfPreCloseTask>
   <unitName>Filters_Output</unitName>
  </hdfPreCloseTask>
  !!]
  subroutine Filters_Output()
    !!{
    Output accumulated filter data to file.
    !!}
    use :: Galacticus_HDF5          , only : galacticusOutputFile
    use :: HDF5_Access              , only : hdf5Access
    use :: IO_HDF5                  , only : hdf5Object
    use :: Numerical_Constants_Units, only : angstromsPerMeter
    implicit none
    type(hdf5Object) :: filtersGroup, dataset

    if (.not.allocated(filterResponses)) return
    !$ call hdf5Access%set()
    filtersGroup=galacticusOutputFile%openGroup('Filters','Properties of filters used.')
    call filtersGroup%writeDataset(filterResponses%name               ,'name'               ,'Filter name.'                                               )
    call filtersGroup%writeDataset(filterResponses%wavelengthEffective,'wavelengthEffective','Effective wavelength of filter [Å].',datasetReturned=dataset)
    call dataset%writeAttribute("Angstroms [Å]"        ,"units"    )
    call dataset%writeAttribute(1.0d0/angstromsPerMeter,"unitsInSI")
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

    allocate(interpolator_)
    call lock%setRead()
    interpolator_=interpolator(filterResponses(filterIndex)%wavelength,filterResponses(filterIndex)%response)
    call lock%unsetRead()
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

  double precision function Filter_Vega_Offset(filterIndex)
    !!{
    Return the Vega to AB magnitude offset for the specified filter.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    integer, intent(in   ) :: filterIndex

    call lock%setRead()
    if (.not.filterResponses(filterIndex)%vegaOffsetAvailable) call Galacticus_Error_Report('Vega offset is not available'//{introspection:location})
    Filter_Vega_Offset=filterResponses(filterIndex)%vegaOffset
    call lock%unsetRead()
   return
  end function Filter_Vega_Offset

  double precision function Filter_Wavelength_Effective(filterIndex)
    !!{
    Return the effective wavelength for the specified filter.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    integer, intent(in   ) :: filterIndex

    call lock%setRead()
    if (.not.filterResponses(filterIndex)%wavelengthEffectiveAvailable) call Galacticus_Error_Report('effective wavelength is not available'//{introspection:location})
    Filter_Wavelength_Effective=filterResponses(filterIndex)%wavelengthEffective
    call lock%unsetRead()
    return
  end function Filter_Wavelength_Effective

end module Instruments_Filters

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of filter response curves.

module Instruments_Filters
  !% Implements calculations of filter response curves.
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Filter_Get_Index, Filter_Response, Filter_Extent, Filter_Vega_Offset, Filter_Name

  type filterType
     !% A structure which holds filter response curves.
     integer                                                        :: nPoints
     double precision                   , allocatable, dimension(:) :: response                       , wavelength
     type            (varying_string   )                            :: name
     logical                                                        :: vegaOffsetAvailable
     double precision                                               :: vegaOffset
     ! Interpolation structures.
     logical                                                        :: reset                   =.true.
     type            (fgsl_interp_accel)                            :: interpolationAccelerator
     type            (fgsl_interp      )                            :: interpolationObject
  end type filterType

  ! Array to hold filter data.
  type(filterType), allocatable, dimension(:) :: filterResponses

contains

  integer function Filter_Get_Index(filterName)
    !% Return the index for the specified filter, loading that filter if necessary.
    implicit none
    type   (varying_string), intent(in   ) :: filterName
    integer                                :: iFilter

    ! See if we already have this filter loaded. If not, load it.
    !$omp critical (Filter_Get_Index_Lock)
    if (.not.allocated(filterResponses)) then
       call Filter_Response_Load(filterName)
       Filter_Get_Index=1
    else
       Filter_Get_Index=0
       do iFilter=1,size(filterResponses)
          if (filterResponses(iFilter)%name == filterName) then
             Filter_Get_Index=iFilter
             exit
          end if
       end do
       if (Filter_Get_Index == 0) then
          call Filter_Response_Load(filterName)
          Filter_Get_Index=size(filterResponses)
       end if
    end if
    !$omp end critical (Filter_Get_Index_Lock)
    return
  end function Filter_Get_Index

  function Filter_Name(filterIndex)
    !% Return the name of the specified filter.
    implicit none
    type   (varying_string)                :: Filter_Name
    integer                , intent(in   ) :: filterIndex

    Filter_Name=filterResponses(filterIndex)%name
    return
  end function Filter_Name

  function Filter_Extent(filterIndex)
    !% Return an array containing the minimum and maximum wavelengths tabulated for this specified filter.
    implicit none
    double precision, dimension(2)  :: Filter_Extent
    integer         , intent(in   ) :: filterIndex

    !$omp critical (Filter_Get_Index_Lock)
    Filter_Extent(1)=filterResponses(filterIndex)%wavelength(1)
    Filter_Extent(2)=filterResponses(filterIndex)%wavelength(filterResponses(filterIndex)%nPoints)
    !$omp end critical (Filter_Get_Index_Lock)
    return
  end function Filter_Extent

  subroutine Filter_Response_Load(filterName)
    !% Load a filter response curve.
    use Memory_Management
    use Galacticus_Error
    use FoX_dom
    use Galacticus_Input_Paths
    use IO_XML
    use String_Handling
    implicit none
    type            (varying_string), intent(in   )               :: filterName
    type            (Node          ), pointer                     :: doc
    type            (filterType    ), allocatable  , dimension(:) :: filterResponsesTemporary
    type            (node          ), pointer                     :: vegaElement
    integer                                                       :: filterIndex             , ioErr
    type            (varying_string)                              :: filterFileName          , errorMessage
    type            (DOMException  )                              :: exception
    logical                                                       :: parseSuccess

    ! Allocate space for this filter.
    if (allocated(filterResponses)) then
       call Move_Alloc(filterResponses,filterResponsesTemporary)
       allocate(filterResponses(size(filterResponsesTemporary)+1))
       filterResponses(1:size(filterResponsesTemporary))=filterResponsesTemporary
       deallocate(filterResponsesTemporary)
       call Memory_Usage_Record(sizeof(filterResponses(1)),blockCount=0)
    else
       allocate(filterResponses(1))
       call Memory_Usage_Record(sizeof(filterResponses))
    end if

    ! Index in array to load into.
    filterIndex=size(filterResponses)

    ! Store the name of the filter.
    filterResponses(filterIndex)%name=filterName

    ! Construct a file name for the filter.
    filterFileName=char(Galacticus_Input_Path())//'data/filters/'//filterName//'.xml'

    ! Parse the XML file.
    !$omp critical (FoX_DOM_Access)
    doc => parseFile(char(filterFileName),iostat=ioErr,ex=exception)
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
         &                              'Filter_Response_Load'                          ,  &
         &                              'unable to read or parse filter response file: '// &
         &                              char(filterFileName)                            // &
         &                                   errorMessage                                  &
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
    ! Destroy the document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    return
  end subroutine Filter_Response_Load

  double precision function Filter_Response(filterIndex,wavelength)
    !% Return the filter response function at the given {\tt wavelength} (specified in Angstroms).  Note that we follow the
    !% convention of \cite{hogg_k_2002} and assume that the filter response gives the fraction of incident photons received by the
    !% detector at a given wavelength, multiplied by the relative photon response (which will be 1 for a photon-counting detector
    !% such as a CCD, or proportional to the photon energy for a bolometer/calorimeter type detector.
    use Numerical_Interpolation
    implicit none
    integer         , intent(in   ) :: filterIndex
    double precision, intent(in   ) :: wavelength

    ! Interpolate in the tabulated response curve.
    !$omp critical (Filter_Get_Index_Lock)
    Filter_Response=Interpolate(filterResponses(filterIndex)%nPoints,filterResponses(filterIndex)%wavelength&
         &,filterResponses(filterIndex)%response,filterResponses(filterIndex)%interpolationObject&
         &,filterResponses(filterIndex)%interpolationAccelerator,wavelength ,reset=filterResponses(filterIndex)%reset)
    !$omp end critical (Filter_Get_Index_Lock)
    return
  end function Filter_Response

  double precision function Filter_Vega_Offset(filterIndex)
    !% Return the Vega to AB magnitude offset for the specified filter.
    use Galacticus_Error
    implicit none
    integer, intent(in   ) :: filterIndex

    if (.not.filterResponses(filterIndex)%vegaOffsetAvailable) call Galacticus_Error_Report('Filter_Vega_Offset','Vega offset is not available')
    Filter_Vega_Offset=filterResponses(filterIndex)%vegaOffset
    return
  end function Filter_Vega_Offset

end module Instruments_Filters

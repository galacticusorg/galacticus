!% Contains a module which implements calculations of filter response curves.

module Instruments_Filters
  !% Implements calculations of filter response curves.
  use ISO_Varying_String
  use Instruments_Filters_Type
  private
  public :: Filter_Get_Index, Filter_Response, Filter_Extent

  ! Array to hold filter data.
  type(filterType), allocatable, dimension(:) :: filterResponses

contains

  integer function Filter_Get_Index(filterName)
    !% Return the index for the specified filter, loading that filter if necessary.
    implicit none
    type(varying_string), intent(in) :: filterName
    integer                          :: iFilter

    ! See if we already have this filter loaded. If not, load it.
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
    return
  end function Filter_Get_Index

  function Filter_Extent(filterIndex)
    !% Return an array containing the minimum and maximum wavelengths tabulated for this specified filter.
    implicit none
    double precision, dimension(2) :: Filter_Extent
    integer,          intent(in)   :: filterIndex
    
    Filter_Extent(1)=filterResponses(filterIndex)%wavelength(1)
    Filter_Extent(2)=filterResponses(filterIndex)%wavelength(filterResponses(filterIndex)%nPoints)
    return
  end function Filter_Extent

  subroutine Filter_Response_Load(filterName)
    !% Load a filter response curve.
    use Memory_Management
    use Galacticus_Error
    use FoX_dom
    implicit none
    type(varying_string), intent(in)                :: filterName
    type(Node),           pointer                   :: doc,datum
    type(NodeList),       pointer                   :: datumList
    type(filterType),     allocatable, dimension(:) :: filterResponsesTemporary
    integer                                         :: iDatum,ioErr,filterIndex
    double precision                                :: datumValues(2)
    type(varying_string)                            :: filterFileName

    ! Allocate space for this filter.
    if (allocated(filterResponses)) then
       call Move_Alloc(filterResponses,filterResponsesTemporary)
       call Alloc_Array(filterResponses,size(filterResponsesTemporary)+1,'filterResponses')
       filterResponses(1:size(filterResponsesTemporary))=filterResponsesTemporary
       call Dealloc_Array(filterResponsesTemporary)
    else
       call Alloc_Array(filterResponses,1,'filterResponses')
    end if

    ! Index in array to load into.
    filterIndex=size(filterResponses)

    ! Store the name of the filter.
    filterResponses(filterIndex)%name=filterName

    ! Construct a file name for the filter.
    filterFileName='data/filters/'//filterName//'.xml'

    ! Parse the XML file.
    !$omp critical (FoX_DOM_Access)
    doc => parseFile(char(filterFileName),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Filter_Response_Load','unable to parse filter response file: '//char(filterFileName))

    ! Get a list of all <datum> elements.
    datumList => getElementsByTagname(doc,"datum")
    filterResponses(filterIndex)%nPoints=getLength(datumList)

    ! Allocate space for the response curve.
    call Alloc_Array(filterResponses(filterIndex)%wavelength,filterResponses(filterIndex)%nPoints,'filterResponses()%wavelength')
    call Alloc_Array(filterResponses(filterIndex)%response  ,filterResponses(filterIndex)%nPoints,'filterResponses()%response'  )

    ! Extract the data from the file.
    do iDatum=0,filterResponses(filterIndex)%nPoints-1
       datum => item(datumList,iDatum)
       call extractDataContent(datum,datumValues)
       filterResponses(filterIndex)%wavelength(iDatum+1)=datumValues(1)
       filterResponses(filterIndex)%response  (iDatum+1)=datumValues(2)
    end do

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
    integer,          intent(in) :: filterIndex
    double precision, intent(in) :: wavelength

    ! Interpolate in the tabulated response curve.
    Filter_Response=Interpolate(filterResponses(filterIndex)%nPoints,filterResponses(filterIndex)%wavelength&
         &,filterResponses(filterIndex)%response,filterResponses(filterIndex)%interpolationObject&
         &,filterResponses(filterIndex)%interpolationAccelerator,wavelength ,reset=filterResponses(filterIndex)%reset)
    return
  end function Filter_Response

end module Instruments_Filters

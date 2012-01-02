!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements calculations of filter response curves.

module Instruments_Filters
  !% Implements calculations of filter response curves.
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Filter_Get_Index, Filter_Response, Filter_Extent

  type filterType
     !% A structure which holds filter response curves.
     integer                                         :: nPoints
     double precision,     allocatable, dimension(:) :: wavelength,response
     type(varying_string)                            :: name
     ! Interpolation structures.
     logical                                         :: reset=.true.
     type(fgsl_interp_accel)                         :: interpolationAccelerator
     type(fgsl_interp)                               :: interpolationObject
  end type filterType

  ! Array to hold filter data.
  type(filterType), allocatable, dimension(:) :: filterResponses

contains

  integer function Filter_Get_Index(filterName)
    !% Return the index for the specified filter, loading that filter if necessary.
    implicit none
    type(varying_string), intent(in) :: filterName
    integer                          :: iFilter

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
    use Galacticus_Input_Paths
    use ISO_Varying_String
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
    doc => parseFile(char(filterFileName),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Filter_Response_Load','unable to parse filter response file: '//char(filterFileName))

    ! Get a list of all <datum> elements.
    datumList => getElementsByTagname(doc,"datum")
    filterResponses(filterIndex)%nPoints=getLength(datumList)

    ! Allocate space for the response curve.
    call Alloc_Array(filterResponses(filterIndex)%wavelength,[filterResponses(filterIndex)%nPoints])
    call Alloc_Array(filterResponses(filterIndex)%response  ,[filterResponses(filterIndex)%nPoints])

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

!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which reads a tabulated transfer function from a file.

module Transfer_Function_File
  !% Implements reading of a tabulated transfer function from a file.
  use ISO_Varying_String
  private
  public :: Transfer_Function_File_Initialize, Transfer_Function_Named_File_Read
  
  ! Flag to indicate if this module has been initialized.
  logical              :: transferFunctionInitialized=.false.

  ! File name for the transfer function data.
  type(varying_string) :: transferFunctionFile

  ! Extrapolation methods.
  integer              :: extrapolateWavenumberLow,extrapolateWavenumberHigh

  ! Number of points per decade to add per decade when extrapolating and the buffer in log wavenumber to use.
  integer,         parameter :: extrapolatePointsPerDecade=10
  double precision           :: extrapolateLogWavenumberBuffer=1.0d0

contains
  
  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_File_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_File_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``transfer function from file'' module.
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: transferFunctionMethod
    procedure(),          pointer, intent(inout) :: Transfer_Function_Tabulate
    
    if (transferFunctionMethod.eq.'file') then
       Transfer_Function_Tabulate => Transfer_Function_File_Read
       !@ <inputParameter>
       !@   <name>transferFunctionFile</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of a file containing a tabulation of the transfer function for the ``file'' transfer function method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('transferFunctionFile',transferFunctionFile)
    end if
    return
  end subroutine Transfer_Function_File_Initialize

  subroutine Transfer_Function_Named_File_Read(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT,fileName)
    !% Read the transfer function from a file specified by {\tt fileName}.
    implicit none
    double precision,                                intent(in)    :: logWavenumber
    double precision,     allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                         intent(out)   :: transferFunctionNumberPoints
    type(varying_string),                            intent(in)    :: fileName

    ! Set the filename to that specified.
    transferFunctionFile=fileName
    ! Flag that module is uninitialized again.
    transferFunctionInitialized=.false.
    ! Call routine to read in the data.
    call Transfer_Function_File_Read(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
         &,transferFunctionLogT)
   return
  end subroutine Transfer_Function_Named_File_Read

  subroutine Transfer_Function_File_Read(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Reads a transfer function from an XML file.
    use FoX_dom
    use Numerical_Comparison
    use Memory_Management
    use Galacticus_Error
    use Galacticus_Display
    use Cosmological_Parameters
    use Numerical_Constants_Math
    use Numerical_Ranges
    use IO_XML
    implicit none
    double precision,                            intent(in)    :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                     intent(out)   :: transferFunctionNumberPoints
    type(Node),       pointer                                  :: doc,datum,thisParameter,nameElement,valueElement&
         &,extrapolationElement,extrapolation
    type(NodeList),   pointer                                  :: datumList,parameterList,wavenumberExtrapolationList
    double precision, allocatable, dimension(:)                :: wavenumberTemporary,transferFunctionTemporary
    integer                                                    :: iDatum,ioErr,iParameter,addCount,iExtrapolation&
         &,extrapolationMethod
    double precision                                           :: datumValues(2),parameterValue
    character(len=32)                                          :: limitType

    ! Read the file if this module has not been initialized.
    if (.not.transferFunctionInitialized) then
       ! Open and parse the data file.
       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(transferFunctionFile),iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Transfer_Function_File_Read','Unable to find transfer function file')
       ! Check that parameters match if any are present.
       parameterList => getElementsByTagname(doc,"parameter")
       do iParameter=0,getLength(parameterList)-1
          thisParameter => item(parameterList,iParameter)
          nameElement => item(getElementsByTagname(thisParameter,"name"),0)
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          call extractDataContent(valueElement,parameterValue)
          select case (getTextContent(nameElement))
          case ("Omega_b")
             if (Values_Differ(parameterValue,Omega_b(),absTol=1.0d-3)) call Galacticus_Display_Message('Omega_b from transfer &
                  & function file does not match internal value')
          case ("Omega_0")
             if (Values_Differ(parameterValue,Omega_0(),absTol=1.0d-3)) call Galacticus_Display_Message('Omega_0 from transfer &
                  & function file does not match internal value')
          case ("Omega_DE")
             if (Values_Differ(parameterValue,Omega_DE(),absTol=1.0d-3)) call Galacticus_Display_Message('Omega_DE from transfer &
                  & function file does not match internal value')
          case ("H_0")
             if (Values_Differ(parameterValue,H_0(),relTol=1.0d-3)) call Galacticus_Display_Message('H_0 from transfer &
                  & function file does not match internal value')
          case ("T_CMB")
             if (Values_Differ(parameterValue,T_CMB(),relTol=1.0d-3)) call Galacticus_Display_Message('T_CMB from transfer &
                  & function file does not match internal value')
          end select
       end do
       ! Get extrapolation methods.
       extrapolationElement        => item(getElementsByTagname(doc                 ,"extrapolation"),0)
       wavenumberExtrapolationList =>      getElementsByTagname(extrapolationElement,"wavenumber"   )
       do iExtrapolation=0,getLength(wavenumberExtrapolationList)-1
          extrapolation => item(wavenumberExtrapolationList,iExtrapolation)
          call XML_Extrapolation_Element_Decode(extrapolation,limitType,extrapolationMethod,allowedMethods=[extrapolateFixed,extrapolatePowerLaw])
          select case (trim(limitType))
          case ('low')
             extrapolateWavenumberLow=extrapolationMethod
          case ('high')
             extrapolateWavenumberHigh=extrapolationMethod
          case default
             call Galacticus_Error_Report('Transfer_Function_File_Read','unrecognized extrapolation limit')
          end select
       end do
       ! Get list of datum elements.
       datumList => getElementsByTagname(doc,"datum")
       ! Allocate transfer function arrays.
       transferFunctionNumberPoints=getLength(datumList)
       ! Deallocate arrays if currently allocated.
       if (allocated(transferFunctionLogWavenumber)) call Dealloc_Array(transferFunctionLogWavenumber)
       if (allocated(transferFunctionLogT))          call Dealloc_Array(transferFunctionLogT         )
       ! Allocate the arrays to current required size.
       call Alloc_Array(transferFunctionLogWavenumber,[transferFunctionNumberPoints])
       call Alloc_Array(transferFunctionLogT         ,[transferFunctionNumberPoints])
       ! Extract data from elements.
       do iDatum=0,getLength(datumList)-1
          datum => item(datumList,iDatum)
          call extractDataContent(datum,datumValues)
          transferFunctionLogWavenumber(iDatum+1)=dlog(datumValues(1))
          transferFunctionLogT         (iDatum+1)=dlog(datumValues(2))
       end do
       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
       ! Flag that transfer function is now initialized.
       transferFunctionInitialized=.true.
    end if
    ! Check that the input wavenumber is within range, extend the range if possible.
    if (logWavenumber < transferFunctionLogWavenumber(1)) then
       ! Determine how many extra points to add.
       addCount=int((transferFunctionLogWavenumber(1)-logWavenumber+extrapolateLogWavenumberBuffer)*dble(extrapolatePointsPerDecade)&
            &/ln10)+1
       ! Allocate temporary arrays.
       call Alloc_Array(wavenumberTemporary      ,[size(transferFunctionLogWavenumber)+addCount])
       call Alloc_Array(transferFunctionTemporary,[size(transferFunctionLogWavenumber)+addCount])
       ! Create additional wavenumber range.
       wavenumberTemporary(1:addCount+1)=Make_Range(logWavenumber-extrapolateLogWavenumberBuffer,transferFunctionLogWavenumber(1)&
            &,addCount+1,rangeType=rangeTypeLinear)
       ! Extrapolate as directed.
       select case (extrapolateWavenumberLow)
       case (extrapolateFixed   )
          transferFunctionTemporary(1:addCount)=transferFunctionLogT(1)
       case (extrapolatePowerLaw)
          transferFunctionTemporary(1:addCount)=transferFunctionLogT(1)+(wavenumberTemporary(1:addCount)&
               &-transferFunctionLogWavenumber(1))*(transferFunctionLogT(2)-transferFunctionLogT(1))&
               &/(transferFunctionLogWavenumber(2)-transferFunctionLogWavenumber(1))
       end select
       ! Copy in original wavenumber and transfer function data.
       wavenumberTemporary      (addCount+1:size(wavenumberTemporary))=transferFunctionLogWavenumber
       transferFunctionTemporary(addCount+1:size(wavenumberTemporary))=transferFunctionLogT
       ! Move the new tabulation into the output arrays.
       call Dealloc_Array(                          transferFunctionLogWavenumber)
       call Dealloc_Array(                          transferFunctionLogT         )
       call Move_Alloc   (wavenumberTemporary      ,transferFunctionLogWavenumber)
       call Move_Alloc   (transferFunctionTemporary,transferFunctionLogT         )
       ! Reset the number of tabulated points.
       transferFunctionNumberPoints=size(transferFunctionLogWavenumber)
    end if
    if (logWavenumber > transferFunctionLogWavenumber(transferFunctionNumberPoints)) then
       ! Determine how many extra points to add.
       addCount=int((logWavenumber-transferFunctionLogWavenumber(transferFunctionNumberPoints)+extrapolateLogWavenumberBuffer)&
            &*dble(extrapolatePointsPerDecade)/ln10)+1
       ! Allocate temporary arrays.
       call Alloc_Array(wavenumberTemporary      ,[size(transferFunctionLogWavenumber)+addCount])
       call Alloc_Array(transferFunctionTemporary,[size(transferFunctionLogWavenumber)+addCount])
       ! Create additional wavenumber range.
       wavenumberTemporary(transferFunctionNumberPoints:size(wavenumberTemporary))&
            &=Make_Range(transferFunctionLogWavenumber(transferFunctionNumberPoints),logWavenumber+extrapolateLogWavenumberBuffer&
            & ,addCount+1,rangeType=rangeTypeLinear)
       ! Extrapolate as directed.
       select case (extrapolateWavenumberLow)
       case (extrapolateFixed   )
          transferFunctionTemporary(transferFunctionNumberPoints+1:size(wavenumberTemporary))=transferFunctionLogT(transferFunctionNumberPoints)
       case (extrapolatePowerLaw)
          transferFunctionTemporary(transferFunctionNumberPoints+1:size(wavenumberTemporary))&
               &=transferFunctionLogT(transferFunctionNumberPoints)+(wavenumberTemporary(transferFunctionNumberPoints&
               &+1:size(wavenumberTemporary))-transferFunctionLogWavenumber(transferFunctionNumberPoints))&
               &*(transferFunctionLogT(transferFunctionNumberPoints)-transferFunctionLogT(transferFunctionNumberPoints-1)) &
               &/(transferFunctionLogWavenumber(transferFunctionNumberPoints)&
               &-transferFunctionLogWavenumber(transferFunctionNumberPoints-1))
       end select
       ! Copy in original wavenumber and transfer function data.
       wavenumberTemporary      (1:transferFunctionNumberPoints)=transferFunctionLogWavenumber
       transferFunctionTemporary(1:transferFunctionNumberPoints)=transferFunctionLogT
       ! Move the new tabulation into the output arrays.
       call Dealloc_Array(                          transferFunctionLogWavenumber)
       call Dealloc_Array(                          transferFunctionLogT         )
       call Move_Alloc   (wavenumberTemporary      ,transferFunctionLogWavenumber)
       call Move_Alloc   (transferFunctionTemporary,transferFunctionLogT         )
       ! Reset the number of tabulated points.
       transferFunctionNumberPoints=size(transferFunctionLogWavenumber)
    end if
    return
  end subroutine Transfer_Function_File_Read
  
end module Transfer_Function_File

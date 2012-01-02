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


!% Contains a module that implements calculations of the intergalactic medium thermal and ionization state read from a file.

module Intergalactic_Medium_State_File
  !% Implements calculations of the intergalactic medium thermal and ionization state read from a file.
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Intergalactic_Medium_State_File_Initialize, Intergalactic_Medium_Temperature_File,&
       & Intergalactic_Medium_Electron_Fraction_File, Intergalactic_Medium_File_Set_File

  ! Name of the file from which to read intergalactic medium state data.
  type(varying_string)                               :: intergalaticMediumStateFileName

  ! Flag indicating whether or not data has been read.
  logical                                            :: intergalaticMediumStateDataRead=.false.

  ! Data tables.
  integer                                            :: redshiftCount
  double precision,        allocatable, dimension(:) :: timeTable,electronFractionTable,temperatureTable
  
  ! Interpolation objects.
  type(fgsl_interp_accel)                            :: interpolationAcceleratorElectronFraction ,interpolationAcceleratorTemperature
  type(fgsl_interp)                                  :: interpolationObjectElectronFraction      ,interpolationObjectTemperature
  logical                                            :: interpolationResetElectronFraction=.true.,interpolationResetTemperature=.true.
  !$omp threadprivate(interpolationAcceleratorElectronFraction,interpolationObjectElectronFraction,interpolationResetElectronFraction)
  !$omp threadprivate(interpolationAcceleratorTemperature     ,interpolationObjectTemperature     ,interpolationResetTemperature     )
  
contains

  !# <intergalaticMediumStateMethod>
  !#  <unitName>Intergalactic_Medium_State_File_Initialize</unitName>
  !# </intergalaticMediumStateMethod>
  subroutine Intergalactic_Medium_State_File_Initialize(intergalaticMediumStateMethod,Intergalactic_Medium_Electron_Fraction_Get&
       &,Intergalactic_Medium_Temperature_Get)
    !% Initializes the ``file'' intergalactic medium state module.
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: intergalaticMediumStateMethod
    procedure(),          pointer, intent(inout) :: Intergalactic_Medium_Electron_Fraction_Get,Intergalactic_Medium_Temperature_Get

    ! Test if our method has been selected.    
    if (intergalaticMediumStateMethod == 'file') then
       ! Set procedure pointers.
       Intergalactic_Medium_Electron_Fraction_Get => Intergalactic_Medium_Electron_Fraction_File
       Intergalactic_Medium_Temperature_Get       => Intergalactic_Medium_Temperature_File
       ! Get the name of the file to read.
       !@ <inputParameter>
       !@   <name>intergalaticMediumStateFileName</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file from which to read intergalactic medium state data.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('intergalaticMediumStateFileName',intergalaticMediumStateFileName)
    end if
    return
  end subroutine Intergalactic_Medium_State_File_Initialize

  subroutine Intergalactic_Medium_File_Set_File(fileName)
    !% Allow an external module to set the filename to be used for intergalatic medium state data.
    type(varying_string), intent(in) :: fileName

    intergalaticMediumStateFileName=fileName
    return
  end subroutine Intergalactic_Medium_File_Set_File

  subroutine Intergalactic_Medium_State_File_Read_Data
    !% Read in data describing the state of the intergalactic medium.
    use Galacticus_Error
    use FoX_dom
    use Memory_Management
    use Cosmology_Functions
    implicit none
    type(Node),      pointer :: doc,thisItem
    type(NodeList),  pointer :: itemList,redshiftList,electronFractionList,temperatureList
    integer                  :: ioErr,iRedshift
    double precision         :: redshift

    ! Check if data has yet to be read.
    if (.not.intergalaticMediumStateDataRead) then

       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(intergalaticMediumStateFileName),iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Intergalactic_Medium_State_File_Read_Data','Unable to parse intergalactic medium state file')
       ! Get the redshift element.
       itemList             => getElementsByTagname(doc     ,"redshift"         )
       thisItem             => item(itemList,0)
       redshiftList         => getElementsByTagname(thisItem,"datum"            )
       itemList             => getElementsByTagname(doc     ,"electronFraction" )
       thisItem             => item(itemList,0)
       electronFractionList => getElementsByTagname(thisItem,"datum"            )
       itemList             => getElementsByTagname(doc     ,"matterTemperature")
       thisItem             => item(itemList,0)
       temperatureList      => getElementsByTagname(thisItem,"datum"            )
       ! Count the number of tabulated redshifts.
       redshiftCount=getLength(redshiftList)
       ! Allocate arrays for table storage.
       call Alloc_Array(timeTable            ,[redshiftCount])
       call Alloc_Array(electronFractionTable,[redshiftCount])
       call Alloc_Array(temperatureTable     ,[redshiftCount])
       ! Extract data.
       do iRedshift=1,redshiftCount
          thisItem => item(redshiftList         ,iRedshift-1)
          call extractDataContent(thisItem,redshift                        )
          timeTable(iRedshift)=Cosmology_Age(Expansion_Factor_from_Redshift(redshift))
          thisItem => item(electronFractionList,iRedshift-1)
          call extractDataContent(thisItem,electronFractionTable(iRedshift))
          thisItem => item(temperatureList     ,iRedshift-1)
          call extractDataContent(thisItem,temperatureTable     (iRedshift))
       end do
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)

       ! Flag that data has now been read.
       intergalaticMediumStateDataRead=.true.
    end if
    return
  end subroutine Intergalactic_Medium_State_File_Read_Data

  double precision function Intergalactic_Medium_Temperature_File(time)
    !% Return the temperature of the intergalactic medium at the specified {\tt time} by interpolating in tabulated data,
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: time

    ! Ensure that data has been read.
    call Intergalactic_Medium_State_File_Read_Data

    ! Interpolate in the tables to get the electron fraction.
    Intergalactic_Medium_Temperature_File=Interpolate(redshiftCount,timeTable,temperatureTable&
         &,interpolationObjectTemperature,interpolationAcceleratorTemperature,time,reset&
         &=interpolationResetTemperature)

    return
  end function Intergalactic_Medium_Temperature_File

  double precision function Intergalactic_Medium_Electron_Fraction_File(time)
    !% Return the electron fraction in the intergalactic medium at the specified {\tt time} by interpolating in tabulated data,
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: time
    
    ! Ensure that data has been read.
    call Intergalactic_Medium_State_File_Read_Data

    ! Interpolate in the tables to get the electron fraction.
    Intergalactic_Medium_Electron_Fraction_File=Interpolate(redshiftCount,timeTable,electronFractionTable&
         &,interpolationObjectElectronFraction ,interpolationAcceleratorElectronFraction,time ,reset&
         &=interpolationResetElectronFraction)

    return
  end function Intergalactic_Medium_Electron_Fraction_File

end module Intergalactic_Medium_State_File

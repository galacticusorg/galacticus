!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!+    Contributions to this file made by:  Luiz Felippe S. Rodrigues.

!% Contains a module that implements calculations of the intergalactic medium thermal and ionization state read from a file.

module Intergalactic_Medium_State_File
  !% Implements calculations of the intergalactic medium thermal and ionization state read from a file.
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Intergalactic_Medium_State_File_Initialize, Intergalactic_Medium_Temperature_File,&
       & Intergalactic_Medium_Electron_Fraction_File, Intergalactic_Medium_File_Set_File, &
       & Intergalactic_Medium_State_File_Current_File_Format_Version

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

  ! Current file format version for intergalactic medium state files.
  integer,                 parameter                 :: fileFormatVersionCurrent=1

contains

  !# <intergalaticMediumStateMethod>
  !#  <unitName>Intergalactic_Medium_State_File_Initialize</unitName>
  !# </intergalaticMediumStateMethod>
  subroutine Intergalactic_Medium_State_File_Initialize(intergalaticMediumStateMethod,Intergalactic_Medium_Electron_Fraction_Get&
       &,Intergalactic_Medium_Temperature_Get)
    !% Initializes the ``file'' intergalactic medium state module.
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in   ) :: intergalaticMediumStateMethod
    procedure(double precision), pointer, intent(inout) :: Intergalactic_Medium_Electron_Fraction_Get,Intergalactic_Medium_Temperature_Get
  
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

  integer function Intergalactic_Medium_State_File_Current_File_Format_Version()
    !% Return the current file format version of intergalactic medium state files.
    implicit none

    Intergalactic_Medium_State_File_Current_File_Format_Version=fileFormatVersionCurrent
    return
  end function Intergalactic_Medium_State_File_Current_File_Format_Version

  subroutine Intergalactic_Medium_State_File_Read_Data
    !% Read in data describing the state of the intergalactic medium.
    use Galacticus_Error
    use FoX_dom
    use Memory_Management
    use Cosmology_Functions
    implicit none
    type(Node),      pointer :: doc,thisItem
    type(NodeList),  pointer :: itemList,redshiftList,electronFractionList,temperatureList
    integer                  :: ioErr,iRedshift,fileFormatVersion
    double precision         :: redshift

    ! Check if data has yet to be read.
    if (.not.intergalaticMediumStateDataRead) then

       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(intergalaticMediumStateFileName),iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Intergalactic_Medium_State_File_Read_Data','Unable to parse intergalactic medium state file')
       ! Check the file format version of the file.
       itemList             => getElementsByTagname(doc,"fileFormat")
       thisItem             => item(itemList,0)
       call extractDataContent(thisItem,fileFormatVersion)
       if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('Intergalactic_Medium_State_File_Read_Data','file format version is out of date')

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

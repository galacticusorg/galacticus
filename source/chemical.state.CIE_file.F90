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


!% Contains a module which reads and interpolates a collisional ionization equilibrium chemical state from a file.

module Chemical_States_CIE_File
  !% Reads and interpolates a collisional ionization equilibrium chemical state from a file.
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Chemical_State_CIE_File_Initialize, Chemical_State_CIE_File_Read, Electron_Density_CIE_File_Interpolate,&
       & Electron_Density_CIE_File_logTemperature_Interpolate, Chemical_Densities_CIE_File_Interpolate,&
       & Chemical_State_CIE_File_Format_Version

  ! Flag to indicate if this module has been initialized.
  logical              :: chemicalStateInitialized         =.false.
  logical              :: chemicalStateChemicalsInitialized=.false.

  ! File name for the chemical state data.
  type(varying_string) :: chemicalStateFile

  ! Ranges of the tabulations.
  double precision     :: metallicityMinimum,metallicityMaximum,temperatureMinimum,temperatureMaximum

  ! Extrapolation methods.
  integer              :: extrapolateTemperatureLow,extrapolateTemperatureHigh,extrapolateMetallicityLow,extrapolateMetallicityHigh

  ! The chemical state tables.
  logical                                       :: logarithmicTable,firstMetallicityIsZero,gotHydrogenAtomic,gotHydrogenCation
  integer                                       :: chemicalStateTemperatureNumberPoints,chemicalStateMetallicityNumberPoints
  double precision                              :: firstNonZeroMetallicity
  double precision, allocatable, dimension(:)   :: chemicalStateMetallicities,chemicalStateTemperatures
  double precision, allocatable, dimension(:,:) :: electronDensityTable,hydrogenAtomicDensityTable,hydrogenCationDensityTable

  ! Interpolation structures.
  logical                                       :: resetTemperature=.true., resetMetallicity=.true.
  type(fgsl_interp_accel)                       :: interpolationAcceleratorTemperature,interpolationAcceleratorMetallicity

  ! Chemical indices.
  integer                                       :: electronChemicalIndex,atomicHydrogenChemicalIndex,atomicHydrogenCationChemicalIndex

  ! Current file format version for CIE cooling files.
  integer         , parameter                   :: fileFormatVersionCurrent=1

contains

  integer function Chemical_State_CIE_File_Format_Version()
    !% Return the current file format version of CIE chemical state files.
    implicit none

    Chemical_State_CIE_File_Format_Version=fileFormatVersionCurrent
    return
  end function Chemical_State_CIE_File_Format_Version

  !# <chemicalStateMethod>
  !#  <unitName>Chemical_State_CIE_File_Initialize</unitName>
  !# </chemicalStateMethod>
  subroutine Chemical_State_CIE_File_Initialize(chemicalStateMethod,Electron_Density_Get&
       &,Electron_Density_Temperature_Log_Slope_Get,Electron_Density_Density_Log_Slope_Get,Chemical_Densities_Get)
    !% Initializes the ``CIE ionization state from file'' module.
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: chemicalStateMethod
    procedure(double precision), pointer, intent(inout) :: Electron_Density_Get,Electron_Density_Temperature_Log_Slope_Get&
         &,Electron_Density_Density_Log_Slope_Get
    procedure(),                 pointer, intent(inout) :: Chemical_Densities_Get

    if (chemicalStateMethod == 'cieFromFile') then
       ! Set the procedure pointer.
       Electron_Density_Get                       => Electron_Density_CIE_File
       Electron_Density_Temperature_Log_Slope_Get => Electron_Density_Temperature_Log_Slope_CIE_File
       Electron_Density_Density_Log_Slope_Get     => Electron_Density_Density_Log_Slope_CIE_File
       Chemical_Densities_Get                    => Chemical_Densities_CIE_FIle

       !@ <inputParameter>
       !@   <name>chemicalStateFile</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file containing a tabulation of the collisional chemical equilibrium chemical state.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('chemicalStateFile',chemicalStateFile)
    end if
    return
  end subroutine Chemical_State_CIE_File_Initialize

  subroutine Chemical_State_CIE_File_Read_Initialize
    !% Ensure that the cooling data file has been read in.
    implicit none

    if (.not.chemicalStateInitialized) then
       !$omp critical (Chemical_State_CIE_File_Initialize)
       if (.not.chemicalStateInitialized) then
          
          ! Call routine to read in the tabulated data.
          call Chemical_State_CIE_File_Read(chemicalStateFile)
          
          ! Get chemical indices.
          call Chemical_State_CIE_Chemicals_Initialize
          
          ! Flag that chemical state is now initialized.
          chemicalStateInitialized=.true.
       end if
       !$omp end critical (Chemical_State_CIE_File_Initialize)
    end if
    return
  end subroutine Chemical_State_CIE_File_Read_Initialize
 
  subroutine Chemical_State_CIE_Chemicals_Initialize
    !% Ensure that chemical indices have been found.
    use Chemical_Abundances_Structure
    implicit none
    
    if (.not.chemicalStateChemicalsInitialized) then
       !$omp critical (Chemical_State_CIE_File_Chemicals_Initialize)
       if (.not.chemicalStateChemicalsInitialized) then
          
          ! Get chemical indices.
          electronChemicalIndex               =Chemicals_Index("Electron"            )
          if (gotHydrogenAtomic) then
             atomicHydrogenChemicalIndex      =Chemicals_Index("AtomicHydrogen"      )
          else
             atomicHydrogenChemicalIndex      =-1
          end if
          if (gotHydrogenCation) then
             atomicHydrogenCationChemicalIndex=Chemicals_Index("AtomicHydrogenCation")
          else
             atomicHydrogenCationChemicalIndex=-1
          end if
          
          ! Flag that chemical state chemical indices are now initialized.
          chemicalStateChemicalsInitialized=.true.
          
       end if
       !$omp end critical (Chemical_State_CIE_File_Chemicals_Initialize)
    end if
    return
  end subroutine Chemical_State_CIE_Chemicals_Initialize

  double precision function Electron_Density_CIE_File(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the electron density by interpolating in tabulated CIE data read from a file.
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Ensure file has been read in.
    call Chemical_State_CIE_File_Read_Initialize

    ! Call routine to interpolate in the tabulated function.
    Electron_Density_CIE_File=Electron_Density_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)

    return
  end function Electron_Density_CIE_File

  double precision function Electron_Density_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)
    !% Compute the chemical state by interpolation in the tabulated data.
    use Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Astronomical
    use IO_XML
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation
    double precision,          save       :: temperaturePrevious=-1.0d0,metallicityPrevious=-1.0d0,electronDensityPrevious
    !$omp threadprivate(temperaturePrevious,metallicityPrevious,electronDensityPrevious)
    integer                               :: iTemperature,iMetallicity
    double precision                      :: temperatureUse,metallicityUse,hTemperature,hMetallicity

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < temperatureMinimum) then
       select case (extrapolateTemperatureLow)
       case (extrapolateZero)
          Electron_Density_CIE_File_Interpolate=0.0d0
          return
       case (extrapolateFixed,extrapolatePowerLaw)
          temperatureUse=temperatureMinimum
       end select
    end if
    if (temperatureUse > temperatureMaximum) then
       select case (extrapolateTemperatureHigh)
       case (extrapolateZero)
          Electron_Density_CIE_File_Interpolate=0.0d0
          return
       case (extrapolateFixed,extrapolatePowerLaw)
          temperatureUse=temperatureMaximum
       end select
    end if

    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(abundances)/metallicitySolar

    if (metallicityUse < metallicityMinimum) then
       select case (extrapolateMetallicityLow)
       case (extrapolateZero)
          Electron_Density_CIE_File_Interpolate=0.0d0
          return
       case (extrapolateFixed)
          metallicityUse=metallicityMinimum
       end select
    end if
    if (metallicityUse > metallicityMaximum) then
       select case (extrapolateMetallicityHigh)
       case (extrapolateZero)
          Electron_Density_CIE_File_Interpolate=0.0d0
          return
       case (extrapolateFixed)
          metallicityUse=metallicityMaximum
       end select
    end if

    ! Check if we need to recompute the chemical state.
    if (temperatureUse /= temperaturePrevious .or. metallicityUse /= metallicityPrevious) then

       ! Get the interpolation.
       call Get_Interpolation(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)

       ! Do the interpolation.
       electronDensityPrevious=Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity,electronDensityTable)

       ! Store the temperature and metallicity for which calculation was performed.
       temperaturePrevious=temperatureUse
       metallicityPrevious=metallicityUse
    end if

    ! Scale to the specified density assuming two-body processes, in which case electron density scales with hydrogen density.
    Electron_Density_CIE_File_Interpolate=electronDensityPrevious*numberDensityHydrogen

    return
  end function Electron_Density_CIE_File_Interpolate

  double precision function Electron_Density_Temperature_Log_Slope_CIE_File(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the logarithmic slope of the electron density with respect to temperature by interpolating in tabulated CIE data
    !% read from a file.
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in)  :: abundances
    type(radiationStructure),  intent(in)  :: radiation

    ! Ensure file has been read in.
    call Chemical_State_CIE_File_Read_Initialize

    ! Call routine to interpolate in the tabulated function.
    Electron_Density_Temperature_Log_Slope_CIE_File=Electron_Density_CIE_File_logTemperature_Interpolate(temperature &
         &,numberDensityHydrogen,abundances,radiation)

    return
  end function Electron_Density_Temperature_Log_Slope_CIE_File

  double precision function Electron_Density_CIE_File_logTemperature_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)
    !% Compute the logarithmic gradient of the electron density with respect to temperature by interpolation in the tabulated data.
    use Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Astronomical
    use IO_XML
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation
    double precision,          save       :: temperaturePrevious=-1.0d0,metallicityPrevious=-1.0d0,electronDensitySlopePrevious
    !$omp threadprivate(temperaturePrevious,metallicityPrevious,electronDensitySlopePrevious)
    integer                               :: iTemperature,iMetallicity
    double precision                      :: temperatureUse,metallicityUse,hTemperature,hMetallicity

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < temperatureMinimum) then
       select case (extrapolateTemperatureLow)
       case (extrapolateZero,extrapolateFixed)
          Electron_Density_CIE_File_logTemperature_Interpolate=0.0d0
          return
       case (extrapolatePowerLaw)
          temperatureUse=temperatureMinimum
       end select
    end if
    if (temperatureUse > temperatureMaximum) then
       select case (extrapolateTemperatureHigh)
       case (extrapolateZero,extrapolateFixed)
          Electron_Density_CIE_File_logTemperature_Interpolate=0.0d0
          return
       case (extrapolatePowerLaw)
          temperatureUse=temperatureMaximum
       end select
    end if

    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(abundances)/metallicitySolar
    if (metallicityUse < metallicityMinimum) then
       select case (extrapolateMetallicityLow)
       case (extrapolateZero)
          Electron_Density_CIE_File_logTemperature_Interpolate=0.0d0
          return
       case (extrapolateFixed)
          metallicityUse=metallicityMinimum
       end select
    end if
    if (metallicityUse > metallicityMaximum) then
       select case (extrapolateMetallicityHigh)
       case (extrapolateZero)
          Electron_Density_CIE_File_logTemperature_Interpolate=0.0d0
          return
       case (extrapolateFixed)
          metallicityUse=metallicityMaximum
       end select
    end if

    ! Check if we need to recompute the cooling function.
    if (temperatureUse /= temperaturePrevious .or. metallicityUse /= metallicityPrevious) then

       ! Get the interpolation.
       call Get_Interpolation(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)

       ! Do the interpolation.
       electronDensitySlopePrevious=( ( electronDensityTable(iTemperature+1,iMetallicity  )                       &
            &                          -electronDensityTable(iTemperature  ,iMetallicity  ))*(1.0d0-hMetallicity) &
            &                        +( electronDensityTable(iTemperature+1,iMetallicity+1)                       & 
            &                          -electronDensityTable(iTemperature  ,iMetallicity+1))*       hMetallicity) &
            & /(chemicalStateTemperatures(iTemperature+1)-chemicalStateTemperatures(iTemperature))

       ! Convert to logarithmic gradient if table was not stored logarithmically.
       if (.not.logarithmicTable) electronDensitySlopePrevious=electronDensitySlopePrevious*temperature&
            &/Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity,electronDensityTable)

       ! Store the temperature and metallicity for which calculation was performed.
       temperaturePrevious=temperatureUse
       metallicityPrevious=metallicityUse

    end if

    ! Return the stored value.
    Electron_Density_CIE_File_logTemperature_Interpolate=electronDensitySlopePrevious

    return
  end function Electron_Density_CIE_File_logTemperature_Interpolate

  double precision function Electron_Density_Density_Log_Slope_CIE_File(temperature,numberDensityHydrogen,abundances&
       &,radiation)
    !% Return the logarithmic slope of the electron density with respect to density assuming atomic CIE as computed by {\sc Cloudy}.
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in)  :: abundances
    type(radiationStructure),  intent(in)  :: radiation

    ! Electron density always scales as total density under CIE conditions.
    Electron_Density_Density_Log_Slope_CIE_File=1.0d0

    return
  end function Electron_Density_Density_Log_Slope_CIE_File

  subroutine Chemical_State_CIE_File_Read(chemicalStateFileToRead,metallicityMaximumTabulated)
    !% Read in data from an chemical state file.
    use Galacticus_Error
    use FoX_dom
    use Memory_Management
    use Numerical_Comparison
    use Galacticus_Display
    use IO_XML
    implicit none
    type(varying_string), intent(in)            :: chemicalStateFileToRead
    double precision,     intent(out), optional :: metallicityMaximumTabulated
    type(Node),           pointer               :: doc,datum,thisChemicalState,metallicityElement,extrapolationElement &
         &,extrapolation,thisTemperature,thisElectronDensity,thisHydrogenAtomicDensity ,thisHydrogenCationDensity,version
    type(NodeList),       pointer               :: temperatureDatumList,electronDatumList,hydrogenAtomicDatumList &
         &,hydrogenCationDatumList,chemicalStateList ,metallicityExtrapolationList,temperatureExtrapolationList,versionList
    integer                                     :: iDatum,ioErr,iChemicalState,iExtrapolation,extrapolationMethod,fileFormatVersion
    double precision                            :: datumValues(1)
    character(len=32)                           :: limitType

    !$omp critical (FoX_DOM_Access)

    ! Parse the XML file.
    call Galacticus_Display_Indent('Parsing file: '//chemicalStateFileToRead,verbosityDebug)
    doc => parseFile(char(chemicalStateFileToRead),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Chemical_State_CIE_File_Read','Unable to find chemical state file')
 
    ! Check the file format version of the file.
    versionList => getElementsByTagname(doc,"fileFormat")
    version     => item(versionList,0)
    call extractDataContent(version,fileFormatVersion)
    if (fileFormatVersion /= fileFormatVersionCurrent) call Galacticus_Error_Report('Chemical_State_CIE_File_Read','file format version is out of date')

    ! Get a list of all <ionizationState> elements.
    chemicalStateList => getElementsByTagname(doc,"ionizationState")
    chemicalStateMetallicityNumberPoints=getLength(chemicalStateList)

    ! Extract data from first chemical state and count number of temperatures present.
    thisChemicalState  => item(chemicalStateList,0)
    thisTemperature      => item(getElementsByTagname(thisChemicalState,"temperature"),0)
    temperatureDatumList => getElementsByTagname(thisTemperature,"datum")
    chemicalStateTemperatureNumberPoints=getLength(temperatureDatumList)

    ! Allocate space for the table.
    if (allocated(chemicalStateMetallicities)) call Dealloc_Array(chemicalStateMetallicities)
    if (allocated(chemicalStateTemperatures )) call Dealloc_Array(chemicalStateTemperatures )
    if (allocated(electronDensityTable      )) call Dealloc_Array(electronDensityTable      )
    if (allocated(hydrogenAtomicDensityTable)) call Dealloc_Array(hydrogenAtomicDensityTable)
    if (allocated(hydrogenCationDensityTable)) call Dealloc_Array(hydrogenCationDensityTable)
    call Alloc_Array(chemicalStateMetallicities                                      ,[chemicalStateMetallicityNumberPoints])
    call Alloc_Array(chemicalStateTemperatures ,[chemicalStateTemperatureNumberPoints                                      ])
    call Alloc_Array(electronDensityTable      ,[chemicalStateTemperatureNumberPoints,chemicalStateMetallicityNumberPoints])

    ! Allocate space for atomic hydrogen density, if such data is included.
    hydrogenAtomicDatumList => getElementsByTagname(doc,"hiDensity")
    if (getLength(hydrogenAtomicDatumList) > 0) then
       call Alloc_Array(hydrogenAtomicDensityTable,[chemicalStateTemperatureNumberPoints,chemicalStateMetallicityNumberPoints])
       gotHydrogenAtomic=.true.
    else
       gotHydrogenAtomic=.false.
    end if
    ! Allocate space for ionized hydrogen density, if such data is included.
    hydrogenCationDatumList => getElementsByTagname(doc,"hiiDensity")
    if (getLength(hydrogenCationDatumList) > 0) then
       call Alloc_Array(hydrogenCationDensityTable ,[chemicalStateTemperatureNumberPoints,chemicalStateMetallicityNumberPoints])
       gotHydrogenCation=.true.
    else
       gotHydrogenCation=.false.
    end if

    ! Extract data from the chemical states and populate metallicity and temperature arrays.
    do iChemicalState=0,getLength(chemicalStateList)-1
       ! Get required chemical state.
       thisChemicalState => item(chemicalStateList,iChemicalState)
       ! Extract the metallicity from the <metallicity> element.
       metallicityElement  => item(getElementsByTagname(thisChemicalState,"metallicity"),0)
       call extractDataContent(metallicityElement,chemicalStateMetallicities(iChemicalState+1))
       ! Extract the data.
       thisTemperature           => item(getElementsByTagname(thisChemicalState           ,"temperature"    ),0)
       temperatureDatumList      =>      getElementsByTagname(thisTemperature             ,"datum"          )
       thisElectronDensity       => item(getElementsByTagname(thisChemicalState           ,"electronDensity"),0)
       electronDatumList         =>      getElementsByTagname(thisElectronDensity         ,"datum"          )
       if (gotHydrogenAtomic) then
          thisHydrogenAtomicDensity => item(getElementsByTagname(thisChemicalState        ,"hiDensity" ),0)
          hydrogenAtomicDatumList   =>      getElementsByTagname(thisHydrogenAtomicDensity,"datum"     )
       end if
       if (gotHydrogenCation) then
          thisHydrogenCationDensity => item(getElementsByTagname(thisChemicalState        ,"hiiDensity"),0)
          hydrogenCationDatumList   =>      getElementsByTagname(thisHydrogenCationDensity,"datum"     )
       end if
       ! Check that number of temperatures is consistent.
       if (getLength(temperatureDatumList) /= chemicalStateTemperatureNumberPoints) call&
            & Galacticus_Error_Report('Chemical_State_CIE_File_Read','sizes of temperatures grids must be the same for all&
            & metallicities')
       ! Check that number of cooling rates matches number of temperatures.
       if (getLength(temperatureDatumList) /= getLength(electronDatumList)) call&
            & Galacticus_Error_Report('Chemical_State_CIE_File_Read','sizes of temperature and electron density arrays must match')
       ! Extract data.
       do iDatum=0,getLength(temperatureDatumList)-1
          datum => item(temperatureDatumList,iDatum)
          call extractDataContent(datum,datumValues)
          if (iChemicalState == 0) then
             ! Store the temperature.
             chemicalStateTemperatures(iDatum+1)=datumValues(1)
          else
             ! Check that temperature grids are aligned.
             if (Values_Differ(chemicalStateTemperatures(iDatum+1),datumValues(1),relTol=1.0d-6)) call&
                  & Galacticus_Error_Report('Chemical_State_CIE_File_Read','temperature grids mismatch')
          end if
          ! Store the chemical state.
          datum => item(electronDatumList,iDatum)
          call extractDataContent(datum,datumValues)
          electronDensityTable(iDatum+1,iChemicalState+1)=datumValues(1)
          if (gotHydrogenAtomic) then
             datum => item(hydrogenAtomicDatumList,iDatum)
             call extractDataContent(datum,datumValues)
             hydrogenAtomicDensityTable(iDatum+1,iChemicalState+1)=datumValues(1)
          end if
          if (gotHydrogenCation) then
             datum => item(hydrogenCationDatumList,iDatum)
             call extractDataContent(datum,datumValues)
             hydrogenCationDensityTable(iDatum+1,iChemicalState+1)=datumValues(1)
          end if
       end do
    end do
    where (chemicalStateMetallicities>-999.0d0)
       chemicalStateMetallicities=10.0d0**chemicalStateMetallicities
    elsewhere
       chemicalStateMetallicities=0.0d0
    end where

    ! Extract extrapolation methods from the file.
    extrapolationElement => item(getElementsByTagname(doc,"extrapolation"),0)
    metallicityExtrapolationList => getElementsByTagname(extrapolationElement,"metallicity")
    do iExtrapolation=0,getLength(metallicityExtrapolationList)-1
       extrapolation => item(metallicityExtrapolationList,iExtrapolation)
       call XML_Extrapolation_Element_Decode(extrapolation,limitType,extrapolationMethod,allowedMethods=[extrapolateZero,extrapolateFixed,extrapolatePowerLaw])
       select case (trim(limitType))
       case ('low')
          extrapolateMetallicityLow=extrapolationMethod
       case ('high')
          extrapolateMetallicityHigh=extrapolationMethod
       case default
          call Galacticus_Error_Report('Chemical_State_CIE_File_Read','unrecognized extrapolation limit')
       end select
    end do
    temperatureExtrapolationList => getElementsByTagname(extrapolationElement,"temperature")
    do iExtrapolation=0,getLength(temperatureExtrapolationList)-1
       extrapolation => item(temperatureExtrapolationList,iExtrapolation)
       call XML_Extrapolation_Element_Decode(extrapolation,limitType,extrapolationMethod,allowedMethods=[extrapolateZero,extrapolateFixed,extrapolatePowerLaw])
       select case (trim(limitType))
       case ('low')
          extrapolateTemperatureLow=extrapolationMethod
       case ('high')
          extrapolateTemperatureHigh=extrapolationMethod
       case default
          call Galacticus_Error_Report('Chemical_State_CIE_File_Read','unrecognized extrapolation limit')
       end select
    end do
    ! Destroy the document.
    call destroy(doc)
    call Galacticus_Display_Unindent('done',verbosityDebug)
    !$omp end critical (FoX_DOM_Access)

    ! Store table ranges for convenience.
    metallicityMinimum=chemicalStateMetallicities(1)
    metallicityMaximum=chemicalStateMetallicities(chemicalStateMetallicityNumberPoints)
    temperatureMinimum=chemicalStateTemperatures(1)
    temperatureMaximum=chemicalStateTemperatures(chemicalStateTemperatureNumberPoints)

    ! Decide whether or not to make the tables logarithmic.
    logarithmicTable= all(electronDensityTable       > 0.0d0)                              &
         &             .and.                                                               &
         &           (all(hydrogenAtomicDensityTable > 0.0d0) .or. .not.gotHydrogenAtomic) &
         &             .and.                                                               &
         &           (all(hydrogenCationDensityTable > 0.0d0) .or. .not.gotHydrogenCation)
    if (logarithmicTable) then
       firstMetallicityIsZero=(chemicalStateMetallicities(1) == 0.0d0)
       if (firstMetallicityIsZero) firstNonZeroMetallicity=chemicalStateMetallicities(2)
       where (chemicalStateMetallicities > 0.0d0)
          chemicalStateMetallicities=dlog(chemicalStateMetallicities)
       elsewhere
          chemicalStateMetallicities=-999.0d0
       end where
       chemicalStateTemperatures=dlog(chemicalStateTemperatures)
       electronDensityTable=dlog(electronDensityTable)
       if (gotHydrogenAtomic) hydrogenAtomicDensityTable=dlog(hydrogenAtomicDensityTable)
       if (gotHydrogenCation) hydrogenCationDensityTable=dlog(hydrogenCationDensityTable)
    else
       if (        extrapolateTemperatureLow  == extrapolatePowerLaw   &
            & .or. extrapolateTemperatureHigh == extrapolatePowerLaw   &
            & .or. extrapolateMetallicityLow  == extrapolatePowerLaw   &
            & .or. extrapolateMetallicityHigh == extrapolatePowerLaw ) &
            & call Galacticus_Error_Report('Chemical_State_CIE_File_Read','power law extrapolation allowed only in loggable tables')
    end if

    ! Force interpolation accelerators to be reset.
    resetTemperature=.true.
    resetMetallicity=.true.

    ! Return the maximum metallicity if requested.
    if (present(metallicityMaximumTabulated)) metallicityMaximumTabulated=metallicityMaximum

    return
  end subroutine Chemical_State_CIE_File_Read

  subroutine Get_Interpolation(temperatureIn,metallicityIn,iTemperature,hTemperature,iMetallicity,hMetallicity)
    !% Determine the interpolating paramters.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in)  :: temperatureIn,metallicityIn
    integer,          intent(out) :: iTemperature,iMetallicity
    double precision, intent(out) :: hTemperature,hMetallicity
    double precision              :: temperatureUse,metallicityUse

    ! Copy the input parameters.
    temperatureUse=temperatureIn
    metallicityUse=max(metallicityIn,0.0d0)
    ! Get interpolation in temperature.
    if (logarithmicTable) temperatureUse=dlog(temperatureUse)
    iTemperature=Interpolate_Locate(chemicalStateTemperatureNumberPoints,chemicalStateTemperatures&
         &,interpolationAcceleratorTemperature,temperatureUse,resetTemperature)
    iTemperature=max(min(iTemperature,chemicalStateTemperatureNumberPoints),1)
    hTemperature=(temperatureUse-chemicalStateTemperatures(iTemperature))/(chemicalStateTemperatures(iTemperature+1)&
         &-chemicalStateTemperatures(iTemperature))

    ! Get interpolation in metallicity.
    if (firstMetallicityIsZero.and.metallicityUse < firstNonZeroMetallicity) then
       iMetallicity=1
       hMetallicity=metallicityUse/firstNonZeroMetallicity
    else
       if (logarithmicTable) metallicityUse=dlog(metallicityUse)
       iMetallicity=Interpolate_Locate(chemicalStateMetallicityNumberPoints,chemicalStateMetallicities&
            &,interpolationAcceleratorMetallicity,metallicityUse,resetMetallicity)
       iMetallicity=max(min(iMetallicity,chemicalStateMetallicityNumberPoints),1)
       hMetallicity=(metallicityUse-chemicalStateMetallicities(iMetallicity))/(chemicalStateMetallicities(iMetallicity+1)&
            &-chemicalStateMetallicities(iMetallicity))
    end if
    return
  end subroutine Get_Interpolation

  double precision function Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity,densityTable)
    !% Perform the interpolation.
    implicit none
    integer,          intent(in)                 :: iTemperature,iMetallicity
    double precision, intent(in)                 :: hTemperature,hMetallicity
    double precision, intent(in), dimension(:,:) :: densityTable

    ! Do the interpolation.
    Do_Interpolation=densityTable(iTemperature  ,iMetallicity  )*(1.0d0-hTemperature)*(1.0d0-hMetallicity)&
         &          +densityTable(iTemperature  ,iMetallicity+1)*(1.0d0-hTemperature)*       hMetallicity &
         &          +densityTable(iTemperature+1,iMetallicity  )*       hTemperature *(1.0d0-hMetallicity)&
         &          +densityTable(iTemperature+1,iMetallicity+1)*       hTemperature *       hMetallicity

    ! Exponentiate the result if the table was stored as the log.
    if (logarithmicTable) Do_Interpolation=dexp(Do_Interpolation)

    return
  end function Do_Interpolation

  subroutine Chemical_Densities_CIE_File(theseAbundances,temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the densities of chemical species at the given temperature and hydrogen density for the specified set of abundances
    !% and radiation field. Units of the returned electron density are cm$^-3$.
    use Abundances_Structure
    use Radiation_Structure
    use Chemical_Abundances_Structure
    implicit none
    type(chemicalAbundancesStructure), intent(inout) :: theseAbundances
    double precision,                  intent(in)    :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)    :: abundances
    type(radiationStructure),          intent(in)    :: radiation

    ! Ensure file has been read in.
    call Chemical_State_CIE_File_Read_Initialize

    ! Call routine to interpolate in the tabulated function.
    call Chemical_Densities_CIE_File_Interpolate(theseAbundances,temperature,numberDensityHydrogen,abundances,radiation)

    return
  end subroutine Chemical_Densities_CIE_File

  subroutine Chemical_Densities_CIE_File_Interpolate(theseChemicals,temperature,numberDensityHydrogen,abundances,radiation)
    !% Compute the chemical state by interpolation in the tabulated data.
    use Abundances_Structure
    use Radiation_Structure
    use Chemical_Abundances_Structure
    use Numerical_Constants_Astronomical
    use IO_XML
    implicit none
    type(chemicalAbundancesStructure), intent(inout) :: theseChemicals
    double precision,                  intent(in)    :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)    :: abundances
    type(radiationStructure),          intent(in)    :: radiation
    double precision,                  save          :: temperaturePrevious=-1.0d0,metallicityPrevious=-1.0d0
    type(chemicalAbundancesStructure), save          :: chemicalDensitiesPrevious
    !$omp threadprivate(temperaturePrevious,metallicityPrevious,chemicalDensitiesPrevious)
    integer                                          :: iTemperature,iMetallicity
    double precision                                 :: temperatureUse,metallicityUse,hTemperature,hMetallicity

    ! Ensure that chemical indices have been determined.
    call Chemical_State_CIE_Chemicals_Initialize

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < temperatureMinimum) then
       select case (extrapolateTemperatureLow)
       case (extrapolateZero)
          call theseChemicals%reset()
          return
       case (extrapolateFixed,extrapolatePowerLaw)
          temperatureUse=temperatureMinimum
       end select
    end if
    if (temperatureUse > temperatureMaximum) then
       select case (extrapolateTemperatureHigh)
       case (extrapolateZero)
          call theseChemicals%reset()
          return
       case (extrapolateFixed,extrapolatePowerLaw)
          temperatureUse=temperatureMaximum
       end select
    end if

    ! Handle out of range metallicities.
    metallicityUse=Abundances_Get_Metallicity(abundances)/metallicitySolar

    if (metallicityUse < metallicityMinimum) then
       select case (extrapolateMetallicityLow)
       case (extrapolateZero)
          call theseChemicals%reset()
          return
       case (extrapolateFixed)
          metallicityUse=metallicityMinimum
       end select
    end if
    if (metallicityUse > metallicityMaximum) then
       select case (extrapolateMetallicityHigh)
       case (extrapolateZero)
          call theseChemicals%reset()
          return
       case (extrapolateFixed)
          metallicityUse=metallicityMaximum
       end select
    end if

    ! Check if we need to recompute the chemical state.
    if (temperatureUse /= temperaturePrevious .or. metallicityUse /= metallicityPrevious) then

       ! Get the interpolation.
       call Get_Interpolation(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)

       ! Reset densities to zero.
       call chemicalDensitiesPrevious%reset()

       ! Do the interpolation.
       if (electronChemicalIndex             > 0) call chemicalDensitiesPrevious%abundanceSet(electronChemicalIndex             &
            & ,Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity,electronDensityTable      ))
       if (atomicHydrogenChemicalIndex       > 0) call chemicalDensitiesPrevious%abundanceSet(atomicHydrogenChemicalIndex       &
            & ,Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity,hydrogenAtomicDensityTable))
       if (atomicHydrogenCationChemicalIndex > 0) call chemicalDensitiesPrevious%abundanceSet(atomicHydrogenCationChemicalIndex &
            & ,Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity,hydrogenCationDensityTable))

       ! Store the temperature and metallicity for which calculation was performed.
       temperaturePrevious=temperatureUse
       metallicityPrevious=metallicityUse
    end if

    ! Scale to the specified density assuming two-body processes, in which case densities scale with hydrogen density.
    theseChemicals=chemicalDensitiesPrevious
    call theseChemicals%multiply(numberDensityHydrogen)

    return
  end subroutine Chemical_Densities_CIE_File_Interpolate

end module Chemical_States_CIE_File

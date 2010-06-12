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








!% Contains a module which reads and interpolates a collisional ionization equilibrium ionization state from a file.

module Ionization_States_CIE_File
  !% Reads and interpolates a collisional ionization equilibrium ionization state from a file.
  use ISO_Varying_String
  use FGSL
  private
  public :: Ionization_State_CIE_File_Initialize, Ionization_State_CIE_File_Read, Electron_Density_CIE_File_Interpolate,&
       & Electron_Density_CIE_File_logTemperature_Interpolate
  
  ! Flag to indicate if this module has been initialized.
  logical              :: ionizationStateInitialized=.false.

  ! File name for the ionization state data.
  type(varying_string) :: ionizationStateFile

  ! Ranges of the tabulations.
  double precision     :: metallicityMinimum,metallicityMaximum,temperatureMinimum,temperatureMaximum

  ! Extrapolation methods.
  integer, parameter :: extrapolateZero=0, extrapolateFixed=1, extrapolatePowerLaw=2
  integer            :: extrapolateTemperatureLow,extrapolateTemperatureHigh,extrapolateMetallicityLow,extrapolateMetallicityHigh

  ! The ionization state tables.
  logical                                       :: logarithmicTable,firstMetallicityIsZero
  integer                                       :: ionizationStateTemperatureNumberPoints,ionizationStateMetallicityNumberPoints
  double precision                              :: firstNonZeroMetallicity
  double precision, allocatable, dimension(:)   :: ionizationStateMetallicities,ionizationStateTemperatures
  double precision, allocatable, dimension(:,:) :: electronDensityTable

  ! Interpolation structures.
  logical                                       :: resetTemperature=.true., resetMetallicity=.true.
  type(fgsl_interp_accel)                       :: interpolationAcceleratorTemperature,interpolationAcceleratorMetallicity

contains
  
  !# <ionizationStateMethod>
  !#  <unitName>Ionization_State_CIE_File_Initialize</unitName>
  !# </ionizationStateMethod>
  subroutine Ionization_State_CIE_File_Initialize(ionizationStateMethod,Electron_Density_Get&
       &,Electron_Density_Temperature_Log_Slope_Get,Electron_Density_Density_Log_Slope_Get)
    !% Initializes the ``CIE ionization state from file'' module.
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: ionizationStateMethod
    procedure(),          pointer, intent(inout) :: Electron_Density_Get,Electron_Density_Temperature_Log_Slope_Get&
         &,Electron_Density_Density_Log_Slope_Get

    if (ionizationStateMethod == 'CIE_from_file') then
       ! Set the procedure pointer.
       Electron_Density_Get                       => Electron_Density_CIE_File
       Electron_Density_Temperature_Log_Slope_Get => Electron_Density_Temperature_Log_Slope_CIE_File
       Electron_Density_Density_Log_Slope_Get     => Electron_Density_Density_Log_Slope_CIE_File

       !@ <inputParameter>
       !@   <name>ionizationStateFile</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file containing a tabulation of the collisional ionization equilibrium ionization state.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('ionizationStateFile',ionizationStateFile)
    end if
    return
  end subroutine Ionization_State_CIE_File_Initialize

  subroutine Ionization_State_CIE_File_Read_Initialize
    !% Ensure that the cooling data file has been read in.
    implicit none
    
    !$omp critical (Ionization_State_CIE_File_Initialize)
    if (.not.ionizationStateInitialized) then
       
       ! Call routine to read in the tabulated data.
       call Ionization_State_CIE_File_Read(ionizationStateFile)
       
       ! Flag that ionization state is now initialized.
       ionizationStateInitialized=.true.
    end if
    !$omp end critical (Ionization_State_CIE_File_Initialize)
    return
  end subroutine Ionization_State_CIE_File_Read_Initialize

  double precision function Electron_Density_CIE_File(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the electron density by interpolating in tabulated CIE data read from a file.
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Ensure file has been read in.
    call Ionization_State_CIE_File_Read_Initialize
    
    ! Call routine to interpolate in the tabulated function.
    Electron_Density_CIE_File=Electron_Density_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)

    return
  end function Electron_Density_CIE_File
  
  double precision function Electron_Density_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)
    !% Compute the ionization state by interpolation in the tabulated data.
    use Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Astronomical
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

    ! Check if we need to recompute the ionization state.
    if (temperatureUse /= temperaturePrevious .or. metallicityUse /= metallicityPrevious) then

       ! Get the interpolation.
       call Get_Interpolation(temperatureUse,metallicityUse,iTemperature,hTemperature,iMetallicity,hMetallicity)
       
       ! Do the interpolation.
       electronDensityPrevious=Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity)
       
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
    call Ionization_State_CIE_File_Read_Initialize
       
    ! Call routine to interpolate in the tabulated function.
    Electron_Density_Temperature_Log_Slope_CIE_File=Electron_Density_CIE_File_logTemperature_Interpolate(temperature &
         &,numberDensityHydrogen,abundances,radiation)
    
    return
  end function Electron_Density_Temperature_Log_Slope_CIE_File

  double precision function Electron_Density_CIE_File_logTemperature_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)
    !% Compute the logarithmic gradient of the electron density with respect to temperature by interpolation in the tabulated data.
    use Abundances_Structure
    use Radiation_Structure
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
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
            & /(ionizationStateTemperatures(iTemperature+1)-ionizationStateTemperatures(iTemperature))
       
       ! Convert to logarithmic gradient if table was not stored logarithmically.
       if (.not.logarithmicTable) electronDensitySlopePrevious=electronDensitySlopePrevious*temperature&
            &/Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity)
       
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
      
  subroutine Ionization_State_CIE_File_Read(ionizationStateFileToRead,metallicityMaximumTabulated)
    !% Read in data from an ionization state file.
    use Galacticus_Error
    use FoX_dom
    use Memory_Management
    use Numerical_Comparison
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)            :: ionizationStateFileToRead
    double precision,     intent(out), optional :: metallicityMaximumTabulated
    type(Node),           pointer               :: doc,datum,thisIonizationState,metallicityElement,extrapolationElement&
         &,extrapolation,limitElement,methodElement,thisTemperature,thisElectronDensity
    type(NodeList),       pointer               :: temperatureDatumList,ionizationDatumList,ionizationStateList&
         &,metallicityExtrapolationList,temperatureExtrapolationList
    integer                                     :: iDatum,ioErr,iIonizationState,iExtrapolation,extrapolationMethod
    double precision                            :: datumValues(1)
    character(len=32)                           :: limitType,methodType

    !$omp critical (FoX_DOM_Access)

    ! Parse the XML file.
    call Galacticus_Display_Indent('Parsing file: '//ionizationStateFileToRead,3)
    doc => parseFile(char(ionizationStateFileToRead),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Ionization_State_CIE_File_Read','Unable to find ionization state file')

    ! Get a list of all <ionizationState> elements.
    ionizationStateList => getElementsByTagname(doc,"ionizationState")
    ionizationStateMetallicityNumberPoints=getLength(ionizationStateList)

    ! Extract data from first ionization state and count number of temperatures present.
    thisIonizationState  => item(ionizationStateList,0)
    thisTemperature      => item(getElementsByTagname(thisIonizationState,"temperature"),0)
    temperatureDatumList => getElementsByTagname(thisIonizationState,"datum")
    ionizationStateTemperatureNumberPoints=getLength(temperatureDatumList)

    ! Allocate space for the table.
    if (allocated(ionizationStateMetallicities)) call Dealloc_Array(ionizationStateMetallicities)
    if (allocated(ionizationStateTemperatures )) call Dealloc_Array(ionizationStateTemperatures )
    if (allocated(electronDensityTable        )) call Dealloc_Array(electronDensityTable        )
    call Alloc_Array(ionizationStateMetallicities                                       ,ionizationStateMetallicityNumberPoints&
         &,'ionizationStateMetallicities')
    call Alloc_Array(ionizationStateTemperatures ,ionizationStateTemperatureNumberPoints                                       &
         &,'ionizationStateTemperature'  )
    call Alloc_Array(electronDensityTable        ,ionizationStateTemperatureNumberPoints,ionizationStateMetallicityNumberPoints&
         &,'electronDensityTable'        )

    ! Extract data from the ionization states and populate metallicity and temperature arrays.
    do iIonizationState=0,getLength(ionizationStateList)-1
       ! Get required ionization state.
       thisIonizationState => item(ionizationStateList,iIonizationState)
       ! Extract the metallicity from the <metallicity> element.
       metallicityElement  => item(getElementsByTagname(thisIonizationState,"metallicity"),0)
       call extractDataContent(metallicityElement,ionizationStateMetallicities(iIonizationState+1))
       ! Extract the data.
       thisTemperature      => item(getElementsByTagname(thisIonizationState,"temperature"),0)
       temperatureDatumList =>      getElementsByTagname(thisTemperature    ,"datum"      )
       thisElectronDensity  => item(getElementsByTagname(thisIonizationState,"coolingRate"),0)
       ionizationDatumList  =>      getElementsByTagname(thisElectronDensity,"datum"      )
       ! Check that number of temperatures is consistent.
       if (getLength(temperatureDatumList) /= ionizationStateTemperatureNumberPoints) call&
            & Galacticus_Error_Report('Ionization_State_CIE_File_Read','sizes of temperatures grids must be the same for all&
            & metallicities')
       ! Check that number of cooling rates matches number of temperatures.
       if (getLength(temperatureDatumList) /= getLength(ionizationDatumList)) call&
            & Galacticus_Error_Report('Ionization_State_CIE_File_Read','sizes of temperature and electron density arrays must match')
       ! Extract data.
       do iDatum=0,getLength(temperatureDatumList)-1
          datum => item(temperatureDatumList,iDatum)
          call extractDataContent(datum,datumValues)
          if (iIonizationState == 0) then
             ! Store the temperature.
             ionizationStateTemperatures(iDatum+1)=datumValues(1)
          else
             ! Check that temperature grids are aligned.
             if (Values_Differ(ionizationStateTemperatures(iDatum+1),datumValues(1),relTol=1.0d-6)) call&
                  & Galacticus_Error_Report('Ionization_State_CIE_File_Read','temperature grids mismatch')
          end if
          ! Store the ionization state.
          datum => item(ionizationDatumList,iDatum)
          call extractDataContent(datum,datumValues)
          electronDensityTable(iDatum+1,iIonizationState+1)=datumValues(1)
       end do
    end do
    where (ionizationStateMetallicities>-999.0d0)
       ionizationStateMetallicities=10.0d0**ionizationStateMetallicities
    elsewhere
       ionizationStateMetallicities=0.0d0
    end where

    ! Extract extrapolation methods from the file.
    extrapolationElement => item(getElementsByTagname(doc,"extrapolation"),0)
    metallicityExtrapolationList => getElementsByTagname(extrapolationElement,"metallicity")
    do iExtrapolation=0,getLength(metallicityExtrapolationList)-1
       extrapolation => item(metallicityExtrapolationList,iExtrapolation)
       limitElement  => item(getElementsByTagname(extrapolation,"limit"),0)
       call extractDataContent(limitElement,limitType)
       methodElement  => item(getElementsByTagname(extrapolation,"method"),0)
       call extractDataContent(methodElement,methodType)
       select case (trim(methodType))
       case ('zero')
          extrapolationMethod=extrapolateZero
       case ('fixed')
          extrapolationMethod=extrapolateFixed
       case ('power law')
          extrapolationMethod=extrapolatePowerLaw
       case default
          call Galacticus_Error_Report('Ionization_State_CIE_File_Read','unrecognized extrapolation method')
       end select
       select case (trim(limitType))
       case ('low')
          extrapolateMetallicityLow=extrapolationMethod
       case ('high')
          extrapolateMetallicityHigh=extrapolationMethod
       case default
         call Galacticus_Error_Report('Ionization_State_CIE_File_Read','unrecognized extrapolation limit')
       end select
    end do
    temperatureExtrapolationList => getElementsByTagname(extrapolationElement,"temperature")
    do iExtrapolation=0,getLength(temperatureExtrapolationList)-1
       extrapolation => item(temperatureExtrapolationList,iExtrapolation)
       limitElement  => item(getElementsByTagname(extrapolation,"limit"),0)
       call extractDataContent(limitElement,limitType)
       methodElement  => item(getElementsByTagname(extrapolation,"method"),0)
       call extractDataContent(methodElement,methodType)
       select case (trim(methodType))
       case ('zero')
          extrapolationMethod=extrapolateZero
       case ('fixed')
          extrapolationMethod=extrapolateFixed
       case ('power law')
          extrapolationMethod=extrapolatePowerLaw
       case default
          call Galacticus_Error_Report('Ionization_State_CIE_File_Read','unrecognized extrapolation method')
       end select
       select case (trim(limitType))
       case ('low')
          extrapolateTemperatureLow=extrapolationMethod
       case ('high')
          extrapolateTemperatureHigh=extrapolationMethod
       case default
          call Galacticus_Error_Report('Ionization_State_CIE_File_Read','unrecognized extrapolation limit')
       end select
    end do
    ! Destroy the document.
    call destroy(doc)
    call Galacticus_Display_Unindent('done',3)
    !$omp end critical (FoX_DOM_Access)
  
    ! Store table ranges for convenience.
    metallicityMinimum=ionizationStateMetallicities(1)
    metallicityMaximum=ionizationStateMetallicities(ionizationStateMetallicityNumberPoints)
    temperatureMinimum=ionizationStateTemperatures(1)
    temperatureMaximum=ionizationStateTemperatures(ionizationStateTemperatureNumberPoints)

    ! Decide whether or not to make the tables logarithmic.
    logarithmicTable=all(electronDensityTable > 0.0d0)
    if (logarithmicTable) then
       firstMetallicityIsZero=(ionizationStateMetallicities(1) == 0.0d0)
       if (firstMetallicityIsZero) firstNonZeroMetallicity=ionizationStateMetallicities(2)
       where (ionizationStateMetallicities > 0.0d0)
          ionizationStateMetallicities=dlog(ionizationStateMetallicities)
       elsewhere
          ionizationStateMetallicities=-999.0d0
       end where
       ionizationStateTemperatures=dlog(ionizationStateTemperatures)
       electronDensityTable=dlog(electronDensityTable)
    else
       if (        extrapolateTemperatureLow  == extrapolatePowerLaw   &
            & .or. extrapolateTemperatureHigh == extrapolatePowerLaw   &
            & .or. extrapolateMetallicityLow  == extrapolatePowerLaw   &
            & .or. extrapolateMetallicityHigh == extrapolatePowerLaw ) &
            & call Galacticus_Error_Report('Ionization_State_CIE_File_Read','power law extrapolation allowed only in loggable tables')
    end if

    ! Force interpolation accelerators to be reset.
    resetTemperature=.true.
    resetMetallicity=.true.

    ! Return the maximum metallicity if requested.
    if (present(metallicityMaximumTabulated)) metallicityMaximumTabulated=metallicityMaximum

    return
  end subroutine Ionization_State_CIE_File_Read

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
    iTemperature=Interpolate_Locate(ionizationStateTemperatureNumberPoints,ionizationStateTemperatures&
         &,interpolationAcceleratorTemperature,temperatureUse,resetTemperature)
    iTemperature=max(min(iTemperature,ionizationStateTemperatureNumberPoints),1)
    hTemperature=(temperatureUse-ionizationStateTemperatures(iTemperature))/(ionizationStateTemperatures(iTemperature+1)&
         &-ionizationStateTemperatures(iTemperature))

    ! Get interpolation in metallicity.
    if (firstMetallicityIsZero.and.metallicityUse < firstNonZeroMetallicity) then
       iMetallicity=1
       hMetallicity=metallicityUse/firstNonZeroMetallicity
    else
       if (logarithmicTable) metallicityUse=dlog(metallicityUse)
       iMetallicity=Interpolate_Locate(ionizationStateMetallicityNumberPoints,ionizationStateMetallicities&
         &,interpolationAcceleratorMetallicity,metallicityUse,resetMetallicity)
       iMetallicity=max(min(iMetallicity,ionizationStateMetallicityNumberPoints),1)
       hMetallicity=(metallicityUse-ionizationStateMetallicities(iMetallicity))/(ionizationStateMetallicities(iMetallicity+1)&
            &-ionizationStateMetallicities(iMetallicity))
    end if
    return
  end subroutine Get_Interpolation

  double precision function Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity)
    !% Perform the interpolation.
    implicit none
    integer,          intent(in) :: iTemperature,iMetallicity
    double precision, intent(in) :: hTemperature,hMetallicity

    ! Do the interpolation.
    Do_Interpolation=electronDensityTable(iTemperature  ,iMetallicity  )*(1.0d0-hTemperature)*(1.0d0-hMetallicity)&
         &          +electronDensityTable(iTemperature  ,iMetallicity+1)*(1.0d0-hTemperature)*       hMetallicity &
         &          +electronDensityTable(iTemperature+1,iMetallicity  )*       hTemperature *(1.0d0-hMetallicity)&
         &          +electronDensityTable(iTemperature+1,iMetallicity+1)*       hTemperature *       hMetallicity

    ! Exponentiate the result if the table was stored as the log.
    if (logarithmicTable) Do_Interpolation=dexp(Do_Interpolation)

    return
  end function Do_Interpolation

end module Ionization_States_CIE_File

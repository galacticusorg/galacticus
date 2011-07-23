!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which reads and interpolates a collisional ionization equilibrium cooling function from a file.

module Cooling_Functions_CIE_File
  !% Reads and interpolates a collisional ionization equilibrium cooling function from a file.
  use ISO_Varying_String
  use FGSL
  private
  public :: Cooling_Function_CIE_File_Initialize, Cooling_Function_CIE_File_Read, Cooling_Function_CIE_File_Interpolate,&
       & Cooling_Function_CIE_File_logTemperature_Interpolate, Cooling_Function_CIE_File,&
       & Cooling_Function_Temperature_Slope_CIE_File, Cooling_Function_Density_Slope_CIE_File
  
  ! Flag indicating whether or not this cooling function is selected.
  logical              :: functionSelected=.false.

  ! Flag to indicate if this module has been initialized.
  logical              :: coolingFunctionInitialized=.false.

  ! File name for the cooling function data.
  type(varying_string) :: coolingFunctionFile

  ! Ranges of the tabulations.
  double precision     :: metallicityMinimum,metallicityMaximum,temperatureMinimum,temperatureMaximum

  ! Extrapolation methods.
  integer              :: extrapolateTemperatureLow,extrapolateTemperatureHigh,extrapolateMetallicityLow,extrapolateMetallicityHigh

  ! The cooling function tables.
  logical                                       :: logarithmicTable,firstMetallicityIsZero
  integer                                       :: coolingFunctionTemperatureNumberPoints,coolingFunctionMetallicityNumberPoints
  double precision                              :: firstNonZeroMetallicity
  double precision, allocatable, dimension(:)   :: coolingFunctionMetallicities,coolingFunctionTemperatures
  double precision, allocatable, dimension(:,:) :: coolingFunctionTable

  ! Interpolation structures.
  logical                                       :: resetTemperature=.true., resetMetallicity=.true.
  type(fgsl_interp_accel)                       :: interpolationAcceleratorTemperature,interpolationAcceleratorMetallicity

contains
  
  !# <coolingFunctionMethods>
  !#  <unitName>Cooling_Function_CIE_File_Initialize</unitName>
  !# </coolingFunctionMethods>
  subroutine Cooling_Function_CIE_File_Initialize(coolingFunctionMethods)
    !% Initializes the ``CIE cooling function from file'' module.
    use Input_Parameters
    implicit none
    type(varying_string), intent(in)    :: coolingFunctionMethods(:)
    
    if (any(coolingFunctionMethods == 'CIE_from_file')) then
       ! Flag that this cooling function has been selected.
       functionSelected=.true.
       !@ <inputParameter>
       !@   <name>coolingFunctionFile</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file containing a tabulation of the collisional ionization equilibrium cooling function.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingFunctionFile',coolingFunctionFile)
    end if
    return
  end subroutine Cooling_Function_CIE_File_Initialize

  subroutine Cooling_Function_CIE_File_Read_Initialize
    !% Ensure that the cooling data file has been read in.
    implicit none
    
    !$omp critical (Cooling_Function_CIE_File_Initialize)
    if (.not.coolingFunctionInitialized) then
       
       ! Call routine to read in the tabulated data.
       call Cooling_Function_CIE_File_Read(coolingFunctionFile)
       
       ! Flag that cooling function is now initialized.
       coolingFunctionInitialized=.true.
    end if
    !$omp end critical (Cooling_Function_CIE_File_Initialize)
    return
  end subroutine Cooling_Function_CIE_File_Read_Initialize

  !# <coolingFunctionCompute>
  !#   <unitName>Cooling_Function_CIE_File</unitName>
  !# </coolingFunctionCompute>
  subroutine Cooling_Function_CIE_File(coolingFunction,temperature,numberDensityHydrogen,abundances,chemicalDensities,radiation)
    !% Return the cooling function by interpolating in tabulated CIE data read from a file.
    use Abundances_Structure
    use Radiation_Structure
    use Chemical_Abundances_Structure
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)  :: abundances
    type(chemicalAbundancesStructure), intent(in) :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunction

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       
       ! Ensure file has been read in.
       call Cooling_Function_CIE_File_Read_Initialize
       
       ! Call routine to interpolate in the tabulated function.
       coolingFunction=Cooling_Function_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)

    else
       
       ! Not selected, return zero.
       coolingFunction=0.0d0
       
    end if

    return
  end subroutine Cooling_Function_CIE_File
  
  !# <coolingFunctionTemperatureSlopeCompute>
  !#   <unitName>Cooling_Function_Temperature_Slope_CIE_File</unitName>
  !# </coolingFunctionTemperatureSlopeCompute>
  subroutine Cooling_Function_Temperature_Slope_CIE_File(coolingFunctionTemperatureSlope,temperature,numberDensityHydrogen&
       &,abundances,chemicalDensities,radiation)
    !% Return the slope of the cooling function with respect to temperature by interpolating in tabulated CIE data
    !% read from a file.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)  :: abundances
    type(chemicalAbundancesStructure), intent(in)  :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunctionTemperatureSlope
    double precision                                :: coolingFunction

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       
       ! Ensure file has been read in.
       call Cooling_Function_CIE_File_Read_Initialize
       
       ! Get the cooling function.
       call Cooling_Function_CIE_File(coolingFunction,temperature,numberDensityHydrogen,abundances,chemicalDensities,radiation)
       
       ! Call routine to interpolate in the tabulated function.
       coolingFunctionTemperatureSlope=Cooling_Function_CIE_File_logTemperature_Interpolate(temperature&
            &,numberDensityHydrogen,abundances,radiation)*coolingFunction/temperature

    else
       
       ! Not selected, return zero.
       coolingFunctionTemperatureSlope=0.0d0
       
    end if
       
    return
  end subroutine Cooling_Function_Temperature_Slope_CIE_File
  
  !# <coolingFunctionDensitySlopeCompute>
  !#   <unitName>Cooling_Function_Density_Slope_CIE_File</unitName>
  !# </coolingFunctionDensitySlopeCompute>
  subroutine Cooling_Function_Density_Slope_CIE_File(coolingFunctionDensitySlope,temperature,numberDensityHydrogen,abundances&
       &,chemicalDensities,radiation)
    !% Return the logarithmic slope of the cooling function with respect to density.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)  :: abundances
    type(chemicalAbundancesStructure), intent(in)  :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunctionDensitySlope
    double precision                               :: coolingFunction

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       
       ! Get the cooling function.
       call Cooling_Function_CIE_File(coolingFunction,temperature,numberDensityHydrogen,abundances,chemicalDensities,radiation)
       
       ! Logarithmic slope is always 2 for a CIE cooling function.
       coolingFunctionDensitySlope=2.0d0*coolingFunction/numberDensityHydrogen
       
    else
       
       ! Not selected, return zero.
       coolingFunctionDensitySlope=0.0d0
       
    end if

    return
  end subroutine Cooling_Function_Density_Slope_CIE_File
  
  double precision function Cooling_Function_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)
    !% Compute the cooling function by interpolation in the tabulated data.
    use Abundances_Structure
    use Radiation_Structure
    use Numerical_Constants_Astronomical
    use IO_XML
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation
    double precision,          save       :: temperaturePrevious=-1.0d0,metallicityPrevious=-1.0d0,coolingFunctionPrevious
    !$omp threadprivate(temperaturePrevious,metallicityPrevious,coolingFunctionPrevious)
    integer                               :: iTemperature,iMetallicity
    double precision                      :: temperatureUse,metallicityUse,hTemperature,hMetallicity

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < temperatureMinimum) then
       select case (extrapolateTemperatureLow)
       case (extrapolateZero)
          Cooling_Function_CIE_File_Interpolate=0.0d0
          return
       case (extrapolateFixed,extrapolatePowerLaw)
          temperatureUse=temperatureMinimum
       end select
    end if
    if (temperatureUse > temperatureMaximum) then
       select case (extrapolateTemperatureHigh)
       case (extrapolateZero)
          Cooling_Function_CIE_File_Interpolate=0.0d0
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
          Cooling_Function_CIE_File_Interpolate=0.0d0
          return
       case (extrapolateFixed)
          metallicityUse=metallicityMinimum
       end select
    end if
    if (metallicityUse > metallicityMaximum) then
       select case (extrapolateMetallicityHigh)
       case (extrapolateZero)
          Cooling_Function_CIE_File_Interpolate=0.0d0
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
       coolingFunctionPrevious=Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity)
       
       ! Store the temperature and metallicity for which calculation was performed.
       temperaturePrevious=temperatureUse
       metallicityPrevious=metallicityUse
    end if

    ! Scale to the specified density assuming all processes are proportional to hydrogen density squared.
    Cooling_Function_CIE_File_Interpolate=coolingFunctionPrevious*(numberDensityHydrogen**2)

    return
  end function Cooling_Function_CIE_File_Interpolate

  double precision function Cooling_Function_CIE_File_logTemperature_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)
    !% Compute the logarithmic gradient of the cooling function with respect to temperature by interpolation in the tabulated data.
    use Abundances_Structure
    use Radiation_Structure
    use Numerical_Interpolation
    use Numerical_Constants_Astronomical
    use IO_XML
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation
    double precision,          save       :: temperaturePrevious=-1.0d0,metallicityPrevious=-1.0d0,coolingFunctionSlopePrevious
    !$omp threadprivate(temperaturePrevious,metallicityPrevious,coolingFunctionSlopePrevious)
    integer                               :: iTemperature,iMetallicity
    double precision                      :: temperatureUse,metallicityUse,hTemperature,hMetallicity

    ! Handle out of range temperatures.
    temperatureUse=temperature
    if (temperatureUse < temperatureMinimum) then
       select case (extrapolateTemperatureLow)
       case (extrapolateZero,extrapolateFixed)
          Cooling_Function_CIE_File_logTemperature_Interpolate=0.0d0
          return
       case (extrapolatePowerLaw)
          temperatureUse=temperatureMinimum
       end select
    end if
    if (temperatureUse > temperatureMaximum) then
       select case (extrapolateTemperatureHigh)
       case (extrapolateZero,extrapolateFixed)
          Cooling_Function_CIE_File_logTemperature_Interpolate=0.0d0
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
          Cooling_Function_CIE_File_logTemperature_Interpolate=0.0d0
          return
       case (extrapolateFixed)
          metallicityUse=metallicityMinimum
       end select
    end if
    if (metallicityUse > metallicityMaximum) then
       select case (extrapolateMetallicityHigh)
       case (extrapolateZero)
          Cooling_Function_CIE_File_logTemperature_Interpolate=0.0d0
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
       coolingFunctionSlopePrevious=( ( coolingFunctionTable(iTemperature+1,iMetallicity  )                       &
            &                          -coolingFunctionTable(iTemperature  ,iMetallicity  ))*(1.0d0-hMetallicity) &
            &                        +( coolingFunctionTable(iTemperature+1,iMetallicity+1)                       & 
            &                          -coolingFunctionTable(iTemperature  ,iMetallicity+1))*       hMetallicity) &
            & /(coolingFunctionTemperatures(iTemperature+1)-coolingFunctionTemperatures(iTemperature))
       
       ! Convert to logarithmic gradient if table was not stored logarithmically.
       if (.not.logarithmicTable) coolingFunctionSlopePrevious=coolingFunctionSlopePrevious*temperature&
            &/Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity)
       
       ! Store the temperature and metallicity for which calculation was performed.
       temperaturePrevious=temperatureUse
       metallicityPrevious=metallicityUse
       
    end if
    
    ! Return the stored value.
    Cooling_Function_CIE_File_logTemperature_Interpolate=coolingFunctionSlopePrevious

    return
  end function Cooling_Function_CIE_File_logTemperature_Interpolate

  subroutine Cooling_Function_CIE_File_Read(coolingFunctionFileToRead,metallicityMaximumTabulated)
    !% Read in data from a cooling function file.
    use Galacticus_Error
    use FoX_dom
    use Memory_Management
    use Numerical_Comparison
    use Galacticus_Display
    use IO_XML
    implicit none
    type(varying_string), intent(in)            :: coolingFunctionFileToRead
    double precision,     intent(out), optional :: metallicityMaximumTabulated
    type(Node),           pointer               :: doc,datum,thisCoolingFunction,metallicityElement,extrapolationElement&
         &,extrapolation,thisTemperature,thisCoolingRate
    type(NodeList),       pointer               :: temperatureDatumList,coolingDatumList,coolingFunctionList&
         &,metallicityExtrapolationList ,temperatureExtrapolationList
    integer                                     :: iDatum,ioErr,iCoolingFunction,iExtrapolation,extrapolationMethod
    double precision                            :: datumValues(1)
    character(len=32)                           :: limitType

    !$omp critical (FoX_DOM_Access)
    ! Parse the XML file.
    call Galacticus_Display_Indent('Parsing file: '//coolingFunctionFileToRead,verbosityWorking)
    call Galacticus_Display_Counter(0,.true.,verbosityWorking)
    doc => parseFile(char(coolingFunctionFileToRead),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Cooling_Function_CIE_File_Read','Unable to find cooling function file')
    call Galacticus_Display_Counter(50,.false.,verbosityWorking)
    ! Get a list of all <coolingFunction> elements.
    coolingFunctionList => getElementsByTagname(doc,"coolingFunction")
    coolingFunctionMetallicityNumberPoints=getLength(coolingFunctionList)
    ! Extract data from first cooling function and count number of temperatures present.
    thisCoolingFunction  => item(coolingFunctionList,0)
    thisTemperature      => item(getElementsByTagname(thisCoolingFunction,"temperature"),0)
    temperatureDatumList => getElementsByTagname(thisTemperature,"datum")
    coolingFunctionTemperatureNumberPoints=getLength(temperatureDatumList)
    ! Allocate space for the table.
    if (allocated(coolingFunctionMetallicities)) call Dealloc_Array(coolingFunctionMetallicities)
    if (allocated(coolingFunctionTemperatures )) call Dealloc_Array(coolingFunctionTemperatures )
    if (allocated(coolingFunctionTable        )) call Dealloc_Array(coolingFunctionTable        )
    call Alloc_Array(coolingFunctionMetallicities                                       ,[coolingFunctionMetallicityNumberPoints])
    call Alloc_Array(coolingFunctionTemperatures ,[coolingFunctionTemperatureNumberPoints                                       ])
    call Alloc_Array(coolingFunctionTable        ,[coolingFunctionTemperatureNumberPoints,coolingFunctionMetallicityNumberPoints])

    ! Extract data from the cooling functions and populate metallicity and temperature arrays.
    do iCoolingFunction=0,getLength(coolingFunctionList)-1
       ! Get required cooling function.
       thisCoolingFunction => item(coolingFunctionList,iCoolingFunction)
       ! Extract the metallicity from the <metallicity> element.
       metallicityElement  => item(getElementsByTagname(thisCoolingFunction,"metallicity"),0)
       call extractDataContent(metallicityElement,coolingFunctionMetallicities(iCoolingFunction+1))
       ! Extract the data.
       thisTemperature      => item(getElementsByTagname(thisCoolingFunction,"temperature"),0)
       temperatureDatumList =>      getElementsByTagname(thisTemperature    ,"datum"      )
       thisCoolingRate      => item(getElementsByTagname(thisCoolingFunction,"coolingRate"),0)
       coolingDatumList     =>      getElementsByTagname(thisCoolingRate    ,"datum"      )
       ! Check that number of temperatures is consistent.
       if (getLength(temperatureDatumList) /= coolingFunctionTemperatureNumberPoints) call&
            & Galacticus_Error_Report('Cooling_Function_CIE_File_Read','sizes of temperatures grids must be the same for all&
            & metallicities')
       ! Check that number of cooling rates matches number of temperatures.
       if (getLength(temperatureDatumList) /= getLength(coolingDatumList)) call&
            & Galacticus_Error_Report('Cooling_Function_CIE_File_Read','sizes of temperature and cooling rate arrays must match')
       ! Extract data.
       do iDatum=0,getLength(temperatureDatumList)-1
          ! Report progress.
          call Galacticus_Display_Counter(                                                                           &
               &                           int(                                                                      &
               &                                50.0d0                                                               &
               &                               +50.0d0                                                               &
               &                               *dble(                                                                &
               &                                      getLength(temperatureDatumList)*iDatum                         &
               &                                     +                                iCoolingFunction               &
               &                                    )                                                                &
               &                               /dble(getLength(temperatureDatumList)*getLength(coolingFunctionList)) &
               &                              )                                                                      &
               &                          ,.false.                                                                   &
               &                          ,verbosityWorking                                                          &
               &                         )
          datum => item(temperatureDatumList,iDatum)
          call extractDataContent(datum,datumValues)
          if (iCoolingFunction == 0) then
             ! Store the temperature.
             coolingFunctionTemperatures(iDatum+1)=datumValues(1)
          else
             ! Check that temperature grids are aligned.
             if (Values_Differ(coolingFunctionTemperatures(iDatum+1),datumValues(1),relTol=1.0d-6)) call&
                  & Galacticus_Error_Report('Cooling_Function_CIE_File_Read','temperature grids mismatch')
          end if
          ! Store the cooling function.
          datum => item(coolingDatumList,iDatum)
          call extractDataContent(datum,datumValues)
          coolingFunctionTable(iDatum+1,iCoolingFunction+1)=datumValues(1)
       end do
    end do
    where (coolingFunctionMetallicities>-999.0d0)
       coolingFunctionMetallicities=10.0d0**coolingFunctionMetallicities
    elsewhere
       coolingFunctionMetallicities=0.0d0
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
         call Galacticus_Error_Report('Cooling_Function_CIE_File_Read','unrecognized extrapolation limit')
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
          call Galacticus_Error_Report('Cooling_Function_CIE_File_Read','unrecognized extrapolation limit')
       end select
    end do
    ! Destroy the document.
    call destroy(doc)
    call Galacticus_Display_Counter_Clear(verbosityWorking)
    call Galacticus_Display_Unindent('done',verbosityWorking)
    !$omp end critical (FoX_DOM_Access)
  
    ! Store table ranges for convenience.
    metallicityMinimum=coolingFunctionMetallicities(1)
    metallicityMaximum=coolingFunctionMetallicities(coolingFunctionMetallicityNumberPoints)
    temperatureMinimum=coolingFunctionTemperatures(1)
    temperatureMaximum=coolingFunctionTemperatures(coolingFunctionTemperatureNumberPoints)

    ! Decide whether or not to make the tables logarithmic.
    logarithmicTable=all(coolingFunctionTable > 0.0d0)
    if (logarithmicTable) then
       firstMetallicityIsZero=(coolingFunctionMetallicities(1) == 0.0d0)
       if (firstMetallicityIsZero) firstNonZeroMetallicity=coolingFunctionMetallicities(2)
       where (coolingFunctionMetallicities > 0.0d0)
          coolingFunctionMetallicities=dlog(coolingFunctionMetallicities)
       elsewhere
          coolingFunctionMetallicities=-999.0d0
       end where
       coolingFunctionTemperatures=dlog(coolingFunctionTemperatures)
       coolingFunctionTable=dlog(coolingFunctionTable)
    else
       if (        extrapolateTemperatureLow  == extrapolatePowerLaw   &
            & .or. extrapolateTemperatureHigh == extrapolatePowerLaw   &
            & .or. extrapolateMetallicityLow  == extrapolatePowerLaw   &
            & .or. extrapolateMetallicityHigh == extrapolatePowerLaw ) &
            & call Galacticus_Error_Report('Cooling_Function_CIE_File_Read','power law extrapolation allowed only in loggable tables')
    end if

    ! Force interpolation accelerators to be reset.
    resetTemperature=.true.
    resetMetallicity=.true.

    ! Return the maximum metallicity if requested.
    if (present(metallicityMaximumTabulated)) metallicityMaximumTabulated=metallicityMaximum

    return
  end subroutine Cooling_Function_CIE_File_Read

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
    iTemperature=Interpolate_Locate(coolingFunctionTemperatureNumberPoints,coolingFunctionTemperatures&
         &,interpolationAcceleratorTemperature,temperatureUse,resetTemperature)
    iTemperature=max(min(iTemperature,coolingFunctionTemperatureNumberPoints),1)
    hTemperature=(temperatureUse-coolingFunctionTemperatures(iTemperature))/(coolingFunctionTemperatures(iTemperature+1)&
         &-coolingFunctionTemperatures(iTemperature))

    ! Get interpolation in metallicity.
    if (firstMetallicityIsZero.and.metallicityUse < firstNonZeroMetallicity) then
       iMetallicity=1
       hMetallicity=metallicityUse/firstNonZeroMetallicity
    else
       if (logarithmicTable) metallicityUse=dlog(metallicityUse)
       iMetallicity=Interpolate_Locate(coolingFunctionMetallicityNumberPoints,coolingFunctionMetallicities&
         &,interpolationAcceleratorMetallicity,metallicityUse,resetMetallicity)
       iMetallicity=max(min(iMetallicity,coolingFunctionMetallicityNumberPoints),1)
       hMetallicity=(metallicityUse-coolingFunctionMetallicities(iMetallicity))/(coolingFunctionMetallicities(iMetallicity+1)&
            &-coolingFunctionMetallicities(iMetallicity))
    end if
    return
  end subroutine Get_Interpolation

  double precision function Do_Interpolation(iTemperature,hTemperature,iMetallicity,hMetallicity)
    !% Perform the interpolation.
    implicit none
    integer,          intent(in) :: iTemperature,iMetallicity
    double precision, intent(in) :: hTemperature,hMetallicity

    ! Do the interpolation.
    Do_Interpolation=coolingFunctionTable(iTemperature  ,iMetallicity  )*(1.0d0-hTemperature)*(1.0d0-hMetallicity)&
         &          +coolingFunctionTable(iTemperature  ,iMetallicity+1)*(1.0d0-hTemperature)*       hMetallicity &
         &          +coolingFunctionTable(iTemperature+1,iMetallicity  )*       hTemperature *(1.0d0-hMetallicity)&
         &          +coolingFunctionTable(iTemperature+1,iMetallicity+1)*       hTemperature *       hMetallicity

    ! Exponentiate the result if the table was stored as the log.
    if (logarithmicTable) Do_Interpolation=dexp(Do_Interpolation)

    return
  end function Do_Interpolation

end module Cooling_Functions_CIE_File

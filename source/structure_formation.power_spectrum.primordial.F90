!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

module Primordial_Power_Spectra
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Primordial_Power_Spectrum, Primordial_Power_Spectrum_Logarithmic_Derivative,&
       & Primordial_Power_Spectrum_State_Retrieve

  ! Flag to indicate if this module has been initialized.
  logical                                                                       :: powerSpectrumInitialized =.false., tablesInitialized         =.false.

  ! Variables to hold the tabulated power spectrum data.
  integer                                                                       :: powerSpectrumNumberPoints=-1
  double precision                                  , allocatable, dimension(:) :: powerSpectrumLogP                , powerSpectrumLogWavenumber
  type            (fgsl_interp                     )                            :: interpolationObject
  type            (fgsl_interp_accel               )                            :: interpolationAccelerator
  logical                                                                       :: resetInterpolation       =.true.

  ! Name of power spectrum method used.
  type            (varying_string                  )                            :: powerSpectrumMethod

  ! Pointer to the subroutine that tabulates the transfer function and template interface for that subroutine.
  procedure       (Power_Spectrum_Tabulate_Template), pointer                   :: Power_Spectrum_Tabulate  =>null()
  abstract interface
     subroutine Power_Spectrum_Tabulate_Template(logWavenumber,powerSpectrumNumberPoints,powerSpectrumLogWavenumber &
          &,powerSpectrumLogP)
    double precision                           , intent(in   ) :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: powerSpectrumLogP        , powerSpectrumLogWavenumber
    integer                                    , intent(  out) :: powerSpectrumNumberPoints
  end subroutine Power_Spectrum_Tabulate_Template
 end interface

contains

  double precision function Primordial_Power_Spectrum(wavenumber)
    !% Return the CDM primordial power spectrum for $k=${\tt wavenumber} [Mpc$^{-1}$].
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   ) :: wavenumber
    double precision                :: logWavenumber

    ! Get logarithm of wavenumber.
    logWavenumber=log(wavenumber)

    !$omp critical(Power_Spectrum_Initialization)
    ! Initialize if necessary.
    if (.not.(powerSpectrumInitialized.and.tablesInitialized)) then
       call Power_Spectrum_Initialize(logWavenumber)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
    end if

    ! If wavenumber is out of range, attempt to remake the table.
    if (logWavenumber<powerSpectrumLogWavenumber(1) .or. logWavenumber>powerSpectrumLogWavenumber(powerSpectrumNumberPoints) )&
         & then
       call Power_Spectrum_Tabulate(logWavenumber,powerSpectrumNumberPoints,powerSpectrumLogWavenumber,powerSpectrumLogP)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
   end if
    !$omp end critical(Power_Spectrum_Initialization)

    ! Interpolate in the tabulated function and return a value.
    Primordial_Power_Spectrum=exp(Interpolate(powerSpectrumNumberPoints,powerSpectrumLogWavenumber,powerSpectrumLogP &
         &,interpolationObject,interpolationAccelerator,logWavenumber,reset=resetInterpolation))

    return
  end function Primordial_Power_Spectrum

  double precision function Primordial_Power_Spectrum_Logarithmic_Derivative(wavenumber)
    !% Return the logarithmic derivative CDM primordial power spectrum for $k=${\tt wavenumber} [Mpc$^{-1}$].
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   ) :: wavenumber
    double precision                :: logWavenumber

    ! Get logarithm of wavenumber.
    logWavenumber=log(wavenumber)

    !$omp critical(Power_Spectrum_Initialization)
    ! Initialize if necessary.
    if (.not.(powerSpectrumInitialized.and.tablesInitialized)) then
       call Power_Spectrum_Initialize(logWavenumber)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
    end if

    ! If wavenumber is out of range, attempt to remake the table.
    if (logWavenumber<powerSpectrumLogWavenumber(1) .or. logWavenumber>powerSpectrumLogWavenumber(powerSpectrumNumberPoints) )&
         & then
       call Power_Spectrum_Tabulate(logWavenumber,powerSpectrumNumberPoints,powerSpectrumLogWavenumber,powerSpectrumLogP)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
    end if
    !$omp end critical(Power_Spectrum_Initialization)

    ! Interpolate in the tabulated function and return a value.
    Primordial_Power_Spectrum_Logarithmic_Derivative=Interpolate_Derivative(powerSpectrumNumberPoints,powerSpectrumLogWavenumber&
         &,powerSpectrumLogP ,interpolationObject,interpolationAccelerator,logWavenumber,reset=resetInterpolation)

    return
  end function Primordial_Power_Spectrum_Logarithmic_Derivative

  subroutine Power_Spectrum_Initialize(logWavenumber)
    !% Initializes the transfer function module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="powerSpectrumMethod" type="moduleUse">
    include 'structure_formation.power_spectrum.primordial.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in   ) :: logWavenumber

    if (.not.powerSpectrumInitialized) then
       ! Get the primordial power spectrum method parameter.
       !@ <inputParameter>
       !@   <name>powerSpectrumMethod</name>
       !@   <defaultValue>powerLaw</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing the primordial power spectrum.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('powerSpectrumMethod',powerSpectrumMethod,defaultValue='powerLaw')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="powerSpectrumMethod" type="functionCall" functionType="void">
       !#  <functionArgs>powerSpectrumMethod,Power_Spectrum_Tabulate</functionArgs>
       include 'structure_formation.power_spectrum.inc'
       !# </include>
       if (.not.associated(Power_Spectrum_Tabulate)) call Galacticus_Error_Report('Power_Spectrum_Initialize','method '&
            &//char(powerSpectrumMethod)//' is unrecognized')
    end if
    ! Call routine to initialize the transfer function.
    call Power_Spectrum_Tabulate(logWavenumber,powerSpectrumNumberPoints,powerSpectrumLogWavenumber,powerSpectrumLogP)
    tablesInitialized=.true.
    ! Flag that the module is now initialized.
    powerSpectrumInitialized=.true.
    return
  end subroutine Power_Spectrum_Initialize

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Primordial_Power_Spectrum_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Primordial_Power_Spectrum_State_Retrieve(stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use Memory_Management
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    powerSpectrumNumberPoints=0
    if (allocated(powerSpectrumLogWavenumber)) call Dealloc_Array(powerSpectrumLogWavenumber)
    if (allocated(powerSpectrumLogP         )) call Dealloc_Array(powerSpectrumLogP         )
    tablesInitialized=.false.
    return
  end subroutine Primordial_Power_Spectrum_State_Retrieve

end module Primordial_Power_Spectra

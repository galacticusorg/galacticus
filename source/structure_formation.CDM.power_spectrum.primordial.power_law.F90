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

!% Contains a module which generates a tabulated power-law primordial power spectrum.

module CDM_Primordial_Power_Spectrum_Power_Law
  !% Implements generation of a tabulated power-law primordial power spectrum. The default power spectrum parameters are taken
  !% from \cite{komatsu_seven-year_2010}.
  implicit none
  private
  public :: CDM_Primordial_Power_Spectrum_Power_Law_Initialize, CDM_Primordial_Power_Spectrum_Power_Law_State_Store,&
       & CDM_Primordial_Power_Spectrum_Power_Law_State_Retrieve
  
  ! Parameters of the power-law.
  double precision            :: powerSpectrumIndex, powerSpectrumRunning, powerSpectrumReferenceWavenumber

  ! Parameters controlling the gridding of the power spectrum and default wavenumber range.
  integer,          parameter :: nPointsPerDecade=1000
  double precision            :: logWavenumberMinimum=dlog(1.0d-5), logWavenumberMaximum=dlog(10.0d0)

contains
  
  !# <powerSpectrumMethod>
  !#  <unitName>CDM_Primordial_Power_Spectrum_Power_Law_Initialize</unitName>
  !# </powerSpectrumMethod>
  subroutine CDM_Primordial_Power_Spectrum_Power_Law_Initialize(powerSpectrumMethod,Power_Spectrum_Tabulate)
    !% Initializes the ``transfer function from CMBFast'' module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: powerSpectrumMethod
    procedure(),          pointer, intent(inout) :: Power_Spectrum_Tabulate
    
    if (powerSpectrumMethod == 'powerLaw') then
       Power_Spectrum_Tabulate => Power_Spectrum_Power_Law_Tabulate
       !@ <inputParameter>
       !@   <name>powerSpectrumIndex</name>
       !@   <defaultValue>0.9538 (\citealt{story_measurement_2012}; CMB$+H_0+$BAO)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The index of the power-law primordial power spectrum.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('powerSpectrumIndex'              ,powerSpectrumIndex              ,defaultValue=0.9538d0)
       !@ <inputParameter>
       !@   <name>powerSpectrumRunning</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The running, $\d n_{\rm s} / \d \ln k$, of the power spectrum index.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('powerSpectrumRunning'            ,powerSpectrumRunning            ,defaultValue=0.0d0  )
       !@ <inputParameter>
       !@   <name>powerSpectrumReferenceWavenumber</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     When a running power spectrum index is used, this is the wavenumber at which the index is equal to {\tt [powerSpectrumIndex]}.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('powerSpectrumReferenceWavenumber',powerSpectrumReferenceWavenumber,defaultValue=1.0d0  )
    end if
    return
  end subroutine CDM_Primordial_Power_Spectrum_Power_Law_Initialize

  subroutine Power_Spectrum_Power_Law_Tabulate(logWavenumber,powerSpectrumNumberPoints,powerSpectrumLogWavenumber&
       &,powerSpectrumLogP)
    !% Tabulate a power-law primordial power spectrum.
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Constants_Math
    implicit none
    double precision,                            intent(in)    :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: powerSpectrumLogWavenumber,powerSpectrumLogP
    integer,                                     intent(out)   :: powerSpectrumNumberPoints
    integer                                                    :: iWavenumber
    double precision                                           :: wavenumber

    ! Determine range of wavenumbers required.
    logWavenumberMinimum=min(logWavenumberMinimum,logWavenumber-ln10)
    logWavenumberMaximum=max(logWavenumberMaximum,logWavenumber+ln10)
    
    ! Determine number of points to tabulate.
    powerSpectrumNumberPoints=int((logWavenumberMaximum-logWavenumberMinimum)*dble(nPointsPerDecade)/ln10)

    ! Deallocate arrays if currently allocated.
    if (allocated(powerSpectrumLogWavenumber)) call Dealloc_Array(powerSpectrumLogWavenumber)
    if (allocated(powerSpectrumLogP))          call Dealloc_Array(powerSpectrumLogP         )
    ! Allocate the arrays to current required size.
    call Alloc_Array(powerSpectrumLogWavenumber,[powerSpectrumNumberPoints])
    call Alloc_Array(powerSpectrumLogP         ,[powerSpectrumNumberPoints])

    ! Tabulate the function.
    powerSpectrumLogWavenumber=Make_Range(logWavenumberMinimum,logWavenumberMaximum,powerSpectrumNumberPoints,rangeTypeLinear)
    do iWavenumber=1,powerSpectrumNumberPoints
       wavenumber=dexp(powerSpectrumLogWavenumber(iWavenumber))
       powerSpectrumLogP(iWavenumber)=(powerSpectrumIndex+0.5d0*powerSpectrumRunning*dlog(wavenumber&
            &/powerSpectrumReferenceWavenumber))*powerSpectrumLogWavenumber(iWavenumber)
    end do
    
    return
  end subroutine Power_Spectrum_Power_Law_Tabulate

  !# <galacticusStateStoreTask>
  !#  <unitName>CDM_Primordial_Power_Spectrum_Power_Law_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine CDM_Primordial_Power_Spectrum_Power_Law_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) logWavenumberMinimum,logWavenumberMaximum
    return
  end subroutine CDM_Primordial_Power_Spectrum_Power_Law_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>CDM_Primordial_Power_Spectrum_Power_Law_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine CDM_Primordial_Power_Spectrum_Power_Law_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) logWavenumberMinimum,logWavenumberMaximum
    return
  end subroutine CDM_Primordial_Power_Spectrum_Power_Law_State_Retrieve
    
end module CDM_Primordial_Power_Spectrum_Power_Law

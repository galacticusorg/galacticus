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


!% Contains a module which generates a tabulated power-law primordial power spectrum.

module CDM_Primordial_Power_Spectrum_Power_Law
  !% Implements generation of a tabulated power-law primordial power spectrum. The default power spectrum parameters are taken
  !% from \cite{komatsu_seven-year_2010}.
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
    
    if (powerSpectrumMethod.eq.'power law') then
       Power_Spectrum_Tabulate => Power_Spectrum_Power_Law_Tabulate
       !@ <inputParameter>
       !@   <name>powerSpectrumIndex</name>
       !@   <defaultValue>0.961 \citep{komatsu_seven-year_2010}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The index of the power-law primordial power spectrum.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('powerSpectrumIndex'              ,powerSpectrumIndex              ,defaultValue=0.961d0)
       !@ <inputParameter>
       !@   <name>powerSpectrumRunning</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The running, $\d n_{\rm s} / \d \ln k$, of the power spectrum index.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('powerSpectrumRunning'            ,powerSpectrumRunning            ,defaultValue=0.0d0  )
       !@ <inputParameter>
       !@   <name>powerSpectrumReferenceWavenumber</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     When a running power spectrum index is used, this is the wavenumber at which the index is equal to {\tt [powerSpectrumIndex]}.
       !@   </description>
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
       powerSpectrumLogP(iWavenumber)=(powerSpectrumIndex+powerSpectrumRunning*dlog(wavenumber&
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

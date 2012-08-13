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


!% Contains a module which implements window functions for computing the variance of the power spectrum.

module Power_Spectrum_Window_Functions
  !% Implements window functions for computing the variance of the power spectrum.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  implicit none
  private
  public :: Power_Spectrum_Window_Function
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: moduleInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: powerSpectrumWindowFunctionMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Power_Spectrum_Window_Function), pointer :: Power_Spectrum_Window_Function_Get => null()

contains

  subroutine Power_Spectrum_Window_Functions_Initialize
    !% Initialize the power spectrum window function module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="powerSpectrumWindowFunctionMethod" type="moduleUse">
    include 'structure_formation.CDM.power_spectrum.variance.window_function.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Power_Spectrum_Window_Functions_Initialization) 
       if (.not.moduleInitialized) then
          ! Get the window function method parameter.
          !@ <inputParameter>
          !@   <name>powerSpectrumWindowFunctionMethod</name>
          !@   <defaultValue>topHat</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing window functions to estimate variance from the power spectrum.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('powerSpectrumWindowFunctionMethod',powerSpectrumWindowFunctionMethod,defaultValue='topHat')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="powerSpectrumWindowFunctionMethod" type="code" action="subroutine">
          !#  <subroutineArgs>powerSpectrumWindowFunctionMethod,Power_Spectrum_Window_Function_Get</subroutineArgs>
          include 'structure_formation.CDM.power_spectrum.variance.window_function.inc'
          !# </include>
          if (.not.associated(Power_Spectrum_Window_Function_Get)) call Galacticus_Error_Report('Power_Spectrum_Window_Functions_Initialize'&
               &,'method '//char(powerSpectrumWindowFunctionMethod)//' is unrecognized')
          moduleInitialized=.true.
       end if
       !$omp end critical(Power_Spectrum_Window_Functions_Initialization) 
    end if
    return
  end subroutine Power_Spectrum_Window_Functions_Initialize

  double precision function Power_Spectrum_Window_Function(wavenumber,smoothingMass)
    !% Returns the window function for power spectrum variance computation at the specified {\tt wavenumber} (in Mpc$^{-1}$) for a
    !% given {\tt smoothingMass} (in $M_\odot$).
    implicit none
    double precision, intent(in) :: wavenumber,smoothingMass

    ! Initialize the module.
    call Power_Spectrum_Window_Functions_Initialize

    ! Call the function that does the work.
    Power_Spectrum_Window_Function=Power_Spectrum_Window_Function_Get(wavenumber,smoothingMass)
    return
  end function Power_Spectrum_Window_Function
  
end module Power_Spectrum_Window_Functions

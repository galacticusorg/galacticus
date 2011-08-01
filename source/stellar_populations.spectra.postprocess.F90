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


!% Contains a module that implements postprocessing of stellar population spectra.

module Stellar_Population_Spectra_Postprocess
  !% Implements postprocessing of stellar population spectra.
  use Abundances_Structure
  use ISO_Varying_String 
  !# <include directive="stellarPopulationSpectraPostprocessMethod" type="moduleUse">
  include 'stellar_populations.spectra.postprocess.modules.inc'
  !# </include>
  implicit none
  private
  public :: Stellar_Population_Spectrum_Postprocess

  ! Flag to indicate if this module has been initialized.  
  logical              :: stellarPopulationSpectraPostprocessInitialized=.false.

  ! Name of cooling time available method used.
  type(varying_string) :: stellarPopulationSpectraPostprocessMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Stellar_Population_Spectrum_Postprocess_Get_Template), pointer :: Stellar_Population_Spectrum_Postprocess_Get => null()
  abstract interface
     double precision function Stellar_Population_Spectrum_Postprocess_Get_Template(wavelength,redshift)
       double precision, intent(in) :: wavelength,redshift
     end function Stellar_Population_Spectrum_Postprocess_Get_Template
  end interface
 
contains

  subroutine Stellar_Population_Spectrum_Postprocess_Initialize
    !% Initialize the stellar population spectra postprocessing module
    use Galacticus_Error
    use Input_Parameters
    implicit none
    
    !$omp critical(Stellar_Population_Spectrum_Postprocess_Initialization) 
    ! Initialize if necessary.
    if (.not.stellarPopulationSpectraPostprocessInitialized) then
       ! Get the cooling function method parameter.
       !@ <inputParameter>
       !@   <name>stellarPopulationSpectraPostprocessMethod</name>
       !@   <defaultValue>Meiksin2006</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for post-processing of stellar population spectra.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationSpectraPostprocessMethod',stellarPopulationSpectraPostprocessMethod,defaultValue='Meiksin2006')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="stellarPopulationSpectraPostprocessMethod" type="code" action="subroutine">
       !#  <subroutineArgs>stellarPopulationSpectraPostprocessMethod,Stellar_Population_Spectrum_Postprocess_Get</subroutineArgs>
       include 'stellar_populations.spectra.postprocess.inc'
       !# </include>
       if (.not.associated(Stellar_Population_Spectrum_Postprocess_Get)) call&
            & Galacticus_Error_Report('Stellar_Population_Spectrum_Postprocess','method ' //char(stellarPopulationSpectraPostprocessMethod)//' is&
            & unrecognized')
       stellarPopulationSpectraPostprocessInitialized=.true.
    end if
    !$omp end critical(Stellar_Population_Spectrum_Postprocess_Initialization) 
    return
  end subroutine Stellar_Population_Spectrum_Postprocess_Initialize

  double precision function Stellar_Population_Spectrum_Postprocess(wavelength,redshift)
    !% Return a multiplicative factor by which a stellar population spectrum should be modified by any postprocessing.
    implicit none
    double precision, intent(in) :: wavelength,redshift
    
    ! Initialize the module.
    call Stellar_Population_Spectrum_Postprocess_Initialize
    
    ! Get the spectrum using the selected method.
    Stellar_Population_Spectrum_Postprocess=Stellar_Population_Spectrum_Postprocess_Get(wavelength,redshift)    
    return
  end function Stellar_Population_Spectrum_Postprocess
  
end module Stellar_Population_Spectra_Postprocess

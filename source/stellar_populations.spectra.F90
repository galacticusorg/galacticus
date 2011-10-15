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


!% Contains a module that implements calculations of stellar population spectra.

module Stellar_Population_Spectra
  !% Implements calculations of stellar population spectra.
  use Abundances_Structure
  use ISO_Varying_String 
  !# <include directive="stellarPopulationSpectraMethod" type="moduleUse">
  include 'stellar_populations.spectra.modules.inc'
  !# </include>
  implicit none
  private
  public :: Stellar_Population_Spectrum, Stellar_Population_Spectrum_Tabulation

  ! Flag to indicate if this module has been initialized.  
  logical              :: stellarPopulationSpectraInitialized=.false.

  ! Name of cooling time available method used.
  type(varying_string) :: stellarPopulationSpectraMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Stellar_Population_Spectrum_Get_Template), pointer :: Stellar_Population_Spectrum_Get => null()
  abstract interface
     double precision function Stellar_Population_Spectrum_Get_Template(abundances,age,wavelength,imfIndex)
       import abundancesStructure
       type(abundancesStructure), intent(in) :: abundances
       double precision,          intent(in) :: age,wavelength
       integer,                   intent(in) :: imfIndex
     end function Stellar_Population_Spectrum_Get_Template
  end interface
  procedure(Stellar_Population_Tabulation_Get_Template), pointer :: Stellar_Population_Spectrum_Tabulation_Get => null()
  abstract interface
     subroutine Stellar_Population_Tabulation_Get_Template(imfIndex,agesCount,metallicitiesCount,ages,metallicity)
       integer,          intent(in)                             :: imfIndex
       integer,          intent(out)                            :: agesCount,metallicitiesCount
       double precision, intent(out), allocatable, dimension(:) :: ages,metallicity
     end subroutine Stellar_Population_Tabulation_Get_Template
  end interface
 
contains

  subroutine Stellar_Population_Spectrum_Initialize
    !% Initialize the stellar population spectra module
    use Galacticus_Error
    use Input_Parameters
    implicit none
    
    !$omp critical(Stellar_Population_Spectrum_Initialization) 
    ! Initialize if necessary.
    if (.not.stellarPopulationSpectraInitialized) then
       ! Get the cooling function method parameter.
       !@ <inputParameter>
       !@   <name>stellarPopulationSpectraMethod</name>
       !@   <defaultValue>Conroy-White-Gunn2010</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of stellar population spectra.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPopulationSpectraMethod',stellarPopulationSpectraMethod,defaultValue='Conroy-White-Gunn2009')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="stellarPopulationSpectraMethod" type="code" action="subroutine">
       !#  <subroutineArgs>stellarPopulationSpectraMethod,Stellar_Population_Spectrum_Get,Stellar_Population_Spectrum_Tabulation_Get</subroutineArgs>
       include 'stellar_populations.spectra.inc'
       !# </include>
       if (.not.(associated(Stellar_Population_Spectrum_Get).and.associated(Stellar_Population_Spectrum_Tabulation_Get))) call&
            & Galacticus_Error_Report('Stellar_Population_Spectrum','method ' //char(stellarPopulationSpectraMethod)//' is&
            & unrecognized')
       stellarPopulationSpectraInitialized=.true.
    end if
    !$omp end critical(Stellar_Population_Spectrum_Initialization) 
    return
  end subroutine Stellar_Population_Spectrum_Initialize

  double precision function Stellar_Population_Spectrum(abundances,age,wavelength,imfIndex)
    !% Return the luminosity (in units of $L_\odot$ Hz$^{-1}$) for a stellar population with composition {\tt abundances}, of the
    !% given {\tt age} (in Gyr) and the specified {\tt wavelength} (in Angstroms).
    implicit none
    type(abundancesStructure), intent(in) :: abundances
    double precision,          intent(in) :: age,wavelength
    integer,                   intent(in) :: imfIndex
    
    ! Initialize the module.
    call Stellar_Population_Spectrum_Initialize
    
    ! Get the spectrum using the selected method.
    Stellar_Population_Spectrum=Stellar_Population_Spectrum_Get(abundances,age,wavelength,imfIndex)
    
    return
  end function Stellar_Population_Spectrum
  
  subroutine Stellar_Population_Spectrum_Tabulation(imfIndex,agesCount,metallicitiesCount,age,metallicity)
    !% Return a tabulation of ages and metallicities at which stellar spectra for the specified IMF should be tabulated.
    implicit none
    integer,          intent(in)                             :: imfIndex
    integer,          intent(out)                            :: agesCount,metallicitiesCount
    double precision, intent(out), allocatable, dimension(:) :: age,metallicity
    
    ! Initialize the module.
    call Stellar_Population_Spectrum_Initialize
    
    ! Get the tabulation using the selected method.
    call Stellar_Population_Spectrum_Tabulation_Get(imfIndex,agesCount,metallicitiesCount,age,metallicity)
    return
  end subroutine Stellar_Population_Spectrum_Tabulation

end module Stellar_Population_Spectra

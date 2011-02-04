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


!% Contains a module that implements calculations of the intergalactic medium thermal and ionization state using RecFast.

module Intergalactic_Medium_State_RecFast
  !% Implements calculations of the intergalactic medium thermal and ionization state using RecFast.
  use ISO_Varying_String
  private
  public :: Intergalactic_Medium_State_RecFast_Initialize
  
contains

  !# <intergalaticMediumStateMethod>
  !#  <unitName>Intergalactic_Medium_State_RecFast_Initialize</unitName>
  !# </intergalaticMediumStateMethod>
  subroutine Intergalactic_Medium_State_RecFast_Initialize(intergalaticMediumStateMethod,Intergalactic_Medium_Electron_Fraction_Get&
       &,Intergalactic_Medium_Temperature_Get)
    !% Initializes the ``RecFast'' intergalactic medium state module.
    use FoX_wxml
    use System_Command
    use Intergalactic_Medium_State_File
    use Cosmology_Functions
    use Cosmological_Parameters
    use Numerical_Constants_Astronomical
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: intergalaticMediumStateMethod
    procedure(),          pointer, intent(inout) :: Intergalactic_Medium_Electron_Fraction_Get,Intergalactic_Medium_Temperature_Get
    character(len=32)                            :: parameterLabel
    type(varying_string)                         :: parameterFile,recfastFile,command
    type(xmlf_t)                                 :: parameterDoc

    ! Test if our method has been selected.    
    if (intergalaticMediumStateMethod == 'RecFast') then
       ! Set procedure pointers (use those from the ``file'' implementation).
       Intergalactic_Medium_Electron_Fraction_Get => Intergalactic_Medium_Electron_Fraction_File
       Intergalactic_Medium_Temperature_Get       => Intergalactic_Medium_Temperature_File

       ! Ensure that the RecFast data file has been generated.
       
       ! Generate the name of the data file and an XML input parameter file.
       recfastFile='data/recFast'
       parameterFile='data/recfast_parameters.xml'
       call xml_OpenFile(char(parameterFile),parameterDoc)
       call xml_NewElement(parameterDoc,"parameters")
       write (parameterLabel,'(f5.3)') Omega_0()
       recfastFile=recfastFile//'_Omega0'//trim(parameterLabel)
       call Write_Parameter(parameterDoc,"Omega_0",parameterLabel)
       write (parameterLabel,'(f5.3)') Omega_DE()
       recfastFile=recfastFile//'_OmegaDE'//trim(parameterLabel)
       call Write_Parameter(parameterDoc,"Omega_DE",parameterLabel)
       write (parameterLabel,'(f6.4)') Omega_b()
       recfastFile=recfastFile//'_Omegab'//trim(parameterLabel)
       call Write_Parameter(parameterDoc,"Omega_b",parameterLabel)
       write (parameterLabel,'(f4.1)') H_0()
       recfastFile=recfastFile//'_H0'//trim(parameterLabel)
       call Write_Parameter(parameterDoc,"H_0",parameterLabel)
       write (parameterLabel,'(f5.3)') T_CMB()
       recfastFile=recfastFile//'_TCMB'//trim(parameterLabel)
       call Write_Parameter(parameterDoc,"T_CMB",parameterLabel)
       write (parameterLabel,'(f4.2)') heliumByMassPrimordial
       recfastFile=recfastFile//'_YHe'//trim(parameterLabel)
       call Write_Parameter(parameterDoc,"Y_He",parameterLabel)
       recfastFile=recfastFile//'.xml'
       call xml_Close(parameterDoc)
       
       ! Run the RecFast driver script to generate the data.
       command='./scripts/aux/RecFast_Driver.pl '//parameterFile//' '//recfastFile
       call System_Command_Do(command)

       ! Set the name of the data file to read in the ``file'' implementation.
       call Intergalactic_Medium_File_Set_File(recfastFile)

    end if
    return
  end subroutine Intergalactic_Medium_State_RecFast_Initialize

end module Intergalactic_Medium_State_RecFast

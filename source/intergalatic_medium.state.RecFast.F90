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

!+    Contributions to this file made by:  Luiz Felippe S. Rodrigues.

!% Contains a module that implements calculations of the intergalactic medium thermal and ionization state using RecFast.

module Intergalactic_Medium_State_RecFast
  !% Implements calculations of the intergalactic medium thermal and ionization state using RecFast.
  use ISO_Varying_String
  implicit none
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
    use Cosmological_Parameters
    use Numerical_Constants_Astronomical
    use Input_Parameters
    use Galacticus_Input_Paths
    implicit none
    type     (varying_string  ), intent(in   )          :: intergalaticMediumStateMethod
    procedure(double precision), intent(inout), pointer :: Intergalactic_Medium_Electron_Fraction_Get
    procedure(double precision), intent(inout), pointer :: Intergalactic_Medium_Temperature_Get
    character(len=32          )                         :: parameterLabel
    type     (varying_string  )                         :: command                                   , parameterFile, &
         &                                                 recfastFile
    type     (xmlf_t          )                         :: parameterDoc

    ! Test if our method has been selected.
    if (intergalaticMediumStateMethod == 'RecFast') then
       ! Set procedure pointers (use those from the ``file'' implementation).
       Intergalactic_Medium_Electron_Fraction_Get => Intergalactic_Medium_Electron_Fraction_File
       Intergalactic_Medium_Temperature_Get       => Intergalactic_Medium_Temperature_File

       ! Ensure that the RecFast data file has been generated.

       ! Generate the name of the data file and an XML input parameter file.
       recfastFile  =char(Galacticus_Input_Path())//'data/intergalacticMedium/recFast'
       parameterFile=char(Galacticus_Input_Path())//'data/intergalacticMedium/recfast_parameters.xml'
       call xml_OpenFile(char(parameterFile),parameterDoc)
       call xml_NewElement(parameterDoc,"parameters")
       write (parameterLabel,'(f5.3)') Omega_Matter()
       recfastFile=recfastFile//'_OmegaM'//trim(parameterLabel)
       call Write_Parameter(parameterDoc,"Omega_Matter",parameterLabel)
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
       write (parameterLabel,'(i1)') Intergalactic_Medium_State_File_Current_File_Format_Version()
       call Write_Parameter(parameterDoc,"fileFormat",parameterLabel)
       call xml_Close(parameterDoc)

       ! Run the RecFast driver script to generate the data.
       command=char(Galacticus_Input_Path())//'scripts/aux/RecFast_Driver.pl '//parameterFile//' '//recfastFile
       call System_Command_Do(command)

       ! Set the name of the data file to read in the ``file'' implementation.
       call Intergalactic_Medium_File_Set_File(recfastFile)

    end if
    return
  end subroutine Intergalactic_Medium_State_RecFast_Initialize

end module Intergalactic_Medium_State_RecFast

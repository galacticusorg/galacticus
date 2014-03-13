!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of the intergalactic medium state class in which state is computed using {\sc RecFast}.

  !# <intergalacticMediumState name="intergalacticMediumStateRecFast">
  !#  <description>The intergalactic medium state is computed using {\sc RecFast}.</description>
  !# </intergalacticMediumState>

  type, extends(intergalacticMediumStateFile) :: intergalacticMediumStateRecFast
     !% An \gls{igm} state class which computes state using {\sc RecFast}.
     private
  end type intergalacticMediumStateRecFast
  
  interface intergalacticMediumStateRecFast
     !% Constructors for the {\sc RecFast} intergalactic medium state class.
     module procedure recFastDefaultConstructor
     module procedure recFastConstructor
  end interface intergalacticMediumStateRecFast

contains

  function recFastDefaultConstructor()
    !% Default constructor for the {\sc RecFast} \gls{igm} state class.
    use Cosmology_Parameters
    implicit none
    type (intergalacticMediumStateRecFast), target  :: recFastDefaultConstructor
    class(cosmologyParametersClass       ), pointer :: cosmologyParameters_

    cosmologyParameters_      => cosmologyParameters(                    )
    recFastDefaultConstructor =  recFastConstructor (cosmologyParameters_)
    return
  end function recFastDefaultConstructor

  function recFastConstructor(thisCosmologyParameters)
    !% Constructor for the {\sc RecFast} \gls{igm} state class.
    use Cosmology_Parameters
    use FoX_wxml
    use System_Command
    use Numerical_Constants_Astronomical
    use Galacticus_Input_Paths
    use Input_Parameters
    implicit none
    type     (intergalacticMediumStateRecFast), target        :: recFastConstructor
    class    (cosmologyParametersClass       ), intent(inout) :: thisCosmologyParameters
    character(len=32                         )                :: parameterLabel
    type     (varying_string                 )                :: command                , parameterFile
    type     (xmlf_t                         )                :: parameterDoc

    ! Generate the name of the data file and an XML input parameter file.
    command='mkdir -p '//char(Galacticus_Input_Path())//'data/intergalacticMedium'
    call System_Command_Do(command)
    recFastConstructor%fileName=char(Galacticus_Input_Path())//'data/intergalacticMedium/recFast'
    parameterFile              =char(Galacticus_Input_Path())//'data/intergalacticMedium/recfast_parameters.xml'
    call xml_OpenFile(char(parameterFile),parameterDoc)
    call xml_NewElement(parameterDoc,"parameters")
    write (parameterLabel,'(f5.3)') thisCosmologyParameters%OmegaMatter()
    recFastConstructor%fileName=recFastConstructor%fileName//'_OmegaM'//trim(parameterLabel)
    call Write_Parameter_XML(parameterDoc,"Omega_Matter",parameterLabel)
    write (parameterLabel,'(f5.3)') thisCosmologyParameters%OmegaDarkEnergy()
    recFastConstructor%fileName=recFastConstructor%fileName//'_OmegaDE'//trim(parameterLabel)
    call Write_Parameter_XML(parameterDoc,"Omega_DE",parameterLabel)
    write (parameterLabel,'(f6.4)') thisCosmologyParameters%OmegaBaryon()
    recFastConstructor%fileName=recFastConstructor%fileName//'_Omegab'//trim(parameterLabel)
    call Write_Parameter_XML(parameterDoc,"Omega_b",parameterLabel)
    write (parameterLabel,'(f4.1)') thisCosmologyParameters%HubbleConstant(unitsStandard)
    recFastConstructor%fileName=recFastConstructor%fileName//'_H0'//trim(parameterLabel)
    call Write_Parameter_XML(parameterDoc,"H_0",parameterLabel)
    write (parameterLabel,'(f5.3)') thisCosmologyParameters%temperatureCMB()
    recFastConstructor%fileName=recFastConstructor%fileName//'_TCMB'//trim(parameterLabel)
    call Write_Parameter_XML(parameterDoc,"T_CMB",parameterLabel)
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    recFastConstructor%fileName=recFastConstructor%fileName//'_YHe'//trim(parameterLabel)
    call Write_Parameter_XML(parameterDoc,"Y_He",parameterLabel)
    recFastConstructor%fileName=recFastConstructor%fileName//'.xml'
    write (parameterLabel,'(i1)') fileFormatVersionCurrent
    call Write_Parameter_XML(parameterDoc,"fileFormat",parameterLabel)
    call xml_Close(parameterDoc)
    ! Run the RecFast driver script to generate the data.
    command=char(Galacticus_Input_Path())//'scripts/aux/RecFast_Driver.pl '//parameterFile//' '//recFastConstructor%fileName
    call System_Command_Do(command)
    return
  end function recFastConstructor

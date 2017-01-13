!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !% An implementation of the intergalactic medium state class in which state is computed using {\normalfont \scshape RecFast}.

  !# <intergalacticMediumState name="intergalacticMediumStateRecFast">
  !#  <description>The intergalactic medium state is computed using {\normalfont \scshape RecFast}.</description>
  !# </intergalacticMediumState>

  type, extends(intergalacticMediumStateFile) :: intergalacticMediumStateRecFast
     !% An \gls{igm} state class which computes state using {\normalfont \scshape RecFast}.
     private
  end type intergalacticMediumStateRecFast
  
  interface intergalacticMediumStateRecFast
     !% Constructors for the {\normalfont \scshape RecFast} intergalactic medium state class.
     module procedure recFastDefaultConstructor
     module procedure recFastConstructor
  end interface intergalacticMediumStateRecFast

contains

  function recFastDefaultConstructor()
    !% Default constructor for the {\normalfont \scshape RecFast} \gls{igm} state class.
    use Cosmology_Parameters
    implicit none
    type (intergalacticMediumStateRecFast), target  :: recFastDefaultConstructor
    class(cosmologyParametersClass       ), pointer :: cosmologyParameters_

    cosmologyParameters_      => cosmologyParameters(                    )
    recFastDefaultConstructor =  recFastConstructor (cosmologyParameters_)
    return
  end function recFastDefaultConstructor

  function recFastConstructor(thisCosmologyParameters)
    !% Constructor for the {\normalfont \scshape RecFast} \gls{igm} state class.
    use Cosmology_Parameters
    use FoX_wxml
    use System_Command
    use Numerical_Constants_Astronomical
    use Galacticus_Input_Paths
    use Input_Parameters
    use Input_Parameters2
    implicit none
    type     (intergalacticMediumStateRecFast), target        :: recFastConstructor
    class    (cosmologyParametersClass       ), intent(inout) :: thisCosmologyParameters
    character(len=32                         )                :: parameterLabel
    type     (varying_string                 )                :: command                , parameterFile
    type     (xmlf_t                         )                :: parameterDoc
    type     (inputParameterList             )                :: parameters

    ! Generate the name of the data file and an XML input parameter file.
    command='mkdir -p '//char(Galacticus_Input_Path())//'data/intergalacticMedium'
    call System_Command_Do(command)
    recFastConstructor%fileName=char(Galacticus_Input_Path())//'data/intergalacticMedium/recFast'
    parameterFile              =char(Galacticus_Input_Path())//'data/intergalacticMedium/recfast_parameters.xml'
    ! Construct a parameter list containing all relevant values.
    parameters=inputParameterList()
    write (parameterLabel,'(f5.3)') thisCosmologyParameters%OmegaMatter    (                   )
    recFastConstructor%fileName=recFastConstructor%fileName//'_OmegaMatter'    //trim(parameterLabel)
    call parameters%add("OmegaMatter"    ,parameterLabel)
    write (parameterLabel,'(f5.3)') thisCosmologyParameters%OmegaDarkEnergy(                   )
    recFastConstructor%fileName=recFastConstructor%fileName//'_OmegaDarkEnergy'//trim(parameterLabel)
    call parameters%add("OmegaDarkEnergy",parameterLabel)
    write (parameterLabel,'(f6.4)') thisCosmologyParameters%OmegaBaryon    (                   )
    recFastConstructor%fileName=recFastConstructor%fileName//'_OmegaBaryon'    //trim(parameterLabel)
    call parameters%add("OmegaBaryon"    ,parameterLabel)
    write (parameterLabel,'(f4.1)') thisCosmologyParameters%HubbleConstant (hubbleUnitsStandard)
    recFastConstructor%fileName=recFastConstructor%fileName//'_HubbleConstant' //trim(parameterLabel)
    call parameters%add("HubbleConstant" ,parameterLabel)
    write (parameterLabel,'(f5.3)') thisCosmologyParameters%temperatureCMB (                   )
    recFastConstructor%fileName=recFastConstructor%fileName//'_temperatureCMB' //trim(parameterLabel)
    call parameters%add("temperatureCMB" ,parameterLabel)
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    recFastConstructor%fileName=recFastConstructor%fileName//'_YHe'            //trim(parameterLabel)
    call parameters%add("Y_He"           ,parameterLabel)
    recFastConstructor%fileName=recFastConstructor%fileName//'.xml'
    write (parameterLabel,'(i1)') fileFormatVersionCurrent
    call parameters%add("fileFormat",parameterLabel)
    ! Generate the parameters XML file.
    call xml_OpenFile(char(parameterFile),parameterDoc)
    call xml_NewElement(parameterDoc,"parameters")
    call parameters%serializeToXML(parameterDoc)
    call xml_Close(parameterDoc)
    ! Run the RecFast driver script to generate the data.
    command=char(Galacticus_Input_Path())//'scripts/aux/RecFast_Driver.pl '//parameterFile//' '//recFastConstructor%fileName
    call System_Command_Do(command)
    return
  end function recFastConstructor

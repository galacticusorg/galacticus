!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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
     module procedure recFastConstructorParameters
     module procedure recFastConstructorInternal
  end interface intergalacticMediumStateRecFast

contains

  function recFastConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \scshape RecFast} \gls{igm} state class.
    use Cosmology_Parameters
    use Input_Parameters
    implicit none
    type (intergalacticMediumStateRecFast)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(cosmologyParametersClass       ), pointer       :: cosmologyParameters_

    ! Check and read parameters.
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    self=intergalacticMediumStateRecFast(cosmologyParameters_)
    !# <objectDestrctor name="cosmologyParameters_"/>
    !# <inputParametersValidate source="parameters"/>
    return
  end function recFastConstructorParameters

  function recFastConstructorInternal(cosmologyParameters_) result(self)
    !% Constructor for the {\normalfont \scshape RecFast} \gls{igm} state class.
    use Cosmology_Parameters
    use FoX_wxml
    use System_Command
    use Numerical_Constants_Astronomical
    use Galacticus_Input_Paths
    use Input_Parameters
    use Input_Parameters
    implicit none
    type     (intergalacticMediumStateRecFast)                :: self
    class    (cosmologyParametersClass       ), intent(inout) :: cosmologyParameters_
    character(len=32                         )                :: parameterLabel
    type     (varying_string                 )                :: command             , parameterFile
    type     (xmlf_t                         )                :: parameterDoc
    type     (inputParameterList             )                :: parameters

    ! Generate the name of the data file and an XML input parameter file.
    command='mkdir -p '//char(Galacticus_Input_Path())//'data/intergalacticMedium'
    call System_Command_Do(command)
    self         %fileName=char(Galacticus_Input_Path())//'data/intergalacticMedium/recFast'
    parameterFile         =char(Galacticus_Input_Path())//'data/intergalacticMedium/recfast_parameters.xml'
    ! Construct a parameter list containing all relevant values.
    parameters=inputParameterList()
    write (parameterLabel,'(f5.3)') cosmologyParameters_%OmegaMatter    (                   )
    self%fileName=self%fileName//'_OmegaMatter'    //trim(parameterLabel)
    call parameters%add("OmegaMatter"    ,parameterLabel)
    write (parameterLabel,'(f5.3)') cosmologyParameters_%OmegaDarkEnergy(                   )
    self%fileName=self%fileName//'_OmegaDarkEnergy'//trim(parameterLabel)
    call parameters%add("OmegaDarkEnergy",parameterLabel)
    write (parameterLabel,'(f6.4)') cosmologyParameters_%OmegaBaryon    (                   )
    self%fileName=self%fileName//'_OmegaBaryon'    //trim(parameterLabel)
    call parameters%add("OmegaBaryon"    ,parameterLabel)
    write (parameterLabel,'(f4.1)') cosmologyParameters_%HubbleConstant (hubbleUnitsStandard)
    self%fileName=self%fileName//'_HubbleConstant' //trim(parameterLabel)
    call parameters%add("HubbleConstant" ,parameterLabel)
    write (parameterLabel,'(f5.3)') cosmologyParameters_%temperatureCMB (                   )
    self%fileName=self%fileName//'_temperatureCMB' //trim(parameterLabel)
    call parameters%add("temperatureCMB" ,parameterLabel)
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    self%fileName=self%fileName//'_YHe'            //trim(parameterLabel)
    call parameters%add("Y_He"           ,parameterLabel)
    self%fileName=self%fileName//'.hdf5'
    write (parameterLabel,'(i1)') fileFormatVersionCurrent
    call parameters%add("fileFormat",parameterLabel)
    ! Generate the parameters XML file.
    call xml_OpenFile(char(parameterFile),parameterDoc)
    call xml_NewElement(parameterDoc,"parameters")
    call parameters%serializeToXML(parameterDoc)
    call xml_Close(parameterDoc)
    ! Run the RecFast driver script to generate the data.
    command=char(Galacticus_Input_Path())//'scripts/aux/RecFast_Driver.pl '//parameterFile//' '//self%fileName
    call System_Command_Do(command)
    return
  end function recFastConstructorInternal

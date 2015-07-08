!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of the accretion disk spectra class using the model of \cite{hopkins_observational_2007}.

  !# <accretionDiskSpectra name="accretionDiskSpectraHopkins2007">
  !#  <description>Accretion disk spectra using the model of \cite{hopkins_observational_2007}.</description>
  !# </accretionDiskSpectra>
  type, extends(accretionDiskSpectraFile) :: accretionDiskSpectraHopkins2007
     !% An accretion disk spectra class which uses the algorithm of \cite{hopkins_observational_2007}.
     private
   contains
     procedure :: descriptor => hopkins2007Descriptor
  end type accretionDiskSpectraHopkins2007
  
  interface accretionDiskSpectraHopkins2007
     !% Constructors for the {\normalfont \ttfamily hopkins2007} accretion disk spectra class.
     module procedure hopkins2007ConstructorParameters
     module procedure hopkins2007ConstructorInternal
  end interface accretionDiskSpectraHopkins2007

contains

  function hopkins2007ConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily hopkins2007} accretion disk spectra class.
    use Input_Parameters2
    implicit none
    type(accretionDiskSpectraHopkins2007)                :: hopkins2007ConstructorParameters
    type(inputParameters                ), intent(in   ) :: parameters

    hopkins2007ConstructorParameters=hopkins2007ConstructorInternal()
    return
  end function hopkins2007ConstructorParameters

  function hopkins2007ConstructorInternal()
    !% Constructor for the {\normalfont \ttfamily hopkins2007} accretion disk spectra class.
    use System_Command
    use Galacticus_Input_Paths
    use String_Handling
    implicit none
    type(accretionDiskSpectraHopkins2007) :: hopkins2007ConstructorInternal
    type(varying_string                 ) :: command
 
    ! Ensure that the data file is generated.
    command=Galacticus_Input_Path()//'scripts/aux/Hopkins2007_AGN_SED_Driver.pl '//fileFormatCurrent
    call System_Command_Do(command)
    ! Load the file.
    call hopkins2007ConstructorInternal%loadFile(char(Galacticus_Input_Path()//"/data/blackHoles/AGN_SEDs_Hopkins2007.hdf5"))
    ! Initialize interpolators.
    hopkins2007ConstructorInternal%resetLuminosity=.true.
    hopkins2007ConstructorInternal%resetWavelength=.true.
    return
  end function hopkins2007ConstructorInternal

  subroutine hopkins2007Descriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    implicit none
    class(accretionDiskSpectraHopkins2007), intent(inout) :: self
    type (inputParameters                ), intent(inout) :: descriptor

    call descriptor%addParameter("accretionDiskSpectraMethod","hopkins2007")
    return
  end subroutine hopkins2007Descriptor

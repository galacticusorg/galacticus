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

  !% An implementation of the accretion disk spectra class for tabulated spectra read from file.

  !# <accretionDiskSpectra name="accretionDiskSpectraHopkins2007">
  !#  <description>Accretion disk spectra using the model of \cite{hopkins_observational_2007}.</description>
  !# </accretionDiskSpectra>

  type, extends(accretionDiskSpectraFile) :: accretionDiskSpectraHopkins2007
     !% An accretion disk spectra class which uses the algorithm of \cite{hopkins_observational_2007}.
     private
   contains
  end type accretionDiskSpectraHopkins2007

  interface accretionDiskSpectraHopkins2007
     !% Constructors for the {\normalfont \ttfamily hopkins2007} accretion disk spectra class.
     module procedure hopkins2007Constructor
  end interface accretionDiskSpectraHopkins2007

contains

  function hopkins2007Constructor()
    !% Constructor for the {\normalfont \ttfamily hopkins2007} accretion disk spectra class.
    use System_Command
    use Galacticus_Input_Paths
    use String_Handling
    implicit none
    type(accretionDiskSpectraHopkins2007), target :: hopkins2007Constructor
    type(varying_string                 )         :: command
 
    ! Ensure that the data file is generated.
    command=Galacticus_Input_Path()//'scripts/aux/Hopkins2007_AGN_SED_Driver.pl '//fileFormatCurrent
    call System_Command_Do(command)
    ! Load the file.
    call hopkins2007Constructor%loadFile(char(Galacticus_Input_Path()//"/data/blackHoles/AGN_SEDs_Hopkins2007.hdf5"))
    ! Initialize interpolators.
    hopkins2007Constructor%resetLuminosity=.true.
    hopkins2007Constructor%resetWavelength=.true.
   return
  end function hopkins2007Constructor

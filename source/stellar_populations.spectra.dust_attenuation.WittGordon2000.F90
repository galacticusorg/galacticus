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

  !% Implements calculations of attenuation of stellar spectra using the model of \cite{witt_multiple_2000}.
  
  !# <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationWittGordon2000">
  !#  <description>Returns the dust attenuation of stellar spectra according to the model of \cite{witt_multiple_2000}.</description>
  !# </stellarSpectraDustAttenuation>

  use Tables

  type, extends(stellarSpectraDustAttenuationTabulated) :: stellarSpectraDustAttenuationWittGordon2000
     !% A class implementing calculations of attenuation of stellar spectra using the model of \cite{witt_multiple_2000}.
     private
   contains
  end type stellarSpectraDustAttenuationWittGordon2000

  interface stellarSpectraDustAttenuationWittGordon2000
     !% Constructors for the ``wittGordon2003'' stellar spectra dust attenuation class.
     module procedure wittGordon2003DefaultConstructor
     module procedure wittGordon2003Constructor
  end interface stellarSpectraDustAttenuationWittGordon2000

  ! Name of the model to use by default.
  type(varying_string) :: wittGordon2003Model

contains

  function wittGordon2003DefaultConstructor()
    !% Default constructor for the ``wittGordon2003'' stellar spectra dust attenuation class.
    use Input_Parameters
    implicit none
    type(stellarSpectraDustAttenuationWittGordon2000) :: wittGordon2003DefaultConstructor

    !@ <inputParameter>
    !@   <name>dustAttenuationWittGordon2000Model</name>
    !@   <defaultValue>SMCbar</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The name of the model from \cite{witt_multiple_2000} to use in dust attenuation calculations.
    !@   </description>
    !@   <type>string</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('dustAttenuationWittGordon2000Model',wittGordon2003Model,defaultValue='MilkyWayShellTau3.0')
    wittGordon2003DefaultConstructor=wittGordon2003Constructor(char(wittGordon2003Model))
    return
  end function wittGordon2003DefaultConstructor

  function wittGordon2003Constructor(model)
    !% Constructor for the ``wittGordon2003'' stellar spectra dust attenuation class.
    use Galacticus_Error
    use Numerical_Constants_Units
    use Numerical_Constants_Astronomical
    implicit none
    type     (stellarSpectraDustAttenuationWittGordon2000)                :: wittGordon2003Constructor
    character(len=*                                      ), intent(in   ) :: model

    ! Initialize fitting function parameters for the chosen model.
    call wittGordon2003Constructor%attenuationTable%create(angstromsPerMicron/[ 1000.0d0,   1142.0d0,   1285.0d0,   1428.0d0,   1571.0d0,  1714.0d0,  1857.0d0,  2000.0d0,   2142.0d0,  2285.0d0,   2428.0d0,   2571.0d0,   2714.0d0,   2857.0d0,   3000.0d0,   3776.0d0,   4754.0d0,   5985.0d0,   7535.0d0,   9487.0d0,  11943.0d0,  15036.0d0,  18929.0d0,  23830.0d0,  30001.0d0],tableCount=1,extrapolationType=extrapolationTypeExtrapolate)
     select case (model)
     case ("MilkyWayShellTau3.0")
        call wittGordon2003Constructor%attenuationTable%populate(magnitudesPerOpticalDepth*[ 15.714d0,  11.754d0,   9.546d0,   8.340d0,   7.752d0,   7.527d0,   7.683d0,   8.529d0,   9.570d0,   8.730d0,   7.416d0,   6.582d0,   6.066d0,   5.715d0,   5.454d0,   4.581d0,   3.597d0,   2.727d0,   2.001d0,   1.320d0,   0.912d0,   0.630d0,   0.435d0,   0.300d0,   0.207d0])
     case ("SMCShellTau3.0"    )
        call wittGordon2003Constructor%attenuationTable%populate(magnitudesPerOpticalDepth*[ 29.025d0,  22.320d0,  18.204d0,  15.501d0,  13.608d0,  12.222d0,  11.100d0,  10.137d0,   9.303d0,   8.571d0,   7.926d0,   7.356d0,   6.846d0,   6.399d0,   6.093d0,   4.830d0,   3.663d0,   2.640d0,   1.890d0,   1.290d0,   0.816d0,   0.498d0,   0.333d0,   0.225d0,   0.099d0])
     case default
        call Galacticus_Error_Report('wittGordon2003Constructor','model unrecognized')
     end select
    return
  end function wittGordon2003Constructor

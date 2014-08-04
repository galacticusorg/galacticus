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

  !% Implements calculations of attenuation of stellar spectra using the model of \cite{gordon_quantitative_2003}.
  
  !# <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationGordon2003">
  !#  <description>Returns the dust attenuation of stellar spectra according to the model of \cite{gordon_quantitative_2003}.</description>
  !# </stellarSpectraDustAttenuation>

  use Tables

  type, extends(stellarSpectraDustAttenuationClass) :: stellarSpectraDustAttenuationGordon2003
     !% A class implementing calculations of attenuation of stellar spectra using the model of \cite{gordon_quantitative_2003}.
     private
     type(table1DGeneric) :: attenuationTable
   contains
     final     ::                gordon2003Destructor
     procedure :: attenuation => gordon2003Attenuation
  end type stellarSpectraDustAttenuationGordon2003

  interface stellarSpectraDustAttenuationGordon2003
     !% Constructors for the ``gordon2003'' stellar spectra dust attenuation class.
     module procedure gordon2003DefaultConstructor
     module procedure gordon2003Constructor
  end interface stellarSpectraDustAttenuationGordon2003

  ! Name of the sample to use by default.
  type(varying_string) :: gordon2003Sample

contains

  function gordon2003DefaultConstructor()
    !% Default constructor for the ``gordon2003'' stellar spectra dust attenuation class.
    use Input_Parameters
    implicit none
    type            (stellarSpectraDustAttenuationGordon2003)            :: gordon2003DefaultConstructor

    !@ <inputParameter>
    !@   <name>dustAttenuationGordon2003Sample</name>
    !@   <defaultValue>SMCbar</defaultValue>
    !@   <attachedTo>module</attachedTo>
    !@   <description>
    !@     The name of the sample from \cite{gordon_quantitative_2003} to use in dust attenuation calculations.
    !@   </description>
    !@   <type>string</type>
    !@   <cardinality>1</cardinality>
    !@ </inputParameter>
    call Get_Input_Parameter('dustAttenuationGordon2003Sample',gordon2003Sample,defaultValue='SMCbar')
    gordon2003DefaultConstructor=gordon2003Constructor(char(gordon2003Sample))
    return
  end function gordon2003DefaultConstructor

  function gordon2003Constructor(sample)
    !% Constructor for the ``gordon2003'' stellar spectra dust attenuation class.
    use Galacticus_Error
    implicit none
    type            (stellarSpectraDustAttenuationGordon2003)                :: gordon2003Constructor
    character       (len=*                                  ), intent(in   ) :: sample

    ! Initialize fitting function parameters for the chosen sample (values from Tables 2 & 3 of Gordon et al.).
     select case (sample)
     case ("SMCbar")
        call gordon2003Constructor%attenuationTable%create([0.455d0,0.606d0,0.800d0,1.235d0,1.538d0,1.818d0,2.273d0,2.703d0,3.375d0,3.625d0,3.875d0,4.125d0,4.375d0,4.625d0,4.875d0,5.125d0,5.375d0,5.625d0,5.875d0,6.125d0,6.375d0,6.625d0,6.875d0,7.125d0,7.375d0,7.625d0,7.875d0,8.125d0,8.375d0,8.625d0],tableCount=1,extrapolationType=extrapolationTypeExtrapolate)
        call gordon2003Constructor%attenuationTable%populate([0.016d0,0.169d0,0.131d0,0.567d0,0.801d0,1.000d0,1.374d0,1.672d0,2.000d0,2.220d0,2.428d0,2.661d0,2.947d0,3.161d0,3.293d0,3.489d0,3.637d0,3.866d0,4.013d0,4.243d0,4.472d0,4.776d0,5.000d0,5.272d0,5.575d0,5.795d0,6.074d0,6.297d0,6.436d0,6.992d0])
     case ("LMC"   )
        call gordon2003Constructor%attenuationTable%create([0.455d0,0.606d0,0.800d0,1.818d0,2.273d0,2.703d0,3.375d0,3.625d0,3.875d0,4.125d0,4.375d0,4.625d0,4.875d0,5.125d0,5.375d0,5.625d0,5.875d0,6.125d0,6.375d0,6.625d0,6.875d0,7.125d0,7.375d0,7.625d0,7.875d0,8.125d0,8.375d0],1,extrapolationType=extrapolationTypeExtrapolate)
        call gordon2003Constructor%attenuationTable%populate([0.030d0,0.186d0,0.257d0,1.000d0,1.293d0,1.518d0,1.786d0,1.969d0,2.149d0,2.391d0,2.771d0,2.967d0,2.846d0,2.646d0,2.565d0,2.566d0,2.598d0,2.607d0,2.668d0,2.787d0,2.874d0,2.983d0,3.118d0,3.231d0,3.374d0,3.366d0,3.467d0])
     case default
        call Galacticus_Error_Report('gordon2003Constructor','sample must be one of "SMCbar", "SMCwing", and "LMC"')
     end select
    return
  end function gordon2003Constructor

  subroutine gordon2003Destructor(self)
    !% Destructor for the ``gordon2003'' stellar spectra dust attenuation class.
    use Input_Parameters
    implicit none
    type(stellarSpectraDustAttenuationGordon2003), intent(inout) :: self

    call self%attenuationTable%destroy()
    return
  end subroutine gordon2003Destructor

  double precision function gordon2003Attenuation(self,wavelength,age,vBandAttenuation)
    !% Return attenuation of stellar spectra according to the model of \cite{gordon_quantitative_2003}.
    use Numerical_Constants_Units
    implicit none
    class           (stellarSpectraDustAttenuationGordon2003), intent(inout) :: self
    double precision                                         , intent(in   ) :: wavelength      , age, &
         &                                                                      vBandAttenuation
    double precision                                                         :: x

    x=1.0d0/(wavelength/angstromsPerMicron)
    gordon2003Attenuation=vBandAttenuation*self%attenuationTable%interpolate(x)
    return
  end function gordon2003Attenuation

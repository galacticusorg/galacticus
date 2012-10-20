!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program which computes various quantities related to the excursion set formalism and stores them in an output file.

program Tests_Excursion_Sets
  !% Computes various quantities related to the excursion set formalism and stores them in an output file.
  use Input_Parameters
  use ISO_Varying_String
  use Memory_Management
  use Cosmology_Functions
  use Cosmological_Parameters
  use Excursion_Sets_Barriers
  use Excursion_Sets_First_Crossings
  use Halo_Mass_Function
  use Galacticus_Display
  use Galacticus_Error
  use Numerical_Ranges
  use CDM_Power_Spectrum
  use Numerical_Constants_Math
  use IO_HDF5
  implicit none
  integer,                             parameter                   :: massCount=200
  double precision,                    parameter                   :: massMinimum=1.0d6,massMaximum=1.0d16
  integer,                             parameter                   :: fileNameLengthMaximum=1024
  double precision,                    allocatable, dimension(:  ) :: haloMass,variance,barrier,firstCrossingProbability,haloMassFunction,wavenumber,powerSpectrum
  double precision,                    allocatable, dimension(:,:) :: firstCrossingRate
  integer                                                          :: iMass,jMass
  double precision                                                 :: time,varianceProgenitor
  character(len=fileNameLengthMaximum)                             :: fileCharacter
  type(varying_string)                                             :: parameterFile,outputFileName
  type(hdf5Object)                                                 :: outputFile

  ! Read in basic code memory usage.
  call Code_Memory_Usage('Excursion_Sets.size')

  ! Get the name of the parameter file from the first command line argument.
  call Get_Command_Argument(1,fileCharacter)
  if (len_trim(fileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Excursion_Sets.exe <parameterFile> <outputFile>")
  parameterFile=fileCharacter
  call Get_Command_Argument(2,fileCharacter)
  if (len_trim(fileCharacter) == 0) call Galacticus_Error_Report(message="Usage: Excursion_Sets.exe <parameterFile> <outputFile>")
  outputFileName=fileCharacter

  ! Open the output file.
  call outputFile%openFile(char(outputFileName),overWrite=.true.,objectsOverwritable=.false.)

  ! Read the parameter file.
  call Input_Parameters_File_Open(parameterFile,outputFile,allowedParametersFile='Excursion_Sets.parameters.xml')

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(1)

  ! Get the current time.
  time=Cosmology_Age(1.0d0)

  ! Allocate arrays.
  call Alloc_Array(haloMass                ,[massCount          ])
  call Alloc_Array(variance                ,[massCount          ])
  call Alloc_Array(barrier                 ,[massCount          ])
  call Alloc_Array(firstCrossingProbability,[massCount          ])
  call Alloc_Array(haloMassFunction        ,[massCount          ])
  call Alloc_Array(wavenumber              ,[massCount          ])
  call Alloc_Array(powerSpectrum           ,[massCount          ])
  call Alloc_Array(firstCrossingRate       ,[massCount,massCount])

  ! Create a grid of masses.
  haloMass=Make_Range(massMinimum,massMaximum,massCount,rangeType=rangeTypeLogarithmic)

  ! Set first crossing rates to unphysical values.
  firstCrossingRate=-1.0d0

  ! Loop over masses.
  do iMass=1,massCount
     wavenumber              (iMass)=(3.0d0*haloMass(iMass)/4.0d0/Pi/Critical_Density()/Omega_Matter())**(-1.0d0/3.0d0)
     powerSpectrum           (iMass)=Power_Spectrum_CDM                       (wavenumber   (iMass))
     variance                (iMass)=sigma_CDM                                (     haloMass(iMass))**2
     barrier                 (iMass)=Excursion_Sets_Barrier                   (variance(iMass),time)
     firstCrossingProbability(iMass)=Excursion_Sets_First_Crossing_Probability(variance(iMass),time)
     haloMassFunction        (iMass)=Halo_Mass_Function_Differential          (time,haloMass(iMass))

     ! Compute halo branching rates.
     do jMass=1,iMass-1
        varianceProgenitor=sigma_CDM(haloMass(jMass))**2
        firstCrossingRate(iMass,jMass)=Excursion_Sets_First_Crossing_Rate(variance(iMass),varianceProgenitor,time)
     end do

  end do

  ! Write results to the output file.
  call outputFile%writeDataset(haloMass                ,'haloMass'                ,'The mass of the halo [M☉]'                       )
  call outputFile%writeDataset(wavenumber              ,'wavenumber'              ,'The wavenumber associated with this mass [Mpc⁻¹]')
  call outputFile%writeDataset(powerSpectrum           ,'powerSpectrum'           ,'The power spectrum at this wavenumber [Mpc³]'    )
  call outputFile%writeDataset(variance                ,'variance'                ,'The variance on this mass scale'                 )
  call outputFile%writeDataset(barrier                 ,'barrier'                 ,'The excursion set barrier for this variance'     )
  call outputFile%writeDataset(firstCrossingProbability,'firstCrossingProbability','The first crossing probability'                  )
  call outputFile%writeDataset(haloMassFunction        ,'haloMassFunction'        ,'The halo mass function [Mpc⁻³ M☉⁻¹]'             )
  call outputFile%writeDataset(firstCrossingRate       ,'firstCrossingRate'       ,'The first crossing rate [Gyr⁻¹]'                 )

  ! Close the parameters file.
  call Input_Parameters_File_Close()

  ! Close the output file.
  call outputFile%close()

end program Tests_Excursion_Sets

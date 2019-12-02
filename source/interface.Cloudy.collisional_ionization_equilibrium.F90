!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which provides an interface to the \gls{cloudy} code for computing tables of cooling functions and chemical state in collisional ionization equilibrium.

module Interfaces_Cloudy_CIE
  !% Provides an interface to the \gls{cloudy} code for computing tables of cooling functions and chemical state in collisional ionization equilibrium.
  use :: File_Utilities, only : lockDescriptor
  private
  public :: Interface_Cloudy_CIE_Tabulate

  ! Global file locks.
  type(lockDescriptor) :: fileLockCoolingFunction, fileLockChemicalState

contains

  subroutine Interface_Cloudy_CIE_Tabulate(metallicityMaximumLogarithmic,fileNameCoolingFunction,fileNameChemicalState,versionFileFormat)
    !% An interface to the \gls{cloudy} code for computing tables of cooling functions and chemical state in collisional ionization equilibrium.
    use :: File_Utilities    , only : File_Exists                , File_Lock                       , File_Remove              , File_Unlock
    use :: Galacticus_Display, only : Galacticus_Display_Counter , Galacticus_Display_Counter_Clear, Galacticus_Display_Indent, Galacticus_Display_Message, &
          &                           Galacticus_Display_Unindent, verbosityWorking
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: IO_HDF5           , only : hdf5Access                 , hdf5Object
    use :: ISO_Varying_String, only : var_str                    , varying_string                  , operator(//)             , char                      , &
         &                            assignment(=)
    use :: Interfaces_Cloudy , only : Interface_Cloudy_Initialize
    use :: Numerical_Ranges  , only : Make_Range                 , rangeTypeLinear
    use :: String_Handling   , only : operator(//)
    use :: System_Command    , only : System_Command_Do
    implicit none
    double precision                , intent(in   )                 :: metallicityMaximumLogarithmic
    type            (varying_string), intent(in   )                 :: fileNameCoolingFunction               , fileNameChemicalState
    integer                         , intent(in   )                 :: versionFileFormat
    integer                         , parameter                     :: versionFileFormatCurrent     = 1
    double precision                , parameter                     :: temperatureMinimumLogarithmic=+2.500d0, temperatureMaximumLogarithmic=+9.000d0, &
         &                                                             temperatureStepLogarithmic   =+0.025d0
    double precision                , parameter                     :: metallicityMinimumLogarithmic=-4.000d0, metallicityZeroLogarithmic   =-9.990d2, &
         &                                                             metallicityStepLogarithmic   =+0.250d0
    double precision                , parameter                     :: heliumAbundancePrimordial    =+0.072d0, heliumAbundanceSolar         =+0.100d0 ! Values as used by Cloudy.
    double precision                , allocatable  , dimension(:  ) :: metallicitiesLogarithmic              , temperaturesLogarithmic
    double precision                , allocatable  , dimension(:,:) :: coolingFunction                       , densityElectron                       , &
         &                                                             densityHydrogenI                      , densityHydrogenII
    logical                                                         :: computeCoolingFunctions               , computeChemicalStates
    type            (hdf5Object    )                                :: outputFile                            , dataset
    integer                                                         :: fileFormatFile                        , metallicityCount                     , &
         &                                                             temperatureCount                      , cloudyScript                         , &
         &                                                             iMetallicity                          , inputFile                            , &
         &                                                             iTemperature                          , status
    type            (varying_string)                                :: cloudyPath                            , cloudyVersion                        , &
         &                                                             fileNameTempCooling                   , fileNameTempOverview
    character       (len=8         )                                :: label
    double precision                                                :: dummy                                 , abundanceHelium

    !$omp critical(cloudyCIEFileLock)
    ! Ensure the requested file format version is compatible.
    if (versionFileFormat /= versionFileFormatCurrent) call Galacticus_Error_Report(var_str("this interface supports file format version ")//versionFileFormatCurrent//" but version "//versionFileFormat//" was requested"//{introspection:location})
    ! Determine if we need to compute cooling functions.
    computeCoolingFunctions=.false.
    if (File_Exists(fileNameCoolingFunction)) then
       call hdf5Access%set()
       call outputFile%openFile(char(fileNameCoolingFunction),readOnly=.true.)
       if (outputFile%hasAttribute('fileFormat')) then
          call outputFile%readAttribute('fileFormat',fileFormatFile)
          if (fileFormatFile /= versionFileFormatCurrent) computeCoolingFunctions=.true.
       end if
       call outputFile%close()
       call hdf5Access%unset()
    else
       computeCoolingFunctions=.true.
    end if
    ! Determine if we need to compute chemical states.
    computeChemicalStates=.false.
    if (File_Exists(fileNameChemicalState)) then
       call hdf5Access%set()
       call outputFile%openFile(char(fileNameChemicalState),readOnly=.true.)
       if (outputFile%hasAttribute('fileFormat')) then
          call outputFile%readAttribute('fileFormat',fileFormatFile)
          if (fileFormatFile /= versionFileFormatCurrent) computeChemicalStates=.true.
       end if
       call outputFile%close()
       call hdf5Access%unset()
    else
       computeChemicalStates=.true.
    end if
    ! Perform calculations if necessary.
    if (computeCoolingFunctions .or. computeChemicalStates) then
       ! Open and lock the cooling function and chemical state files.
       call File_Lock(char(fileNameCoolingFunction),fileLockCoolingFunction)
       call File_Lock(char(fileNameChemicalState  ),fileLockChemicalState  )
       ! Generate metallicity and temperature arrays.
       metallicityCount        =int((+metallicityMaximumLogarithmic-metallicityMinimumLogarithmic)/metallicityStepLogarithmic+1.5d0)
       temperatureCount        =int((+temperatureMaximumLogarithmic-temperatureMinimumLogarithmic)/temperatureStepLogarithmic+1.5d0)
       allocate(metallicitiesLogarithmic(metallicityCount+1                 ))
       allocate(temperaturesLogarithmic (                   temperatureCount))
       allocate(coolingFunction         (metallicityCount+1,temperatureCount))
       allocate(densityElectron         (metallicityCount+1,temperatureCount))
       allocate(densityHydrogenI        (metallicityCount+1,temperatureCount))
       allocate(densityHydrogenII       (metallicityCount+1,temperatureCount))
       metallicitiesLogarithmic(1                   )=metallicityZeroLogarithmic
       metallicitiesLogarithmic(2:metallicityCount+1)=Make_Range(metallicityMinimumLogarithmic,metallicityMaximumLogarithmic,metallicityCount,rangeTypeLinear)
       temperaturesLogarithmic                       =Make_Range(temperatureMinimumLogarithmic,temperatureMaximumLogarithmic,temperatureCount,rangeTypeLinear)
       ! Initialize Cloudy.
       call Interface_Cloudy_Initialize(cloudyPath,cloudyVersion)
       ! Specify file names for temporary Cloudy data.
       fileNameTempCooling ="cloudy_cooling.tmp"
       fileNameTempOverview="cloudy_overview.tmp"
       ! Begin iterating over metallicities.
       call Galacticus_Display_Indent("Computing cooling functions and chemical states using Cloudy (this may take a long time)....",verbosityWorking)
       do iMetallicity=1,metallicityCount+1
          if (metallicitiesLogarithmic(iMetallicity) <= metallicityZeroLogarithmic) then
             write (label,'(a)'   ) '  -∞'
          else
             write (label,'(f8.3)') metallicitiesLogarithmic(iMetallicity)
          end if
          call Galacticus_Display_Message("Computing for log₁₀(Z/Z☉)="//trim(label),verbosityWorking)
          call Galacticus_Display_Counter(int(100.0d0*dble(iMetallicity-1)/dble(metallicityCount+1)),iMetallicity==0,verbosityWorking)
          ! Generate an input file for Cloudy.
          open(newUnit=cloudyScript,file=char(cloudyPath//'/source/input.in'),status='unknown',form='formatted')
          write (cloudyScript,'(a)') 'print off'
          write (cloudyScript,'(a)') 'background, z=0'            ! Use a very low level incident continuum.
          write (cloudyScript,'(a)') 'cosmic rays background'     ! Include cosmic ray background ionization rate.
          write (cloudyScript,'(a)') 'stop zone 1'                ! Stop after a single zone.
          write (cloudyScript,'(a)') 'no photoionization'         ! Do three iterations to ensure convergence is reached.
          write (cloudyScript,'(a)') 'hden 0.0'
          if (metallicitiesLogarithmic(iMetallicity) <= metallicityZeroLogarithmic) then
             write (cloudyScript,'(a)') 'abundances primordial'
          else
             write (cloudyScript,'(a,f6.3)') 'metals _log ',metallicitiesLogarithmic(iMetallicity)
             ! Assume a linear growth of helium abundance with metallicity.
             abundanceHelium=heliumAbundancePrimordial+(heliumAbundanceSolar-heliumAbundancePrimordial)*(10.0d0**metallicitiesLogarithmic(iMetallicity))
             write (cloudyScript,'(a,f6.3)') 'element abundance linear helium ',abundanceHelium
          end if
          write (cloudyScript,'(a,f6.3,a)') 'constant temper ',temperatureMinimumLogarithmic,' vary'
          write (cloudyScript,'(a,f6.3,a,f6.3,a,f6.3)') 'grid ',temperatureMinimumLogarithmic,' to ',temperatureMaximumLogarithmic,' step ',temperatureStepLogarithmic
          write (cloudyScript,'(a)') 'no molecules'
          write (cloudyScript,'(a)') 'set trim -20'
          write (cloudyScript,'(a)') 'punch cooling "' //char(fileNameTempCooling )//'"'
          write (cloudyScript,'(a)') 'punch overview "'//char(fileNameTempOverview)//'"'
          close(cloudyScript)
          call System_Command_Do("cd "//cloudyPath//"/source; cloudy.exe -r input",status);
          if (status /= 0) call Galacticus_Error_Report('Cloudy failed'//{introspection:location})
          ! Extract the cooling rate.
          open(newUnit=inputFile,file=char(cloudyPath//"/source/"//fileNameTempCooling),status='old')
          read (inputFile,*) ! Skip the header line.
          do iTemperature=1,temperatureCount
             read (inputFile,*) dummy,dummy,dummy,coolingFunction(iMetallicity,iTemperature)
             read (inputFile,*)
          end do
          close(inputFile)
          call File_Remove(fileNameTempCooling)
          ! Extract the electron and hydrogen density.
          open(newUnit=inputFile,file=char(cloudyPath//"/source/"//fileNameTempOverview),status='old')
          read (inputFile,*) ! Skip the header line.
          do iTemperature=1,temperatureCount
             read (inputFile,*) dummy,dummy,dummy,dummy,densityElectron(iMetallicity,iTemperature),dummy,densityHydrogenI(iMetallicity,iTemperature),densityHydrogenII(iMetallicity,iTemperature)
             read (inputFile,*)
          end do
          close(inputFile)
          call File_Remove(fileNameTempOverview)
       end do
       call Galacticus_Display_Counter_Clear(verbosityWorking)
       ! Output cooling functions to an HDF5 file.
       if (computeCoolingFunctions) then
          call hdf5Access%set()
          call outputFile%openFile      (char(fileNameCoolingFunction))
          ! Store data.
          call outputFile%writeDataset  (metallicitiesLogarithmic                                  ,'metallicity'    ,datasetReturned=dataset)
          call dataset   %writeAttribute('fix'                                                     ,'extrapolateLow'                         )
          call dataset   %writeAttribute('fix'                                                     ,'extrapolateHigh'                        )
          call dataset   %writeAttribute('K'                                                       ,'units'                                  )
          call dataset   %writeAttribute(1.0d0                                                     ,'unitsInSI'                              )
          call dataset   %close         (                                                                                                    )
          call outputFile%writeDataset  (10.0d0**temperaturesLogarithmic                           ,'temperature'    ,datasetReturned=dataset)
          call dataset   %writeAttribute('powerLaw'                                                ,'extrapolateLow'                         )
          call dataset   %writeAttribute('powerLaw'                                                ,'extrapolateHigh'                        )
          call dataset   %close         (                                                                                                    )
          call outputFile%writeDataset  (coolingFunction                                           ,'coolingRate'    ,datasetReturned=dataset)
          call dataset   %close         (                                                                                                    )
          ! Add attributes.
          call outputFile%writeAttribute("CIE cooling functions computed by Cloudy "//cloudyVersion,'description'                            )
          call outputFile%writeAttribute(versionFileFormatCurrent                                  ,'fileFormat'                             )
          call outputFile%close         (                                                                                                    )
          call hdf5Access%unset()
       end if
       ! Output chemical states to an HDF5 file.
       if (computeChemicalStates) then
          call hdf5Access%set()
          call outputFile%openFile      (char(fileNameChemicalState))
          ! Store data.
          call outputFile%writeDataset  (metallicitiesLogarithmic                                  ,'metallicity'    ,datasetReturned=dataset)
          call dataset   %writeAttribute('fix'                                                     ,'extrapolateLow'                         )
          call dataset   %writeAttribute('fix'                                                     ,'extrapolateHigh'                        )
          call dataset   %writeAttribute('K'                                                       ,'units'                                  )
          call dataset   %writeAttribute(1.0d0                                                     ,'unitsInSI'                              )
          call dataset   %close         (                                                                                                    )
          call outputFile%writeDataset  (10.0d0**temperaturesLogarithmic                           ,'temperature'    ,datasetReturned=dataset)
          call dataset   %writeAttribute('powerLaw'                                                ,'extrapolateLow'                         )
          call dataset   %writeAttribute('powerLaw'                                                ,'extrapolateHigh'                        )
          call dataset   %close         (                                                                                                    )
          call outputFile%writeDataset  (densityElectron                                           ,'electronDensity',datasetReturned=dataset)
          call dataset   %close         (                                                                                                    )
          call outputFile%writeDataset  (densityHydrogenI                                          ,'hiDensity'      ,datasetReturned=dataset)
          call dataset   %close         (                                                                                                    )
          call outputFile%writeDataset  (densityHydrogenII                                         ,'hiiDensity'     ,datasetReturned=dataset)
          call dataset   %close         (                                                                                                    )
          ! Add attributes.
          call outputFile%writeAttribute("CIE ionization states computed by Cloudy "//cloudyVersion,'description'                            )
          call outputFile%writeAttribute(versionFileFormatCurrent                                  ,'fileFormat'                             )
          call outputFile%close         (                                                                                                    )
          call hdf5Access%unset()
       end if
       call File_Unlock(fileLockChemicalState  )
       call File_Unlock(fileLockCoolingFunction)
       ! Write message.
       call Galacticus_Display_Unindent("...done",verbosityWorking)
    end if
    !$omp end critical(cloudyCIEFileLock)
    return
  end subroutine Interface_Cloudy_CIE_Tabulate

end module Interfaces_Cloudy_CIE

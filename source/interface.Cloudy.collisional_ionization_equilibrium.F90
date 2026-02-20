!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module which provides an interface to the \gls{cloudy} code for computing tables of cooling functions and chemical state in collisional ionization equilibrium.
!!}

module Interfaces_Cloudy_CIE
  !!{
  Provides an interface to the \gls{cloudy} code for computing tables of cooling functions and chemical state in collisional ionization equilibrium.
  !!}
  use :: File_Utilities, only : lockDescriptor
  private
  public :: Interface_Cloudy_CIE_Tabulate

  ! Global file locks.
  type(lockDescriptor) :: fileLockCoolingFunction, fileLockChemicalState

contains

  subroutine Interface_Cloudy_CIE_Tabulate(metallicityMaximumLogarithmic,fileNameCoolingFunction,fileNameChemicalState,versionFileFormat,includeContinuum)
    !!{
    An interface to the \gls{cloudy} code for computing tables of cooling functions and chemical state in collisional ionization equilibrium.
    !!}
    use :: Display                         , only : displayCounter                     , displayCounterClear           , displayIndent       , displayMessage, &
          &                                         displayUnindent                    , verbosityLevelWorking
    use :: File_Utilities                  , only : File_Exists                        , File_Lock                     , File_Remove         , File_Unlock   , &
         &                                          Directory_Make                     , File_Path
    use :: Error                           , only : Error_Report
    use :: Hashes_Cryptographic            , only : Hash_MD5
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: ISO_Varying_String              , only : assignment(=)                      , char                          , operator(//)        , var_str       , &
          &                                         varying_string
    use :: Input_Paths                     , only : inputPath                          , pathTypeDataDynamic
    use :: Interfaces_Cloudy               , only : Interface_Cloudy_Initialize
    use :: Numerical_Constants_Astronomical, only : heliumToHydrogenAbundancePrimordial, heliumToHydrogenAbundanceSolar
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Units       , only : electronVolt
    use :: Numerical_Ranges                , only : Make_Range                         , rangeTypeLinear               , rangeTypeLogarithmic
    use :: String_Handling                 , only : operator(//)
    use :: System_Command                  , only : System_Command_Do
    implicit none
    double precision                , intent(in   )                   :: metallicityMaximumLogarithmic
    type            (varying_string), intent(in   )                   :: fileNameCoolingFunction                , fileNameChemicalState
    integer                         , intent(in   )                   :: versionFileFormat
    logical                         , intent(in   ), optional         :: includeContinuum
    integer                         , parameter                       :: versionFileFormatCurrent     = 1
    double precision                , parameter                       :: temperatureMinimumLogarithmic=+2.500d+0, temperatureMaximumLogarithmic=+9.000d0, &
         &                                                               temperatureStepLogarithmic   =+0.025d+0
    double precision                , parameter                       :: metallicityMinimumLogarithmic=-4.000d+0, metallicityZeroLogarithmic   =-9.990d2, &
         &                                                               metallicityStepLogarithmic   =+0.250d+0
    double precision                , parameter                       :: energyMinimum                =+1.000d-3, energyMaximum                =+1.000d2
    integer                         , parameter                       :: energyBinsPerDecade          =50
    double precision                , allocatable  , dimension(:    ) :: metallicitiesLogarithmic               , temperaturesLogarithmic               , &
         &                                                               energyContinuum
    double precision                , allocatable  , dimension(:,:  ) :: coolingFunction                        , densityElectron                       , &
         &                                                               densityHydrogenI                       , densityHydrogenII
    double precision                , allocatable  , dimension(:,:,:) :: powerEmittedFractionalCumulative
    logical                                                           :: computeCoolingFunctions                , computeChemicalStates
    type            (hdf5Object    )                                  :: outputFile                             , dataset
    integer                                                           :: fileFormatFile                         , metallicityCount                     , &
         &                                                               temperatureCount                       , cloudyScript                         , &
         &                                                               iMetallicity                           , inputFile                            , &
         &                                                               iTemperature                           , status                               , &
         &                                                               energyCount                            , iEnergy                              , &
         &                                                               i
    type            (varying_string)                                  :: cloudyPath                             , cloudyVersion                        , &
         &                                                               fileNameTempCooling                    , fileNameTempOverview                 , &
         &                                                               fileNameTempContinuum                  , fileNameCoolingFunctionLock          , &
         &                                                               fileNameChemicalStateLock
    character       (len=   8      )                                  :: label
    character       (len=1024      )                                  :: line
    double precision                                                  :: dummy                                  , abundanceHelium                      , &
         &                                                               energy                                 , intensity                            , &
         &                                                               powerTotal
    !![
    <optionalArgument name="includeContinuum" defaultsTo=".false."/>
    !!]
    
    ! Ensure the requested file format version is compatible.
    if (versionFileFormat /= versionFileFormatCurrent) call Error_Report(var_str("this interface supports file format version ")//versionFileFormatCurrent//" but version "//versionFileFormat//" was requested"//{introspection:location})
    ! Determine if we need to compute cooling functions.
    call Directory_Make(File_Path(fileNameCoolingFunction))
    call Directory_Make(File_Path(fileNameChemicalState  ))
    call Directory_Make(inputPath(pathTypeDataDynamic)//'aux')
    fileNameCoolingFunctionLock=inputPath(pathTypeDataDynamic)//'aux/cloudyCIE_coolingFunction_'//Hash_MD5(fileNameCoolingFunction)
    fileNameChemicalStateLock  =inputPath(pathTypeDataDynamic)//'aux/cloudyCIE_chemicalState_'  //Hash_MD5(fileNameChemicalState  )
    computeCoolingFunctions=.false.
    computeChemicalStates  =.false.
    do i=1,2
       call File_Lock(fileNameCoolingFunction,fileLockCoolingFunction,lockIsShared=i == 1 .or. .not.computeCoolingFunctions,fileNameLock=fileNameCoolingFunctionLock)
       call File_Lock(fileNameChemicalState  ,fileLockChemicalState  ,lockIsShared=i == 1 .or. .not.computeChemicalStates  ,fileNameLock=fileNameChemicalStateLock  )
       computeCoolingFunctions=.false.
       computeChemicalStates  =.false.
       if (File_Exists(fileNameCoolingFunction)) then
          !$ call hdf5Access%set     (                                             )
          call    outputFile%openFile(char(fileNameCoolingFunction),readOnly=.true.)
          if (outputFile%hasAttribute('fileFormat')) then
             call outputFile%readAttribute('fileFormat',fileFormatFile)
             if (fileFormatFile /= versionFileFormatCurrent) computeCoolingFunctions=.true.
          end if
          call    outputFile%close()
          !$ call hdf5Access%unset()
       else
          computeCoolingFunctions=.true.
       end if
       ! Determine if we need to compute chemical states.
       if (File_Exists(fileNameChemicalState)) then
          !$ call hdf5Access%set     (                                           )
          call    outputFile%openFile(char(fileNameChemicalState),readOnly=.true.)
          if (outputFile%hasAttribute('fileFormat')) then
             call outputFile%readAttribute('fileFormat',fileFormatFile)
             if (fileFormatFile /= versionFileFormatCurrent) computeChemicalStates=.true.
          end if
          call    outputFile%close()
          !$ call hdf5Access%unset()
       else
          computeChemicalStates=.true.
       end if
       if ((computeCoolingFunctions .or. computeChemicalStates) .and. i == 1) then
          call File_Unlock(fileLockCoolingFunction,sync=.false.)
          call File_Unlock(fileLockChemicalState  ,sync=.false.)
          cycle
       end if
       ! Perform calculations if necessary.
       if (computeCoolingFunctions .or. computeChemicalStates) then
          ! Generate metallicity and temperature arrays.
          metallicityCount=int(     (+metallicityMaximumLogarithmic-metallicityMinimumLogarithmic)/     metallicityStepLogarithmic+1.5d0)
          temperatureCount=int(     (+temperatureMaximumLogarithmic-temperatureMinimumLogarithmic)/     temperatureStepLogarithmic+1.5d0)
          energyCount     =int(log10(+energyMaximum                /energyMinimum                )*dble(energyBinsPerDecade      )+1    )
          allocate(metallicitiesLogarithmic(metallicityCount+1                 ))
          allocate(temperaturesLogarithmic (                   temperatureCount))
          allocate(coolingFunction         (metallicityCount+1,temperatureCount))
          allocate(densityElectron         (metallicityCount+1,temperatureCount))
          allocate(densityHydrogenI        (metallicityCount+1,temperatureCount))
          allocate(densityHydrogenII       (metallicityCount+1,temperatureCount))
          metallicitiesLogarithmic(1                   )=metallicityZeroLogarithmic
          metallicitiesLogarithmic(2:metallicityCount+1)=Make_Range(metallicityMinimumLogarithmic,metallicityMaximumLogarithmic,metallicityCount,rangeTypeLinear)
          temperaturesLogarithmic                       =Make_Range(temperatureMinimumLogarithmic,temperatureMaximumLogarithmic,temperatureCount,rangeTypeLinear)
          if (includeContinuum_) then
             allocate(energyContinuum                 (                                    energyCount))
             allocate(powerEmittedFractionalCumulative(metallicityCount+1,temperatureCount,energyCount))
             energyContinuum                 =Make_Range(energyMinimum,energyMaximum,energyCount,rangeTypeLogarithmic)
             powerEmittedFractionalCumulative=0.0d0
          end if
          ! Initialize Cloudy.
          call Interface_Cloudy_Initialize(cloudyPath,cloudyVersion)
          ! Specify file names for temporary Cloudy data.
          fileNameTempCooling  ="cloudy_cooling.tmp"
          fileNameTempOverview ="cloudy_overview.tmp"
          fileNameTempContinuum="cloudy_continuum.tmp"
          ! Begin iterating over metallicities.
          call displayIndent("Computing cooling functions and chemical states using Cloudy (this may take a long time)....",verbosityLevelWorking)
          do iMetallicity=1,metallicityCount+1
             if (metallicitiesLogarithmic(iMetallicity) <= metallicityZeroLogarithmic) then
                write (label,'(a)'   ) '  -∞'
             else
                write (label,'(f8.3)') metallicitiesLogarithmic(iMetallicity)
             end if
             call displayMessage("Computing for log₁₀(Z/Z☉)="//trim(label),verbosityLevelWorking)
             call displayCounter(int(100.0d0*dble(iMetallicity-1)/dble(metallicityCount+1)),iMetallicity==0,verbosityLevelWorking)
             ! Generate an input file for Cloudy.
             open(newUnit=cloudyScript,file=char(cloudyPath//'/source/input.in'),status='unknown',form='formatted')
             write (cloudyScript,'(a)') 'print off'
             write (cloudyScript,'(a)') 'background, z=0'            ! Use a very low level incident continuum.
             write (cloudyScript,'(a)') 'cosmic rays background'     ! Include cosmic ray background ionization rate.
             write (cloudyScript,'(a)') 'stop zone 1'                ! Stop after a single zone.
             write (cloudyScript,'(a)') 'no photoionization'         ! Set no photoionization.
             write (cloudyScript,'(a)') 'hden 0.0'
             if (metallicitiesLogarithmic(iMetallicity) <= metallicityZeroLogarithmic) then
                write (cloudyScript,'(a)') 'abundances primordial'
             else
                write (cloudyScript,'(a,f6.3)') 'metals _log ',metallicitiesLogarithmic(iMetallicity)
                ! Assume a linear growth of helium abundance with metallicity.
                abundanceHelium=heliumToHydrogenAbundancePrimordial+(heliumToHydrogenAbundanceSolar-heliumToHydrogenAbundancePrimordial)*(10.0d0**metallicitiesLogarithmic(iMetallicity))
                write (cloudyScript,'(a,f6.3)') 'element abundance linear helium ',abundanceHelium
             end if
             write        (cloudyScript,'(a,f6.3,a)') 'constant temper ',temperatureMinimumLogarithmic,' vary'
             write        (cloudyScript,'(a,f6.3,a,f6.3,a,f6.3)') 'grid ',temperatureMinimumLogarithmic,' to ',temperatureMaximumLogarithmic,' step ',temperatureStepLogarithmic
             write        (cloudyScript,'(a)') 'no molecules'
             write        (cloudyScript,'(a)') 'set trim -20'
             write        (cloudyScript,'(a)') 'save cooling "'                     //char(fileNameTempCooling  )//'"'
             write        (cloudyScript,'(a)') 'save overview "'                    //char(fileNameTempOverview )//'"'
             if (includeContinuum_) &
                  & write (cloudyScript,'(a)') 'save emitted continuum units _keV "'//char(fileNameTempContinuum)//'"'
             close(cloudyScript)
             call System_Command_Do("cd "//cloudyPath//"/source; ./cloudy.exe -r input",status);
             if (status /= 0) call Error_Report('Cloudy failed'//{introspection:location})
             ! Extract the cooling rate.
             open(newUnit=inputFile,file=char(cloudyPath//"/source/"//fileNameTempCooling),status='old')
             read (inputFile,*) ! Skip the header line.
             do iTemperature=1,temperatureCount
                read (inputFile,*) dummy,dummy,dummy,coolingFunction(iMetallicity,iTemperature)
                read (inputFile,*)
             end do
             close(inputFile)
             call File_Remove(cloudyPath//"/source/"//fileNameTempCooling)
             ! Extract the electron and hydrogen density.
             open(newUnit=inputFile,file=char(cloudyPath//"/source/"//fileNameTempOverview),status='old')
             read (inputFile,*) ! Skip the header line.
             do iTemperature=1,temperatureCount
                read (inputFile,*) dummy,dummy,dummy,dummy,densityElectron(iMetallicity,iTemperature),dummy,densityHydrogenI(iMetallicity,iTemperature),densityHydrogenII(iMetallicity,iTemperature)
                read (inputFile,*)
             end do
             close(inputFile)
             call File_Remove(cloudyPath//"/source/"//fileNameTempOverview)
             ! Extract the emitted continuum. The continuum is given in units of ν f_ν [ergs cm⁻² s⁻¹], and frequencies are
             ! spaced uniformly in their logarithm. Therefore, to find the total power emitted we can simply sum the column
             ! values. We sum these into bins of energy, while also accumulating the total power. After reading all lines we
             ! normalize by the total power and cumulate to get the fraction of power emitted up to a given energy.
             if (includeContinuum_) then
                open(newUnit=inputFile,file=char(cloudyPath//"/source/"//fileNameTempContinuum),status='old')
                read (inputFile,*,iostat=status) ! Skip the header line.
                iTemperature=1
                iEnergy     =1
                powerTotal  =0.0d0
                do while (status == 0)
                   read (inputFile,'(a)',iostat=status) line
                   if (status /= 0) exit
                   if (line(1:1) == "#") then
                      ! New iteration reached.
                      !! Normalize and cumulate the power
                      powerEmittedFractionalCumulative(iMetallicity,iTemperature,:)=+powerEmittedFractionalCumulative(iMetallicity,iTemperature,:) &
                           &                                                        /powerTotal
                      do iEnergy=2,energyCount
                         powerEmittedFractionalCumulative(iMetallicity,iTemperature,iEnergy)=+powerEmittedFractionalCumulative(iMetallicity,iTemperature,iEnergy  ) &
                              &                                                              +powerEmittedFractionalCumulative(iMetallicity,iTemperature,iEnergy-1)
                      end do
                      !! Move to the next temperature.
                      iTemperature=iTemperature+1
                      iEnergy     =             1
                      powerTotal  =             0.0d0
                   else
                      read (line,*) energy,dummy,intensity
                      do while (energy > energyContinuum(iEnergy))
                         iEnergy=iEnergy+1
                         if (iEnergy > energyCount) exit
                      end do
                      if (iEnergy < energyCount) &
                           powerEmittedFractionalCumulative(iMetallicity,iTemperature,iEnergy)=powerEmittedFractionalCumulative(iMetallicity,iTemperature,iEnergy)+intensity
                      powerTotal=powerTotal+intensity
                   end if
                end do
                close(inputFile)
                call File_Remove(cloudyPath//"/source/"//fileNameTempContinuum)
             end if
          end do
          call displayCounterClear(verbosityLevelWorking)
          ! Output cooling functions to an HDF5 file.
          if (computeCoolingFunctions) then
             !$ call hdf5Access%set           (                                                                                                                     )
             call    outputFile%openFile      (char(fileNameCoolingFunction)                                                                                        )
             ! Store data.
             call    outputFile%writeDataset  (metallicitiesLogarithmic                                  ,'metallicity'                     ,datasetReturned=dataset)
             call    dataset   %writeAttribute('fix'                                                     ,'extrapolateLow'                                          )
             call    dataset   %writeAttribute('fix'                                                     ,'extrapolateHigh'                                         )
             call    dataset   %writeAttribute('K'                                                       ,'units'                                                   )
             call    dataset   %writeAttribute(1.0d0                                                     ,'unitsInSI'                                               )
             call    dataset   %close         (                                                                                                                     )
             call    outputFile%writeDataset  (10.0d0**temperaturesLogarithmic                           ,'temperature'                     ,datasetReturned=dataset)
             call    dataset   %writeAttribute('extrapolate'                                             ,'extrapolateLow'                                          )
             call    dataset   %writeAttribute('extrapolate'                                             ,'extrapolateHigh'                                         )
             call    dataset   %close         (                                                                                                                     )
             call    outputFile%writeDataset  (coolingFunction                                           ,'coolingRate'                     ,datasetReturned=dataset)
             call    dataset   %close         (                                                                                                                     )
             if (includeContinuum_) then
                call outputFile%writeDataset  (energyContinuum                                           ,'energyContinuum'                 ,datasetReturned=dataset)
                call dataset   %writeAttribute('extrapolate'                                             ,'extrapolateLow'                                          )
                call dataset   %writeAttribute('extrapolate'                                             ,'extrapolateHigh'                                         )
                call dataset   %writeAttribute('keV'                                                     ,'units'                                                   )
                call dataset   %writeAttribute(kilo*electronVolt                                         ,'unitsInSI'                                               )
                call dataset   %close         (                                                                                                                     )
                call outputFile%writeDataset  (powerEmittedFractionalCumulative                          ,'powerEmittedFractionalCumulative',datasetReturned=dataset)
                call dataset   %close         (                                                                                                                     )
             end if
             ! Add attributes.
             call    outputFile%writeAttribute("CIE cooling functions computed by Cloudy "//cloudyVersion,'description'                                             )
             call    outputFile%writeAttribute(versionFileFormatCurrent                                  ,'fileFormat'                                              )
             call    outputFile%close         (                                                                                                                     )
             !$ call hdf5Access%unset         (                                                                                                                     )
          end if
          ! Output chemical states to an HDF5 file.
          if (computeChemicalStates) then
             !$ call hdf5Access%set           (                                                                                                    )
             call    outputFile%openFile      (char(fileNameChemicalState)                                                                         )
             ! Store data.
             call    outputFile%writeDataset  (metallicitiesLogarithmic                                  ,'metallicity'    ,datasetReturned=dataset)
             call    dataset   %writeAttribute('fix'                                                     ,'extrapolateLow'                         )
             call    dataset   %writeAttribute('fix'                                                     ,'extrapolateHigh'                        )
             call    dataset   %writeAttribute('K'                                                       ,'units'                                  )
             call    dataset   %writeAttribute(1.0d0                                                     ,'unitsInSI'                              )
             call    dataset   %close         (                                                                                                    )
             call    outputFile%writeDataset  (10.0d0**temperaturesLogarithmic                           ,'temperature'    ,datasetReturned=dataset)
             call    dataset   %writeAttribute('extrapolate'                                             ,'extrapolateLow'                         )
             call    dataset   %writeAttribute('extrapolate'                                             ,'extrapolateHigh'                        )
             call    dataset   %close         (                                                                                                    )
             call    outputFile%writeDataset  (densityElectron                                           ,'electronDensity',datasetReturned=dataset)
             call    dataset   %close         (                                                                                                    )
             call    outputFile%writeDataset  (densityHydrogenI                                          ,'hiDensity'      ,datasetReturned=dataset)
             call    dataset   %close         (                                                                                                    )
             call    outputFile%writeDataset  (densityHydrogenII                                         ,'hiiDensity'     ,datasetReturned=dataset)
             call    dataset   %close         (                                                                                                    )
             ! Add attributes.
             call    outputFile%writeAttribute("CIE ionization states computed by Cloudy "//cloudyVersion,'description'                            )
             call    outputFile%writeAttribute(versionFileFormatCurrent                                  ,'fileFormat'                             )
             call    outputFile%close         (                                                                                                    )
             !$ call hdf5Access%unset         (                                                                                                    )
          end if
          ! Write message.
          call displayUnindent("...done",verbosityLevelWorking)
       end if
       call File_Unlock(fileLockChemicalState  )
       call File_Unlock(fileLockCoolingFunction)
    end do
    return
  end subroutine Interface_Cloudy_CIE_Tabulate

end module Interfaces_Cloudy_CIE

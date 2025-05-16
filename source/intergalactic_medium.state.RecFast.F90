!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  An implementation of the intergalactic medium state class in which state is computed using {\normalfont \scshape RecFast}.
  !!}

  use :: File_Utilities, only : lockDescriptor

  !![
  <intergalacticMediumState name="intergalacticMediumStateRecFast">
   <description>
    An intergalactic medium state class which computes the state of the intergalactic medium using the
    \href{https://www.astro.ubc.ca/people/scott/recfast.html}{{\normalfont \scshape RecFast}} code
    \cite{seager_how_2000,wong_how_2008}. The {\normalfont \scshape RecFast} code will be downloaded and run to compute the
    intergalactic medium state as needed, which will then be stored for future use.
   </description>
  </intergalacticMediumState>
  !!]
  type, extends(intergalacticMediumStateFile) :: intergalacticMediumStateRecFast
     !!{
     An \gls{igm} state class which computes state using {\normalfont \scshape RecFast}.
     !!}
     private
     type (lockDescriptor) :: fileLock
  end type intergalacticMediumStateRecFast

  interface intergalacticMediumStateRecFast
     !!{
     Constructors for the \refClass{intergalacticMediumStateRecFast} intergalactic medium state class.
     !!}
     module procedure recFastConstructorParameters
     module procedure recFastConstructorInternal
  end interface intergalacticMediumStateRecFast

contains

  function recFastConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the \refClass{intergalacticMediumStateRecFast} \gls{igm} state class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (intergalacticMediumStateRecFast)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass       ), pointer       :: cosmologyParameters_

    ! Check and read parameters.
    !![
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=intergalacticMediumStateRecFast(cosmologyFunctions_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function recFastConstructorParameters

  function recFastConstructorInternal(cosmologyFunctions_,cosmologyParameters_) result(self)
    !!{
    Constructor for the \refClass{intergalacticMediumStateRecFast} \gls{igm} state class.
    !!}
    use :: Cosmology_Parameters            , only : cosmologyParametersClass    , hubbleUnitsStandard
    use :: Dates_and_Times                 , only : Formatted_Date_and_Time
    use :: File_Utilities                  , only : Count_Lines_in_File         , Directory_Make     , File_Exists, File_Lock, &
          &                                         File_Unlock                 , File_Name_Temporary, File_Remove
    use :: Input_Paths                     , only : inputPath                   , pathTypeDataDynamic
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Interfaces_RecFast              , only : Interface_RecFast_Initialize
    use :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial
    use :: System_Command                  , only : System_Command_Do
    implicit none
    type            (intergalacticMediumStateRecFast)                              :: self
    class           (cosmologyFunctionsClass        ), intent(in   ), target       :: cosmologyFunctions_
    class           (cosmologyParametersClass       ), intent(in   ), target       :: cosmologyParameters_
    double precision                                 , allocatable  , dimension(:) :: redshift            , electronFraction , &
         &                                                                            hIonizedFraction    , heIonizedFraction, &
         &                                                                            matterTemperature
    character       (len=32                         )                              :: parameterLabel
    type            (varying_string                 )                              :: parameterFile       , recFastFile      , &
         &                                                                            recfastPath         , recfastVersion
    double precision                                                               :: omegaDarkMatter
    integer                                                                        :: fileFormatVersion   , i                , &
         &                                                                            countRedshift       , parametersUnit   , &
         &                                                                            recFastUnit
    type            (hdf5Object                     )                              :: outputFile          , dataset          , &
         &                                                                            provenance          , recFastProvenance
    logical                                                                        :: buildFile
    type            (lockDescriptor                 )                              :: fileLock
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_"/>
    !!]

    ! Compute dark matter density.
    omegaDarkMatter=self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon()
    ! Construct the file name.
    self%fileName=char(inputPath(pathTypeDataDynamic))//'intergalacticMedium/recFast'
    write (parameterLabel,'(f6.4)') self%cosmologyParameters_%OmegaMatter    (                   )
    self%fileName=self%fileName//'_OmegaMatter'    //trim(parameterLabel)
    write (parameterLabel,'(f6.4)') self%cosmologyParameters_%OmegaDarkEnergy(                   )
    self%fileName=self%fileName//'_OmegaDarkEnergy'//trim(parameterLabel)
    write (parameterLabel,'(f6.4)') self%cosmologyParameters_%OmegaBaryon    (                   )
    self%fileName=self%fileName//'_OmegaBaryon'    //trim(parameterLabel)
    write (parameterLabel,'(f4.1)') self%cosmologyParameters_%HubbleConstant (hubbleUnitsStandard)
    self%fileName=self%fileName//'_HubbleConstant' //trim(parameterLabel)
    write (parameterLabel,'(f6.4)') self%cosmologyParameters_%temperatureCMB (                   )
    self%fileName=self%fileName//'_temperatureCMB' //trim(parameterLabel)
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    self%fileName=self%fileName//'_YHe'            //trim(parameterLabel)
    self%fileName=self%fileName//'.hdf5'
    ! Create directory for output.
    call Directory_Make(inputPath(pathTypeDataDynamic)//'intergalacticMedium')
    ! Lock file.
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock(char(self%fileName),self%fileLock,lockIsShared=.false.)
    ! Check existence of file.
    buildFile=.false.
    if (File_Exists(char(self%fileName))) then
       ! Check file version number.
       !$ call hdf5Access%set()
       call outputFile%openFile     (char(self%fileName),overwrite=.false.          ,readOnly=.true.)
       call outputFile%readAttribute('fileFormat'       ,          fileFormatVersion                )
       call outputFile%close        (                                                               )
       !$ call hdf5Access%unset()
       buildFile=fileFormatVersion /= fileFormatVersionCurrent
    else
       buildFile=.true.
    end if
    ! Build file if necessary.
    if (buildFile) then
       ! Initialize RecFast.
       call Interface_RecFast_Initialize(recfastPath,recfastVersion)
       ! Build RecFast parameter file.
       parameterFile=File_Name_Temporary("recFastParameters")
       recFastFile  =parameterFile//".out"
       open(newUnit=parametersUnit,file=char(parameterFile),status='unknown',form='formatted')
       write (parametersUnit,'(a)')  char(recFastFile)
       write (parametersUnit,*    )  self%cosmologyParameters_%OmegaBaryon   (),omegaDarkMatter                           ,self%cosmologyParameters_%OmegaDarkEnergy()
       write (parametersUnit,*    )  self%cosmologyParameters_%HubbleConstant(),self%cosmologyParameters_%temperatureCMB(),heliumByMassPrimordial
       write (parametersUnit,*    )  1
       write (parametersUnit,*    )  6
       close(parametersUnit)
       ! Run RecFast.
       call File_Lock(char(recfastPath//"recfast.exe"),fileLock,lockIsShared=.false.)
       call System_Command_Do(recfastPath//"recfast.exe < "//parameterFile)
       call File_Unlock(fileLock)
       ! Parse the output file.
       countRedshift=Count_Lines_in_File(recFastFile)-1
       allocate(redshift         (countRedshift))
       allocate(electronFraction (countRedshift))
       allocate(hIonizedFraction (countRedshift))
       allocate(heIonizedFraction(countRedshift))
       allocate(matterTemperature(countRedshift))
       open(newUnit=recFastUnit,file=char(recFastFile),status='old',form='formatted')
       read (recFastUnit,*) ! Skip header line.
       do i=1,countRedshift
          read (recFastUnit,*) redshift(i),electronFraction(i),hIonizedFraction(i),heIonizedFraction(i),matterTemperature(i)
       end do
       close(recFastUnit)
       call File_Remove(char(recFastFile  ))
       call File_Remove(char(parameterFile))
       ! Create the output file.
       !$ call hdf5Access%set()
       call outputFile%openFile      (char(self%fileName),overwrite=.true.)
       call outputFile%writeDataset  (redshift           ,'redshift'         ,'Redshift'                                            )
       call outputFile%writeDataset  (electronFraction   ,'electronFraction' ,'Electron fraction'                                   )
       call outputFile%writeDataset  (hIonizedFraction   ,'hIonizedFraction' ,'Fraction of ionized hydrogen'                        )
       call outputFile%writeDataset  (heIonizedFraction  ,'heIonizedFraction','Fraction of ionized helium'                          )
       call outputFile%writeDataset  (matterTemperature  ,'matterTemperature','Temperature of matter'       ,datasetReturned=dataset)
       call dataset   %writeAttribute('Kelvin'           ,'units'                                                                   )
       call dataset   %writeAttribute(1.0d0              ,'unitsInSI'                                                               )
       call dataset   %close         (                                                                                              )
       ! Add description and provenance to output structure.
       call outputFile%writeAttribute('IGM ionization/thermal state computed using RecFast','description'         )
       call outputFile%writeAttribute(fileFormatVersionCurrent                             ,'fileFormat'          )
       call outputFile%writeAttribute(1                                                    ,'extrapolationAllowed')
       provenance=outputFile%openGroup('provenance')
       call provenance%writeAttribute(char(Formatted_Date_and_Time())                      ,'date'                )
       call provenance%writeAttribute('Galacticus via RecFast'                             ,'source'              )
       recFastProvenance=provenance%openGroup('recFast'                                                                                            )
       call recFastProvenance%writeAttribute(trim(recFastVersion)                                                                        ,'version')
       call recFastProvenance%writeAttribute('Includes modification of H recombination. Includes all modifications for HeI recombination','notes'  )
       call recFastProvenance%close         (                                                                                                      )
       call provenance       %close         (                                                                                                      )
       call outputFile       %close         (                                                                                                      )
       !$ call hdf5Access%unset()
    end if
    call File_Unlock(self%fileLock)
    return
  end function recFastConstructorInternal

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
  An implementation of the intergalactic medium state class in which state is read from file.
  !!}

  use :: Numerical_Interpolation, only : interpolator
  use :: Table_Labels           , only : enumerationExtrapolationTypeType

  ! Current file format version for intergalactic medium state files.
  integer, parameter :: fileFormatVersionCurrent=1

  !![
  <intergalacticMediumState name="intergalacticMediumStateFile">
   <description>
    An intergalactic medium state class which reads the state of the intergalactic medium from a file and interpolates in the
    tabulated results. The HDF5 file containing the table should have the following form:
    \begin{verbatim}
    HDF5 "igmState.hdf5" {
    GROUP "/" {
       ATTRIBUTE "fileFormat" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SCALAR
          DATA {
          (0): 1
          }
       }
       GROUP "Parameters" {
          ATTRIBUTE "HubbleConstant" {
             DATATYPE  H5T_IEEE_F64LE
             DATASPACE  SCALAR
             DATA {
             (0): 67.8
             }
          }
          ATTRIBUTE "OmegaBaryon" {
             DATATYPE  H5T_IEEE_F64LE
             DATASPACE  SCALAR
             DATA {
             (0): 0.0484
             }
          }
          ATTRIBUTE "OmegaDarkEnergy" {
             DATATYPE  H5T_IEEE_F64LE
             DATASPACE  SCALAR
             DATA {
             (0): 0.692
             }
          }
          ATTRIBUTE "OmegaMatter" {
             DATATYPE  H5T_IEEE_F64LE
             DATASPACE  SCALAR
             DATA {
             (0): 0.308
             }
          }
          ATTRIBUTE "Y_He" {
             DATATYPE  H5T_IEEE_F64LE
             DATASPACE  SCALAR
             DATA {
             (0): 0.22
             }
          }
          ATTRIBUTE "temperatureCMB" {
             DATATYPE  H5T_IEEE_F64LE
             DATASPACE  SCALAR
             DATA {
             (0): 2.725
             }
          }
       }
       DATASET "electronFraction" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 10000 ) / ( 10000 ) }
       }
       DATASET "hIonizedFraction" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 10000 ) / ( 10000 ) }
       }
       DATASET "heIonizedFraction" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 10000 ) / ( 10000 ) }
       }
       DATASET "matterTemperature" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 10000 ) / ( 10000 ) }
          ATTRIBUTE "units" {
             DATATYPE  H5T_STRING {
                STRSIZE 6;
                STRPAD H5T_STR_NULLTERM;
                CSET H5T_CSET_ASCII;
                CTYPE H5T_C_S1;
             }
             DATASPACE  SCALAR
             DATA {
             (0): "Kelvin"
             }
          }
          ATTRIBUTE "unitsInSI" {
             DATATYPE  H5T_IEEE_F64LE
             DATASPACE  SCALAR
             DATA {
             (0): 1
             }
          }
       }
       DATASET "redshift" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 10000 ) / ( 10000 ) }
       }
    }
    }
    \end{verbatim}
     The {\normalfont \ttfamily electronFraction}, {\normalfont \ttfamily hIonizedFraction}, {\normalfont \ttfamily
    heIonizedFraction}, and {\normalfont \ttfamily matterTemperature} datasets contain the relevant quantity for each redshift
    in the {\normalfont \ttfamily redshift} dataset.
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </intergalacticMediumState>
  !!]
  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateFile
     !!{
     An \gls{igm} state class which reads state from file.
     !!}
     private
     ! Name of the file from which to read intergalactic medium state data.
     type            (varying_string                  )                            :: fileName
     ! Flag indicating whether or not data has been read.
     logical                                                                       :: dataRead                           =.false.
     ! Data tables.
     integer                                                                       :: redshiftCount
     double precision                                  , allocatable, dimension(:) :: electronFractionTable                      , temperatureTable                 , &
          &                                                                           timeTable                                  , ionizedHydrogenFractionTable     , &
          &                                                                           ionizedHeliumFractionTable
     type            (interpolator                    )                            :: interpolatorElectronFraction               , interpolatorTemperature          , &
          &                                                                           interpolatorIonizedHydrogenFraction        , interpolatorIonizedHeliumFraction
     type            (enumerationExtrapolationTypeType)                            :: extrapolationType
   contains
     final     ::                                fileDestructor
     procedure :: electronFraction            => fileElectronFraction
     procedure :: neutralHydrogenFraction     => fileNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => fileNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => fileSinglyIonizedHeliumFraction
     procedure :: temperature                 => fileTemperature
  end type intergalacticMediumStateFile

  interface intergalacticMediumStateFile
     !!{
     Constructors for the file intergalactic medium state class.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface intergalacticMediumStateFile

contains

  function fileConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the file \gls{igm} state class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (intergalacticMediumStateFile)                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass    ), pointer       :: cosmologyParameters_
    type (varying_string              )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <variable>fileName</variable>
      <description>The name of the file from which to read intergalactic medium state data.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=intergalacticMediumStateFile(fileName,cosmologyFunctions_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName,cosmologyFunctions_,cosmologyParameters_) result(self)
    !!{
    Constructor for the file \gls{igm} state class.
    !!}
    use :: File_Utilities, only : File_Name_Expand
    implicit none
    type (intergalacticMediumStateFile)                        :: self
    type (varying_string              ), intent(in   )         :: fileName
    class(cosmologyFunctionsClass     ), intent(inout), target :: cosmologyFunctions_
    class(cosmologyParametersClass    ), intent(inout), target :: cosmologyParameters_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_"/>
    !!]

    self%fileName=File_Name_Expand(char(fileName))
    return
  end function fileConstructorInternal

  subroutine fileDestructor(self)
    !!{
    Destructor for the file \gls{igm} state class.
    !!}
    implicit none
    type(intergalacticMediumStateFile), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine fileDestructor

  subroutine fileReadData(self)
    !!{
    Read in data describing the state of the intergalactic medium.
    !!}
    use :: File_Utilities, only : File_Exists
    use :: Error         , only : Error_Report
    use :: HDF5_Access   , only : hdf5Access
    use :: IO_HDF5       , only : hdf5Object
    use :: Table_Labels  , only : extrapolationTypeAbort, extrapolationTypeExtrapolate
    implicit none
    class  (intergalacticMediumStateFile), intent(inout) :: self
    integer                                              :: fileFormatVersion   , iRedshift, &
         &                                                  extrapolationAllowed

    ! Check if data has yet to be read.
    if (self%dataRead) return
    block
      type(hdf5Object) :: file
      
      if (.not.File_Exists(self%fileName)) call Error_Report('Unable to find intergalactic medium state file "' //char(self%fileName)//'"'//{introspection:location})
      !$ call hdf5Access%set()
      ! Open the file.
      file=hdf5Object(self%fileName,readOnly=.true.)
      ! Check the file format version of the file.
      call file%readAttribute('fileFormat',fileFormatVersion)
      if (fileFormatVersion /= fileFormatVersionCurrent) call Error_Report('file format version is out of date'//{introspection:location})
      ! Check if extrapolation is allowed.
      self%extrapolationType=extrapolationTypeAbort
      if (file%hasAttribute('extrapolationAllowed')) then
         call file%readAttribute('extrapolationAllowed',extrapolationAllowed)
         if (extrapolationAllowed /= 0) self%extrapolationType=extrapolationTypeExtrapolate
      end if
      call file%readAttribute('fileFormat',fileFormatVersion)
      ! Read the data.
      call file%readDataset('redshift'         ,self%timeTable                   )
      call file%readDataset('electronFraction' ,self%electronFractionTable       )
      call file%readDataset('hIonizedFraction' ,self%ionizedHydrogenFractionTable)
      call file%readDataset('heIonizedFraction',self%ionizedHeliumFractionTable  )
      call file%readDataset('matterTemperature',self%temperatureTable            )
      !$ call hdf5Access%unset()
      self%redshiftCount=size(self%timeTable)
      ! Convert redshifts to times.
      do iRedshift=1,self%redshiftCount
         self%timeTable(iRedshift)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%timeTable(iRedshift)))
      end do
      ! Build interpolators.
      self%interpolatorElectronFraction       =interpolator(self%timeTable,self%electronFractionTable       ,extrapolationType=self%extrapolationType)
      self%interpolatorTemperature            =interpolator(self%timeTable,self%temperatureTable            ,extrapolationType=self%extrapolationType)
      self%interpolatorIonizedHydrogenFraction=interpolator(self%timeTable,self%ionizedHydrogenFractionTable,extrapolationType=self%extrapolationType)
      self%interpolatorIonizedHeliumFraction  =interpolator(self%timeTable,self%ionizedHeliumFractionTable  ,extrapolationType=self%extrapolationType)
      ! Flag that data has now been read.
      self%dataRead=.true.
    end block
    return
  end subroutine fileReadData

  double precision function fileTemperature(self,time)
    !!{
    Return the temperature of the intergalactic medium at the specified {\normalfont \ttfamily time} by interpolating in tabulated data.
    !!}
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the temperature.
    fileTemperature=max(0.0d0,self%interpolatorTemperature%interpolate(time))
    return
  end function fileTemperature

  double precision function fileElectronFraction(self,time)
    !!{
    Return the electron fraction in the intergalactic medium at the specified {\normalfont \ttfamily time} by interpolating in tabulated data,
    !!}
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileElectronFraction=min(1.0d0,max(0.0d0,self%interpolatorElectronFraction%interpolate(time)))
    return
  end function fileElectronFraction

  double precision function fileNeutralHydrogenFraction(self,time)
    !!{
    Return the neutral hydrogen fraction in the intergalactic medium at the specified {\normalfont \ttfamily time} by interpolating in tabulated data,
    !!}
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileNeutralHydrogenFraction=+1.0d0                                                                            &
         &                      -min(1.0d0,max(0.0d0,self%interpolatorIonizedHydrogenFraction%interpolate(time)))
    return
  end function fileNeutralHydrogenFraction

  double precision function fileNeutralHeliumFraction(self,time)
    !!{
    Return the neutral helium fraction in the intergalactic medium at the specified {\normalfont \ttfamily time} by interpolating in tabulated data,
    !!}
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileNeutralHeliumFraction=+1.0d0                                                                          &
         &                    -min(1.0d0,max(0.0d0,self%interpolatorIonizedHeliumFraction%interpolate(time)))
    return
  end function fileNeutralHeliumFraction

  double precision function fileSinglyIonizedHeliumFraction(self,time)
    !!{
    Return the neutral helium fraction in the intergalactic medium at the specified {\normalfont \ttfamily time} by interpolating in tabulated data,
    !!}
    implicit none
    class           (intergalacticMediumStateFile), intent(inout) :: self
    double precision                              , intent(in   ) :: time

    ! Ensure that data has been read.
    call fileReadData(self)
    ! Interpolate in the tables to get the electron fraction.
    fileSinglyIonizedHeliumFraction=min(1.0d0,max(0.0d0,self%interpolatorIonizedHeliumFraction%interpolate(time)))
    return
  end function fileSinglyIonizedHeliumFraction

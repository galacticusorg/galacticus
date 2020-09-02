!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a file-based transfer function class.

  use :: Cosmology_Functions , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass
  use :: Tables              , only : table1DGeneric

  !# <transferFunction name="transferFunctionFile">
  !#  <description>
  !# Provides a transfer function from a tabulation given in an HDF5 file with the following structure:
  !# \begin{verbatim}
  !# HDF5 "transferFunction.hdf5" {
  !# GROUP "/" {
  !#    ATTRIBUTE "description" {
  !#       DATATYPE  H5T_STRING {
  !#          STRSIZE 71;
  !#          STRPAD H5T_STR_NULLTERM;
  !#          CSET H5T_CSET_ASCII;
  !#          CTYPE H5T_C_S1;
  !#       }
  !#       DATASPACE  SCALAR
  !#    }
  !#    ATTRIBUTE "fileFormat" {
  !#       DATATYPE  H5T_STD_I32LE
  !#       DATASPACE  SCALAR
  !#    }
  !#    ATTRIBUTE "redshift" {
  !#       DATATYPE  H5T_STD_I32LE
  !#       DATASPACE  SCALAR
  !#    }
  !#    GROUP "extrapolation" {
  !#       GROUP "wavenumber" {
  !#          ATTRIBUTE "high" {
  !#             DATATYPE  H5T_STRING {
  !#                STRSIZE 11;
  !#                STRPAD H5T_STR_NULLTERM;
  !#                CSET H5T_CSET_ASCII;
  !#                CTYPE H5T_C_S1;
  !#             }
  !#             DATASPACE  SCALAR
  !#          }
  !#          ATTRIBUTE "low" {
  !#             DATATYPE  H5T_STRING {
  !#                STRSIZE 3;
  !#                STRPAD H5T_STR_NULLTERM;
  !#                CSET H5T_CSET_ASCII;
  !#                CTYPE H5T_C_S1;
  !#             }
  !#             DATASPACE  SCALAR
  !#          }
  !#       }
  !#    }
  !#    GROUP "parameters" {
  !#       ATTRIBUTE "HubbleConstant" {
  !#          DATATYPE  H5T_STRING {
  !#             STRSIZE 4;
  !#             STRPAD H5T_STR_NULLTERM;
  !#             CSET H5T_CSET_ASCII;
  !#             CTYPE H5T_C_S1;
  !#          }
  !#          DATASPACE  SCALAR
  !#       }
  !#       ATTRIBUTE "OmegaBaryon" {
  !#          DATATYPE  H5T_STRING {
  !#             STRSIZE 6;
  !#             STRPAD H5T_STR_NULLTERM;
  !#             CSET H5T_CSET_ASCII;
  !#             CTYPE H5T_C_S1;
  !#          }
  !#          DATASPACE  SCALAR
  !#       }
  !#       ATTRIBUTE "OmegaDarkEnergy" {
  !#          DATATYPE  H5T_STRING {
  !#             STRSIZE 5;
  !#             STRPAD H5T_STR_NULLTERM;
  !#             CSET H5T_CSET_ASCII;
  !#             CTYPE H5T_C_S1;
  !#          }
  !#          DATASPACE  SCALAR
  !#       }
  !#       ATTRIBUTE "OmegaMatter" {
  !#          DATATYPE  H5T_STRING {
  !#             STRSIZE 5;
  !#             STRPAD H5T_STR_NULLTERM;
  !#             CSET H5T_CSET_ASCII;
  !#             CTYPE H5T_C_S1;
  !#          }
  !#          DATASPACE  SCALAR
  !#       }
  !#    }
  !#    DATASET "transferFunction" {
  !#       DATATYPE  H5T_IEEE_F64LE
  !#       DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ) }
  !#    }
  !#    DATASET "wavenumber" {
  !#       DATATYPE  H5T_IEEE_F64LE
  !#       DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ) }
  !#    }
  !# }
  !# }
  !# \end{verbatim}
  !# </description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionFile
     !% A transfer function class which interpolates a transfer function given in a file.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     type            (varying_string          )          :: fileName
     type            (table1DGeneric          )          :: transfer
     double precision                                    :: time                          , redshift
   contains
     !@ <objectMethods>
     !@   <object>transferFunctionFile</object>
     !@   <objectMethod>
     !@     <method>readFile</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless char(len=*)\textgreater} fileName\argin</arguments>
     !@     <description>Read the named transfer function file.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                          fileDestructor
     procedure :: readFile              => fileReadFile
     procedure :: value                 => fileValue
     procedure :: logarithmicDerivative => fileLogarithmicDerivative
     procedure :: halfModeMass          => fileHalfModeMass
     procedure :: epochTime             => fileEpochTime
  end type transferFunctionFile

  interface transferFunctionFile
     !% Constructors for the file transfer function class.
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface transferFunctionFile

  ! Current file format version for transfer function files. Note that this file format matches that used by the CAMB interface
  ! module.
  integer, parameter :: fileFormatVersionCurrent=2

contains

  function fileConstructorParameters(parameters) result(self)
    !% Constructor for the file transfer function class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (transferFunctionFile    )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    type            (varying_string          )                :: fileName
    double precision                                          :: redshift

    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of the file from which to read a tabulated transfer function.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>redshift</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The redshift of the transfer function to read.</description>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    self=transferFunctionFile(char(fileName),redshift,cosmologyParameters_,cosmologyFunctions_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="cosmologyFunctions_" />
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName,redshift,cosmologyParameters_,cosmologyFunctions_) result(self)
    !% Internal constructor for the file transfer function class.
    implicit none
    type            (transferFunctionFile    )                        :: self
    character       (len=*                   ), intent(in   )         :: fileName
    double precision                          , intent(in   )         :: redshift
    class           (cosmologyParametersClass), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="fileName, redshift, *cosmologyParameters_, *cosmologyFunctions_"/>

    self%time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    call self%readFile(fileName)
    return
  end function fileConstructorInternal

  subroutine fileReadFile(self,fileName)
    !% Internal constructor for the file transfer function class.
    use :: Cosmology_Parameters   , only : cosmologyParametersSimple
    use :: File_Utilities         , only : File_Name_Expand
    use :: Galacticus_Display     , only : Galacticus_Display_Message
    use :: Galacticus_Error       , only : Galacticus_Error_Report
    use :: IO_HDF5                , only : hdf5Access                        , hdf5Object
    use :: Numerical_Comparison   , only : Values_Differ
    use :: Numerical_Interpolation, only : GSL_Interp_cSpline
    use :: Table_Labels           , only : enumerationExtrapolationTypeEncode
    implicit none
    class           (transferFunctionFile    ), intent(inout)             :: self
    character       (len=*                   ), intent(in   )             :: fileName
    double precision                          , allocatable, dimension(:) :: transfer                       , wavenumber               , &
         &                                                                   transferLogarithmic            , wavenumberLogarithmic
    class           (cosmologyParametersClass), pointer                   :: cosmologyParametersFile
    double precision                          , parameter                 :: toleranceUniformity     =1.0d-6
    double precision                                                      :: HubbleConstant                 , OmegaBaryon              , &
         &                                                                   OmegaMatter                    , OmegaDarkEnergy          , &
         &                                                                   temperatureCMB
    integer                                                               :: extrapolateWavenumberLow       , extrapolateWavenumberHigh, &
         &                                                                   versionNumber
    character       (len=32                  )                            :: datasetName
    type            (varying_string          )                            :: limitTypeVar
    type            (hdf5Object              )                            :: fileObject                     , parametersObject         , &
         &                                                                   extrapolationObject            , wavenumberObject         , &
         &                                                                   darkMatterGroup

    ! Open and read the HDF5 data file.
    call hdf5Access%set()
    call fileObject%openFile(char(File_Name_Expand(fileName)),readOnly=.true.)
    ! Check that the file has the correct format version number.
    call fileObject%readAttribute('fileFormat',versionNumber,allowPseudoScalar=.true.)
    if (versionNumber /= fileFormatVersionCurrent) call Galacticus_Error_Report('file has the incorrect version number'//{introspection:location})
    ! Check that parameters match if any are present.
    parametersObject=fileObject%openGroup('parameters')
    allocate(cosmologyParametersSimple :: cosmologyParametersFile)
    select type (cosmologyParametersFile)
    type is (cosmologyParametersSimple)
       call parametersObject%readAttribute('OmegaMatter'    ,OmegaMatter    )
       call parametersObject%readAttribute('OmegaDarkEnergy',OmegaDarkEnergy)
       call parametersObject%readAttribute('OmegaBaryon'    ,OmegaBaryon    )
       call parametersObject%readAttribute('HubbleConstant' ,HubbleConstant )
       call parametersObject%readAttribute('temperatureCMB' ,temperatureCMB )
       cosmologyParametersFile=cosmologyParametersSimple(OmegaMatter,OmegaBaryon,OmegaDarkEnergy,temperatureCMB,HubbleConstant)
       if (Values_Differ(cosmologyParametersFile%OmegaBaryon    (),self%cosmologyParameters_%OmegaBaryon    (),absTol=1.0d-3)) &
            & call Galacticus_Display_Message('OmegaBaryon from transfer function file does not match internal value'    )
       if (Values_Differ(cosmologyParametersFile%OmegaMatter    (),self%cosmologyParameters_%OmegaMatter    (),absTol=1.0d-3)) &
            & call Galacticus_Display_Message('OmegaMatter from transfer function file does not match internal value'    )
       if (Values_Differ(cosmologyParametersFile%OmegaDarkEnergy(),self%cosmologyParameters_%OmegaDarkEnergy(),absTol=1.0d-3)) &
            & call Galacticus_Display_Message('OmegaDarkEnergy from transfer function file does not match internal value')
       if (Values_Differ(cosmologyParametersFile%HubbleConstant (),self%cosmologyParameters_%HubbleConstant (),relTol=1.0d-3)) &
            & call Galacticus_Display_Message('HubbleConstant from transfer function file does not match internal value' )
       if (Values_Differ(cosmologyParametersFile%temperatureCMB (),self%cosmologyParameters_%temperatureCMB (),relTol=1.0d-3)) &
            & call Galacticus_Display_Message('temperatureCMB from transfer function file does not match internal value' )
    end select
    deallocate(cosmologyParametersFile)
    call parametersObject%close()
    ! Get extrapolation methods.
    extrapolationObject=fileObject         %openGroup('extrapolation')
    wavenumberObject   =extrapolationObject%openGroup('wavenumber'   )
    call wavenumberObject%readAttribute('low' ,limitTypeVar)
    extrapolateWavenumberLow =enumerationExtrapolationTypeEncode(char(limitTypeVar),includesPrefix=.false.)
    call wavenumberObject%readAttribute('high',limitTypeVar)
    extrapolateWavenumberHigh=enumerationExtrapolationTypeEncode(char(limitTypeVar),includesPrefix=.false.)
    call wavenumberObject   %close()
    call extrapolationObject%close()
    ! Read the transfer function from file.
    darkMatterGroup=fileObject%openGroup('darkMatter')
    call fileObject     %readDataset('wavenumber'                                   ,wavenumber)
    write (datasetName,'(f9.4)') self%redshift
    call darkMatterGroup%readDataset('transferFunctionZ'//trim(adjustl(datasetName)),transfer  )
    call darkMatterGroup%close      (                                                          )
    ! Close the file.
    call fileObject%close()
    call hdf5Access%unset()
    ! Construct the tabulated transfer function.
    call self%transfer%destroy()
    wavenumberLogarithmic=log(wavenumber)
    transferLogarithmic  =log(transfer  )
    ! Create the table.
    call self%transfer%create  (                                                    &
         &                      wavenumberLogarithmic                             , &
         &                      extrapolationType=[                                 &
         &                                         extrapolateWavenumberLow       , &
         &                                         extrapolateWavenumberHigh        &
         &                                        ]                               , &
         &                      interpolationType= GSL_Interp_cSpline               &
         &                     )
    call self%transfer%populate(                                                    &
         &                      transferLogarithmic                                 &
         &                     )
    return
  end subroutine fileReadFile

  subroutine fileDestructor(self)
    !% Destructor for the file transfer function class.
    implicit none
    type(transferFunctionFile), intent(inout) :: self

    call self%transfer%destroy()
    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    return
  end subroutine fileDestructor

  double precision function fileValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionFile), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber

    fileValue=exp(self%transfer%interpolate(log(wavenumber)))
    return
  end function fileValue

  double precision function fileLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionFile), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber

    fileLogarithmicDerivative=+self%transfer%interpolateGradient(log(wavenumber))
    return
  end function fileLogarithmicDerivative

  double precision function fileHalfModeMass(self,status)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is
    !% suppressed by a factor of two relative to a \gls{cdm} transfer function. Not supported in
    !% this implementation.
    use :: Galacticus_Error, only : Galacticus_Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionFile), intent(inout)           :: self
    integer                      , intent(  out), optional :: status
    !$GLC attributes unused :: self

    fileHalfModeMass=0.0d0
    if (present(status)) then
       status=errorStatusFail
    else
       call Galacticus_Error_Report('not supported by this implementation'//{introspection:location})
    end if
    return
  end function fileHalfModeMass

  double precision function fileEpochTime(self)
    !% Return the cosmic time at the epoch at which this transfer function is defined.
    implicit none
    class(transferFunctionFile), intent(inout) :: self

    fileEpochTime=self%time
    return
  end function fileEpochTime

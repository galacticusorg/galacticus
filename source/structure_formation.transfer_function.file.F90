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

  !!{
  Implements a file-based transfer function class.
  !!}

  use :: Cosmology_Functions , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass
  use :: Tables              , only : table1DGeneric

  !![
  <transferFunction name="transferFunctionFile">
   <description>
  Provides a transfer function from a tabulation given in an HDF5 file with the following structure:
  \begin{verbatim}
  HDF5 "transferFunction.hdf5" {
  GROUP "/" {
     ATTRIBUTE "description" {
        DATATYPE  H5T_STRING {
           STRSIZE 71;
           STRPAD H5T_STR_NULLTERM;
           CSET H5T_CSET_ASCII;
           CTYPE H5T_C_S1;
        }
        DATASPACE  SCALAR
     }
     ATTRIBUTE "fileFormat" {
        DATATYPE  H5T_STD_I32LE
        DATASPACE  SCALAR
     }
     ATTRIBUTE "redshift" {
        DATATYPE  H5T_STD_I32LE
        DATASPACE  SCALAR
     }
     GROUP "extrapolation" {
        GROUP "wavenumber" {
           ATTRIBUTE "high" {
              DATATYPE  H5T_STRING {
                 STRSIZE 11;
                 STRPAD H5T_STR_NULLTERM;
                 CSET H5T_CSET_ASCII;
                 CTYPE H5T_C_S1;
              }
              DATASPACE  SCALAR
           }
           ATTRIBUTE "low" {
              DATATYPE  H5T_STRING {
                 STRSIZE 3;
                 STRPAD H5T_STR_NULLTERM;
                 CSET H5T_CSET_ASCII;
                 CTYPE H5T_C_S1;
              }
              DATASPACE  SCALAR
           }
        }
     }
     GROUP "parameters" {
        ATTRIBUTE "HubbleConstant" {
           DATATYPE  H5T_STRING {
              STRSIZE 4;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
        ATTRIBUTE "OmegaBaryon" {
           DATATYPE  H5T_STRING {
              STRSIZE 6;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
        ATTRIBUTE "OmegaDarkEnergy" {
           DATATYPE  H5T_STRING {
              STRSIZE 5;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
        ATTRIBUTE "OmegaMatter" {
           DATATYPE  H5T_STRING {
              STRSIZE 5;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
     }
     GROUP "darkMatter" {
        DATASET "transferFunctionZ0.0000" {
           DATATYPE  H5T_IEEE_F64LE
           DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ) }
        }
     }
     DATASET "wavenumber" {
        DATATYPE  H5T_IEEE_F64LE
        DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ) }
     }
  }
  }
  \end{verbatim}

  If an optional \refClass{transferFunctionClass} object named {\normalfont \ttfamily transferFunctionReference} is supplied then
  that transfer function is multiplied by the tabulated transfer function. In this case half and quarter-mode masses relative to
  {\normalfont \ttfamily transferFunctionReference} are also computed.
   </description>
   <runTimeFileDependencies paths="fileName"/>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionFile
     !!{
     A transfer function class which interpolates a transfer function given in a file.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer                   :: cosmologyFunctions_                => null()
     class           (transferFunctionClass   ), pointer                   :: transferFunctionReference          => null()
     type            (varying_string          )                            :: fileName
     type            (table1DGeneric          )                            :: transfer
     logical                                                               :: massHalfModeAvailable                       , massQuarterModeAvailable, &
          &                                                                   transferFunctionReferenceAvailable          , acceptNegativeValues
     double precision                                                      :: time                                        , redshift                , &
          &                                                                   massHalfMode                                , massQuarterMode
     double precision                          , allocatable, dimension(:) :: wavenumbersLocalMinima_
   contains
     !![
     <methods>
       <method description="Read the named transfer function file." method="readFile" />
     </methods>
     !!]
     final     ::                           fileDestructor
     procedure :: readFile               => fileReadFile
     procedure :: value                  => fileValue
     procedure :: logarithmicDerivative  => fileLogarithmicDerivative
     procedure :: wavenumbersLocalMinima => fileWavenumbersLocalMinima
     procedure :: halfModeMass           => fileHalfModeMass
     procedure :: quarterModeMass        => fileQuarterModeMass
     procedure :: fractionModeMass       => fileFractionModeMass
     procedure :: epochTime              => fileEpochTime
  end type transferFunctionFile

  interface transferFunctionFile
     !!{
     Constructors for the file transfer function class.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface transferFunctionFile

  ! Current file format version for transfer function files. Note that this file format matches that used by the CAMB interface
  ! module.
  integer, parameter :: fileFormatVersionCurrent=2

contains

  function fileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the file transfer function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (transferFunctionFile    )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    class           (transferFunctionClass   ), pointer       :: transferFunctionReference
    type            (varying_string          )                :: fileName
    double precision                                          :: redshift
    logical                                                   :: acceptNegativeValues

    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file from which to read a tabulated transfer function.</description>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift of the transfer function to read.</description>
    </inputParameter>
    <inputParameter>
      <name>acceptNegativeValues</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, negative values in the transfer function are allowed (and the absolute value is taken prior to interpolation). Otherwise, negative values result in an error.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    if (parameters%isPresent('transferFunctionReference')) then
       !![
       <objectBuilder class="transferFunction" name="transferFunctionReference" source="parameters" parameterName="transferFunctionReference"/>
       !!]
    else
       transferFunctionReference => null()
    end if
    !![
    <conditionalCall>
      <call>self=transferFunctionFile(char(fileName),redshift,acceptNegativeValues,cosmologyParameters_,cosmologyFunctions_{conditions})</call>
      <argument name="transferFunctionReference" value="transferFunctionReference" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    if (parameters%isPresent('transferFunctionReference')) then
       !![
       <objectDestructor name="transferFunctionReference"/>
       !!]
    end if
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName,redshift,acceptNegativeValues,cosmologyParameters_,cosmologyFunctions_,transferFunctionReference) result(self)
    !!{
    Internal constructor for the file transfer function class.
    !!}
    implicit none
    type            (transferFunctionFile    )                                  :: self
    character       (len=*                   ), intent(in   )                   :: fileName
    double precision                          , intent(in   )                   :: redshift
    logical                                   , intent(in   )                   :: acceptNegativeValues
    class           (cosmologyParametersClass), intent(in   ), target           :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), intent(in   ), target           :: cosmologyFunctions_
    class           (transferFunctionClass   ), intent(in   ), target, optional :: transferFunctionReference
    integer                                                                     :: status
    !![
    <constructorAssign variables="fileName, redshift, acceptNegativeValues, *cosmologyParameters_, *cosmologyFunctions_, *transferFunctionReference"/>
    !!]

    self%time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    call self%readFile(fileName)
    ! Compute half and quarter-mode masses.
    self%transferFunctionReferenceAvailable=present(transferFunctionReference)
    self%massHalfModeAvailable             =.false.
    self%massQuarterModeAvailable          =.false.
    if (self%transferFunctionReferenceAvailable) then
       self%massHalfMode            =self%fractionModeMass(fraction=0.50d0,status=status)
       self%massHalfModeAvailable   =status == errorStatusSuccess
       self%massQuarterMode         =self%fractionModeMass(fraction=0.25d0,status=status)
       self%massQuarterModeAvailable=status == errorStatusSuccess
    end if
    return
  end function fileConstructorInternal

  subroutine fileReadFile(self,fileName)
    !!{
    Internal constructor for the file transfer function class.
    !!}
    use :: Cosmology_Parameters   , only : cosmologyParametersSimple
    use :: Display                , only : displayMessage                  , displayMagenta                    , displayReset, displayGreen, &
         &                                 displayYellow                   , displayBlue
    use :: File_Utilities         , only : File_Name_Expand
    use :: Error                  , only : Error_Report
    use :: HDF5_Access            , only : hdf5Access
    use :: IO_HDF5                , only : hdf5Object
    use :: Numerical_Comparison   , only : Values_Differ
    use :: Numerical_Interpolation, only : GSL_Interp_cSpline
    use :: Table_Labels           , only : enumerationExtrapolationTypeType, enumerationExtrapolationTypeEncode
    implicit none
    class           (transferFunctionFile            ), intent(inout)             :: self
    character       (len=*                           ), intent(in   )             :: fileName
    double precision                                  , allocatable, dimension(:) :: transfer                       , wavenumber               , &
         &                                                                           transferLogarithmic            , wavenumberLogarithmic
    class           (cosmologyParametersClass        ), pointer                   :: cosmologyParametersFile
    double precision                                  , parameter                 :: toleranceUniformity     =1.0d-6
    double precision                                                              :: HubbleConstant                 , OmegaBaryon              , &
         &                                                                           OmegaMatter                    , OmegaDarkEnergy          , &
         &                                                                           temperatureCMB
    type            (enumerationExtrapolationTypeType)                            :: extrapolateWavenumberLow       , extrapolateWavenumberHigh
    integer                                                                       :: versionNumber                  , i                        , &
         &                                                                           countLocalMinima
    character       (len=32                          )                            :: datasetName
    type            (varying_string                  )                            :: limitTypeVar
    type            (hdf5Object                      )                            :: fileObject

    ! Open and read the HDF5 data file.
    block
      type(hdf5Object) :: darkMatterGroup    , parametersObject, &
           &              extrapolationObject, wavenumberObject
      !$ call hdf5Access%set()
      fileObject=hdf5Object(char(File_Name_Expand(fileName)),readOnly=.true.)
      ! Check that the file has the correct format version number.
      call fileObject%readAttribute('fileFormat',versionNumber,allowPseudoScalar=.true.)
      if (versionNumber /= fileFormatVersionCurrent) call Error_Report('file has the incorrect version number'//{introspection:location})
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
              & call displayMessage(displayMagenta()//'WARNING: '//displayReset()//'OmegaBaryon from transfer function file does not match internal value'    )
         if (Values_Differ(cosmologyParametersFile%OmegaMatter    (),self%cosmologyParameters_%OmegaMatter    (),absTol=1.0d-3)) &
              & call displayMessage(displayMagenta()//'WARNING: '//displayReset()//'OmegaMatter from transfer function file does not match internal value'    )
         if (Values_Differ(cosmologyParametersFile%OmegaDarkEnergy(),self%cosmologyParameters_%OmegaDarkEnergy(),absTol=1.0d-3)) &
              & call displayMessage(displayMagenta()//'WARNING: '//displayReset()//'OmegaDarkEnergy from transfer function file does not match internal value')
         if (Values_Differ(cosmologyParametersFile%HubbleConstant (),self%cosmologyParameters_%HubbleConstant (),relTol=1.0d-3)) &
              & call displayMessage(displayMagenta()//'WARNING: '//displayReset()//'HubbleConstant from transfer function file does not match internal value' )
         if (Values_Differ(cosmologyParametersFile%temperatureCMB (),self%cosmologyParameters_%temperatureCMB (),relTol=1.0d-3)) &
              & call displayMessage(displayMagenta()//'WARNING: '//displayReset()//'temperatureCMB from transfer function file does not match internal value' )
      end select
      deallocate(cosmologyParametersFile)
      ! Get extrapolation methods.
      extrapolationObject=fileObject         %openGroup('extrapolation')
      wavenumberObject   =extrapolationObject%openGroup('wavenumber'   )
      call wavenumberObject%readAttribute('low' ,limitTypeVar)
      extrapolateWavenumberLow =enumerationExtrapolationTypeEncode(char(limitTypeVar),includesPrefix=.false.)
      call wavenumberObject%readAttribute('high',limitTypeVar)
      extrapolateWavenumberHigh=enumerationExtrapolationTypeEncode(char(limitTypeVar),includesPrefix=.false.)
      ! Read the transfer function from file.
      darkMatterGroup=fileObject%openGroup('darkMatter')
      call fileObject     %readDataset('wavenumber'                                   ,wavenumber)
      write (datasetName,'(f9.4)') self%redshift
      call darkMatterGroup%readDataset('transferFunctionZ'//trim(adjustl(datasetName)),transfer  )
      !$ call hdf5Access%unset()
    end block
    ! Validate the transfer function.
    if (any(transfer == 0.0d0)) call Error_Report('tabulated transfer function contains points at which T(k) = 0 - all points must be non-zero'//{introspection:location})
    if (any(transfer <  0.0d0)) then
       if (self%acceptNegativeValues) then
          transfer=abs(transfer)
       else
          call Error_Report(                                                                                                                                                                  &
          &                 'tabulated transfer function contains points at which T(k) < 0 - all points must be positive'                                           //char(10)             // &
          &                 displayGreen()//"HELP: "//displayReset()//'set '                                                                                                               // &
          &                 '<'//displayBlue()//'acceptNegativeValues'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"true"'//displayReset()//'/> '// &
          &                 'to interpolate in |T(k)| such that negative values are acceptable'                                                                                            // &
          &                 {introspection:location}                                                                                                                                          &
          &                 )
       end if
    end if
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
    ! Determine local minima of the transfer function.
    if (size(transferLogarithmic) > 2) then
       countLocalMinima=0
       do i=2,size(transferLogarithmic)-1
          if     (                                                   &
               &   transferLogarithmic(i) < transferLogarithmic(i-1) &
               &  .and.                                              &
               &   transferLogarithmic(i) < transferLogarithmic(i+1) &
               & ) countLocalMinima=countLocalMinima+1
       end do
       if (allocated(self%wavenumbersLocalMinima_)) deallocate(self%wavenumbersLocalMinima_)
       allocate(self%wavenumbersLocalMinima_(countLocalMinima))
       countLocalMinima=0
       do i=2,size(transferLogarithmic)-1
          if     (                                                   &
               &   transferLogarithmic(i) < transferLogarithmic(i-1) &
               &  .and.                                              &
               &   transferLogarithmic(i) < transferLogarithmic(i+1) &
               & ) then
             countLocalMinima=countLocalMinima+1
             self%wavenumbersLocalMinima_(countLocalMinima)=wavenumber(i)
          end if
       end do
    else
       allocate(self%wavenumbersLocalMinima_(               0))
    end if
    return
  end subroutine fileReadFile

  subroutine fileDestructor(self)
    !!{
    Destructor for the file transfer function class.
    !!}
    implicit none
    type(transferFunctionFile), intent(inout) :: self

    call self%transfer%destroy()
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine fileDestructor

  double precision function fileValue(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    use :: Error             , only : errorStatusOutOfRange               , errorStatusSuccess, Error_Report
    use :: Table_Labels      , only : enumerationExtrapolationTypeDescribe
    use :: ISO_Varying_String, only : var_str
    implicit none
    class           (transferFunctionFile), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber
    integer                                               :: status

    fileValue=exp(self%transfer%interpolate(log(wavenumber),status=status))
    select case (status)
    case (errorStatusSuccess)
       ! Success - nothing to do.
    case (errorStatusOutOfRange)
       call Error_Report(                                                                                                                                                                   &
            &            var_str('wavenumber is outside tabulated range:')                                                                                                     //char(10)// &
            &            'an extrapolation option "abort" was chosen - either extend the range of your tabulated transfer function or choose a different interpolation method:'          // &
            &            enumerationExtrapolationTypeDescribe()                                                                                                                          // &
            &            {introspection:location}                                                                                                                                           &
            &            )
    case default
       call Error_Report('transfer function interpolation failed for unknown reason'//{introspection:location})
    end select
    if (self%transferFunctionReferenceAvailable)                           &
         & fileValue=+                               fileValue             &
         &           *self%transferFunctionReference%    value(wavenumber)
    return
  end function fileValue

  double precision function fileLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    use :: Error             , only : errorStatusOutOfRange               , errorStatusSuccess, Error_Report
    use :: Table_Labels      , only : enumerationExtrapolationTypeDescribe
    use :: ISO_Varying_String, only : var_str
    implicit none
    class           (transferFunctionFile), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber
    integer                                               :: status

    fileLogarithmicDerivative=+self%transfer%interpolateGradient(log(wavenumber),status=status)
    select case (status)
    case (errorStatusSuccess)
       ! Success - nothing to do.
    case (errorStatusOutOfRange)
       call Error_Report(                                                                                                                                                                   &
            &            var_str('wavenumber is outside tabulated range:')                                                                                                     //char(10)// &
            &            'an extrapolation option "abort" was chosen - either extend the range of your tabulated transfer function or choose a different interpolation method:'          // &
            &            enumerationExtrapolationTypeDescribe()                                                                                                                          // &
            &            {introspection:location}                                                                                                                                           &
            &            )
    case default
       call Error_Report('transfer function interpolation failed for unknown reason'//{introspection:location})
    end select
    if (self%transferFunctionReferenceAvailable) &
         & fileLogarithmicDerivative=+                               fileLogarithmicDerivative             &
         &                           +self%transferFunctionReference%    logarithmicDerivative(wavenumber)
    return
  end function fileLogarithmicDerivative
  
  subroutine fileWavenumbersLocalMinima(self,wavenumbers)
    !!{
    Return a list of wavenumbers corresponding to local minima in the transfer function.
    !!}
    implicit none
    class           (transferFunctionFile), intent(inout)                            :: self
    double precision                      , intent(  out), allocatable, dimension(:) :: wavenumbers

    wavenumbers=self%wavenumbersLocalMinima_
    return
  end subroutine fileWavenumbersLocalMinima
  
  double precision function fileHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of two relative to a \gls{cdm} transfer function.
    !!}
    use :: Error, only : Error_Report, errorStatusFail, errorStatusSuccess
    implicit none
    class  (transferFunctionFile), intent(inout), target   :: self
    integer                      , intent(  out), optional :: status

    if (self%massHalfModeAvailable) then
       fileHalfModeMass=self%massHalfMode
       if (present(status)) status=errorStatusSuccess
    else
       fileHalfModeMass=0.0d0
       if (present(status)) then
          status=errorStatusFail
       else
          call Error_Report('not supported by this implementation'//{introspection:location})
       end if
    end if
    return
  end function fileHalfModeMass

  double precision function fileQuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of four relative to a \gls{cdm} transfer function.
    !!}
    use :: Error, only : Error_Report, errorStatusFail
    implicit none
    class  (transferFunctionFile), intent(inout), target   :: self
    integer                      , intent(  out), optional :: status
    
    if (self%massQuarterModeAvailable) then
       fileQuarterModeMass=self%massQuarterMode
       if (present(status)) status=errorStatusSuccess
    else
       fileQuarterModeMass=0.0d0
       if (present(status)) then
          status=errorStatusFail
       else
          call Error_Report('not supported by this implementation'//{introspection:location})
       end if
    end if
    return
  end function fileQuarterModeMass

  double precision function fileFractionModeMass(self,fraction,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    reduced by {\normalfont \ttfamily fraction} relative to a \gls{cdm} transfer function.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Error                   , only : Error_Report, errorStatusFail
    implicit none
    class           (transferFunctionFile), intent(inout), target   :: self
    double precision                      , intent(in   )           :: fraction
    integer                               , intent(  out), optional :: status
    double precision                                                :: fractionCurrent, fractionPrevious, &
         &                                                             wavenumber     , matterDensity
    integer                                                         :: j
    logical                                                         :: modeFound
    
    if (present(status)) status=errorStatusSuccess
    if (self%transferFunctionReferenceAvailable) then
       matterDensity          =+self%cosmologyParameters_%OmegaMatter    () &
            &                  *self%cosmologyParameters_%densityCritical()
       ! Step from large to small scales until the target fraction is reached.
       fractionCurrent =0.0d0
       fractionPrevious=0.0d0
       modeFound       =.false.
       do j=1,self%transfer%size()
          fractionPrevious=+fractionCurrent
          fractionCurrent =+exp(self%transfer%y(j))
          if (fractionCurrent < fraction) then
             if (j > 1) then            
                ! Interpolate to find the exact wavenumber.
                wavenumber=+exp(                                                                      &
                     &          +(+              fraction          -              fractionPrevious  ) &
                     &          /(+              fractionCurrent   -              fractionPrevious  ) &
                     &          *(+self%transfer%x              (j)-self%transfer%x            (j-1)) &
                     &          +                                   self%transfer%x            (j-1)  &
                     &         )
                ! Compute the mode mass.
                modeFound           =.true.
                fileFractionModeMass=+4.0d0         &
                     &               *Pi            &
                     &               /3.0d0         &
                     &               *matterDensity &
                     &               *(             &
                     &                 +Pi          &
                     &                 /wavenumber  &
                     &                )**3
             end if
             exit
          end if
       end do
       if (.not.modeFound) then
          fileFractionModeMass=0.0d0
          if (present(status)) then
             status=errorStatusFail
          else
             call Error_Report('not supported by this implementation'//{introspection:location})
          end if
       end if
    else
       fileFractionModeMass=0.0d0
       if (present(status)) then
          status=errorStatusFail
       else
          call Error_Report('not supported by this implementation'//{introspection:location})
       end if
    end if
    return
  end function fileFractionModeMass

  double precision function fileEpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionFile), intent(inout) :: self

    fileEpochTime=self%time
    return
  end function fileEpochTime

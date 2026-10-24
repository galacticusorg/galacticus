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

  !!{RST
  A transferred primordial power spectrum class which reads the power spectrum from file and interpolates in it.
  !!}
  
  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParameters, cosmologyParametersClass
  use :: Numerical_Interpolation, only : interpolator
  
  !![
  <powerSpectrumPrimordialTransferred name="powerSpectrumPrimordialTransferredFile" docformat="rst">
    <description>
    A transferred primordial spectrum class which reads the power spectrum from file and interpolates in it. The HDF5 file containing the data should have the following structure:

    .. code-block:: none

        HDF5 "powerSpectrum.hdf5" {
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
           ATTRIBUTE "extrapolationWavenumber" {
              DATATYPE  H5T_STRING {
                 STRSIZE 11;
                 STRPAD H5T_STR_NULLTERM;
                 CSET H5T_CSET_ASCII;
                 CTYPE H5T_C_S1;
              }
              DATASPACE  SCALAR
           }
           ATTRIBUTE "extrapolationRedshift" {
              DATATYPE  H5T_STRING {
                 STRSIZE 11;
                 STRPAD H5T_STR_NULLTERM;
                 CSET H5T_CSET_ASCII;
                 CTYPE H5T_C_S1;
              }
              DATASPACE  SCALAR
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
           DATASET "wavenumber" {
              DATATYPE  H5T_IEEE_F64LE
              DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ) }
           }
           DATASET "redshift" {
              DATATYPE  H5T_IEEE_F64LE
              DATASPACE  SIMPLE { ( 10 ) / ( 10 ) }
           }
           DATASET "power" {
              DATATYPE  H5T_IEEE_F64LE
              DATASPACE  SIMPLE { ( 10 ) / ( 10 ), ( 1000 ) / ( 1000 ) }
           }
        }
        }

    The `power` dataset should tabulate the transferred power spectrum as a function of :math:`(z,k)` where :math:`z` is redshift and :math:`k` is wavenumber.

    The tabulated power spectrum is multiplied by a normalization determined by a :galacticus-class:`powerSpectrumNormalizerClass`
    object. If the tabulated power spectrum is already normalized use a :galacticus-class:`powerSpectrumNormalizerIdentity`.
    </description>
    <deepCopy>
      <ignore variables="powerSpectrumNormalizer_"/>
    </deepCopy>
    <assignment>
      <functionClass variables="powerSpectrumNormalizer_"/>
    </assignment>
    <allowedParameters>
      <object variable="powerSpectrumNormalizer_" name="powerSpectrumNormalizer" source="parameters" function="Power_Spectrum_Normalizer_Allowed_Parameters_" module="Functions_Global"/>
    </allowedParameters>
    <runTimeFileDependencies paths="fileName"/>
  </powerSpectrumPrimordialTransferred>
  !!]
  type, extends(powerSpectrumPrimordialTransferredClass) :: powerSpectrumPrimordialTransferredFile
     !!{RST
     Implements a file transferred primordial power spectrum.
     !!}
     private
     class           (cosmologyParametersClass), pointer                     :: cosmologyParameters_     => null()
     class           (cosmologyFunctionsClass ), pointer                     :: cosmologyFunctions_      => null()
     ! The normalizer is held as an unlimited polymorphic object, and operated on through global functions, to avoid a circular
     ! module dependency (see Power_Spectrum_Normalizer_Utilities).
     class           (*                       ), pointer                     :: powerSpectrumNormalizer_ => null()
     double precision                                                        :: normalization_
     type            (varying_string          )                              :: fileName
     double precision                          , allocatable, dimension(:  ) :: timeLogarithmic                   , wavenumberLogarithmic
     double precision                          , allocatable, dimension(:,:) :: powerLogarithmic
     type            (interpolator            )                              :: interpolatorWavenumber            , interpolatorTime
   contains
     !![
     <methods docformat="rst">
       <method description="Read the named power spectrum file."                                                                     method="readFile"         />
       <method description="Return the unnormalized transferred primordial power spectrum at the given ``wavenumber`` and ``time``." method="powerUnnormalized"/>
       <method description="Return the normalization of the power spectrum, computing it if necessary."                              method="normalization"    />
     </methods>
     !!]
     final     ::                                fileDestructor
     procedure :: readFile                    => fileReadFile
     procedure :: power                       => filePower
     procedure :: powerUnnormalized           => filePowerUnnormalized
     procedure :: normalization               => fileNormalization
     procedure :: logarithmicDerivative       => fileLogarithmicDerivative
     procedure :: growthIsWavenumberDependent => fileGrowthIsWavenumberDependent
     procedure :: descriptor                  => fileDescriptor
     procedure :: deepCopy                    => fileDeepCopy
     procedure :: deepCopyReset               => fileDeepCopyReset
     procedure :: deepCopyFinalize            => fileDeepCopyFinalize
  end type powerSpectrumPrimordialTransferredFile

  interface powerSpectrumPrimordialTransferredFile
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumPrimordialTransferredFile` transferred primordial power spectrum class.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface powerSpectrumPrimordialTransferredFile

  ! Current file format version for power spectrum files.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function fileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumPrimordialTransferredFile` transferred primordial power spectrum class which takes a parameter set as input.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Construct_, Power_Spectrum_Normalizer_Release_
    use :: Input_Parameters, only : inputParameter                      , inputParameters
    implicit none
    type (powerSpectrumPrimordialTransferredFile)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(cosmologyParametersClass              ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class(*                                     ), pointer       :: powerSpectrumNormalizer_
    type (varying_string                        )                :: fileName
    
    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <source>parameters</source>
      <description>
      The name of the file from which to read a tabulated transfer function.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    powerSpectrumNormalizer_ => Power_Spectrum_Normalizer_Construct_(parameters)
    self=powerSpectrumPrimordialTransferredFile(fileName,cosmologyParameters_,cosmologyFunctions_,powerSpectrumNormalizer_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    call Power_Spectrum_Normalizer_Release_(powerSpectrumNormalizer_)
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName,cosmologyParameters_,cosmologyFunctions_,powerSpectrumNormalizer_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumPrimordialTransferredFile` transferred primordial power
    spectrum class. The ``powerSpectrumNormalizer_`` argument must be of :galacticus-class:`powerSpectrumNormalizerClass` class.
    !!}
    implicit none
    type (powerSpectrumPrimordialTransferredFile)                        :: self
    type (varying_string                        ), intent(in   )         :: fileName
    class(cosmologyParametersClass              ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class(*                                     ), intent(in   ), target :: powerSpectrumNormalizer_
    !![
    <constructorAssign variables="fileName, *cosmologyParameters_, *cosmologyFunctions_, *powerSpectrumNormalizer_"/>
    !!]

    call self%readFile(fileName)
    self%normalization_=-1.0d0
    return
  end function fileConstructorInternal

  subroutine fileDestructor(self)
    !!{RST
    Destructor for the file transfer function class.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Release_
    implicit none
    type(powerSpectrumPrimordialTransferredFile), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    if (associated(self%powerSpectrumNormalizer_)) call Power_Spectrum_Normalizer_Release_(self%powerSpectrumNormalizer_)
    return
  end subroutine fileDestructor

  subroutine fileReadFile(self,fileName)
    !!{RST
    Read a tabulated power spectrum from file.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Cosmology_Parameters   , only : cosmologyParametersSimple
    use            :: Display                , only : displayMessage                  , displayMagenta                    , displayGreen                , displayReset
    use            :: File_Utilities         , only : File_Exists
    use            :: Error                  , only : Error_Report
    use            :: HDF5_Access            , only : hdf5Access
    use            :: IO_HDF5                , only : hdf5File                        , hdf5Group
    use            :: Numerical_Comparison   , only : Values_Differ
    use            :: Numerical_Interpolation, only : GSL_Interp_cSpline
    use            :: Table_Labels           , only : enumerationExtrapolationTypeType, enumerationExtrapolationTypeEncode, extrapolationTypeExtrapolate
    use            :: ISO_Varying_String     , only : var_str                         , operator(//)
    use            :: Sorting                , only : sortIndex
    use            :: String_Handling        , only : operator(//)
    implicit none
    class           (powerSpectrumPrimordialTransferredFile), intent(inout)               :: self
    type            (varying_string                        ), intent(in   )               :: fileName
    class           (cosmologyParametersClass              ), pointer                     :: cosmologyParametersFile
    double precision                                        , dimension(:  ), allocatable :: wavenumber             , redshift
    double precision                                        , dimension(:,:), allocatable :: power
    integer         (c_size_t                              ), dimension(:  ), allocatable :: order
    double precision                                                                      :: HubbleConstant         , OmegaBaryon        , &
         &                                                                                   OmegaMatter            , OmegaDarkEnergy    , &
         &                                                                                   temperatureCMB
    type            (enumerationExtrapolationTypeType      )                              :: extrapolateWavenumber  , extrapolateRedshift
    integer                                                                               :: versionNumber          , i
    type            (hdf5File                              )                              :: fileObject
    type            (hdf5Group                             )                              :: parametersObject
    type            (varying_string                        )                              :: limitTypeVar           , message

    ! Check that the file exists.
    if (.not.File_Exists(fileName)) call Error_Report("file '"//char(fileName)//"' does not exist"//{introspection:location})
    ! Open and read the HDF5 data file.
    !$ call hdf5Access%set()
    fileObject=hdf5File(fileName,readOnly=.true.)
    ! Check that the file has the correct format version number.
    call fileObject%readAttribute('fileFormat',versionNumber,allowPseudoScalar=.true.)
    if (versionNumber /= fileFormatVersionCurrent) call Error_Report(var_str('file has the incorrect format version number (expected fileFormat=1, found fileFormat=')//versionNumber//')'//{introspection:location})
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
            & call displayMessage(displayMagenta()//"WARNING:"//displayReset()//' OmegaBaryon from transfer function file does not match internal value'    )
       if (Values_Differ(cosmologyParametersFile%OmegaMatter    (),self%cosmologyParameters_%OmegaMatter    (),absTol=1.0d-3)) &
            & call displayMessage(displayMagenta()//"WARNING:"//displayReset()//' OmegaMatter from transfer function file does not match internal value'    )
       if (Values_Differ(cosmologyParametersFile%OmegaDarkEnergy(),self%cosmologyParameters_%OmegaDarkEnergy(),absTol=1.0d-3)) &
            & call displayMessage(displayMagenta()//"WARNING:"//displayReset()//' OmegaDarkEnergy from transfer function file does not match internal value')
       if (Values_Differ(cosmologyParametersFile%HubbleConstant (),self%cosmologyParameters_%HubbleConstant (),relTol=1.0d-3)) &
            & call displayMessage(displayMagenta()//"WARNING:"//displayReset()//' HubbleConstant from transfer function file does not match internal value' )
       if (Values_Differ(cosmologyParametersFile%temperatureCMB (),self%cosmologyParameters_%temperatureCMB (),relTol=1.0d-3)) &
            & call displayMessage(displayMagenta()//"WARNING:"//displayReset()//' temperatureCMB from transfer function file does not match internal value' )
    end select
    deallocate(cosmologyParametersFile)
    ! Get extrapolation methods.
    call fileObject%readAttribute('extrapolationWavenumber',limitTypeVar)
    extrapolateWavenumber=enumerationExtrapolationTypeEncode(char(limitTypeVar),includesPrefix=.false.)
    call fileObject%readAttribute('extrapolationRedshift'  ,limitTypeVar)
    extrapolateRedshift  =enumerationExtrapolationTypeEncode(char(limitTypeVar),includesPrefix=.false.)
    ! Read the power spectrum from file.
    call fileObject%readDataset('wavenumber',wavenumber)
    call fileObject%readDataset('redshift'  ,redshift  )
    call fileObject%readDataset('power'     ,power     )
    !$ call hdf5Access%unset()
    ! Validate data.
    if (any(shape(power) /= [size(wavenumber),size(redshift)])) then
       message='dimensions of `power` are mismatched with those of `redshift` and `wavenumber` arrays'
       if (all(shape(power) == [size(redshift),size(wavenumber)])) message=message//' (it appears that the axes of your `power` array may be transposed - they should be [redshift,wavenumber] in the HDF5 dataset)'
       call Error_Report(message//{introspection:location})
    end if
    ! Construct the tabulated power spectrum and interpolators. Note that the tabulated power must be in order of increasing time, so sort on redshift and index in reverse.
    order=sortIndex(redshift)
    allocate(self%wavenumberLogarithmic(size(wavenumber)               ))
    allocate(self%      timeLogarithmic(                 size(redshift)))
    allocate(self%     powerLogarithmic(size(wavenumber),size(redshift)))
    self%wavenumberLogarithmic=log(wavenumber)
    do i=1,size(redshift)
       self%timeLogarithmic (  i)=log(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift(  order(size(redshift)+1-i)))))
       self%powerLogarithmic(:,i)=log(                                                                                            power(:,order(size(redshift)+1-i))  )
    end do
    self%interpolatorWavenumber=interpolator(self%wavenumberLogarithmic,extrapolationType=extrapolateWavenumber)
    self%interpolatorTime      =interpolator(self%timeLogarithmic      ,extrapolationType=extrapolateRedshift  )
    ! Warn about potential problems extrapolating power to large wavenumber.
    if     (                                                                                                                                                                                                                                           &
         &            extrapolateWavenumber                   ==      extrapolationTypeExtrapolate                                                                                                                                                     &
         &  .and.                                                                                                                                                                                                                                      &
         &   any(self%powerLogarithmic     (:,size(redshift)) >  self%powerLogarithmic            (:,size(redshift)-1))                                                                                                                                &
         & ) call Warn(                                                                                                                                                                                                                                &
         &             displayMagenta()// "WARNING:"//displayReset()//" You have requested extrapolation of the power spectrum in wavenumber, but your power spectrum is increasing with wavenumber at the largest wavenumbers tabulated."//char(10)// &
         &                               "         "                //" This can lead to excessive power when extrapolated to very large wavenumbers."                                                                                    //char(10)// &
         &             displayGreen  ()//"    HELP:"//displayReset()//" Either ensure that your power spectrum is decreasing at the highest wavenumbers or, if your power spectrum is strongly truncated at large wavenumber, consider "            // &
         &                                                            "setting the `extrapolationWavenumber` attribute to `zero` to simply truncate the power spectrum beyond the maximum tabulated wavenumber."                                       &
         &            )
    return
  end subroutine fileReadFile
  
  double precision function filePower(self,wavenumber,time)
    !!{RST
    Return the transferred primordial power spectrum at the given ``wavenumber``.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredFile), intent(inout) :: self
    double precision                                        , intent(in   ) :: wavenumber, time

    filePower=+self%normalization    (                ) &
         &    *self%powerUnnormalized(wavenumber,time)
    return
  end function filePower

  double precision function fileNormalization(self)
    !!{RST
    Return the normalization of the power spectrum, computing it (using the normalizer) on the first call.
    !!}
    implicit none
    class(powerSpectrumPrimordialTransferredFile), intent(inout) :: self

    ! The computation is kept in a separate procedure, as the internal procedure which it passes to the normalizer requires
    ! a trampoline, which would otherwise be set up on every call to this function (and so on every evaluation of the power
    ! spectrum), not only when the normalization must be computed.
    if (self%normalization_ < 0.0d0) call fileNormalizationCompute(self)
    fileNormalization=self%normalization_
    return
  end function fileNormalization

  subroutine fileNormalizationCompute(self)
    !!{RST
    Compute the normalization of the power spectrum using the normalizer.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Normalization_
    implicit none
    class(powerSpectrumPrimordialTransferredFile), intent(inout), target :: self

    self%normalization_=Power_Spectrum_Normalizer_Normalization_(self%powerSpectrumNormalizer_,powerUnnormalized)
    return

  contains

    double precision function powerUnnormalized(wavenumber,time)
      !!{RST
      Return the unnormalized transferred primordial power spectrum, for use by the normalizer.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber, time

      powerUnnormalized=self%powerUnnormalized(wavenumber,time)
      return
    end function powerUnnormalized

  end subroutine fileNormalizationCompute

  double precision function filePowerUnnormalized(self,wavenumber,time) result(power)
    !!{RST
    Return the unnormalized (i.e. as tabulated) transferred primordial power spectrum at the given ``wavenumber``.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredFile), intent(inout)  :: self
    double precision                                        , intent(in   )  :: wavenumber , time
    integer         (c_size_t                              )                 :: iWavenumber, iTime, &
         &                                                                      jWavenumber, jTime
    double precision                                        , dimension(0:1) :: hWavenumber, hTime

    call self%interpolatorWavenumber%linearFactors(log(wavenumber),iWavenumber,hWavenumber) 
    call self%interpolatorTime      %linearFactors(log(time      ),iTime      ,hTime      )
    power=0.0d0
    do jWavenumber=0,1
       do jTime   =0,1
          power  =+     power                                                 &
               &  +self%powerLogarithmic(iWavenumber+jWavenumber,iTime+jTime) &
               &  *     hWavenumber     (            jWavenumber            ) &
               &  *     hTime           (                              jTime)
       end do
    end do
    power=exp(power)
    return
  end function filePowerUnnormalized

  double precision function fileLogarithmicDerivative(self,wavenumber,time)
    !!{RST
    Return the logarithmic derivative of the transferred primordial power spectrum at the given ``wavenumber``.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredFile), intent(inout)  :: self
    double precision                                        , intent(in   )  :: wavenumber , time
    integer         (c_size_t                              )                 :: iWavenumber, iTime, &
         &                                                                      jTime
    double precision                                        , dimension(0:1) :: hWavenumber, hTime
    
    call self%interpolatorWavenumber%linearFactors(log(wavenumber),iWavenumber,hWavenumber) 
    call self%interpolatorTime      %linearFactors(log(time      ),iTime      ,hTime      )
    fileLogarithmicDerivative=0.0d0
    do jTime=0,1
       fileLogarithmicDerivative=+            fileLogarithmicDerivative                                     &
            &                    +(                                                                         &
            &                      +self%     powerLogarithmic         (iWavenumber+1_c_size_t,iTime+jTime) &
            &                      -self%     powerLogarithmic         (iWavenumber+0_c_size_t,iTime+jTime) &
            &                     )                                                                         &
            &                    /(                                                                         &
            &                      +self%wavenumberLogarithmic         (iWavenumber+1_c_size_t            ) &
            &                      -self%wavenumberLogarithmic         (iWavenumber+0_c_size_t            ) &
            &                     )                                                                         &
            &                    *       hTime                         (                             jTime)
    end do
    return
  end function fileLogarithmicDerivative

  logical function fileGrowthIsWavenumberDependent(self)
    !!{RST
    Return true if the growth of the power spectrum is wavenumber-dependent.
    !!}
    implicit none
    class(powerSpectrumPrimordialTransferredFile), intent(inout) :: self

    fileGrowthIsWavenumberDependent=.true.
    return
  end function fileGrowthIsWavenumberDependent

  subroutine fileDescriptor(self,descriptor,includeClass,includeFileModificationTimes,parameterName)
    !!{RST
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Error             , only : Error_Report                         , errorStatusSuccess     , errorStatusNotExist
    use :: File_Utilities    , only : File_Modification_Time
    use :: Functions_Global  , only : Power_Spectrum_Normalizer_Descriptor_
    use :: Input_Parameters  , only : inputParameters                      , descriptorParameterName
    use :: ISO_Varying_String, only : char                                 , operator(//)
    implicit none
    class    (powerSpectrumPrimordialTransferredFile), intent(inout)           :: self
    type     (inputParameters                       ), intent(inout)           :: descriptor
    logical                                          , intent(in   ), optional :: includeClass    , includeFileModificationTimes
    character(len=*                                 ), intent(in   ), optional :: parameterName
    type     (inputParameters                       )                          :: parameters
    character(len=30                                )                          :: timeModification
    integer                                                                    :: status

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter(descriptorParameterName('powerSpectrumPrimordialTransferred',parameterName),'file')
    parameters=descriptor%subparameters(descriptorParameterName('powerSpectrumPrimordialTransferred',parameterName))
    call parameters%addParameter('fileName',char(self%fileName))
    call self%cosmologyParameters_%descriptor(parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    call self%cosmologyFunctions_ %descriptor(parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    call Power_Spectrum_Normalizer_Descriptor_(self%powerSpectrumNormalizer_,parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
    if (present(includeFileModificationTimes)) then
       if (includeFileModificationTimes) then
          timeModification=File_Modification_Time(self%fileName,status)
          if      (status == errorStatusSuccess ) then
             call descriptor%addParameter('runTimeFileDependency1',char(self%fileName//": "//trim(timeModification)))
          else if (status /= errorStatusNotExist) then
             call Error_Report('unable to get file modification time'//{introspection:location})
          end if
       end if
    end if
    return
  end subroutine fileDescriptor

  subroutine fileDeepCopyReset(self)
    !!{RST
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Deep_Copy_Reset_
    implicit none
    class(powerSpectrumPrimordialTransferredFile), intent(inout) :: self

    self%copiedSelf => null()
    if (associated(self%cosmologyParameters_    )) call self%cosmologyParameters_%deepCopyReset()
    if (associated(self%cosmologyFunctions_     )) call self%cosmologyFunctions_ %deepCopyReset()
    if (associated(self%powerSpectrumNormalizer_)) call Power_Spectrum_Normalizer_Deep_Copy_Reset_(self%powerSpectrumNormalizer_)
    return
  end subroutine fileDeepCopyReset

  subroutine fileDeepCopyFinalize(self)
    !!{RST
    Finalize a deep copy of the object.
    !!}
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Deep_Copy_Finalize_
    implicit none
    class(powerSpectrumPrimordialTransferredFile), intent(inout) :: self

    if (associated(self%cosmologyParameters_    )) call self%cosmologyParameters_%deepCopyFinalize()
    if (associated(self%cosmologyFunctions_     )) call self%cosmologyFunctions_ %deepCopyFinalize()
    if (associated(self%powerSpectrumNormalizer_)) call Power_Spectrum_Normalizer_Deep_Copy_Finalize_(self%powerSpectrumNormalizer_)
    return
  end subroutine fileDeepCopyFinalize

  subroutine fileDeepCopy(self,destination)
    !!{RST
    Perform a deep copy of the object.
    !!}
    use :: Error           , only : Error_Report
    use :: Function_Classes, only : functionClass
    use :: Functions_Global, only : Power_Spectrum_Normalizer_Deep_Copy_
    implicit none
    class  (powerSpectrumPrimordialTransferredFile ), intent(inout), target :: self
    class  (powerSpectrumPrimordialTransferredClass), intent(inout)         :: destination
    integer                                                                 :: referenceCount

    call self%deepCopy_(destination)
    select type (destination)
    type is (powerSpectrumPrimordialTransferredFile)
       ! The (shallow) assignment of self to destination in deepCopy_() acquired a reference to our normalizer. Release that
       ! reference, as the destination is given its own (deep) copy of the normalizer.
       if (associated(destination%powerSpectrumNormalizer_)) then
          select type (powerSpectrumNormalizer_ => destination%powerSpectrumNormalizer_)
          class is (functionClass)
             referenceCount=powerSpectrumNormalizer_%referenceCountDecrement()
          end select
       end if
       nullify(destination%powerSpectrumNormalizer_)
       if (associated(self%powerSpectrumNormalizer_)) call Power_Spectrum_Normalizer_Deep_Copy_(self%powerSpectrumNormalizer_,destination%powerSpectrumNormalizer_)
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine fileDeepCopy

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
  A transferred primordial power spectrum class which reads the power spectrum from file and interpolates in it.
  !!}
  
  use :: Cosmology_Functions    , only : cosmologyFunctions , cosmologyFunctionsClass
  use :: Cosmology_Parameters   , only : cosmologyParameters, cosmologyParametersClass
  use :: Numerical_Interpolation, only : interpolator
  
  !![
  <powerSpectrumPrimordialTransferred name="powerSpectrumPrimordialTransferredFile">
    <description>
      A transferred primordial spectrum class which reads the power spectrum from file and interpolates in it. The HDF5 file
      containing the data should have the following structure:
       \begin{verbatim}
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
             DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ) }
          }
          DATASET "power" {
             DATATYPE  H5T_IEEE_F64LE
             DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ), ( 1000 ) / ( 1000 ) }
          }
       }
       }
       \end{verbatim}
      The `power` dataset should tabulate the transferred power spectrum as a function of $(k,z)$ where $k$ is wavenumber and $z$
      is redshift.
    </description>
    <runTimeFileDependencies paths="fileName"/>
  </powerSpectrumPrimordialTransferred>
  !!]
  type, extends(powerSpectrumPrimordialTransferredClass) :: powerSpectrumPrimordialTransferredFile
     !!{
     Implements a file transferred primordial power spectrum.
     !!}
     private
     class           (cosmologyParametersClass), pointer                     :: cosmologyParameters_   => null()
     class           (cosmologyFunctionsClass ), pointer                     :: cosmologyFunctions_    => null()
     type            (varying_string          )                              :: fileName
     double precision                          , allocatable, dimension(:  ) :: timeLogarithmic                 , wavenumberLogarithmic
     double precision                          , allocatable, dimension(:,:) :: powerLogarithmic
     type            (interpolator            )                              :: interpolatorWavenumber          , interpolatorTime
   contains
     !![
     <methods>
       <method description="Read the named power spectrum file." method="readFile" />
     </methods>
     !!]
     final     ::                                fileDestructor
     procedure :: readFile                    => fileReadFile
     procedure :: power                       => filePower
     procedure :: logarithmicDerivative       => fileLogarithmicDerivative
     procedure :: growthIsWavenumberDependent => fileGrowthIsWavenumberDependent
  end type powerSpectrumPrimordialTransferredFile

  interface powerSpectrumPrimordialTransferredFile
     !!{
     Constructors for the \refClass{powerSpectrumPrimordialTransferredFile} transferred primordial power spectrum class.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface powerSpectrumPrimordialTransferredFile

  ! Current file format version for power spectrum files.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function fileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{powerSpectrumPrimordialTransferredFile} transferred primordial power spectrum class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (powerSpectrumPrimordialTransferredFile)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(cosmologyParametersClass              ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    type (varying_string                        )                :: fileName
    
    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file from which to read a tabulated transfer function.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !!]
    self=powerSpectrumPrimordialTransferredFile(fileName,cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumPrimordialTransferredFile} transferred primordial power spectrum class.
    !!}
    implicit none
    type (powerSpectrumPrimordialTransferredFile)                        :: self
    type (varying_string                        ), intent(in   )         :: fileName
    class(cosmologyParametersClass              ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="fileName, *cosmologyParameters_, *cosmologyFunctions_"/>
    !!]

    call self%readFile(fileName)
    return
  end function fileConstructorInternal

  subroutine fileDestructor(self)
    !!{
    Destructor for the file transfer function class.
    !!}
    implicit none
    type(powerSpectrumPrimordialTransferredFile), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine fileDestructor

  subroutine fileReadFile(self,fileName)
    !!{
    Read a tabulated power spectrum from file.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Cosmology_Parameters   , only : cosmologyParametersSimple
    use            :: Display                , only : displayMessage                  , displayMagenta                    , displayGreen                , displayReset
    use            :: File_Utilities         , only : File_Name_Expand                , File_Exists
    use            :: Error                  , only : Error_Report
    use            :: HDF5_Access            , only : hdf5Access
    use            :: IO_HDF5                , only : hdf5Object
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
    type            (hdf5Object                            )                              :: fileObject             , parametersObject
    type            (varying_string                        )                              :: limitTypeVar

    ! Check that the file exists.
    if (.not.File_Exists(File_Name_Expand(char(fileName)))) call Error_Report("file '"//char(fileName)//"' does not exist"//{introspection:location})
    ! Open and read the HDF5 data file.
    !$ call hdf5Access%set()
    fileObject=hdf5Object(char(File_Name_Expand(char(fileName))),readOnly=.true.)
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
    !!{
    Return the transferred primordial power spectrum at the given {\normalfont \ttfamily
    wavenumber}.
    !!}
    implicit none
    class           (powerSpectrumPrimordialTransferredFile), intent(inout)  :: self
    double precision                                        , intent(in   )  :: wavenumber , time
    integer         (c_size_t                              )                 :: iWavenumber, iTime, &
         &                                                                      jWavenumber, jTime
    double precision                                        , dimension(0:1) :: hWavenumber, hTime

    call self%interpolatorWavenumber%linearFactors(log(wavenumber),iWavenumber,hWavenumber) 
    call self%interpolatorTime      %linearFactors(log(time      ),iTime      ,hTime      )
    filePower=0.0d0
    do jWavenumber=0,1
       do jTime   =0,1
          filePower=+     filePower                                             &
               &    +self%powerLogarithmic(iWavenumber+jWavenumber,iTime+jTime) &
               &    *     hWavenumber     (            jWavenumber            ) &
               &    *     hTime           (                              jTime)
       end do
    end do
    filePower=exp(filePower)
    return
  end function filePower

  double precision function fileLogarithmicDerivative(self,wavenumber,time)
    !!{
    Return the logarithmic derivative of the transferred primordial power spectrum at the
    given {\normalfont \ttfamily wavenumber}.
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
    !!{
    Return true if the growth of the power spectrum is wavenumber-dependent.
    !!}
    implicit none
    class(powerSpectrumPrimordialTransferredFile), intent(inout) :: self

    fileGrowthIsWavenumberDependent=.true.
    return
  end function fileGrowthIsWavenumberDependent

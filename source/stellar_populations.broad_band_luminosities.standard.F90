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

  !% Implements the standard stellar populations broad band luminosities class.
  
  use, intrinsic :: ISO_C_Binding                         , only : c_size_t
  use            :: Locks                                 , only : ompReadWriteLock
  use            :: Numerical_Interpolation               , only : interpolator
  use            :: Stellar_Population_Spectra_Postprocess, only : stellarPopulationSpectraPostprocessorClass
  use            :: Stellar_Population_Spectra            , only : stellarPopulationSpectraClass

  type luminosityTable
     !% Structure for holding tables of simple stellar population luminosities.
     integer                                                       :: agesCount           , metallicitiesCount
     integer         (c_size_t    )                                :: isTabulatedMaximum=0
     logical                       , allocatable, dimension(:    ) :: isTabulated
     double precision              , allocatable, dimension(:    ) :: age                 , metallicity
     double precision              , allocatable, dimension(:,:,:) :: luminosity
     type            (interpolator)                                :: interpolatorAge     , interpolatorMetallicity
  end type luminosityTable

  !# <stellarPopulationBroadBandLuminosities name="stellarPopulationBroadBandLuminositiesStandard">
  !#  <description>The standard stellar populations broad band luminosities class.</description>
  !# </stellarPopulationBroadBandLuminosities>
  type, extends(stellarPopulationBroadBandLuminositiesClass) :: stellarPopulationBroadBandLuminositiesStandard
     !% The standard stellar population broad band luminosities class.
     private
     type            (luminosityTable ), allocatable, dimension(:) :: luminosityTables
     double precision                                              :: integrationToleranceRelative
     logical                                                       :: integrationToleranceDegrade
     logical                                                       :: storeToFile
     type            (varying_string  )                            :: storeDirectory
     logical                                                       :: maximumAgeExceededIsFatal
     type            (ompReadWriteLock)                            :: luminosityTableLock
   contains
     !@ <objectMethods>
     !@   <object>stellarPopulationBroadBandLuminositiesStandard</object>
     !@   <objectMethod>
     !@     <method>tabulate</method>
     !@     <type>\void</type>
     !@     <arguments>\intone\ luminosityIndex\argin,\intone\ filterIndex\argin, \textcolor{red}{\textless type(stellarLuminosities)(:)\textgreater} stellarPopulationSpectraPostprocessor_\argin, \textcolor{red}{\textless type(stellarLuminosities)\textgreater} stellarPopulation_\arginout, \doublezero\ redshift\argin</arguments>
     !@     <description>Tabulate broad band stellar luminosities.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: luminosities     => standardLuminosities
     procedure :: luminosityTracks => standardLuminosityTracks
     procedure :: tabulate         => standardTabulate
  end type stellarPopulationBroadBandLuminositiesStandard

  interface stellarPopulationBroadBandLuminositiesStandard
     !% Constructors for the {\normalfont \ttfamily standard} stellar population broad band luminosities class.
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface stellarPopulationBroadBandLuminositiesStandard

  ! Module scope variables used in integrations.
  double precision                                                      :: standardAge                                  , standardRedshift
  integer                                                               :: standardFilterIndex
  integer         (c_size_t                                  )          :: standardPopulationID
  type            (abundances                                )          :: standardAbundances
  class           (stellarPopulationSpectraPostprocessorClass), pointer :: standardStellarPopulationSpectraPostprocessor
  class           (stellarPopulationSpectraClass             ), pointer :: standardStellarPopulationSpectra
  !$omp threadprivate(standardAge,standardRedshift,standardAbundances,standardFilterIndex,standardPopulationID,standardStellarPopulationSpectraPostprocessor,standardStellarPopulationSpectra)

contains

  function standardConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily standard} stellar population broad band luminosities class which takes a
    !% parameter set as input.
    use :: Input_Parameters, only : inputParameters
    use :: Galacticus_Paths, only : galacticusPath , pathTypeDataDynamic
    implicit none
    type            (stellarPopulationBroadBandLuminositiesStandard)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    double precision                                                                :: integrationToleranceRelative
    logical                                                                         :: storeToFile                 , integrationToleranceDegrade, &
         &                                                                             maximumAgeExceededIsFatal
    type            (varying_string                                )                :: storeDirectory

    !# <inputParameter>
    !#   <name>integrationToleranceRelative</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>4.0d-3</defaultValue>
    !#   <description>The relative tolerance used when integrating the flux of stellar populations through filters.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>integrationToleranceDegrade</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If {\normalfont \ttfamily true}, automatically degrade the relative tolerance used when integrating the flux of stellar populations through filters to ensure convergence.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>storeToFile</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>Specifies whether or not stellar populations luminosities (integrated under a filter) should be stored to file for rapid reuse.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>storeDirectory</name>
    !#   <defaultValue>galacticusPath(pathTypeDataDynamic)//'stellarPopulations'</defaultValue>
    !#   <description>Specifies the directory to which stellar populations luminosities (integrated under a filter) should be stored to file for rapid reuse.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>maximumAgeExceededIsFatal</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>Specifies whether or not exceeding the maximum available age of the stellar population is fatal.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    self=stellarPopulationBroadBandLuminositiesStandard(integrationToleranceRelative,integrationToleranceDegrade,maximumAgeExceededIsFatal,storeToFile,storeDirectory)
    !# <inputParametersValidate source="parameters"/>
    return
  end function standardConstructorParameters

  function standardConstructorInternal(integrationToleranceRelative,integrationToleranceDegrade,maximumAgeExceededIsFatal,storeToFile,storeDirectory) result(self)
    !% Internal constructor for the {\normalfont \ttfamily standard} stellar population broad band luminosities class.
    implicit none
    type            (stellarPopulationBroadBandLuminositiesStandard)                :: self
    double precision                                                , intent(in   ) :: integrationToleranceRelative
    logical                                                         , intent(in   ) :: storeToFile                 , integrationToleranceDegrade, &
         &                                                                             maximumAgeExceededIsFatal
    type            (varying_string                                ), intent(in   ) :: storeDirectory
    !# <constructorAssign variables="integrationToleranceRelative, integrationToleranceDegrade, maximumAgeExceededIsFatal,storeToFile,storeDirectory"/>
    
    self%luminosityTableLock=ompReadWriteLock()
    return
  end function standardConstructorInternal

  function standardLuminosities(self,luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,abundancesStellar,age,redshift)
    !% Returns the luminosity for a $1 M_\odot$ simple {\normalfont \ttfamily stellarPopulation\_} of given {\normalfont \ttfamily
    !% abundances} and {\normalfont \ttfamily age} and observed through the filter specified by {\normalfont \ttfamily
    !% filterIndex}.
    use            :: Abundances_Structure   , only : Abundances_Get_Metallicity, logMetallicityZero, metallicityTypeLogarithmicByMassSolar
    use            :: Galacticus_Error       , only : Galacticus_Error_Report
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Stellar_Populations    , only : stellarPopulationClass
    implicit none
    class           (stellarPopulationBroadBandLuminositiesStandard)                                  , intent(inout) :: self
    integer                                                         , dimension( :                   ), intent(in   ) :: luminosityIndex                       , filterIndex
    double precision                                                , dimension( :                   ), intent(in   ) :: age                                   , redshift
    type            (abundances                                    )                                  , intent(in   ) :: abundancesStellar
    type            (stellarPopulationSpectraPostprocessorList     ), dimension( :                   ), intent(in   ) :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                        )                                  , intent(inout) :: stellarPopulation_
    double precision                                                , dimension(size(luminosityIndex))                :: standardLuminosities
    double precision                                                , dimension(0:1                  )                :: hAge                                  , hMetallicity
    integer         (c_size_t                                      )                                                  :: iAge                                  , iMetallicity   , &
         &                                                                                                               iLuminosityStart                      , iLuminosityEnd , &
         &                                                                                                               populationID
    double precision                                                                                                  :: metallicity

    ! Tabulate the luminosities.
    call self%tabulate(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    ! Obtain a read lock on the luminosity tables.
    call self%luminosityTableLock%setRead()
    ! Get interpolation in metallicity.
    populationID=stellarPopulation_%uniqueID()
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < self%luminosityTables(populationID)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > self%luminosityTables(populationID)%metallicity(self%luminosityTables(populationID)%metallicitiesCount)) then
       iMetallicity=self%luminosityTables(populationID)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       call self%luminosityTables(populationID)%interpolatorMetallicity%linearFactors(metallicity,iMetallicity,hMetallicity)
    end if
    ! Do the interpolations.
    iLuminosityStart=0
    do while (iLuminosityStart < size(luminosityIndex))
       iLuminosityStart=iLuminosityStart+1
       ! Find all contiguous luminosities at this age.
       iLuminosityEnd=iLuminosityStart
       do while (age(iLuminosityEnd) == age(iLuminosityStart) .or. (age(iLuminosityStart) < 0.0d0 .and. age(iLuminosityEnd) < 0.0d0))
          iLuminosityEnd=iLuminosityEnd+1
          if (iLuminosityEnd > size(luminosityIndex)) exit
       end do
       iLuminosityEnd=iLuminosityEnd-1
       ! Only compute luminosities for entries with positive age (negative age implies that the luminosity required is for a
       ! population observed prior to the formation of this population).
       if (age(iLuminosityStart) < 0.0d0) then
          standardLuminosities(iLuminosityStart:iLuminosityEnd)=0.0d0
       else
          ! Check for out of range age.
          if (age(iLuminosityStart) > self%luminosityTables(populationID)%age(self%luminosityTables(populationID)%agesCount)) then
             if (self%maximumAgeExceededIsFatal) then
                call Galacticus_Error_Report('age exceeds the maximum tabulated'//{introspection:location})
             else
                iAge=self%luminosityTables(populationID)%agesCount-1
                hAge=[0.0d0,1.0d0]
             end if
          else
             call self%luminosityTables(populationID)%interpolatorAge%linearFactors(age(iLuminosityStart),iAge,hAge)
          end if
          standardLuminosities(iLuminosityStart:iLuminosityEnd)=                                                                                                 &
               & +self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosityStart:iLuminosityEnd),iAge+0,iMetallicity+0)*hAge(0)*hMetallicity(0) &
               & +self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosityStart:iLuminosityEnd),iAge+0,iMetallicity+1)*hAge(0)*hMetallicity(1) &
               & +self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosityStart:iLuminosityEnd),iAge+1,iMetallicity+0)*hAge(1)*hMetallicity(0) &
               & +self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosityStart:iLuminosityEnd),iAge+1,iMetallicity+1)*hAge(1)*hMetallicity(1)
       end if
       iLuminosityStart=iLuminosityEnd
    end do
    call self%luminosityTableLock%unsetRead()
    ! Prevent interpolation from returning negative fluxes.
    standardLuminosities=max(standardLuminosities,0.0d0)
    return
  end function standardLuminosities

  subroutine standardLuminosityTracks(self,luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,abundancesStellar,redshift,ages,luminosities)
    !% Returns the luminosity for a $1 M_\odot$ simple stellar population of given {\normalfont \ttfamily abundances} drawn from
    !% the given {\normalfont \ttfamily stellarPopulation} and observed through the filter specified by {\normalfont \ttfamily
    !% filterIndex}, for all available ages.
    use            :: Abundances_Structure, only : Abundances_Get_Metallicity, logMetallicityZero, metallicityTypeLogarithmicByMassSolar
    use, intrinsic :: ISO_C_Binding       , only : c_size_t
    use            :: Memory_Management   , only : allocateArray
    implicit none
    class           (stellarPopulationBroadBandLuminositiesStandard), intent(inout) :: self
    integer                                                         , intent(in   ), dimension( :   )              :: luminosityIndex                       , filterIndex
    double precision                                                , intent(in   ), dimension( :   )              :: redshift
    double precision                                                , intent(  out), dimension( :   ), allocatable :: ages
    double precision                                                , intent(  out), dimension( : ,:), allocatable :: luminosities
    type            (stellarPopulationSpectraPostprocessorList     ), intent(in   ), dimension( :   )              :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                        ), intent(inout)                                :: stellarPopulation_
    type            (abundances                                    ), intent(in   )                                :: abundancesStellar
    double precision                                                               , dimension(0:1  )              :: hMetallicity
    integer         (c_size_t                                      )                                               :: iLuminosity                           , iMetallicity, &
         &                                                                                                            jMetallicity                          , populationID
    double precision                                                                                               :: metallicity

    ! Obtain a read lock on the luminosity tables.
    call self%luminosityTableLock%setRead()
    ! Tabulate the luminosities.
    call self%tabulate(luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    ! Get interpolation in metallicity.
    populationID=stellarPopulation_%uniqueID()
    metallicity=Abundances_Get_Metallicity(abundancesStellar,metallicityType=metallicityTypeLogarithmicByMassSolar)
    if (metallicity == logMetallicityZero .or. metallicity < self%luminosityTables(populationID)%metallicity(1)) then
       iMetallicity=1
       hMetallicity=[1.0d0,0.0d0]
    else if (metallicity > self%luminosityTables(populationID)%metallicity(self%luminosityTables(populationID)%metallicitiesCount)) then
       iMetallicity=self%luminosityTables(populationID)%metallicitiesCount-1
       hMetallicity=[0.0d0,1.0d0]
    else
       call self%luminosityTables(populationID)%interpolatorMetallicity%linearFactors(metallicity,iMetallicity,hMetallicity)
    end if
    ! Allocate arrays for ages and luminosities.
    call allocateArray(ages        ,[self%luminosityTables(populationID)%agesCount                      ])
    call allocateArray(luminosities,[self%luminosityTables(populationID)%agesCount,size(luminosityIndex)])
    ! Assign ages.
    ages=self%luminosityTables(populationID)%age
    ! Do the interpolation.
    luminosities(:,:)=0.0d0
    do iLuminosity=1,size(luminosityIndex)
       do jMetallicity=0,1
          luminosities                                          (:,                iLuminosity                             )= &
               & +luminosities                                  (:,                iLuminosity                             )  &
               & +self%luminosityTables(populationID)%luminosity(  luminosityIndex(iLuminosity),:,iMetallicity+jMetallicity)  &
               & *hMetallicity                                  (                                              jMetallicity)
       end do
    end do
    ! Release the read lock on the luminosity tables.
    call self%luminosityTableLock%setRead()
    ! Prevent interpolation from returning negative fluxes.
    luminosities=max(luminosities,0.0d0)
    return
  end subroutine standardLuminosityTracks

  subroutine standardTabulate(self,luminosityIndex,filterIndex,stellarPopulationSpectraPostprocessor_,stellarPopulation_,redshift)
    !% Tabulate stellar population luminosity in the given filters.
    use            :: Abundances_Structure            , only : logMetallicityZero        , metallicityTypeLogarithmicByMassSolar
    use            :: File_Utilities                  , only : File_Exists               , File_Lock                            , File_Unlock               , lockDescriptor
    use            :: Galacticus_Display              , only : Galacticus_Display_Counter, Galacticus_Display_Counter_Clear     , Galacticus_Display_Indent , Galacticus_Display_Unindent, &
         &                                                     verbosityWorking
    use            :: Galacticus_Error                , only : Galacticus_Error_Report   , Galacticus_Warn                      , errorStatusFail           , errorStatusSuccess
    use :: Input_Parameters, only : inputParameters
    use            :: IO_HDF5                         , only : hdf5Access                , hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: ISO_Varying_String              , only : assignment(=)             , char                                 , operator(//)              , var_str
    use            :: Instruments_Filters             , only : Filter_Extent             , Filter_Name
    use            :: Memory_Management               , only : Memory_Usage_Record       , allocateArray                        , deallocateArray
    use            :: Numerical_Constants_Astronomical, only : metallicitySolar
    use            :: Numerical_Integration           , only : integrator                , GSL_Integ_Gauss15
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class           (stellarPopulationBroadBandLuminositiesStandard), intent(inout)                   :: self
    integer                                                         , intent(in   ), dimension(:    ) :: filterIndex                                   , luminosityIndex
    double precision                                                , intent(in   ), dimension(:    ) :: redshift
    type            (stellarPopulationSpectraPostprocessorList     ), intent(in   ), dimension(:    ) :: stellarPopulationSpectraPostprocessor_
    class           (stellarPopulationClass                        ), intent(inout)                   :: stellarPopulation_
    type            (luminosityTable                               ), allocatable  , dimension(:    ) :: luminosityTablesTemporary
    double precision                                                , allocatable  , dimension(:,:,:) :: luminosityTemporary
    logical                                                         , allocatable  , dimension(:    ) :: isTabulatedTemporary
    double precision                                                               , dimension(2    ) :: wavelengthRange
    type            (lockDescriptor                                )                                  :: lockFileDescriptor
    class           (stellarPopulationSpectraClass                 ), pointer                         :: stellarPopulationSpectra_
    class           (stellarPopulationSpectraPostprocessorClass    ), pointer                         :: stellarPopulationSpectraPostprocessorPrevious_
    type            (integrator                                    ), allocatable                     :: integrator_                                   , integratorAB_
    integer         (c_size_t                                      )                                  :: iAge                                          , iLuminosity                          , &
         &                                                                                               iMetallicity                                  , jLuminosity                          , &
         &                                                                                               populationID
    integer                                                                                           :: loopCountMaximum                              , loopCount                            , &
         &                                                                                               errorStatus                                   , luminosityIndexMaximum
    logical                                                                                           :: computeTable                                  , calculateLuminosity                  , &
         &                                                                                               stellarPopulationHashedDescriptorComputed
    double precision                                                                                  :: toleranceRelative                             , normalization
    type            (varying_string                                )                                  :: message                                       , luminositiesFileName                 , &
         &                                                                                               descriptorString                              , stellarPopulationHashedDescriptor    , &
         &                                                                                               postprocessorHashedDescriptor
    character       (len=16                                        )                                  :: datasetName                                   , redshiftLabel                        , &
         &                                                                                               label
    type            (hdf5Object                                    )                                  :: luminositiesFile
    type            (inputParameters                               )                                  :: descriptor

    ! Obtain a read lock on the luminosity tables.
    call self%luminosityTableLock%setRead()
    ! Allocate table storage. First test if the tables must be resized (or allocated).
    populationID=stellarPopulation_%uniqueID()
    if (.not.allocated(self%luminosityTables).or.size(self%luminosityTables) < populationID) then
       ! A resize or allocation is required. Obtain a write lock on the tables and then re-test if resizing is required (as
       ! another thread could potentially have alreaedy done this while we waited for the write lock).
       call self%luminosityTableLock%setWrite(haveReadLock=.true.)
       if (allocated(self%luminosityTables)) then
          if (size(self%luminosityTables) < populationID) then
             call Move_Alloc(self%luminosityTables,luminosityTablesTemporary)
             allocate(self%luminosityTables(populationID))
             self%luminosityTables(1:size(luminosityTablesTemporary))=luminosityTablesTemporary
             self%luminosityTables(size(luminosityTablesTemporary)+1:populationID)%isTabulatedMaximum=0
             deallocate(luminosityTablesTemporary)
             call Memory_Usage_Record(sizeof(self%luminosityTables(1)),blockCount=0)
          end if
       else
          allocate(self%luminosityTables(populationID))
          self%luminosityTables%isTabulatedMaximum=0
          call Memory_Usage_Record(sizeof(self%luminosityTables))
       end if
       call self%luminosityTableLock%unsetWrite(haveReadLock=.true.)
    end if
    ! Determine if we have tabulated luminosities for this luminosityIndex in this population yet.
    luminosityIndexMaximum=maxval(luminosityIndex)
    if (.not.allocated(self%luminosityTables(populationID)%isTabulated) .or. self%luminosityTables(populationID)%isTabulatedMaximum < luminosityIndexMaximum) then
       ! Tabulation is required. Obtain a write lock on the luminosity tables and then retest if tabulation is required (as
       ! another thread may have tabulated while we waited for the write lock).
       call self%luminosityTableLock%setWrite(haveReadLock=.true.)
       luminosityIndexMaximum=maxval(luminosityIndex)
       if (.not.allocated(self%luminosityTables(populationID)%isTabulated) .or. self%luminosityTables(populationID)%isTabulatedMaximum < luminosityIndexMaximum) then
          stellarPopulationSpectra_                      => stellarPopulation_%spectra()
          luminosityIndexMaximum                         =  maxval(luminosityIndex)
          stellarPopulationHashedDescriptorComputed      =  .false.
          stellarPopulationSpectraPostprocessorPrevious_ => null()
          do iLuminosity=1,size(luminosityIndex)
             !$omp critical(broadBandLuminositiesStandardComputeTable)
             if (allocated(self%luminosityTables(populationID)%isTabulated)) then
                if (size(self%luminosityTables(populationID)%isTabulated) >= luminosityIndex(iLuminosity)) then
                   computeTable=.not.self%luminosityTables(populationID)%isTabulated(luminosityIndex(iLuminosity))
                else
                   call Move_Alloc (self%luminosityTables(populationID)%isTabulated,isTabulatedTemporary)
                   call Move_Alloc (self%luminosityTables(populationID)%luminosity ,luminosityTemporary )
                   call allocateArray(self%luminosityTables(populationID)%isTabulated,[luminosityIndexMaximum])
                   call allocateArray(self%luminosityTables(populationID)%luminosity ,[luminosityIndexMaximum&
                        &,self%luminosityTables(populationID)%agesCount,self%luminosityTables(populationID)%metallicitiesCount])
                   self%luminosityTables(populationID)%isTabulated(1:size(isTabulatedTemporary)    )=isTabulatedTemporary
                   self%luminosityTables(populationID)%isTabulated(  size(isTabulatedTemporary)+1:luminosityIndexMaximum)=.false.
                   self%luminosityTables(populationID)%luminosity (1:size(isTabulatedTemporary),:,:)=luminosityTemporary
                   call deallocateArray(isTabulatedTemporary)
                   call deallocateArray(luminosityTemporary)
                   computeTable=.true.
                end if
             else
                call allocateArray(self%luminosityTables(populationID)%isTabulated,[luminosityIndexMaximum])
                self%luminosityTables(populationID)%isTabulated=.false.
                ! Since we have not yet tabulated any luminosities yet for this population, we need to get a list of suitable
                ! metallicities and ages at which to tabulate.
                call stellarPopulationSpectra_%tabulation(self%luminosityTables(populationID)%agesCount,self%luminosityTables(populationID)%metallicitiesCount,self%luminosityTables(populationID)%age,self%luminosityTables(populationID)%metallicity)
                where (self%luminosityTables(populationID)%metallicity > 0.0d0)
                   self%luminosityTables(populationID)%metallicity=log10(self%luminosityTables(populationID)%metallicity/metallicitySolar)
                elsewhere
                   self%luminosityTables(populationID)%metallicity=logMetallicityZero
                end where
                call allocateArray(self%luminosityTables(populationID)%luminosity,[luminosityIndexMaximum,self%luminosityTables(populationID)%agesCount ,self%luminosityTables(populationID)%metallicitiesCount])
                self%luminosityTables(populationID)%interpolatorAge        =interpolator(self%luminosityTables(populationID)%age        )
                self%luminosityTables(populationID)%interpolatorMetallicity=interpolator(self%luminosityTables(populationID)%metallicity)
                computeTable=.true.
             end if
             ! If we haven't, do so now.
             if (computeTable) then
                ! Determine if we can read the required luminosity from file.
                calculateLuminosity=.true.
                if (self%storeToFile) then
                   ! Construct name of the file to which this would be stored.
                   if (.not.stellarPopulationHashedDescriptorComputed) then
                      stellarPopulationHashedDescriptor        =stellarPopulation_%hashedDescriptor(includeSourceDigest=.true.)
                      stellarPopulationHashedDescriptorComputed=.true.
                   end if
                   if (.not.associated(stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_,stellarPopulationSpectraPostprocessorPrevious_)) then
                      postprocessorHashedDescriptor                  =  stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_%hashedDescriptor(includeSourceDigest=.true.)
                      stellarPopulationSpectraPostprocessorPrevious_ => stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_
                   end if
                   luminositiesFileName=self%storeDirectory                  // &
                        &               "/stellarLuminosities::filter:"      // &
                        &               Filter_Name(filterIndex(iLuminosity))// &
                        &               "::postprocessor:"                   // &
                        &               postprocessorHashedDescriptor        // &
                        &               "::population:"                      // &
                        &               stellarPopulationHashedDescriptor    // &
                        &               ".hdf5"
                   if (File_Exists(luminositiesFileName)) then
                      ! Construct the dataset name.
                      write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                      datasetName="redshift"//adjustl(trim(redshiftLabel))
                      ! Open the file and check for the required dataset.
                      ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
                      call File_Lock(char(luminositiesFileName),lockFileDescriptor,lockIsShared=.true.)
                      !$ call hdf5Access%set()
                      call luminositiesFile%openFile(char(luminositiesFileName),readOnly=.true.)
                      if (luminositiesFile%hasDataset(trim(datasetName))) then
                         ! Read the dataset.
                         call luminositiesFile%readDatasetStatic(trim(datasetName),self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:))
                         ! We do not need to calculate this luminosity.
                         calculateLuminosity=.false.
                      end if
                      call luminositiesFile%close()
                      !$ call hdf5Access%unset()
                      call File_Unlock(lockFileDescriptor)
                   end if
                end if
                ! Compute the luminosity if necessary.
                if (calculateLuminosity) then
                   ! Display a message and counter.
                   message=var_str('Tabulating stellar luminosities for stellar population #')//populationID//', luminosity '
                   write (redshiftLabel,'(f6.3)') redshift(iLuminosity)
                   message=message                                               // &
                        &  Filter_Name                 (filterIndex(iLuminosity))// &
                        &  ":z"                                                  // &
                        &  trim(adjustl(redshiftLabel))                          // &
                        &  " "                                                   // &
                        &                                           iLuminosity  // &
                        &  " of "                                                // &
                        &  size(luminosityIndex)
                   call Galacticus_Display_Indent (message,verbosityWorking)
                   call Galacticus_Display_Counter(0,.true.,verbosityWorking)
                   ! Get wavelength extent of the filter.
                   wavelengthRange=Filter_Extent(filterIndex(iLuminosity))
                   ! Integrate over the wavelength range.
                   standardFilterIndex = filterIndex     (iLuminosity)
                   standardRedshift    = redshift        (iLuminosity)
                   standardPopulationID= populationID
                   loopCountMaximum    =+self%luminosityTables(populationID)%metallicitiesCount &
                        &               *self%luminosityTables(populationID)%agesCount
                   loopCount           = 0
                   !$omp parallel private(iAge,iMetallicity,integrator_,toleranceRelative,errorStatus) copyin(standardFilterIndex,standardRedshift,standardPopulationID)
                   allocate(integrator_)
                   integrator_=integrator(integrandFilteredLuminosity,toleranceRelative=self%integrationToleranceRelative,integrationRule=GSL_Integ_Gauss15,intervalsMaximum=10000_c_size_t)
                   allocate(standardStellarPopulationSpectra             ,mold=stellarPopulationSpectra_                                                                 )
                   allocate(standardStellarPopulationSpectraPostprocessor,mold=stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_)
                   !$omp critical(broadBandLuminositiesDeepCopy)
                   !# <deepCopyReset variables="stellarPopulationSpectra_ stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_"/>
                   !# <deepCopy source="stellarPopulationSpectra_"                                                                  destination="standardStellarPopulationSpectra"             />
                   !# <deepCopy source="stellarPopulationSpectraPostprocessor_(iLuminosity)%stellarPopulationSpectraPostprocessor_" destination="standardStellarPopulationSpectraPostprocessor"/>
                   !$omp end critical(broadBandLuminositiesDeepCopy)
                   !$omp do
                   do iAge=1,self%luminosityTables(populationID)%agesCount
                      standardAge=self%luminosityTables(populationID)%age(iAge)
                      do iMetallicity=1,self%luminosityTables(populationID)%metallicitiesCount
                         ! Update the counter.
                         !$omp atomic
                         loopCount=loopCount+1
                         call Galacticus_Display_Counter(int(100.0d0*dble(loopCount)/dble(loopCountMaximum)),.false.,verbosityWorking)
                         call standardAbundances%metallicitySet(self%luminosityTables(populationID)%metallicity(iMetallicity) &
                              &,metallicityType=metallicityTypeLogarithmicByMassSolar)
                         toleranceRelative=self%integrationToleranceRelative
                         errorStatus      =errorStatusFail
                         do while (errorStatus /= errorStatusSuccess)
                            call integrator_%toleranceSet(toleranceRelative=toleranceRelative)
                            self%luminosityTables(populationID)%luminosity(                              &
                                 &                                         luminosityIndex(iLuminosity), &
                                 &                                         iAge                        , &
                                 &                                         iMetallicity                  &
                                 &                                        )                              &
                                 & =integrator_%integrate(                           &
                                 &                               wavelengthRange(1), &
                                 &                               wavelengthRange(2), &
                                 &                        status=errorStatus         &
                                 &                       )
                            if (errorStatus /= errorStatusSuccess) then
                               if (self%integrationToleranceDegrade.and.toleranceRelative < 1.0d0) then
                                  toleranceRelative=2.0d0*toleranceRelative
                                  write (label,'(e9.3)') 2.0d0*self%integrationToleranceRelative
                                  message=         "WARNING: increasing relative tolerance for stellar population luminosities to"    //char(10)
                                  message=message//trim(adjustl(label))//" and retrying integral"//{introspection:location}
                                  call Galacticus_Warn(message)
                               else if (self%integrationToleranceDegrade) then
                                  message="integration of stellar populations failed"
                                  call Galacticus_Error_Report(message//{introspection:location})
                               else
                                  write (label,'(e9.3)') 2.0d0*self%integrationToleranceRelative
                                  message=         "integration of stellar populations failed"                                        //char(10)
                                  message=message//"HELP: consider increasing the [self%integrationToleranceRelative]"                //char(10)
                                  message=message//"      parameter to "//trim(adjustl(label))//" to reduce the integration tolerance"//char(10)
                                  message=message//"      required if you can accept this lower accuracy."
                                  call Galacticus_Error_Report(message//{introspection:location})
                               end if
                            end if
                         end do
                      end do
                   end do
                   !$omp end do
                   deallocate(integrator_)
                   !# <objectDestructor name="standardStellarPopulationSpectra"             />
                   !# <objectDestructor name="standardStellarPopulationSpectraPostprocessor"/>
                   !$omp end parallel
                   ! Clear the counter and write a completion message.
                   call Galacticus_Display_Counter_Clear(           verbosityWorking)
                   call Galacticus_Display_Unindent     ('finished',verbosityWorking)
                   ! Get the normalization by integrating a zeroth magnitude (AB) source through the filter.
                   allocate(integratorAB_)
                   integratorAB_=integrator(integrandFilteredLuminosityAB,toleranceRelative=self%integrationToleranceRelative)
                   normalization=integratorAB_%integrate(wavelengthRange(1),wavelengthRange(2))
                   deallocate(integratorAB_)
                   ! Normalize the luminosity.
                   self         %luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:) &
                        & =+self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:) &
                        &  /normalization
                   ! Store the luminosities to file.
                   if (self%storeToFile) then
                      ! Construct the dataset name.
                      write (redshiftLabel,'(f7.4)') redshift(iLuminosity)
                      datasetName="redshift"//adjustl(trim(redshiftLabel))
                      ! Open the file.
                      descriptor=inputParameters()
                      call stellarPopulation_%descriptor(descriptor,includeMethod=.true.)
                      descriptorString=descriptor%serializeToString()
                      call descriptor%destroy()
                      ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
                      call File_Lock(char(luminositiesFileName),lockFileDescriptor,lockIsShared=.false.)
                      !$ call hdf5Access%set()
                      call luminositiesFile%openFile      (char(luminositiesFileName)             )
                      if (.not.luminositiesFile%hasAttribute('parameters')) call luminositiesFile%writeAttribute(char(descriptorString),'parameters')
                      ! Write the dataset.
                      if (.not.luminositiesFile%hasDataset(trim(datasetName))) &
                           & call luminositiesFile%writeDataset(self%luminosityTables(populationID)%luminosity(luminosityIndex(iLuminosity),:,:),datasetName=trim(datasetName),commentText="Tabulated luminosities at redshift z="//adjustl(trim(redshiftLabel)))
                      ! Close the file.
                      call luminositiesFile%close()
                      !$ call hdf5Access%unset()
                      call File_Unlock(lockFileDescriptor)
                   end if
                end if
                ! Flag that calculations have been performed for this filter.
                self%luminosityTables(populationID)%isTabulated(luminosityIndex(iLuminosity))=.true.
                if (luminosityIndex(iLuminosity) > self%luminosityTables(populationID)%isTabulatedMaximum) then
                   jLuminosity=self%luminosityTables(populationID)%isTabulatedMaximum
                   do while (jLuminosity < size(self%luminosityTables(populationID)%isTabulated) .and. self%luminosityTables(populationID)%isTabulated(min(jLuminosity+1,size(self%luminosityTables(populationID)%isTabulated))))
                      jLuminosity=jLuminosity+1
                   end do
                   self%luminosityTables(populationID)%isTabulatedMaximum=jLuminosity
                end if
             end if
             !$omp end critical(broadBandLuminositiesStandardComputeTable)
          end do
       end if
       call self%luminosityTableLock%unsetWrite(haveReadLock=.true.)
    end if
    ! Release our read lock on the luminosity tables.
    call self%luminosityTableLock%unsetRead()
    return

  contains

    double precision function integrandFilteredLuminosity(wavelength)
      !% Integrand for the luminosity through a given filter.
      use :: Instruments_Filters, only : Filter_Response
      implicit none
      double precision, intent(in   ) :: wavelength
      double precision                :: wavelengthRedshifted

      ! If this luminosity is for a redshifted spectrum, then we shift wavelength at which we sample the stellar population
      ! spectrum to be a factor of (1+z) smaller. We therefore integrate over the stellar SED at shorter wavelengths, since these
      ! will be shifted into the filter by z=0. Factor of 1/wavelength appears since we want to integrate F_ν (dν / ν) and dν =
      ! -c/λ² dλ. Note that we follow the convention of Hogg et al. (2002) and assume that the filter response gives the fraction
      ! of incident photons received by the detector at a given wavelength, multiplied by the relative photon response (which will
      ! be 1 for a photon-counting detector such as a CCD, or proportional to the photon energy for a bolometer/calorimeter type
      ! detector).
      wavelengthRedshifted=wavelength/(1.0d0+standardRedshift)
      integrandFilteredLuminosity=+Filter_Response                                         (standardFilterIndex             ,wavelength          ) &
           &                      *standardStellarPopulationSpectra             %luminosity(standardAbundances  ,standardAge,wavelengthRedshifted) &
           &                      *standardStellarPopulationSpectraPostprocessor%multiplier(wavelengthRedshifted,standardAge,standardRedshift    ) &
           &                      /                                                                                          wavelength
      return
    end function integrandFilteredLuminosity

    double precision function integrandFilteredLuminosityAB(wavelength)
      !% Integrand for the luminosity of a zeroth magnitude (AB) source through a given filter.
      use :: Instruments_Filters             , only : Filter_Response
      use :: Numerical_Constants_Astronomical, only : luminositySolar, luminosityZeroPointAB
      implicit none
      double precision, intent(in   ) :: wavelength
      ! Luminosity of a zeroth magintude (AB) source in Solar luminosities per Hz.
      double precision, parameter     :: luminosityZeroPointABSolar=luminosityZeroPointAB/luminositySolar

      integrandFilteredLuminosityAB=Filter_Response(standardFilterIndex,wavelength)*luminosityZeroPointABSolar/wavelength
      return
    end function integrandFilteredLuminosityAB

  end subroutine standardTabulate


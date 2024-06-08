!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  !+    Contributions to this file made by: Sachi Weerasooriya

  !!{
  Contains a module which implements a property extractor class for the emission line luminosity of a component.
  !!}
  use :: Galactic_Structure_Options     , only : enumerationComponentTypeType
  use :: Output_Times                   , only : outputTimesClass
  use :: Star_Formation_Histories       , only : starFormationHistoryClass
  use :: HII_Region_Luminosity_Functions, only : hiiRegionLuminosityFunctionClass
  use :: Star_Formation_Histories       , only : starFormationHistoryClass
  
  type:: emissionLineLuminosityTemplate
     !!{
     Type used to store luminosity templates for emission lines.
     !!}
     private
     integer         (c_size_t)                                :: countLines            =-1_c_size_t
     double precision          , allocatable, dimension(:,:,:) :: emissionLineLuminosity
  end type emissionLineLuminosityTemplate

  !![
  <nodePropertyExtractor name="nodePropertyExtractorLuminosityEmissionLine">
    <description>An emission line luminosity property extractor class. The luminosity of the named emission line (given by the {\normalfont
      \ttfamily lineNames} parameter: if multiple lines are named, the sum of their luminosities) is computed.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLuminosityEmissionLine
     !!{
     A property extractor class for the emission line luminosity of a component.
     !!}
     private
     class           (starFormationHistoryClass       ), pointer                       :: starFormationHistory_                => null()
     class           (outputTimesClass                ), pointer                       :: outputTimes_                         => null()
     class           (hiiRegionLuminosityFunctionClass), pointer                       :: hiiRegionLuminosityFunction_         => null()
     type            (enumerationComponentTypeType    )                                :: component
     integer                                                                           :: countWavelengths                             , countLines
     type            (varying_string                  ), allocatable, dimension(:    ) :: lineNames                                    , names_                      , &
          &                                                                               descriptions_
     double precision                                  , allocatable, dimension(:    ) :: metallicityBoundaries                        , metallicities               , &
          &                                                                               ages
     double precision                                  , allocatable, dimension(:,:,:) :: luminositiesReduced
     double precision                                  , allocatable, dimension(:,:  ) :: ionizingLuminosityHydrogenNormalized
     type            (emissionLineLuminosityTemplate  ), allocatable, dimension(:    ) :: templates
     double precision                                                                  :: metallicityPopulationMinimum                 , metallicityPopulationMaximum, &
          &                                                                               agePopulationMaximum                         , resolution                  , &
          &                                                                               factorWavelength                             , toleranceRelative           , &
          &                                                                               ionizingLuminosityHydrogenMean
     logical                                                                           :: useluminosityTemplates
   contains
     !![
     <methods>
       <method description="Return a hashed descriptor of the object which incorporates the time and metallicity binning of the star formation history." method="historyHashedDescriptor"/>
       <method description="Compute the mean luminosity of the stellar population in the given bin of the star formation history."                       method="luminosityMean"         />
       <method description="Return the index of the template time to use."                                                                               method="indexTemplateTime"      />
       <method description="Return the index of the template luminosities to use."                                                                       method="indexTemplateNode"      />
     </methods>
     !!]
     final     ::                            emissionLineLuminosityDestructor
     procedure :: historyHashedDescriptor => emissionLineLuminosityHistoryHashedDescriptor
     procedure :: elementCount            => emissionLineLuminosityElementCount
     procedure :: extract                 => emissionLineLuminosityExtract
     procedure :: names                   => emissionLineLuminosityNames
     procedure :: descriptions            => emissionLineLuminosityDescriptions
     procedure :: unitsInSI               => emissionLineLuminosityUnitsInSI
     procedure :: luminosityMean          => emissionLineLuminosityMean
     procedure :: indexTemplateTime       => emissionLineLuminosityIndexTemplateTime
     procedure :: indexTemplateNode       => emissionLineLuminosityIndexTemplateNode 
  end type nodePropertyExtractorLuminosityEmissionLine
  
  interface nodePropertyExtractorLuminosityEmissionLine
     !!{
     Constructors for the ``emissionLineLuminosity'' output analysis class.
     !!}
     module procedure emissionLineLuminosityConstructorParameters
     module procedure emissionLineLuminosityConstructorInternal
  end interface nodePropertyExtractorLuminosityEmissionLine
      
contains

  function emissionLineLuminosityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily emission line luminosity} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type            (nodePropertyExtractorLuminosityEmissionLine)                              :: self
    type            (inputParameters                            ), intent(inout)               :: parameters
    class           (starFormationHistoryClass                  ), pointer                     :: starFormationHistory_
    class           (outputTimesClass                           ), pointer                     :: outputTimes_
    class           (hiiRegionLuminosityFunctionClass           ), pointer                     :: hiiRegionLuminosityFunction_
    type            (varying_string                             ), allocatable  , dimension(:) :: lineNames
    type            (varying_string                             )                              :: component
    double precision                                                                           :: toleranceRelative
    
    allocate(lineNames(parameters%count('lineNames')))
    !![
    <inputParameter>
      <name>lineNames</name>
      <source>parameters</source>
      <description>The emission lines to extract.</description>
    </inputParameter>
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation rate.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelative</name>
      <source>parameters</source>
      <defaultValue>1.0d-3</defaultValue>
      <description>The relative tolerance used in integration over stellar population spectra.</description>
    </inputParameter>
    <objectBuilder class="starFormationHistory"        name="starFormationHistory_"        source="parameters"/>
    <objectBuilder class="outputTimes"                 name="outputTimes_"                 source="parameters"/>
    <objectBuilder class="hiiRegionLuminosityFunction" name="hiiRegionLuminosityFunction_" source="parameters"/>
    !!]
    self=nodePropertyExtractorLuminosityEmissionLine(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),lineNames,toleranceRelative,starFormationHistory_,outputTimes_,hiiRegionLuminosityFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"       />
    <objectDestructor name="outputTimes_"                />
    <objectDestructor name="hiiRegionLuminosityFunction_"/>
    !!]
    return
  end function emissionLineLuminosityConstructorParameters

  function emissionLineLuminosityConstructorInternal(component,lineNames,toleranceRelative,starFormationHistory_,outputTimes_,hiiRegionLuminosityFunction_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sed} property extractor class.
    !!}
    use :: Galactic_Structure_Options      , only : componentTypeDisk, componentTypeSpheroid
    use :: Galacticus_Nodes                , only : nodeComponentDisk, nodeComponentSpheroid
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Error                           , only : Error_Report
    use :: Input_Paths                     , only : inputPath        , pathTypeDataStatic
    implicit none
    type            (nodePropertyExtractorLuminosityEmissionLine      )                             :: self
    type            (enumerationComponentTypeType              ), intent(in   )                     :: component
    type            (varying_string                            ), intent(in   ), dimension(:      ) :: lineNames
    class           (starFormationHistoryClass                 ), intent(in   ), target             :: starFormationHistory_
    class           (outputTimesClass                          ), intent(in   ), target             :: outputTimes_
    class           (hiiRegionLuminosityFunctionClass          ), intent(in   ), target             :: hiiRegionLuminosityFunction_
    double precision                                            , intent(in   )                     :: toleranceRelative    
    double precision                                            ,                                   :: deltaIonizingFluxHydrogen         , rateHydrogenIonizingPhotonsMinimum, &
         &                                                                                             rateHydrogenIonizingPhotonsMaximum 
    double precision                                            , allocatable  , dimension(:      ) :: ionizingFluxHydrogen
    double precision                                            , allocatable  , dimension(:,:,:,:) :: luminosities
    type            (hdf5Object                                )                                    :: emissionLinesFile                 , lines
    integer                                                                                         :: i
    !![
    <constructorAssign variables="lineNames, component, toleranceRelative, *starFormationHistory_, *outputTimes_,*hiiRegionLuminosityFunction_"/>
    !!]
    
    if     (                                                                                                    &
         &   component /= componentTypeDisk                                                                     &
         &  .and.                                                                                               &
         &   component /= componentTypeSpheroid                                                                 &
         & ) call Error_Report("only 'disk' and 'spheroid' components are supported"//{introspection:location})
    ! Get details of the star formation rate tabulation.
    self%metallicityBoundaries =self%starFormationHistory_%metallicityBoundaries      ()
    self%useluminosityTemplates=self%starFormationHistory_%perOutputTabulationIsStatic()
    ! Read the table of emission line luminosities.
    !$ call hdf5Access%set()
    call emissionLinesFile%openFile(char(inputPath(pathTypeDataStatic))//"hiiRegions/cloudyTableBC2003.hdf5",readOnly=.true.)
    lines=emissionLinesFile%openGroup('lines')
    do i=1,size(lineNames)
       if (.not.lines%hasDataset(char(self%lineNames(i)))) call Error_Report('line "'//char(self%lineNames(i))//'" not found'//{introspection:location})
    end do
    call emissionLinesFile%readDataset('metallicity'                         ,self%metallicities                       )
    call emissionLinesFile%readDataset('age'                                 ,self%ages                                )
    call emissionLinesFile%readDataset('ionizingLuminosityHydrogen'          ,     ionizingFluxHydrogen                )
    call emissionLinesFile%readDataset('ionizingLuminosityHydrogenNormalized',self%ionizingLuminosityHydrogenNormalized)
    self%metallicityPopulationMinimum=minval(self%metallicities)
    self%metallicityPopulationMaximum=maxval(self%metallicities)
    self%agePopulationMaximum        =maxval(self%ages         )
    allocate(                                  &
         &        luminosities                 &
         &   (                                 &
         &    size(self%ages                ), &
         &    size(self%metallicities       ), &
         &    size(     ionizingFluxHydrogen), &
         &    size(     lineNames           )  &
         &   )                                 &
         &  )
    allocate(                                  &
         &   self%luminositiesReduced          &
         &   (                                 &
         &    size(self%ages                ), &
         &    size(self%metallicities       ), &
         &    size(     lineNames           )  &
         &   )                                 &
         &  )
    self%luminositiesReduced=0.0d0
    do i=1,size(lineNames)
       call lines%readDatasetStatic(char(lineNames(i)),luminosities(:,:,:,i))      
    end do
    call lines            %close()
    call emissionLinesFile%close()
    !$ call hdf5Access%unset()
    ! Calculate emission line luminosities as a function of age and metallicity by averaging over the distribution of HII region
    ! luminosities.
    deltaIonizingFluxHydrogen=+ionizingFluxHydrogen(2) &
         &                    /ionizingFluxHydrogen(1)
    do i=1,size(ionizingFluxHydrogen)
       rateHydrogenIonizingPhotonsMinimum=ionizingFluxHydrogen(i)/sqrt(deltaIonizingFluxHydrogen) 
       rateHydrogenIonizingPhotonsMaximum=ionizingFluxHydrogen(i)*sqrt(deltaIonizingFluxHydrogen)
       ! Accumulate the luminosity weighted by the cumulative fraction of HII regions in this luminosity interval.
       self%luminositiesReduced=+self%luminositiesReduced                                                                                                                &
            &                   +self%hiiRegionLuminosityFunction_%cumulativeDistributionFunction(rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum) &  
            &                   *                                  luminosities                  (:,:,i,:)
    end do
    ! Normalize reduced luminosities to the total fraction of HII regions in the luminosity interval spanned by the table. Also,
    ! find the mean ionizing luminosity of HII regions in this luminosity interval.
    rateHydrogenIonizingPhotonsMinimum =+ionizingFluxHydrogen(size(ionizingFluxHydrogen))/sqrt(deltaIonizingFluxHydrogen) 
    rateHydrogenIonizingPhotonsMaximum =+ionizingFluxHydrogen(size(ionizingFluxHydrogen))*sqrt(deltaIonizingFluxHydrogen)
    self%luminositiesReduced           =+self%luminositiesReduced                                                                                                                &
         &                              /self%hiiRegionLuminosityFunction_%cumulativeDistributionFunction(rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum)
    self%ionizingLuminosityHydrogenMean=+self%hiiRegionLuminosityFunction_%cumulativeLuminosity          (rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum)
    ! Construct property names and descriptions.
    allocate(self%names_       (size(lineNames)))
    allocate(self%descriptions_(size(lineNames)))
    do i=1,size(lineNames)
       select case (self%component%ID)
       case (componentTypeDisk    %ID)
          self%names_       (i)="luminosityEmissionLineDisk:"    //lineNames(i)
       case (componentTypeSpheroid%ID)
          self%names_       (i)="luminosityEmissionLineSpheroid:"//lineNames(i)
       end select
       self%descriptions_(i)="Luminosity of the "                //lineNames(i)//" emission line [ergs/s]"
    end do
    self%countLines=size(lineNames)
    return    
  end function emissionLineLuminosityConstructorInternal

  subroutine emissionLineLuminosityDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily emission line luminosity} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationHistory_"       />
    <objectDestructor name="self%outputTimes_"                />
    <objectDestructor name="self%hiiRegionLuminosityFunction_"/>
    !!]
    return
  end subroutine emissionLineLuminosityDestructor
  
  integer function emissionLineLuminosityElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily emissionLineLuminosity} property extractors.
    !!}
    implicit none
    class     (nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    emissionLineLuminosityElementCount=self%countLines
    return
  end function emissionLineLuminosityElementCount

  function emissionLineLuminosityExtract(self,node,time,instance) result(luminosity)
    !!{
    Implement a {\normalfont \ttfamily luminosityEmissionLine} property extractor.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentDisk, nodeComponentSpheroid
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid
    use :: Histories                 , only : history
    implicit none
    double precision                                             , dimension(:  )            , allocatable :: luminosity
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)   , target                :: self
    type            (treeNode                                   ), intent(inout)   , target                :: node
    double precision                                             , intent(in   )                           :: time
    type            (multiCounter                               ), intent(inout)   , optional              :: instance
    class           (nodeComponentDisk                          )                  , pointer               :: disk
    class           (nodeComponentSpheroid                      )                  , pointer               :: spheroid
    double precision                                             , dimension(:,:,:), pointer               :: luminosityTemplate_
    double precision                                             , dimension(:,:,:), target  , allocatable :: luminosityTemplate
    type            (history                                    )                                          :: starFormationHistory
    integer                                                                                                :: indexTemplate       , iLines
    !$GLC attributes unused :: instance

    allocate(luminosity(self%countLines))
    luminosity=0.0d0
    ! Get the relevant star formation history.
    select case (self%component%ID)
    case (componentTypeDisk    %ID)
       disk                 => node    %disk                ()
       starFormationHistory =  disk    %starFormationHistory()
    case (componentTypeSpheroid%ID)
       spheroid             => node    %spheroid            ()
       starFormationHistory =  spheroid%starFormationHistory()
    end select
    if (.not.starFormationHistory%exists()) return
    ! Get the index of the template to use.
    indexTemplate=self%indexTemplateNode(node,starFormationHistory)
    if (indexTemplate > 0) then
       ! Stored templates can be used, so point to the relevant set.
       luminosityTemplate_ => self%templates(indexTemplate)%emissionLineLuminosity
    else
       ! Stored templates can not be used, get the templates for this specific case, and point to them.
       luminosityTemplate  =  self%luminosityMean(time,starFormationHistory)
       luminosityTemplate_ => luminosityTemplate
    end if
    do iLines=1,size(luminosity,dim=1)
       luminosity(iLines)=sum(luminosityTemplate_(iLines,:,:)*starFormationHistory%data(:,:))
    end do
    return
  end function emissionLineLuminosityExtract

  subroutine emissionLineLuminosityNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily emissionLines}.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)                            :: self
    double precision                                             , intent(in   )                            :: time
    type            (varying_string                             ), intent(inout), dimension(:), allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%countLines))
    names=self%names_
    return
  end subroutine emissionLineLuminosityNames
  
  subroutine emissionLineLuminosityDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily emission line luminosity} property.
    !!}
    implicit none
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)                            :: self
    double precision                                             , intent(in   )                            :: time
    type            (varying_string                             ), intent(inout), dimension(:), allocatable :: descriptions
    !$GLC attributes unused :: self, time
    
    allocate(descriptions(self%countLines))
    descriptions=self%descriptions_
    return
  end subroutine emissionLineLuminosityDescriptions
  
  function emissionLineLuminosityUnitsInSI(self,time) result(unitsInSI)
  !!{
    Return the units of the {\normalfont \ttfamily emissionLineLuminosity} properties in the SI system.
    !!}
    use :: Numerical_Constants_Units, only : ergs
    implicit none
    double precision                                             , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)               :: self
    double precision                                             , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(unitsInSI(self%countLines))
    unitsInSI=ergs
    return
  end function emissionLineLuminosityUnitsInSI

  integer function emissionLineLuminosityIndexTemplateTime(self,time) result(index)
    !!{
    Find the index of the template emission lines to use.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self
    double precision                                             , intent(in   ) :: time
    integer         (c_size_t                                   )                :: indexOutput

    index=-1
    if (.not.self%useluminosityTemplates) return
    ! Check that the time is an output time.
    indexOutput=self%outputTimes_%index(time,findClosest=.true.)
    if (.not.Values_Agree(time,self%outputTimes_%time(indexOutput),relTol=1.0d-6)) return
    index=int(indexOutput)
    return
  end function emissionLineLuminosityIndexTemplateTime

  integer function emissionLineLuminosityIndexTemplateNode(self,node,starFormationHistory) result(index)
    !!{
    Find the index of the template emission line luminosity to use, and also compute the template.
    !!}
    use :: Display             , only : displayMessage    , verbosityLevelWorking
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Histories           , only : history
    use :: ISO_Varying_String  , only : var_str
    use :: HDF5_Access         , only : hdf5Access
    use :: IO_HDF5             , only : hdf5Object
    use :: Numerical_Comparison, only : Values_Agree
    use :: File_Utilities      , only : File_Exists       , File_Lock            , File_Unlock, lockDescriptor
    use :: String_Handling     , only : operator(//)
    use :: Input_Paths         , only : inputPath         , pathTypeDataDynamic
    implicit none
    class    (nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self
    type     (treeNode                                   ), intent(inout) :: node
    type     (history                                    ), intent(in   ) :: starFormationHistory
    class    (nodeComponentBasic                         ), pointer       :: basic
    integer  (c_size_t                                   )                :: indexOutput
    integer                                                               :: i
    type     (lockDescriptor                             )                :: fileLock
    type     (hdf5Object                                 )                :: file
    type     (varying_string                             )                :: fileName
    character(len=16                                     )                :: label
    
    ! Return a negative index if templates are not being used.
    index=-1
    if (.not.self%useluminosityTemplates) return
    ! Check that the node exists at at output time.
    basic       => node             %basic(                               )
    indexOutput =  self%outputTimes_%index(basic%time(),findClosest=.true.)
    if (.not.Values_Agree(basic%time(),self%outputTimes_%time(indexOutput),relTol=1.0d-6)) return
    index=int(indexOutput)
    do i=1,(self%countLines)
       ! Ensure that the templates have been built for this index.
       if (.not.allocated(self%templates)) allocate(self%templates(self%outputTimes_%count()))
       if (.not.allocated(self%templates(index)%emissionLineLuminosity)) then
          ! Construct the file name.
          fileName=inputPath(pathTypeDataDynamic)                         // &
               &        'stellarPopulations/'                             // &
               &        self%objectType             (                    )// &
               &        '_'                                               // &
               &        self%historyHashedDescriptor(starFormationHistory)// &
               &        '_'                                               // &
               &        self%lineNames(i)                                 // &
               &        '_'                                               // &
               &        indexOutput                                       // &
               &        '.hdf5'
          ! Check if the templates can be retrieved from file.
          !! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
          call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
          if (File_Exists(fileName)) then
             !$ call hdf5Access%set()
             call file%openFile(char(fileName))
             if (file%hasDataset('luminosityTemplate')) then
                !$omp critical(gfortranInternalIO)
                write (label,'(f12.8)') self%outputTimes_%time(indexOutput)
                !$omp end critical(gfortranInternalIO)
                call displayMessage("reading SED tabulation for time "//trim(adjustl(label))//" Gyr from file '"//fileName//"'",verbosityLevelWorking)
                call file%readDataset('luminosityTemplate',self%templates(index)%emissionLineLuminosity)
             end if
             call file%close()
             !$ call hdf5Access%unset()
          end if
          if (.not.allocated(self%templates(index)%emissionLineLuminosity)) then
             self%templates(index)%emissionLineLuminosity=self%luminosityMean(self%outputTimes_%time(indexOutput),starFormationHistory,parallelize=.true.)
             !$omp critical(gfortranInternalIO)
             write (label,'(f12.8)') self%outputTimes_%time(indexOutput)
             !$omp end critical(gfortranInternalIO)
             call displayMessage("storing Emission Line Luminosity tabulation for time "//trim(adjustl(label))//" Gyr to file '"//fileName//"'",verbosityLevelWorking)
             !$ call hdf5Access%set()
             call file%openFile(char(fileName),overWrite=.false.,readOnly=.false.)
             call file%writeDataset(self%templates(index)%emissionLineLuminosity,'luminosityTemplate')
             call file%close()
             !$ call hdf5Access%unset()
          end if
          call File_Unlock(fileLock)
       end if
    end do
    return
  end function emissionLineLuminosityIndexTemplateNode

  function emissionLineLuminosityMean(self,time,starFormationHistory,parallelize) result(luminosityMean)
    !!{
    Compute the mean luminosity of the stellar population in each bin of the star formation history.
    !!}
    use :: Display                , only : displayIndent        , displayUnindent, displayCounter, displayCounterClear, &
         &                                 verbosityLevelWorking
    use :: Error                  , only : Error_Report
    use :: Histories              , only : history
    use :: Numerical_Integration  , only : integrator
    use :: Multi_Counters         , only : multiCounter
    use :: Locks                  , only : ompLock
    use :: Numerical_Interpolation, only : interpolator
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num
    implicit none
    double precision                                             , dimension(:,:,:), allocatable :: luminosityMean
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)                 :: self
    double precision                                             , intent(in   )                 :: time
    type            (history                                    ), intent(in   )                 :: starFormationHistory
    logical                                                      , intent(in   )   , optional    :: parallelize
    type            (integrator                                 ), save            , allocatable :: integratorTime            , integratorMetallicity
    type            (interpolator                               ), save            , allocatable :: interpolatorTime          , interpolatorMetallicity
    integer         (c_size_t                                   )                                :: iTime                     , iMetallicity           , &
         &                                                                                          counter                   , counterMaximum         , &
         &                                                                                          iterator
    integer         (c_size_t                                   ), save                          :: iLine
    double precision                                                                             :: metallicityMinimum        , metallicityMaximum
    double precision                                             , save                          :: timeMinimum               , timeMaximum            , &
         &                                                                                          age                       , metallicity_
    character       (len=12                                     )                                :: label
    type            (multiCounter                               )                                :: state
    type            (ompLock                                    )                                :: stateLock
    !$omp threadprivate(iLine,integratorTime,integratorMetallicity,interpolatorTime,interpolatorMetallicity,timeMinimum,timeMaximum,age,metallicity_)
    !![
    <optionalArgument name="parallelize" defaultsTo=".false." />
    !!]
    
    allocate(luminosityMean(self%countLines,size(starFormationHistory%data,dim=1),size(starFormationHistory%data,dim=2)))
    counter       =-1
    counterMaximum=product     ([size(luminosityMean,dim=1              ),size(starFormationHistory%data,dim=1              ),size(starFormationHistory%data,dim=2              )])
    state         =multiCounter([size(luminosityMean,dim=1,kind=c_size_t),size(starFormationHistory%data,dim=1,kind=c_size_t),size(starFormationHistory%data,dim=2,kind=c_size_t)])
    stateLock     =ompLock ()
    !$omp parallel private (iTime,iMetallicity,metallicityMinimum,metallicityMaximum)
    allocate(integratorTime       )
    allocate(integratorMetallicity)
    integratorTime       =integrator(emissionLineLuminosityIntegrandTime       ,toleranceRelative=self%toleranceRelative)
    integratorMetallicity=integrator(emissionLineLuminosityIntegrandMetallicity,toleranceRelative=self%toleranceRelative)
    allocate(interpolatorTime       )
    allocate(interpolatorMetallicity)
    interpolatorTime       =interpolator(self%ages         )
    interpolatorMetallicity=interpolator(self%metallicities)
    !$omp master
    if (parallelize_) then
       !$omp critical(gfortranInternalIO)
       write (label,'(f12.8)') time
       !$omp end critical(gfortranInternalIO)
       call displayIndent("computing template emission line luminosities for time "//trim(adjustl(label))//" Gyr",verbosityLevelWorking)
    end if
    !$omp end master
    ! Iterate over (time,metallicity).
    !$omp do
    do iterator=0,counterMaximum-1
       call stateLock%  set()
       if (state%increment()) then
          iLine       =state%state(1_c_size_t)
          iTime       =state%state(2_c_size_t)
          iMetallicity=state%state(3_c_size_t)
       else
          iLine       =0_c_size_t
          iTime       =0_c_size_t
          iMetallicity=0_c_size_t
          call Error_Report('unable to increment counter'//{introspection:location})
       end if
       call stateLock%unset()
       if (parallelize_) then
          !$omp atomic
          counter=counter+1
          call displayCounter(percentageComplete=int(100.0d0*dble(counter)/dble(counterMaximum)),isNew=counter==0,verbosity=verbosityLevelWorking)
       end if       
       ! Determine times.
       if (iTime == 1) then
          timeMinimum=                         0.0d0
       else
          timeMinimum=    starFormationHistory%time(iTime-1)
       end if
       timeMaximum   =min(starFormationHistory%time(iTime  ),time)
       if (timeMaximum <= timeMinimum) cycle
       ! Determine metallicities.
       if (iMetallicity == 1) then
          metallicityMinimum=                                                   self%metallicityPopulationMinimum
       else
          metallicityMinimum=min(max(self%metallicityBoundaries(iMetallicity-1),self%metallicityPopulationMinimum),self%metallicityPopulationMaximum)
       end if
       metallicityMaximum   =min(max(self%metallicityBoundaries(iMetallicity  ),self%metallicityPopulationMinimum),self%metallicityPopulationMaximum)
       if (metallicityMaximum > metallicityMinimum) then
          luminosityMean(iLine,iTime,iMetallicity)=+integratorMetallicity%integrate(metallicityMinimum,metallicityMaximum) &
               &                                            /                      (timeMaximum       -timeMinimum       ) &
               &                                            /                      (metallicityMaximum-metallicityMinimum)
       else
          metallicity_                            =metallicityMinimum
          luminosityMean(iLine,iTime,iMetallicity)=+integratorTime       %integrate(timeMinimum       ,timeMaximum       ) &
               &                                            /                      (timeMaximum       -timeMinimum       )
       end if
    end do
    !$omp end do
    !$omp master
    if (parallelize_) then
       call displayCounterClear(       verbosityLevelWorking)
       call displayUnindent    ("done",verbosityLevelWorking)
    end if
    !$omp end master
    deallocate(integratorTime         )
    deallocate(integratorMetallicity  )
    deallocate(interpolatorTime       )
    deallocate(interpolatorMetallicity)
    !$omp end parallel
    return

  contains

    double precision function emissionLineLuminosityIntegrandMetallicity(metallicity) result(integrand)
      !!{
      Integrand over metallicity of the stellar population.
      !!}
      implicit none
      double precision, intent(in   ) :: metallicity

      metallicity_=metallicity
      integrand   =integratorTime%integrate(timeMinimum,timeMaximum)
      return
    end function emissionLineLuminosityIntegrandMetallicity

    double precision function emissionLineLuminosityIntegrandTime(timeBirth) result(integrand)
      !!{
      Integrand over birth time of the stellar population.
      !!}
      implicit none
      double precision          , intent(in   )  :: timeBirth
      integer                                    :: iTime                , iMetallicity
      integer         (c_size_t), dimension(0:1) :: interpolateIndexTime , interpolateIndexMetallicity
      double precision          , dimension(0:1) :: interpolateFactorTime, interpolateFactorMetallicity
      
      age    =min(                            &
           &      +     time                  &
           &      -     timeBirth           , &
           &      +self%agePopulationMaximum  &
           &     )
      call interpolatorTime       %linearFactors(age         ,interpolateIndexTime       (0),interpolateFactorTime       )
      call interpolatorMetallicity%linearFactors(metallicity_,interpolateIndexMetallicity(0),interpolateFactorMetallicity)
      interpolateIndexTime       (1)=interpolateIndexTime       (0)+1
      interpolateIndexMetallicity(1)=interpolateIndexMetallicity(0)+1
      integrand                     =0.0d0
      do iTime=0,1
         do iMetallicity=0,1
           integrand=+integrand                                                                                                              &
                 &   +self%luminositiesReduced                 (interpolateIndexTime(iTime),interpolateIndexMetallicity(iMetallicity),iLine) &
                 &   *self%ionizingLuminosityHydrogenNormalized(interpolateIndexTime(iTime),interpolateIndexMetallicity(iMetallicity)      ) &
                 &   /self%ionizingLuminosityHydrogenMean                                                                                    &
                 &   *     interpolateFactorTime                                    (iTime)                                                  &
                 &   *     interpolateFactorMetallicity                                                                (iMetallicity)
         end do
      end do
      return
    end function emissionLineLuminosityIntegrandTime

  end function emissionLineLuminosityMean

  function emissionLineLuminosityHistoryHashedDescriptor(self,starFormationHistory) result(hashedDescriptor)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters    , only : inputParameters
    use :: String_Handling     , only : String_C_To_Fortran
    use :: Hashes_Cryptographic, only : Hash_MD5
    use :: FoX_DOM             , only : setLiveNodeLists
    use :: Histories           , only : history
    implicit none
    type     (varying_string                             )                :: hashedDescriptor
    class    (nodePropertyExtractorLuminosityEmissionLine), intent(in   ) :: self
    type     (history                                    ), intent(in   ) :: starFormationHistory
    character(len=18                                     )                :: parameterLabel
    type     (inputParameters                            )                :: descriptor
    type     (varying_string                             )                :: descriptorString    , values
    integer                                                               :: i
    !![
    <workaround type="gfortran" PR="102845" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=102845">
      <description>
	Memory leak possibly due to OpenMP parallelism, or some failing of gfortran.
      </description>
    </workaround>
    !!]
    ! type     (varying_string          ), save          :: descriptorStringPrevious  , hashedDescriptorPrevious
    ! !$omp threadprivate(descriptorStringPrevious,hashedDescriptorPrevious)
    
    descriptor=inputParameters()
    call setLiveNodeLists(descriptor%document,.false.)
    call self%starFormationHistory_       %descriptor(descriptor)
    call self%outputTimes_                %descriptor(descriptor)
    call self%hiiRegionLuminosityFunction_%descriptor(descriptor)
    values=""
    do i=1,size(self%metallicityBoundaries)
       !$omp critical(gfortranInternalIO)
       write (parameterLabel,'(e17.10)') self                %metallicityBoundaries(i)
       !$omp end critical(gfortranInternalIO)
       values=values//trim(adjustl(parameterLabel))
       if (i < size(self%metallicityBoundaries)) values=values//":"
    end do
    call descriptor%addParameter('metallicity',char(values))
    values=""
    do i=1,size(starFormationHistory%time)
       !$omp critical(gfortranInternalIO)
       write (parameterLabel,'(e17.10)') starFormationHistory%time                 (i)
       !$omp end critical(gfortranInternalIO)
       values=values//trim(adjustl(parameterLabel))
       if (i < size(starFormationHistory%time)) values=values//":"
    end do
    call descriptor%addParameter('time'       ,char(values))
    descriptorString=descriptor%serializeToString()
    call descriptor%destroy()
    descriptorString=descriptorString//":sourceDigest{"//String_C_To_Fortran(nodePropertyExtractorLuminosityEmissionLine5)//"}"
    !![
    <workaround type="gfortran" PR="102845" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=102845">
     <description>
      Memory leak possibly due to OpenMP parallelism, or some failing of gfortran.
     </description>
    </workaround>
    !!]
    hashedDescriptor=Hash_MD5(descriptorString)
    return
  end function emissionLineLuminosityHistoryHashedDescriptor



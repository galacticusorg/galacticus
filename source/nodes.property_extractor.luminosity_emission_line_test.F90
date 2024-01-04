!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Contains a module which implements a property extractor class for the emission line luminosity of a component.
  !!}
  use :: Cosmology_Functions                   , only : cosmologyFunctionsClass
  use :: Galactic_Structure_Options            , only : enumerationComponentTypeType
  use :: Output_Times                          , only : outputTimesClass
  use :: Star_Formation_Histories              , only : starFormationHistoryClass
  use :: hii_Region_Luminosity_Functions       , only : hiiRegionLuminosityFunctionClass
  use :: Stellar_Luminosities_Structure        , only : enumerationFrameType, enumerationFrameDecode
  type:: emissionLineLuminosityTemplate
     !!{
     Type used to store luminosity templates for emission lines.
     !!}
     private
     integer         (c_size_t)                                :: countlines=-1_c_size_t
     double precision          , allocatable, dimension(:    ) :: lines
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
     class           (starFormationHistoryClass                 ), pointer                   :: starFormationHistory_                  => null()
     class           (outputTimesClass                          ), pointer                   :: outputTimes_                           => null()
     class           (cosmologyFunctionsClass                   ), pointer                   :: cosmologyFunctions_                    => null()
     class           (hiiRegionLuminosityFunctionClass          ), pointer                   :: hiiRegionLuminosityFunction_           => null()
     type            (varying_string                            )                            :: name_                                   , description_
     type            (enumerationComponentTypeType              )                            :: component
     type            (enumerationFrameType                      )                            :: frame
     integer                                                                                 :: countLines
     double precision                                            , allocatable, dimension(:) :: lines_                                    , metallicityBoundaries
     type            (emissionLineLuminosityTemplate            ), allocatable, dimension(:) :: templates
     double precision                                                                        :: metallicityPopulationMinimum, metallicityPopulationMaximum, &
                                                                          &                     agePopulationMaximum, resolution, factorWavelength
     integer                                                                                 :: abundanceIndex
     logical                                                                                 :: useluminosityTemplates
   contains
     !![
     <methods>
       <method description="Return a hashed descriptor of the object which incorporates the time and metallicity binning of the star formation history." method="luminosityHistoryHashedDescriptor"/>
       <method description="Compute the mean luminosity of the stellar population in the given bin of the star formation history."                       method="luminosityMean"         />
       <method description="Return the index of the template time to use."                                                                               method="indexTemplateTime"      />
       <method description="Return the index of the template luminosities to use."                                                                       method="indexTemplateNode"      />
     </methods>
     !!]
     final     ::                            emissionLineLuminosityDestructor
     procedure :: luminosityHistoryHashedDescriptor => emissionLineLuminosityHistoryHashedDescriptor
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
    Constructor for the {\normalfont \ttfamily sed} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type            (nodePropertyExtractorLuminosityEmissionLine                  )                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (stellarPopulationSpectraClass             ), pointer       :: stellarPopulationSpectra_
    class           (stellarPopulationSpectraPostprocessorClass), pointer       :: stellarPopulationSpectraPostprocessor_
    class           (starFormationHistoryClass                 ), pointer       :: starFormationHistory_
    class           (outputTimesClass                          ), pointer       :: outputTimes_
    class           (cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_
    class           (hiiRegionLuminosityFunctionClass          ), pointer       :: hiiRegionLuminosityFunction_
    type            (varying_string                            )                :: component                             
    ! Read the table of emission line luminosities.
    !$ call hdf5Access%set()
    call emissionLinesFile%openFile(char(inputPath(pathTypeDataStatic))//"hiiRegions/emissionLines.hdf5",readOnly=.true.)
    lines=emissionLinesFile%openGroup('lines')
    do i=1,size(lineNames)
       if (.not.lines%hasDataset(char(self%lineNames(i)))) call Error_Report('line "'//char(self%lineNames(i))//'" not found'//{introspection:location})
    end do
    call emissionLinesFile%readDataset('metallicity'                  ,self%metallicity                 )
    call emissionLinesFile%readDataset('ionizingFluxHydrogen'         ,self%ionizingFluxHydrogen        )
    allocate(                                          &
         &   self%luminosity                           &
         &   (                                         &
         &    size(self%ionizingFluxHydrogen        ), &
         &    size(self%metallicity                 ), &
         &    size(self%lineNames                   )  &
         &   )                                         &
         &  )
    allocate(                                          &
         &   self%wavelength                           &
         &   (                                         &
         &    size(self%lineNames                   )  &
         &   )                                         &
         &  )
    do i=1,size(lineNames)
       call lines      %readDatasetStatic(char(self%lineNames(i)),self%luminosity(:,:,:,:,:,i))
       lineDataset=lines%openDataset(char(self%lineNames(i)))
       call lineDataset%readAttribute('wavelength',self%wavelength(i))
       call lineDataset%close        (                               )
    end do
    call lines            %close      (                                                                 )
    call emissionLinesFile%close      (                                                                 )
    !$ call hdf5Access%unset()
    ! Convert parameters and luminosities to log form.
    self%metallicity                 =log10(self%metallicity                 )
    self%densityHydrogen             =log10(self%densityHydrogen             )
    self%ionizingFluxHydrogen        =log10(self%ionizingFluxHydrogen        )
    self%luminosity                  =log10(self%luminosity                  )

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation rate.</description>
    </inputParameter>
    <objectBuilder class="starFormationHistory"                  name="starFormationHistory_"                  source="parameters"/>
    <objectBuilder class="outputTimes"                           name="outputTimes_"                           source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                    name="cosmologyFunctions_"                    source="parameters"/>
    <objectBuilder class="hiiRegionLuminosityFunction"           name="hiiRegionLuminosityFunction_"           source="parameters"/>
    !!]
    self=nodePropertyExtractorLuminosityEmissionLine(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),enumerationFrameEncode(char(frame),includesPrefix=.false.),toleranceRelative,starFormationHistory_,outputTimes_,cosmologyFunctions_,hiiRegionLuminosityFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"                 />
    <objectDestructor name="outputTimes_"                          />
    <objectDestructor name="cosmologyFunctions_"                   />
    <objectDestructor name="hiiRegionLuminosityFunction_"   />
    !!]
    return
  end function emissionLineLuminosityConstructorParameters

  function emissionLineLuminosityConstructorInternal(component,frame,starFormationHistory_,outputTimes_,cosmologyFunctions_,hiiRegionLuminosityFunction_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sed} property extractor class.
    !!}
    use :: Atomic_Data                     , only : Abundance_Pattern_Lookup
    use :: Galactic_Structure_Options      , only : componentTypeDisk       , componentTypeSpheroid
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    type            (nodePropertyExtractorLuminosityEmissionLine              )               :: self
    type            (enumerationComponentTypeType              ), intent(in   )               :: component
    type            (enumerationFrameType                      ), intent(in   )               :: frame
    class           (starFormationHistoryClass                 ), intent(in   ), target       :: starFormationHistory_
    class           (outputTimesClass                          ), intent(in   ), target       :: outputTimes_
    class           (cosmologyFunctionsClass                   ), intent(in   ), target       :: cosmologyFunctions_
    class           (hiiRegionLuminosityFunctionClass          ), intent(in   ), target       :: hiiRegionLuminosityFunction_
    double precision                                            ,                             :: exponent
    double precision                                            , allocatable  , dimension(:) :: ages                                  , metallicities
    integer                                                                                   :: agesCount                             , metallicitiesCount
    
    !![
    <constructorAssign variables="component, frame, *starFormationHistory_, *outputTimes_, *cosmologyFunctions_,*hiiRegionLuminosityFunction_"/>
    !!]
    
    if     (                                                                                                               &
         &   component /= componentTypeDisk                                                                                &
         &  .and.                                                                                                          &
         &   component /= componentTypeSpheroid                                                                            &
         & ) call Error_Report("only 'disk' and 'spheroid' components are supported"//{introspection:location})
    call self%hiiRegionLuminosityFunction_%lines(self%countLines                   ,self%lines_              )
    call self%hiiRegionLuminosityFunction_%tabulation (     agesCount       ,metallicitiesCount,     ages        ,metallicities)    
    self%metallicityBoundaries       =self%starFormationHistory_%metallicityBoundaries()
    self%agePopulationMaximum        =ages         (agesCount         )
    self%metallicityPopulationMaximum=metallicities(metallicitiesCount)/metallicitySolar
    self%metallicityPopulationMinimum=metallicities(                 1)/metallicitySolar
    self%abundanceIndex              =Abundance_Pattern_Lookup(abundanceName="solar")
    self%useluminosityTemplates      =self%starFormationHistory_%perOutputTabulationIsStatic()
    
    ! Read the table of emission line luminosities.
    !$ call hdf5Access%set()
    call emissionLinesFile%openFile(char(inputPath(pathTypeDataStatic))//"hiiRegions/emissionLines.hdf5",readOnly=.true.)
    lines=emissionLinesFile%openGroup('lines')
    do i=1,size(lineNames)
       if (.not.lines%hasDataset(char(self%lineNames(i)))) call Error_Report('line "'//char(self%lineNames(i))//'" not found'//{introspection:location})
    end do
    call emissionLinesFile%readDataset('metallicity'                  ,self%metallicity                 )
    call emissionLinesFile%readDataset('age'                          ,self%age                         )
    call emissionLinesFile%readDataset('ionizingFluxHydrogen'         ,self%ionizingFluxHydrogen        )


    allocate(                                          &
         &   self%luminosity                           &
         &   (                                         &
         &    size(self%lineNames                   )  &
         &    size(self%age                         )  &
         &    size(self%metallicity                 ), &
         &    size(self%ionizingFluxHydrogen        ), &
         &   )     

    do i=1,size(lineNames)
       call lines      %readDatasetStatic(char(self%lineNames(i)),self%luminosity(:,:,:,:,:,i))
       lineDataset=lines%openDataset(char(self%lineNames(i)))
       call lineDataset%close        (                               )
    end do
    call lines            %close      (                                                                 )
    call emissionLinesFile%close      (                                                                 )
    !$ call hdf5Access%unset()
    ! Convert parameters and luminosities to log form.
    self%metallicity                 =log10(self%metallicity                 )    
    self%ionizingFluxHydrogen        =log10(self%ionizingFluxHydrogen        )
    self%luminosity                  =log10(self%luminosity                  )

    ! Find indices of ionizing continuua filters.
    allocate(self%ionizingContinuumIndex(self%outputTimes_%count(),3_c_size_t))
    do i=1,self%outputTimes_%count()
       if (present(outputMask).and..not.outputMask(i)) then
          self%ionizingContinuumIndex(i,:                        )=-1
       else
          self%ionizingContinuumIndex(i,ionizingContinuumHydrogen%ID)=unitStellarLuminosities%index('Lyc'            ,'rest',self%outputTimes_%redshift(i))
          self%ionizingContinuumIndex(i,ionizingContinuumHelium  %ID)=unitStellarLuminosities%index('HeliumContinuum','rest',self%outputTimes_%redshift(i))
          self%ionizingContinuumIndex(i,ionizingContinuumOxygen  %ID)=unitStellarLuminosities%index('OxygenContinuum','rest',self%outputTimes_%redshift(i))
       end if
    end do
    ! Read wavelength intervals of ionizing continuum filters.
    self%filterExtent(:,ionizingContinuumHydrogen%ID)=Filter_Extent(Filter_Get_Index(var_str('Lyc'            )))
    self%filterExtent(:,ionizingContinuumHelium  %ID)=Filter_Extent(Filter_Get_Index(var_str('HeliumContinuum')))
    self%filterExtent(:,ionizingContinuumOxygen  %ID)=Filter_Extent(Filter_Get_Index(var_str('OxygenContinuum')))

    !Calculate emission line luminosity for some age and metallicity
    do line=1,size(self%luminosity,dim=6)
        do k=1,size(self%ionizingFluxHydrogen)
           rateHydrogenIonizingPhotonsMinimum=10.0d0**(self%ionizingFluxHydrogen(k)-0.5d0) 
           rateHydrogenIonizingPhotonsMaximum=10.0d0**(self%ionizingFluxHydrogen(k)+0.5d0)
           ! luminosity at metallicity and age 
           self%emissionLineLuminosity_age_z=self%emissionLineLuminosity_age_z                             &
                     &             *self%hiiRegionLuminosityFunctionis_%powerLawCumulativeLuminosity(rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum) &                      &                *10.0d0**self%luminosity(          &
                     &             *10.0d0**self%luminosity(             &
                     &                                      line,        &
                     &                                      :,           &
                     &                                      :,           &
                     &                                      k,           &
                     &                                     )
	end do
    end do

     ! Construct name and description.
    self%name_       ="luminosityEmissionLine:"//String_Join(lineNames,"+")
    self%description_="Luminosity of the "     //String_Join(lineNames,"+")//" emission line"
    if (size(lineNames) > 1) self%description_=self%description_//"s"
    self%description_=self%description_//" [ergs/s]"	
    
  end function emissionLineLuminosityConstructorInternal

  subroutine emissionLineLuminosityDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily emission line luminosity} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self
    !![
    <objectDestructor name="self%starFormationHistory_"                 />
    <objectDestructor name="self%outputTimes_"                          />
    <objectDestructor name="self%cosmologyFunctions_"                   />
    <objectDestructor name="self%hiiRegionLuminosityFunction_"          />
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

    emissionLineLuminosityElementCount=1
    return
  end function emissionLineLuminosityElementCount


  function emissionLineLuminosityExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily luminosityEmissionLine} property extractor.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentDisk, nodeComponentSpheroid
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid
    use :: Histories                 , only : history
    implicit none
    double precision                          , dimension(:,:  )          , allocatable :: emissionLineLuminosityExtract
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)   , target                :: self
    type            (treeNode                ), intent(inout)   , target                :: node
    double precision                          , intent(in   )                           :: time
    type            (multiCounter            ), intent(inout)   , optional              :: instance
    class           (nodeComponentDisk       )                  , pointer               :: disk
    class           (nodeComponentSpheroid   )                  , pointer               :: spheroid
    double precision                          , dimension(:,:,:), pointer               :: luminosityTemplate_
    double precision                          , dimension(:,:,:), target  , allocatable :: luminosityTemplate
    type            (history                 )                                          :: starFormationHistory
    integer                                                                             :: indexTemplate       , iLines
    !$GLC attributes unused :: instance

    allocate(emissionLineLuminosityExtract(self%size(time),1))
    emissionLineLuminosityExtract=0.0d0
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
       luminosityTemplate_ => self%templates(indexTemplate)%luminosityEmissionLine
    else
       ! Stored templates can not be used, get the templates for this specific case, and point to them.
       luminosityTemplate  =  self%luminosityMean(time,starFormationHistory)
       luminosityTemplate_ => luminosityTemplate
    end if    
    do iLines=1,size(emissionLineExtract,dim=1)
       emissionLineExtract(iLines,1)=sum(luminosityTemplate_(iLines,:,:)*starFormationHistory%data(:,:))
    end do
    return
  end function emissionLineLuminosityExtract

  subroutine emissionLineLuminosityNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily emissionLines}.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self
    double precision                          , intent(in   )     :: time
    type            (varying_string            ), intent(inout), dimension(:) , allocatable :: names
allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(1))
    names(1)=enumerationComponentTypeDecode(self%component,includePrefix=.false.)//"StellarEmissionLineLuminosity:"//self%lines_%objectType(short=.true.)
    return
  end subroutine emissionLineLuminosityNames
  
  subroutine emissionLineLuminosityDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily emission line luminosity} property.
    !!}
    implicit none
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)            :: self
    double precision                                       , intent(in   )                  :: time
    type            (varying_string            ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(1))
     descriptions(1)=var_str('Emission line luminosity of the galaxy [ergs].'   )
    return
  end subroutine emissionLineLuminosityDescriptions
  
  function emissionLineLuminosityUnitsInSI(self,time)
  !!{
    Return the units of the {\normalfont \ttfamily emissionLineLuminosity} properties in the SI
system.
    !!}
    use :: Numerical_Constants_Units, only : ergs
    implicit none
    double precision                          , allocatable  , dimension(:)      ::
emissionLineLuminosityUnitsInSI
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout) ::
self
    double precision                          , intent(in   )                    :: time
    type            (varying_string          ), intent(inout), dimension(:) ,
allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(emissionLineLuminosityUnitsInSI(1))
    emissionLineLuminosityUnitsInSI(1)=ergs
    return
  end function emissionLineLuminosityUnitsInSI

  integer function emissionLineLuminosityIndexTemplateTime(self,time)
    !!{
    Find the index of the template emission lines to use.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self
    double precision                          , intent(in   ) :: time
    integer         (c_size_t                )                :: indexOutput

    emsiionLineLuminosityIndexTemplateTime=-1
    if (.not.self%useluminosityTemplates) return
    ! Check that the time is an output time.
    indexOutput=self%outputTimes_%index(time,findClosest=.true.)
    if (.not.Values_Agree(time,self%outputTimes_%time(indexOutput),relTol=1.0d-6)) return
    emissionLineLuminosityIndexTemplateTime=int(indexOutput)
    return
  end function emissionLineLuminosityIndexTemplateTime

  integer function emissionLineLuminosityIndexTemplateNode(self,node,starFormationHistory)
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
    type     (treeNode                ), intent(inout) :: node
    type     (history                 ), intent(in   ) :: starFormationHistory
    class    (nodeComponentBasic      ), pointer       :: basic
    integer  (c_size_t                )                :: indexOutput
    type     (lockDescriptor          )                :: fileLock
    type     (hdf5Object              )                :: file
    type     (varying_string          )                :: fileName
    character(len=16                  )                :: label

    ! Return a negative index if templates are not being used.
    emissionLineIndexTemplateNode=-1
    if (.not.self%useluminosityTemplates) return
    ! Check that the node exists at at output time.
    basic       => node             %basic(                               )
    indexOutput =  self%outputTimes_%index(basic%time(),findClosest=.true.)
    if (.not.Values_Agree(basic%time(),self%outputTimes_%time(indexOutput),relTol=1.0d-6)) return
    emissionLineIndexTemplateNode=int(indexOutput)
    ! Ensure that the templates have been built for this index.
    if (.not.allocated(self%templates)) allocate(self%templates(self%outputTimes_%count()))
    if (.not.allocated(self%templates(emissionLineIndexTemplateNode)%emissionLineLuminosity)) then
       ! Construct the file name.
       fileName=inputPath(pathTypeDataDynamic)                         // &
            &        'stellarPopulations/'                             // &
            &        self%objectType             (                    )// &
            &        '_'                                               // &
            &        self%elHistoryHashedDescriptor(starFormationHistory)// &
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
             call file%readDataset('luminosityTemplate',self%templates(emissionLineIndexTemplateNode)%sed)
          end if
          call file%close()
          !$ call hdf5Access%unset()
       end if
       if (.not.allocated(self%templates(emissionLineIndexTemplateNode)%emissionLineLuminosity)) then
          self%templates(emissionLineIndexTemplateNode)%emissionLineLuminosity=self%luminosityMean(self%outputTimes_%time(indexOutput),starFormationHistory,parallelize=.true.)
          !$omp critical(gfortranInternalIO)
          write (label,'(f12.8)') self%outputTimes_%time(indexOutput)
          !$omp end critical(gfortranInternalIO)
          call displayMessage("storing Emission Line Luminosity tabulation for time "//trim(adjustl(label))//" Gyr to file '"//fileName//"'",verbosityLevelWorking)
          !$ call hdf5Access%set()
          call file%openFile(char(fileName),overWrite=.false.,readOnly=.false.)
          call file%writeDataset(self%templates(emissionLineIndexTemplateNode)%emissionLineLuminosity,'luminosityTemplate')
          call file%close()
          !$ call hdf5Access%unset()
       end if
       call File_Unlock(fileLock)
    end if
    return
  end function emissionLineLuminosityIndexTemplateNode

  double precision function emissionLineLuminosityMean(self,time,starFormationHistory,parallelize)
    !!{
    Compute the mean luminosity of the stellar population in each bin of the star formation history.
    !!}
    use    :: Abundances_Structure , only : abundances             , metallicityTypeLinearByMassSolar, adjustElementsReset
    use    :: Display              , only : displayIndent          , displayUnindent                 , displayCounter     , displayCounterClear, &
         &                                  verbosityLevelWorking
    use    :: Error                , only : Error_Report
    use    :: Histories            , only : history
    use    :: Numerical_Integration, only : integrator
    use    :: Multi_Counters       , only : multiCounter
    use    :: Locks                , only : ompLock
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num
    implicit none
    double precision                                            , dimension(:,:,:), allocatable :: emissionLineLuminosityMean
    class           (nodePropertyExtractorLuminosityEmissionLine                  ), intent(inout)                 :: self
    double precision                                            , intent(in   )                 :: time
    type            (history                                   ), intent(in   )                 :: starFormationHistory
    logical                                                     , intent(in   )   , optional    :: parallelize
    class           (hiiRegionLuminosityFunctionClass         ), pointer         , save        :: hiiRegionLuminosityFunction_
    class           (cosmologyFunctionsClass                   ), pointer         , save        :: cosmologyFunctions_
    type            (integrator                                ), allocatable     , save        :: integratorTime                        , integratorMetallicity
    integer         (c_size_t                                  ), dimension(:    ), allocatable :: jLines
    double precision                                            , dimension(:    ), allocatable :: wavelengthMinima                      , wavelengthMaxima
    integer         (c_size_t                                  )                                :: iTime                                 , iMetallicity,          &                                                                                                                                                         , kWavelength          , &
         &                                                                                         counter                               , counterMaximum       , &
         &                                                                                         iterator
    double precision                                                                            :: metallicityMinimum                    , metallicityMaximum   , &
         &                                                                                         expansionFactor
    double precision                                                              , save        :: timeMinimum                           , timeMaximum          , &
         &                                                                                         lines                                 , age                  , &
         &                                                                                         redshift
    type            (abundances                                )                  , save        :: abundancesStellar
    character       (len=12                                    )                                :: label
    type            (multiCounter                              )                                :: state
    type            (ompLock                                   )                                :: stateLock
    !$omp threadprivate(stellarPopulationSpectra_,stellarPopulationSpectraPostprocessor_,cosmologyFunctions_,integratorTime,integratorWavelength,integratorMetallicity,abundancesStellar,lines,timeMinimum,timeMaximum,age,redshift)
    !![
    <optionalArgument name="parallelize" defaultsTo=".false." />
    !!]
    
    allocate(emissionLineLuminosityMean(self%size(time),size(starFormationHistory%data,dim=1),size(starFormationHistory%data,dim=2)))
    select case (self%frame%ID)
    case (frameRest    %ID)
       expansionFactor=1.0d0
    case (frameObserved%ID)
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
    case default
       expansionFactor=0.0d0
       call Error_Report('unknown frame'//{introspection:location})
    end select

    counter       =-1
    counterMaximum=product     ([size(emissionLineLuminosityMean,dim=1              ),size(starFormationHistory%data,dim=1              ),size(starFormationHistory%data,dim=2              )])
    state         =multiCounter([size(emissionLineLuminosityMean,dim=1,kind=c_size_t),size(starFormationHistory%data,dim=1,kind=c_size_t),size(starFormationHistory%data,dim=2,kind=c_size_t)])
    stateLock     =ompLock     (                                                                                                                                                     )
    !$omp parallel private (iTime,iMetallicity,metallicityMinimum,metallicityMaximum)
    allocate(integratorTime       )
    allocate(integratorMetallicity)
    integratorTime       =integrator(emissionLineLuminosityIntegrandTime       ,toleranceRelative=self%toleranceRelative)
    integratorMetallicity=integrator(emissionLineLuminosityIntegrandMetallicity,toleranceRelative=self%toleranceRelative)
    if (parallelize_) then
       allocate(starFormationHistory_,mold=self%starFormationHistory_)
       allocate(cosmologyFunctions_                   ,mold=self%cosmologyFunctions_                   )
       !$omp critical(nodePropertyExtractLuminosityEmissionLineDeepCopy)
       !![
       <deepCopyReset variables="hiiRegionLuminosityFunction_ self%cosmologyFunctions_ starFormationHistory_"/>
       !!]
       !![
       <deepCopy source="self%starFormationHistory_"                  destination="starFormationHistory_"             />
       <deepCopy source="self%cosmologyFunctions_"                    destination="cosmologyFunctions_"               />
       <deepCopy source="hiiRegionLuminosityFunction_"                destination="hiiRegionLuminosityFunction_"      />
       !!]
       !![
       <deepCopyFinalize variables="hiiRegionLuminosityFunction_ cosmologyFunctions_ starFormationHistory_"/>
       !!]
       !$omp end critical(nodePropertyExtractLuminosityEmissionLineDeepCopy)
    end if
    !$omp master
    if (parallelize_) then
       !$omp critical(gfortranInternalIO)
       write (label,'(f12.8)') time
       !$omp end critical(gfortranInternalIO)
       call displayIndent("computing template SEDs for time "//trim(adjustl(label))//" Gyr",verbosityLevelWorking)
    end if
    !$omp end master
    ! Iterate over (lines,time,metallicity).
    !$omp do
    do iterator=0,counterMaximum-1
       call stateLock%  set()
       if (state%increment()) then
          iLines      =state%state(1_c_size_t)
          iTime       =state%state(2_c_size_t)
          iMetallicity=state%state(3_c_size_t)
       else
          iLines.     =0_c_size_t
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
       ! Determine the line.
       
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
          emissionLineLuminosityMean(iLines,iTime,iMetallicity)=+integratorMetallicity%integrate(metallicityMinimum,metallicityMaximum) &
               &                                            /                               (timeMaximum       -timeMinimum       ) &
               &                                            /                               (metallicityMaximum-metallicityMinimum)
       else
          call abundancesStellar%metallicitySet(                                                       &
               &                                metallicity    =     metallicityMinimum              , &
               &                                metallicityType=     metallicityTypeLinearByMassSolar, &
               &                                adjustElements =     adjustElementsReset             , &
               &                                abundanceIndex =self%abundanceIndex                    &
               &                               )
          emissionLineLuminosityMean(iLines,iTime,iMetallicity)=+integratorTime       %integrate(timeMinimum       ,timeMaximum       ) &
               &                                            /                               (timeMaximum       -timeMinimum       )
       end if
    end do
    !$omp end do
    !$omp master
    if (parallelize_) then
       call displayCounterClear(       verbosityLevelWorking)
       call displayUnindent    ("done",verbosityLevelWorking)
    end if
    !$omp end master
    !![
    <objectDestructor name="cosmologyFunctions_"                   />
    !!]
    deallocate(integratorTime       )
    deallocate(integratorMetallicity)
    !$omp end parallel
    return

  contains

    double precision function emissionLineLuminosityIntegrandMetallicity(metallicity)
      !!{
      Integrand over metallicity of the stellar population.
      !!}
      implicit none
      double precision, intent(in   ) :: metallicity

      call abundancesStellar%metallicitySet(                                                       &
           &                                metallicity    =     metallicity                     , &
           &                                metallicityType=     metallicityTypeLinearByMassSolar, &
           &                                adjustElements =     adjustElementsReset             , &
           &                                abundanceIndex =self%abundanceIndex                    &
           &                               )
      emissionLineLuminosityIntegrandMetallicity=integratorTime%integrate(timeMinimum,timeMaximum)
      return
    end function emissionLineLuminosityIntegrandMetallicity

    double precision function emissionLineLuminosityIntegrandTime(timeBirth)
      !!{
      Integrand over birth time of the stellar population.
      !!}
      implicit none
      double precision, intent(in   ) :: timeBirth

      age             =min(                            &
           &               +     time                  &
           &               -     timeBirth           , &
           &               +self%agePopulationMaximum  &
           &              )
      if (parallelize_) then
         redshift=     cosmologyFunctions_%redshiftFromExpansionFactor(     cosmologyFunctions_%expansionFactor(time))
      else
         redshift=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(time))
      end if
      if (parallelize_) then
      	emissionLineLuminosityIntegrandTime=+self%starFormationHistory_             *self%emissionLineLuminosity_age_z
      end if         
      return
    end function emissionLineLuminosityIntegrandTime

    
  end function emissionLineLuminosityMean

  function emissionLineLuminosityHistoryHashedDescriptor(self,starFormationHistory)  
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters    , only : inputParameters
    use :: String_Handling     , only : String_C_To_Fortran
    use :: Hashes_Cryptographic, only : Hash_MD5
    use :: FoX_DOM             , only : setLiveNodeLists
    use :: Histories           , only : history
    implicit none
    type     (varying_string          )                :: emissionLineLuminosityHistoryHashedDescriptor
    class    (nodePropertyExtractorLuminosityEmissionLine), intent(in   ) :: self
    type     (history                 ), intent(in   ) :: starFormationHistory
    character(len=18                  )                :: parameterLabel
    type     (inputParameters         )                :: descriptor
    type     (varying_string          )                :: descriptorString          , values
    integer                                            :: i
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
    call descriptor%addParameter('frame'            ,char(enumerationFrameDecode(self%frame,includePrefix=.false.)))
    call self%starFormationHistory_                 %descriptor(descriptor)
    call self%outputTimes_                          %descriptor(descriptor)
    call self%cosmologyFunctions_                   %descriptor(descriptor)
    call self%hiiRegionLuminosityFunction_          %descriptor(descriptor)
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
    emissionLineLuminosityHistoryHashedDescriptor=Hash_MD5(descriptorString)
    return
  end function emissionLineLuminosityHistoryHashedDescriptor

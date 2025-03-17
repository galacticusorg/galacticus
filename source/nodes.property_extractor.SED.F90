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
  Implements a property extractor class for the SED of a component.
  !!}
  use :: Cosmology_Functions                   , only : cosmologyFunctionsClass
  use :: Galactic_Structure_Options            , only : enumerationComponentTypeType
  use :: Output_Times                          , only : outputTimesClass
  use :: Stellar_Population_Spectra            , only : stellarPopulationSpectraClass
  use :: Stellar_Population_Spectra_Postprocess, only : stellarPopulationSpectraPostprocessorClass
  use :: Star_Formation_Histories              , only : starFormationHistoryClass
  use :: Stellar_Luminosities_Structure        , only : enumerationFrameType

  type :: sedTemplate
     !!{
     Type used to store SED templates.
     !!}
     private
     integer         (c_size_t)                                :: countWavelengths=-1_c_size_t
     double precision          , allocatable, dimension(:    ) :: wavelength
     double precision          , allocatable, dimension(:,:,:) :: sed
  end type sedTemplate
     
  !![
  <nodePropertyExtractor name="nodePropertyExtractorSED">
    <description>A property extractor class for the SED of a component.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorSED
     !!{
     A property extractor class for the SED of a component.
     !!}
     private
     class           (stellarPopulationSpectraClass             ), pointer                   :: stellarPopulationSpectra_              => null()
     class           (stellarPopulationSpectraPostprocessorClass), pointer                   :: stellarPopulationSpectraPostprocessor_ => null()
     class           (starFormationHistoryClass                 ), pointer                   :: starFormationHistory_                  => null()
     class           (outputTimesClass                          ), pointer                   :: outputTimes_                           => null()
     class           (cosmologyFunctionsClass                   ), pointer                   :: cosmologyFunctions_                    => null()
     type            (enumerationComponentTypeType              )                            :: component
     integer                                                                                 :: countWavelengths
     double precision                                            , allocatable, dimension(:) :: wavelengths_                                    , metallicityBoundaries
     type            (sedTemplate                               ), allocatable, dimension(:) :: templates
     double precision                                                                        :: metallicityPopulationMinimum                    , metallicityPopulationMaximum, &
          &                                                                                     wavelengthMinimum                               , wavelengthMaximum           , &
          &                                                                                     agePopulationMaximum                            , resolution                  , &
          &                                                                                     factorWavelength                                , toleranceRelative
     integer                                                                                 :: abundanceIndex
     type            (enumerationFrameType                      )                            :: frame
   contains
     !![
     <methods>
       <method description="Return the index of the template SEDs to use."                                                                               method="indexTemplateTime"      />
       <method description="Return the index of the template SEDs to use."                                                                               method="indexTemplateNode"      />
       <method description="Compute the mean luminosity of the stellar population in the given bin of the star formation history."                       method="luminosityMean"         />
       <method description="Return a hashed descriptor of the object which incorporates the time and metallicity binning of the star formation history." method="historyHashedDescriptor"/>
       <method description="Return an array of the wavelengths at which the SED is computed."                                                            method="wavelengths"            />
     </methods>
     !!]
     final     ::                            sedDestructor
     procedure :: historyHashedDescriptor => sedHistoryHashedDescriptor
     procedure :: wavelengths             => sedWavelengths
     procedure :: columnDescriptions      => sedColumnDescriptions
     procedure :: size                    => sedSize
     procedure :: elementCount            => sedElementCount
     procedure :: extract                 => sedExtract
     procedure :: names                   => sedNames
     procedure :: descriptions            => sedDescriptions
     procedure :: unitsInSI               => sedUnitsInSI
     procedure :: luminosityMean          => sedLuminosityMean
     procedure :: indexTemplateTime       => sedIndexTemplateTime
     procedure :: indexTemplateNode       => sedIndexTemplateNode
  end type nodePropertyExtractorSED
  
  interface nodePropertyExtractorSED
     !!{
     Constructors for the ``sed'' output analysis class.
     !!}
     module procedure sedConstructorParameters
     module procedure sedConstructorInternal
  end interface nodePropertyExtractorSED
      
contains

  function sedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sed} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters              , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options    , only : enumerationComponentTypeEncode
    use :: Stellar_Luminosities_Structure, only : enumerationFrameEncode
    implicit none
    type            (nodePropertyExtractorSED                  )                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (stellarPopulationSpectraClass             ), pointer       :: stellarPopulationSpectra_
    class           (stellarPopulationSpectraPostprocessorClass), pointer       :: stellarPopulationSpectraPostprocessor_
    class           (starFormationHistoryClass                 ), pointer       :: starFormationHistory_
    class           (outputTimesClass                          ), pointer       :: outputTimes_
    class           (cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_
    type            (varying_string                            )                :: component                             , frame
    double precision                                                            :: wavelengthMinimum                     , wavelengthMaximum, &
         &                                                                         resolution                            , toleranceRelative
    
    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation rate.</description>
    </inputParameter>
    <inputParameter>
      <name>frame</name>
      <source>parameters</source>
      <defaultValue>var_str('rest')</defaultValue>
      <description>The frame ({\normalfont \ttfamily rest} or {\normalfont \ttfamily observed}) for which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum wavelength at which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMaximum</name>
      <source>parameters</source>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>The maximum wavelength at which to compute the SED.</description>
    </inputParameter>
    <inputParameter>
      <name>resolution</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>The resolution, $\lambda/\Delta\lambda$, at which to compute the SED. If a negative value is given the SED will be computed at the full resolution provided by the stellar population spectra class.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelative</name>
      <source>parameters</source>
      <defaultValue>1.0d-3</defaultValue>
      <description>The relative tolerance used in integration over stellar population spectra.</description>
    </inputParameter>
    <objectBuilder class="stellarPopulationSpectraPostprocessor" name="stellarPopulationSpectraPostprocessor_" source="parameters"/>
    <objectBuilder class="stellarPopulationSpectra"              name="stellarPopulationSpectra_"              source="parameters"/>
    <objectBuilder class="starFormationHistory"                  name="starFormationHistory_"                  source="parameters"/>
    <objectBuilder class="outputTimes"                           name="outputTimes_"                           source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                    name="cosmologyFunctions_"                    source="parameters"/>
    !!]
    self=nodePropertyExtractorSED(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),enumerationFrameEncode(char(frame),includesPrefix=.false.),wavelengthMinimum,wavelengthMaximum,resolution,toleranceRelative,stellarPopulationSpectra_,stellarPopulationSpectraPostprocessor_,starFormationHistory_,outputTimes_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarPopulationSpectra_"             />
    <objectDestructor name="stellarPopulationSpectraPostprocessor_"/>
    <objectDestructor name="starFormationHistory_"                 />
    <objectDestructor name="outputTimes_"                          />
    <objectDestructor name="cosmologyFunctions_"                   />
    !!]
    return
  end function sedConstructorParameters

  function sedConstructorInternal(component,frame,wavelengthMinimum,wavelengthMaximum,resolution,toleranceRelative,stellarPopulationSpectra_,stellarPopulationSpectraPostprocessor_,starFormationHistory_,outputTimes_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sed} property extractor class.
    !!}
    use :: Atomic_Data                     , only : Abundance_Pattern_Lookup
    use :: Galactic_Structure_Options      , only : componentTypeDisk       , componentTypeSpheroid, componentTypeNuclearStarCluster
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    type            (nodePropertyExtractorSED                  )                              :: self
    type            (enumerationComponentTypeType              ), intent(in   )               :: component
    type            (enumerationFrameType                      ), intent(in   )               :: frame
    class           (stellarPopulationSpectraClass             ), intent(in   ), target       :: stellarPopulationSpectra_
    class           (stellarPopulationSpectraPostprocessorClass), intent(in   ), target       :: stellarPopulationSpectraPostprocessor_
    class           (starFormationHistoryClass                 ), intent(in   ), target       :: starFormationHistory_
    class           (outputTimesClass                          ), intent(in   ), target       :: outputTimes_
    class           (cosmologyFunctionsClass                   ), intent(in   ), target       :: cosmologyFunctions_
    double precision                                            , intent(in   )               :: wavelengthMinimum                     , wavelengthMaximum, &
         &                                                                                       resolution                            , toleranceRelative
    double precision                                            , allocatable  , dimension(:) :: ages                                  , metallicities
    integer                                                                                   :: agesCount                             , metallicitiesCount
    !![
    <constructorAssign variables="component, frame, wavelengthMinimum, wavelengthMaximum, resolution, toleranceRelative, *stellarPopulationSpectra_, *stellarPopulationSpectraPostprocessor_, *starFormationHistory_, *outputTimes_, *cosmologyFunctions_"/>
    !!]
    
    if     (                                                                                                                          &
         &   component /= componentTypeDisk                                                                                           &
         &  .and.                                                                                                                     &
         &   component /= componentTypeSpheroid                                                                                       &
         &  .and.                                                                                                                     &
         &   component /= componentTypeNuclearStarCluster                                                                             &
         & ) call Error_Report("only 'disk', 'spheroid' and 'nuclearStarCluster' components are supported"//{introspection:location})
    call self%stellarPopulationSpectra_%wavelengths(self%countWavelengths                   ,self%wavelengths_              )
    call self%stellarPopulationSpectra_%tabulation (     agesCount       ,metallicitiesCount,     ages        ,metallicities)    
    self%metallicityBoundaries       =self%starFormationHistory_%metallicityBoundaries()
    self%agePopulationMaximum        =ages         (agesCount         )
    self%metallicityPopulationMaximum=metallicities(metallicitiesCount)/metallicitySolar
    self%metallicityPopulationMinimum=metallicities(                 1)/metallicitySolar
    self%abundanceIndex              =Abundance_Pattern_Lookup(abundanceName="solar")
    ! Compute the factor by which the minimum/maximum wavelength in a resolution element differ from the central wavelength.
    if (resolution > 0.0d0) self%factorWavelength=(1.0d0+sqrt(1.0d0+4.0d0*resolution**2))/2.0d0/resolution
    return
  end function sedConstructorInternal

  subroutine sedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily sed} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSED), intent(inout) :: self
    
    !![
    <objectDestructor name="self%stellarPopulationSpectra_"             />
    <objectDestructor name="self%stellarPopulationSpectraPostprocessor_"/>
    <objectDestructor name="self%starFormationHistory_"                 />
    <objectDestructor name="self%outputTimes_"                          />
    <objectDestructor name="self%cosmologyFunctions_"                   />
    !!]
    return
  end subroutine sedDestructor

  integer function sedElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily sed} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorSED), intent(inout) :: self
    double precision                          , intent(in   ) :: time
    !$GLC attributes unused :: time

    sedElementCount=1
    return
  end function sedElementCount

  function sedSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily sed} property extractors.
    !!}
    use :: Error                         , only : Error_Report
    use :: Stellar_Luminosities_Structure, only : frameRest   , frameObserved 
    implicit none
    integer         (c_size_t                )                              :: sedSize
    class           (nodePropertyExtractorSED), intent(inout)               :: self
    double precision                          , intent(in   )               :: time
    logical                                   , allocatable  , dimension(:) :: selection
    integer         (c_size_t                )                              :: indexTemplate  , countTemplates, &
         &                                                                     i
    double precision                                                        :: expansionFactor

    if (self%resolution < 0.0d0) then
       ! Full resolution, so the number of wavelengths is simply the total number available within the wavelength range.
       allocate(selection(self%countWavelengths))
       select case (self%frame%ID)
       case (frameRest    %ID)
          expansionFactor=1.0d0
          selection      =                                                             &
               &           self%wavelengths_                 >= self%wavelengthMinimum &
               &          .and.                                                        &
               &           self%wavelengths_                 <= self%wavelengthMaximum
       case (frameObserved%ID)
          expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
          selection      =                                                             &
               &           self%wavelengths_/expansionFactor >= self%wavelengthMinimum &
               &          .and.                                                        &
               &           self%wavelengths_/expansionFactor <= self%wavelengthMaximum
       case default
          expansionFactor=1.0d0
          sedSize        =0_c_size_t
          call Error_Report('unknown frame'//{introspection:location})
       end select
       sedSize      =count(selection)
       indexTemplate=self%indexTemplateTime(time,countTemplates)
       if (indexTemplate > 0) then
          if (.not.allocated(self%templates)) allocate(self%templates(countTemplates))
          if (self%templates(indexTemplate)%countWavelengths < 0_c_size_t) then
             self%templates(indexTemplate)%countWavelengths=sedSize
             allocate(self%templates(indexTemplate)%wavelength(sedSize))
             self%templates(indexTemplate)%wavelength=+pack(self%wavelengths_   ,mask=selection) &
                  &                                   /          expansionFactor
          end if
       end if
    else
       ! The number of wavelengths must be computed, or retrieved from a template.
       indexTemplate=self%indexTemplateTime(time,countTemplates)
       if (indexTemplate > 0) then
          ! A template can be used. If the result is already computed for this template, use it.
          if (.not.allocated(self%templates)) allocate(self%templates(countTemplates))
          if (self%templates(indexTemplate)%countWavelengths > -1_c_size_t) then
             sedSize=self%templates(indexTemplate)%countWavelengths
             return
          end if
       end if
       ! Compute the number of wavelengths.
       sedSize=int(log(self%wavelengthMaximum/self%wavelengthMinimum)/log(self%factorWavelength)/2.0d0)+1
       ! Store the result for future re-use if possible.
       if (indexTemplate > 0) then
          self%templates(indexTemplate)%countWavelengths=sedSize
          allocate(self%templates(indexTemplate)%wavelength(sedSize))
          self%templates(indexTemplate)%wavelength(1)=+self%wavelengthMinimum &
               &                                      *self%factorWavelength
          if (sedSize > 1_c_size_t) then
             do i=2_c_size_t,sedSize
                self%templates(indexTemplate)%wavelength(i)=+self%templates(indexTemplate)%wavelength      (i-1)    &
                     &                                      *self                         %factorWavelength     **2
             end do
          end if
       end if
    end if
    return
  end function sedSize

  function sedExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily sed} property extractor.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, componentTypeNuclearStarCluster
    use :: Histories                 , only : history
    implicit none
    double precision                          , dimension(:,:  )          , allocatable :: sedExtract
    class           (nodePropertyExtractorSED), intent(inout)   , target                :: self
    type            (treeNode                ), intent(inout)   , target                :: node
    double precision                          , intent(in   )                           :: time
    type            (multiCounter            ), intent(inout)   , optional              :: instance
    class           (nodeComponentDisk       )                  , pointer               :: disk
    class           (nodeComponentSpheroid   )                  , pointer               :: spheroid
    class           (nodeComponentNSC        )                  , pointer               :: nuclearStarCluster
    double precision                          , dimension(:,:,:), pointer               :: sedTemplate_
    double precision                          , dimension(:,:,:), target  , allocatable :: sedTemplate
    double precision                          , dimension(  :,:)          , allocatable :: masses
    type            (history                 )                                          :: starFormationHistory
    integer         (c_size_t                )                                          :: countTemplates
    integer                                                                             :: indexTemplate       , iWavelength
    !$GLC attributes unused :: instance

    allocate(sedExtract(self%size(time),1))
    sedExtract=0.0d0
    ! Get the relevant star formation history.
    select case (self%component%ID)
    case (componentTypeDisk               %ID)
       disk                 => node              %disk                ()
       starFormationHistory =  disk              %starFormationHistory()
    case (componentTypeSpheroid           %ID)
       spheroid             => node              %spheroid            ()
       starFormationHistory =  spheroid          %starFormationHistory()
     case (componentTypeNuclearStarCluster%ID)
       nuclearStarCluster   => node              %NSC                 ()
       starFormationHistory =  nuclearStarCluster%starFormationHistory()
    end select
    if (.not.starFormationHistory%exists()) return
    ! Get the index of the template to use.
    indexTemplate=self%indexTemplateNode(node,starFormationHistory,countTemplates)
    if (indexTemplate > 0) then
       ! Stored templates can be used, so point to the relevant set.
       sedTemplate_ => self%templates(indexTemplate)%sed
    else
       ! Stored templates can not be used, get the templates for this specific case, and point to them.
       sedTemplate  =  self%luminosityMean(time,node,starFormationHistory)
       sedTemplate_ => sedTemplate
    end if
    masses=self%starFormationHistory_%masses(node,starFormationHistory,allowTruncation=.false.)
    do iWavelength=1,size(sedExtract,dim=1)
       sedExtract(iWavelength,1)=sum(sedTemplate_(iWavelength,:,:)*masses(:,:))
    end do
    return
  end function sedExtract

  subroutine sedNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily sed} properties.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class           (nodePropertyExtractorSED), intent(inout)                             :: self
    double precision                          , intent(in   ), optional                   :: time
    type            (varying_string          ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(1))
    names(1)=enumerationComponentTypeDecode(self%component,includePrefix=.false.)//"StellarSED:"//self%stellarPopulationSpectraPostprocessor_%objectType(short=.true.)
    return
  end subroutine sedNames

  subroutine sedDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily sed} property.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    class           (nodePropertyExtractorSED), intent(inout)                             :: self
    double precision                          , intent(in   ), optional                   :: time
    type            (varying_string          ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(1))
    descriptions(1)="Spectral energy density (SED), dL/dν for the "//enumerationComponentTypeDecode(self%component,includePrefix=.false.)//" [L☉ Hz⁻¹]."
    return
  end subroutine sedDescriptions

  function sedWavelengths(self,time)
    !!{
    Return wavelengths at which the SED is tabulated.
    !!}
    use :: Error                         , only : Error_Report
    use :: Stellar_Luminosities_Structure, only : frameRest   , frameObserved
    implicit none
    double precision                          , dimension(:) , allocatable :: sedWavelengths
    class           (nodePropertyExtractorSED), intent(inout)              :: self
    double precision                          , intent(in   )              :: time
    integer         (c_size_t                )                             :: i             , j              , &
         &                                                                    indexTemplate , countTemplates
    double precision                                                       :: wavelength    , expansionFactor

    allocate(sedWavelengths(self%size(time)))
    indexTemplate  =self%indexTemplateTime(time,countTemplates)
    j              =0
    select case (self%frame%ID)
    case (frameRest    %ID)
       expansionFactor=1.0d0
    case (frameObserved%ID)
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
    case default
       expansionFactor=0.0d0
       call Error_Report('unknown frame'//{introspection:location})
    end select
    do i=1,size(sedWavelengths)
       if (indexTemplate > 0) then
          ! Use wavelength from the pre-computed template.
          wavelength=self%templates(indexTemplate)%wavelength(i)
       else if (self%resolution < 0.0d0) then
          ! Full resolution SED.
          j=j+1
          do while (self%wavelengths_(j)/expansionFactor < self%wavelengthMinimum)
             j=j+1
          end do
          wavelength=self%wavelengths_(j)/expansionFactor
       else
          ! Finite resolution SED.
          if (i == 1) then
             wavelength=+self%wavelengthMinimum    &
                  &     *self%factorWavelength
          else
             wavelength=+     wavelength           &
                  &     *self%factorWavelength **2
          end if
       end if
       sedWavelengths(i)=wavelength
    end do
    return
  end function sedWavelengths

  subroutine sedColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily sed} property.
    !!}
    use :: Numerical_Constants_Units, only : metersToAngstroms
    implicit none
    class           (nodePropertyExtractorSED), intent(inout)                            :: self
    double precision                          , intent(in   ), optional                  :: time
    type            (varying_string          ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                          , intent(inout), dimension(:), allocatable :: values 
    type            (varying_string          ), intent(  out)                            :: valuesDescription
    double precision                          , intent(  out)                            :: valuesUnitsInSI
    integer         (c_size_t                )                                           :: i
    character       (len=18                  )                                           :: label
    
    allocate(descriptions(self%size(time)))
    allocate(values      (self%size(time)))
    values=self%wavelengths(time)
    do i=1,size(descriptions)      
       write (label,'(a2,1x,e12.6,1x,a1)') "λ=",values(i),"Å"
       descriptions(i)=trim(label)
    end do
    valuesDescription=var_str('Wavelengths at which the SED is tabulated [in units of Å].')
    valuesUnitsInSI  =1.0d0/metersToAngstroms
    return
  end subroutine sedColumnDescriptions

  function sedUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily sed} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    implicit none
    double precision                          , allocatable  , dimension(:) :: sedUnitsInSI
    class           (nodePropertyExtractorSED), intent(inout)               :: self
    double precision                          , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(sedUnitsInSI(1))
    sedUnitsInSI(1)=luminositySolar
    return
  end function sedUnitsInSI

  integer function sedIndexTemplateTime(self,time,countTemplates) result(indexTemplate)
    !!{
    Find the index of the template SEDs to use.
    !!}
    use :: Numerical_Comparison          , only : Values_Agree
    use :: Star_Formation_Histories      , only : starFormationHistoryAgesFixed, starFormationHistoryAgesFixedPerOutput
    use :: Stellar_Luminosities_Structure, only : frameRest
    implicit none
    class           (nodePropertyExtractorSED), intent(inout) :: self
    double precision                          , intent(in   ) :: time
    integer         (c_size_t                ), intent(  out) :: countTemplates 
    integer         (c_size_t                )                :: indexOutput

    if      (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed         ) then
       ! Ages are fixed. A single template can be used if we are computing rest-frame SEDs, and postprocessing is redshift independent.
       if (self%frame == frameRest .and. .not.self%stellarPopulationSpectraPostprocessor_%isRedshiftDependent()) then
          indexTemplate =+1
          countTemplates=+1
       else
          indexTemplate =-1
          countTemplates=-1
       end if
    else if (self%starFormationHistory_%ageDistribution() ==starFormationHistoryAgesFixedPerOutput) then
       ! Check that the time is an output time.
       indexOutput=self%outputTimes_%index(time,findClosest=.true.)
       if (Values_Agree(time,self%outputTimes_%time(indexOutput),relTol=1.0d-6)) then
          ! The time corresponds to an output time - use a template.
          indexTemplate=int(indexOutput)
       else
          ! The time does not correspond to an output time - a template can not be used.
          indexTemplate=-1
       end if
       countTemplates=self%outputTimes_%count()
    else
       ! Templates can not be used.
       indexTemplate =-1
       countTemplates=-1
    end if
    return
  end function sedIndexTemplateTime

  integer function sedIndexTemplateNode(self,node,starFormationHistory,countTemplates) result(indexTemplate)
    !!{
    Find the index of the template SEDs to use, and also compute the template.
    !!}
    use :: Display                       , only : displayMessage               , verbosityLevelWorking
    use :: Galacticus_Nodes              , only : nodeComponentBasic
    use :: Histories                     , only : history
    use :: ISO_Varying_String            , only : var_str
    use :: HDF5_Access                   , only : hdf5Access
    use :: IO_HDF5                       , only : hdf5Object
    use :: Numerical_Comparison          , only : Values_Agree
    use :: File_Utilities                , only : File_Exists                  , File_Lock                             , File_Unlock, lockDescriptor
    use :: String_Handling               , only : operator(//)
    use :: Star_Formation_Histories      , only : starFormationHistoryAgesFixed, starFormationHistoryAgesFixedPerOutput
    use :: Stellar_Luminosities_Structure, only : frameRest
    use :: Input_Paths                   , only : inputPath                    , pathTypeDataDynamic
    implicit none
    class    (nodePropertyExtractorSED), intent(inout) :: self
    type     (treeNode                ), intent(inout) :: node
    type     (history                 ), intent(in   ) :: starFormationHistory
    integer  (c_size_t                ), intent(  out) :: countTemplates 
    class    (nodeComponentBasic      ), pointer       :: basic
    integer  (c_size_t                )                :: indexOutput
    type     (lockDescriptor          )                :: fileLock
    type     (hdf5Object              )                :: file
    type     (varying_string          )                :: fileName
    character(len=16                  )                :: label

    if      (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed         ) then
       ! Ages are fixed. A single template can be used if we are computing rest-frame SEDs, and postprocessing is redshift independent.
       if (self%frame == frameRest .and. .not.self%stellarPopulationSpectraPostprocessor_%isRedshiftDependent()) then
          indexTemplate =+1
          countTemplates=+1
       else
          indexTemplate =-1
          countTemplates=-1
       end if
    else if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixedPerOutput) then
       ! Check that the node exists at at output time.
       basic       => node             %basic(                               )
       indexOutput =  self%outputTimes_%index(basic%time(),findClosest=.true.)
       if (Values_Agree(basic%time(),self%outputTimes_%time(indexOutput),relTol=1.0d-6)) then
          ! The time corresponds to an output time - use a template.
          indexTemplate=int(indexOutput)
       else
          ! The time does not correspond to an output time - a template can not be used.
          indexTemplate=-1
       end if
       countTemplates=self%outputTimes_%count()
    else
       ! Templates can not be used.
       indexTemplate =-1
       countTemplates=-1
    end if
    ! Return if a template can not be used.
    if (indexTemplate < 0) return
    ! Ensure that the templates have been built for this index.
    if (.not.allocated(self%templates)) allocate(self%templates(countTemplates))
    if (.not.allocated(self%templates(indexTemplate)%sed)) then
       ! Construct the file name.
       fileName=inputPath(pathTypeDataDynamic)                              // &
            &        'stellarPopulations/'                                  // &
            &        self%objectType             (                         )// &
            &        '_'                                                    // &
            &        self%historyHashedDescriptor(node,starFormationHistory)// &
            &        '_'                                                    // &
            &        indexTemplate                                          // &
            &        '.hdf5'
       ! Check if the templates can be retrieved from file.
       !! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
       if (File_Exists(fileName)) then
          !$ call hdf5Access%set()
          call file%openFile(char(fileName))
          if (file%hasDataset('sedTemplate')) then
             if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed) then
                call displayMessage("reading SED tabulation from file '"                                        //fileName//"'",verbosityLevelWorking)
             else
                !$omp critical(gfortranInternalIO)
                write (label,'(f12.8)') self%outputTimes_%time(indexOutput)
                !$omp end critical(gfortranInternalIO)
                call displayMessage("reading SED tabulation for time "//trim(adjustl(label))//" Gyr from file '"//fileName//"'",verbosityLevelWorking)
             end if
             call file%readDataset('sedTemplate',self%templates(indexTemplate)%sed)
          end if
          call file%close()
          !$ call hdf5Access%unset()
       end if
       if (.not.allocated(self%templates(indexTemplate)%sed)) then
          basic                              => node%basic         (                                                         )
          self %templates(indexTemplate)%sed =  self%luminosityMean(basic%time(),node,starFormationHistory,parallelize=.true.)
          if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed) then
             call displayMessage("storing SED tabulation to file '"                                        //fileName//"'",verbosityLevelWorking)
          else
             !$omp critical(gfortranInternalIO)
             write (label,'(f12.8)') self%outputTimes_%time(indexOutput)
             !$omp end critical(gfortranInternalIO)
             call displayMessage("storing SED tabulation for time "//trim(adjustl(label))//" Gyr to file '"//fileName//"'",verbosityLevelWorking)
          end if
          !$ call hdf5Access%set()
          call file%openFile(char(fileName),overWrite=.false.,readOnly=.false.)
          call file%writeDataset(self%templates(indexTemplate)%sed,'sedTemplate')
          call file%close()
          !$ call hdf5Access%unset()
       end if
       call File_Unlock(fileLock)
    end if
    return
  end function sedIndexTemplateNode

  double precision function sedLuminosityMean(self,time,node,starFormationHistory,parallelize)
    !!{
    Compute the mean luminosity of the stellar population in each bin of the star formation history.
    !!}
    use :: Abundances_Structure          , only : abundances                   , metallicityTypeLinearByMassSolar, adjustElementsReset
    use :: Display                       , only : displayIndent                , displayUnindent                 , displayCounter     , displayCounterClear, &
         &                                        verbosityLevelWorking        , displayMessage
    use :: Error                         , only : Error_Report
    use :: Histories                     , only : history
    use :: Numerical_Integration         , only : integrator
    use :: Numerical_Constants_Prefixes  , only : siFormat
    use :: Multi_Counters                , only : multiCounter
    use :: Locks                         , only : ompLock
    use :: Stellar_Luminosities_Structure, only : frameRest                    , frameObserved
    use :: Star_Formation_Histories      , only : starFormationHistoryAgesFixed
    use :: Timers                        , only : timer
    implicit none
    double precision                                            , dimension(:,:,:), allocatable :: sedLuminosityMean
    class           (nodePropertyExtractorSED                  ), intent(inout)                 :: self
    double precision                                            , intent(in   )                 :: time
    type            (treeNode                                  ), intent(inout)                 :: node
    type            (history                                   ), intent(in   )                 :: starFormationHistory
    logical                                                     , intent(in   )   , optional    :: parallelize
    class           (stellarPopulationSpectraClass             ), pointer         , save        :: stellarPopulationSpectra_
    class           (stellarPopulationSpectraPostprocessorClass), pointer         , save        :: stellarPopulationSpectraPostprocessor_
    class           (cosmologyFunctionsClass                   ), pointer         , save        :: cosmologyFunctions_
    type            (integrator                                ), allocatable     , save        :: integratorTime                        , integratorMetallicity, &
         &                                                                                         integratorWavelength
    integer         (c_size_t                                  ), dimension(:    ), allocatable :: jWavelength
    double precision                                            , dimension(:    ), allocatable :: wavelengthMinima                      , wavelengthMaxima     , &
         &                                                                                         times
    double precision                                            , dimension(  :,:), allocatable :: masses
    integer         (c_size_t                                  )                                :: iWavelength                           , iTime                , &
         &                                                                                         iMetallicity                          , kWavelength          , &
         &                                                                                         counter                               , counterMaximum       , &
         &                                                                                         iterator
    double precision                                                                            :: metallicityMinimum                    , metallicityMaximum   , &
         &                                                                                         expansionFactor                       , timeStart
    double precision                                                              , save        :: timeMinimum                           , timeMaximum          , &
         &                                                                                         wavelength                            , wavelengthMinimum    , &
         &                                                                                         wavelengthMaximum                     , age                  , &
         &                                                                                         redshift
    type            (abundances                                )                  , save        :: abundancesStellar
    character       (len=12                                    )                                :: label
    type            (multiCounter                              )                                :: state
    type            (ompLock                                   )                                :: stateLock
    type            (timer                                     )                                :: timer_
    !$omp threadprivate(stellarPopulationSpectra_,stellarPopulationSpectraPostprocessor_,cosmologyFunctions_,integratorTime,integratorWavelength,integratorMetallicity,abundancesStellar,wavelength,wavelengthMinimum,wavelengthMaximum,timeMinimum,timeMaximum,age,redshift)
    !$GLC attributes initialized :: masses
    !![
    <optionalArgument name="parallelize" defaultsTo=".false." />
    !!]

    times =self%starFormationHistory_%times (node=node,starFormationHistory=starFormationHistory,allowTruncation=.false.,timeStart=timeStart)
    masses=self%starFormationHistory_%masses(node=node,starFormationHistory=starFormationHistory,allowTruncation=.false.                    )
    allocate(sedLuminosityMean(self%size(time),size(masses,dim=1),size(masses,dim=2)))
    select case (self%frame%ID)
    case (frameRest    %ID)
       expansionFactor=1.0d0
    case (frameObserved%ID)
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
    case default
       expansionFactor=0.0d0
       call Error_Report('unknown frame'//{introspection:location})
    end select
    if (self%resolution < 0.0d0) then
       allocate(jWavelength(size(sedLuminosityMean,dim=1)))
       kWavelength=0
       do iWavelength=1,self%countWavelengths
          if     (                                                                          &
               &   self%wavelengths_(iWavelength)/expansionFactor >= self%wavelengthMinimum &
               &  .and.                                                                     &
               &   self%wavelengths_(iWavelength)/expansionFactor <= self%wavelengthMaximum &
               & ) then
             kWavelength             =kWavelength+1
             jWavelength(kWavelength)=iWavelength
          end if
       end do
    else
       allocate(wavelengthMinima(size(sedLuminosityMean,dim=1)))
       allocate(wavelengthMaxima(size(sedLuminosityMean,dim=1)))
       wavelengthMinima(1)=self%wavelengthMinimum
       wavelengthMaxima(1)=wavelengthMinima(1)*self%factorWavelength**2
       do iWavelength=2,size(wavelengthMinima)
          wavelengthMinima(iWavelength)=wavelengthMaxima(iWavelength-1)
          wavelengthMaxima(iWavelength)=wavelengthMinima(iWavelength  )*self%factorWavelength**2
       end do
    end if
    counter       =-1
    counterMaximum=product     ([size(sedLuminosityMean,dim=1              ),size(masses,dim=1              ),size(masses,dim=2              )])
    state         =multiCounter([size(sedLuminosityMean,dim=1,kind=c_size_t),size(masses,dim=1,kind=c_size_t),size(masses,dim=2,kind=c_size_t)])
    stateLock     =ompLock     (                                                                                                               )
    !$omp parallel private (iWavelength,iTime,iMetallicity,metallicityMinimum,metallicityMaximum)
    allocate(integratorTime       )
    allocate(integratorMetallicity)
    allocate(integratorWavelength )
    integratorTime       =integrator(sedIntegrandTime       ,toleranceRelative=self%toleranceRelative)
    integratorMetallicity=integrator(sedIntegrandMetallicity,toleranceRelative=self%toleranceRelative)
    integratorWavelength =integrator(sedIntegrandWavelength ,toleranceRelative=self%toleranceRelative)
    if (parallelize_) then
       allocate(stellarPopulationSpectra_             ,mold=self%stellarPopulationSpectra_             )
       allocate(stellarPopulationSpectraPostprocessor_,mold=self%stellarPopulationSpectraPostprocessor_)
       allocate(cosmologyFunctions_                   ,mold=self%cosmologyFunctions_                   )
       !$omp critical(nodePropertyExtractSEDDeepCopy)
       !![
       <deepCopyReset variables="self%stellarPopulationSpectra_ self%stellarPopulationSpectraPostprocessor_ self%cosmologyFunctions_"/>
       !!]
       !![
       <deepCopy source="self%stellarPopulationSpectra_"              destination="stellarPopulationSpectra_"             />
       <deepCopy source="self%stellarPopulationSpectraPostprocessor_" destination="stellarPopulationSpectraPostprocessor_"/>
       <deepCopy source="self%cosmologyFunctions_"                    destination="cosmologyFunctions_"                   />
       !!]
       !![
       <deepCopyFinalize variables="stellarPopulationSpectra_ stellarPopulationSpectraPostprocessor_ cosmologyFunctions_"/>
       !!]
       !$omp end critical(nodePropertyExtractSEDDeepCopy)
    end if
    !$omp master
    if (parallelize_) then
       if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed) then
          call displayIndent("computing template SEDs"                                        ,verbosityLevelWorking)
       else
          !$omp critical(gfortranInternalIO)
          write (label,'(f12.8)') time
          !$omp end critical(gfortranInternalIO)
          call displayIndent("computing template SEDs for time "//trim(adjustl(label))//" Gyr",verbosityLevelWorking)
       end if
       call timer_%start()
    end if
    !$omp end master
    ! Iterate over (wavelength,time,metallicity).
    !$omp do
    do iterator=0,counterMaximum-1
       call stateLock%  set()
       if (state%increment()) then
          iWavelength =state%state(1_c_size_t)
          iTime       =state%state(2_c_size_t)
          iMetallicity=state%state(3_c_size_t)
       else
          iWavelength =0_c_size_t
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
       ! Determine the wavelength.
       if (self%resolution < 0.0d0) then
          ! Full resolution SED.
          wavelength=self%wavelengths_(jWavelength(iWavelength))
       else
          ! Finite resolution SED.
          wavelengthMinimum=wavelengthMinima(iWavelength)*expansionFactor
          wavelengthMaximum=wavelengthMaxima(iWavelength)*expansionFactor
       end if
       ! Determine times.
       if (iTime == 1) then
          timeMinimum=    timeStart
       else
          timeMinimum=    times(iTime-1)
       end if
       timeMaximum   =min(times(iTime  ),time)
       if (timeMaximum <= timeMinimum) then
          sedLuminosityMean(iWavelength,iTime,iMetallicity)=0.0d0
          cycle
       end if
       ! Determine metallicities.
       if (iMetallicity == 1) then
          metallicityMinimum=                                                   self%metallicityPopulationMinimum
       else
          metallicityMinimum=min(max(self%metallicityBoundaries(iMetallicity-1),self%metallicityPopulationMinimum),self%metallicityPopulationMaximum)
       end if
       metallicityMaximum   =min(max(self%metallicityBoundaries(iMetallicity  ),self%metallicityPopulationMinimum),self%metallicityPopulationMaximum)
       if (metallicityMaximum > metallicityMinimum) then
          sedLuminosityMean(iWavelength,iTime,iMetallicity)=+integratorMetallicity%integrate(metallicityMinimum,metallicityMaximum) &
               &                                            /                               (timeMaximum       -timeMinimum       ) &
               &                                            /                               (metallicityMaximum-metallicityMinimum)
       else
          call abundancesStellar%metallicitySet(                                                       &
               &                                metallicity    =     metallicityMinimum              , &
               &                                metallicityType=     metallicityTypeLinearByMassSolar, &
               &                                adjustElements =     adjustElementsReset             , &
               &                                abundanceIndex =self%abundanceIndex                    &
               &                               )
          sedLuminosityMean(iWavelength,iTime,iMetallicity)=+integratorTime       %integrate(timeMinimum       ,timeMaximum       ) &
               &                                            /                               (timeMaximum       -timeMinimum       )
       end if
    end do
    !$omp end do
    !$omp master
    if (parallelize_) then
       call timer_%stop()
       call displayCounterClear(                                                                                          verbosityLevelWorking)
       call displayMessage     ("table is "//trim(adjustl(siFormat(dble(sizeof(sedLuminosityMean)),'f9.3')))//"B in size",verbosityLevelWorking)
       call displayUnindent    ("done in " //trim(adjustl(timer_%reportText()                             ))             ,verbosityLevelWorking)
    end if
    !$omp end master
    !![
    <objectDestructor name="stellarPopulationSpectra_"             />
    <objectDestructor name="stellarPopulationSpectraPostprocessor_"/>
    <objectDestructor name="cosmologyFunctions_"                   />
    !!]
    deallocate(integratorTime       )
    deallocate(integratorMetallicity)
    deallocate(integratorWavelength )
    !$omp end parallel
    return

  contains

    double precision function sedIntegrandMetallicity(metallicity)
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
      sedIntegrandMetallicity=integratorTime%integrate(timeMinimum,timeMaximum)
      return
    end function sedIntegrandMetallicity

    double precision function sedIntegrandTime(timeBirth)
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
      if (self%resolution < 0.0d0) then
         ! Full resolution - evaluate at the given wavelength.
         if (parallelize_) then
            sedIntegrandTime=+     stellarPopulationSpectra_             %luminosity(abundancesStellar,age,wavelength) &
                 &           *     stellarPopulationSpectraPostprocessor_%multiplier(wavelength       ,age,redshift  )
         else
            sedIntegrandTime=+self%stellarPopulationSpectra_             %luminosity(abundancesStellar,age,wavelength) &
                 &           *self%stellarPopulationSpectraPostprocessor_%multiplier(wavelength       ,age,redshift  )
         end if
      else
         ! Finite resolution - integrate over wavelength.
         sedIntegrandTime=integratorWavelength%integrate(wavelengthMinimum,wavelengthMaximum)
      end if
      return
    end function sedIntegrandTime

    double precision function sedIntegrandWavelength(wavelength)
      !!{
      Integrand over wavelength of the stellar population. The assumption made here is that we have a photon-counting detector
      (such as a CCD). So, we integrate the photon rate over the wavelength range, and then multiply by the energy at the
      central wavelength $\bar{\lambda} = \mathrm{c}/\bar{\nu}$, to find the luminosity in the wavelength range. We then divide
      by the width of the range in frequency to get our SED. Specifically,
      \begin{equation}
        \langle L_\nu \rangle = \frac{\mathrm{h} \bar{\nu}}{\Delta \nu} \int_{\lambda_\mathrm{min}}^{\lambda_\mathrm{max}} \mathrm{d}\nu \frac{L_\nu}{\mathrm{h}\nu}.
      \end{equation}
      Using the fact that $\Delta\nu = \nu_1-\nu_2 = (\mathrm{c}/\bar{\lambda})(f-f^{-1})$ this can be written as
      \begin{equation}
        \langle L_\nu \rangle = (f-f^{-1})^{-1} \int_{\lambda_\mathrm{min}}^{\lambda_\mathrm{max}} \frac{\mathrm{d}\lambda}{\lambda} L_\nu.
      \end{equation}
      !!}
      implicit none
      double precision, intent(in   ) :: wavelength

      if (parallelize_) then
         sedIntegrandWavelength=+     stellarPopulationSpectra_             %luminosity(abundancesStellar,age,wavelength) &
                 &              *     stellarPopulationSpectraPostprocessor_%multiplier(wavelength       ,age,redshift  )
      else
         sedIntegrandWavelength=+self%stellarPopulationSpectra_             %luminosity(abundancesStellar,age,wavelength) &
                 &              *self%stellarPopulationSpectraPostprocessor_%multiplier(wavelength       ,age,redshift  )
      end if
      sedIntegrandWavelength=+sedIntegrandWavelength        &
           &                 /wavelength                    &
           &                 /(                             &
           &                   +1.0d0*self%factorWavelength &
           &                   -1.0d0/self%factorWavelength &
           &                 )
      return
    end function sedIntegrandWavelength

  end function sedLuminosityMean

  function sedHistoryHashedDescriptor(self,node,starFormationHistory)  
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Galacticus_Nodes              , only : treeNode
    use :: Input_Parameters              , only : inputParameters
    use :: String_Handling               , only : String_C_To_Fortran
    use :: Hashes_Cryptographic          , only : Hash_MD5
    use :: FoX_DOM                       , only : setLiveNodeLists
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : enumerationFrameDecode
    use :: Star_Formation_Histories      , only : starFormationHistoryAgesFixed
    implicit none
    type            (varying_string          )                              :: sedHistoryHashedDescriptor
    class           (nodePropertyExtractorSED), intent(in   )               :: self
    type            (treeNode                ), intent(inout)               :: node
    type            (history                 ), intent(in   )               :: starFormationHistory
    double precision                          , allocatable  , dimension(:) :: times
    character       (len=18                  )                              :: parameterLabel
    type            (inputParameters         )                              :: descriptor
    type            (varying_string          )                              :: descriptorString          , values
    integer                                                                 :: i
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
    !$omp critical(gfortranInternalIO)
    write (parameterLabel,'(e17.10)') self%wavelengthMinimum
    !$omp end critical(gfortranInternalIO)
    call descriptor%addParameter('wavelengthMinimum',trim(adjustl(parameterLabel)))
    !$omp critical(gfortranInternalIO)
    write (parameterLabel,'(e17.10)') self%wavelengthMaximum
    !$omp end critical(gfortranInternalIO)
    call descriptor%addParameter('wavelengthMaximum',trim(adjustl(parameterLabel)))
    !$omp critical(gfortranInternalIO)
    write (parameterLabel,'(e17.10)') self%resolution
    !$omp end critical(gfortranInternalIO)
    call descriptor%addParameter('resolution'       ,trim(adjustl(parameterLabel)))
    call self%stellarPopulationSpectra_             %descriptor(descriptor)
    call self%stellarPopulationSpectraPostprocessor_%descriptor(descriptor)
    call self%starFormationHistory_                 %descriptor(descriptor)
    call self%outputTimes_                          %descriptor(descriptor)
    call self%cosmologyFunctions_                   %descriptor(descriptor)
    values=""
    do i=1,size(self%metallicityBoundaries)
       !$omp critical(gfortranInternalIO)
       write (parameterLabel,'(e17.10)') self                %metallicityBoundaries(i)
       !$omp end critical(gfortranInternalIO)
       values=values//trim(adjustl(parameterLabel))
       if (i < size(self%metallicityBoundaries)) values=values//":"
    end do
    call descriptor%addParameter('metallicity',char(values))
    ! Times are only added if ages are not fixed. For fixed ages, the history is the same (for our purposes) always.
    if (self%starFormationHistory_%ageDistribution() /= starFormationHistoryAgesFixed) then
       values=""
       times =self %starFormationHistory_%times(node=node,starFormationHistory=starFormationHistory)
       do i=1,size(times)
          !$omp critical(gfortranInternalIO)
          write (parameterLabel,'(e17.10)') times(i)
          !$omp end critical(gfortranInternalIO)
          values=values//trim(adjustl(parameterLabel))
          if (i < size(times)) values=values//":"
       end do
       call descriptor%addParameter('time',char(values))
    end if
    descriptorString=descriptor%serializeToString()
    call descriptor%destroy()
    descriptorString=descriptorString//":sourceDigest{"//String_C_To_Fortran(nodePropertyExtractorSED5)//"}"
    !![
    <workaround type="gfortran" PR="102845" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=102845">
     <description>
      Memory leak possibly due to OpenMP parallelism, or some failing of gfortran.
     </description>
    </workaround>
    !!]
    !if (descriptorString /= descriptorStringPrevious) then
    !   descriptorStringPrevious=         descriptorString
    !   hashedDescriptorPrevious=Hash_MD5(descriptorString)
    !end if
    !sedHistoryHashedDescriptor=hashedDescriptorPrevious
    sedHistoryHashedDescriptor=Hash_MD5(descriptorString)
    return
  end function sedHistoryHashedDescriptor

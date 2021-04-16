!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% Contains a module which implements a property extractor class for the SED of a component.
  use :: Output_Times              , only : outputTimesClass
  use :: Stellar_Population_Spectra, only : stellarPopulationSpectraClass
  use :: Star_Formation_Histories  , only : starFormationHistoryClass

  type :: sedTemplate
     !% Type used to store SED templates.
     private
     double precision, allocatable, dimension(:,:,:) :: sed
  end type sedTemplate
  
  !# <nodePropertyExtractor name="nodePropertyExtractorSED">
  !#  <description>A property extractor class for the SED of a component.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorSED
     !% A property extractor class for the SED of a component.
     private
     class           (stellarPopulationSpectraClass), pointer                   :: stellarPopulationSpectra_    => null()
     class           (starFormationHistoryClass    ), pointer                   :: starFormationHistory_        => null()
     class           (outputTimesClass             ), pointer                   :: outputTimes_                 => null()
     integer                                                                    :: countWavelengths                      , component
     double precision                               , allocatable, dimension(:) :: wavelengths                           , metallicityBoundaries
     type            (sedTemplate                  ), allocatable, dimension(:) :: templates
     double precision                                                           :: metallicityPopulationMaximum          , agePopulationMaximum , &
          &                                                                        metallicityPopulationMinimum
     integer                                                                    :: abundanceIndex
     logical                                                                    :: useSEDTemplates
     type            (varying_string               )                            :: fileName
   contains
     !# <methods>
     !#  <method description="Return the index of the template SEDs to use."                                                         method="indexTemplate" />
     !#  <method description="Compute the mean luminosity of the stellar population in the given bin of the star formation history." method="luminosityMean"/>
     !# </methods>
     final     ::                       sedDestructor
     procedure :: columnDescriptions => sedColumnDescriptions
     procedure :: size               => sedSize
     procedure :: elementCount       => sedElementCount
     procedure :: extract            => sedExtract
     procedure :: names              => sedNames
     procedure :: descriptions       => sedDescriptions
     procedure :: unitsInSI          => sedUnitsInSI
     procedure :: type               => sedType
     procedure :: indexTemplate      => sedIndexTemplate
     procedure :: luminosityMean     => sedLuminosityMean
  end type nodePropertyExtractorSED

  interface nodePropertyExtractorSED
     !% Constructors for the ``sed'' output analysis class.
     module procedure sedConstructorParameters
     module procedure sedConstructorInternal
  end interface nodePropertyExtractorSED

contains

  function sedConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sed} property extractor class which takes a parameter set as input.
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    implicit none
    type (nodePropertyExtractorSED     )                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(stellarPopulationSpectraClass), pointer       :: stellarPopulationSpectra_
    class(starFormationHistoryClass    ), pointer       :: starFormationHistory_
    class(outputTimesClass             ), pointer       :: outputTimes_
    type (varying_string               )                :: component

    !# <inputParameter>
    !#   <name>component</name>
    !#   <source>parameters</source>
    !#   <description>The component from which to extract star formation rate.</description>
    !# </inputParameter>
    !# <objectBuilder class="stellarPopulationSpectra" name="stellarPopulationSpectra_" source="parameters"/>
    !# <objectBuilder class="starFormationHistory"     name="starFormationHistory_"     source="parameters"/>
    !# <objectBuilder class="outputTimes"              name="outputTimes_"              source="parameters"/>
    self=nodePropertyExtractorSED(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),stellarPopulationSpectra_,starFormationHistory_,outputTimes_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="stellarPopulationSpectra_"/>
    !# <objectDestructor name="starFormationHistory_"    />
    !# <objectDestructor name="outputTimes_"             />
    return
  end function sedConstructorParameters

  function sedConstructorInternal(component,stellarPopulationSpectra_,starFormationHistory_,outputTimes_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sed} property extractor class.
    use :: Atomic_Data                     , only : Abundance_Pattern_Lookup
    use :: Galacticus_Paths                , only : galacticusPath          , pathTypeDataDynamic
    use :: Galactic_Structure_Options      , only : componentTypeDisk       , componentTypeSpheroid
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    implicit none
    type            (nodePropertyExtractorSED     )                              :: self
    integer                                        , intent(in   )               :: component
    class           (stellarPopulationSpectraClass), intent(in   ), target       :: stellarPopulationSpectra_
    class           (starFormationHistoryClass    ), intent(in   ), target       :: starFormationHistory_
    class           (outputTimesClass             ), intent(in   ), target       :: outputTimes_
    double precision                               , allocatable  , dimension(:) :: ages                     , metallicities
    integer                                                                      :: agesCount                , metallicitiesCount
    !# <constructorAssign variables="component, *stellarPopulationSpectra_, *starFormationHistory_, *outputTimes_"/>

    if     (                                                                                                               &
         &   component /= componentTypeDisk                                                                                &
         &  .and.                                                                                                          &
         &   component /= componentTypeSpheroid                                                                            &
         & ) call Galacticus_Error_Report("only 'disk' and 'spheroid' components are supported"//{introspection:location})
    call self%stellarPopulationSpectra_%wavelengths(self%countWavelengths                   ,self%wavelengths              )
    call self%stellarPopulationSpectra_%tabulation (     agesCount       ,metallicitiesCount,     ages       ,metallicities)
    self%metallicityBoundaries       =self%starFormationHistory_%metallicityBoundaries()
    self%agePopulationMaximum        =ages         (agesCount         )
    self%metallicityPopulationMaximum=metallicities(metallicitiesCount)/metallicitySolar
    self%metallicityPopulationMinimum=metallicities(                 1)/metallicitySolar
    self%abundanceIndex              =Abundance_Pattern_Lookup(abundanceName="solar")
    self%useSEDTemplates             =self%starFormationHistory_%perOutputTabualtionIsStatic()
    self%fileName                    =galacticusPath(pathTypeDataDynamic)              // &
         &                            'stellarPopulations/'                            // &
         &                            self%objectType      (                          )// &
         &                            '_'                                              // &
         &                            self%hashedDescriptor(includeSourceDigest=.true.)// &
         &                            '.hdf5'
    return
  end function sedConstructorInternal

  subroutine sedDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sed} property extractor class.
    implicit none
    type(nodePropertyExtractorSED), intent(inout) :: self

    !# <objectDestructor name="self%stellarPopulationSpectra_"/>
    !# <objectDestructor name="self%starFormationHistory_"    />
    !# <objectDestructor name="self%outputTimes_"             />
    return
  end subroutine sedDestructor

  integer function sedElementCount(self,time)
    !% Return the number of elements in the {\normalfont \ttfamily sed} property extractors.
    implicit none
    class           (nodePropertyExtractorSED), intent(inout) :: self
    double precision                          , intent(in   ) :: time
    !$GLC attributes unused :: time

    sedElementCount=1
    return
  end function sedElementCount

  function sedSize(self,time)
    !% Return the number of array alements in the {\normalfont \ttfamily sed} property extractors.
    implicit none
    integer         (c_size_t                )                :: sedSize
    class           (nodePropertyExtractorSED), intent(inout) :: self
    double precision                          , intent(in   ) :: time
    !$GLC attributes unused :: time

    sedSize=int(self%countWavelengths,kind=c_size_t)
    return
  end function sedSize

  function sedExtract(self,node,time,instance)
    !% Implement a {\normalfont \ttfamily sed} property extractor.
    use :: Galacticus_Nodes          , only : nodeComponentDisk, nodeComponentSpheroid
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid
    use :: Histories                 , only : history
    implicit none
    double precision                          , dimension(:,:  )          , allocatable :: sedExtract
    class           (nodePropertyExtractorSED), intent(inout)   , target                :: self
    type            (treeNode                ), intent(inout)   , target                :: node
    double precision                          , intent(in   )                           :: time
    type            (multiCounter            ), intent(inout)   , optional              :: instance
    class           (nodeComponentDisk       )                  , pointer               :: disk
    class           (nodeComponentSpheroid   )                  , pointer               :: spheroid
    double precision                          , dimension(:,:,:), pointer               :: sedTemplate_
    double precision                          , dimension(:,:,:), target  , allocatable :: sedTemplate
    type            (history                 )                                          :: starFormationHistory
    integer                                                                             :: indexTemplate       , iWavelength
    !$GLC attributes unused :: instance

    allocate(sedExtract(self%countWavelengths,1))
    sedExtract =  0.0d0
    ! Get the relevant star formation history.
    select case (self%component)
    case (componentTypeDisk)
       disk                 => node    %disk                ()
       starFormationHistory =  disk    %starFormationHistory()
    case (componentTypeSpheroid)
       spheroid             => node    %spheroid            ()
       starFormationHistory =  spheroid%starFormationHistory()
    end select
    if (.not.starFormationHistory%exists()) return
    ! Get the index of the template to use.
    indexTemplate=self%indexTemplate(node,starFormationHistory)
    if (indexTemplate > 0) then
       ! Stored templates can be used, so point to the relevant set.
       sedTemplate_ => self%templates(indexTemplate)%sed
    else
       ! Stored templates can not be used, get the templates for this specific case, and point to them.
       sedTemplate  =  self%luminosityMean(time,starFormationHistory)
       sedTemplate_ => sedTemplate
    end if
    do iWavelength=1,self%countWavelengths
       sedExtract(iWavelength,1)=sum(sedTemplate_(iWavelength,:,:)*starFormationHistory%data(:,:))
    end do
    return
  end function sedExtract

  function sedNames(self,time)
    !% Return the names of the {\normalfont \ttfamily sed} properties.
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    type            (varying_string          ), dimension(:) , allocatable :: sedNames
    class           (nodePropertyExtractorSED), intent(inout)              :: self
    double precision                          , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(sedNames(1))
    sedNames(1)=enumerationComponentTypeDecode(self%component,includePrefix=.false.)//"StellarSED"
    return
  end function sedNames

  function sedDescriptions(self,time)
    !% Return descriptions of the {\normalfont \ttfamily sed} property.
    use :: Galactic_Structure_Options, only : enumerationComponentTypeDecode
    implicit none
    type            (varying_string          ), dimension(:) , allocatable :: sedDescriptions
    class           (nodePropertyExtractorSED), intent(inout)              :: self
    double precision                          , intent(in   )              :: time
    !$GLC attributes unused :: time

    allocate(sedDescriptions(1))
    sedDescriptions(1)="Spectral energy density (SED) for the "//enumerationComponentTypeDecode(self%component,includePrefix=.false.)//" [L☉/Hz⁻¹]."
    return
  end function sedDescriptions

  function sedColumnDescriptions(self,time)
    !% Return column descriptions of the {\normalfont \ttfamily sed} property.
    implicit none
    type            (varying_string          ), dimension(:) , allocatable :: sedColumnDescriptions
    class           (nodePropertyExtractorSED), intent(inout)              :: self
    double precision                          , intent(in   )              :: time
    integer                                                                :: i
    character      (len=18                   )                             :: label
    !$GLC attributes unused :: time

    allocate(sedColumnDescriptions(self%countWavelengths))
    do i=1,self%countWavelengths
       write (label,'(a2,1x,e12.6,1x,a1)') "λ=",self%wavelengths(i),"Å"
       sedColumnDescriptions(i)=trim(label)
    end do
    return
  end function sedColumnDescriptions

  function sedUnitsInSI(self,time)
    !% Return the units of the {\normalfont \ttfamily sed} properties in the SI system.
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    implicit none
    double precision                          , allocatable  , dimension(:) :: sedUnitsInSI
    class           (nodePropertyExtractorSED), intent(inout)               :: self
    double precision                          , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(sedUnitsInSI(1))
    sedUnitsInSI(1)=luminositySolar
    return
  end function sedUnitsInSI

  integer function sedType(self)
    !% Return the type of the {\normalfont \ttfamily sed} properties.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorSED), intent(inout) :: self
    !$GLC attributes unused :: self

    sedType=outputAnalysisPropertyTypeLinear
    return
  end function sedType

  integer function sedIndexTemplate(self,node,starFormationHistory)
    !% Find the index of the template SEDs to use.
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    use :: Histories           , only : history
    use :: ISO_Varying_String  , only : var_str
    use :: IO_HDF5             , only : hdf5Access        , hdf5Object
    use :: Numerical_Comparison, only : Values_Agree
    use :: File_Utilities      , only : File_Exists       , File_Lock , File_Unlock, lockDescriptor
    use :: String_Handling     , only : operator(//)
    implicit none
    class  (nodePropertyExtractorSED), intent(inout) :: self
    type   (treeNode                ), intent(inout) :: node
    type   (history                 ), intent(in   ) :: starFormationHistory
    class  (nodeComponentBasic      ), pointer       :: basic
    integer(c_size_t                )                :: indexOutput
    type   (lockDescriptor          )                :: fileLock
    type   (hdf5Object              )                :: file
    type   (varying_string          )                :: datasetName

    ! Return a negative index if templates are not being used.
    sedIndexTemplate=-1
    if (.not.self%useSEDTemplates) return
    ! Check that the node exists at at output time.
    basic       => node             %basic(                               )
    indexOutput =  self%outputTimes_%index(basic%time(),findClosest=.true.)
    if (.not.Values_Agree(basic%time(),self%outputTimes_%time(indexOutput),relTol=1.0d-6)) return
    sedIndexTemplate=int(indexOutput)
    ! Ensure that the templates have been built for this index.
    if (.not.allocated(self%templates)) allocate(self%templates(self%outputTimes_%count()))
    if (.not.allocated(self%templates(sedIndexTemplate)%sed)) then
       ! Check if the templates can be retrieved from file.
       datasetName=var_str("sedTemplates")//indexOutput
       !! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
       if (File_Exists(self%fileName)) then
          !$ call hdf5Access%set()
          call file%openFile(char(self%fileName))
          if (file%hasDataset(char(datasetName))) call file%readDataset(char(datasetName),self%templates(sedIndexTemplate)%sed)
          call file%close()
          !$ call hdf5Access%unset()
       end if
       if (.not.allocated(self%templates(sedIndexTemplate)%sed)) then
          self%templates(sedIndexTemplate)%sed=self%luminosityMean(self%outputTimes_%time(indexOutput),starFormationHistory,parallelize=.true.)
          !$ call hdf5Access%set()
          call file%openFile(char(self%fileName),overWrite=.false.,readOnly=.false.)
          call file%writeDataset(self%templates(sedIndexTemplate)%sed,char(datasetName))
          call file%close()
          !$ call hdf5Access%unset()
       end if
       call File_Unlock(fileLock)
    end if
    return
  end function sedIndexTemplate
  
  double precision function sedLuminosityMean(self,time,starFormationHistory,parallelize)
    !% Compute the mean luminosity of the stellar population in each bin of the star formation history.
    use    :: Abundances_Structure , only : abundances           , metallicityTypeLinearByMassSolar, adjustElementsReset
    use    :: Display              , only : displayIndent        , displayUnindent                 , displayCounter     , displayCounterClear, &
         &                                  verbosityLevelWorking
    use    :: Histories            , only : history
    use    :: Numerical_Integration, only : integrator
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num
    implicit none
    double precision                               , dimension(:,:,:), allocatable :: sedLuminosityMean
    class           (nodePropertyExtractorSED     ), intent(inout)                 :: self
    double precision                               , intent(in   )                 :: time
    type            (history                      ), intent(in   )                 :: starFormationHistory
    logical                                        , intent(in   )   , optional    :: parallelize
    class           (stellarPopulationSpectraClass), pointer         , save        :: stellarPopulationSpectra_
    type            (integrator                   ), allocatable     , save        :: integratorTime           , integratorMetallicity
    integer                                                                        :: iWavelength              , iTime                , &
         &                                                                            iMetallicity
    double precision                                                               :: metallicityMinimum       , metallicityMaximum
    double precision                                                 , save        :: timeMinimum              , timeMaximum          , &
         &                                                                            wavelength
    type            (abundances                   )                  , save        :: abundancesStellar
    character       (len=12                       )                                :: label
    integer                                                                        :: counter
    !$omp threadprivate(stellarPopulationSpectra_,integratorTime,integratorMetallicity,abundancesStellar,wavelength,timeMinimum,timeMaximum)
    !# <optionalArgument name="parallelize" defaultsTo=".false." />

    allocate(sedLuminosityMean(self%countWavelengths,size(starFormationHistory%data,dim=1),size(starFormationHistory%data,dim=2)))
    counter=-1
    !$omp parallel private (iWavelength,iTime,iMetallicity,metallicityMinimum,metallicityMaximum)
    allocate(integratorTime       )
    allocate(integratorMetallicity)
    integratorTime       =integrator(sedIntegrandTime       ,toleranceRelative=1.0d-3)
    integratorMetallicity=integrator(sedIntegrandMetallicity,toleranceRelative=1.0d-3)
    if (parallelize_) then
       allocate(stellarPopulationSpectra_,mold=self%stellarPopulationSpectra_)
       !$omp critical(nodePropertyExtractSEDDeepCopy)
       !# <deepCopyReset variables="self%stellarPopulationSpectra_"/>
       !# <deepCopy source="self%stellarPopulationSpectra_" destination="stellarPopulationSpectra_"/>
       !# <deepCopyFinalize variables="stellarPopulationSpectra_"/>
       !$omp end critical(nodePropertyExtractSEDDeepCopy)
    end if
    !$omp master
    if (parallelize_) then
       write (label,'(f12.8)') time
       call displayIndent("computing template SEDs for time "//trim(adjustl(label))//" Gyr",verbosityLevelWorking)
    end if
    !$omp end master
    !$omp do
    do iWavelength=1,self%countWavelengths
       if (parallelize_) then
          !$omp atomic
          counter=counter+1
          call displayCounter(percentageComplete=int(100.0d0*dble(counter)/dble(self%countWavelengths)),isNew=counter==0,verbosity=verbosityLevelWorking)
       end if
       wavelength=self%wavelengths(iWavelength)
       do iTime=1,size(starFormationHistory%data,dim=1)
          if (iTime == 1) then
             timeMinimum=                         0.0d0
          else
             timeMinimum=    starFormationHistory%time(iTime-1)
          end if
          timeMaximum   =min(starFormationHistory%time(iTime  ),time)
          if (timeMaximum <= timeMinimum) cycle
          do iMetallicity=1,size(starFormationHistory%data,dim=2)
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
       end do
    end do
    !$omp end do
    !$omp master
    if (parallelize_) then
       call displayCounterClear(       verbosityLevelWorking)
       call displayUnindent    ("done",verbosityLevelWorking)
    end if
    !$omp end master
    !# <objectDestructor name="stellarPopulationSpectra_"/>
    deallocate(integratorTime       )
    deallocate(integratorMetallicity)
    !$omp end parallel
    return

  contains
    
    double precision function sedIntegrandMetallicity(metallicity)
      !% Integrand over metallicity of the stellar population.
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
      !% Integrand over birth time of the stellar population.
      implicit none
      double precision, intent(in   ) :: timeBirth
      double precision                :: age

      age             =min(                            &
           &               +     time                  &
           &               -     timeBirth           , &
           &               +self%agePopulationMaximum  &
           &              )
      if (parallelize) then
         sedIntegrandTime=     stellarPopulationSpectra_%luminosity(abundancesStellar,age,wavelength)
      else
         sedIntegrandTime=self%stellarPopulationSpectra_%luminosity(abundancesStellar,age,wavelength)
      end if
      return
    end function sedIntegrandTime

  end function sedLuminosityMean
  

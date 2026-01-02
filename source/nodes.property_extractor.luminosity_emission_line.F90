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

  !+    Contributions to this file made by: Sachi Weerasooriya

  !!{
  Implements a property extractor class for the emission line luminosity of a component.
  !!}
  use, intrinsic :: ISO_C_Binding                   , only : c_size_t
  use            :: Galactic_Structure_Options      , only : enumerationComponentTypeType
  use            :: Output_Times                    , only : outputTimesClass
  use            :: Star_Formation_Histories        , only : starFormationHistoryClass
  use            :: HII_Region_Luminosity_Functions , only : hiiRegionLuminosityFunctionClass
  use            :: HII_Region_Mass_Functions       , only : hiiRegionMassFunctionClass
  use            :: Star_Formation_Histories        , only : starFormationHistoryClass
  use            :: HII_Region_Density_Distributions, only : hiiRegionDensityDistributionClass  
  use            :: HII_Region_Escape_Fraction      , only : hiiRegionEscapeFractionClass

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
   <description>
    An emission line luminosity property extractor class. The luminosity of the named emission line (given by the {\normalfont
    \ttfamily lineNames} parameter: if multiple lines are named, the sum of their luminosities) is computed.
   </description>
   <runTimeFileDependencies paths="cloudyTableFileName"/>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorLuminosityEmissionLine
     !!{
     A property extractor class for the emission line luminosity of a component.
     !!}
     private
     class           (starFormationHistoryClass        ), pointer                       :: starFormationHistory_                => null()
     class           (outputTimesClass                 ), pointer                       :: outputTimes_                         => null()
     class           (hiiRegionLuminosityFunctionClass ), pointer                       :: hiiRegionLuminosityFunction_         => null()
     class           (hiiRegionMassFunctionClass       ), pointer                       :: hiiRegionMassFunction_               => null()
     class           (hiiRegionDensityDistributionClass), pointer                       :: hiiRegionDensityDistribution_        => null()
     class           (hiiRegionEscapeFractionClass     ), pointer                       :: hiiRegionEscapeFraction_             => null()
     type            (enumerationComponentTypeType     )                                :: component
     type            (varying_string                   )                                :: parametersGroupPath
     integer                                                                            :: countLines
     integer         (c_size_t                         )                                :: indexAge                                     , indexMetallicity            , &
          &                                                                                indexIonizingLuminosityHydrogen              , indexDensityHydrogen        , &
          &                                                                                indexMassStellar
     type            (varying_string                   ), allocatable, dimension(:    ) :: lineNames                                    , names_                      , &
          &                                                                                descriptions_
     double precision                                   , allocatable, dimension(:    ) :: metallicityBoundaries                        , metallicities               , &
          &                                                                                ages                                         , wavelengths
     double precision                                   , allocatable, dimension(:,:,:) :: luminositiesReduced
     double precision                                   , allocatable, dimension(:,:  ) :: ionizingLuminosityHydrogenNormalized
     type            (emissionLineLuminosityTemplate   ), allocatable, dimension(:    ) :: templates
     double precision                                                                   :: metallicityPopulationMinimum                 , metallicityPopulationMaximum, &
          &                                                                                agePopulationMaximum                         , resolution                  , &
          &                                                                                factorWavelength                             , toleranceRelative           , &
          &                                                                                ionizingLuminosityHydrogenMean               , massStellarMean
     type            (varying_string                   )                                :: cloudyTableFileName
     logical                                                                            :: tabulatedByMass
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
     Constructors for the \refClass{nodePropertyExtractorLuminosityEmissionLine} output analysis class.
     !!}
     module procedure emissionLineLuminosityConstructorParameters
     module procedure emissionLineLuminosityConstructorInternal
  end interface nodePropertyExtractorLuminosityEmissionLine
      
contains

  function emissionLineLuminosityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorLuminosityEmissionLine} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    use :: IO_HDF5                   , only : hdf5Object
    implicit none
    type            (nodePropertyExtractorLuminosityEmissionLine)                              :: self
    type            (inputParameters                            ), intent(inout)               :: parameters
    class           (starFormationHistoryClass                  ), pointer                     :: starFormationHistory_
    class           (outputTimesClass                           ), pointer                     :: outputTimes_
    class           (hiiRegionLuminosityFunctionClass           ), pointer                     :: hiiRegionLuminosityFunction_
    class           (hiiRegionMassFunctionClass                 ), pointer                     :: hiiRegionMassFunction_
    class           (hiiRegionDensityDistributionClass          ), pointer                     :: hiiRegionDensityDistribution_
    class           (hiiRegionEscapeFractionClass               ), pointer                     :: hiiRegionEscapeFraction_
    type            (varying_string                             ), allocatable  , dimension(:) :: lineNames
    type            (varying_string                             )                              :: component                   , cloudyTableFileName
    type            (hdf5Object                                 )                              :: parametersGroup
    double precision                                                                           :: toleranceRelative
    
    allocate(lineNames(parameters%count('lineNames')))
    !![
    <inputParameter>
      <name>cloudyTableFileName</name>
      <defaultValue>var_str('%DATASTATICPATH%/hiiRegions/emissionLineLuminosities_BC2003_highResolution_imfChabrier.hdf5')</defaultValue>
      <source>parameters</source>
      <description>The file of emission line luminosities to use.</description>
    </inputParameter>
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
    <objectBuilder class="starFormationHistory"         name="starFormationHistory_"         source="parameters"/>
    <objectBuilder class="outputTimes"                  name="outputTimes_"                  source="parameters"/>
    <objectBuilder class="hiiRegionLuminosityFunction"  name="hiiRegionLuminosityFunction_"  source="parameters"/>
    <objectBuilder class="hiiRegionMassFunction"        name="hiiRegionMassFunction_"        source="parameters"/>
    <objectBuilder class="hiiRegionDensityDistribution" name="hiiRegionDensityDistribution_" source="parameters"/>
    <objectBuilder class="hiiRegionEscapeFraction"      name="hiiRegionEscapeFraction_"      source="parameters"/>
    !!]
    self=nodePropertyExtractorLuminosityEmissionLine(cloudyTableFileName,enumerationComponentTypeEncode(char(component),includesPrefix=.false.),lineNames,toleranceRelative,starFormationHistory_,outputTimes_,hiiRegionLuminosityFunction_,hiiRegionMassFunction_,hiiRegionDensityDistribution_,hiiRegionEscapeFraction_)
    parametersGroup=parameters%parametersGroup()
    self%parametersGroupPath=parametersGroup%pathTo(includeFileName=.false.)
    call parametersGroup%close()
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationHistory_"        />
    <objectDestructor name="outputTimes_"                 />
    <objectDestructor name="hiiRegionLuminosityFunction_" />
    <objectDestructor name="hiiRegionMassFunction_"       />
    <objectDestructor name="hiiRegionDensityDistribution_"/>
    <objectDestructor name="hiiRegionEscapeFraction_"     />
    !!]
    return
  end function emissionLineLuminosityConstructorParameters

  function emissionLineLuminosityConstructorInternal(cloudyTableFileName,component,lineNames,toleranceRelative,starFormationHistory_,outputTimes_,hiiRegionLuminosityFunction_,hiiRegionMassFunction_,hiiRegionDensityDistribution_,hiiRegionEscapeFraction_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorLuminosityEmissionLine} property extractor class.
    !!}
    use :: Array_Utilities                 , only : slice5Dto2D
    use :: Galactic_Structure_Options      , only : componentTypeDisk, componentTypeSpheroid, componentTypeAll
    use :: Galacticus_Nodes                , only : nodeComponentDisk, nodeComponentSpheroid
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : metallicitySolar
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Error                           , only : Error_Report
    use :: Input_Paths                     , only : inputPath        , pathTypeDataStatic
    implicit none
    type            (nodePropertyExtractorLuminosityEmissionLine)                                      :: self
    type            (varying_string                             ), intent(in   )                       :: cloudyTableFileName
    type            (enumerationComponentTypeType               ), intent(in   )                       :: component
    type            (varying_string                             ), intent(in   ), dimension(:        ) :: lineNames
    class           (starFormationHistoryClass                  ), intent(in   ), target               :: starFormationHistory_
    class           (outputTimesClass                           ), intent(in   ), target               :: outputTimes_
    class           (hiiRegionLuminosityFunctionClass           ), intent(in   ), target               :: hiiRegionLuminosityFunction_
    class           (hiiRegionMassFunctionClass                 ), intent(in   ), target               :: hiiRegionMassFunction_
    class           (hiiRegionDensityDistributionClass          ), intent(in   ), target               :: hiiRegionDensityDistribution_
    class           (hiiRegionEscapeFractionClass               ), intent(in   ), target               :: hiiRegionEscapeFraction_
    double precision                                             , intent(in   )                       :: toleranceRelative    
    double precision                                             ,                                     :: deltaIonizingLuminosityHydrogen   , deltaMassStellar                  , &
         &                                                                                                rateHydrogenIonizingPhotonsMinimum, rateHydrogenIonizingPhotonsMaximum, &
         &                                                                                                massStellarMinimum                , massStellarMaximum                , &
         &                                                                                                densityHydrogenMinimum            , densityHydrogenMaximum            , &
         &                                                                                                deltaDensityHydrogen              , fractionNormalization
    double precision                                             , allocatable  , dimension(:        ) :: ionizingLuminosityHydrogen        , densityHydrogen                   , &
         &                                                                                                massStellar
    integer                                                                     , dimension(5        ) :: shapeLines
    double precision                                             , allocatable  , dimension(:,:,:,:,:) :: luminosities
    integer         (c_size_t                                   )               , dimension(2        ) :: permutation
    type            (hdf5Object                                 )                                      :: emissionLinesFile                 , lines                             , &
         &                                                                                                dataset
    integer         (c_size_t                                   )                                      :: i                                 , k                                 , &
         &                                                                                                iAge                              , indexNormalization                , &
         &                                                                                                sizeNormalization
    !![
    <constructorAssign variables="cloudyTableFileName, lineNames, component, toleranceRelative, *starFormationHistory_, *outputTimes_, *hiiRegionLuminosityFunction_, *hiiRegionMassFunction_, *hiiRegionDensityDistribution_, *hiiRegionEscapeFraction_"/>
    !!]
    if     (                                                                                                    &
         &   component /= componentTypeDisk                                                                     &
         &  .and.                                                                                               &
         &   component /= componentTypeSpheroid                                                                 &
         &  .and.                                                                                               &
         &   component /= componentTypeAll                                                                      &
         & ) call Error_Report("only 'disk' and 'spheroid' components are supported"//{introspection:location})
    ! Initialize no output Parameters group.
    self%parametersGroupPath=""
    ! Get details of the star formation rate tabulation.
    self%metallicityBoundaries=self%starFormationHistory_%metallicityBoundaries      ()
    !$ call hdf5Access%set()
    call emissionLinesFile%openFile(self%cloudyTableFileName,readOnly=.true.)
    lines=emissionLinesFile%openGroup('lines')
    do i=1,size(lineNames)
       if (.not.lines%hasDataset(char(self%lineNames(i)))) call Error_Report('line "'//char(self%lineNames(i))//'" not found'//{introspection:location})
    end do
    ! Detect if the tabulation is by Qₕ or M⭑.
    if      (emissionLinesFile%hasDataset('ionizingLuminosityHydrogen')) then
       self%tabulatedByMass=.false.
    else if (emissionLinesFile%hasDataset('massStellar'               )) then
       self%tabulatedByMass=.true.
    else
       self%tabulatedByMass=.false.
       call Error_Report('expected to find either `ionizingLuminosityHydrogen` or `massStellar` in the tabulation'//{introspection:location})
    end if
    call    emissionLinesFile%readDataset('ionizingLuminosityHydrogenNormalized',self%ionizingLuminosityHydrogenNormalized)
    call    emissionLinesFile%readDataset('metallicity'                         ,self%metallicities                       )
    call    emissionLinesFile%readDataset('age'                                 ,self%ages                                )
    call    emissionLinesFile%readDataset('densityHydrogen'                     ,     densityHydrogen                     )
    if (self%tabulatedByMass) then
       call emissionLinesFile%readDataset('massStellar'                         ,     massStellar                         )
    else
       call emissionLinesFile%readDataset('ionizingLuminosityHydrogen'          ,     ionizingLuminosityHydrogen          )
    end if
    ! Extract indexing into the lines arrays.
    dataset   =emissionLinesFile%openDataset('metallicity'               )
    call    dataset%readAttribute('index',self%indexMetallicity               )
    call    dataset%close        (                     )
    dataset   =emissionLinesFile%openDataset('age'                       )
    call    dataset%readAttribute('index',self%indexAge                       )
    call    dataset%close        (                     )
    dataset=   emissionLinesFile%openDataset('densityHydrogen'           )
    call    dataset%readAttribute('index',self%indexDensityHydrogen           )
    call    dataset%close        (                     )
    if (self%tabulatedByMass) then
       dataset=emissionLinesFile%openDataset('massStellar'               )
       call dataset%readAttribute('index',self%indexMassStellar               )
       call dataset%close        (                     )
       self%indexIonizingLuminosityHydrogen=-huge(0_c_size_t)
    else
       dataset=emissionLinesFile%openDataset('ionizingLuminosityHydrogen')
       call dataset%readAttribute('index',self%indexIonizingLuminosityHydrogen)
       call dataset%close        (                     )
       self%indexMassStellar               =-huge(0_c_size_t)
    end if
    ! Offset indexing to Fortran standard (i.e. starting from 1 instead of 0).
    self%indexMetallicity               =self%indexMetallicity               +1
    self%indexAge                       =self%indexAge                       +1
    self%indexIonizingLuminosityHydrogen=self%indexIonizingLuminosityHydrogen+1
    self%indexMassStellar               =self%indexMassStellar               +1
    self%indexDensityHydrogen           =self%indexDensityHydrogen           +1
    ! Establish arrays.
    self%metallicityPopulationMinimum=minval(self%metallicities)
    self%metallicityPopulationMaximum=maxval(self%metallicities)
    self%agePopulationMaximum        =maxval(self%ages         )
    shapeLines   (self%indexMetallicity               )=size(self%metallicities             )
    shapeLines   (self%indexAge                       )=size(self%ages                      )
    if (self%tabulatedByMass) then
       shapeLines(self%indexMassStellar               )=size(     massStellar               )
    else
       shapeLines(self%indexIonizingLuminosityHydrogen)=size(     ionizingLuminosityHydrogen)
    end if
    shapeLines   (self%indexDensityHydrogen           )=size(     densityHydrogen           )
    shapeLines   (     5                              )=size(     lineNames                 )
    !![
    <allocate variable="luminosities" shape="shapeLines"/>
    !!]
    allocate(                                        &
         &   self%luminositiesReduced                &
         &   (                                       &
         &    size(self%ages                      ), &
         &    size(self%metallicities             ), &
         &    size(     lineNames                 )  &
         &   )                                       &
         &  )
    allocate(                                        &
         &   self%wavelengths                        &
         &   (                                       &
         &    size(     lineNames                 )  &
         &   )                                       &
         &  )
    self%luminositiesReduced=0.0d0
    do i=1,size(lineNames)
       dataset=lines%openDataset(char(lineNames(i)))
       call dataset%readAttribute    ('wavelength'      ,self%wavelengths (        i))
       call dataset%close            (                                               )
       call lines  %readDatasetStatic(char(lineNames(i)),     luminosities(:,:,:,:,i))
    end do
    call lines            %close()
    call emissionLinesFile%close()
    !$ call hdf5Access%unset()
    ! Calculate emission line luminosities as a function of age and metallicity by averaging over the distribution of HII region
    ! luminosities/masses. Account for any needed permutation in the indexing needed to get our final array to be
    ! (age,metallicity,line) ordered.    
    !! Find the width of Qₕ/M⭑ and nₕ bins.
    if (self%tabulatedByMass) then
       deltaMassStellar               =+massStellar               (2) &
            &                          /massStellar               (1)
       deltaIonizingLuminosityHydrogen=-huge(0.0d0)
    else
       deltaIonizingLuminosityHydrogen=+ionizingLuminosityHydrogen(2) &
            &                          /ionizingLuminosityHydrogen(1)
       deltaMassStellar               =-huge(0.0d0)
    end if
    deltaDensityHydrogen              =+densityHydrogen           (2) &
         &                             /densityHydrogen           (1)
    !! Construct the permutation needed to convert the table read from file into the ordering needed internally.
    permutation(1)=self%indexMetallicity
    permutation(2)=     2
    if (self%tabulatedByMass) then
       if (self%indexMassStellar                < self%indexMetallicity) permutation(1)=permutation(1)-1_c_size_t
       indexNormalization=self%indexMassStellar
       sizeNormalization =size(massStellar               )
    else
       if (self%indexIonizingLuminosityHydrogen < self%indexMetallicity) permutation(1)=permutation(1)-1_c_size_t
       indexNormalization=self%indexIonizingLuminosityHydrogen
       sizeNormalization =size(ionizingLuminosityHydrogen)
    end if
    if    (self%indexDensityHydrogen            < self%indexMetallicity) permutation(1)=permutation(1)-1_c_size_t
    if    (self%indexAge                        < self%indexMetallicity) permutation(1)=permutation(1)-1_c_size_t
    do i=1,sizeNormalization
       if (self%tabulatedByMass) then
          massStellarMinimum                =massStellar               (i)/sqrt(deltaMassStellar               )
          massStellarMaximum                =massStellar               (i)*sqrt(deltaMassStellar               )
       else
          rateHydrogenIonizingPhotonsMinimum=ionizingLuminosityHydrogen(i)/sqrt(deltaIonizingLuminosityHydrogen)
          rateHydrogenIonizingPhotonsMaximum=ionizingLuminosityHydrogen(i)*sqrt(deltaIonizingLuminosityHydrogen)
       end if
       do k=1,size(densityHydrogen)
          densityHydrogenMinimum         =densityHydrogen           (k)/sqrt(deltaDensityHydrogen           )
          densityHydrogenMaximum         =densityHydrogen           (k)*sqrt(deltaDensityHydrogen           )
          do iAge=1,size(self%ages)
             ! Accumulate the luminosity weighted by the cumulative fraction of HII regions in this luminosity interval.
             if (self%tabulatedByMass) then
                fractionNormalization=self%hiiRegionMassFunction_      %cumulativeDistributionFunction(massStellarMinimum                ,massStellarMaximum                )
             else
                fractionNormalization=self%hiiRegionLuminosityFunction_%cumulativeDistributionFunction(rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum)
             end if
             self%luminositiesReduced(iAge,:,:)=+        self%luminositiesReduced                                        (          iAge        ,:,:                   ) &
                  &                             +(+1.0d0-self%hiiRegionEscapeFraction_     %escapeFraction               (self%ages(iAge))                             ) & 
                  &                             *                                           fractionNormalization                                                        &
                  &                             *        self%hiiRegionDensityDistribution_%cumulativeDensityDistribution(densityHydrogenMinimum,densityHydrogenMaximum) &
                  &                             *reshape(                                                                                                                &
                  &                                      slice5Dto2D(                                                                                                    &
                  &                                                   luminosities                                                           ,                           &
                  &                                                   [indexNormalization,self%indexDensityHydrogen,self%indexAge],                                      &
                  &                                                   [i                 ,     k                   ,     iAge    ]                                       &
                  &                                                  )                                                                       ,                           &
                  &                                                   [size(luminosities,dim=self%indexMetallicity),size(luminosities,dim=5)],                           &
                  &                                       order=permutation                                                                                              &
                  &                                      )
          end do
       end do
    end do
    ! Normalize reduced luminosities to the total fraction of HII regions in the luminosity/mass interval spanned by the table. Also,
    ! find the mean ionizing luminosity of HII regions in this luminosity/mass interval.
    if (self%tabulatedByMass) then
       massStellarMinimum                 =+massStellar               (                              1 )/sqrt(deltaMassStellar               )
       massStellarMaximum                 =+massStellar               (size(massStellar               ))*sqrt(deltaMassStellar               )
       fractionNormalization              =+self%hiiRegionMassFunction_      %cumulativeDistributionFunction(                massStellarMinimum,                massStellarMaximum)
    else
       rateHydrogenIonizingPhotonsMinimum =+ionizingLuminosityHydrogen(                              1 )/sqrt(deltaIonizingLuminosityHydrogen)
       rateHydrogenIonizingPhotonsMaximum =+ionizingLuminosityHydrogen(size(ionizingLuminosityHydrogen))*sqrt(deltaIonizingLuminosityHydrogen)
       fractionNormalization              =+self%hiiRegionLuminosityFunction_%cumulativeDistributionFunction(rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum)
    end if
    densityHydrogenMinimum                =+densityHydrogen           (                              1 )/sqrt(deltaDensityHydrogen           )
    densityHydrogenMaximum                =+densityHydrogen           (size(           densityHydrogen))*sqrt(deltaDensityHydrogen           )
    self%luminositiesReduced              =+self%                              luminositiesReduced                                                                                   &
         &                                 /                                   fractionNormalization                                                                                 &
         &                                 /self%hiiRegionDensityDistribution_%cumulativeDensityDistribution (            densityHydrogenMinimum,            densityHydrogenMaximum)
    if (self%tabulatedByMass) then
       self%massStellarMean               =+self%hiiRegionMassFunction_       %cumulativeMass                (                massStellarMinimum,                massStellarMaximum)
    else
       self%ionizingLuminosityHydrogenMean=+self%hiiRegionLuminosityFunction_ %cumulativeLuminosity          (rateHydrogenIonizingPhotonsMinimum,rateHydrogenIonizingPhotonsMaximum)
    end if
    ! Construct property names and descriptions.
    allocate(self%names_       (size(lineNames)))
    allocate(self%descriptions_(size(lineNames)))
    do i=1,size(lineNames)
       select case (self%component%ID)
       case (componentTypeDisk    %ID)
          self%names_       (i)="luminosityEmissionLineDisk:"    //lineNames(i)
       case (componentTypeSpheroid%ID)
          self%names_       (i)="luminosityEmissionLineSpheroid:"//lineNames(i)
       case (componentTypeAll%ID)
          self%names_       (i)="luminosityEmissionLineTotal:"   //lineNames(i)
       end select
       self   %descriptions_(i)="Luminosity of the "             //lineNames(i)//" emission line [ergs/s]"
    end do
    self%countLines=size(lineNames)
    return    
  end function emissionLineLuminosityConstructorInternal

  subroutine emissionLineLuminosityDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorLuminosityEmissionLine} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationHistory_"        />
    <objectDestructor name="self%outputTimes_"                 />
    <objectDestructor name="self%hiiRegionLuminosityFunction_" />
    <objectDestructor name="self%hiiRegionMassFunction_"       />
    <objectDestructor name="self%hiiRegionDensityDistribution_"/>
    <objectDestructor name="self%hiiRegionEscapeFraction_"     />
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
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, componentTypeAll
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
    double precision                                             , dimension(  :,:)          , allocatable :: masses
    type            (history                                    )                                          :: starFormationHistory
    integer(c_size_t)                                                                                      :: indexTemplate      
    logical                                                                                                :: parallelize
    integer                                                                                                :: iLines              , i
    !$GLC attributes unused :: instance

    allocate(luminosity(self%countLines))
    luminosity=0.0d0
    do i=1,2
       ! Get the relevant star formation history.
       select case (i)
       case (1)
          if (self%component == componentTypeSpheroid) cycle
          disk                 => node    %disk                ()
          starFormationHistory =  disk    %starFormationHistory()
       case (2)
          if (self%component == componentTypeDisk) cycle
          spheroid             => node    %spheroid            ()
          starFormationHistory =  spheroid%starFormationHistory()
       end select
       if (.not.starFormationHistory%exists()) cycle
       ! Get the index of the template to use.
       indexTemplate=self%indexTemplateNode(node,starFormationHistory)
       if (indexTemplate > 0) then
          ! Stored templates can be used, so point to the relevant set.
          luminosityTemplate_ => self%templates(indexTemplate)%emissionLineLuminosity
       else
          ! Stored templates can not be used, get the templates for this specific case, and point to them.
          luminosityTemplate  =  self%luminosityMean(time,node,indexTemplate,starFormationHistory,parallelize)
          luminosityTemplate_ => luminosityTemplate
       end if
       masses=self%starFormationHistory_%masses(node,starFormationHistory,allowTruncation=.false.)
       do iLines=1,size(luminosity,dim=1)
          luminosity(iLines)=luminosity(iLines)+sum(luminosityTemplate_(iLines,:,:)*masses(:,:))
       end do
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

  integer function emissionLineLuminosityIndexTemplateTime(self,time) result(indexTemplate)
    !!{
    Find the index of the template emission lines to use.
    !!}
    use :: Numerical_Comparison    , only : Values_Agree
    use :: Star_Formation_Histories, only : starFormationHistoryAgesFixed, starFormationHistoryAgesFixedPerOutput
    implicit none
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout) :: self
    double precision                                             , intent(in   ) :: time
    integer         (c_size_t                                   )                :: indexOutput

    if      (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed         ) then
       ! Ages are fixed - a single template can be used.
       indexTemplate =+1
    else if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixedPerOutput) then
       ! Check that the time is an output time.
       indexOutput=self%outputTimes_%index(time,findClosest=.true.)
       if (Values_Agree(time,self%outputTimes_%time(indexOutput),relTol=1.0d-6)) then
          ! The time corresponds to an output time - use a template.
          indexTemplate=int(indexOutput)
       else
          ! The time does not correspond to an output time - a template can not be used.
          indexTemplate=-1
       end if
    else
       ! Templates can not be used.
       indexTemplate =-1
    end if
    return
  end function emissionLineLuminosityIndexTemplateTime

  function emissionLineLuminosityIndexTemplateNode(self,node,starFormationHistory) result(indexTemplate)
    !!{
    Find the index of the template emission line luminosity to use, and also compute the template.
    !!}
    use :: Display                 , only : displayMessage               , verbosityLevelWorking
    use :: Galacticus_Nodes        , only : nodeComponentBasic
    use :: Histories               , only : history
    use :: ISO_Varying_String      , only : var_str
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: Output_HDF5             , only : outputFile
    use :: Numerical_Comparison    , only : Values_Agree
    use :: File_Utilities          , only : File_Exists                  , File_Lock                             , File_Unlock, lockDescriptor, &
         &                                  Directory_Make               , File_Path
    use :: String_Handling         , only : operator(//)
    use :: Input_Paths             , only : inputPath                    , pathTypeDataDynamic
    use :: Star_Formation_Histories, only : starFormationHistoryAgesFixed, starFormationHistoryAgesFixedPerOutput
    implicit none
    integer         (c_size_t                                   )                              :: indexTemplate
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)               :: self
    type            (treeNode                                   ), intent(inout)               :: node
    type            (history                                    ), intent(in   )               :: starFormationHistory
    class           (nodeComponentBasic                         ), pointer                     :: basic
    double precision                                             , allocatable  , dimension(:) :: times
    type            (hdf5Object                                 ), allocatable  , dimension(:) :: parametersGroups
    integer         (c_size_t                                   )                              :: indexOutput         , countTemplates, &
         &                                                                                        i
    type            (lockDescriptor                             )                              :: fileLock
    type            (hdf5Object                                 )                              :: file
    type            (varying_string                             )                              :: fileName
    character       (len=16                                     )                              :: label

    if      (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed         ) then
       ! Ages are fixed - a single template can be used.
       indexTemplate =1
       countTemplates=1
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
    if (.not.allocated(self%templates                                      )) allocate(self%templates(countTemplates))
    if (.not.allocated(self%templates(indexTemplate)%emissionLineLuminosity)) then
       ! Construct the file name.
       fileName=inputPath(pathTypeDataDynamic)                                          // &
            &        'stellarPopulations/'                                              // &
            &        self%objectType             (                                     )// &
            &        '_'                                                                // &
            &        self%historyHashedDescriptor(node,indexOutput,starFormationHistory)// &
            &        '_'                                                                // &
            &        indexTemplate                                                      // &
            &        '.hdf5'
       ! Store the file name used to the output file parameters group for this object.
       !$ call hdf5Access%set()
       if (self%parametersGroupPath /= "") then
          parametersGroups=outputFile%openGroupPath(char(self%parametersGroupPath))
          call parametersGroups(size(parametersGroups))%writeAttribute(fileName,char(var_str('meta:luminosityEmissionLineMatrixFileName')//indexTemplate))
          do i=1,size(parametersGroups)
             call parametersGroups(i)%close()
          end do
       end if
       !$ call hdf5Access%unset()
       ! Check if the templates can be retrieved from file.
       !! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call Directory_Make(File_Path(fileName))
       call File_Lock(fileName,fileLock,lockIsShared=.false.)
       if (File_Exists(fileName)) then
          !$ call hdf5Access%set()
          call file%openFile(char(fileName))
          if (file%hasDataset('luminosityTemplate')) then
             if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed) then
                call displayMessage("reading emission line luminosity tabulation from file '"                                        //fileName//"'",verbosityLevelWorking)
             else
                !$omp critical(gfortranInternalIO)
                write (label,'(f12.8)') self%outputTimes_%time(indexOutput)
                !$omp end critical(gfortranInternalIO)
                call displayMessage("reading emission line luminosity tabulation for time "//trim(adjustl(label))//" Gyr from file '"//fileName//"'",verbosityLevelWorking)
             end if
          call file%readDataset('luminosityTemplate',self%templates(indexTemplate)%emissionLineLuminosity)
          end if
          call file%close()
          !$ call hdf5Access%unset()
       end if
       if (.not.allocated(self%templates(indexTemplate)%emissionLineLuminosity)) then
          basic                                                => node%basic         (                                                                                    )
          self%templates(indexTemplate)%emissionLineLuminosity =  self%luminosityMean(basic%time(),node,indexTemplate,starFormationHistory,parallelize=.true.,times_=times)
          if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed) then
             call displayMessage("storing emission line luminosity tabulation to file '"                                        //fileName//"'",verbosityLevelWorking)
          else
             !$omp critical(gfortranInternalIO)
             write (label,'(f12.8)') self%outputTimes_%time(indexOutput)
             !$omp end critical(gfortranInternalIO)
             call displayMessage("storing emission line luminosity tabulation for time "//trim(adjustl(label))//" Gyr to file '"//fileName//"'",verbosityLevelWorking)
          end if
          !$ call hdf5Access%set()
          call    file%openFile(char(fileName),overWrite=.false.,readOnly=.false.)
          call    file%writeDataset(self %templates             (indexTemplate)%emissionLineLuminosity      ,'luminosityTemplate','A matrix mapping star formation history to emission line luminosities.' )
          call    file%writeDataset(self %lineNames                                                         ,'lineNames'         ,'The names of the emission lines'                                        )
          call    file%writeDataset(self %wavelengths                                                       ,'wavelengths'       ,'The wavelengths of the emission lines [Å]'                              )
          call    file%writeDataset(self %metallicityBoundaries                                             ,'metallicity'       ,'The metallicities at which the star formation history is tabulated [Z☉]')
          if (self%starFormationHistory_%ageDistribution() == starFormationHistoryAgesFixed) then
             call file%writeDataset(basic%time                  (             )                       -times,'ages'              ,'The ages at which the star formation history is tabulated [Gyr]'        )
          else
             call file%writeDataset(      times                                                             ,'time'              ,'The times at which the star formation history is tabulated [Gyr]'       )
          end if
          call file%close()
          !$ call hdf5Access%unset()
       end if
       call File_Unlock(fileLock)
    end if
    return
  end function emissionLineLuminosityIndexTemplateNode

  function emissionLineLuminosityMean(self,time,node,indexOutput,starFormationHistory,parallelize,times_) result(luminosityMean)
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
    double precision                                             , dimension(:,:,:)                            , allocatable :: luminosityMean
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(inout)                                             :: self
    double precision                                             , intent(in   )                                             :: time
    type            (treeNode                                   ), intent(inout)                                             :: node
    integer         (c_size_t                                   ), intent(in   )                                             :: indexOutput
    type            (history                                    ), intent(in   )                                             :: starFormationHistory
    logical                                                      , intent(in   )   , optional                                :: parallelize
    double precision                                             , intent(  out)   , optional, dimension(:    ), allocatable :: times_
    double precision                                             , dimension(:    )                            , allocatable :: times
    double precision                                             , dimension(  :,:)                            , allocatable :: masses
    type            (integrator                                 ), save                                        , allocatable :: integratorTime            , integratorMetallicity
    type            (interpolator                               ), save                                        , allocatable :: interpolatorTime          , interpolatorMetallicity
    integer         (c_size_t                                   )                                                            :: iTime                     , iMetallicity           , &
         &                                                                                                                      counter                   , counterMaximum         , &
         &                                                                                                                      iterator
    integer         (c_size_t                                   ), save                                                      :: iLine
    double precision                                                                                                         :: metallicityMinimum        , metallicityMaximum     , &
         &                                                                                                                      timeStart
    double precision                                             , save                                                      :: timeMinimum               , timeMaximum            , &
         &                                                                                                                      age                       , metallicity_
    character       (len=12                                     )                                                            :: label
    type            (multiCounter                               )                                                            :: state
    type            (ompLock                                    )                                                            :: stateLock
    !$omp threadprivate(iLine,integratorTime,integratorMetallicity,interpolatorTime,interpolatorMetallicity,timeMinimum,timeMaximum,age,metallicity_)
    !$GLC attributes initialized :: masses
    !![
    <optionalArgument name="parallelize" defaultsTo=".false." />
    !!]
    
    times =self%starFormationHistory_%times (node=node,indexOutput=indexOutput,starFormationHistory=starFormationHistory,allowTruncation=.false.,timeStart=timeStart)
    masses=self%starFormationHistory_%masses(node=node                        ,starFormationHistory=starFormationHistory,allowTruncation=.false.                    )
    if (present(times_)) times_=times
    allocate(luminosityMean(self%countLines,size(masses,dim=1),size(masses,dim=2)))
    counter       =-1
    counterMaximum=product     ([size(luminosityMean,dim=1              ),size(masses,dim=1              ),size(masses,dim=2              )])
    state         =multiCounter([size(luminosityMean,dim=1,kind=c_size_t),size(masses,dim=1,kind=c_size_t),size(masses,dim=2,kind=c_size_t)])
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
          timeMinimum=    timeStart
       else
          timeMinimum=    times(iTime-1)
       end if
       timeMaximum   =min(times(iTime  ),time)
       if (timeMaximum <= timeMinimum) then
          luminosityMean(iLine,iTime,iMetallicity)=0.0d0
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
      integer                                    :: iTime                   , iMetallicity
      integer         (c_size_t), dimension(0:1) :: interpolateIndexTime    , interpolateIndexMetallicity
      double precision          , dimension(0:1) :: interpolateFactorTime   , interpolateFactorMetallicity
      double precision                           :: massStellarNormalization
      
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
            if (self%tabulatedByMass) then
               massStellarNormalization=+1.0d0                                                                                                            &
                    &                   /self%massStellarMean
            else
               massStellarNormalization=+self%ionizingLuminosityHydrogenNormalized(interpolateIndexTime(iTime),interpolateIndexMetallicity(iMetallicity)) &
                    &                   /self%ionizingLuminosityHydrogenMean
            end if
            integrand=+integrand                                                                                                     &
                 &   +self%luminositiesReduced         (interpolateIndexTime(iTime),interpolateIndexMetallicity(iMetallicity),iLine) &
                 &   *     massStellarNormalization                                                                                  &
                 &   *     interpolateFactorTime                            (iTime)                                                  &
                 &   *     interpolateFactorMetallicity                                                        (iMetallicity)
         end do
      end do
      return
    end function emissionLineLuminosityIntegrandTime

  end function emissionLineLuminosityMean

  function emissionLineLuminosityHistoryHashedDescriptor(self,node,indexOutput,starFormationHistory) result(hashedDescriptor)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters        , only : inputParameters
    use :: String_Handling         , only : String_C_To_Fortran
    use :: Hashes_Cryptographic    , only : Hash_MD5
    use :: FoX_DOM                 , only : setLiveNodeLists
    use :: Histories               , only : history
    use :: File_Utilities          , only : File_Modification_Time
    use :: String_Handling         , only : String_Join
    use :: Star_Formation_Histories, only : starFormationHistoryAgesFixed
    implicit none
    type            (varying_string                             )                              :: hashedDescriptor
    class           (nodePropertyExtractorLuminosityEmissionLine), intent(in   )               :: self
    type            (treeNode                                   ), intent(inout)               :: node
    integer         (c_size_t                                   ), intent(in   )               :: indexOutput
    type            (history                                    ), intent(in   )               :: starFormationHistory
    double precision                                             , allocatable  , dimension(:) :: times
    character       (len=18                                     )                              :: parameterLabel
    character       (len=30                                     )                              :: timeModification
    type            (inputParameters                            )                              :: descriptor
    type            (varying_string                             )                              :: descriptorString    , values
    integer                                                                                    :: i                   , status
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
    ! Add composited object descriptors.
    call self%starFormationHistory_        %descriptor(descriptor)
    call self%outputTimes_                 %descriptor(descriptor)
    call self%hiiRegionLuminosityFunction_ %descriptor(descriptor)
    call self%hiiRegionDensityDistribution_%descriptor(descriptor)
    call self%hiiRegionEscapeFraction_     %descriptor(descriptor)
    ! Add descriptors for star formation history tabulation.
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
       times =self %starFormationHistory_%times(node=node,indexOutput=indexOutput,starFormationHistory=starFormationHistory)
       do i=1,size(times)
          !$omp critical(gfortranInternalIO)
          write (parameterLabel,'(e17.10)') times(i)
          !$omp end critical(gfortranInternalIO)
          values=values//trim(adjustl(parameterLabel))
          if (i < size(times)) values=values//":"
       end do
       call descriptor%addParameter('time'       ,char(values))
    end if
    ! Add line names.
    call descriptor%addParameter('lineNames',char(String_Join(self%lineNames," ")))
    ! Add tolerance.
    !$omp critical(gfortranInternalIO)
    write (parameterLabel,'(e17.10)') self%toleranceRelative
    !$omp end critical(gfortranInternalIO)
    call descriptor%addParameter('toleranceRelative',parameterLabel)
    ! Add descriptor for Cloudy table.
    timeModification=File_Modification_Time(self%cloudyTableFileName,status)
    if (status == errorStatusSuccess) then
       call descriptor%addParameter('runTimeFileDependency1',char(self%cloudyTableFileName//": "//trim(timeModification)))
    else if (status /= errorStatusNotExist) then
       call Error_Report('unable to get file modification time'//{introspection:location})
    end if
    ! Construct the hash.
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



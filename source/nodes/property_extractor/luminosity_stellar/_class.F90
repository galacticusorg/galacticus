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
Implements a stellar mass output analysis property extractor class.
!!}

  use :: Galactic_Structure_Options, only : enumerationComponentTypeType
  use :: ISO_Varying_String        , only : varying_string              , var_str
  use :: Output_Times              , only : outputTimesClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorLuminosityStellar" docformat="rst">
   <description>
   A property extractor that returns the stellar luminosity of a node in a specified broadband filter, in units of the AB zero-point. The ``filterName`` and ``filterType`` parameters select the photometric band and whether to use rest-frame or observer-frame luminosities. The ``component`` parameter selects the galactic component whose luminosity is returned (``all``, the default, sums over all components; ``disk``, ``spheroid``, and ``nuclearStarCluster`` select an individual component---the latter are needed when applying dust attenuation, which cannot meaningfully be applied to a summed luminosity). The optional ``redshiftBand`` shifts the band to a fixed redshift (for K-corrections), and ``postprocessChain`` applies a named spectral postprocessing chain (e.g.\ :term:`IGM` attenuation) before the photometric integration. Luminosity indices are pre-computed per output time for efficiency.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorLuminosityStellar
     !!{RST
     A stellar luminosity output analysis property extractor class.
     !!}
     private
     type            (varying_string              )                              :: filterName                   , filterType, &
          &                                                                         postprocessChain             , name_     , &
          &                                                                         description_
     type            (enumerationComponentTypeType)                              :: component
     double precision                                                            :: redshiftBand
     integer                                       , allocatable, dimension(:  ) :: luminosityIndex
     ! Names of the postprocessing chains used to resolve this luminosity into bins of stellar population age, and the
     ! index of the corresponding luminosity at each output time. Empty unless age resolution has been requested.
     type            (varying_string              ), allocatable, dimension(:  ) :: postprocessChains
     integer                                       , allocatable, dimension(:,:) :: luminosityIndexChain
     class           (outputTimesClass            ), pointer                     :: outputTimes_         => null()
   contains
     final     ::                luminosityStellarDestructor
     procedure :: extract             => luminosityStellarExtract
     procedure :: supportsAttenuation => luminosityStellarSupportsAttenuation
     procedure :: decompose           => luminosityStellarDecompose
     !![
     <methods docformat="rst">
       <method method="luminosityValue" description="Return the luminosity of this extractor's component for a given entry in the stellar luminosities structure."/>
     </methods>
     !!]
     procedure :: luminosityValue => luminosityStellarLuminosityValue
     procedure :: quantity        => luminosityStellarQuantity
     procedure :: name            => luminosityStellarName
     procedure :: description     => luminosityStellarDescription
     procedure :: unitsInSI       => luminosityStellarUnitsInSI
     procedure :: units           => luminosityStellarUnits
  end type nodePropertyExtractorLuminosityStellar

  interface nodePropertyExtractorLuminosityStellar
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorLuminosityStellar` property extractor class.
     !!}
     module procedure luminosityStellarConstructorParameters
     module procedure luminosityStellarConstructorInternal
  end interface nodePropertyExtractorLuminosityStellar

contains

  function luminosityStellarConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellar` property extractor class which takes a parameter set as input.
    !!}
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    implicit none
    type            (nodePropertyExtractorLuminosityStellar)                              :: self
    type            (inputParameters                       ), intent(inout)               :: parameters
    class           (outputTimesClass                      ), pointer                     :: outputTimes_
    type            (varying_string                        )                              :: filterName           , filterType               , &
         &                                                                                   postprocessChain     , component
    type            (varying_string                        ), allocatable  , dimension(:) :: postprocessChains
    double precision                                                                      :: redshiftBand
    logical                                                                               :: redshiftBandIsPresent, postprocessChainIsPresent
    integer                                                                               :: countChains

    redshiftBandIsPresent    =parameters%isPresent('redshiftBand'    )
    postprocessChainIsPresent=parameters%isPresent('postprocessChain')
    !![
    <inputParameter docformat="rst">
      <name>filterName</name>
      <source>parameters</source>
      <description>
      The filter to select.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>filterType</name>
      <source>parameters</source>
      <description>
      The filter type (rest or observed) to select.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>component</name>
      <defaultValue>var_str('all')</defaultValue>
      <source>parameters</source>
      <description>
      The galactic component whose luminosity is to be extracted---one of ``all``, ``disk``, ``spheroid``, or
      ``nuclearStarCluster``. Note that dust attenuation can not be applied to the ``all`` case, since the components
      are attenuated differently and so must be attenuated separately before being summed.
      </description>
    </inputParameter>
    !!]
    if (redshiftBandIsPresent) then
       !![
       <inputParameter docformat="rst">
         <name>redshiftBand</name>
         <source>parameters</source>
         <description>
         The redshift of the band (if not the output redshift).
         </description>
       </inputParameter>
       !!]
    end if
    if (postprocessChainIsPresent) then
       !![
       <inputParameter docformat="rst">
         <name>postprocessChain</name>
         <source>parameters</source>
         <description>
         The postprocessing chain to use.
         </description>
       </inputParameter>
       !!]
    end if
    ! The age-resolving chains are read only if given: a deferred-shape parameter can not be read into an unallocated
    ! array, so its extent is established first.
    countChains=parameters%count('postprocessChains',zeroIfNotPresent=.true.)
    allocate(postprocessChains(countChains))
    if (countChains > 0) then
       !![
       <inputParameter docformat="rst">
         <name>postprocessChains</name>
         <source>parameters</source>
         <description>
         The names of the postprocessing chains to use when resolving this luminosity into bins of stellar population
         age, as required when applying a dust model which attenuates young and old populations differently. Each
         named chain must already be defined by the ``[stellarPopulationSpectraPostprocessorBuilder]``, and a
         luminosity using it must be tracked for this filter.

         The age range spanned by each chain is obtained from the chain itself, not declared here, so it can not drift
         out of step with the chain's configuration. Leave this unset---the default---for a dust model which does not
         distinguish populations by age.
         </description>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    if (redshiftBandIsPresent) then
       if (postprocessChainIsPresent) then
          self=nodePropertyExtractorLuminosityStellar(char(filterName),char(filterType),enumerationComponentTypeEncode(char(component),includesPrefix=.false.),outputTimes_,postprocessChains=postprocessChains,redshiftBand=redshiftBand,postprocessChain=char(postprocessChain))
       else
          self=nodePropertyExtractorLuminosityStellar(char(filterName),char(filterType),enumerationComponentTypeEncode(char(component),includesPrefix=.false.),outputTimes_,postprocessChains=postprocessChains,redshiftBand=redshiftBand                                        )
       end if
    else
       if (postprocessChainIsPresent) then
          self=nodePropertyExtractorLuminosityStellar(char(filterName),char(filterType),enumerationComponentTypeEncode(char(component),includesPrefix=.false.),outputTimes_,postprocessChains=postprocessChains,                          postprocessChain=char(postprocessChain))
       else
          self=nodePropertyExtractorLuminosityStellar(char(filterName),char(filterType),enumerationComponentTypeEncode(char(component),includesPrefix=.false.),outputTimes_,postprocessChains=postprocessChains                                                                  )
       end if
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function luminosityStellarConstructorParameters

  function luminosityStellarConstructorInternal(filterName,filterType,component,outputTimes_,redshiftBand,postprocessChain,postprocessChains,outputMask) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellar` property extractor class.
    !!}
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Error                         , only : Error_Report
    use            :: Galactic_Structure_Options    , only : componentTypeAll               , componentTypeDisk             , componentTypeSpheroid, &
         &                                                   componentTypeNuclearStarCluster, enumerationComponentTypeDecode
    use            :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    type            (nodePropertyExtractorLuminosityStellar)                                        :: self
    character       (len=*                                 ), intent(in   )                         :: filterName      , filterType
    type            (enumerationComponentTypeType          ), intent(in   )                         :: component
    class           (outputTimesClass                      ), intent(in   ), target                 :: outputTimes_
    character       (len=*                                 ), intent(in   ), optional               :: postprocessChain
    double precision                                        , intent(in   ), optional               :: redshiftBand
    logical                                                 , intent(in   ), dimension(:), optional :: outputMask
    type            (varying_string                        ), intent(in   ), dimension(:), optional :: postprocessChains
    integer         (c_size_t                              )                                        :: i
    integer                                                                                         :: j                , countChains
    character       (len=7                                 )                                        :: label
    !![
    <constructorAssign variables="filterName, filterType, component, redshiftBand, postprocessChain, *outputTimes_"/>
    !!]

    ! Validate the component. Only components which can host a stellar population are meaningful here.
    select case (component%ID)
    case (componentTypeAll%ID,componentTypeDisk%ID,componentTypeSpheroid%ID,componentTypeNuclearStarCluster%ID)
       ! These are supported.
    case default
       call Error_Report(                                                                            &
            &            'component "'                                                            // &
            &            enumerationComponentTypeDecode(component,includePrefix=.false.)          // &
            &            '" can not host a stellar population - use "all", "disk", "spheroid", or'// &
            &            ' "nuclearStarCluster"'                                                  // &
            &            {introspection:location}                                                    &
            &           )
    end select

    allocate(self%luminosityIndex(self%outputTimes_%count()))
    do i=1,self%outputTimes_%count()
       if (present(outputMask).and..not.outputMask(i)) then
          self%luminosityIndex(i)=-1
       else
          self%luminosityIndex(i)=unitStellarLuminosities%index(filterName,filterType,self%outputTimes_%redshift(i),redshiftBand,postprocessChain)
       end if
    end do
    ! Resolve the luminosity index of each age-resolving chain, at each output time. A chain named here must already
    ! be tracked as a luminosity for this filter; `unitStellarLuminosities%index` reports helpfully if it is not.
    countChains=0
    if (present(postprocessChains)) then
       do j=1,size(postprocessChains)
          if (len_trim(postprocessChains(j)) > 0) countChains=countChains+1
       end do
    end if
    allocate(self%postprocessChains   (countChains                          ))
    allocate(self%luminosityIndexChain(countChains,self%outputTimes_%count()))
    if (countChains > 0) then
       countChains=0
       do j=1,size(postprocessChains)
          if (len_trim(postprocessChains(j)) == 0) cycle
          countChains                        =countChains+1
          self%postprocessChains(countChains)=postprocessChains(j)
          do i=1,self%outputTimes_%count()
             if (present(outputMask).and..not.outputMask(i)) then
                self%luminosityIndexChain(countChains,i)=-1
             else
                self%luminosityIndexChain(countChains,i)=unitStellarLuminosities%index(filterName,filterType,self%outputTimes_%redshift(i),redshiftBand,char(postprocessChains(j)))
             end if
          end do
       end do
    end if
    self%name_       ="luminosityStellar:"//filterName//":"//filterType
    if (component == componentTypeAll) then
       self%description_="Total stellar luminosity in the "//filterType//"-frame "//filterName//" filter"
    else
       self%name_       =self%name_       //":"                                            // &
            &            enumerationComponentTypeDecode(component,includePrefix=.false.)
       self%description_="Stellar luminosity of the "                                      // &
            &            enumerationComponentTypeDecode(component,includePrefix=.false.)   // &
            &            " component in the "//filterType//"-frame "//filterName//" filter"
    end if
    if (present(redshiftBand)) then
       write (label,'(f7.3)') redshiftBand
       self%name_       =self%name_        //":z"            //trim(adjustl(label))
       self%description_=self%description_//" shifted to z="//trim(adjustl(label))
    end if
    if (present(postprocessChain)) then
       self%name_       =self%name_                                //postprocessChain
       self%description_=self%description_//" postprocessed with '"//postprocessChain//"'"
    end if
    return
  end function luminosityStellarConstructorInternal

  subroutine luminosityStellarDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorLuminosityStellar` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLuminosityStellar), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine luminosityStellarDestructor

  double precision function luminosityStellarExtract(self,node,instance)
    !!{RST
    Implement a stellar luminosity output analysis property extractor.
    !!}
    use            :: Galactic_Structure_Options, only : massTypeStellar      , weightByLuminosity
    use            :: Galacticus_Nodes          , only : nodeComponentBasic   , treeNode
    use            :: Mass_Distributions        , only : massDistributionClass
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    implicit none
    class  (nodePropertyExtractorLuminosityStellar), intent(inout), target   :: self
    type   (treeNode                              ), intent(inout), target   :: node
    type   (multiCounter                          ), intent(inout), optional :: instance
    class  (massDistributionClass                 )               , pointer  :: massDistribution_
    class  (nodeComponentBasic                    ), pointer                 :: basic
    integer(c_size_t                              )                          :: i
    !$GLC attributes unused :: instance

    basic                    => node             %basic           (                                                                                                                     )
    i                        =  self%outputTimes_%index           (basic%time(),findClosest=.true.                                                                                      )
    massDistribution_        => node             %massDistribution(componentType=self%component,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=self%luminosityIndex(i))
    luminosityStellarExtract =  massDistribution_%massTotal       (                                                                                                                     )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function luminosityStellarExtract


  function luminosityStellarQuantity(self)
    !!{RST
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityLuminosity
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: luminosityStellarQuantity
    class(nodePropertyExtractorLuminosityStellar       ), intent(inout) :: self
    !$GLC attributes unused :: self

    luminosityStellarQuantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function luminosityStellarQuantity

  function luminosityStellarName(self)
    !!{RST
    Return the name of the luminosityStellar property.
    !!}
    implicit none
    type (varying_string                        )                :: luminosityStellarName
    class(nodePropertyExtractorLuminosityStellar), intent(inout) :: self

    luminosityStellarName=self%name_
    return
  end function luminosityStellarName

  function luminosityStellarDescription(self)
    !!{RST
    Return a description of the luminosityStellar property.
    !!}
    implicit none
    type (varying_string                        )                :: luminosityStellarDescription
    class(nodePropertyExtractorLuminosityStellar), intent(inout) :: self

    luminosityStellarDescription=self%description_
    return
  end function luminosityStellarDescription

  double precision function luminosityStellarUnitsInSI(self)
    !!{RST
    Return the units of the luminosityStellar property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminosityZeroPointAB
    implicit none
    class(nodePropertyExtractorLuminosityStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    luminosityStellarUnitsInSI=luminosityZeroPointAB
    return
  end function luminosityStellarUnitsInSI

  function luminosityStellarUnits(self) result(units)
    !!{RST
    Return the units of the luminosityStellar property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType                              )                :: units
    class(nodePropertyExtractorLuminosityStellar), intent(inout) :: self

    units=unitType(self%unitsInSI(),description='AB-magnitude zero point',quantity='4.465920e17 W/Hz')
    return
  end function luminosityStellarUnits

  logical function luminosityStellarSupportsAttenuation(self) result(supportsAttenuation)
    !!{RST
    Return true: a stellar luminosity can be decomposed for attenuation by dust.
    !!}
    implicit none
    class(nodePropertyExtractorLuminosityStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    supportsAttenuation=.true.
    return
  end function luminosityStellarSupportsAttenuation

  function luminosityStellarDecompose(self,node,time,request) result(decomposition)
    !!{RST
    Decompose this luminosity into parcels of emission which may be attenuated separately.

    When no age resolution is requested the whole luminosity is a single parcel. When it is, the luminosity of each
    age bin is recovered from the luminosities computed with the postprocessing chains named by
    ``[postprocessChains]``, whose age windows are read from the chains themselves. Two arrangements of those chains
    are understood:

    *Cumulative*
       Every chain admits ages from zero, each to a different upper limit---the usual case, in which a ``recent`` chain
       restricted to young populations sits alongside an unrestricted ``default`` chain. The luminosity of a bin is
       then the difference between successive chains.

    *Disjoint*
       Each chain admits a range beginning where the previous one ended, so the luminosity of a bin is that of its
       chain directly.

    Anything else, a chain whose age window is not sharp, or a requested age boundary which no chain provides, is
    reported as an error rather than approximated: silently splitting a luminosity at the wrong age would produce
    plausible but wrong colours.
    !!}
    use :: Dust_Attenuation_Descriptors  , only : emissionSourceStellar
    use :: Error                         , only : Error_Report
    use :: Galactic_Structure_Options    , only : componentTypeAll       , enumerationComponentTypeDecode
    use :: Galacticus_Nodes              , only : nodeComponentBasic
    use :: ISO_Varying_String            , only : operator(//)           , var_str
    use :: Numerical_Comparison          , only : Values_Agree
    use :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (luminosityDecomposition               )                              :: decomposition
    class           (nodePropertyExtractorLuminosityStellar), intent(inout), target       :: self
    type            (treeNode                              ), intent(inout), target       :: node
    double precision                                        , intent(in   )               :: time
    type            (decompositionRequest                  ), intent(in   )               :: request
    class           (nodeComponentBasic                    )               , pointer      :: basic
    double precision                                        , allocatable  , dimension(:) :: ageMinimum   , ageMaximum   , &
         &                                                                                   luminosity
    integer                                                 , allocatable  , dimension(:) :: order
    integer         (c_size_t                              )                              :: indexTime
    integer                                                                               :: countChains  , i            , &
         &                                                                                   j            , indexSwap    , &
         &                                                                                   countBins
    logical                                                                               :: isSharp      , isCumulative , &
         &                                                                                   isDisjoint   , found
    double precision                                                                      :: wavelength   , ageLower     , &
         &                                                                                   luminosityBin
    character       (len=12                                )                              :: labelAge

    ! Dust can not be applied to a luminosity summed over components, since components are attenuated differently.
    if (self%component == componentTypeAll)                                                                       &
         & call Error_Report(                                                                                     &
         &                   'this luminosity is summed over all components, and so can not be attenuated; set'// &
         &                   ' [component] to an individual component'                                         // &
         &                   {introspection:location}                                                             &
         &                  )
    basic     => node             %basic(                               )
    indexTime =  self%outputTimes_%index(basic%time(),findClosest=.true.)
    ! With no age resolution requested, the luminosity is a single parcel spanning all ages.
    if (request%countAgeBins() == 1) then
       call decomposition%initialize(1,1)
       decomposition%luminosities(1)              =self%luminosityValue(node,self%luminosityIndex(indexTime))
       decomposition%elementIndex(1)              =1
       decomposition%descriptors (1)%wavelength   =unitStellarLuminosities%wavelengthRestFrame(self%luminosityIndex(indexTime))
       decomposition%descriptors (1)%componentType=self%component
       decomposition%descriptors (1)%sourceType   =emissionSourceStellar
       return
    end if
    ! Age resolution has been requested, so the age-resolving chains are needed.
    countChains=size(self%postprocessChains)
    if (countChains == 0)                                                                                              &
         & call Error_Report(                                                                                          &
         &                   'the dust model requires this luminosity to be resolved by stellar population age, but'// &
         &                   ' no age-resolving postprocessing chains are given; set [postprocessChains]'           // &
         &                   {introspection:location}                                                                  &
         &                  )
    allocate(ageMinimum(countChains))
    allocate(ageMaximum(countChains))
    allocate(luminosity(countChains))
    allocate(order     (countChains))
    do i=1,countChains
       call unitStellarLuminosities%ageWindow(self%luminosityIndexChain(i,indexTime),ageMinimum(i),ageMaximum(i),isSharp)
       if (.not.isSharp)                                                                                                  &
            & call Error_Report(                                                                                          &
            &                   var_str('postprocessing chain "')//self%postprocessChains(i)//'" does not have a sharp'// &
            &                   ' age window, so it can not be used to isolate the light of a range of ages'           // &
            &                   {introspection:location}                                                                  &
            &                  )
       if (ageMaximum(i) <= ageMinimum(i))                                                                                &
            & call Error_Report(                                                                                          &
            &                   var_str('postprocessing chain "')//self%postprocessChains(i)//'" admits no ages at all'// &
            &                   {introspection:location}                                                                  &
            &                  )
       luminosity(i)=self%luminosityValue(node,self%luminosityIndexChain(i,indexTime))
       order     (i)=i
    end do
    ! Order the chains by increasing upper age limit.
    do i=2,countChains
       indexSwap=order(i)
       j        =i-1
       do while (j >= 1)
          if (ageMaximum(order(j)) <= ageMaximum(indexSwap)) exit
          order(j+1)=order(j)
          j         =j-1
       end do
       order(j+1)=indexSwap
    end do
    ! Classify the arrangement of the chains.
    isCumulative=.true.
    isDisjoint  =.true.
    do i=1,countChains
       if (ageMinimum(order(i)) > 0.0d0) isCumulative=.false.
       if (i == 1) then
          if (ageMinimum(order(i)) > 0.0d0) isDisjoint=.false.
       else if (.not.Values_Agree(ageMinimum(order(i)),ageMaximum(order(i-1)),relTol=1.0d-6,absTol=1.0d-12)) then
          isDisjoint=.false.
       end if
    end do
    if (.not.(isCumulative.or.isDisjoint))                                                                                 &
         & call Error_Report(                                                                                              &
         &                   'the age windows of [postprocessChains] are neither cumulative (all beginning at zero age)'// &
         &                   ' nor disjoint (each beginning where the previous ends), so the luminosity of each age bin'// &
         &                   ' can not be determined'                                                                   // &
         &                   {introspection:location}                                                                      &
         &                  )
    ! The highest chain must extend to all ages, or light older than its limit is simply missing.
    if (ageMaximum(order(countChains)) < huge(0.0d0))                                                                     &
         & call Error_Report(                                                                                             &
         &                   'no chain in [postprocessChains] admits the oldest populations, so part of the luminosity'// &
         &                   ' would be discarded; include an unrestricted chain'                                      // &
         &                   {introspection:location}                                                                     &
         &                  )
    ! Every age boundary the dust model asks for must be one that the chains provide.
    do i=1,request%countAgeBins()-1
       found=.false.
       do j=1,countChains-1
          if (Values_Agree(request%ageBoundaries(i),ageMaximum(order(j)),relTol=1.0d-6,absTol=1.0d-12)) then
             found=.true.
             exit
          end if
       end do
       if (.not.found) then
          write (labelAge,'(e12.6)') request%ageBoundaries(i)
          call Error_Report(                                                                              &
               &            var_str('the dust model requires the luminosity to be split at an age of ')// &
               &            trim(adjustl(labelAge))                                                    // &
               &            ' Gyr, but no chain in [postprocessChains] has that upper age limit'       // &
               &            {introspection:location}                                                      &
               &           )
       end if
    end do
    ! Build one parcel per age bin. All parcels contribute to the single output element of this scalar extractor.
    countBins =countChains
    wavelength=unitStellarLuminosities%wavelengthRestFrame(self%luminosityIndexChain(order(countChains),indexTime))
    call decomposition%initialize(countBins,1)
    ageLower=0.0d0
    do i=1,countBins
       if (isCumulative) then
          if (i == 1) then
             luminosityBin=luminosity(order(i))
          else
             luminosityBin=luminosity(order(i))-luminosity(order(i-1))
          end if
       else
          luminosityBin=luminosity(order(i))
       end if
       decomposition%luminosities(i)              =luminosityBin
       decomposition%elementIndex(i)              =1
       decomposition%descriptors (i)%wavelength   =wavelength
       decomposition%descriptors (i)%componentType=self%component
       decomposition%descriptors (i)%sourceType   =emissionSourceStellar
       decomposition%descriptors (i)%ageMinimum   =ageLower
       decomposition%descriptors (i)%ageMaximum   =ageMaximum(order(i))
       ageLower                                   =ageMaximum(order(i))
    end do
    return
  end function luminosityStellarDecompose

  double precision function luminosityStellarLuminosityValue(self,node,index) result(luminosity)
    !!{RST
    Return the luminosity of this extractor's component for the given entry in the stellar luminosities structure.
    !!}
    use :: Galactic_Structure_Options, only : massTypeStellar      , weightByLuminosity
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class  (nodePropertyExtractorLuminosityStellar), intent(inout) :: self
    type   (treeNode                              ), intent(inout) :: node
    integer                                        , intent(in   ) :: index
    class  (massDistributionClass                 ), pointer       :: massDistribution_

    if (index < 0) then
       luminosity=0.0d0
       return
    end if
    massDistribution_ => node             %massDistribution(componentType=self%component,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=index)
    luminosity        =  massDistribution_%massTotal       (                                                                                                   )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function luminosityStellarLuminosityValue

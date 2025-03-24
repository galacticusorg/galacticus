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
Implements a stellar mass output analysis property extractor class.
!!}

  use :: ISO_Varying_String, only : varying_string
  use :: Output_Times      , only : outputTimesClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorLuminosityStellar">
   <description>A stellar luminosity output analysis property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorLuminosityStellar
     !!{
     A stellar luminosity output analysis property extractor class.
     !!}
     private
     type            (varying_string  )                            :: filterName                , filterType, &
          &                                                           postprocessChain          , name_     , &
          &                                                           description_
     double precision                                              :: redshiftBand
     integer                           , allocatable, dimension(:) :: luminosityIndex
     class           (outputTimesClass), pointer                   :: outputTimes_     => null()
   contains
     final     ::                luminosityStellarDestructor
     procedure :: extract     => luminosityStellarExtract
     procedure :: quantity    => luminosityStellarQuantity
     procedure :: name        => luminosityStellarName
     procedure :: description => luminosityStellarDescription
     procedure :: unitsInSI   => luminosityStellarUnitsInSI
  end type nodePropertyExtractorLuminosityStellar

  interface nodePropertyExtractorLuminosityStellar
     !!{
     Constructors for the {\normalfont \ttfamily luminosityStellar} output analysis class.
     !!}
     module procedure luminosityStellarConstructorParameters
     module procedure luminosityStellarConstructorInternal
  end interface nodePropertyExtractorLuminosityStellar

contains

  function luminosityStellarConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily luminosityStellar} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorLuminosityStellar)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (outputTimesClass                      ), pointer       :: outputTimes_
    type            (varying_string                        )                :: filterName           , filterType               , &
         &                                                                     postprocessChain
    double precision                                                        :: redshiftBand
    logical                                                                 :: redshiftBandIsPresent, postprocessChainIsPresent

    redshiftBandIsPresent    =parameters%isPresent('redshiftBand'    )
    postprocessChainIsPresent=parameters%isPresent('postprocessChain')
    !![
    <inputParameter>
      <name>filterName</name>
      <source>parameters</source>
      <description>The filter to select.</description>
    </inputParameter>
    <inputParameter>
      <name>filterType</name>
      <source>parameters</source>
      <description>The filter type (rest or observed) to select.</description>
    </inputParameter>
    !!]
    if (redshiftBandIsPresent) then
       !![
       <inputParameter>
         <name>redshiftBand</name>
         <source>parameters</source>
         <description>The redshift of the band (if not the output redshift).</description>
       </inputParameter>
       !!]
    end if
    if (postprocessChainIsPresent) then
       !![
       <inputParameter>
         <name>postprocessChain</name>
         <source>parameters</source>
         <description>The postprocessing chain to use.</description>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    if (redshiftBandIsPresent) then
       if (postprocessChainIsPresent) then
          self=nodePropertyExtractorLuminosityStellar(char(filterName),char(filterType),outputTimes_,redshiftBand=redshiftBand,postprocessChain=char(postprocessChain))
       else
          self=nodePropertyExtractorLuminosityStellar(char(filterName),char(filterType),outputTimes_,redshiftBand=redshiftBand                                        )
       end if
    else
       if (postprocessChainIsPresent) then
          self=nodePropertyExtractorLuminosityStellar(char(filterName),char(filterType),outputTimes_,                          postprocessChain=char(postprocessChain))
       else
          self=nodePropertyExtractorLuminosityStellar(char(filterName),char(filterType),outputTimes_                                                                  )
       end if
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function luminosityStellarConstructorParameters

  function luminosityStellarConstructorInternal(filterName,filterType,outputTimes_,redshiftBand,postprocessChain,outputMask) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily luminosityStellar} output analysis property extractor class.
    !!}
    use, intrinsic :: ISO_C_Binding                 , only : c_size_t
    use            :: Stellar_Luminosities_Structure, only : unitStellarLuminosities
    implicit none
    type            (nodePropertyExtractorLuminosityStellar)                                        :: self
    character       (len=*                                 ), intent(in   )                         :: filterName      , filterType
    class           (outputTimesClass                      ), intent(in   ), target                 :: outputTimes_
    character       (len=*                                 ), intent(in   ), optional               :: postprocessChain
    double precision                                        , intent(in   ), optional               :: redshiftBand
    logical                                                 , intent(in   ), dimension(:), optional :: outputMask
    integer         (c_size_t                              )                                        :: i
    character       (len=7                                 )                                        :: label
    !![
    <constructorAssign variables="filterName, filterType, redshiftBand, postprocessChain, *outputTimes_"/>
    !!]

    allocate(self%luminosityIndex(self%outputTimes_%count()))
    do i=1,self%outputTimes_%count()
       if (present(outputMask).and..not.outputMask(i)) then
          self%luminosityIndex(i)=-1
       else
          self%luminosityIndex(i)=unitStellarLuminosities%index(filterName,filterType,self%outputTimes_%redshift(i),redshiftBand,postprocessChain)
       end if
    end do
    self%name_       ="luminosityStellar:"//filterName//":"//filterType
    self%description_="Total stellar luminosity luminosity in the "//filterType//"-frame "//filterName//" filter"
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
    !!{
    Destructor for the {\normalfont \ttfamily luminosityStellar} output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorLuminosityStellar), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine luminosityStellarDestructor

  double precision function luminosityStellarExtract(self,node,instance)
    !!{
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

    basic                    => node             %basic           (                                                                                        )
    i                        =  self%outputTimes_%index           (basic%time(),findClosest=.true.                                                         )
    massDistribution_        => node             %massDistribution(massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=self%luminosityIndex(i))
    luminosityStellarExtract =  massDistribution_%massTotal       (                                                                                        )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function luminosityStellarExtract


  function luminosityStellarQuantity(self)
    !!{
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
    !!{
    Return the name of the luminosityStellar property.
    !!}
    implicit none
    type (varying_string                        )                :: luminosityStellarName
    class(nodePropertyExtractorLuminosityStellar), intent(inout) :: self

    luminosityStellarName=self%name_
    return
  end function luminosityStellarName

  function luminosityStellarDescription(self)
    !!{
    Return a description of the luminosityStellar property.
    !!}
    implicit none
    type (varying_string                        )                :: luminosityStellarDescription
    class(nodePropertyExtractorLuminosityStellar), intent(inout) :: self

    luminosityStellarDescription=self%description_
    return
  end function luminosityStellarDescription

  double precision function luminosityStellarUnitsInSI(self)
    !!{
    Return the units of the luminosityStellar property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : luminosityZeroPointAB
    implicit none
    class(nodePropertyExtractorLuminosityStellar), intent(inout) :: self
    !$GLC attributes unused :: self

    luminosityStellarUnitsInSI=luminosityZeroPointAB
    return
  end function luminosityStellarUnitsInSI

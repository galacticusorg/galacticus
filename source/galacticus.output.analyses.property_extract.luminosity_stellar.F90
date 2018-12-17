!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a stellar mass output analysis property extractor class.

  use ISO_Varying_String
  use Output_Times

  !# <outputAnalysisPropertyExtractor name="outputAnalysisPropertyExtractorLuminosityStellar" defaultThreadPrivate="yes">
  !#  <description>A stellar luminosity output analysis property extractor class.</description>
  !# </outputAnalysisPropertyExtractor>
  type, extends(outputAnalysisPropertyExtractorClass) :: outputAnalysisPropertyExtractorLuminosityStellar
     !% A stellar luminosity output analysis property extractor class.
     private
     type            (varying_string  )                            :: filterName      , filterType, &
          &                                                           postprocessChain
     double precision                                              :: redshiftBand
     integer                           , allocatable, dimension(:) :: luminosityIndex
     class           (outputTimesClass), pointer                   :: outputTimes_
   contains
     final     ::             luminosityStellarDestructor
     procedure :: extract  => luminosityStellarExtract
     procedure :: type     => luminosityStellarType
     procedure :: quantity => luminosityStellarQuantity
  end type outputAnalysisPropertyExtractorLuminosityStellar

  interface outputAnalysisPropertyExtractorLuminosityStellar
     !% Constructors for the ``luminosityStellar'' output analysis class.
     module procedure luminosityStellarConstructorParameters
     module procedure luminosityStellarConstructorInternal
  end interface outputAnalysisPropertyExtractorLuminosityStellar

contains

  function luminosityStellarConstructorParameters(parameters) result(self)
    !% Constructor for the ``luminosityStellar'' output analysis property extractor class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisPropertyExtractorLuminosityStellar)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (outputTimesClass                                ), pointer       :: outputTimes_
    type            (varying_string                                  )                :: filterName           , filterType               , &
         &                                                                               postprocessChain
    double precision                                                                  :: redshiftBand
    logical                                                                           :: redshiftBandIsPresent, postprocessChainIsPresent

    redshiftBandIsPresent    =parameters%isPresent('redshiftBand'    )
    postprocessChainIsPresent=parameters%isPresent('postprocessChain')
    !# <inputParameter>
    !#   <name>filterName</name>
    !#   <source>parameters</source>
    !#   <description>The filter to select.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>filterType</name>
    !#   <source>parameters</source>
    !#   <description>The filter type (rest or observed) to select.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    if (redshiftBandIsPresent) then
       !# <inputParameter>
       !#   <name>redshiftBand</name>
       !#   <source>parameters</source>
       !#   <description>The redshift of the band (if not the output redshift).</description>
       !#   <type>float</type>
       !#   <cardinality>0..1</cardinality>
       !# </inputParameter>
    end if
    if (postprocessChainIsPresent) then
       !# <inputParameter>
       !#   <name>postprocessChain</name>
       !#   <source>parameters</source>
       !#   <description>The postprocessing chain to use.</description>
       !#   <type>string</type>
       !#   <cardinality>0..1</cardinality>
       !# </inputParameter>
    end if
    !# <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    if (redshiftBandIsPresent) then
       if (postprocessChainIsPresent) then
          self=outputAnalysisPropertyExtractorLuminosityStellar(char(filterName),char(filterType),outputTimes_,redshiftBand=redshiftBand,postprocessChain=char(postprocessChain))
       else
          self=outputAnalysisPropertyExtractorLuminosityStellar(char(filterName),char(filterType),outputTimes_,redshiftBand=redshiftBand                                        )
       end if
    else
       if (postprocessChainIsPresent) then
          self=outputAnalysisPropertyExtractorLuminosityStellar(char(filterName),char(filterType),outputTimes_,                          postprocessChain=char(postprocessChain))
       else
          self=outputAnalysisPropertyExtractorLuminosityStellar(char(filterName),char(filterType),outputTimes_                                                                  )
       end if
    end if
    !# <inputParametersValidate source="parameters"/>
    return
  end function luminosityStellarConstructorParameters
  
  function luminosityStellarConstructorInternal(filterName,filterType,outputTimes_,redshiftBand,postprocessChain,outputMask) result(self)
    !% Internal constructor for the ``luminosityStellar'' output analysis property extractor class.
    use, intrinsic :: ISO_C_Binding
    use               Stellar_Luminosities_Structure
    use               Memory_Management
    implicit none
    type            (outputAnalysisPropertyExtractorLuminosityStellar)                                        :: self
    character       (len=*                                           ), intent(in   )                         :: filterName      , filterType
    class           (outputTimesClass                                ), intent(in   ), target                 :: outputTimes_
    character       (len=*                                           ), intent(in   ), optional               :: postprocessChain
    double precision                                                  , intent(in   ), optional               :: redshiftBand
    logical                                                           , intent(in   ), dimension(:), optional :: outputMask
    integer         (c_size_t                                        )                                        :: i
    !# <constructorAssign variables="filterName, filterType, redshiftBand, postprocessChain, *outputTimes_"/>

    call allocateArray(self%luminosityIndex,[self%outputTimes_%count()])
    do i=1,self%outputTimes_%count()
       if (present(outputMask).and..not.outputMask(i)) then
          self%luminosityIndex(i)=-1
       else
          self%luminosityIndex(i)=unitStellarLuminosities%index(filterName,filterType,self%outputTimes_%redshift(i),redshiftBand,postprocessChain)
       end if
    end do
    return
  end function luminosityStellarConstructorInternal
  
  subroutine luminosityStellarDestructor(self)
    !% Destructor for the ``luminosityStellar'' output analysis property extractor class.
    use               Memory_Management
    implicit none
    type(outputAnalysisPropertyExtractorLuminosityStellar), intent(inout) :: self

    !# <objectDestructor name="self%outputTimes_"/>
    return
  end subroutine luminosityStellarDestructor

  double precision function luminosityStellarExtract(self,node)
    !% Implement a stellar luminosity output analysis property extractor.
    use, intrinsic :: ISO_C_Binding
    use               Galactic_Structure_Enclosed_Masses
    use               Galactic_Structure_Options
    use               Galacticus_Nodes                  , only : nodeComponentBasic
    implicit none
    class  (outputAnalysisPropertyExtractorLuminosityStellar), intent(inout) :: self
    type   (treeNode                                        ), intent(inout) :: node
    class  (nodeComponentBasic                              ), pointer       :: basic
    integer(c_size_t                                        )                :: i
    
    basic                    =>                                  node %basic()
    i                        =  self%outputTimes_%index         (basic%time ()                                                                                                     )
    luminosityStellarExtract =  Galactic_Structure_Enclosed_Mass(node         ,radiusLarge,massType=massTypeStellar,weightBy=weightByLuminosity,weightIndex=self%luminosityIndex(i))
    return
  end function luminosityStellarExtract
  
  integer function luminosityStellarType(self)
    !% Return the type of the stellar luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorLuminosityStellar), intent(inout) :: self
    !GCC$ attributes unused :: self

    luminosityStellarType=outputAnalysisPropertyTypeLinear
    return
  end function luminosityStellarType

  integer function luminosityStellarQuantity(self)
    !% Return the class of the stellar luminosity property.
    use Output_Analyses_Options
    implicit none
    class(outputAnalysisPropertyExtractorLuminosityStellar), intent(inout) :: self
    !GCC$ attributes unused :: self

    luminosityStellarQuantity=outputAnalysisPropertyQuantityLuminosity
    return
  end function luminosityStellarQuantity

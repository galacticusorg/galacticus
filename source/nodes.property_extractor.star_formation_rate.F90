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
Implements a star formation rate property extractor class.
!!}

  use :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorStarFormationRate">
   <description>
    A node property extractor which extracts the star formation rate in a galaxy. The type of star formation rate is controlled by
    the {\normalfont \ttfamily [component]} parameter, which can be either ``{\normalfont \ttfamily disk}'', ``{\normalfont
    \ttfamily spheroid}'', ``{\normalfont \ttfamily nsc}'' or ``{\normalfont \ttfamily total}''. The corresponding star formation
    rate is extracted as {\normalfont \ttfamily \textless\ component\textgreater\ StarFormationRate} in units of $M_\odot$/Gyr.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorStarFormationRate
     !!{
     A star formation rate property extractor class.
     !!}
     private
     class(starFormationRateDisksClass              ), pointer :: starFormationRateDisks_               => null()
     class(starFormationRateSpheroidsClass          ), pointer :: starFormationRateSpheroids_           => null()
     class(starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     type (varying_string                           )          :: name_                                          , description_, &
          &                                                       component
   contains
     final     ::                starFormationRateDestructor
     procedure :: extract     => starFormationRateExtract
     procedure :: quantity    => starFormationRateQuantity
     procedure :: name        => starFormationRateName
     procedure :: description => starFormationRateDescription
     procedure :: unitsInSI   => starFormationRateUnitsInSI
  end type nodePropertyExtractorStarFormationRate

  interface nodePropertyExtractorStarFormationRate
     !!{
     Constructors for the \refClass{nodePropertyExtractorStarFormationRate} output analysis class.
     !!}
     module procedure starFormationRateConstructorParameters
     module procedure starFormationRateConstructorInternal
  end interface nodePropertyExtractorStarFormationRate

contains

  function starFormationRateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorStarFormationRate} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorStarFormationRate   )                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(starFormationRateDisksClass              ), pointer       :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass          ), pointer       :: starFormationRateSpheroids_
    class(starFormationRateNuclearStarClustersClass), pointer       :: starFormationRateNuclearStarClusters_
    type (varying_string                           )                :: component
    type (enumerationGalacticComponentType         )                :: component_

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component from which to extract star formation rate.</description>
    </inputParameter>
    !!]
    component_=enumerationGalacticComponentEncode(char(component),includesPrefix=.false.)
    select case (component_%ID)
    case (galacticComponentDisk              %ID)
       !![
       <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
       !!]
       starFormationRateSpheroids_           => null()
       starFormationRateNuclearStarClusters_ => null()
    case (galacticComponentSpheroid          %ID)
       !![
       <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
       !!]
       starFormationRateDisks_               => null()
       starFormationRateNuclearStarClusters_ => null()
    case (galacticComponentNuclearStarCluster%ID)
       !![
       <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
       !!]
       starFormationRateDisks_              => null()
       starFormationRateSpheroids_          => null()
    case (galacticComponentTotal            %ID)
       !![
       <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
       <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
       <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
       !!]
    end select
    self=nodePropertyExtractorStarFormationRate(starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function starFormationRateConstructorParameters

  function starFormationRateConstructorInternal(starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorStarFormationRate} property extractor class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorStarFormationRate   )                        :: self
    class(starFormationRateDisksClass              ), intent(in   ), target :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass          ), intent(in   ), target :: starFormationRateSpheroids_
    class(starFormationRateNuclearStarClustersClass), intent(in   ), target :: starFormationRateNuclearStarClusters_
    !![
    <constructorAssign variables="*starFormationRateDisks_, *starFormationRateSpheroids_, *starFormationRateNuclearStarClusters_"/>
    !!]

    if      (associated(self%starFormationRateDisks_).and.associated(self%starFormationRateSpheroids_).and.associated(self%starFormationRateNuclearStarClusters_)) then
       self%name_       ="totalStarFormationRate"
       self%description_="Total (disk + spheroid + nuclear star cluster) star formation rate [M☉ Gyr⁻¹]."
       self%component   ="total"
    else if (associated(self%starFormationRateDisks_)                                                                                                            ) then
       self%name_       ="diskStarFormationRate"
       self%description_="Disk star formation rate [M☉ Gyr⁻¹]."
       self%component   ="disk"
    else if (                                             associated(self%starFormationRateSpheroids_)                                                           ) then
       self%name_       ="spheroidStarFormationRate"
       self%description_="Spheroid star formation rate [M☉ Gyr⁻¹]."
       self%component   ="spheroid"
    else if (                                                                                              associated(self%starFormationRateNuclearStarClusters_)) then
       self%name_       ="nuclearStarClusterStarFormationRate"
       self%description_="Nuclear star cluster star formation rate [M☉ Gyr⁻¹]."
       self%component   ="NSC"
    else
       call Error_Report('No star formation rate specified.'//{introspection:location})
    end if
    return
  end function starFormationRateConstructorInternal

  subroutine starFormationRateDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorStarFormationRate} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorStarFormationRate), intent(inout) :: self
  
    !![
    <objectDestructor name="self%starFormationRateDisks_"              />
    <objectDestructor name="self%starFormationRateSpheroids_"          />
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine starFormationRateDestructor

  double precision function starFormationRateExtract(self,node,instance)
    !!{
    Implement a star formation rate output analysis property extractor.
    !!}
    implicit none
    class(nodePropertyExtractorStarFormationRate), intent(inout), target   :: self
    type (treeNode                              ), intent(inout), target   :: node
    type (multiCounter                          ), intent(inout), optional :: instance
    !$GLC attributes unused :: instance

    starFormationRateExtract=0.0d0
    if (associated(self%starFormationRateDisks_              )) starFormationRateExtract=starFormationRateExtract+self%starFormationRateDisks_              %rate(node)
    if (associated(self%starFormationRateSpheroids_          )) starFormationRateExtract=starFormationRateExtract+self%starFormationRateSpheroids_          %rate(node)
    if (associated(self%starFormationRateNuclearStarClusters_)) starFormationRateExtract=starFormationRateExtract+self%starFormationRateNuclearStarClusters_%rate(node)
    return
  end function starFormationRateExtract

  function starFormationRateQuantity(self)
    !!{
    Return the class of the stellar luminosity property.
    !!}
    use :: Output_Analyses_Options, only : outputAnalysisPropertyQuantityStarFormationRate
    implicit none
    type (enumerationOutputAnalysisPropertyQuantityType)                :: starFormationRateQuantity
    class(nodePropertyExtractorStarFormationRate       ), intent(inout) :: self
    !$GLC attributes unused :: self

    starFormationRateQuantity=outputAnalysisPropertyQuantityStarFormationRate
    return
  end function starFormationRateQuantity

  function starFormationRateName(self)
    !!{
    Return the name of the starFormationRate property.
    !!}
    implicit none
    type (varying_string                        )                :: starFormationRateName
    class(nodePropertyExtractorStarFormationRate), intent(inout) :: self

    starFormationRateName=self%name_
    return
  end function starFormationRateName

  function starFormationRateDescription(self)
    !!{
    Return a description of the starFormationRate property.
    !!}
    implicit none
    type (varying_string                        )                :: starFormationRateDescription
    class(nodePropertyExtractorStarFormationRate), intent(inout) :: self

    starFormationRateDescription=self%description_
    return
  end function starFormationRateDescription

  double precision function starFormationRateUnitsInSI(self)
    !!{
    Return the units of the starFormationRate property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    class(nodePropertyExtractorStarFormationRate), intent(inout) :: self
    !$GLC attributes unused :: self

    starFormationRateUnitsInSI=massSolar/gigaYear
    return
  end function starFormationRateUnitsInSI

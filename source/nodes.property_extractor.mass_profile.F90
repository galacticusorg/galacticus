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
  Contains a module which implements a property extractor class for the enclosed mass at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale   , darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier
  use :: Galactic_Structure                  , only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassProfile">
   <description>A property extractor class for the enclosed mass at a set of radii.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorMassProfile
     !!{
     A property extractor class for the enclosed mass at a set of radii.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_          => null()
     class  (galacticStructureClass  ), pointer                   :: galacticStructure_            => null()
     integer                                                      :: radiiCount                             , elementCount_
     logical                                                      :: includeRadii
     type   (varying_string          ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusSpecifier         ), allocatable, dimension(:) :: radii
     logical                                                      :: darkMatterScaleRadiusIsNeeded          , diskIsNeeded        , &
          &                                                          spheroidIsNeeded                       , virialRadiusIsNeeded
   contains
     final     ::                       massProfileDestructor
     procedure :: columnDescriptions => massProfileColumnDescriptions
     procedure :: size               => massProfileSize
     procedure :: elementCount       => massProfileElementCount
     procedure :: extract            => massProfileExtract
     procedure :: names              => massProfileNames
     procedure :: descriptions       => massProfileDescriptions
     procedure :: unitsInSI          => massProfileUnitsInSI
  end type nodePropertyExtractorMassProfile

  interface nodePropertyExtractorMassProfile
     !!{
     Constructors for the ``massProfile'' output analysis class.
     !!}
     module procedure massProfileConstructorParameters
     module procedure massProfileConstructorInternal
  end interface nodePropertyExtractorMassProfile

contains

  function massProfileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily massProfile} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorMassProfile)                              :: self
    type   (inputParameters                 ), intent(inout)               :: parameters
    type   (varying_string                  ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass        ), pointer                     :: darkMatterHaloScale_
    class  (galacticStructureClass          ), pointer                     :: galacticStructure_
    logical                                                                :: includeRadii

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter>
      <name>radiusSpecifiers</name>
      <description>A list of radius specifiers at which to output the enclosed mass profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the radii at which enclosed mass data are output should also be included in the output file.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="galacticStructure"   name="galacticStructure_"   source="parameters"/>
    !!]
    self=nodePropertyExtractorMassProfile(radiusSpecifiers,includeRadii,darkMatterHaloScale_,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="galacticStructure_"  />
    !!]
    return
  end function massProfileConstructorParameters

  function massProfileConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily massProfile} property extractor class.
    !!}
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorMassProfile)                              :: self
    type   (varying_string                  ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass        ), intent(in   ), target       :: darkMatterHaloScale_
    class  (galacticStructureClass          ), intent(in   ), target       :: galacticStructure_
    logical                                  , intent(in   )               :: includeRadii
    !![
    <constructorAssign variables="radiusSpecifiers, includeRadii, *darkMatterHaloScale_, *galacticStructure_"/>
    !!]

    if (includeRadii) then
       self%elementCount_=2
    else
       self%elementCount_=1
    end if
    self%radiiCount      =size(radiusSpecifiers)
    call Galactic_Structure_Radii_Definition_Decode(                                    &
         &                                          radiusSpecifiers                  , &
         &                                          self%radii                        , &
         &                                          self%diskIsNeeded                 , &
         &                                          self%spheroidIsNeeded             , &
         &                                          self%virialRadiusIsNeeded         , &
         &                                          self%darkMatterScaleRadiusIsNeeded  &
         &                                         )
    return
  end function massProfileConstructorInternal

  subroutine massProfileDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily massProfile} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMassProfile), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%galacticStructure_"  />
    !!]
    return
  end subroutine massProfileDestructor

  integer function massProfileElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily massProfile} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorMassProfile), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    massProfileElementCount=self%elementCount_
    return
  end function massProfileElementCount

  function massProfileSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily massProfile} property extractors.
    !!}
    implicit none
    integer         (c_size_t                        )                :: massProfileSize
    class           (nodePropertyExtractorMassProfile), intent(inout) :: self
    double precision                                  , intent(in   ) :: time
    !$GLC attributes unused :: time

    massProfileSize=self%radiiCount
    return
  end function massProfileSize

  function massProfileExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily massProfile} property extractor.
    !!}
    use :: Galactic_Structure_Options          , only : componentTypeAll               , massTypeGalactic            , massTypeStellar
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius, radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius            , radiusTypeGalacticLightFraction, &
          &                                             radiusTypeGalacticMassFraction , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius, radiusTypeSpheroidRadius       , &
          &                                             radiusTypeStellarMassFraction  , radiusTypeVirialRadius
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile , nodeComponentDisk           , nodeComponentSpheroid           , treeNode
    implicit none
    double precision                                  , dimension(:,:), allocatable :: massProfileExtract
    class           (nodePropertyExtractorMassProfile), intent(inout) , target      :: self
    type            (treeNode                        ), intent(inout) , target      :: node
    double precision                                  , intent(in   )               :: time
    type            (multiCounter                    ), intent(inout) , optional    :: instance
    class           (nodeComponentDisk               ), pointer                     :: disk
    class           (nodeComponentSpheroid           ), pointer                     :: spheroid
    class           (nodeComponentDarkMatterProfile  ), pointer                     :: darkMatterProfile
    integer                                                                         :: i
    double precision                                                                :: radius             , radiusVirial
    !$GLC attributes unused :: time, instance

    allocate(massProfileExtract(self%radiiCount,self%elementCount_))
    radiusVirial                                              =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial      =  self%darkMatterHaloScale_%radiusVirial(node                    )
    if (self%                 diskIsNeeded) disk              =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid          =>                                        node%spheroid         ()
    if (self%darkMatterScaleRadiusIsNeeded) darkMatterProfile =>                                        node%darkMatterProfile()
    do i=1,self%radiiCount
       radius=self%radii(i)%value
       select case (self%radii(i)%type%ID)
       case   (radiusTypeRadius                %ID)
          ! Nothing to do.
       case   (radiusTypeVirialRadius          %ID)
          radius=+radius*radiusVirial
       case   (radiusTypeDarkMatterScaleRadius %ID)
          radius=+radius*darkMatterProfile%         scale()
       case   (radiusTypeDiskRadius            %ID)
          radius=+radius*disk             %        radius()
       case   (radiusTypeSpheroidRadius        %ID)
          radius=+radius*spheroid         %        radius()
       case   (radiusTypeDiskHalfMassRadius    %ID)
          radius=+radius*disk             %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius%ID)
          radius=+radius*spheroid         %halfMassRadius()
       case   (radiusTypeGalacticMassFraction  %ID,  &
            &  radiusTypeGalacticLightFraction %ID)
          radius=+radius                                           &
               & *self%galacticStructure_%radiusEnclosingMass      &
               &  (                                                &
               &   node                                         ,  &
               &   massFractional=self%radii(i)%fraction        ,  &
               &   massType      =              massTypeGalactic,  &
               &   componentType =              componentTypeAll,  &
               &   weightBy      =self%radii(i)%weightBy        ,  &
               &   weightIndex   =self%radii(i)%weightByIndex      &
               &  )
        case   (radiusTypeStellarMassFraction  %ID)
          radius=+radius                                           &
               & *self%galacticStructure_%radiusEnclosingMass      &
               &  (                                                &
               &   node                                         ,  &
               &   massFractional=self%radii(i)%fraction        ,  &
               &   massType      =              massTypeStellar ,  &
               &   componentType =              componentTypeAll,  &
               &   weightBy      =self%radii(i)%weightBy        ,  &
               &   weightIndex   =self%radii(i)%weightByIndex      &
               &  )
       end select
       massProfileExtract       (i,1)=self%galacticStructure_%massEnclosed(                                       &
            &                                                              node                                 , &
            &                                                              radius                               , &
            &                                                              componentType=self%radii(i)%component, &
            &                                                              massType     =self%radii(i)%mass       &
            &                                                             )
       if (self%includeRadii)                                                                                     &
            & massProfileExtract(i,2)=                                     radius
    end do
    return
  end function massProfileExtract

  subroutine massProfileNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily massProfile} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorMassProfile), intent(inout)                             :: self
    double precision                                  , intent(in   ), optional                   :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names(1)="massProfile"
    if (self%includeRadii) names(2)="massProfileRadius"
    return
  end subroutine massProfileNames

  subroutine massProfileDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily massProfile} property.
    !!}
    implicit none
    class           (nodePropertyExtractorMassProfile), intent(inout)                             :: self
    double precision                                  , intent(in   ), optional                   :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions       (1)="Enclosed mass at a given radius [M☉]."
    if (self%includeRadii)                                                  &
         & descriptions(2)="Radius at which enclosed mass is output [Mpc]."
    return
  end subroutine massProfileDescriptions

  subroutine massProfileColumnDescriptions(self,descriptions,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily massProfile} property.
    !!}
    implicit none
    class           (nodePropertyExtractorMassProfile), intent(inout)                             :: self
    double precision                                  , intent(in   ), optional                   :: time
    type            (varying_string                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%radiiCount))
    descriptions=self%radii%name
    return
  end subroutine massProfileColumnDescriptions

  function massProfileUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily massProfile} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                  , allocatable  , dimension(:) :: massProfileUnitsInSI
    class           (nodePropertyExtractorMassProfile), intent(inout)               :: self
    double precision                                  , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(massProfileUnitsInSI(self%elementCount_))
    massProfileUnitsInSI       (1)=massSolar
    if (self%includeRadii)                    &
         & massProfileUnitsInSI(2)=megaParsec
    return
  end function massProfileUnitsInSI


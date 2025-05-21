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
  Implements a property extractor class for the enclosed mass at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier
  use :: Cosmology_Parameters                , only : cosmologyParametersClass

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
     class  (cosmologyParametersClass), pointer                   :: cosmologyParameters_          => null()
     integer                                                      :: radiiCount                             , elementCount_
     logical                                                      :: includeRadii
     type   (varying_string          ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusSpecifier         ), allocatable, dimension(:) :: radii

     logical                                                      :: darkMatterScaleRadiusIsNeeded          , diskIsNeeded        , &
          &                                                          spheroidIsNeeded                       , virialRadiusIsNeeded, &
          &                                                          nuclearStarClusterIsNeeded             , satelliteIsNeeded   , &
          &                                                          hotHaloIsNeeded
     double precision                                             :: fractionDarkMatter
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
     Constructors for the \refClass{nodePropertyExtractorMassProfile} output analysis class.
     !!}
     module procedure massProfileConstructorParameters
     module procedure massProfileConstructorInternal
  end interface nodePropertyExtractorMassProfile

contains

  function massProfileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMassProfile} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorMassProfile)                              :: self
    type   (inputParameters                 ), intent(inout)               :: parameters
    type   (varying_string                  ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass        ), pointer                     :: darkMatterHaloScale_
    class  (cosmologyParametersClass        ), pointer                     :: cosmologyParameters_
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
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=nodePropertyExtractorMassProfile(radiusSpecifiers,includeRadii,darkMatterHaloScale_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function massProfileConstructorParameters

  function massProfileConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMassProfile} property extractor class.
    !!}
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorMassProfile)                              :: self
    type   (varying_string                  ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass        ), intent(in   ), target       :: darkMatterHaloScale_
    class  (cosmologyParametersClass        ), intent(in   ), target       :: cosmologyParameters_
    logical                                  , intent(in   )               :: includeRadii
    !![
    <constructorAssign variables="radiusSpecifiers, includeRadii, *darkMatterHaloScale_, *cosmologyParameters_"/>
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
         &                                          self%hotHaloIsNeeded              , &
         &                                          self%diskIsNeeded                 , &
         &                                          self%spheroidIsNeeded             , &
         &                                          self%nuclearStarClusterIsNeeded   , &
         &                                          self%satelliteIsNeeded            , &
         &                                          self%virialRadiusIsNeeded         , &
         &                                          self%darkMatterScaleRadiusIsNeeded  &
         &                                         )
    self%fractionDarkMatter=+(                                         &
         &                    +self%cosmologyParameters_%OmegaMatter() &
         &                    -self%cosmologyParameters_%OmegaBaryon() &
         &                   )                                         &
         &                  /  self%cosmologyParameters_%OmegaMatter()
    return
  end function massProfileConstructorInternal

  subroutine massProfileDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorMassProfile} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMassProfile), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%cosmologyParameters_"/>
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

    use :: Galactic_Structure_Options          , only : componentTypeAll                          , massTypeGalactic            , massTypeStellar                     , massTypeDark
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius           , radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius                , radiusTypeGalacticLightFraction   , &
         &                                              radiusTypeGalacticMassFraction            , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius    , radiusTypeSpheroidRadius          , &
         &                                              radiusTypeStellarMassFraction             , radiusTypeVirialRadius      , radiusTypeSatelliteBoundMassFraction, radiusTypeNuclearStarClusterRadius, &
         &                                              radiusTypeNuclearStarClusterHalfMassRadius, radiusTypeHotHaloOuterRadius
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile            , nodeComponentDisk           , nodeComponentSpheroid               , nodeComponentNSC                  , & 
         &                                              nodeComponentSatellite                    , nodeComponentHotHalo        , treeNode
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Error                               , only : Error_Report
    implicit none
    double precision                                  , dimension(:,:), allocatable :: massProfileExtract
    class           (nodePropertyExtractorMassProfile), intent(inout) , target      :: self
    type            (treeNode                        ), intent(inout) , target      :: node
    double precision                                  , intent(in   )               :: time
    type            (multiCounter                    ), intent(inout) , optional    :: instance
    class           (nodeComponentHotHalo            ), pointer                     :: hotHalo
    class           (nodeComponentDisk               ), pointer                     :: disk
    class           (nodeComponentSpheroid           ), pointer                     :: spheroid
    class           (nodeComponentNSC                ), pointer                     :: nuclearStarCluster
    class           (nodeComponentDarkMatterProfile  ), pointer                     :: darkMatterProfile
    class           (nodeComponentSatellite          ), pointer                     :: satellite
    class           (massDistributionClass           ), pointer                     :: massDistribution_
    integer                                                                         :: i
    double precision                                                                :: radius             , radiusVirial, &
         &                                                                             mass
    !$GLC attributes unused :: time, instance

    allocate(massProfileExtract(self%radiiCount,self%elementCount_))
    radiusVirial                                               =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial       =  self%darkMatterHaloScale_%radiusVirial(node                    )
    if (self%              hotHaloIsNeeded) hotHalo            =>                                        node%hotHalo          ()
    if (self%                 diskIsNeeded) disk               =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid           =>                                        node%spheroid         ()
    if (self%   nuclearStarClusterIsNeeded) nuclearStarCluster =>                                        node%NSC              ()
    if (self%            satelliteIsNeeded) satellite          =>                                        node%satellite        ()
    if (self%darkMatterScaleRadiusIsNeeded) darkMatterProfile  =>                                        node%darkMatterProfile()
    do i=1,self%radiiCount
       radius=self%radii(i)%value
       select case (self%radii(i)%type%ID)
       case   (radiusTypeRadius                          %ID)
          ! Nothing to do.
       case   (radiusTypeVirialRadius                    %ID)
          radius=+radius*radiusVirial
       case   (radiusTypeDarkMatterScaleRadius           %ID)
          radius=+radius*darkMatterProfile %         scale()
       case   (radiusTypeHotHaloOuterRadius              %ID)
          radius=+radius*hotHalo           %   outerRadius()
       case   (radiusTypeDiskRadius                      %ID)
          radius=+radius*disk              %        radius()
       case   (radiusTypeSpheroidRadius                  %ID)
          radius=+radius*spheroid          %        radius()
       case  (radiusTypeNuclearStarClusterRadius         %ID)
          radius=+radius*nuclearStarCluster%        radius()
       case   (radiusTypeDiskHalfMassRadius              %ID)
          radius=+radius*disk              %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius%ID)
          radius=+radius*spheroid          %halfMassRadius()
       case   (radiusTypeNuclearStarClusterHalfMassRadius%ID)
          radius=+radius*nuclearStarCluster%halfMassRadius()
       case   (radiusTypeSatelliteBoundMassFraction      %ID)
          mass              =  +satellite        %boundMass          (                                                &
               &                                                     )                                                &
               &               *self             %fractionDarkMatter
          massDistribution_ =>  node             %massDistribution   (                                                &
               &                                                      massType      =              massTypeDark    ,  &
               &                                                      componentType =              componentTypeAll,  &
               &                                                      weightBy      =self%radii(i)%weightBy        ,  &
               &                                                      weightIndex   =self%radii(i)%weightByIndex      &
               &                                                     )
          radius            =  +radius                                                                                &
               &               *massDistribution_%radiusEnclosingMass(                                                &
               &                                                      mass          =              mass               &
               &                                                     )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       case   (radiusTypeGalacticMassFraction  %ID,  &
            &  radiusTypeGalacticLightFraction %ID)
          massDistribution_ =>  node             %massDistribution   (                                                &
               &                                                      massType      =              massTypeStellar ,  &
               &                                                      componentType =              componentTypeAll,  &
               &                                                      weightBy      =self%radii(i)%weightBy        ,  &
               &                                                      weightIndex   =self%radii(i)%weightByIndex      &
               &                                                     )
          radius            =  +radius                                                                                &
               &               *massDistribution_%radiusEnclosingMass(                                                &
               &                                                      massFractional=self%radii(i)%fraction           &
               &                                                     )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       case   (radiusTypeStellarMassFraction  %ID)
           massDistribution_ =>  node             %massDistribution  (                                                &
               &                                                      massType      =              massTypeStellar ,  &
               &                                                      componentType =              componentTypeAll,  &
               &                                                      weightBy      =self%radii(i)%weightBy        ,  &
               &                                                      weightIndex   =self%radii(i)%weightByIndex      &
               &                                                     )
          radius            =  +radius                                                                                &
               &               *massDistribution_%radiusEnclosingMass(                                                &
               &                                                      massFractional=self%radii(i)%fraction           &
               &                                                     )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       case default
          call Error_Report('unrecognized radius type'//{introspection:location})
       end select
       massDistribution_              => node             %massDistribution    (                                       &
            &                                                                   componentType=self%radii(i)%component, &
            &                                                                   massType     =self%radii(i)%mass       &
            &                                                                  )
       massProfileExtract       (i,1) =  massDistribution_%massEnclosedBySphere(                                       &
            &                                                                   radius       =              radius     &
            &                                                                  )
       if (self%includeRadii)                                                                                          &
            & massProfileExtract(i,2) =                                                                     radius
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
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

  subroutine massProfileColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily massProfile} property.
    !!}
    implicit none
    class           (nodePropertyExtractorMassProfile), intent(inout)                            :: self
    double precision                                  , intent(in   ), optional                  :: time
    type            (varying_string                  ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                  , intent(inout), dimension(:), allocatable :: values 
    type            (varying_string                  ), intent(  out)                            :: valuesDescription
    double precision                                  , intent(  out)                            :: valuesUnitsInSI
    !$GLC attributes unused :: time

    allocate(descriptions(self%radiiCount))
    allocate(values      (              0))
    valuesDescription=var_str('')
    valuesUnitsInSI  =0.0d0
    descriptions     =self%radii%name
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


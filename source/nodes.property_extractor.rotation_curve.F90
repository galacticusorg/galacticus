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
  Contains a module which implements a property extractor class for the rotation curve at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale, darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRotationCurve">
   <description>
    A property extractor class for the rotation curve at a set of radii. The radii and types of rotation curve to output
    are specified by the {\normalfont \ttfamily radiusSpecifiers} parameter. This parameter's value can contain multiple
    entries, each of which should be a valid
    \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#sec.radiusSpecifiers}{radius
    specifier}.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorRotationCurve
     !!{
     A property extractor class for the rotation curve at a set of radii.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_          => null()
     integer                                                      :: radiiCount                             , elementCount_
     logical                                                      :: includeRadii
     type   (varying_string          ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusSpecifier         ), allocatable, dimension(:) :: radii
     logical                                                      :: darkMatterScaleRadiusIsNeeded          , diskIsNeeded        , &
          &                                                          spheroidIsNeeded                       , virialRadiusIsNeeded, & 
          &                                                          nuclearStarClusterIsNeeded             , satelliteIsNeeded
   contains
     final     ::                       rotationCurveDestructor
     procedure :: columnDescriptions => rotationCurveColumnDescriptions
     procedure :: size               => rotationCurveSize
     procedure :: elementCount       => rotationCurveElementCount
     procedure :: extract            => rotationCurveExtract
     procedure :: names              => rotationCurveNames
     procedure :: descriptions       => rotationCurveDescriptions
     procedure :: unitsInSI          => rotationCurveUnitsInSI
  end type nodePropertyExtractorRotationCurve

  interface nodePropertyExtractorRotationCurve
     !!{
     Constructors for the ``rotationCurve'' output analysis class.
     !!}
     module procedure rotationCurveConstructorParameters
     module procedure rotationCurveConstructorInternal
  end interface nodePropertyExtractorRotationCurve

contains

  function rotationCurveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily rotationCurve} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorRotationCurve)                              :: self
    type   (inputParameters                   ), intent(inout)               :: parameters
    type   (varying_string                    ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass          ), pointer                     :: darkMatterHaloScale_
    logical                                                                  :: includeRadii

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter>
      <name>radiusSpecifiers</name>
      <description>A list of radius specifiers at which to output the rotation curve.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the radii at which rotation curve data are output should also be included in the output file.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRotationCurve(radiusSpecifiers,includeRadii,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function rotationCurveConstructorParameters

  function rotationCurveConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily rotationCurve} property extractor class.
    !!}
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorRotationCurve)                              :: self
    type   (varying_string                    ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass          ), intent(in   ), target       :: darkMatterHaloScale_
    logical                                    , intent(in   )               :: includeRadii
    !![
    <constructorAssign variables="radiusSpecifiers, includeRadii, *darkMatterHaloScale_"/>
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
         &                                          self%nuclearStarClusterIsNeeded   , &
         &                                          self%satelliteIsNeeded            , &
         &                                          self%virialRadiusIsNeeded         , &
         &                                          self%darkMatterScaleRadiusIsNeeded  &
         &                                         )
    return
  end function rotationCurveConstructorInternal

  subroutine rotationCurveDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily rotationCurve} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorRotationCurve), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine rotationCurveDestructor

  integer function rotationCurveElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily rotationCurve} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorRotationCurve), intent(inout) :: self
    double precision                                    , intent(in   ) :: time
    !$GLC attributes unused :: time

    rotationCurveElementCount=self%elementCount_
    return
  end function rotationCurveElementCount

  function rotationCurveSize(self,time)
    !!{
    Return the number of array elements in the {\normalfont \ttfamily rotationCurve} property extractors.
    !!}
    implicit none
    integer         (c_size_t                          )                :: rotationCurveSize
    class           (nodePropertyExtractorRotationCurve), intent(inout) :: self
    double precision                                    , intent(in   ) :: time
    !$GLC attributes unused :: time

    rotationCurveSize=self%radiiCount
    return
  end function rotationCurveSize

  function rotationCurveExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily rotationCurve} property extractor.
    !!}
    use :: Galactic_Structure_Options          , only : componentTypeAll               , massTypeGalactic                  , massTypeStellar
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius, radiusTypeDiskHalfMassRadius      , radiusTypeDiskRadius                      , radiusTypeGalacticLightFraction, &
          &                                             radiusTypeGalacticMassFraction , radiusTypeRadius                  , radiusTypeSpheroidHalfMassRadius          , radiusTypeSpheroidRadius       , &
          &                                             radiustypestellarmassfraction  , radiusTypeNuclearStarClusterRadius, radiusTypeNuclearStarClusterHalfMassRadius, radiusTypeVirialRadius
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile , nodeComponentDisk                 , nodeComponentSpheroid                     , nodeComponentNSC               , &
          &                                             treeNode
    use :: Error                               , only : Error_Report
    use :: Mass_Distributions                  , only : massDistributionClass
    implicit none
    double precision                                    , dimension(:,:), allocatable :: rotationCurveExtract
    class           (nodePropertyExtractorRotationCurve), intent(inout) , target      :: self
    type            (treeNode                          ), intent(inout) , target      :: node
    double precision                                    , intent(in   )               :: time
    type            (multiCounter                      ), intent(inout) , optional    :: instance
    class           (nodeComponentDisk                 ), pointer                     :: disk
    class           (nodeComponentSpheroid             ), pointer                     :: spheroid
    class           (nodeComponentNSC                  ), pointer                     :: nuclearStarCluster
    class           (nodeComponentDarkMatterProfile    ), pointer                     :: darkMatterProfile
    class           (massDistributionClass             ), pointer                     :: massDistribution_
    integer                                                                           :: i
    double precision                                                                  :: radius                , radiusVirial
    !$GLC attributes unused :: time, instance

    allocate(rotationCurveExtract(self%radiiCount,self%elementCount_))
    radiusVirial                                         =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial       =  self%darkMatterHaloScale_%radiusVirial(node                    )
    if (self%                 diskIsNeeded) disk               =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid           =>                                        node%spheroid         ()
    if (self%   nuclearStarClusterIsNeeded) nuclearStarCluster =>                                        node%NSC              ()
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
       case   (radiusTypeDiskRadius                      %ID)
          radius=+radius*disk              %        radius()
       case   (radiusTypeSpheroidRadius                  %ID)
          radius=+radius*spheroid          %        radius()
       case   (radiusTypeNuclearStarClusterRadius        %ID)
          radius=+radius*nuclearStarCluster%        radius()
       case   (radiusTypeDiskHalfMassRadius              %ID)
          radius=+radius*disk              %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius          %ID)
          radius=+radius*spheroid          %halfMassRadius()
       case   (radiusTypeNuclearStarClusterHalfMassRadius%ID)
          radius=+radius*nuclearStarCluster%halfMassRadius()
       case   (radiusTypeGalacticMassFraction            %ID,  &
            &  radiusTypeGalacticLightFraction           %ID)
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
       massDistribution_                => node             %massDistribution(                                       &
               &                                                              componentType=self%radii(i)%component, &
               &                                                              massType     =self%radii(i)%mass       &
               &                                                             )
       rotationCurveExtract       (i,1) =  massDistribution_%rotationCurve   (                                       &
            &                                                                                             radius     &
            &                                                                )
       if (self%includeRadii)                                                                                        &
            & rotationCurveExtract(i,2) =                                                                 radius
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end do
    return
  end function rotationCurveExtract

  subroutine rotationCurveNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily rotationCurve} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorRotationCurve), intent(inout)                             :: self
    double precision                                    , intent(in   ), optional                   :: time
    type            (varying_string                    ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names       (1)="rotationCurve"
    if (self%includeRadii)                             &
         & names(2)="rotationCurveRadius"
    return
  end subroutine rotationCurveNames

  subroutine rotationCurveDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily rotationCurve} property.
    !!}
    implicit none
    class           (nodePropertyExtractorRotationCurve), intent(inout)                             :: self
    double precision                                    , intent(in   ), optional                   :: time
    type            (varying_string                    ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time
    
    allocate(descriptions(self%elementCount_))
    descriptions       (1)="Rotation curve at a given radius [km s⁻¹]."
    if (self%includeRadii)                                                                &
         & descriptions(2)="Radius at which rotation curve is output [Mpc]."
    return
  end subroutine rotationCurveDescriptions

  subroutine rotationCurveColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnitsInSI,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily rotationCurve} property.
    !!}
    implicit none
    class           (nodePropertyExtractorRotationCurve), intent(inout)                            :: self
    double precision                                    , intent(in   ), optional                  :: time
    type            (varying_string                    ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                    , intent(inout), dimension(:), allocatable :: values
    type            (varying_string                    ), intent(  out)                            :: valuesDescription
    double precision                                    , intent(  out)                            :: valuesUnitsInSI
    !$GLC attributes unused :: time

    allocate(descriptions(self%radiiCount))
    allocate(values      (              0))
    valuesDescription=var_str('')
    valuesUnitsInSI  =0.0d0
    descriptions     =self%radii%name
    return
  end subroutine rotationCurveColumnDescriptions

  function rotationCurveUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily rotationCurve} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                    , allocatable  , dimension(:) :: rotationCurveUnitsInSI
    class           (nodePropertyExtractorRotationCurve), intent(inout)               :: self
    double precision                                    , intent(in   ), optional     :: time
    !$GLC attributes unused :: time
    
    allocate(rotationCurveUnitsInSI(self%elementCount_))
    rotationCurveUnitsInSI       (1)=kilo
    if (self%includeRadii)                      &
         & rotationCurveUnitsInSI(2)=megaParsec
    return
  end function rotationCurveUnitsInSI


!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Contains a module which implements a property extractor class for the projected density at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale   , darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier
  use :: Galactic_Structure                  , only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorProjectedDensity">
   <description>
    A property extractor class for the projected density at a set of radii. The radii and types of projected density to output
    is specified by the {\normalfont \ttfamily radiusSpecifiers} parameter. This parameter's value can contain multiple
    entries, each of which should be a valid
    \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#sec.radiusSpecifiers}{radius
    specifier}.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorProjectedDensity
     !!{
     A property extractor class for the projected density at a set of radii.
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
     final     ::                       projectedDensityDestructor
     procedure :: columnDescriptions => projectedDensityColumnDescriptions
     procedure :: size               => projectedDensitySize
     procedure :: elementCount       => projectedDensityElementCount
     procedure :: extract            => projectedDensityExtract
     procedure :: names              => projectedDensityNames
     procedure :: descriptions       => projectedDensityDescriptions
     procedure :: unitsInSI          => projectedDensityUnitsInSI
  end type nodePropertyExtractorProjectedDensity

  interface nodePropertyExtractorProjectedDensity
     !!{
     Constructors for the ``projectedDensity'' output analysis class.
     !!}
     module procedure projectedDensityConstructorParameters
     module procedure projectedDensityConstructorInternal
  end interface nodePropertyExtractorProjectedDensity

  ! Module-scope variables used in integrands.
  double precision :: projectedDensityRadius
  !$omp threadprivate(projectedDensityRadius)

contains

  function projectedDensityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily projectedDensity} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorProjectedDensity)                              :: self
    type   (inputParameters                      ), intent(inout)               :: parameters
    type   (varying_string                       ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass             ), pointer                     :: darkMatterHaloScale_
    class  (galacticStructureClass               ), pointer                     :: galacticStructure_
    logical                                                                     :: includeRadii

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter>
      <name>radiusSpecifiers</name>
      <description>A list of radius specifiers at which to output the projected density profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the radii at which projected density data are output should also be included in the output file.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="galacticStructure"   name="galacticStructure_"   source="parameters"/>
    !!]
    self=nodePropertyExtractorProjectedDensity(radiusSpecifiers,includeRadii,darkMatterHaloScale_,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="galacticStructure_"  />
    !!]
    return
  end function projectedDensityConstructorParameters

  function projectedDensityConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily projectedDensity} property extractor class.
    !!}
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorProjectedDensity)                              :: self
    type   (varying_string                       ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass             ), intent(in   ), target       :: darkMatterHaloScale_
    class  (galacticStructureClass               ), intent(in   ), target       :: galacticStructure_
    logical                                       , intent(in   )               :: includeRadii
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
  end function projectedDensityConstructorInternal

  subroutine projectedDensityDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily projectedDensity} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorProjectedDensity), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%galacticStructure_"  />
    !!]
    return
  end subroutine projectedDensityDestructor

  integer function projectedDensityElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily projectedDensity} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    projectedDensityElementCount=self%elementCount_
    return
  end function projectedDensityElementCount

  function projectedDensitySize(self,time)
    !!{
    Return the number of array alements in the {\normalfont \ttfamily projectedDensity} property extractors.
    !!}
    implicit none
    integer         (c_size_t                             )                :: projectedDensitySize
    class           (nodePropertyExtractorProjectedDensity), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    projectedDensitySize=self%radiiCount
    return
  end function projectedDensitySize

  function projectedDensityExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily projectedDensity} property extractor.
    !!}
    use :: Galactic_Structure_Options          , only : componentTypeAll               , massTypeGalactic            , massTypeStellar
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius, radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius            , radiusTypeGalacticLightFraction, &
          &                                             radiusTypeGalacticMassFraction , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius, radiusTypeSpheroidRadius       , &
          &                                             radiusTypeStellarMassFraction  , radiusTypeVirialRadius
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile , nodeComponentDisk           , nodeComponentSpheroid           , treeNode
    use :: Numerical_Integration               , only : integrator, GSL_Integ_Gauss15
    implicit none
    double precision                                       , dimension(:,:), allocatable :: projectedDensityExtract
    class           (nodePropertyExtractorProjectedDensity), intent(inout) , target      :: self
    type            (treeNode                             ), intent(inout) , target      :: node
    double precision                                       , intent(in   )               :: time
    type            (multiCounter                         ), intent(inout) , optional    :: instance
    class           (nodeComponentDisk                    ), pointer                     :: disk
    class           (nodeComponentSpheroid                ), pointer                     :: spheroid
    class           (nodeComponentDarkMatterProfile       ), pointer                     :: darkMatterProfile
    double precision                                       , parameter                   :: epsilonSingularity     =1.0d-3
    type            (integrator                           )                              :: integrator_
    integer                                                                              :: i
    double precision                                                                     :: radiusVirial                  , radiusOuter, &
         &                                                                                  radiusSingularity
    !$GLC attributes unused :: time, instance

    allocate(projectedDensityExtract(self%radiiCount,self%elementCount_))
    radiusVirial                                         =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial      =  self%darkMatterHaloScale_%radiusVirial(node                    )
    if (self%                 diskIsNeeded) disk              =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid          =>                                        node%spheroid         ()
    if (self%darkMatterScaleRadiusIsNeeded) darkMatterProfile =>                                        node%darkMatterProfile()
    integrator_=integrator(projectedDensityIntegrand,toleranceRelative=1.0d-3,hasSingularities=.true.,integrationRule=GSL_Integ_Gauss15)
    do i=1,self%radiiCount
       projectedDensityRadius=self%radii(i)%value
       select case (self%radii(i)%type%ID)
       case   (radiusTypeRadius                %ID)
          ! Nothing to do.
       case   (radiusTypeVirialRadius          %ID)
          projectedDensityRadius=+projectedDensityRadius*radiusVirial
       case   (radiusTypeDarkMatterScaleRadius %ID)
          projectedDensityRadius=+projectedDensityRadius*darkMatterProfile%         scale()
       case   (radiusTypeDiskRadius            %ID)
          projectedDensityRadius=+projectedDensityRadius*disk             %        radius()
       case   (radiusTypeSpheroidRadius        %ID)
          projectedDensityRadius=+projectedDensityRadius*spheroid         %        radius()
       case   (radiusTypeDiskHalfMassRadius    %ID)
          projectedDensityRadius=+projectedDensityRadius*disk             %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius%ID)
          projectedDensityRadius=+projectedDensityRadius*spheroid         %halfMassRadius()
       case   (radiusTypeGalacticMassFraction  %ID,  &
            &  radiusTypeGalacticLightFraction %ID)
          projectedDensityRadius=+projectedDensityRadius           &
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
          projectedDensityRadius=+projectedDensityRadius           &
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
       radiusOuter                        =self       %darkMatterHaloScale_%radiusVirial(node                              )
       ! Cut out a small region round the coordinate singularity at the inner radius. This region will be integrated analytically
       ! assuming a constant density over this region. The region outside of this cut-out will be integrated numerically.
       radiusSingularity=min(projectedDensityRadius*(1.0d0+epsilonSingularity),radiusOuter)
       !! Analytic integral within the cut-out.
       projectedDensityExtract       (i,1)=+2.0d0                                                                  &
            &                              *sqrt(                                                                  &
            &                                    +radiusSingularity     **2                                        &
            &                                    -projectedDensityRadius**2                                        &
            &                                   )                                                                  &
            &                              *self%galacticStructure_%density(                                       &
            &                                                               node                                 , &
            &                                                               [                                      &
            &                                                                projectedDensityRadius              , &
            &                                                                0.0d0                               , &
            &                                                                0.0d0                                 &
            &                                                               ]                                    , &
            &                                                               componentType=self%radii(i)%component, &
            &                                                               massType     =self%radii(i)%mass       &
            &                                                              )
       !! Numerical integral outside of the cut-out.
       if (radiusSingularity > radiusOuter)                                                                &       
            & projectedDensityExtract(i,1)=+projectedDensityExtract(i,1)                                   &
            &                              +integrator_%integrate(log(radiusSingularity),log(radiusOuter))
       if (self%includeRadii)                                                                              &
            & projectedDensityExtract(i,2)=projectedDensityRadius
    end do
    return

  contains

    double precision function projectedDensityIntegrand(radiusLogarithmic)
      !!{
      Integrand function used for computing projected densities.
      !!}
      implicit none
      double precision, intent(in   ) :: radiusLogarithmic
      double precision                :: radius

      radius=exp(radiusLogarithmic)
      if (radius <= projectedDensityRadius) then
         projectedDensityIntegrand=+0.0d0
      else
         projectedDensityIntegrand=+2.0d0                                                                  &
              &                    *radius                      **2                                        &
              &                    /sqrt(                                                                  &
              &                          +radius                **2                                        &
              &                          -projectedDensityRadius**2                                        &
              &                    )                                                                       &
              &                    *self%galacticStructure_%density(                                       &
              &                                                     node                                 , &
              &                                                     [                                      &
              &                                                      radius                              , &
              &                                                      0.0d0                               , &
              &                                                      0.0d0                                 &
              &                                                     ]                                    , &
              &                                                     componentType=self%radii(i)%component, &
              &                                                     massType     =self%radii(i)%mass       &
              &                                                    )
      end if
      return
    end function projectedDensityIntegrand

  end function projectedDensityExtract

  subroutine projectedDensityNames(self,names,time)
    !!{
    Return the names of the {\normalfont \ttfamily projectedDensity} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names       (1)="projectedDensity"
    if (self%includeRadii)                   &
         & names(2)="projectedDensityRadius"
    return
  end subroutine projectedDensityNames

  subroutine projectedDensityDescriptions(self,descriptions,time)
    !!{
    Return descriptions of the {\normalfont \ttfamily projectedDensity} property.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions       (1)="Projected density at a given radius [M☉/Mpc⁻²]."
    if (self%includeRadii)                                                      &
         & descriptions(2)="Radius at which projected density is output [Mpc]."
    return
  end subroutine projectedDensityDescriptions

  subroutine projectedDensityColumnDescriptions(self,descriptions,time)
    !!{
    Return column descriptions of the {\normalfont \ttfamily projectedDensity} property.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%radiiCount))
    descriptions=self%radii%name
    return
  end subroutine projectedDensityColumnDescriptions

  function projectedDensityUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily projectedDensity} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                       , allocatable  , dimension(:) :: projectedDensityUnitsInSI
    class           (nodePropertyExtractorProjectedDensity), intent(inout)               :: self
    double precision                                       , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(projectedDensityUnitsInSI(self%elementCount_))
    projectedDensityUnitsInSI       (1)=massSolar/megaParsec**2
    if (self%includeRadii)                                       &
         & projectedDensityUnitsInSI(2)=          megaParsec
    return
  end function projectedDensityUnitsInSI


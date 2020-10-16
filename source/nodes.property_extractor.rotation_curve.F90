!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Contains a module which implements a property extractor class for the rotation curve at a set of radii.
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale, darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier

  !# <nodePropertyExtractor name="nodePropertyExtractorRotationCurve">
  !#  <description>
  !#   A property extractor class for the rotation curve at a set of radii. The radii and types of rotation curve to output
  !#   are specified by the {\normalfont \ttfamily radiusSpecifiers} parameter. This parameter's value can contain multiple
  !#   entries, each of which should be a valid
  !#   \href{https://github.com/galacticusorg/galacticus/releases/download/masterRelease/Galacticus_Physics.pdf\#sec.radiusSpecifiers}{radius
  !#   specifier}.
  !#  </description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorRotationCurve
     !% A property extractor class for the rotation curve at a set of radii.
     private
     class  (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_
     integer                                                      :: radiiCount                   , elementCount_       , &
          &                                                          step
     logical                                                      :: includeRadii
     type   (varying_string          ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusSpecifier         ), allocatable, dimension(:) :: radii
     logical                                                      :: darkMatterScaleRadiusIsNeeded, diskIsNeeded        , &
          &                                                          spheroidIsNeeded             , virialRadiusIsNeeded
   contains
     final     ::                 rotationCurveDestructor
     procedure :: elementCount => rotationCurveElementCount
     procedure :: extract      => rotationCurveExtract
     procedure :: names        => rotationCurveNames
     procedure :: descriptions => rotationCurveDescriptions
     procedure :: unitsInSI    => rotationCurveUnitsInSI
     procedure :: type         => rotationCurveType
  end type nodePropertyExtractorRotationCurve

  interface nodePropertyExtractorRotationCurve
     !% Constructors for the ``rotationCurve'' output analysis class.
     module procedure rotationCurveConstructorParameters
     module procedure rotationCurveConstructorInternal
  end interface nodePropertyExtractorRotationCurve

contains

  function rotationCurveConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily rotationCurve} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorRotationCurve)                              :: self
    type   (inputParameters                   ), intent(inout)               :: parameters
    type   (varying_string                    ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass          ), pointer                     :: darkMatterHaloScale_
    logical                                                                  :: includeRadii

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !# <inputParameter>
    !#   <name>radiusSpecifiers</name>
    !#   <description>A list of radius specifiers at which to output the rotation curve.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>includeRadii</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not the radii at which rotation curve data are output should also be included in the output file.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=nodePropertyExtractorRotationCurve(radiusSpecifiers,includeRadii,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_" />
    return
  end function rotationCurveConstructorParameters

  function rotationCurveConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily rotationCurve} property extractor class.
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorRotationCurve)                              :: self
    type   (varying_string                    ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass          ), intent(in   ), target       :: darkMatterHaloScale_
    logical                                    , intent(in   )               :: includeRadii
    !# <constructorAssign variables="radiusSpecifiers, includeRadii, *darkMatterHaloScale_"/>

    if (includeRadii) then
       self%step=2
    else
       self%step=1
    end if
    self%radiiCount   =size(radiusSpecifiers)
    self%elementCount_=self%step*self%radiiCount
    call Galactic_Structure_Radii_Definition_Decode(                                    &
         &                                          radiusSpecifiers                  , &
         &                                          self%radii                        , &
         &                                          self%diskIsNeeded                 , &
         &                                          self%spheroidIsNeeded             , &
         &                                          self%virialRadiusIsNeeded         , &
         &                                          self%darkMatterScaleRadiusIsNeeded  &
         &                                         )
    return
  end function rotationCurveConstructorInternal

  subroutine rotationCurveDestructor(self)
    !% Destructor for the {\normalfont \ttfamily rotationCurve} property extractor class.
    implicit none
    type(nodePropertyExtractorRotationCurve), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine rotationCurveDestructor

  integer function rotationCurveElementCount(self,time)
    !% Return the number of elements in the {\normalfont \ttfamily rotationCurve} property extractors.
    implicit none
    class           (nodePropertyExtractorRotationCurve), intent(inout) :: self
    double precision                                    , intent(in   ) :: time
    !$GLC attributes unused :: time

    rotationCurveElementCount=self%elementCount_
    return
  end function rotationCurveElementCount

  function rotationCurveExtract(self,node,time,instance)
    !% Implement a {\normalfont \ttfamily rotationCurve} property extractor.
    use :: Galactic_Structure_Enclosed_Masses  , only : Galactic_Structure_Radius_Enclosing_Mass
    use :: Galactic_Structure_Options          , only : componentTypeAll                        , massTypeGalactic
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius         , radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius            , radiusTypeGalacticLightFraction, &
          &                                             radiusTypeGalacticMassFraction          , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius, radiusTypeSpheroidRadius       , &
          &                                             radiusTypeVirialRadius
    use :: Galactic_Structure_Rotation_Curves  , only : Galactic_Structure_Rotation_Curve
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile          , nodeComponentDisk           , nodeComponentSpheroid           , treeNode
    implicit none
    double precision                                     , dimension(:) , allocatable :: rotationCurveExtract
    class           (nodePropertyExtractorRotationCurve), intent(inout), target      :: self
    type            (treeNode                           ), intent(inout), target      :: node
    double precision                                     , intent(in   )              :: time
    type            (multiCounter                       ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                  ), pointer                    :: disk
    class           (nodeComponentSpheroid              ), pointer                    :: spheroid
    class           (nodeComponentDarkMatterProfile     ), pointer                    :: darkMatterProfile
    integer                                                                           :: i
    double precision                                                                  :: radius                , radiusVirial
    !$GLC attributes unused :: time, instance

    allocate(rotationCurveExtract(self%elementCount_))
    radiusVirial                                         =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial      =  self%darkMatterHaloScale_%virialRadius(node                    )
    if (self%                 diskIsNeeded) disk              =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid          =>                                        node%spheroid         ()
    if (self%darkMatterScaleRadiusIsNeeded) darkMatterProfile =>                                        node%darkMatterProfile()
    do i=1,self%radiiCount
       radius=self%radii(i)%value
       select case (self%radii(i)%type)
       case   (radiusTypeRadius                )
          ! Nothing to do.
       case   (radiusTypeVirialRadius          )
          radius=+radius*radiusVirial
       case   (radiusTypeDarkMatterScaleRadius )
          radius=+radius*darkMatterProfile%         scale()
       case   (radiusTypeDiskRadius            )
          radius=+radius*disk             %        radius()
       case   (radiusTypeSpheroidRadius        )
          radius=+radius*spheroid         %        radius()
       case   (radiusTypeDiskHalfMassRadius    )
          radius=+radius*disk             %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius)
          radius=+radius*spheroid         %halfMassRadius()
       case   (radiusTypeGalacticMassFraction ,  &
            &  radiusTypeGalacticLightFraction )
          radius=+radius                                           &
               & *Galactic_Structure_Radius_Enclosing_Mass         &
               &  (                                                &
               &   node                                         ,  &
               &   fractionalMass=self%radii(i)%fraction        ,  &
               &   massType      =              massTypeGalactic,  &
               &   componentType =              componentTypeAll,  &
               &   weightBy      =self%radii(i)%weightBy        ,  &
               &   weightIndex   =self%radii(i)%weightByIndex      &
               &  )
       end select
       rotationCurveExtract       ((i-1)*self%step+1)=Galactic_Structure_Rotation_Curve(                                       &
               &                                                                        node                                 , &
               &                                                                        radius                               , &
               &                                                                        componentType=self%radii(i)%component, &
               &                                                                        massType     =self%radii(i)%mass       &
               &                                                                       )
       if (self%includeRadii)                                                                                                  &
            & rotationCurveExtract((i-1)*self%step+2)=                                  radius
    end do
    return
  end function rotationCurveExtract

  function rotationCurveNames(self,time)
    !% Return the names of the {\normalfont \ttfamily rotationCurve} properties.
    implicit none
    type            (varying_string                    ), dimension(:) , allocatable :: rotationCurveNames
    class           (nodePropertyExtractorRotationCurve), intent(inout)              :: self
    double precision                                    , intent(in   )              :: time
    integer                                                                          :: i
    !$GLC attributes unused :: time

    allocate(rotationCurveNames(self%elementCount_))
    do i=1,size(self%radii)
       rotationCurveNames       ((i-1)*self%step+1)="rotationCurve:"      //char(self%radii(i)%name)
       if (self%includeRadii)                                                                          &
            & rotationCurveNames((i-1)*self%step+2)="rotationCurveRadius:"//char(self%radii(i)%name)
    end do
    return
  end function rotationCurveNames

  function rotationCurveDescriptions(self,time)
    !% Return descriptions of the {\normalfont \ttfamily rotationCurve} property.
    implicit none
    type            (varying_string                    ), dimension(:) , allocatable :: rotationCurveDescriptions
    class           (nodePropertyExtractorRotationCurve), intent(inout)              :: self
    double precision                                    , intent(in   )              :: time
    integer                                                                          :: i
    !$GLC attributes unused :: time

    allocate(rotationCurveDescriptions(self%elementCount_))
    do i=1,size(self%radii)
       rotationCurveDescriptions       ((i-1)*self%step+1)="Rotation curve at a given radius [km s⁻¹]."
       if (self%includeRadii)                                                                                &
            & rotationCurveDescriptions((i-1)*self%step+2)="Radius at which rotation curve is output [Mpc]."
    end do
    return
  end function rotationCurveDescriptions

  function rotationCurveUnitsInSI(self,time)
    !% Return the units of the {\normalfont \ttfamily rotationCurve} properties in the SI system.
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                    , allocatable  , dimension(:) :: rotationCurveUnitsInSI
    class           (nodePropertyExtractorRotationCurve), intent(inout)               :: self
    double precision                                    , intent(in   )               :: time
    integer                                                                           :: i
    !$GLC attributes unused :: time

    allocate(rotationCurveUnitsInSI(self%elementCount_))
    do i=1,size(self%radii)
       rotationCurveUnitsInSI       ((i-1)*self%step+1)=kilo
       if (self%includeRadii)                                      &
            & rotationCurveUnitsInSI((i-1)*self%step+2)=megaParsec
    end do
    return
  end function rotationCurveUnitsInSI

  integer function rotationCurveType(self)
    !% Return the type of the {\normalfont \ttfamily rotationCurve} properties.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorRotationCurve), intent(inout) :: self
    !$GLC attributes unused :: self

    rotationCurveType=outputAnalysisPropertyTypeLinear
    return
  end function rotationCurveType

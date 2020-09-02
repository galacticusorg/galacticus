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

  !% Contains a module which implements a property extractor class for the density at a set of radii.
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale, darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier

  !# <nodePropertyExtractor name="nodePropertyExtractorDensityProfile">
  !#  <description>A property extractor class for the density at a set of radii.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorDensityProfile
     !% A property extractor class for the density at a set fo radii.
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
     final     ::                 densityProfileDestructor
     procedure :: elementCount => densityProfileElementCount
     procedure :: extract      => densityProfileExtract
     procedure :: names        => densityProfileNames
     procedure :: descriptions => densityProfileDescriptions
     procedure :: unitsInSI    => densityProfileUnitsInSI
     procedure :: type         => densityProfileType
  end type nodePropertyExtractorDensityProfile

  interface nodePropertyExtractorDensityProfile
     !% Constructors for the ``densityProfile'' output analysis class.
     module procedure densityProfileConstructorParameters
     module procedure densityProfileConstructorInternal
  end interface nodePropertyExtractorDensityProfile

contains

  function densityProfileConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily densityProfile} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorDensityProfile)                              :: self
    type   (inputParameters                    ), intent(inout)               :: parameters
    type   (varying_string                     ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass           ), pointer                     :: darkMatterHaloScale_
    logical                                                                   :: includeRadii

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !# <inputParameter>
    !#   <name>radiusSpecifiers</name>
    !#   <description>A list of radius specifiers at which to output the density profile.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>includeRadii</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not the radii at which density data are output should also be included in the output file.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=nodePropertyExtractorDensityProfile(radiusSpecifiers,includeRadii,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_" />
    return
  end function densityProfileConstructorParameters

  function densityProfileConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily densityProfile} property extractor class.
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorDensityProfile)                              :: self
    type   (varying_string                     ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass           ), intent(in   ), target       :: darkMatterHaloScale_
    logical                                     , intent(in   )               :: includeRadii
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
  end function densityProfileConstructorInternal

  subroutine densityProfileDestructor(self)
    !% Destructor for the {\normalfont \ttfamily densityProfile} property extractor class.
    implicit none
    type(nodePropertyExtractorDensityProfile), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine densityProfileDestructor

  integer function densityProfileElementCount(self,time)
    !% Return the number of elements in the {\normalfont \ttfamily densityProfile} property extractors.
    implicit none
    class           (nodePropertyExtractorDensityProfile), intent(inout) :: self
    double precision                                     , intent(in   ) :: time
    !$GLC attributes unused :: time

    densityProfileElementCount=self%elementCount_
    return
  end function densityProfileElementCount

  function densityProfileExtract(self,node,time,instance)
    !% Implement a {\normalfont \ttfamily densityProfile} property extractor.
    use :: Galactic_Structure_Densities        , only : Galactic_Structure_Density
    use :: Galactic_Structure_Enclosed_Masses  , only : Galactic_Structure_Radius_Enclosing_Mass
    use :: Galactic_Structure_Options          , only : componentTypeAll                        , massTypeGalactic
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius         , radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius            , radiusTypeGalacticLightFraction, &
          &                                             radiusTypeGalacticMassFraction          , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius, radiusTypeSpheroidRadius       , &
          &                                             radiusTypeVirialRadius
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile          , nodeComponentDisk           , nodeComponentSpheroid           , treeNode
    implicit none
    double precision                                     , dimension(:) , allocatable :: densityProfileExtract
    class           (nodePropertyExtractorDensityProfile), intent(inout), target      :: self
    type            (treeNode                           ), intent(inout), target      :: node
    double precision                                     , intent(in   )              :: time
    type            (multiCounter                       ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                  ), pointer                    :: disk
    class           (nodeComponentSpheroid              ), pointer                    :: spheroid
    class           (nodeComponentDarkMatterProfile     ), pointer                    :: darkMatterProfile
    integer                                                                           :: i
    double precision                                                                  :: radius                , radiusVirial
    !$GLC attributes unused :: time, instance

    allocate(densityProfileExtract(self%elementCount_))
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
       densityProfileExtract       ((i-1)*self%step+1)=Galactic_Structure_Density(                                       &
            &                                                                     node                                 , &
            &                                                                     [                                      &
            &                                                                      radius                              , &
            &                                                                      0.0d0                               , &
            &                                                                      0.0d0                                 &
            &                                                                     ]                                    , &
            &                                                                     componentType=self%radii(i)%component, &
            &                                                                     massType     =self%radii(i)%mass       &
            &                                                                    )
       if (self%includeRadii)                                                                                            &
            & densityProfileExtract((i-1)*self%step+2)=                            radius
    end do
    return
  end function densityProfileExtract

  function densityProfileNames(self,time)
    !% Return the names of the {\normalfont \ttfamily densityProfile} properties.
    implicit none
    type            (varying_string                     ), dimension(:) , allocatable :: densityProfileNames
    class           (nodePropertyExtractorDensityProfile), intent(inout)              :: self
    double precision                                     , intent(in   )              :: time
    integer                                                                           :: i
    !$GLC attributes unused :: time

    allocate(densityProfileNames(self%elementCount_))
    do i=1,size(self%radii)
       densityProfileNames       ((i-1)*self%step+1)="densityProfile:"      //char(self%radii(i)%name)
       if (self%includeRadii)                                                                          &
            & densityProfileNames((i-1)*self%step+2)="densityProfileRadius:"//char(self%radii(i)%name)
    end do
    return
  end function densityProfileNames

  function densityProfileDescriptions(self,time)
    !% Return descriptions of the {\normalfont \ttfamily densityProfile} property.
    implicit none
    type            (varying_string                     ), dimension(:) , allocatable :: densityProfileDescriptions
    class           (nodePropertyExtractorDensityProfile), intent(inout)              :: self
    double precision                                     , intent(in   )              :: time
    integer                                                                           :: i
    !$GLC attributes unused :: time

    allocate(densityProfileDescriptions(self%elementCount_))
    do i=1,size(self%radii)
       densityProfileDescriptions       ((i-1)*self%step+1)="Density at a given radius [M☉/Mpc⁻³]."
       if (self%includeRadii)                                                                          &
            & densityProfileDescriptions((i-1)*self%step+2)="Radius at which density is output [Mpc]."
    end do
    return
  end function densityProfileDescriptions

  function densityProfileUnitsInSI(self,time)
    !% Return the units of the {\normalfont \ttfamily densityProfile} properties in the SI system.
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                     , allocatable  , dimension(:) :: densityProfileUnitsInSI
    class           (nodePropertyExtractorDensityProfile), intent(inout)               :: self
    double precision                                     , intent(in   )               :: time
    integer                                                                            :: i
    !$GLC attributes unused :: time

    allocate(densityProfileUnitsInSI(self%elementCount_))
    do i=1,size(self%radii)
       densityProfileUnitsInSI       ((i-1)*self%step+1)=massSolar/megaParsec**3
       if (self%includeRadii)                                                       &
            & densityProfileUnitsInSI((i-1)*self%step+2)=          megaParsec
    end do
    return
  end function densityProfileUnitsInSI

  integer function densityProfileType(self)
    !% Return the type of the {\normalfont \ttfamily densityProfile} properties.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorDensityProfile), intent(inout) :: self
    !$GLC attributes unused :: self

    densityProfileType=outputAnalysisPropertyTypeLinear
    return
  end function densityProfileType

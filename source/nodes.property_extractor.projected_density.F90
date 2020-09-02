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

  !% Contains a module which implements a property extractor class for the projected density at a set of radii.
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale, darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier

  !# <nodePropertyExtractor name="nodePropertyExtractorProjectedDensity">
  !#  <description>A property extractor class for the projected density at a set of radii.</description>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorProjectedDensity
     !% A property extractor class for the projected density at a set of radii.
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
     final     ::                 projectedDensityDestructor
     procedure :: elementCount => projectedDensityElementCount
     procedure :: extract      => projectedDensityExtract
     procedure :: names        => projectedDensityNames
     procedure :: descriptions => projectedDensityDescriptions
     procedure :: unitsInSI    => projectedDensityUnitsInSI
     procedure :: type         => projectedDensityType
  end type nodePropertyExtractorProjectedDensity

  interface nodePropertyExtractorProjectedDensity
     !% Constructors for the ``projectedDensity'' output analysis class.
     module procedure projectedDensityConstructorParameters
     module procedure projectedDensityConstructorInternal
  end interface nodePropertyExtractorProjectedDensity

  ! Module-scope variables used in integrands.
  double precision :: projectedDensityRadius
  !$omp threadprivate(projectedDensityRadius)

contains

  function projectedDensityConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily projectedDensity} property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorProjectedDensity)                              :: self
    type   (inputParameters                      ), intent(inout)               :: parameters
    type   (varying_string                       ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass             ), pointer                     :: darkMatterHaloScale_
    logical                                                                     :: includeRadii

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !# <inputParameter>
    !#   <name>radiusSpecifiers</name>
    !#   <description>A list of radius specifiers at which to output the projected density profile.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>includeRadii</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether or not the radii at which projected density data are output should also be included in the output file.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=nodePropertyExtractorProjectedDensity(radiusSpecifiers,includeRadii,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_" />
    return
  end function projectedDensityConstructorParameters

  function projectedDensityConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily projectedDensity} property extractor class.
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorProjectedDensity)                              :: self
    type   (varying_string                       ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass             ), intent(in   ), target       :: darkMatterHaloScale_
    logical                                       , intent(in   )               :: includeRadii
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
  end function projectedDensityConstructorInternal

  subroutine projectedDensityDestructor(self)
    !% Destructor for the {\normalfont \ttfamily projectedDensity} property extractor class.
    implicit none
    type(nodePropertyExtractorProjectedDensity), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    return
  end subroutine projectedDensityDestructor

  integer function projectedDensityElementCount(self,time)
    !% Return the number of elements in the {\normalfont \ttfamily projectedDensity} property extractors.
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    projectedDensityElementCount=self%elementCount_
    return
  end function projectedDensityElementCount

  function projectedDensityExtract(self,node,time,instance)
    !% Implement a {\normalfont \ttfamily projectedDensity} property extractor.
    use :: Galactic_Structure_Densities        , only : Galactic_Structure_Density
    use :: Galactic_Structure_Enclosed_Masses  , only : Galactic_Structure_Radius_Enclosing_Mass
    use :: Galactic_Structure_Options          , only : componentTypeAll                        , massTypeGalactic
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius         , radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius            , radiusTypeGalacticLightFraction, &
          &                                             radiusTypeGalacticMassFraction          , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius, radiusTypeSpheroidRadius       , &
          &                                             radiusTypeVirialRadius
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile          , nodeComponentDisk           , nodeComponentSpheroid           , treeNode
    use :: Numerical_Integration               , only : integrator
    implicit none
    double precision                                       , dimension(:) , allocatable :: projectedDensityExtract
    class           (nodePropertyExtractorProjectedDensity), intent(inout), target      :: self
    type            (treeNode                             ), intent(inout), target      :: node
    double precision                                       , intent(in   )              :: time
    type            (multiCounter                         ), intent(inout), optional    :: instance
    class           (nodeComponentDisk                    ), pointer                    :: disk
    class           (nodeComponentSpheroid                ), pointer                    :: spheroid
    class           (nodeComponentDarkMatterProfile       ), pointer                    :: darkMatterProfile
    type            (integrator                           )                             :: integrator_
    integer                                                                             :: i
    double precision                                                                    :: radiusVirial           , radiusOuter
    !$GLC attributes unused :: time, instance

    allocate(projectedDensityExtract(self%elementCount_))
    radiusVirial                                         =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial      =  self%darkMatterHaloScale_%virialRadius(node                    )
    if (self%                 diskIsNeeded) disk              =>                                        node%disk             ()
    if (self%             spheroidIsNeeded) spheroid          =>                                        node%spheroid         ()
    if (self%darkMatterScaleRadiusIsNeeded) darkMatterProfile =>                                        node%darkMatterProfile()
    integrator_=integrator(projectedDensityIntegrand,toleranceRelative=1.0d-3)
    do i=1,self%radiiCount
       projectedDensityRadius=self%radii(i)%value
       select case (self%radii(i)%type)
       case   (radiusTypeRadius                )
          ! Nothing to do.
       case   (radiusTypeVirialRadius          )
          projectedDensityRadius=+projectedDensityRadius*radiusVirial
       case   (radiusTypeDarkMatterScaleRadius )
          projectedDensityRadius=+projectedDensityRadius*darkMatterProfile%         scale()
       case   (radiusTypeDiskRadius            )
          projectedDensityRadius=+projectedDensityRadius*disk             %        radius()
       case   (radiusTypeSpheroidRadius        )
          projectedDensityRadius=+projectedDensityRadius*spheroid         %        radius()
       case   (radiusTypeDiskHalfMassRadius    )
          projectedDensityRadius=+projectedDensityRadius*disk             %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius)
          projectedDensityRadius=+projectedDensityRadius*spheroid         %halfMassRadius()
       case   (radiusTypeGalacticMassFraction ,  &
            &  radiusTypeGalacticLightFraction )
          projectedDensityRadius=+projectedDensityRadius           &
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
       radiusOuter                                      =self       %darkMatterHaloScale_%virialRadius(node                              )       
       projectedDensityExtract       ((i-1)*self%step+1)=integrator_                     %integrate   (projectedDensityRadius,radiusOuter)
       if (self%includeRadii)                                                                                                              &
            & projectedDensityExtract((i-1)*self%step+2)=                                              projectedDensityRadius
    end do
    return

  contains

    double precision function projectedDensityIntegrand(radius)
      !% Integrand function used for computing projected densities.
      use :: Galactic_Structure_Densities, only : Galactic_Structure_Density
      implicit none
      double precision, intent(in   ) :: radius

      if (radius <= projectedDensityRadius) then
         projectedDensityIntegrand=+0.0d0
      else
         projectedDensityIntegrand=+2.0d0                                                             &
              &                    *radius                                                            &
              &                    /sqrt(                                                             &
              &                          +radius                **2                                   &
              &                          -projectedDensityRadius**2                                   &
              &                    )                                                                  &
              &                    *Galactic_Structure_Density(                                       &
              &                                                node                                 , &
              &                                                [                                      &
              &                                                 radius                              , &
              &                                                 0.0d0                               , &
              &                                                 0.0d0                                 &
              &                                                ]                                    , &
              &                                                componentType=self%radii(i)%component, &
              &                                                massType     =self%radii(i)%mass       &
              &                                               )
      end if
      return
    end function projectedDensityIntegrand

  end function projectedDensityExtract

  function projectedDensityNames(self,time)
    !% Return the names of the {\normalfont \ttfamily projectedDensity} properties.
    implicit none
    type            (varying_string                     ), dimension(:) , allocatable :: projectedDensityNames
    class           (nodePropertyExtractorProjectedDensity), intent(inout)              :: self
    double precision                                     , intent(in   )              :: time
    integer                                                                           :: i
    !$GLC attributes unused :: time

    allocate(projectedDensityNames(self%elementCount_))
    do i=1,size(self%radii)
       projectedDensityNames       ((i-1)*self%step+1)="projectedDensity:"      //char(self%radii(i)%name)
       if (self%includeRadii)                                                                          &
            & projectedDensityNames((i-1)*self%step+2)="projectedDensityRadius:"//char(self%radii(i)%name)
    end do
    return
  end function projectedDensityNames

  function projectedDensityDescriptions(self,time)
    !% Return descriptions of the {\normalfont \ttfamily projectedDensity} property.
    implicit none
    type            (varying_string                       ), dimension(:) , allocatable :: projectedDensityDescriptions
    class           (nodePropertyExtractorProjectedDensity), intent(inout)              :: self
    double precision                                       , intent(in   )              :: time
    integer                                                                             :: i
    !$GLC attributes unused :: time

    allocate(projectedDensityDescriptions(self%elementCount_))
    do i=1,size(self%radii)
       projectedDensityDescriptions       ((i-1)*self%step+1)="Projected density at a given radius [M☉/Mpc⁻²]."
       if (self%includeRadii)                                                                                      &
            & projectedDensityDescriptions((i-1)*self%step+2)="Radius at which projected density is output [Mpc]."
    end do
    return
  end function projectedDensityDescriptions

  function projectedDensityUnitsInSI(self,time)
    !% Return the units of the {\normalfont \ttfamily projectedDensity} properties in the SI system.
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                       , allocatable  , dimension(:) :: projectedDensityUnitsInSI
    class           (nodePropertyExtractorProjectedDensity), intent(inout)               :: self
    double precision                                       , intent(in   )               :: time
    integer                                                                              :: i
    !$GLC attributes unused :: time

    allocate(projectedDensityUnitsInSI(self%elementCount_))
    do i=1,size(self%radii)
       projectedDensityUnitsInSI       ((i-1)*self%step+1)=massSolar/megaParsec**2
       if (self%includeRadii)                                                       &
            & projectedDensityUnitsInSI((i-1)*self%step+2)=          megaParsec
    end do
    return
  end function projectedDensityUnitsInSI

  integer function projectedDensityType(self)
    !% Return the type of the {\normalfont \ttfamily projectedDensity} properties.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorProjectedDensity), intent(inout) :: self
    !$GLC attributes unused :: self

    projectedDensityType=outputAnalysisPropertyTypeLinear
    return
  end function projectedDensityType

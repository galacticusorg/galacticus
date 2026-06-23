!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{RST
  Implements a property extractor class for the dark matter only density at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO            , only : darkMatterProfileDMOClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusSpecifier

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDensityDMOProfile" docformat="rst">
   <description>
   A property extractor class for the dark matter only density at a set of radii.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorDensityDMOProfile
     !!{RST
     A property extractor class for the dark matter only density at a set of radii.
     !!}
     private
     class  (darkMatterHaloScaleClass ), pointer                   :: darkMatterHaloScale_          => null()
     class  (darkMatterProfileDMOClass), pointer                   :: darkMatterProfileDMO_         => null()
     integer                                                       :: radiiCount                             , elementCount_
     logical                                                       :: includeRadii
     type   (varying_string           ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusSpecifier          ), allocatable, dimension(:) :: radii
     logical                                                       :: darkMatterScaleRadiusIsNeeded          , diskIsNeeded        , &
          &                                                           spheroidIsNeeded                       , virialRadiusIsNeeded, &
          &                                                           nuclearStarClusterIsNeeded             , satelliteIsNeeded   , &
          &                                                           hotHaloIsNeeded
   contains
     final     ::                       densityDMOProfileDestructor
     procedure :: columnDescriptions => densityDMOProfileColumnDescriptions
     procedure :: size               => densityDMOProfileSize
     procedure :: elementCount       => densityDMOProfileElementCount
     procedure :: extract            => densityDMOProfileExtract
     procedure :: names              => densityDMOProfileNames
     procedure :: descriptions       => densityDMOProfileDescriptions
     procedure :: unitsInSI          => densityDMOProfileUnitsInSI
     procedure :: units              => densityDMOProfileUnits
  end type nodePropertyExtractorDensityDMOProfile

  interface nodePropertyExtractorDensityDMOProfile
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorDensityDMOProfile` output analysis class.
     !!}
     module procedure densityDMOProfileConstructorParameters
     module procedure densityDMOProfileConstructorInternal
  end interface nodePropertyExtractorDensityDMOProfile

contains

  function densityDMOProfileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorDensityDMOProfile` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorDensityDMOProfile)                              :: self
    type   (inputParameters                       ), intent(inout)               :: parameters
    type   (varying_string                        ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass              ), pointer                     :: darkMatterHaloScale_
    class  (darkMatterProfileDMOClass             ), pointer                     :: darkMatterProfileDMO_
    logical                                                                      :: includeRadii

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter docformat="rst">
      <name>radiusSpecifiers</name>
      <description>
      A list of radius specifiers at which to output the density profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>
      Specifies whether or not the radii at which density data are output should also be included in the output file.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=nodePropertyExtractorDensityDMOProfile(radiusSpecifiers,includeRadii,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function densityDMOProfileConstructorParameters

  function densityDMOProfileConstructorInternal(radiusSpecifiers,includeRadii,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorDensityDMOProfile` property extractor class.
    !!}
    use :: Error                               , only : Error_Report
    use :: Galactic_Structure_Options          , only : componentTypeDarkMatterOnly               , massTypeDark, massTypeAll    
    use :: Galactic_Structure_Radii_Definitions, only : Galactic_Structure_Radii_Definition_Decode
    implicit none
    type   (nodePropertyExtractorDensityDMOProfile)                              :: self
    type   (varying_string                        ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass              ), intent(in   ), target       :: darkMatterHaloScale_
    class  (darkMatterProfileDMOClass             ), intent(in   ), target       :: darkMatterProfileDMO_
    logical                                        , intent(in   )               :: includeRadii
    !![
    <constructorAssign variables="radiusSpecifiers, includeRadii, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
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
    if (any(self%radii%component /= componentTypeDarkMatterOnly                                    )) call Error_Report('only the dark halo component can be output'//{introspection:location})
    if (any(self%radii%mass      /= massTypeDark                .and.self%radii%mass /= massTypeAll)) call Error_Report('only the dark matter can be output'        //{introspection:location})
    return
  end function densityDMOProfileConstructorInternal

  subroutine densityDMOProfileDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorDensityDMOProfile` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDensityDMOProfile), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine densityDMOProfileDestructor

  integer function densityDMOProfileElementCount(self,time)
    !!{RST
    Return the number of elements in the ``densityDMOProfile`` property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityDMOProfile), intent(inout) :: self
    double precision                                        , intent(in   ) :: time
    !$GLC attributes unused :: time

    densityDMOProfileElementCount=self%elementCount_
    return
  end function densityDMOProfileElementCount

  function densityDMOProfileSize(self,time)
    !!{RST
    Return the number of array elements in the ``densityDMOProfile`` property extractors.
    !!}
    implicit none
    integer         (c_size_t                              )                :: densityDMOProfileSize
    class           (nodePropertyExtractorDensityDMOProfile), intent(inout) :: self
    double precision                                        , intent(in   ) :: time
    !$GLC attributes unused :: time

    densityDMOProfileSize=self%radiiCount
    return
  end function densityDMOProfileSize

  function densityDMOProfileExtract(self,node,time,instance)
    !!{RST
    Implement a ``densityDMOProfile`` property extractor.
    !!}
    use :: Galactic_Structure_Options          , only : componentTypeAll               , massTypeGalactic            , massTypeStellar
    use :: Galactic_Structure_Radii_Definitions, only : radiusTypeDarkMatterScaleRadius, radiusTypeDiskHalfMassRadius, radiusTypeDiskRadius                      , radiusTypeGalacticLightFraction   , &
          &                                             radiusTypeGalacticMassFraction , radiusTypeRadius            , radiusTypeSpheroidHalfMassRadius          , radiusTypeSpheroidRadius          , &
          &                                             radiusTypeStellarMassFraction  , radiusTypeVirialRadius      , radiusTypeNuclearStarClusterHalfMassRadius, radiusTypeNuclearStarClusterRadius, &
          &                                             radiusTypeHotHaloOuterRadius
    use :: Galacticus_Nodes                    , only : nodeComponentDarkMatterProfile , nodeComponentDisk           , nodeComponentSpheroid                     , nodeComponentNSC                  , &
         &                                              nodeComponentHotHalo           , treeNode
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Coordinates                         , only : coordinateSpherical            , assignment(=)
    use :: Numerical_Constants_Math            , only : Pi
    use :: Error                               , only : Error_Report
    implicit none
    double precision                                        , dimension(:,:), allocatable :: densityDMOProfileExtract
    class           (nodePropertyExtractorDensityDMOProfile), intent(inout) , target      :: self
    type            (treeNode                              ), intent(inout) , target      :: node
    double precision                                        , intent(in   )               :: time
    type            (multiCounter                          ), intent(inout) , optional    :: instance
    class           (nodeComponentHotHalo                  ), pointer                     :: hotHalo
    class           (nodeComponentDisk                     ), pointer                     :: disk
    class           (nodeComponentSpheroid                 ), pointer                     :: spheroid
    class           (nodeComponentNSC                      ), pointer                     :: nuclearStarCluster
    class           (nodeComponentDarkMatterProfile        ), pointer                     :: darkMatterProfile
    class           (massDistributionClass                 ), pointer                     :: massDistribution_
    type            (coordinateSpherical                   )                              :: coordinates
    integer                                                                               :: i
    double precision                                                                      :: radius                , radiusVirial
    !$GLC attributes unused :: time, instance

    allocate(densityDMOProfileExtract(self%radiiCount,self%elementCount_))
    radiusVirial                                               =  0.0d0
    if (self%         virialRadiusIsNeeded) radiusVirial       =  self%darkMatterHaloScale_%radiusVirial(node                    )
    if (self%              hotHaloIsNeeded) hotHalo            =>                                        node%hotHalo          ()
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
       case   (radiusTypeHotHaloOuterRadius              %ID)
          radius=+radius*hotHalo           %   outerRadius()
       case   (radiusTypeDiskRadius                      %ID)
          radius=+radius*disk              %        radius()
       case   (radiusTypeSpheroidRadius                  %ID)
          radius=+radius*spheroid          %        radius()
       case   (radiusTypeNuclearStarClusterRadius        %ID)
          radius=+radius*nuclearStarCluster%        radius()
       case   (radiusTypeDiskHalfMassRadius              %ID)
          radius=+radius*disk             %halfMassRadius()
       case   (radiusTypeSpheroidHalfMassRadius          %ID)
          radius=+radius*spheroid         %halfMassRadius()
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
       coordinates                          =  [radius,Pi/2.0d0,0.0d0]
       massDistribution_                    => self             %darkMatterProfileDMO_%get    (node       )
       densityDMOProfileExtract       (i,1) =  massDistribution_                      %density(coordinates)
       if (self%includeRadii)                                                                               &
            & densityDMOProfileExtract(i,2) =                                                  radius
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end do
    return
  end function densityDMOProfileExtract

  subroutine densityDMOProfileNames(self,names,time)
    !!{RST
    Return the names of the ``densityDMOProfile`` properties.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityDMOProfile), intent(inout)                             :: self
    double precision                                        , intent(in   ), optional                   :: time
    type            (varying_string                        ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names(1)="densityDMOProfile"
    if (self%includeRadii) names(2)="densityDMOProfileRadius"
    return
  end subroutine densityDMOProfileNames

  subroutine densityDMOProfileDescriptions(self,descriptions,time)
    !!{RST
    Return descriptions of the ``densityDMOProfile`` property.
    !!}
    implicit none
    class           (nodePropertyExtractorDensityDMOProfile), intent(inout)                             :: self
    double precision                                        , intent(in   ), optional                   :: time
    type            (varying_string                        ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions       (1)="Dark matter only density at a given radius [M☉/Mpc³]."
    if (self%includeRadii)                                                             &
         & descriptions(2)="Radius at which dark matter only density is output [Mpc]."
    return
  end subroutine densityDMOProfileDescriptions

  subroutine densityDMOProfileColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnits,time)
    !!{RST
    Return column descriptions of the ``densityDMOProfile`` property.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    class           (nodePropertyExtractorDensityDMOProfile), intent(inout)                            :: self
    double precision                                        , intent(in   ), optional                  :: time
    type            (varying_string                        ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                        , intent(inout), dimension(:), allocatable :: values
    type            (varying_string                        ), intent(  out)                            :: valuesDescription
    type            (unitType                              ), intent(  out)                            :: valuesUnits
    !$GLC attributes unused :: time

    allocate(descriptions(self%radiiCount))
    allocate(values      (              0))
    valuesDescription=var_str('')
    valuesUnits      =unitType(1.0d0)
    descriptions     =self%radii%name
    return
  end subroutine densityDMOProfileColumnDescriptions

  function densityDMOProfileUnitsInSI(self,time)
    !!{RST
    Return the units of the ``densityDMOProfile`` properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                        , allocatable  , dimension(:) :: densityDMOProfileUnitsInSI
    class           (nodePropertyExtractorDensityDMOProfile), intent(inout)               :: self
    double precision                                        , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(densityDMOProfileUnitsInSI(self%elementCount_))
    densityDMOProfileUnitsInSI       (1)=massSolar/megaParsec**3
    if (self%includeRadii)                                       &
         & densityDMOProfileUnitsInSI(2)=          megaParsec
    return
  end function densityDMOProfileUnitsInSI

  function densityDMOProfileUnits(self,time) result(units)
    !!{RST
    Return the units of the ``densityDMOProfile`` properties.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    use :: Units_MetaData                  , only : unitType
    implicit none
    type            (unitType                              ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorDensityDMOProfile), intent(inout)             :: self
    double precision                                        , intent(in   ), optional   :: time
    !$GLC attributes unused :: time

    allocate(units(self%elementCount_))
    units       (1)=unitType(massSolar/megaParsec**3,description='M☉/Mpc³',quantity='solMass/Mpc^3')
    if (self%includeRadii)                                                                           &
         & units(2)=unitType(          megaParsec   ,description='Mpc'    ,quantity='Mpc'          )
    return
  end function densityDMOProfileUnits


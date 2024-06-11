!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
Contains a module which implements the standard galactic structure functions.
!!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles      , only : darkMatterProfileClass
  use :: Root_Finder               , only : rootFinder
  use :: Kind_Numbers              , only : kind_int8
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType, enumerationMassTypeType, enumerationWeightByType, enumerationStructureErrorCodeType
  !![
  <galacticStructure name="galacticStructureStandard">
   <description>
     The standard galactic structure functions.
   </description>
  </galacticStructure>
  !!]
  type, extends(galacticStructureClass) :: galacticStructureStandard
     !!{
     The standard galactic structure functions.
     !!}
     private
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_                   => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_                  => null()
     type            (rootFinder              )          :: finderSurfaceDensity
     double precision                                    :: radiusEnclosingSurfaceDensityPrevious          , potentialOffset
     integer         (kind_int8               )          :: uniqueIDPrevious
     logical                                             :: potentialOffsetComputed
   contains
     !![
     <methods>
       <method method="defaults"         description="Set default options."        />
       <method method="restore"          description="Restore defaults from stack."/>
       <method method="calculationReset" description="Reset memoized calculations."/>
     </methods>
     !!]
     final     ::                                  standardDestructor
     procedure :: autoHook                      => standardAutoHook
     procedure :: density                       => standardDensity
     procedure :: densitySphericalAverage       => standardDensitySphericalAverage
     procedure :: massEnclosed                  => standardMassEnclosed
     procedure :: radiusEnclosingMass           => standardRadiusEnclosingMass
     procedure :: velocityRotation              => standardVelocityRotation
     procedure :: velocityRotationGradient      => standardVelocityRotationGradient
     procedure :: potential                     => standardPotential
     procedure :: surfaceDensity                => standardSurfaceDensity
     procedure :: radiusEnclosingSurfaceDensity => standardRadiusEnclosingSurfaceDensity
     procedure :: acceleration                  => standardAcceleration
     procedure :: tidalTensor                   => standardTidalTensor
     procedure :: chandrasekharIntegral         => standardChandrasekharIntegral
     procedure :: velocityDispersion            => standardVelocityDispersion
     procedure :: defaults                      => standardDefaults
     procedure :: restore                       => standardRestore
     procedure :: calculationReset              => standardCalculationReset
  end type galacticStructureStandard

  interface galacticStructureStandard
     !!{
     Constructors for the ``standard'' galactic structure class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface galacticStructureStandard

  ! Types used to store state to allow recursive calling of these functions.
  type :: galacticStructureState
     type            (enumerationComponentTypeType) :: componentType_
     type            (enumerationMassTypeType     ) :: massType_
     type            (enumerationWeightByType     ) :: weightBy_
     integer                                        :: weightIndex_
     double precision                               :: radius_
  end type galacticStructureState

  type :: galacticStructureStateList     
     type(galacticStructureState    )          :: state
     type(galacticStructureStateList), pointer :: next  => null(), previous => null()
  end type galacticStructureStateList

  ! State stack used to allow recursive calling of these functions.
  type            (galacticStructureStateList       ), pointer      :: galacticStructureState_, galacticStructureStateHead_
  !$omp threadprivate(galacticStructureState_,galacticStructureStateHead_)

  ! Submodule-scope variables used in callback functions and root-finding.
  type            (enumerationStructureErrorCodeType),              :: status_
  double precision                                   , dimension(3) :: positionSpherical_     , positionCartesian_ , &
       &                                                               positionCylindrical_   , velocityCartesian_
  !$omp threadprivate(status_,positionSpherical_,positionCartesian_,positionCylindrical_,velocityCartesian_)

  ! Submodule-scope variables used in root finding.
  double precision                                                  :: surfaceDensityTarget
  type            (treeNode                         ), pointer      :: node_
  class           (galacticStructureStandard        ), pointer      :: self_
  !$omp threadprivate(self_,node_,surfaceDensityTarget)

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``standard'' galactic structure class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (galacticStructureStandard)                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_
    class(darkMatterHaloScaleClass ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=galacticStructureStandard(cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
   return
  end function standardConstructorParameters

  function standardConstructorInternal(cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the ``standard'' galactic structure class.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type(galacticStructureStandard)                        :: self
    class(cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    class(darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]

    self%potentialOffsetComputed              =.false.
    self%radiusEnclosingSurfaceDensityPrevious=-huge(0.0d0)
    self%uniqueIDPrevious                     =-1_kind_int8
    self%finderSurfaceDensity                 =rootFinder(                                                             &
            &                                             rootFunction                 =surfaceDensityRoot           , &
            &                                             rangeExpandDownward          =0.5d0                        , &
            &                                             rangeExpandUpward            =2.0d0                        , &
            &                                             rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                                             rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                                             rangeExpandType              =rangeExpandMultiplicative    , &
            &                                             toleranceAbsolute            =0.0d+0                       , &
            &                                             toleranceRelative            =1.0d-6                         &
            &                                            )
    return
  end function standardConstructorInternal

  subroutine standardAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(galacticStructureStandard), intent(inout) :: self

    call calculationResetEvent%attach(self,standardCalculationReset,openMPThreadBindingAllLevels,label='galacticStructureStandard')
    return
  end subroutine standardAutoHook

  subroutine standardDestructor(self)
    !!{
    Destructor for the ``standard'' galactic structure class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(galacticStructureStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    if (calculationResetEvent%isAttached(self,standardCalculationReset)) call calculationResetEvent%detach(self,standardCalculationReset)
    return
  end subroutine standardDestructor

  subroutine standardCalculationReset(self,node,uniqueID)
    !!{
    Reset calculations for galactic structure potentials.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    class  (galacticStructureStandard), intent(inout) :: self
    type   (treeNode                 ), intent(in   ) :: node
    integer(kind_int8                ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%potentialOffsetComputed              =.false.
    self%radiusEnclosingSurfaceDensityPrevious=-huge(0.0d0)
    self%uniqueIDPrevious                     =uniqueID
    return
  end subroutine standardCalculationReset

  double precision function standardDensity(self,node,position,coordinateSystem,componentType,massType,weightBy,weightIndex) result(density)
    !!{
    Compute the density (of given {\normalfont \ttfamily massType}) at the specified {\normalfont \ttfamily position}. Assumes that galactic structure has already
    been computed.
    !!}
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Coordinate_Systems        , only : Coordinates_Cartesian_To_Spherical, Coordinates_Cylindrical_To_Spherical
    use :: Galactic_Structure_Options, only : coordinateSystemCartesian         , coordinateSystemCylindrical         , coordinateSystemSpherical, enumerationCoordinateSystemType, &
         &                                    massTypeAll                       , componentTypeAll
    use :: Coordinates               , only : assignment(=)                     , coordinateSpherical
    use :: Error                     , only : Error_Report
    implicit none
    class           (galacticStructureStandard      ), intent(inout)               :: self
    type            (treeNode                       ), intent(inout)               :: node
    type            (enumerationComponentTypeType   ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional     :: massType
    type            (enumerationWeightByType        ), intent(in   ), optional     :: weightBy
    integer                                          , intent(in   ), optional     :: weightIndex
    type            (enumerationCoordinateSystemType), intent(in   ), optional     :: coordinateSystem
    double precision                                 , intent(in   ), dimension(3) :: position
    class           (massDistributionClass          ), pointer                     :: massDistribution_
    type            (enumerationCoordinateSystemType)                              :: coordinateSystemActual
    type            (coordinateSpherical            )                              :: position_

    ! Determine position in spherical coordinate system to use.
    if (present(coordinateSystem)) then
       coordinateSystemActual=coordinateSystem
    else
       coordinateSystemActual=coordinateSystemSpherical
    end if
    select case (coordinateSystemActual%ID)
    case (coordinateSystemSpherical  %ID)
       positionSpherical_=position
    case (coordinateSystemCylindrical%ID)
       positionSpherical_=Coordinates_Cylindrical_To_Spherical(position)
    case (coordinateSystemCartesian  %ID)
       positionSpherical_=Coordinates_Cartesian_To_Spherical  (position)
    case default
       call Error_Report('unknown coordinate system type'//{introspection:location})
    end select
    call self%defaults(componentType=componentType,massType=massType,weightBy=weightBy,weightIndex=weightIndex)
    ! Evaluate the density.
    position_         =  positionSpherical_
    massDistribution_ => node              %massDistribution(galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    density           =  massDistribution_ %density         (position_                                                                                                                                                              )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    call self%restore()
    return
  end function standardDensity

  double precision function standardDensitySphericalAverage(self,node,radius,componentType,massType,weightBy,weightIndex) result(density)
    !!{
    Compute the density (of given {\normalfont \ttfamily massType}) at the specified {\normalfont \ttfamily position}. Assumes that galactic structure has already
    been computed.
    !!}
    use :: Galactic_Structure_Options, only : massTypeAll          , componentTypeAll
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (galacticStructureStandard       ), intent(inout)           :: self
    type            (treeNode                        ), intent(inout)           :: node
    double precision                                  , intent(in   )           :: radius
    type            (enumerationComponentTypeType    ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType         ), intent(in   ), optional :: massType
    type            (enumerationWeightByType         ), intent(in   ), optional :: weightBy
    integer                                           , intent(in   ), optional :: weightIndex
    class           (massDistributionClass           ), pointer                 :: massDistribution_
    
    call self%defaults(radius=radius,componentType=componentType,massType=massType,weightBy=weightBy,weightIndex=weightIndex)
    ! Compute the spherically-averaged density.
    massDistribution_ => node             %massDistribution       (galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    density           =  massDistribution_%densitySphericalAverage(radius                                                                                                                                                                 )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    call self%restore()
    return
  end function standardDensitySphericalAverage
  
  double precision function standardMassEnclosed(self,node,radius,componentType,massType,weightBy,weightIndex) result(massEnclosed)
    !!{
    Compute the mass within a given radius, or the total mass if no radius is specified.
    !!}
    use :: Galactic_Structure_Options, only : massTypeAll          , componentTypeAll
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (galacticStructureStandard   ), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    double precision                              , intent(in   ), optional :: radius
    class           (massDistributionClass       ), pointer                 :: massDistribution_

    call self%defaults(radius,componentType,massType,weightBy,weightIndex)
    ! Compute the contribution from components directly, by mapping a function over all components.
    massDistribution_ => node             %massDistribution    (galacticStructureState_%state%componentType_  ,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    massEnclosed      =  massDistribution_%massEnclosedBySphere(galacticStructureState_%state%radius_)
    ! Call routines to supply the masses for all components.
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    call self%restore()
    return
  end function standardMassEnclosed

  double precision function standardRadiusEnclosingMass(self,node,mass,massFractional,componentType,massType,weightBy,weightIndex)
    !!{
    Return the radius enclosing a given mass (or fractional mass) in {\normalfont \ttfamily node}.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (galacticStructureStandard   ), intent(inout), target   :: self
    type            (treeNode                    ), intent(inout), target   :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    double precision                              , intent(in   ), optional :: massFractional, mass
    class           (massDistributionClass       ), pointer                 :: massDistribution_

    call self%defaults(componentType=componentType,massType=massType,weightBy=weightBy,weightIndex=weightIndex)
    massDistribution_           => node             %massDistribution   (galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    standardRadiusEnclosingMass =  massDistribution_%radiusEnclosingMass(mass,massFractional                                                                                                                                                    )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    call self%restore()
    return
  end function standardRadiusEnclosingMass

  double precision function standardVelocityRotation(self,node,radius,componentType,massType) result(velocityRotation)
    !!{
    Compute the rotation curve a given radius.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (galacticStructureStandard   ), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                              , intent(in   )           :: radius
    class           (massDistributionClass       ), pointer                 :: massDistribution_
    double precision                                                        :: rotationCurveSquared

    call self%defaults(radius=radius,componentType=componentType,massType=massType)
    massDistribution_    => node             %massDistribution(galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    rotationCurveSquared =  massDistribution_%rotationCurve   (galacticStructureState_%state%radius_)**2
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! We've added velocities in quadrature, so now take the square root.
    velocityRotation=sqrt(rotationCurveSquared)
    call self%restore()
    return
  end function standardVelocityRotation

  double precision function standardVelocityRotationGradient(self,node,radius,componentType,massType) result(velocityRotationGradient)
    !!{
    Solve for the rotation curve gradient at a given radius. Assumes the galactic structure has already been computed.
    !!}
    use :: Error             , only : Error_Report
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (galacticStructureStandard   ), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                              , intent(in   )           :: radius
    class           (massDistributionClass       )               , pointer  :: massDistribution_
    double precision                                                        :: velocityRotation

    call self%defaults(radius=radius,componentType=componentType,massType=massType)
    massDistribution_        => node             %massDistribution     (galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    velocityRotationGradient =  massDistribution_%rotationCurveGradient(galacticStructureState_%state%radius_)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Convert the summed dV²/dr to dV/dr.
    velocityRotation=+self%velocityRotation(node,radius,componentType,massType)
    if (velocityRotation > 0.0d0) then
       velocityRotationGradient=+0.5d0                    &
            &                   *velocityRotationGradient &
            &                   /velocityRotation
    else if (velocityRotationGradient /= 0.0d0) then
       call Error_Report('rotation curve is zero, but gradient is non-zero'//{introspection:location})
    end if
    call self%restore()
    return
  end function standardVelocityRotationGradient

  double precision function standardPotential(self,node,radius,componentType,massType,status) result(potential)
    !!{
    Solve for the gravitational potential at a given radius. Assumes the galactic structure has already been computed.
    !!}
    use :: Coordinates               , only : assignment(=)        , coordinateCylindrical
    use :: Galactic_Structure_Options, only : componentTypeAll     , massTypeAll          , structureErrorCodeSuccess
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (galacticStructureStandard        ), intent(inout)           :: self
    type            (treeNode                         ), intent(inout)           :: node
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType
    double precision                                   , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    class           (massDistributionClass            ), pointer                 :: massDistribution_
    type            (coordinateCylindrical            )                          :: position

    ! Initialize status.
    if (present(status)) status=structureErrorCodeSuccess
    ! Reset calculations if this is a new node.
    if (node%uniqueID() /= self%uniqueIDPrevious) call self%calculationReset(node,node%uniqueID())
    ! Evaluate the potential at the halo virial radius.
    if (.not.self%potentialOffsetComputed) then
       call self%defaults(componentType=componentTypeAll,massType=massTypeAll,radius=self%darkMatterHaloScale_%radiusVirial(node))
       massDistribution_ => node%massDistribution(galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
       position =[galacticStructureState_%state%radius_,0.0d0,0.0d0]
       potential=massDistribution_%potential(position,status_)
       if (status_ /= structureErrorCodeSuccess) status=status_
       call self%restore()
       ! Compute the potential offset such that the total gravitational potential at the virial radius is -V² where V is the
       ! virial velocity.
       self%potentialOffset        =-potential-self%darkMatterHaloScale_%velocityVirial(node)**2
       self%potentialOffsetComputed=.true.
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end if
    call self%defaults(radius=radius,componentType=componentType,massType=massType)
   ! Determine which component type to use.
    if (present(componentType)) then
       galacticStructureState_%state%componentType_=componentType
    else
       galacticStructureState_%state%componentType_=componentTypeAll
    end if
    massDistribution_ => node%massDistribution(galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    position =[galacticStructureState_%state%radius_,0.0d0,0.0d0]
    potential=+massDistribution_%potential(position,status_) &
         &    +self%potentialOffset
    if (status_ /= structureErrorCodeSuccess) status=status_
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    call self%restore()
    return
  end function standardPotential

  double precision function standardSurfaceDensity(self,node,position,coordinateSystem,componentType,massType,weightBy,weightIndex) result(surfaceDensity)
    !!{
    Compute the surface density of given {\normalfont \ttfamily massType}) at the specified {\normalfont \ttfamily position}.
    !!}
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Coordinate_Systems        , only : Coordinates_Cartesian_To_Cylindrical, Coordinates_Spherical_To_Cylindrical
    use :: Galactic_Structure_Options, only : coordinateSystemCartesian           , coordinateSystemCylindrical         , coordinateSystemSpherical, enumerationCoordinateSystemType, &
         &                                     massTypeAll                        , componentTypeAll
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : treeNode
    use :: Coordinates               , only : assignment(=)                       , coordinateCylindrical
    implicit none
    class           (galacticStructureStandard      ), intent(inout)               :: self
    type            (treeNode                       ), intent(inout)               :: node
    type            (enumerationComponentTypeType   ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional     :: massType
    type            (enumerationWeightByType        ), intent(in   ), optional     :: weightBy
    integer                                          , intent(in   ), optional     :: weightIndex
    type            (enumerationCoordinateSystemType), intent(in   ), optional     :: coordinateSystem
    double precision                                 , intent(in   ), dimension(3) :: position
    type            (coordinateCylindrical           )                             :: position_
    class           (massDistributionClass          ), pointer                     :: massDistribution_
    !![
    <optionalArgument name="coordinateSystem" defaultsTo="coordinateSystemCylindrical" />
    !!]
       
    call self%defaults(componentType=componentType,massType=massType,weightBy=weightBy,weightIndex=weightIndex)
    ! Determine position in cylindrical coordinate system to use.
    select case (coordinateSystem_%ID)
    case (coordinateSystemSpherical  %ID)
       positionCylindrical_=Coordinates_Spherical_To_Cylindrical (position)
    case (coordinateSystemCylindrical%ID)
       positionCylindrical_=position
    case (coordinateSystemCartesian  %ID)
       positionCylindrical_=Coordinates_Cartesian_To_Cylindrical(position)
    case default
       call Error_Report('unknown coordinate system type'//{introspection:location})
    end select
    ! Compute the surface density.
    position_         =  positionCylindrical_
    massDistribution_ => node                %massDistribution(galacticStructureState_%state%componentType_  ,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    surfaceDensity    =  massDistribution_   %surfaceDensity  (position_)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    call self%restore()
    return
  end function standardSurfaceDensity

  double precision function standardRadiusEnclosingSurfaceDensity(self,node,surfaceDensity,componentType,massType,weightBy,weightIndex) result(radius)
    !!{
    Return the radius enclosing a given surface density in {\normalfont \ttfamily node}.
    !!}
    use :: Kind_Numbers, only : kind_int8
    use :: Root_Finder , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (galacticStructureStandard   ), intent(inout), target   :: self
    type            (treeNode                    ), intent(inout), target   :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    double precision                              , intent(in   )           :: surfaceDensity
    double precision                                                        :: radiusGuess

    call self%defaults(componentType=componentType,massType=massType,weightBy=weightBy,weightIndex=weightIndex)
    self_                => self
    node_                => node
    surfaceDensityTarget =  surfaceDensity
    if (self%radiusEnclosingSurfaceDensityPrevious >= 0.0d0) then
       radiusGuess=self                     %radiusEnclosingSurfaceDensityPrevious
    else
       radiusGuess=self%darkMatterHaloScale_%radiusVirial                         (node)
    end if
    self%radiusEnclosingSurfaceDensityPrevious=self%finderSurfaceDensity%find(rootGuess=radiusGuess)
    radius=self%radiusEnclosingSurfaceDensityPrevious
    call self%restore()
    return
  end function standardRadiusEnclosingSurfaceDensity

  double precision function surfaceDensityRoot(radius)
    !!{
    Root function used in solving for the radius that encloses a given surface density.
    !!}
    use :: Galactic_Structure_Options, only : coordinateSystemCylindrical
    implicit none
    double precision, intent(in   ) :: radius

    ! Evaluate the root function.
    surfaceDensityRoot=+self_%surfaceDensity      (node_,[radius,0.0d0,0.0d0],coordinateSystemCylindrical,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_) &
         &             -      surfaceDensityTarget
    return
  end function surfaceDensityRoot
  
  function standardAcceleration(self,node,positionCartesian,componentType,massType) result(acceleration)
    !!{
    Compute the gravitational acceleration at a given position.
    !!}
    use :: Coordinates       , only : assignment(=)        , coordinateCartesian
    use :: Galacticus_Nodes  , only : treeNode
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (galacticStructureStandard   ), intent(inout)               :: self
    double precision                                             , dimension(3) :: acceleration
    type            (treeNode                    ), intent(inout)               :: node
    double precision                              , intent(in   ), dimension(3) :: positionCartesian
    type            (enumerationComponentTypeType), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional     :: massType
    class           (massDistributionClass       ), pointer                     :: massDistribution_
    type            (coordinateCartesian         )                              :: position

    call self%defaults(componentType=componentType,massType=massType)
    positionCartesian_ =  positionCartesian
    position           =  positionCartesian
    massDistribution_  => node             %massDistribution(galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    acceleration       =  massDistribution_%acceleration    (position)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    call self%restore()
    return
  end function standardAcceleration

  function standardTidalTensor(self,node,positionCartesian,componentType,massType) result(tidalTensor)
    !!{
    Compute the gravitational tidal tensor at a given position.
    !!}
    use :: Galacticus_Nodes  , only : treeNode
    use :: Mass_Distributions, only : massDistributionClass
    use :: Coordinates       , only : assignment(=)                 , coordinateCartesian
    use :: Tensors           , only : tensorRank2Dimension3Symmetric
    implicit none
    class           (galacticStructureStandard     ), intent(inout)               :: self
    type            (tensorRank2Dimension3Symmetric)                              :: tidalTensor
    type            (treeNode                      ), intent(inout)               :: node
    double precision                                , intent(in   ), dimension(3) :: positionCartesian
    type            (enumerationComponentTypeType  ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType       ), intent(in   ), optional     :: massType
    class           (massDistributionClass         ), pointer                     :: massDistribution_
    type            (coordinateCartesian           )                              :: position
    
    call self%defaults(componentType=componentType,massType=massType)
    positionCartesian_ =  positionCartesian
    position           =  positionCartesian
    massDistribution_  => node             %massDistribution(galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    tidalTensor        =  massDistribution_%tidalTensor     (position)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    call self%restore()
    return
  end function standardTidalTensor

  function standardChandrasekharIntegral(self,node,nodeSatellite,positionCartesian,velocityCartesian,massPerturber,componentType,massType) result(chandrasekharIntegral)
    !!{
    Compute the integral appearing in the \cite{chandrasekhar_dynamical_1943} dynamical friction model:
    \begin{equation}
      \rho(\boldsymbol{x}_\mathrm{s}) \int \mathrm{d}\boldsymbol{v} f(\boldsymbol{v}) {\boldsymbol{v}-\boldsymbol{v}_\mathrm{s} \over |\boldsymbol{v}-\boldsymbol{v}_\mathrm{s}|^3},
    \end{equation}  
    where $\rho(\boldsymbol{x}_\mathrm{s})$ is the density at the position of the perturber, $\boldsymbol{x}_\mathrm{s}$,
    $f(\boldsymbol{v})$ is the velocity distribution function at velocity $\boldsymbol{v}$, and $\boldsymbol{v}_\mathrm{s}$ is
    the velocity of the perturber.
    !!}
    use :: Galacticus_Nodes          , only : treeNode
    use :: Galactic_Structure_Options, only : componentTypeAll     , massTypeAll
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Coordinates               , only : assignment(=)        , coordinateCartesian
    implicit none
    double precision                                               , dimension(3) :: chandrasekharIntegral
    class           (galacticStructureStandard     ), intent(inout)               :: self
    type            (treeNode                      ), intent(inout), target       :: node                       , nodeSatellite
    double precision                                , intent(in   ), dimension(3) :: positionCartesian          , velocityCartesian
    double precision                                , intent(in   )               :: massPerturber
    type            (enumerationComponentTypeType  ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType       ), intent(in   ), optional     :: massType
    integer                                         , parameter                   :: chandrasekharIntegralSize=3
    class           (massDistributionClass         ), pointer                     :: massDistribution_          , massDistributionPerturber
    type            (coordinateCartesian           )                              :: position                   , velocity

    call self%defaults(componentType=componentType,massType=massType)
    positionCartesian_ =  positionCartesian
    velocityCartesian_ =  velocityCartesian
    position           =  positionCartesian
    velocity           =  velocityCartesian
    ! Evaluate the density.
    massDistribution_         => node             %massDistribution     (galacticStructureState_%state%componentType_  ,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    massDistributionPerturber => nodeSatellite    %massDistribution     (componentTypeAll,                              massTypeAll,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    chandrasekharIntegral     =  massDistribution_%chandrasekharIntegral(massDistribution_,massDistributionPerturber,massPerturber,position,velocity)
    !![
    <objectDestructor name="massDistribution_"        />
    <objectDestructor name="massDistributionPerturber"/>
    !!]
    call self%restore()
    return
  end function standardChandrasekharIntegral

  double precision function standardVelocityDispersion(self,node,radius,radiusOuter,componentType,massType) result(velocityDispersion)
    !!{
    Returns the velocity dispersion of the specified {\normalfont \ttfamily componentType} in {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius}.
    !!}
    use :: Galactic_Structure_Options, only : radiusLarge
    use :: Numerical_Integration     , only : integrator
    use :: Numerical_Constants_Math  , only : Pi
    class           (galacticStructureStandard   ), intent(inout), target   :: self
    type            (treeNode                    ), intent(inout), target   :: node
    double precision                              , intent(in   )           :: radius                 , radiusOuter
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                                                        :: densitySphericalAverage, densityVelocityVariance, &
         &                                                                     massTotal
    type            (integrator                  )                          :: integrator_

    self_         => self
    node_         => node
    call self%defaults(componentType=componentType,massType=massType)
    massTotal=self%massEnclosed(node,radiusLarge,componentType=componentType,massType=massType)
    ! Return with zero dispersion if the component is massless.
    if (massTotal <= 0.0d0) then
       velocityDispersion=0.0d0
       return
    end if
    ! Integrate the Jeans equation.
    integrator_            =integrator           (integrandVelocityDispersion,toleranceRelative=1.0d-3)
    densityVelocityVariance=integrator_%integrate(radius                     ,radiusOuter             )
    ! Get the density at this radius.
    densitySphericalAverage=self_%densitySphericalAverage(                                                            &
         &                                                              node_                                       , &
         &                                                              radius                                      , & 
         &                                                componentType=galacticStructureState_%state%componentType_, &
         &                                                massType     =galacticStructureState_%state%massType_       &
         &                                               )
    ! Check for zero density.
    if (densitySphericalAverage <= 0.0d0) then
       velocityDispersion=0.0d0
    else
       velocityDispersion=sqrt(max(densityVelocityVariance,0.0d0)/densitySphericalAverage)
    end if
    call self%restore()
    return
  end function standardVelocityDispersion

  double precision function integrandVelocityDispersion(radius)
    !!{
    Integrand function used for finding velocity dispersions using Jeans equation.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    double precision, intent(in   ) :: radius

    if (radius == 0.0d0) then
       integrandVelocityDispersion=0.0d0
    else
       integrandVelocityDispersion=+gravitationalConstantGalacticus                                                           &
            &                      *self_%massEnclosed           (                                                            &
            &                                                                    node_                                      , &
            &                                                                    radius                                       &
            &                                                    )                                                            &
            &                      /radius**2                                                                                 &
            &                      *self_%densitySphericalAverage(                                                            &
            &                                                                   node_                                       , &
            &                                                                   radius                                      , & 
            &                                                     componentType=galacticStructureState_%state%componentType_, &
            &                                                     massType     =galacticStructureState_%state%massType_       &
            &                                                    )
    end if
    return
  end function integrandVelocityDispersion

  subroutine standardDefaults(self,radius,componentType,massType,weightBy,weightIndex)
    !!{
    Set the default values for options in the enclosed mass functions.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll, massTypeAll, weightByLuminosity, weightByMass, &
         &                                    radiusLarge
    use :: Error                     , only : Error_Report
    implicit none
    class           (galacticStructureStandard   ), intent(inout)              :: self
    double precision                              , intent(in   ), optional    :: radius
    type            (enumerationComponentTypeType), intent(in   ), optional    :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional    :: massType
    type            (enumerationWeightByType     ), intent(in   ), optional    :: weightBy
    integer                                       , intent(in   ), optional    :: weightIndex
    !![
    <optionalArgument name="componentType" defaultsTo="componentTypeAll" />  
    <optionalArgument name="massType"      defaultsTo="massTypeAll"      />
    <optionalArgument name="weightBy"      defaultsTo="weightByMass"     />
    <optionalArgument name="weightIndex"   defaultsTo="-1"               />
    <optionalArgument name="radius"        defaultsTo="radiusLarge"      />
    !!]

    ! Expand the state stack if necessary.
    if (.not.associated(galacticStructureState_)) then
       if (.not.associated(galacticStructureStateHead_)) allocate(galacticStructureStateHead_)
       galacticStructureState_ => galacticStructureStateHead_
    else
       if (.not.associated(galacticStructureState_%next)) then
          allocate(galacticStructureState_%next)
          galacticStructureState_%next%previous => galacticStructureState_
       end if
       galacticStructureState_ => galacticStructureState_%next
    end if
    ! Set defaults.    
    galacticStructureState_%state%radius_        =radius_
    galacticStructureState_%state%massType_      =massType_
    galacticStructureState_%state%componentType_ =componentType_
    galacticStructureState_%state%weightBy_      =weightBy_
    select case (weightBy_%ID)
    case (weightByLuminosity%ID)
       if (.not.present(weightIndex)) call Error_Report('weightIndex should be specified for luminosity weighting'//{introspection:location})
       galacticStructureState_%state%weightIndex_=weightIndex_
    end select
    return
  end subroutine standardDefaults

  subroutine standardRestore(self)
    !!{
    Restore the previous state from the stack.
    !!}
    implicit none
    class(galacticStructureStandard), intent(inout) :: self
    !$GLC attributes unused :: self

    galacticStructureState_ => galacticStructureState_%previous
    return
  end subroutine standardRestore

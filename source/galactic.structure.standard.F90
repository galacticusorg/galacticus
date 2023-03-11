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
     class           (darkMatterProfileClass  ), pointer :: darkMatterProfile_                    => null()
     type            (rootFinder              )          :: finderMass                                     , finderSurfaceDensity
     double precision                                    :: radiusEnclosingMassPrevious                    , potentialOffset     , &
          &                                                 radiusEnclosingSurfaceDensityPrevious
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
  type            (galacticStructureStateList       ), pointer :: galacticStructureState_, galacticStructureStateHead_
  !$omp threadprivate(galacticStructureState_,galacticStructureStateHead_)

  ! Submodule-scope variables used in callback functions and root-finding.
  type            (enumerationStructureErrorCodeType),                            :: status_
  double precision                                   , dimension(3)               :: positionSpherical_           , positionCartesian_ , &
       &                                                                             positionCylindrical_         , velocityCartesian_
  !$omp threadprivate(status_,positionSpherical_,positionCartesian_,positionCylindrical_,velocityCartesian_)

  ! Submodule-scope variables used in root finding.
  double precision                                                                :: massTarget                   , surfaceDensityTarget
  type            (treeNode                         ), pointer                    :: node_                        , nodeSatellite_
  class           (galacticStructureStandard        ), pointer                    :: self_
  !$omp threadprivate(self_,node_,nodeSatellite_,massTarget,surfaceDensityTarget)

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
    class(darkMatterProfileClass   ), pointer       :: darkMatterProfile_

    !![
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="darkMatterProfile"   name="darkMatterProfile_"   source="parameters"/>
    !!]
    self=galacticStructureStandard(cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfile_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="darkMatterProfile_"  />
    !!]
   return
  end function standardConstructorParameters

  function standardConstructorInternal(cosmologyFunctions_,darkMatterHaloScale_,darkMatterProfile_) result(self)
    !!{
    Internal constructor for the ``standard'' galactic structure class.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type(galacticStructureStandard)                        :: self
    class(cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    class(darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    class(darkMatterProfileClass  ), intent(in   ), target :: darkMatterProfile_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *darkMatterHaloScale_, *darkMatterProfile_"/>
    !!]

    self%potentialOffsetComputed              =.false.
    self%radiusEnclosingMassPrevious          =-huge(0.0d0)
    self%radiusEnclosingSurfaceDensityPrevious=-huge(0.0d0)
    self%uniqueIDPrevious                     =-1_kind_int8
    self%finderMass                           =rootFinder(                                                             &
               &                                          rootFunction                 =massEnclosedRoot             , &
               &                                          rangeExpandDownward          =0.5d0                        , &
               &                                          rangeExpandUpward            =2.0d0                        , &
               &                                          rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                                          rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                                          rangeExpandType              =rangeExpandMultiplicative    , &
               &                                          toleranceAbsolute            =0.0d+0                       , &
               &                                          toleranceRelative            =1.0d-6                         &
               &                                         )
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
    <objectDestructor name="self%darkMatterProfile_"  />
    !!]
    if (calculationResetEvent%isAttached(self,standardCalculationReset)) call calculationResetEvent%detach(self,standardCalculationReset)
    return
  end subroutine standardDestructor

  subroutine standardCalculationReset(self,node)
    !!{
    Reset calculations for galactic structure potentials.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(galacticStructureStandard), intent(inout) :: self
    type (treeNode                 ), intent(in   ) :: node

    self%potentialOffsetComputed              =.false.
    self%radiusEnclosingMassPrevious          =-huge(0.0d0)
    self%radiusEnclosingSurfaceDensityPrevious=-huge(0.0d0)
    self%uniqueIDPrevious                     =node%uniqueID()
    return
  end subroutine standardCalculationReset

  double precision function standardDensity(self,node,position,coordinateSystem,componentType,massType,weightBy,weightIndex) result(density)
    !!{
    Compute the density (of given {\normalfont \ttfamily massType}) at the specified {\normalfont \ttfamily position}. Assumes that galactic structure has already
    been computed.
    !!}
    use :: Coordinate_Systems        , only : Coordinates_Cartesian_To_Spherical, Coordinates_Cylindrical_To_Spherical
    use :: Galactic_Structure_Options, only : componentTypeAll                  , coordinateSystemCartesian           , coordinateSystemCylindrical, coordinateSystemSpherical      , &
          &                                   massTypeAll                       , weightByLuminosity                  , weightByMass               , enumerationCoordinateSystemType
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : optimizeForDensitySummation       , reductionSummation                  , treeNode
    !![
    <include directive="densityTask" type="moduleUse">
    !!]
    include 'galactic_structure.density.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureStandard      ), intent(inout)               :: self
    type            (treeNode                       ), intent(inout)               :: node
    type            (enumerationComponentTypeType   ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional     :: massType
    type            (enumerationWeightByType        ), intent(in   ), optional     :: weightBy
    integer                                          , intent(in   ), optional     :: weightIndex
    type            (enumerationCoordinateSystemType), intent(in   ), optional     :: coordinateSystem
    double precision                                 , intent(in   ), dimension(3) :: position
    procedure       (densityComponent               ), pointer                     :: densityComponent_
    type            (enumerationCoordinateSystemType)                              :: coordinateSystemActual
    double precision                                                               :: densityComponent__

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
    ! Call routines to supply the densities for all components.
    densityComponent_ => densityComponent
    density           =  node%mapDouble0(densityComponent_,reductionSummation,optimizeFor=optimizeForDensitySummation)
    !![
    <include directive="densityTask" type="functionCall" functionType="function" returnParameter="densityComponent__">
     <functionArgs>node,positionSpherical_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_</functionArgs>
     <onReturn>density=density+densityComponent__</onReturn>
    !!]
    include 'galactic_structure.density.tasks.inc'
    !![
    </include>
    !!]
    call self%restore()
    return
  end function standardDensity

  double precision function densityComponent(component)
    !!{
    Unary function returning the density in a component. Suitable for mapping over components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    densityComponent=component%density(positionSpherical_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    return
  end function densityComponent

  double precision function standardDensitySphericalAverage(self,node,radius,componentType,massType,weightBy,weightIndex) result(density)
    !!{
    Compute the density (of given {\normalfont \ttfamily massType}) at the specified {\normalfont \ttfamily position}. Assumes that galactic structure has already
    been computed.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll                           , massTypeAll       , weightByLuminosity, weightByMass
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : optimizeForDensitySphericalAverageSummation, reductionSummation, treeNode
    !![
    <include directive="densitySphericalAverageTask" type="moduleUse">
    !!]
    include 'galactic_structure.density_spherical_average.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureStandard       ), intent(inout)           :: self
    type            (treeNode                        ), intent(inout)           :: node
    double precision                                  , intent(in   )           :: radius
    type            (enumerationComponentTypeType    ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType         ), intent(in   ), optional :: massType
    type            (enumerationWeightByType         ), intent(in   ), optional :: weightBy
    integer                                           , intent(in   ), optional :: weightIndex
    procedure       (densitySphericalAverageComponent), pointer                 :: densitySphericalAverageComponent_
    double precision                                                            :: densitySphericalAverageComponent__

    call self%defaults(radius=radius,componentType=componentType,massType=massType,weightBy=weightBy,weightIndex=weightIndex)
    ! Call routines to supply the densities for all components.
    densitySphericalAverageComponent_ => densitySphericalAverageComponent
    density                           =  node%mapDouble0(densitySphericalAverageComponent_,reductionSummation,optimizeFor=optimizeForDensitySphericalAverageSummation)
    !![
    <include directive="densitySphericalAverageTask" type="functionCall" functionType="function" returnParameter="densitySphericalAverageComponent__">
     <functionArgs>node,radius,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_</functionArgs>
     <onReturn>density=density+densitySphericalAverageComponent__</onReturn>
    !!]
    include 'galactic_structure.density_spherical_average.tasks.inc'
    !![
    </include>
    !!]
    call self%restore()
    return
  end function standardDensitySphericalAverage
  
  double precision function densitySphericalAverageComponent(component)
    !!{
    Unary function returning the spherically-averaged density in a component. Suitable for mapping over components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    densitySphericalAverageComponent=component%densitySphericalAverage(galacticStructureState_%state%radius_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    return
  end function densitySphericalAverageComponent

  double precision function standardMassEnclosed(self,node,radius,componentType,massType,weightBy,weightIndex) result(massEnclosed)
    !!{
    Compute the mass within a given radius, or the total mass if no radius is specified.
    !!}
    use :: Galacticus_Nodes, only : optimizeForEnclosedMassSummation, reductionSummation
    !![
    <include directive="enclosedMassTask" type="moduleUse">
    !!]
    include 'galactic_structure.enclosed_mass.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureStandard   ), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    double precision                              , intent(in   ), optional :: radius
    procedure       (massComponentEnclosed       ), pointer                 :: massComponentEnclosed_
    double precision                                                        :: massComponent

    call self%defaults(radius,componentType,massType,weightBy,weightIndex)
    ! Compute the contribution from components directly, by mapping a function over all components.
    massComponentEnclosed_ => massComponentEnclosed
    massEnclosed           =  node                 %mapDouble0(massComponentEnclosed_,reductionSummation,optimizeFor=optimizeForEnclosedMassSummation)
    ! Call routines to supply the masses for all components.
    !![
    <include directive="enclosedMassTask" type="functionCall" functionType="function" returnParameter="massComponent">
     <functionArgs>node,galacticStructureState_%state%radius_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_</functionArgs>
     <onReturn>massEnclosed=massEnclosed+massComponent</onReturn>
    !!]
    include 'galactic_structure.enclosed_mass.tasks.inc'
    !![
    </include>
    !!]
    call self%restore()
    return
  end function standardMassEnclosed

  double precision function massComponentEnclosed(component)
    !!{
    Unary function returning the enclosed mass in a component. Suitable for mapping over components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    massComponentEnclosed=component%enclosedMass(galacticStructureState_%state%radius_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    return
  end function massComponentEnclosed

  double precision function standardRadiusEnclosingMass(self,node,mass,massFractional,componentType,massType,weightBy,weightIndex)
    !!{
    Return the radius enclosing a given mass (or fractional mass) in {\normalfont \ttfamily node}.
    !!}
    use :: Dark_Matter_Profile_Structure_Tasks, only : Dark_Matter_Profile_Radius_Enclosing_Mass
    use :: Display                            , only : displayMessage                           , verbosityLevelWarn
    use :: Galactic_Structure_Options         , only : componentTypeDarkHalo                    , massTypeDark
    use :: Error                              , only : Error_Report
    use :: ISO_Varying_String                 , only : assignment(=)                            , operator(//)      , varying_string
    use :: String_Handling                    , only : operator(//)
    implicit none
    class           (galacticStructureStandard   ), intent(inout), target   :: self
    type            (treeNode                    ), intent(inout), target   :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type            (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    double precision                              , intent(in   ), optional :: massFractional, mass
    double precision                                                        :: radiusGuess
    type            (varying_string              )                          :: message
    character       (len=11                      )                          :: massLabel

    call self%defaults(componentType=componentType,massType=massType,weightBy=weightBy,weightIndex=weightIndex)
    ! Determine what mass to use.
    if (present(mass)) then
       if (present(massFractional)) call Error_Report('only one mass or massFractional can be specified'//{introspection:location})
       massTarget=mass
    else if (present(massFractional)) then
       massTarget=massFractional*self%massEnclosed(node,componentType=galacticStructureState_%state%componentType_,massType=galacticStructureState_%state%massType_,weightBy=galacticStructureState_%state%weightBy_,weightIndex=galacticStructureState_%state%weightIndex_)
    else
       call Error_Report('either mass or massFractional must be specified'//{introspection:location})
    end if
    if (massTarget <= 0.0d0) then
       standardRadiusEnclosingMass=0.0d0
       call self%restore()
      return
    end if
    self_ => self
    node_ => node
    ! If the dark matter component is queried and its density profile is unaffected by baryons, compute the radius from the dark
    ! matter profile. Otherwise, find the radius numerically.
    if     (                                                                       &
         &   galacticStructureState_%state%componentType_ == componentTypeDarkHalo &
         &  .or.                                                                   &
         &   galacticStructureState_%state%massType_      == massTypeDark          &
         & ) then
       if (.not.associated(self%darkMatterProfile_)) call Error_Report('object is not expecting dark matter requests'//{introspection:location})       
       ! Use the function provided by the dark matter profile structure tasks module here. This ensures precise consistency
       ! between calculations here and in the enclosed mass function.
       standardRadiusEnclosingMass=Dark_Matter_Profile_Radius_Enclosing_Mass(node_,massTarget)
    else
       ! Solve for the radius.
       if (massEnclosedRoot(0.0d0) >= 0.0d0) then
          message='Enclosed mass in galaxy (ID='
          write (massLabel,'(e10.4)') self%massEnclosed(node_,0.0d0,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
          message=message//node%index()//') seems to be finite ('//trim(massLabel)
          write (massLabel,'(e10.4)') massTarget
          message=message//') at zero radius (was seeking '      //trim(massLabel)
          message=message//') - returning zero radius.'
          call displayMessage(message,verbosityLevelWarn)
          standardRadiusEnclosingMass=0.0d0
          call self%restore()
          return
       end if
       if (self%radiusEnclosingMassPrevious >= 0.0d0) then
          radiusGuess=self                     %radiusEnclosingMassPrevious
       else
          radiusGuess=self%darkMatterHaloScale_%radiusVirial               (node)
       end if
       self%radiusEnclosingMassPrevious=self%finderMass%find                       (rootGuess=radiusGuess)
       self%uniqueIDPrevious           =node           %uniqueID                   (                     )
       standardRadiusEnclosingMass     =self           %radiusEnclosingMassPrevious
    end if
    call self%restore()
    return
  end function standardRadiusEnclosingMass

  double precision function massEnclosedRoot(radius)
    !!{
    Root function used in solving for the radius that encloses a given mass.
    !!}
    double precision, intent(in   ) :: radius

    massEnclosedRoot=+self_%massEnclosed(node_,radius,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_) &
         &           -      massTarget
    return
  end function massEnclosedRoot

  double precision function standardVelocityRotation(self,node,radius,componentType,massType) result(velocityRotation)
    !!{
    Compute the rotation curve a given radius.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll     , massTypeAll
    use :: Galacticus_Nodes          , only : treeNode
    use :: Mass_Distributions        , only : massDistributionClass
    !![
    <include directive="rotationCurveTask" type="moduleUse">
    !!]
    include 'galactic_structure.rotation_curve.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureStandard   ), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                              , intent(in   )           :: radius
    class           (massDistributionClass       ), pointer                 :: massDistribution_
    double precision                                                        :: velocityRotationComponent__, rotationCurveSquared

    call self%defaults(radius=radius,componentType=componentType,massType=massType)
    massDistribution_    => node             %massDistribution(                                                                                                                          )
    rotationCurveSquared =  massDistribution_%rotationCurve   (galacticStructureState_%state%radius_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)**2
    !![
    <objectDestructor name="massDistribution_"/>
    <include directive="rotationCurveTask" type="functionCall" functionType="function" returnParameter="velocityRotationComponent__">
     <functionArgs>node,galacticStructureState_%state%radius_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_</functionArgs>
     <onReturn>rotationCurveSquared=rotationCurveSquared+velocityRotationComponent__**2</onReturn>
    !!]
    include 'galactic_structure.rotation_curve.tasks.inc'
    !![
    </include>
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
    use :: Galacticus_Nodes  , only : treeNode
    use :: Mass_Distributions, only : massDistributionClass
    !![
    <include directive="rotationCurveGradientTask" type="moduleUse">
    !!]
    include 'galactic_structure.rotation_curve.gradient.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureStandard   ), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    double precision                              , intent(in   )           :: radius
    class           (massDistributionClass       )               , pointer  :: massDistribution_
    double precision                                                        :: velocityRotationGradientComponent__, velocityRotation

    call self%defaults(radius=radius,componentType=componentType,massType=massType)
    massDistribution_        => node             %massDistribution     (                                                                                                                          )
    velocityRotationGradient =  massDistribution_%rotationCurveGradient(galacticStructureState_%state%radius_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    !![
    <objectDestructor name="massDistribution_"/>
    <include directive="rotationCurveGradientTask" type="functionCall" functionType="function" returnParameter="velocityRotationGradientComponent__">
     <functionArgs>node,galacticStructureState_%state%radius_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_</functionArgs>
     <onReturn>velocityRotationGradient=velocityRotationGradient+velocityRotationGradientComponent__</onReturn>
    !!]
    include 'galactic_structure.rotation_curve.gradient.tasks.inc'
    !![
    </include>
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
    !![
    <include directive="potentialTask" type="moduleUse">
    !!]
    include 'galactic_structure.potential.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureStandard        ), intent(inout)           :: self
    type            (treeNode                         ), intent(inout)           :: node
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType
    double precision                                   , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    class           (massDistributionClass            ), pointer                 :: massDistribution_
    type            (coordinateCylindrical            )                          :: position
    double precision                                                             :: potentialComponent__

    ! Initialize status.
    if (present(status)) status=structureErrorCodeSuccess
    ! Reset calculations if this is a new node.
    if (node%uniqueID() /= self%uniqueIDPrevious) call self%calculationReset(node)
    ! Get the mass distribution for this calculation.
    massDistribution_ => node%massDistribution()
    ! Evaluate the potential at the halo virial radius.
    if (.not.self%potentialOffsetComputed) then
       call self%defaults(componentType=componentTypeAll,massType=massTypeAll,radius=self%darkMatterHaloScale_%radiusVirial(node))
       position =[galacticStructureState_%state%radius_,0.0d0,0.0d0]
       potential=massDistribution_%potential(position,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,status_)
       if (status_ /= structureErrorCodeSuccess) status=status_
       !![
       <include directive="potentialTask" type="functionCall" functionType="function" returnParameter="potentialComponent__">
        <functionArgs>node,galacticStructureState_%state%radius_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,status</functionArgs>
        <onReturn>potential=potential+potentialComponent__</onReturn>
       !!]
       include 'galactic_structure.potential.tasks.inc'
       !![
       </include>
       !!]
       call self%restore()
       ! Compute the potential offset such that the total gravitational potential at the virial radius is -V² where V is the
       ! virial velocity.
       self%potentialOffset        =-potential-self%darkMatterHaloScale_%velocityVirial(node)**2
       self%potentialOffsetComputed=.true.
    end if
    call self%defaults(radius=radius,componentType=componentType,massType=massType)
    ! Determine which component type to use.
    if (present(componentType)) then
       galacticStructureState_%state%componentType_=componentType
    else
       galacticStructureState_%state%componentType_=componentTypeAll
    end if
    position =[galacticStructureState_%state%radius_,0.0d0,0.0d0]
    potential=+massDistribution_%potential(position,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,status_) &
         &    +self%potentialOffset
    if (status_ /= structureErrorCodeSuccess) status=status_
    include 'galactic_structure.potential.tasks.inc'
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
    use :: Coordinate_Systems        , only : Coordinates_Cartesian_To_Cylindrical, Coordinates_Spherical_To_Cylindrical
    use :: Galactic_Structure_Options, only : componentTypeAll                    , coordinateSystemCartesian           , coordinateSystemCylindrical, coordinateSystemSpherical      , &
          &                                   massTypeAll                         , weightByMass                        , weightIndexNull            , enumerationCoordinateSystemType
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : optimizeForSurfaceDensitySummation  , optimizeforsurfacedensitysummation  , reductionSummation         , reductionsummation             , &
          &                                   treeNode
    implicit none
    class           (galacticStructureStandard      ), intent(inout)               :: self
    type            (treeNode                       ), intent(inout)               :: node
    type            (enumerationComponentTypeType   ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional     :: massType
    type            (enumerationWeightByType        ), intent(in   ), optional     :: weightBy
    integer                                          , intent(in   ), optional     :: weightIndex
    type            (enumerationCoordinateSystemType), intent(in   ), optional     :: coordinateSystem
    double precision                                 , intent(in   ), dimension(3) :: position
    procedure       (surfaceDensityComponent        ), pointer                     :: surfaceDensityComponent_
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
    ! Call routines to supply the densities for all components.
    surfaceDensityComponent_ => surfaceDensityComponent
    surfaceDensity           =  node%mapDouble0(surfaceDensityComponent_,reductionSummation,optimizeFor=optimizeForSurfaceDensitySummation)
    call self%restore()
    return
  end function standardSurfaceDensity

  double precision function surfaceDensityComponent(component)
    !!{
    Unary function returning the surface density in a component. Suitable for mapping over components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    class(nodeComponent), intent(inout) :: component

    surfaceDensityComponent=component%surfaceDensity(positionCylindrical_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_,galacticStructureState_%state%weightBy_,galacticStructureState_%state%weightIndex_)
    return
  end function surfaceDensityComponent

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
    !![
    <include directive="accelerationTask" type="moduleUse">
    !!]
    include 'galactic_structure.acceleration.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureStandard   ), intent(inout)               :: self
    double precision                                             , dimension(3) :: acceleration
    type            (treeNode                    ), intent(inout)               :: node
    double precision                              , intent(in   ), dimension(3) :: positionCartesian
    type            (enumerationComponentTypeType), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional     :: massType
    double precision                                             , dimension(3) :: accelerationComponent__
    class           (massDistributionClass       ), pointer                     :: massDistribution_
    type            (coordinateCartesian         )                              :: position

    call self%defaults(componentType=componentType,massType=massType)
    positionCartesian_ =  positionCartesian
    position           =  positionCartesian
    massDistribution_  => node             %massDistribution(                                                                                             )
    acceleration       =  massDistribution_%acceleration    (position,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    !![
    <objectDestructor name="massDistribution_"/>
    <include directive="accelerationTask" type="functionCall" functionType="function" returnParameter="accelerationComponent__">
     <functionArgs>node,positionCartesian_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_</functionArgs>
     <onReturn>acceleration=acceleration+accelerationComponent__</onReturn>
    !!]
    include 'galactic_structure.acceleration.tasks.inc'
    !![
    </include>
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
    !![
    <include directive="tidalTensorTask" type="moduleUse">
    !!]
    include 'galactic_structure.tidal_tensor.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureStandard     ), intent(inout)               :: self
    type            (tensorRank2Dimension3Symmetric)                              :: tidalTensor
    type            (treeNode                      ), intent(inout)               :: node
    double precision                                , intent(in   ), dimension(3) :: positionCartesian
    type            (enumerationComponentTypeType  ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType       ), intent(in   ), optional     :: massType
    type            (tensorRank2Dimension3Symmetric)                              :: tidalTensorComponent__
    class           (massDistributionClass         ), pointer                     :: massDistribution_
    type            (coordinateCartesian           )                              :: position
    
    call self%defaults(componentType=componentType,massType=massType)
    positionCartesian_ =  positionCartesian
    position           =  positionCartesian
    massDistribution_  => node             %massDistribution(                                                                                             )
    tidalTensor        =  massDistribution_%tidalTensor     (position,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    !![
    <objectDestructor name="massDistribution_"/>
    <include directive="tidalTensorTask" type="functionCall" functionType="function" returnParameter="tidalTensorComponent__">
     <functionArgs>node,positionCartesian_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_</functionArgs>
     <onReturn>tidalTensor=tidalTensor+tidalTensorComponent__</onReturn>
    !!]
    include 'galactic_structure.tidal_tensor.tasks.inc'
    !![
    </include>
    !!]
    call self%restore()
    return
  end function standardTidalTensor

  function standardChandrasekharIntegral(self,node,nodeSatellite,positionCartesian,velocityCartesian,componentType,massType) result(chandrasekharIntegral)
    !!{
    Compute the integral appearing in the \cite{chandrasekhar_dynamical_1943} dynamical friction model:
    \begin{equation}
      \rho(\boldsymbol{x}_\mathrm{s}) \int \mathrm{d}\boldsymbol{v} f(\boldsymbol{v}) {\boldsymbol{v}-\boldsymbol{v}_\mathrm{s} \over |\boldsymbol{v}-\boldsymbol{v}_\mathrm{s}|^3},
    \end{equation}  
    where $\rho(\boldsymbol{x}_\mathrm{s})$ is the density at the position of the perturber, $\boldsymbol{x}_\mathrm{s}$,
    $f(\boldsymbol{v})$ is the velocity distribution function at velocity $\boldsymbol{v}$, and $\boldsymbol{v}_\mathrm{s}$ is
    the velocity of the perturber.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll                         , massTypeAll
    use :: Galacticus_Nodes          , only : optimizeForChandrasekharIntegralSummation, reductionSummation, treeNode
    !![
    <include directive="chandrasekharIntegralTask" type="moduleUse">
    !!]
    include 'galactic_structure.chandrasekharIntegral.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    double precision                                               , dimension(3) :: chandrasekharIntegral
    class           (galacticStructureStandard     ), intent(inout)               :: self
    type            (treeNode                      ), intent(inout), target       :: node                              , nodeSatellite
    double precision                                , intent(in   ), dimension(3) :: positionCartesian                 , velocityCartesian
    type            (enumerationComponentTypeType  ), intent(in   ), optional     :: componentType
    type            (enumerationMassTypeType       ), intent(in   ), optional     :: massType
    integer                                         , parameter                   :: chandrasekharIntegralSize       =3
    procedure       (chandrasekharIntegralComponent), pointer                     :: chandrasekharIntegralComponent_
    double precision                                               , dimension(3) :: chandrasekharIntegralComponent__

    call self%defaults(componentType=componentType,massType=massType)
    positionCartesian_              =  positionCartesian
    velocityCartesian_              =  velocityCartesian
    nodeSatellite_                  => nodeSatellite
    chandrasekharIntegralComponent_ => chandrasekharIntegralComponent
    chandrasekharIntegral           =  node%mapDouble1(chandrasekharIntegralComponent_,chandrasekharIntegralSize,reductionSummation,optimizeFor=optimizeForChandrasekharIntegralSummation)
    !![
    <include directive="chandrasekharIntegralTask" type="functionCall" functionType="function" returnParameter="chandrasekharIntegralComponent__">
     <functionArgs>node,nodeSatellite,positionCartesian_,velocityCartesian_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_</functionArgs>
     <onReturn>chandrasekharIntegral=chandrasekharIntegral+chandrasekharIntegralComponent__</onReturn>
    !!]
    include 'galactic_structure.chandrasekharIntegral.tasks.inc'
    !![
    </include>
    !!]
    call self%restore()
    return
  end function standardChandrasekharIntegral

  function chandrasekharIntegralComponent(component,resultSize)
    !!{
    Function returning the Chandrasekhar integral in a component. Suitable for mapping over components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    integer                        , intent(in   )         :: resultSize
    class           (nodeComponent), intent(inout)         :: component
    double precision               , dimension(resultSize) :: chandrasekharIntegralComponent

    chandrasekharIntegralComponent=component%chandrasekharIntegral(nodeSatellite_,positionCartesian_,velocityCartesian_,galacticStructureState_%state%componentType_,galacticStructureState_%state%massType_)
    return
  end function chandrasekharIntegralComponent

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

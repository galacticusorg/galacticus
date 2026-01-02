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

  !!{
  Implementation of a ``fixed'' solver for galactic structure (no self-gravity of baryons, and size simply scales in
  proportion to specific angular momentum).
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass
  
  !![
  <enumeration>
   <name>radiusFixed</name>
   <description>Enumerates the possible definitions of radius used by the ``fixed'' galactic structure solver.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="virial"    />
   <entry label="turnaround"/>
  </enumeration>
  !!]

  !![
  <galacticStructureSolver name="galacticStructureSolverFixed">
   <description>
    A galactic structure solver that determines the sizes of galactic components by assuming that radius equals \begin{equation} r
    = f_\mathrm{r} \lambda r_0 \end{equation} where $r_0$ is the virial or turnaround radius of the \gls{node} if {\normalfont
    \ttfamily [radiusFixed]}$=${\normalfont \ttfamily virialRadius} or {\normalfont \ttfamily turnaround} respectively, $\lambda$
    is its spin parameter and $f_\mathrm{r}=${\normalfont \ttfamily [factor]} is a parameter. Optionally, different values of
    $f_\mathrm{r}$ can be specified for disks and spheroids using the {\normalfont \ttfamily [factorDisk]} and {\normalfont
    \ttfamily [factorSpheroid]} parameters respectively---if either or both are not provided the value of {\normalfont \ttfamily
    [factor]} will be used for the corresponding component.
   </description>
  </galacticStructureSolver>
  !!]
  type, extends(galacticStructureSolverClass) :: galacticStructureSolverFixed
     !!{
     Implementation of a ``fixed'' solver for galactic structure (no self-gravity of baryons, and size simply scales in
     proportion to specific angular momentum).
     !!}
     private
     double precision                                      :: factor                          , factorDisk        , &
          &                                                   factorSpheroid
     type            (enumerationRadiusFixedType)          :: radiusFixed
     class           (darkMatterHaloScaleClass  ), pointer :: darkMatterHaloScale_   => null()
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
   contains
     final     ::             fixedDestructor
     procedure :: solve    => fixedSolve
     procedure :: revert   => fixedRevert
     procedure :: autoHook => fixedAutoHook
  end type galacticStructureSolverFixed

  interface galacticStructureSolverFixed
     !!{
     Constructors for the \refClass{galacticStructureSolverFixed} galactic structure solver class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface galacticStructureSolverFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticStructureSolverFixed} galactic structure solver class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticStructureSolverFixed)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_
    class           (virialDensityContrastClass  ), pointer       :: virialDensityContrast_
    double precision                                              :: factor                , factorDisk, &
         &                                                           factorSpheroid
    type            (varying_string              )                :: radiusFixed

    !![
    <inputParameter>
      <name>factor</name>
      <defaultSource>\citep{mo_formation_1998}</defaultSource>
      <defaultValue>sqrt(0.5d0)</defaultValue>
      <description>The ratio of galaxy radius to $\lambda r_\mathrm{vir}$ in the ``fixed'' galactic structure radius solver algorithm. This will be applied to any component for which no component-specific value is provided.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>factorDisk</name>
      <defaultSource>\citep{mo_formation_1998}</defaultSource>
      <defaultValue>sqrt(0.5d0)</defaultValue>
      <description>The ratio of galaxy radius to $\lambda r_\mathrm{vir}$ in the ``fixed'' galactic structure radius solver algorithm for disks. This will override the generic value supplied by {\normalfont \ttfamily [factor]} for disks.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>factorSpheroid</name>
      <defaultSource>\citep{mo_formation_1998}</defaultSource>
      <defaultValue>sqrt(0.5d0)</defaultValue>
      <description>The ratio of galaxy radius to $\lambda r_\mathrm{vir}$ in the ``fixed'' galactic structure radius solver algorithm for spheroids. This will override the generic value supplied by {\normalfont \ttfamily [factor]} for spheroids.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusFixed</name>
      <defaultValue>var_str('virial')</defaultValue>
      <description>The radius to use in the ``fixed'' galactic structure radius solver algorithm. Allowed options are ``virial'' and ``turnaround''.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    <conditionalCall>
      <call>self=galacticStructureSolverFixed(darkMatterHaloScale_,virialDensityContrast_,enumerationRadiusFixedEncode(char(radiusFixed),includesPrefix=.false.),factor{conditions})</call>
      <argument name="factorDisk"     value="factorDisk"     parameterPresent="parameters"/>
      <argument name="factorSpheroid" value="factorSpheroid" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(darkMatterHaloScale_,virialDensityContrast_,radiusFixed,factor,factorDisk,factorSpheroid) result(self)
    !!{
    Internal constructor for the \refClass{galacticStructureSolverFixed} galactic structure solver class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (galacticStructureSolverFixed)                          :: self
    class           (darkMatterHaloScaleClass    ), intent(in   ), target   :: darkMatterHaloScale_
    class           (virialDensityContrastClass  ), intent(in   ), target   :: virialDensityContrast_
    type            (enumerationRadiusFixedType  ), intent(in   )           :: radiusFixed
    double precision                              , intent(in   )           :: factor
    double precision                              , intent(in   ), optional :: factorDisk            , factorSpheroid
    !![
    <constructorAssign variables="factor, radiusFixed, *darkMatterHaloScale_, *virialDensityContrast_, factorDisk, factorSpheroid"/>
    !!]

    if (.not.enumerationRadiusFixedIsValid(radiusFixed)) call Error_Report('invalid radiusFixed'//{introspection:location})
    if (.not.present(factorDisk    )) self%factorDisk    =factor
    if (.not.present(factorSpheroid)) self%factorSpheroid=factor
    return
  end function fixedConstructorInternal

  subroutine fixedAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent  , openMPThreadBindingAtLevel, postEvolveEvent, preDerivativeEvent, &
          &                     satelliteMergerEvent, dependencyDirectionAfter  , dependencyRegEx
    implicit none
    class(galacticStructureSolverFixed), intent(inout) :: self
    type (dependencyRegEx             ), dimension(1)  :: dependencies

    dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^nodeComponent')
    call   preDerivativeEvent%attach(self,fixedSolvePreDeriativeHook,openMPThreadBindingAtLevel,label='structureSolverFixed',dependencies=dependencies)
    call      postEvolveEvent%attach(self,fixedSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverFixed',dependencies=dependencies)
    call satelliteMergerEvent%attach(self,fixedSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverFixed',dependencies=dependencies)
    call   nodePromotionEvent%attach(self,fixedSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverFixed',dependencies=dependencies)
    return
  end subroutine fixedAutoHook

  subroutine fixedDestructor(self)
    !!{
    Destructor for the \refClass{galacticStructureSolverFixed} galactic structure solver class.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent, postEvolveEvent, preDerivativeEvent, satelliteMergerEvent
    implicit none
    type(galacticStructureSolverFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"  />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    if (  preDerivativeEvent%isAttached(self,fixedSolvePreDeriativeHook)) call   preDerivativeEvent%detach(self,fixedSolvePreDeriativeHook)
    if (     postEvolveEvent%isAttached(self,fixedSolveHook            )) call      postEvolveEvent%detach(self,fixedSolveHook            )
    if (satelliteMergerEvent%isAttached(self,fixedSolveHook            )) call satelliteMergerEvent%detach(self,fixedSolveHook            )
    if (  nodePromotionEvent%isAttached(self,fixedSolveHook            )) call   nodePromotionEvent%detach(self,fixedSolveHook            )
    return
  end subroutine fixedDestructor

  subroutine fixedSolveHook(self,node)
    !!{
    Hookable wrapper around the solver.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node

    select type (self)
    type is (galacticStructureSolverFixed)
       call self%solve(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine fixedSolveHook

  subroutine fixedSolvePreDeriativeHook(self,node,propertyType)
    !!{
    Hookable wrapper around the solver.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (*       ), intent(inout)         :: self
    type   (treeNode), intent(inout), target :: node
    integer          , intent(in   )         :: propertyType
    !$GLC attributes unused :: propertyType

    select type (self)
    type is (galacticStructureSolverFixed)
       call self%solve(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine fixedSolvePreDeriativeHook

  subroutine fixedSolve(self,node,plausibilityOnly)
    !!{
    Solve for the structure of galactic components assuming no self-gravity of baryons, and that size simply scales in
    proportion to specific angular momentum.
    !!}
    use :: Calculations_Resets       , only : Calculations_Reset
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    implicit none
    class           (galacticStructureSolverFixed), intent(inout)           :: self
    type            (treeNode                    ), intent(inout), target   :: node
    logical                                       , intent(in   ), optional :: plausibilityOnly
    logical                                       , parameter               :: specificAngularMomentumRequired=.false.
    procedure       (solverGet                   ), pointer                 :: radiusGet                              , velocityGet
    procedure       (solverSet                   ), pointer                 :: radiusSet                              , velocitySet
    logical                                                                 :: componentActive
    double precision                                                        :: specificAngularMomentum
    type            (enumerationComponentTypeType)                          :: component
    !![
    <optionalArgument name="plausibilityOnly" defaultsTo=".false."/>
    !!]

    ! Check that the galaxy is physical plausible. In this fixed solver, we don't act on this.
    node%isPhysicallyPlausible=.true.
    node%isSolvable           =.true.
    include 'galactic_structure.radius_solver.plausible.inc'
    if (.not.node%isPhysicallyPlausible .or. plausibilityOnly_) return
    call Calculations_Reset(node)
    include 'galactic_structure.radius_solver.tasks.inc'
    return

  contains

    subroutine radiusSolve(node,component,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet)
      !!{
      Solve for the equilibrium radius of the given component.
      !!}
      use :: Dark_Matter_Halo_Spins    , only : Dark_Matter_Halo_Angular_Momentum_Scale
      use :: Galacticus_Nodes          , only : nodeComponentBasic                     , nodeComponentSpin, treeNode
      use :: Mass_Distributions        , only : massDistributionClass
      use :: Galactic_Structure_Options, only : componentTypeDarkMatterOnly            , massTypeDark     , componentTypeDisk, componentTypeSpheroid
      implicit none
      type            (treeNode                    ), intent(inout)          :: node
      type            (enumerationComponentTypeType), intent(in   )          :: component
      double precision                              , intent(in   )          :: specificAngularMomentum
      procedure       (solverGet                   ), intent(in   ), pointer :: radiusGet              , velocityGet
      procedure       (solverSet                   ), intent(in   ), pointer :: radiusSet              , velocitySet
      class           (nodeComponentSpin           )               , pointer :: spin
      class           (nodeComponentBasic          )               , pointer :: basic
      class           (massDistributionClass       )               , pointer :: massDistribution_
      double precision                                                       :: radius                 , velocity   , &
           &                                                                    factor
      !$GLC attributes unused :: radiusGet, velocityGet, specificAngularMomentum

      ! Determine the factor to use.
      select case (component%ID)
      case (componentTypeDisk    %ID)
         factor=self%factorDisk
      case (componentTypeSpheroid%ID)
         factor=self%factorSpheroid
      case default
         factor=self%factor
      end select
      ! Find the radius of the component, assuming radius is a fixed fraction of radius times spin parameter.
      spin => node%spin()
      select case (self%radiusFixed%ID)
      case (radiusFixedVirial    %ID)
         velocity          =  +self             %darkMatterHaloScale_  %velocityVirial              (     node                                                          )
         radius            =  +self             %darkMatterHaloScale_  %radiusVirial                (     node                                                          ) &
              &               *                                         factor                                                                                            &
              &               *spin                                    %angularMomentum             (                                                                   ) &
              &               /Dark_Matter_Halo_Angular_Momentum_Scale                              (     node                       ,     self %darkMatterHaloScale_   )
      case (radiusFixedTurnaround%ID)
         massDistribution_ =>  node                                    %massDistribution            (     componentTypeDarkMatterOnly,           massTypeDark           )
         basic             =>  node                                    %basic                       (                                                                   )
         velocity          =  +massDistribution_                       %velocityRotationCurveMaximum(                                                                   )
         radius            =  +self             %darkMatterHaloScale_  %radiusVirial                (     node                                                          ) &
              &               *self             %virialDensityContrast_%turnAroundOverVirialRadii   (mass=basic%mass()               ,time=basic%timeLastIsolated     ()) &
              &               *                                         factor                                                                                            &
              &               *spin                                    %angularMomentum             (                                                                   ) &
              &               /Dark_Matter_Halo_Angular_Momentum_Scale                              (     node                       ,     self %darkMatterHaloScale_   )
         !![
	 <objectDestructor name="massDistribution_"/>
         !!]
      end select
      ! Set the component size to new radius and velocity.
      call radiusSet  (node,radius  )
      call velocitySet(node,velocity)
      return
    end subroutine radiusSolve

  end subroutine fixedSolve

  subroutine fixedRevert(self,node)
    !!{
    Revert radii for the fixed galactic structure solve. Not necessary for this algorithm.
    !!}
    implicit none
    class(galacticStructureSolverFixed), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    return
  end subroutine fixedRevert

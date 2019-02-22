!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a galactic structure initial radius class using the adiabatic contraction algorithm of
  !% \cite{gnedin_response_2004}.

  use Kind_Numbers
  use Galactic_Structure_Options
  use Math_Exponentiation
  use Cosmology_Parameters
  use Dark_Matter_Profiles
  use Dark_Matter_Halo_Scales

  ! Number of previous radius solutions to store.
  integer, parameter :: gnedin2004StoreCount=10

  !# <galacticStructureRadiiInitial name="galacticStructureRadiiInitialGnedin2004">
  !#  <description>A galactic structure initial radius class using the adiabatic contraction algorithm of \cite{gnedin_response_2004}.</description>
  !# </galacticStructureRadiiInitial>
  type, extends(galacticStructureRadiiInitialClass) :: galacticStructureRadiiInitialGnedin2004
     !% A galactic structure initial radius class using the adiabatic contraction algorithm of \cite{gnedin_response_2004}.
     private
     class           (cosmologyParametersClass), pointer                         :: cosmologyParameters_ => null()
     class           (darkMatterHaloScaleClass), pointer                         :: darkMatterHaloScale_ => null()
     class           (darkMatterProfileClass  ), pointer                         :: darkMatterProfile_ => null()
     ! Parameters of the adiabatic contraction algorithm.
     double precision                                                            :: A                              , omega
     ! Stored solutions for reuse.
     integer         (kind=kind_int8          )                                  :: uniqueIDPrevious
     integer                                                                     :: radiusPreviousIndex            , radiusPreviousIndexMaximum
     double precision                          , dimension(gnedin2004StoreCount) :: radiusPrevious                 , radiusInitialPrevious
     type            (fastExponentiator       )                                  :: radiusExponentiator
     ! Quantities used in solving the initial radius root function.
     double precision                                                            :: baryonicFinalTerm              , baryonicFinalTermDerivative, &
          &                                                                         darkMatterFraction             , initialMassFraction        , &
          &                                                                         radiusFinal                    , radiusFinalMean            , &
          &                                                                         radiusFinalMeanSelfDerivative  , radiusInitial              , &
          &                                                                         radiusInitialMeanSelfDerivative, radiusShared               , &
          &                                                                         radiusVirial
     logical                                                                     :: massesComputed
   contains
     !@ <objectMethods>
     !@   <object>galacticStructureRadiiInitialGnedin2004</object>
     !@   <objectMethod>
     !@     <method>computeFactors</method>
     !@     <description>Compute factors needed for solving adiabatic contraction.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} node\arginout, \doublezero\ radius\argin, \logicalzero\ [computeGradientFactors]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusOrbitalMean</method>
     !@     <description>Compute the orbit-averaged radius for dark matter.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>radiusOrbitalMeanDerivative</method>
     !@     <description>Compute the derivative of the orbit-averaged radius for dark matter.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin</arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                gnedin2004Destructor
     procedure :: radius                      => gnedin2004Radius
     procedure :: radiusDerivative            => gnedin2004RadiusDerivative
     procedure :: calculationReset            => gnedin2004CalculationReset
     procedure :: computeFactors              => gnedin2004ComputeFactors
     procedure :: radiusOrbitalMean           => gnedin2004RadiusOrbitalMean
     procedure :: radiusOrbitalMeanDerivative => gnedin2004RadiusOrbitalMeanDerivative
  end type galacticStructureRadiiInitialGnedin2004

  interface galacticStructureRadiiInitialGnedin2004
     !% Constructors for the {\normalfont \ttfamily gnedin2004} galactic structure initial radius class.
     module procedure gnedin2004ConstructorParameters
     module procedure gnedin2004ConstructorInternal
  end interface galacticStructureRadiiInitialGnedin2004

  ! Module-scope quantities used in solving the initial radius root function.
  integer                                         , parameter :: gnedin2004ComponentType=componentTypeAll, gnedin2004MassType   =massTypeBaryonic, &
       &                                                         gnedin2004WeightBy     =weightByMass    , gnedin2004WeightIndex=weightIndexNull
  logical                                         , parameter :: gnedin2004HaloLoaded   =.false.
  type   (treeNode                               ), pointer   :: gnedin2004Node
  class  (galacticStructureRadiiInitialGnedin2004), pointer   :: gnedin2004Self
  !$omp threadprivate(gnedin2004Self,gnedin2004Node)
  
contains

  function gnedin2004ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily gnedin2004} galactic structure initial radius class which takes a parameter list as
    !% input.
    use Input_Parameters
    implicit none
    type            (galacticStructureRadiiInitialGnedin2004)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileClass                 ), pointer       :: darkMatterProfile_
    double precision                                                         :: A                   , omega

    !# <inputParameter>
    !#   <name>A</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultSource>(\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultSource>
    !#   <defaultValue>0.80d0</defaultValue>
    !#   <description>The parameter $A$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>omega</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultSource>(\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultSource>
    !#   <defaultValue>0.77d0</defaultValue>
    !#   <description>The parameter $\omega$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !# <objectBuilder class="darkMatterProfile"   name="darkMatterProfile_"   source="parameters"/>
    self=galacticStructureRadiiInitialGnedin2004(A,omega,cosmologyParameters_,darkMatterHaloScale_,darkMatterProfile_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    !# <objectDestructor name="darkMatterProfile_"  />
    return
  end function gnedin2004ConstructorParameters

  function gnedin2004ConstructorInternal(A,omega,cosmologyParameters_,darkMatterHaloScale_,darkMatterProfile_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily gnedin2004} galactic structure initial radius class
    implicit none
    type            (galacticStructureRadiiInitialGnedin2004)                        :: self
    double precision                                         , intent(in   )         :: A                   , omega
    class           (cosmologyParametersClass               ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileClass                 ), intent(in   ), target :: darkMatterProfile_
    !# <constructorAssign variables="A, omega, *cosmologyParameters_, *darkMatterHaloScale_, *darkMatterProfile_"/>
    
    self%uniqueIDPrevious   =-1_kind_int8
    self%massesComputed=.false.
    self%radiusExponentiator=fastExponentiator(1.0d-3,1.0d0,omega,1.0d4,.false.)
    return
  end function gnedin2004ConstructorInternal

  subroutine gnedin2004Destructor(self)
    !% Destructor for the {\normalfont \ttfamily gnedin2004} galactic structure initial radius class.
    implicit none
    type(galacticStructureRadiiInitialGnedin2004), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    !# <objectDestructor name="self%darkMatterProfile_"  />
    return
  end subroutine gnedin2004Destructor

  double precision function gnedin2004Radius(self,node,radius)
    !% Compute the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    !% \cite{gnedin_response_2004}.
    use Root_Finder
    implicit none
    class           (galacticStructureRadiiInitialGnedin2004), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    double precision                                         , intent(in   ) :: radius
    double precision                                         , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-2
    type            (rootFinder                             ), save          :: finder
    !$omp threadprivate(finder)
    integer                                                                  :: i                      , j                       , &
         &                                                                      iMod
    double precision                                                         :: radiusUpperBound
    
    ! Reset stored solutions if the node has changed.
    if (node%uniqueID() /= self%uniqueIDPrevious) call self%calculationReset(node)
    ! Check for a previously computed solution.
    if (self%radiusPreviousIndexMaximum > 0 .and. any(self%radiusPrevious(1:self%radiusPreviousIndexMaximum) == radius)) then
       gnedin2004Radius=0.0d0
       do i=1,self%radiusPreviousIndexMaximum
          if (self%radiusPrevious(i) == radius) then
             gnedin2004Radius=self%radiusInitialPrevious(i)
             exit
          end if
       end do
    else
       ! Get the virial radius of the node.
       self%radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
       ! Return radius unchanged if larger than the virial radius.
       if (radius >= self%radiusVirial) then
          gnedin2004Radius=radius
       else
          ! Compute the various factors needed by this calculation.
          call self%computeFactors(node,radius,computeGradientFactors=.false.)
          if (gnedin2004Solver(self%radiusVirial) < 0.0d0) then ! Check that solution is within bounds.
             gnedin2004Radius=self%radiusVirial
          else
             j=-1
             if (self%radiusPreviousIndexMaximum > 0) then
                ! No exact match exists, look for approximate matches.
                do i=1,self%radiusPreviousIndexMaximum
                   iMod=modulo(self%radiusPreviousIndex-i,gnedin2004StoreCount)+1
                   if (abs(radius-self%radiusPrevious(iMod))/self%radiusPrevious(iMod) < toleranceRelative) then
                      j=iMod
                      exit
                   end if
                end do
             end if
             ! Initialize our root finder.
             if (.not.finder%isInitialized()) then
                call finder%rootFunction(gnedin2004Solver                   )
                call finder%tolerance   (toleranceAbsolute,toleranceRelative)
             end if
             ! Find the solution for initial radius.
             if (j == -1) then
                ! No previous solution to use as an initial guess. Instead, we make an estimate of the initial radius under the
                ! assumption that the mass of dark matter (in the initial profile) enclosed within the mean initial radius is the
                ! same as enclosed within the mean final radius. Since the initial and final radii are typically not too
                ! different, and since the mean radius is a weak (ѡ<1) function of the radius this is a useful
                ! approximation. Furthermore, since it will underestimate the actual mass within the initial mean radius it gives
                ! an overestimate of the initial radius. This means that we have a bracketing of the initial radius which we can
                ! use in the solver.
                radiusUpperBound   =  +(                                                                                     &
                     &                  +self%baryonicFinalTerm                                                              &
                     &                  /self%darkMatterProfile_%enclosedMass(node,self%radiusOrbitalMean(self%radiusFinal)) &
                     &                  +self%darkMatterFraction                                                             &
                     &                  *self%radiusFinal                                                                    &
                     &                 )                                                                                     &
                     &                /  self%initialMassFraction
                if (radiusUpperBound < radius) radiusUpperBound=radius
                call finder%rangeExpand(                                                             &
                     &                  rangeExpandUpward            =1.1d0                        , &
                     &                  rangeExpandDownward          =0.9d0                        , &
                     &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
                     &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
                     &                  rangeExpandType              =rangeExpandMultiplicative      &
                     &                 )
                gnedin2004Radius=finder%find(rootRange=[radius,radiusUpperBound])
             else
               ! Use previous solution as an initial guess.
               call finder%rangeExpand(                                                                   &
                    &                  rangeExpandDownward          =1.0d0/sqrt(1.0d0+toleranceRelative), &
                    &                  rangeExpandUpward            =1.0d0*sqrt(1.0d0+toleranceRelative), &
                    &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative      , &
                    &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive      , &
                    &                  rangeExpandType              =rangeExpandMultiplicative            &
                    &                 )
               gnedin2004Radius=finder%find(                                                                        &
                    &                       rootRange=[                                                             &
                    &                                  self%radiusInitialPrevious(j)/sqrt(1.0d0+toleranceRelative), &
                    &                                  self%radiusInitialPrevious(j)*sqrt(1.0d0+toleranceRelative)  &
                    &                                 ]                                                             &
                    &                      )
            end if
          end if
       end if
       ! Store this solution.       
       self%radiusPreviousIndex                                 =modulo(self%radiusPreviousIndex         ,gnedin2004StoreCount)+1
       self%radiusPreviousIndexMaximum                          =min   (self%radiusPreviousIndexMaximum+1,gnedin2004StoreCount)
       self%radiusPrevious            (self%radiusPreviousIndex)=radius
       self%radiusInitialPrevious     (self%radiusPreviousIndex)=gnedin2004Radius
    end if
    return
  end function gnedin2004Radius

  double precision function gnedin2004RadiusDerivative(self,node,radius)
    !% Compute the derivative of the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    !% \cite{gnedin_response_2004}.
    use Root_Finder
    implicit none
    class           (galacticStructureRadiiInitialGnedin2004), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    double precision                                         , intent(in   ) :: radius
    double precision                                         , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    type            (rootFinder                             ), save          :: finder
    !$omp threadprivate(finder)
    
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rangeExpand (                                                             &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeDownwardLimit           =0.0d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
       call finder%rootFunction(gnedin2004DerivativeSolver                  )
       call finder%tolerance   (toleranceAbsolute         ,toleranceRelative)
    end if
    ! Reset stored solutions if the node has changed.
    if (node%uniqueID() /= self%uniqueIDPrevious) call self%calculationReset(node)
    ! Compute the various factors needed by this calculation.
    call self%computeFactors(node,radius,computeGradientFactors=.true.)
    ! Return unit derivative if radius is larger than the virial radius.
    if (radius >= self%radiusVirial) then
       gnedin2004RadiusDerivative=1.0d0
       return
    end if
    ! Compute initial radius, and derivatives of initial and final mean radii.
    self%radiusFinal                    =                                           radius
    self%radiusInitial                  =self%radius                     (node,     radius       )
    self%radiusInitialMeanSelfDerivative=self%radiusOrbitalMeanDerivative(     self%radiusInitial)
    self%radiusFinalMeanSelfDerivative  =self%radiusOrbitalMeanDerivative(          radius       )
    ! Find the solution for initial radius.
    gnedin2004RadiusDerivative=finder%find(rootGuess=1.0d0)
    return
  end function gnedin2004RadiusDerivative

  subroutine gnedin2004ComputeFactors(self,node,radius,computeGradientFactors)
    !% Compute various factors needed when solving for the initial radius in the dark matter halo using the adiabatic contraction
    !% algorithm of \cite{gnedin_response_2004}.
    use Numerical_Constants_Physical
    use Galacticus_Nodes            , only : treeNode          , nodeComponentBasic               , optimizeForRotationCurveGradientSummation, optimizeForEnclosedMassSummation, &
         &                                   reductionSummation, optimizeForRotationCurveSummation
    !# <include directive="rotationCurveTask" name="radiusSolverRotationCurveTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'galactic_structure.radius_solver.initial_radii.adiabatic.rotation_curve.tasks.modules.inc'
    !# </include>
    !# <include directive="rotationCurveGradientTask" name="radiusSolverRotationCurveGradientTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'galactic_structure.radius_solver.initial_radii.adiabatic.rotation_curve_gradient.tasks.modules.inc'
    !# </include>
    !# <include directive="enclosedMassTask" name="radiusSolverEnclosedMassTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'galactic_structure.radius_solver.initial_radii.adiabatic.enclosed_mass.tasks.modules.inc'
    !# </include>
    implicit none
    class           (galacticStructureRadiiInitialGnedin2004), intent(inout), target  :: self
    type            (treeNode                               ), intent(inout), target  :: node
    double precision                                         , intent(in   )          :: radius
    logical                                                  , intent(in   )          :: computeGradientFactors
    type            (treeNode                               )               , pointer :: nodeCurrent                           , nodeHost
    class           (nodeComponentBasic                     )               , pointer :: basic
    double precision                                         , parameter              :: toleranceAbsolute               =0.0d0, toleranceRelative   =1.0d-3
    procedure       (gnedin2004MassEnclosed                 )               , pointer :: mapFunction
    double precision                                                                  :: massBaryonicSelfTotal                 , massBaryonicTotal          , &
         &                                                                               massComponent                         , velocityComponent          , &
         &                                                                               velocityComponentSquaredGradient      , velocityCircularSquared    , &
         &                                                                               velocityCircularSquaredGradient
    
    ! Set module-scope pointers to node and self.
    gnedin2004Node => node
    gnedin2004Self => self
    ! Store the final radius and its orbit-averaged mean.
    self%radiusFinal    =                       radius
    self%radiusFinalMean=self%radiusOrbitalMean(radius)
    self%radiusShared   =self%radiusFinalMean
    ! Compute the baryonic contribution to the rotation curve.
    mapFunction             => gnedin2004VelocityCircularSquared
    velocityCircularSquared =  node%mapDouble0(mapFunction,reductionSummation,optimizeFor=optimizeForRotationCurveSummation)
    !# <include directive="rotationCurveTask" name="radiusSolverRotationCurveTask" type="functionCall" functionType="function" returnParameter="velocityComponent">
    !#  <exclude>Dark_Matter_Profile_Rotation_Curve_Task</exclude>
    !#  <functionArgs>node,self%radiusFinalMean,gnedin2004ComponentType,gnedin2004MassType,gnedin2004HaloLoaded</functionArgs>
    !#  <onReturn>velocityCircularSquared=velocityCircularSquared+velocityComponent**2</onReturn>
    include 'galactic_structure.radius_solver.initial_radii.adiabatic.rotation_curve.tasks.inc'
    !# </include>
    self%baryonicFinalTerm=velocityCircularSquared*self%radiusFinalMean*self%radiusFinal/gravitationalConstantGalacticus
    ! Compute the baryonic contribution to the rotation curve.
    if (computeGradientFactors) then
       mapFunction                     => gnedin2004VelocityCircularSquaredGradient
       velocityCircularSquaredGradient =  node%mapDouble0(mapFunction,reductionSummation,optimizeFor=optimizeForRotationCurveGradientSummation)
       !# <include directive="rotationCurveGradientTask" name="radiusSolverRotationCurveGradientTask" type="functionCall" functionType="function" returnParameter="velocityComponentSquaredGradient">
       !#  <exclude>Dark_Matter_Profile_Rotation_Curve_Gradient_Task</exclude>
       !#  <functionArgs>node,self%radiusFinalMean,gnedin2004ComponentType,gnedin2004MassType,gnedin2004HaloLoaded</functionArgs>
       !#  <onReturn>velocityCircularSquaredGradient=velocityCircularSquaredGradient+velocityComponentSquaredGradient</onReturn>
       include 'galactic_structure.radius_solver.initial_radii.adiabatic.rotation_curve_gradient.tasks.inc'
       !# </include>
       self%baryonicFinalTermDerivative=velocityCircularSquaredGradient*self%radiusOrbitalMeanDerivative(self%radiusFinal)*self%radiusFinalMean*self%radiusFinal/gravitationalConstantGalacticus
    end if
    ! Compute the initial baryonic contribution from this halo, and any satellites.
    if (.not.self%massesComputed) then
       massBaryonicTotal     =  0.0d0
       massBaryonicSelfTotal =  0.0d0
       nodeCurrent           => node
       nodeHost              => node
       do while (associated(nodeCurrent))
          mapFunction       => gnedin2004MassEnclosed
          massBaryonicTotal =  massBaryonicTotal+nodeCurrent%mapDouble0(mapFunction,reductionSummation,optimizeFor=optimizeForEnclosedMassSummation)
          !# <include directive="enclosedMassTask" name="radiusSolverEnclosedMassTask" type="functionCall" functionType="function" returnParameter="massComponent">
          !#  <exclude>Dark_Matter_Profile_Enclosed_Mass_Task</exclude>
          !#  <functionArgs>nodeCurrent,radiusLarge,gnedin2004ComponentType,gnedin2004MassType,gnedin2004WeightBy,gnedin2004WeightIndex,gnedin2004HaloLoaded</functionArgs>
          !#  <onReturn>massBaryonicTotal=massBaryonicTotal+massComponent</onReturn>
          include 'galactic_structure.radius_solver.initial_radii.adiabatic.enclosed_mass.tasks.inc'
          !# </include>
          if (associated(nodeCurrent,nodeHost)) then
             massBaryonicSelfTotal=massBaryonicTotal
             do while (associated(nodeCurrent%firstSatellite))
                nodeCurrent => nodeCurrent%firstSatellite
             end do
             if (associated(nodeCurrent,nodeHost)) nodeCurrent => null()
          else
             if (associated(nodeCurrent%sibling)) then
                nodeCurrent => nodeCurrent%sibling
                do while (associated(nodeCurrent%firstSatellite))
                   nodeCurrent => nodeCurrent%firstSatellite
                end do
             else
                nodeCurrent => nodeCurrent%parent
                if (associated(nodeCurrent,node)) nodeCurrent => null()
             end if
          end if
       end do
       ! Limit masses to physical values.
       massBaryonicSelfTotal=max(massBaryonicSelfTotal,0.0d0)
       massBaryonicTotal    =max(massBaryonicTotal    ,0.0d0)
       ! Compute the dark matter fraction.
       basic => node%basic()
       self%darkMatterFraction =min((self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter()+(massBaryonicTotal-massBaryonicSelfTotal)/basic%mass(),1.0d0)
       ! Compute the initial mass fraction.
       self%initialMassFraction=min((self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter()+ massBaryonicTotal                       /basic%mass(),1.0d0)
       ! Record that masses (and mass fractions) have been computed.
       self%massesComputed=.true.
    end if
    return
  end subroutine gnedin2004ComputeFactors

  double precision function gnedin2004Solver(radiusInitial)
    !% Root function used in finding the initial radius in the dark matter halo when solving for adiabatic contraction.
    implicit none
    double precision, intent(in   ) :: radiusInitial
    double precision                :: massDarkMatterInitial, radiusInitialMean

    ! Find the initial mean orbital radius.
    radiusInitialMean    =gnedin2004Self                   %radiusOrbitalMean(               radiusInitial    )
    ! Get the mass of dark matter inside the initial radius.
    massDarkMatterInitial=gnedin2004Self%darkMatterProfile_%enclosedMass     (gnedin2004Node,radiusInitialMean)
    ! Compute the root function.
    gnedin2004Solver=+massDarkMatterInitial                                             &
         &           *(                                                                 &
         &             +gnedin2004Self%initialMassFraction*               radiusInitial &
         &             -gnedin2004Self%darkMatterFraction *gnedin2004Self%radiusFinal   &
         &            )                                                                 &
         &           -gnedin2004Self%baryonicFinalTerm
    return
  end function gnedin2004Solver

  double precision function gnedin2004DerivativeSolver(radiusInitialDerivative)
    !% Root function used in finding the derivative of the initial radius in the dark matter halo when solving for adiabatic
    !% contraction.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: radiusInitialDerivative
    double precision                :: densityDarkMatterInitial, massDarkMatterInitial, &
         &                             radiusInitialMean

    ! Find the initial mean orbital radius.
    radiusInitialMean       =gnedin2004Self                   %radiusOrbitalMean(               gnedin2004Self%radiusInitial    )
    ! Get the mass of dark matter inside the initial radius.
    massDarkMatterInitial   =gnedin2004Self%darkMatterProfile_%enclosedMass     (gnedin2004Node,               radiusInitialMean)
    ! Get the mass of dark matter inside the initial radius.
    densityDarkMatterInitial=gnedin2004Self%darkMatterProfile_%density          (gnedin2004Node,               radiusInitialMean)
    ! Compute the root function.
    gnedin2004DerivativeSolver=+massDarkMatterInitial                                      &
         & *(                                                                              &
         &   +gnedin2004Self%initialMassFraction*               radiusInitialDerivative    &
         &   -gnedin2004Self%darkMatterFraction                                            &
         &  )                                                                              &
         & +(                                                                              &
         &   +gnedin2004Self%initialMassFraction*gnedin2004Self%radiusInitial              &
         &   -gnedin2004Self%darkMatterFraction *gnedin2004Self%radiusFinal                &
         &  )                                                                              &
         &  *4.0d0                                                                         &
         &  *Pi                                                                            &
         &  *radiusInitialMean**2                                                          &
         &  *densityDarkMatterInitial                                                      &
         &  *gnedin2004Self%radiusInitialMeanSelfDerivative                                &
         &  *               radiusInitialDerivative                                        &
         &  -gnedin2004Self%baryonicFinalTerm                                              &
         &  *(                                                                             &
         &    +1.0d0                                       /gnedin2004Self%radiusFinal     &
         &    +gnedin2004Self%radiusFinalMeanSelfDerivative/gnedin2004Self%radiusFinalMean &
         &   )                                                                             &
         &  -gnedin2004Self%baryonicFinalTermDerivative
    return
  end function gnedin2004DerivativeSolver

  double precision function gnedin2004RadiusOrbitalMean(self,radius)
    !% Returns the orbit averaged radius for dark matter corresponding the given {\normalfont \ttfamily radius} using the model of
    !% \cite{gnedin_response_2004}.
    implicit none
    class           (galacticStructureRadiiInitialGnedin2004), intent(inout) :: self
    double precision                                         , intent(in   ) :: radius

    gnedin2004RadiusOrbitalMean=+self%A                                                          &
         &                      *self%radiusVirial                                               &
         &                      *self%radiusExponentiator%exponentiate(radius/self%radiusVirial)
    return
  end function gnedin2004RadiusOrbitalMean

  double precision function gnedin2004RadiusOrbitalMeanDerivative(self,radius)
    !% Returns the derivative of the orbit averaged radius for dark matter corresponding the given {\normalfont \ttfamily radius} using the model of
    !% \cite{gnedin_response_2004}.
    implicit none
    class           (galacticStructureRadiiInitialGnedin2004), intent(inout) :: self
    double precision                                         , intent(in   ) :: radius

    gnedin2004RadiusOrbitalMeanDerivative=+self%A                &
         &                                *(                     &
         &                                  +     radius         &
         &                                  /self%radiusVirial   &
         &                                 )**(self%omega-1.0d0)
    return
  end function gnedin2004RadiusOrbitalMeanDerivative

  double precision function gnedin2004MassEnclosed(component)
    !% Unary function returning the enclosed mass in a component. Suitable for mapping over components. Ignores the dark matter
    !% profile.
    use Galacticus_Nodes, only : nodeComponent, nodeComponentDarkMatterProfile
    implicit none
    class(nodeComponent), intent(inout) :: component

    select type (component)
    class is (nodeComponentDarkMatterProfile)
       gnedin2004MassEnclosed=0.0d0
    class default
       gnedin2004MassEnclosed=component%enclosedMass(radiusLarge,gnedin2004ComponentType,gnedin2004MassType,gnedin2004WeightBy,gnedin2004WeightIndex,gnedin2004HaloLoaded)
    end select
    return
  end function gnedin2004MassEnclosed

  double precision function gnedin2004VelocityCircularSquared(component)
    !% Unary function returning the squared rotation curve in a component. Suitable for mapping over components.
    use Galacticus_Nodes, only : nodeComponent, nodeComponentDarkMatterProfile
    implicit none
    class(nodeComponent), intent(inout) :: component

    select type (component)
    class is (nodeComponentDarkMatterProfile)
       gnedin2004VelocityCircularSquared=0.0d0
    class default
       gnedin2004VelocityCircularSquared=component%rotationCurve(gnedin2004Self%radiusShared,gnedin2004ComponentType,gnedin2004MassType,gnedin2004HaloLoaded)**2
    end select
    return
  end function gnedin2004VelocityCircularSquared

  double precision function gnedin2004VelocityCircularSquaredGradient(component)
    !% Unary function returning the squared rotation curve gradient in a component. Suitable for mapping over components.
    use Galacticus_Nodes, only : nodeComponent, nodeComponentDarkMatterProfile
    implicit none
    class(nodeComponent), intent(inout) :: component

    select type (component)
    class is (nodeComponentDarkMatterProfile)
       gnedin2004VelocityCircularSquaredGradient=0.0d0
    class default
       gnedin2004VelocityCircularSquaredGradient=component%rotationCurveGradient(gnedin2004Self%radiusShared,gnedin2004ComponentType,gnedin2004MassType,gnedin2004HaloLoaded)
    end select
    return
  end function gnedin2004VelocityCircularSquaredGradient

  subroutine gnedin2004CalculationReset(self,node)
    !% Remove memory of stored computed values as we're about to begin computing derivatives anew.
    use Galacticus_Nodes, only : treeNode
    implicit none
    class(galacticStructureRadiiInitialGnedin2004), intent(inout) :: self
    type (treeNode                               ), intent(inout) :: node

    self%uniqueIDPrevious          =node%uniqueID()
    self%radiusPreviousIndex       = 0
    self%radiusPreviousIndexMaximum= 0
    self%radiusPrevious            =-1.0d0
    self%massesComputed            =.false.
    return
  end subroutine gnedin2004CalculationReset

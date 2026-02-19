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
  An implementation of adiabaticGnedin2004 dark matter halo profiles.
  !!}

  use :: Cosmology_Parameters      , only : cosmologyParameters              , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale              , darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMO             , darkMatterProfileDMOClass
  use :: Mass_Distributions        , only : enumerationNonAnalyticSolversType

  !![
  <darkMatterProfile name="darkMatterProfileAdiabaticGnedin2004">
   <description>
    A dark matter profile class which builds \refClass{massDistributionSphericalAdiabaticGnedin2004} objects to apply adiabatic
    contraction to other dark matter profiles.
   </description>
  </darkMatterProfile>
  !!]
  type, extends(darkMatterProfileClass) :: darkMatterProfileAdiabaticGnedin2004
     !!{
     A dark matter halo profile class implementing adiabaticGnedin2004 dark matter halos.
     !!}
     private
     class           (cosmologyParametersClass         ), pointer  :: cosmologyParameters_   => null()
     class           (darkMatterProfileDMOClass        ), pointer  :: darkMatterProfileDMO_  => null()
     class           (darkMatterHaloScaleClass         ), pointer  :: darkMatterHaloScale_   => null()
     type            (enumerationNonAnalyticSolversType)           :: nonAnalyticSolver
     double precision                                              :: A                               , omega               , &
          &                                                           radiusFractionalPivot           , toleranceRelative   , &
          &                                                           darkMatterFraction
   contains    
     final     ::        adiabaticGnedin2004Destructor
     procedure :: get => adiabaticGnedin2004Get
  end type darkMatterProfileAdiabaticGnedin2004

  interface darkMatterProfileAdiabaticGnedin2004
     !!{
     Constructors for the \refClass{darkMatterProfileAdiabaticGnedin2004} dark matter halo profile class.
     !!}
     module procedure adiabaticGnedin2004ConstructorParameters
     module procedure adiabaticGnedin2004ConstructorInternal
  end interface darkMatterProfileAdiabaticGnedin2004

contains

  function adiabaticGnedin2004ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily adiabaticGnedin2004} dark matter halo profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    use :: Input_Parameters  , only : inputParameters
    implicit none
    type            (darkMatterProfileAdiabaticGnedin2004)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (cosmologyParametersClass            ), pointer       :: cosmologyParameters_
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass           ), pointer       :: darkMatterProfileDMO_
    type            (varying_string                      )                :: nonAnalyticSolver
    double precision                                                      :: A                    , omega            , &
          &                                                                  radiusFractionalPivot, toleranceRelative
    
    !![
    <inputParameter>
      <name>A</name>
      <defaultSource>(\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultSource>
      <defaultValue>0.80d0</defaultValue>
      <description>The parameter $A$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>omega</name>
      <defaultSource>(\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultSource>
      <defaultValue>0.77d0</defaultValue>
      <description>The parameter $\omega$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusFractionalPivot</name>
      <defaultSource>\citep{gnedin_response_2004}</defaultSource>
      <defaultValue>1.0d0</defaultValue>
      <description>The pivot radius (in units of the virial radius), $r_0$, appearing in equation~(\ref{eq:adiabaticContractionGnedinPowerLaw}).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelative</name>
      <defaultValue>1.0d-2</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in solving for the initial radius in the adiabatically-contracted dark matter profile.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring adiabatic contraction by baryons is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=darkMatterProfileAdiabaticGnedin2004(A,omega,radiusFractionalPivot,toleranceRelative,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_" />
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function adiabaticGnedin2004ConstructorParameters

  function adiabaticGnedin2004ConstructorInternal(A,omega,radiusFractionalPivot,toleranceRelative,nonAnalyticSolver,cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileAdiabaticGnedin2004} dark matter profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversIsValid
    use :: Error             , only : Error_Report
    implicit none
    type            (darkMatterProfileAdiabaticGnedin2004)                        :: self
    double precision                                      , intent(in   )         :: A                    , omega            , &
         &                                                                           radiusFractionalPivot, toleranceRelative
    class           (cosmologyParametersClass            ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterProfileDMOClass           ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    type            (enumerationNonAnalyticSolversType   ), intent(in   )         :: nonAnalyticSolver
    !![
    <constructorAssign variables="A, omega, radiusFractionalPivot, toleranceRelative, nonAnalyticSolver, *cosmologyParameters_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]
    
    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    ! Evaluate the dark matter fraction.
    self%darkMatterFraction=+1.0d0                                   &
         &                  -self%cosmologyParameters_%OmegaBaryon() &
         &                  /self%cosmologyParameters_%OmegaMatter()
    return
  end function adiabaticGnedin2004ConstructorInternal

  subroutine adiabaticGnedin2004Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileAdiabaticGnedin2004} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileAdiabaticGnedin2004), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine adiabaticGnedin2004Destructor

  function adiabaticGnedin2004Get(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                       , massTypeDark                       , massTypeBaryonic             , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalAdiabaticGnedin2004, kinematicsDistributionCollisionless, massDistributionSpherical    , kinematicsDistributionClass, &
         &                                    sphericalAdiabaticGnedin2004Initializor     , kinematicsDistributionUndecorator  , nonAnalyticSolversFallThrough
    implicit none
    class           (massDistributionClass                  ), pointer                 :: massDistribution_
    class           (kinematicsDistributionClass            ), pointer                 :: kinematicsDistribution_      , kinematicsDistribution__
    class           (darkMatterProfileAdiabaticGnedin2004   ), intent(inout), target   :: self
    type            (treeNode                               ), intent(inout), target   :: node
    type            (enumerationWeightByType                ), intent(in   ), optional :: weightBy
    integer                                                  , intent(in   ), optional :: weightIndex
    class           (massDistributionClass                  ), pointer                 :: massDistributionDecorated    , massDistributionBaryonic
    double precision                                                                   :: massBaryonicSelfTotal        , massBaryonicTotal        , &
         &                                                                                darkMatterDistributedFraction, initialMassFraction
    procedure       (sphericalAdiabaticGnedin2004Initializor), pointer                 :: initializationFunction
    class           (*                                      ), pointer                 :: initializationSelf           , initializationArgument
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]
    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Set the baryonic component to zero - we will compute this later during initialization.
    massDistributionBaryonic      => null()
    massBaryonicTotal             =  0.0d0
    massBaryonicSelfTotal         =  0.0d0
    darkMatterDistributedFraction =  0.0d0
    initialMassFraction           =  0.0d0
    ! Create the mass distribution.
    allocate(massDistributionSphericalAdiabaticGnedin2004 :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalAdiabaticGnedin2004)
       massDistributionDecorated => self%darkMatterProfileDMO_%get(node,weightBy,weightIndex)
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          initializationFunction => adiabaticGnedin2004Initialize
          initializationSelf     => self
          initializationArgument => node
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalAdiabaticGnedin2004(                                                                                             &amp;
	      &amp;                                        A                            =self                     %A                                  , &amp;
	      &amp;                                        omega                        =self                     %omega                              , &amp;
	      &amp;                                        radiusVirial                 =self%darkMatterHaloScale_%radiusVirial                 (node), &amp;
	      &amp;                                        radiusFractionalPivot        =self                     %radiusFractionalPivot              , &amp;
	      &amp;                                        darkMatterFraction           =self                     %darkMatterFraction                 , &amp;
	      &amp;                                        darkMatterDistributedFraction=                          darkMatterDistributedFraction      , &amp;
	      &amp;                                        massFractionInitial          =                          initialMassFraction                , &amp;
	      &amp;                                        nonAnalyticSolver            =self                     %nonAnalyticSolver                  , &amp;
	      &amp;                                        toleranceRelative            =self                     %toleranceRelative                  , &amp;
	      &amp;                                        massDistribution_            =                          massDistributionDecorated          , &amp;
	      &amp;                                        massDistributionBaryonic     =                          massDistributionBaryonic           , &amp;
	      &amp;                                        initializationFunction       =                          initializationFunction             , &amp;
	      &amp;                                        initializationSelf           =                          initializationSelf                 , &amp;
	      &amp;                                        initializationArgument       =                          initializationArgument             , &amp;
              &amp;                                        componentType                =                          componentTypeDarkHalo              , &amp;
              &amp;                                        massType                     =                          massTypeDark                         &amp;
              &amp;                                       )
	    </constructor>
          </referenceConstruct>
          !!]
          kinematicsDistribution__ => massDistributionDecorated%kinematicsDistribution()
          if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
             allocate(kinematicsDistributionUndecorator :: kinematicsDistribution_)
             select type (kinematicsDistribution_)
             type is (kinematicsDistributionUndecorator  )
                !![
		<referenceConstruct object="kinematicsDistribution_">
		  <constructor>
		    kinematicsDistributionUndecorator(kinematicsDistribution__)
		  </constructor>
		</referenceConstruct>
             !!]
             end select
          else
             allocate(kinematicsDistributionCollisionless :: kinematicsDistribution_)
             select type (kinematicsDistribution_)
             type is (kinematicsDistributionCollisionless)
                !![
		<referenceConstruct object="kinematicsDistribution_">
		  <constructor>
		    kinematicsDistributionCollisionless(kinematicsDistribution__)
		  </constructor>
		</referenceConstruct>
                !!]
             end select
          end if
          call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
          !![
	  <objectDestructor name="kinematicsDistribution_"  />
	  <objectDestructor name="kinematicsDistribution__" />
          <objectDestructor name="massDistributionDecorated"/>
          !!]
       class default
          call Error_Report('expected a spherical mass distribution'//{introspection:location})
       end select
    end select
    return
  end function adiabaticGnedin2004Get

  subroutine adiabaticGnedin2004Initialize(self,node,massDistributionBaryonic,darkMatterDistributedFraction,massFractionInitial)
    !!{
    Initialize the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : massTypeBaryonic
    use :: Mass_Distributions        , only : massDistributionSphericalAdiabaticGnedin2004
    use :: Error                     , only : Error_Report
    implicit none
    class           (*                    ), intent(inout), target  :: self                         , node
    class           (massDistributionClass), intent(  out), pointer :: massDistributionBaryonic
    double precision                       , intent(  out)          :: darkMatterDistributedFraction, massFractionInitial
    type            (treeNode             )               , pointer :: nodeCurrent
    class           (nodeComponentBasic   )               , pointer :: basic
    double precision                                                :: massBaryonicSelf             , massBaryonicTotal  , &
         &                                                             massBaryonicSubhalos
    
    select type (self)
    type is (darkMatterProfileAdiabaticGnedin2004)
       select type (node)
       type is (treeNode)
          ! Compute the initial baryonic contribution from this halo, and any satellites.
          massDistributionBaryonic => node%massDistribution(massType=massTypeBaryonic)
          massBaryonicSelf         =  node%massBaryonic    (                         )
          ! Compute baryonic mass in subhalos.
          massBaryonicSubhalos     =  0.0d0
          nodeCurrent              => node
          do while (associated(nodeCurrent%firstSatellite))
             nodeCurrent => nodeCurrent%firstSatellite
          end do
          if (associated(nodeCurrent,node)) nodeCurrent => null()
          do while (associated(nodeCurrent))
             massBaryonicSubhalos=+            massBaryonicSubhalos   &
                  &               +nodeCurrent%massBaryonic        ()
             if (associated(nodeCurrent%sibling)) then
                nodeCurrent => nodeCurrent%sibling
                do while (associated(nodeCurrent%firstSatellite))
                   nodeCurrent => nodeCurrent%firstSatellite
                end do
             else
                nodeCurrent => nodeCurrent%parent
                if (associated(nodeCurrent,node)) nodeCurrent => null()
             end if
          end do
          ! Compute the total baryonic mass.
          massBaryonicTotal=+massBaryonicSubhalos &
               &            +massBaryonicSelf
          ! Limit masses to physical values.
          massBaryonicSelf =max(massBaryonicSelf ,0.0d0)
          massBaryonicTotal=max(massBaryonicTotal,0.0d0)
          ! Compute the fraction of matter assumed to be distributed like the dark matter.
          basic                         => node%basic()
          darkMatterDistributedFraction =  min(self%darkMatterFraction+(massBaryonicTotal-massBaryonicSelf)/basic%mass(),1.0d0)
          ! Compute the initial mass fraction.
          massFractionInitial           =  min(self%darkMatterFraction+ massBaryonicTotal                  /basic%mass(),1.0d0)
        class default
          call Error_Report('unexpected class'//{introspection:location})
       end select
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine adiabaticGnedin2004Initialize

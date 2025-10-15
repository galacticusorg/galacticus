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
  An implementation of dark matter halo profile concentrations using the \cite{ludlow_mass-concentration-redshift_2016}
  algorithm.
  !!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Root_Finder             , only : rootFinder

  !![
  <darkMatterProfileScaleRadius name="darkMatterProfileScaleRadiusLudlow2016">
   <description>
    Dark matter halo scale radii are computed using the algorithm of \cite{ludlow_mass-concentration-redshift_2016}. While
    \cite{ludlow_mass-concentration-redshift_2016} used $\Delta = 200 \rho_\mathrm{crit}$ to define halos, their model actually
    predicts the scale radius, $r_{-2}$, rather than the concentration. Therefore, here we report that the
    \cite{ludlow_mass-concentration-redshift_2016} concentrations are defined using the model's own virial density contrast
    definition --- this ensures that the predicted scale radii are applied directly to model halos.
   </description>
  </darkMatterProfileScaleRadius>
  !!]
  type, extends(darkMatterProfileScaleRadiusClass) :: darkMatterProfileScaleRadiusLudlow2016
     !!{
     A dark matter halo profile concentration class implementing the algorithm of
     \cite{ludlow_mass-concentration-redshift_2016}.
     !!}
     private
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_           => null()
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_          => null()
     class           (darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
     class           (virialDensityContrastClass       ), pointer :: virialDensityContrast_        => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_         => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_          => null()
     double precision                                             :: C                                      , f              , &
          &                                                          timeFormationSeekDelta                 , densityContrast
   contains
     !![
     <methods>
       <method description="Evaluate a function which goes to zero at the formation time of the tree."          method="formationTimeRoot"           />
       <method description="Initialize a root finder object for use in finding the formation time of the tree." method="formationTimeRootFunctionSet"/>
       <method description="Specifies if a branch history is required for application of the algorithm."        method="requireBranchHistory"        />
     </methods>
     !!]
     final             ::                                 ludlow2016Destructor
     procedure         :: radius                       => ludlow2016Radius
     procedure, nopass :: formationTimeRoot            => ludlow2016FormationTimeRoot
     procedure         :: formationTimeRootFunctionSet => ludlow2016FormationTimeRootFunctionSet
     procedure         :: requireBranchHistory         => ludlow2016RequireBranchHistory
  end type darkMatterProfileScaleRadiusLudlow2016

  interface darkMatterProfileScaleRadiusLudlow2016
     !!{
     Constructors for the \refClass{darkMatterProfileScaleRadiusLudlow2016} dark matter halo profile concentration class.
     !!}
     module procedure ludlow2016ConstructorParameters
     module procedure ludlow2016ConstructorInternal
  end interface darkMatterProfileScaleRadiusLudlow2016

  ! Array used to store state.
  type :: ludlow2016State
     class           (darkMatterProfileScaleRadiusLudlow2016), pointer     :: self                   => null()
     class           (cosmologyFunctionsClass               ), pointer     :: cosmologyFunctions_    => null()
     type            (treeNode                              ), pointer     :: node                   => null()
     type            (rootFinder                            ), allocatable :: finder
     double precision                                                      :: massHaloCharacteristic          , massLimit      , &
          &                                                                   hubbleParameterPresent          , densityContrast, &
          &                                                                   timePrevious
  end type ludlow2016State
  type   (ludlow2016State), allocatable, dimension(:) :: states
  integer                 , parameter                 :: statesIncrement=10
  integer                                             :: stateCount     = 0
  !$omp threadprivate(states, stateCount)

contains

  function ludlow2016ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily ludlow2016} dark matter halo profile concentration class.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileScaleRadiusLudlow2016)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass              ), pointer       :: cosmologyParameters_
    class           (darkMatterProfileScaleRadiusClass     ), pointer       :: darkMatterProfileScaleRadius_
    class           (virialDensityContrastClass            ), pointer       :: virialDensityContrast_
    class           (darkMatterProfileDMOClass             ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass              ), pointer       :: darkMatterHaloScale_
    double precision                                                        :: C                            , f, &
         &                                                                     timeFormationSeekDelta

    if (.not.parameters%isPresent('darkMatterProfileScaleRadius')) call Error_Report('a fallback scale radius method must be specified'//{introspection:location})
    !![
    <inputParameter>
      <name>C</name>
      <source>parameters</source>
      <defaultValue>400.0d0</defaultValue>
      <description>The parameter $C$ appearing in the halo concentration algorithm of \cite{ludlow_mass-concentration-redshift_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>f</name>
      <source>parameters</source>
      <defaultValue>0.02d0</defaultValue>
      <description>The parameter $f$ appearing in the halo concentration algorithm of \cite{ludlow_mass-concentration-redshift_2016}.</description>
    </inputParameter>
    <inputParameter>
      <name>timeFormationSeekDelta</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The parameter $\Delta \log t$ by which the logarithm of the trial formation time is incremented when stepping through the formation history of a node to find the formation time. If set to zero (or a negative value) the cumulative mass histories of nodes are assumed to be monotonic functions of time, and the formation time is instead found by a root finding algorithm,</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"           name="cosmologyFunctions_"           source="parameters"/>
    <objectBuilder class="cosmologyParameters"          name="cosmologyParameters_"          source="parameters"/>
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    <objectBuilder class="virialDensityContrast"        name="virialDensityContrast_"        source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"         name="darkMatterProfileDMO_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"          name="darkMatterHaloScale_"          source="parameters"/>
    !!]
    self=darkMatterProfileScaleRadiusLudlow2016(C,f,timeFormationSeekDelta,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileScaleRadius_,virialDensityContrast_,darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"          />
    <objectDestructor name="cosmologyParameters_"         />
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    <objectDestructor name="virialDensityContrast_"       />
    <objectDestructor name="darkMatterProfileDMO_"        />
    <objectDestructor name="darkMatterHaloScale_"         />
    !!]
    return
  end function ludlow2016ConstructorParameters

  function ludlow2016ConstructorInternal(C,f,timeFormationSeekDelta,cosmologyFunctions_,cosmologyParameters_,darkMatterProfileScaleRadius_,virialDensityContrast_,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileScaleRadiusLudlow2016} dark matter halo profile concentration class.
    !!}
    implicit none
    type            (darkMatterProfileScaleRadiusLudlow2016)                        :: self
    double precision                                        , intent(in   )         :: C                            , f, &
         &                                                                             timeFormationSeekDelta
    class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass              ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterProfileScaleRadiusClass     ), intent(in   ), target :: darkMatterProfileScaleRadius_
    class           (virialDensityContrastClass            ), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterProfileDMOClass             ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass              ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="C, f, timeFormationSeekDelta, *cosmologyFunctions_, *cosmologyParameters_, *darkMatterProfileScaleRadius_, *virialDensityContrast_, *darkMatterProfileDMO_, *darkMatterHaloScale_"/>
    !!]

    ! Find the density contrast as used to define masses by Ludlow et al. (2016).
    self%densityContrast=200.0d0/self%cosmologyParameters_%omegaMatter()
    return
  end function ludlow2016ConstructorInternal

  subroutine ludlow2016Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileScaleRadiusLudlow2016} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileScaleRadiusLudlow2016), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"          />
    <objectDestructor name="self%cosmologyParameters_"         />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    <objectDestructor name="self%virialDensityContrast_"       />
    <objectDestructor name="self%darkMatterProfileDMO_"        />
    <objectDestructor name="self%darkMatterHaloScale_"         />
    !!]
    return
  end subroutine ludlow2016Destructor

  double precision function ludlow2016Radius(self,node)
    !!{
    Return the scale radius of the dark matter halo profile of {\normalfont \ttfamily node} using the
    \cite{ludlow_mass-concentration-redshift_2016} algorithm.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Calculations_Resets                 , only : Calculations_Reset
    use :: Display                             , only : displayGreen                       , displayReset
    use :: Error                               , only : Error_Report                       , errorStatusSuccess
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , nodeComponentDarkMatterProfile, treeNode
    use :: Galactic_Structure_Options          , only : componentTypeDarkMatterOnly        , massTypeDark
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Merger_Tree_Walkers                 , only : mergerTreeWalkerIsolatedNodesBranch
    use :: Numerical_Comparison                , only : Values_Agree
    use :: Numerical_Constants_Math            , only : Pi
    use :: Root_Finder                         , only : rangeExpandMultiplicative          , rangeExpandSignExpectNegative , rangeExpandSignExpectPositive
    use :: String_Handling                     , only : stringXMLFormat
    implicit none
    class           (darkMatterProfileScaleRadiusLudlow2016), intent(inout), target       :: self
    type            (treeNode                              ), intent(inout), target       :: node
    type            (treeNode                              )               , pointer      :: nodeBranch
    class           (nodeComponentBasic                    )               , pointer      :: basic                        , basicBranch
    class           (nodeComponentDarkMatterProfile        )               , pointer      :: darkMatterProfile_           , darkMatterProfileChild_
    class           (massDistributionClass                 ), pointer                     :: massDistribution_
    integer                                                 , parameter                   :: iterationCountMaximum =100
    double precision                                        , parameter                   :: concentrationGuess    =10.0d0
    type            (ludlow2016State                       ), allocatable  , dimension(:) :: statesTmp
    type            (mergerTreeWalkerIsolatedNodesBranch   )                              :: treeWalker
    double precision                                                                      :: massHaloCharacteristic       , timeBranchEarliest     , &
         &                                                                                   densityMeanScaleRadius       , timeFormation          , &
         &                                                                                   radiusScale                  , radiusScalePrevious    , &
         &                                                                                   timeFormationTrial           , timeFormationEarliest  , &
         &                                                                                   timeFormationLatest          , timeFormationPrevious  , &
         &                                                                                   radiusScalePrevious2nd       , rangeExpandFactor
    integer                                                                               :: iterationCount               , i                      , &
         &                                                                                   status

    ! For halos with no progenitors, simply keep the fall-back result. Otherwise, perform our calculation
    if (self%requireBranchHistory() .and. .not.associated(node%firstChild)) then
       ludlow2016Radius=self%darkMatterProfileScaleRadius_%radius(node)
    else
       ! Increment the state counter. This is necessary to ensure that this function can be called recursively.
       if (allocated(states)) then
          if (stateCount == size(states)) then
             call move_alloc(states,statesTmp)
             allocate(states(size(statesTmp)+statesIncrement))
             states(1:size(statesTmp))=statesTmp
             do i=1,size(statesTmp)
                nullify(statesTmp(i)%self               )
                nullify(statesTmp(i)%node               )
                nullify(statesTmp(i)%cosmologyFunctions_)
             end do
             deallocate(statesTmp)
          end if
       else
          allocate(states(statesIncrement))
       end if
       stateCount=stateCount+1
       allocate(states(stateCount)%cosmologyFunctions_,mold=self%cosmologyFunctions_)
       !$omp critical(darkMatterProfilesStructureScaleLudlow2016DeepCopy)
       !![
       <deepCopyReset variables="self%cosmologyFunctions_"/>
       <deepCopy source="self%cosmologyFunctions_" destination="states(stateCount)%cosmologyFunctions_"/>
       <deepCopyFinalize variables="states(stateCount)%cosmologyFunctions_"/>
       !!]
       !$omp end critical(darkMatterProfilesStructureScaleLudlow2016DeepCopy)
       allocate(states(stateCount)%finder)
       call self%formationTimeRootFunctionSet(states(stateCount)%finder)
       ! Get the dark matter profile component of the node.
       darkMatterProfile_ => node%darkMatterProfile()
       ! Set an initial guess to the scale radius. We use the scale radius of the primary progenitor if available - under the
       ! assumption that the scale radius should change only slowly this should be a reasonable guess. If this is not available,
       ! use a fraction of the virial radius.
       if (associated(node%firstChild)) then
          darkMatterProfileChild_ =>   node                   %firstChild          %darkMatterProfile (    )
          radiusScalePrevious     =   +darkMatterProfileChild_                     %scale             (    )
       else
          radiusScalePrevious     =   +self                   %darkMatterHaloScale_%radiusVirial      (node) &
               &                      /                                            concentrationGuess
       end if
       radiusScalePrevious2nd     =   -huge(0.0d0)
       timeFormation              =   -huge(0.0d0)
       timeFormationPrevious      =   -huge(0.0d0)
       call darkMatterProfile_%scaleSet(radiusScalePrevious)
       call Calculations_Reset(node)
       ! Begin iteratively seeking a solution for the scale radius.
       iterationCount       =0
       timeFormationLatest  =0.0d0
       timeFormationEarliest=0.0d0
       do while (iterationCount < iterationCountMaximum)
          iterationCount=iterationCount+1
          ! Compute the characteristic halo mass, M₋₂.
          if (iterationCount == 1) then
             basic                                     =>  node                    %basic                     (                                        )
             states(stateCount)%self                   =>  self
             states(stateCount)%node                   =>  node
             states(stateCount)%hubbleParameterPresent =   self%cosmologyFunctions_%hubbleParameterEpochal    (expansionFactor=1.0d0                   )
             states(stateCount)%timePrevious           =  -1.0d0
             states(stateCount)%densityContrast        =  -huge(0.0d0)
          end if
          massDistribution_                            =>  self                     %darkMatterProfileDMO_ %get(node                                    )
          massHaloCharacteristic                       =  +massDistribution_        %massEnclosedBySphere      (darkMatterProfile_%scale()              )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
          states(stateCount)%massHaloCharacteristic    =   massHaloCharacteristic
          states(stateCount)%massLimit                 =  +self%f                                                                                                                                     &
               &                                                              *Dark_Matter_Profile_Mass_Definition(                                                                                   &
               &                                                                                                                          node                                                      , &
               &                                                                                                                          ludlow2016DensityContrast(states(stateCount),basic%time()), &
               &                                                                                                   cosmologyParameters_  =self%cosmologyParameters_                                 , &
               &                                                                                                   cosmologyFunctions_   =self%cosmologyFunctions_                                  , &
               &                                                                                                   virialDensityContrast_=self%virialDensityContrast_                               , &
               &                                                                                                   darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                  &
               &                                                                                                  )
          ! Determine the ranges of formation time.
          if (self%requireBranchHistory()) then
             ! Find the earliest time in the branch. Also estimate the earliest and latest times between which the formation time will lie.
             if (iterationCount == 1) then
                timeBranchEarliest   =huge(0.0d0)
                timeFormationEarliest=huge(0.0d0)
                timeFormationLatest  =basic%time()
                treeWalker           =mergerTreeWalkerIsolatedNodesBranch(node)
                do while (treeWalker%next(nodeBranch))
                   basicBranch =>nodeBranch%basic()
                   timeBranchEarliest                                                                       =min(timeBranchEarliest   ,basicBranch%time())
                   if (basicBranch%mass() > states(stateCount)%massLimit             ) timeFormationEarliest=min(timeFormationEarliest,basicBranch%time())
                   if (basicBranch%mass() > states(stateCount)%massHaloCharacteristic) timeFormationLatest  =min(timeFormationLatest  ,basicBranch%time())
                end do
                timeFormationLatest=max(timeFormationLatest,timeFormationEarliest)
             end if
          end if
          ! Test if the formation time is before the earliest time in the branch.
          if (self%requireBranchHistory() .and. self%formationTimeRoot(timeBranchEarliest) > 0.0d0) then
             ! The characteristic halo mass is never resolved in this branch - fall though to an alternative concentration calculation.
             radiusScale        =+self%darkMatterProfileScaleRadius_%radius(node)
             ! Force convergence by setting the previous scale radius to that which we just set - since we're choosing a (possibly
             ! random) concentration from a distribution we no longer need to iterate to find a solution.
             radiusScalePrevious=+radiusScale
          else if (self%formationTimeRoot(basic%time()) < 0.0d0) then
             ! The characteristic mass exceeds that at the present time. This is possible in early iterations if the progenitor
             ! halo was assigned a very large scale radius (usually from the fall-through method). Fall through to the alternative
             ! concentration calculation, but do not force convergence.
             radiusScale        =+self%darkMatterProfileScaleRadius_%radius(node)
          else
             ! Find the time at which the mass in progenitors equals this characteristic halo mass.
             if (iterationCount > 2) then
                if (radiusScalePrevious2nd > radiusScalePrevious) then
                   rangeExpandFactor=radiusScalePrevious2nd/radiusScalePrevious
                else
                   rangeExpandFactor=radiusScalePrevious   /radiusScalePrevious2nd
                end if
             else
                rangeExpandFactor=1.1d0
             end if
             if (self%requireBranchHistory()) then
                call states(stateCount)%finder%rangeExpand(                                                             &
                     &                                     rangeExpandUpward            =1.0d0*rangeExpandFactor      , &
                     &                                     rangeExpandDownward          =1.0d0/rangeExpandFactor      , &
                     &                                     rangeExpandType              =rangeExpandMultiplicative    , &
                     &                                     rangeUpwardLimit             =basic%time()                 , &
                     &                                     rangeDownwardLimit           =timeBranchEarliest           , &
                     &                                     rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
                     &                                     rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &
                     &                                    )
                if (iterationCount == 1) then
                   timeFormation=states(stateCount)%finder%find(rootRange=[timeFormationEarliest,timeFormationLatest],status=status)
                else
                   timeFormation=states(stateCount)%finder%find(rootGuess= timeFormationPrevious                     ,status=status)
                end if
             else
                call states(stateCount)%finder%rangeExpand(                                                             &
                     &                                     rangeExpandUpward            =1.0d0*rangeExpandFactor      , &
                     &                                     rangeExpandDownward          =1.0d0/rangeExpandFactor      , &
                     &                                     rangeExpandType              =rangeExpandMultiplicative    , &
                     &                                     rangeUpwardLimit             =basic%time()                 , &
                     &                                     rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
                     &                                     rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &
                     &                                    )
                if (timeFormationPrevious <= 0.0d0) then
                   timeFormation=states(stateCount)%finder%find(rootGuess= basic%time()                              ,status=status)
                else
                   timeFormation=states(stateCount)%finder%find(rootGuess= timeFormationPrevious                     ,status=status)
                end if
             end if
             if (status /= errorStatusSuccess)                                                                                                                                                                                &
                  & call Error_Report(                                                                                                                                                                                        &
                  &                   'solving for formation time failed'//char        (10)                                                                                                            //                     &
                  &                   displayGreen()//' HELP:'           //displayReset(  )                                                                                                            //                     &
                  &                   ' if you are using '//stringXMLFormat('<darkMatterProfileScaleRadius value="concentration"/>')                                                                   //char(10)//           &
                  &                   ' as the fall back method for setting scale radii, consider using mean concentrations in the fall-back method'                                                   //char(10)//           &
                  &                   ' (as scatter in the concentration-mass relation can lead to poor convergence here) by setting the highlighted'                                                  //char(10)//           &
                  &                   ' option in your parameter file as shown below:'                                                                                                                 //char(10)//char(10)// &
                  &                    stringXMLFormat('<darkMatterProfileScaleRadius value="concentration">**B<useMeanConcentration value="true"/>**C</darkMatterProfileScaleRadius>',indentInitial=3)//                     &
                  &                   {introspection:location}                                                                                                                                                                &
                  &                  )
             ! If requested, check for possible earlier formation times by simply stepping through trial times and finding the
             ! earliest at which the required mass threshold is reached. This is used for cases where the cumulative mass history
             ! is not monotonic.
             if (self%timeFormationSeekDelta > 0.0d0) then
                timeFormationTrial=timeBranchEarliest
                do while (timeFormationTrial < timeFormation)
                   if (self%formationTimeRoot(timeFormationTrial) > 0.0d0) then
                      timeFormation=timeFormationTrial
                      exit
                   end if
                   timeFormationTrial=+timeFormationTrial               &
                        &             *exp(self%timeFormationSeekDelta)
                end do
             end if
             ! Find the mean density with the scale radius, ⟨ρ₋₂⟩.
             densityMeanScaleRadius=+  self%C                                                                          &
                  &                 *  self%cosmologyParameters_%densityCritical       (                             ) &
                  &                 *(                                                                                 &
                  &                   +self%cosmologyFunctions_ %hubbleParameterEpochal(time           =timeFormation) &
                  &                   /self%cosmologyFunctions_ %hubbleParameterEpochal(expansionFactor=1.0d0        ) &
                  &                  )**2
             ! Find the corresponding scale radius.
             radiusScale=(                        &
                  &       +3.0d0                  &
                  &       *massHaloCharacteristic &
                  &       /densityMeanScaleRadius &
                  &       /4.0d0                  &
                  &       /Pi                     &
                  &      )**(1.0d0/3.0d0)
          end if
          ! Test for convergence.
          if (Values_Agree(radiusScale,radiusScalePrevious,relTol=1.0d-3)) exit
          ! Convergence was not attained - record current results and perform another iteration.
          call darkMatterProfile_%scaleSet(radiusScale)
          call Calculations_Reset(node)
          radiusScalePrevious2nd=radiusScalePrevious
          radiusScalePrevious   =radiusScale
          timeFormationPrevious =timeFormation
       end do
       ! Note that we do not check if convergence was actually reached. Due to the nature of the algorithm it is possible that the
       ! function for which we are seeking the root is discontinuous. As the scale radius changes, so does M₀ and therefore the
       ! mass limit f·M₀. As a result some progenitors will discontinuously pass/fail the mass limit check, making the root
       ! function discontinuous. Oscillating solutions can't be avoided in such situations, so we simply take the final iteration
       ! as our best estimate of the scale radius.
       ludlow2016Radius    =radiusScale
       ! Decrement the state index.
       !![
       <objectDestructor name="states(stateCount)%cosmologyFunctions_"/>
       !!]
       deallocate(states(stateCount)%finder)
       stateCount=stateCount-1
    end if
    return
  end function ludlow2016Radius

  subroutine ludlow2016FormationTimeRootFunctionSet(self,finder)
    !!{
    Initialize the finder object to compute the relevant formation history.
    !!}
    use :: Root_Finder, only : rootFinder
    implicit none
    class           (darkMatterProfileScaleRadiusLudlow2016), intent(inout) :: self
    type            (rootFinder                            ), intent(inout) :: finder
    double precision                                        , parameter     :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-4
    !$GLC attributes unused :: self
    
    finder=rootFinder(                                               &
         &            rootFunction     =ludlow2016FormationTimeRoot, &
         &            toleranceAbsolute=toleranceAbsolute          , &
         &            toleranceRelative=toleranceRelative            &
         &           )
    return
  end subroutine ludlow2016FormationTimeRootFunctionSet

  double precision function ludlow2016FormationTimeRoot(timeFormation)
    !!{
    Function used to find the formation time of a halo in the {\normalfont \ttfamily ludlow2016} concentration algorithm.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Merger_Tree_Walkers                 , only : mergerTreeWalkerIsolatedNodesBranch
    implicit none
    double precision                                     , intent(in   ) :: timeFormation
    type            (treeNode                           ), pointer       :: nodeBranch   , nodeChild        , &
         &                                                                  nodeSibling
    class           (nodeComponentBasic                 ), pointer       :: basicBranch  , basicChild       , &
         &                                                                  basicSibling
    type            (mergerTreeWalkerIsolatedNodesBranch)                :: treeWalker
    double precision                                                     :: massBranch   , massAccretionRate, &
         &                                                                  massSiblings , massProgenitor

    treeWalker=mergerTreeWalkerIsolatedNodesBranch(states(stateCount)%node,timeEarliest=timeFormation)
    massBranch=0.0d0
    do while (treeWalker%next(nodeBranch))
       basicBranch => nodeBranch%basic()
       if (associated(nodeBranch%firstChild).and.basicBranch%time() >= timeFormation) then
          nodeChild => nodeBranch%firstChild
          do while (associated(nodeChild))
             basicChild => nodeChild%basic()
             if (basicChild%time() < timeFormation) then
                ! Find the mass of the primary progenitor.
                massProgenitor=Dark_Matter_Profile_Mass_Definition(                                                                                        &
                     &                                                                    nodeChild                                                      , &
                     &                                                                    ludlow2016DensityContrast(states(stateCount),basicChild%time()), &
                     &                                             cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                   , &
                     &                                             cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                    , &
                     &                                             virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                 , &
                     &                                             darkMatterProfileDMO_ =states(stateCount)%self%darkMatterProfileDMO_                   &
                     &                                            )
                if (nodeChild%isPrimaryProgenitor()) then
                   ! Interpolate in mass for primary progenitors.
                   massSiblings =  massProgenitor
                   nodeSibling  => nodeChild%sibling
                   do while (associated(nodeSibling))
                      basicSibling =>   nodeSibling %basic  ()
                      massSiblings =   +massSiblings                                                                                                                  &
                           &           +Dark_Matter_Profile_Mass_Definition(                                                                                          &
                           &                                                                       nodeSibling                                                      , &
                           &                                                                       ludlow2016DensityContrast(states(stateCount),basicSibling%time()), &
                           &                                                cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                     , &
                           &                                                cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                      , &
                           &                                                virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                   , &
                           &                                                darkMatterProfileDMO_ =states(stateCount)%self%darkMatterProfileDMO_                     &
                           &                                               )
                      nodeSibling  =>   nodeSibling %sibling
                   end do
                   massAccretionRate=+(                                                                                                                               &
                        &              +Dark_Matter_Profile_Mass_Definition(                                                                                          &
                        &                                                                          nodeBranch                                                       , &
                        &                                                                          ludlow2016DensityContrast(states(stateCount),basicBranch %time()), &
                        &                                                   cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                     , &
                        &                                                   cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                      , &
                        &                                                   virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                   , &
                        &                                                   darkMatterProfileDMO_ =states(stateCount)%self%darkMatterProfileDMO_                      &
                        &                                                  )                                                                                          &
                        &              -massSiblings                                                                                                                  &
                        &             )                                                                                                                               &
                        &            /(                                                                                                                               &
                        &              +basicBranch%time()                                                                                                            &
                        &              -basicChild %time()                                                                                                            &
                        &             )
                   massProgenitor   =+massProgenitor               &
                        &            +massAccretionRate            &
                        &            *(                            &
                        &              +           timeFormation   &
                        &              -basicChild%time         () &
                        &             )
                   if (massProgenitor >= states(stateCount)%massLimit) &
                        & massBranch=+massBranch                       &
                        &            +massProgenitor
                else
                   ! No interpolation in mass for non-primary progenitors.
                   if (massProgenitor >= states(stateCount)%massLimit) &
                        & massBranch=+massBranch                       &
                        &            +massProgenitor
                end if
             end if
             nodeChild => nodeChild%sibling
          end do
       else if (.not.associated(nodeBranch%firstChild).and.basicBranch%time() == timeFormation) then
          massProgenitor=Dark_Matter_Profile_Mass_Definition(                                                                                        &
               &                                                                   nodeBranch                                                      , &
               &                                                                   ludlow2016DensityContrast(states(stateCount),basicBranch%time()), &
               &                                             cosmologyParameters_  =states(stateCount)%self%cosmologyParameters_                   , &
               &                                             cosmologyFunctions_   =states(stateCount)%self%cosmologyFunctions_                    , &
               &                                             virialDensityContrast_=states(stateCount)%self%virialDensityContrast_                 , &
               &                                             darkMatterProfileDMO_ =states(stateCount)%self%darkMatterProfileDMO_                    &
               &                                            )                                           
          if (massProgenitor >= states(stateCount)%massLimit) &
               & massBranch=+massBranch                       &
               &            +massProgenitor
       end if
    end do
    ludlow2016FormationTimeRoot=+massBranch                                &
         &                      -states(stateCount)%massHaloCharacteristic
    return
  end function ludlow2016FormationTimeRoot

  double precision function ludlow2016DensityContrast(state,time)
    !!{
    Compute the density contrast for this epoch.
    !!}
    implicit none
    type            (ludlow2016State), intent(inout) :: state
    double precision                 , intent(in   ) :: time

    if (time /= state%timePrevious) then
       state%timePrevious   =+time
       state%densityContrast=+  state%self                    %densityContrast                       &
            &                *(                                                                      &
            &                  +state%self%cosmologyFunctions_%hubbleParameterEpochal(time=time )    &
            &                  /state%                         hubbleParameterPresent                &
            &                 )                                                                  **2 &
            &                *  state%self%cosmologyFunctions_%expansionFactor       (     time )**3
    end if
    ludlow2016DensityContrast=state%densityContrast
    return
  end function ludlow2016DensityContrast

  logical function ludlow2016RequireBranchHistory(self) result(requireBranchHistory)
    !!{
    Specify if the branch history is required for the scale radius calculation.
    !!}
    implicit none
    class(darkMatterProfileScaleRadiusLudlow2016), intent(inout) :: self
    !$GLC attributes unused :: self

    requireBranchHistory=.true.
    return
  end function ludlow2016RequireBranchHistory

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
  Implements a node operator class that evaluates the properties of prompt cusps following the model of :cite:t:`delos_cusp-halo_2025`.
  !!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Power_Spectra           , only : powerSpectrumClass
  use :: Linear_Growth           , only : linearGrowthClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass

  !![
  <nodeOperator name="nodeOperatorDarkMatterProfilePromptCusps" docformat="rst">
   <description>
   A node operator class that evaluates the properties of prompt cusps following the model of :cite:t:`delos_cusp-halo_2025`, with
   a log-normal scatter of :math:`\mu \exp(-1/\sigma_0)` dex added to the cusp amplitude, where :math:`\mu=`\
   ``[coefficientScatter]``.

   Prompt cusp formation is, by definition, the formation event of a halo. A halo for which the cusp collapse epoch would lie in
   its own future can therefore not have formed: such halos are assigned no cusp and are labeled ``promptCuspUnformed``. **This
   class must be used together with a merger tree operator which prunes those labeled halos from the tree**, for example

   .. code-block:: xml

      &lt;mergerTreeOperator value="pruneByFilter"&gt;
        &lt;preservePrimaryProgenitor value="false"/&gt;
        &lt;galacticFilter value="labelled"&gt;
          &lt;label value="promptCuspUnformed"/&gt;
        &lt;/galacticFilter&gt;
      &lt;/mergerTreeOperator&gt;

   placed ahead of any other merger tree operator which depends on tree structure. Without it, halos which can not have formed are
   retained and evolved, and the branches above them are assigned cusp properties inconsistent with those halos.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfilePromptCusps
     !!{RST
     A node operator class that evaluates the properties of prompt cusps following the model of :cite:t:`delos_cusp-halo_2025`.
     !!}
     private
     class           (powerSpectrumClass        ), pointer                   :: powerSpectrum_              => null()
     class           (linearGrowthClass         ), pointer                   :: linearGrowth_               => null()
     class           (cosmologyParametersClass  ), pointer                   :: cosmologyParameters_        => null()
     class           (cosmologyFunctionsClass   ), pointer                   :: cosmologyFunctions_         => null()
     class           (darkMatterHaloScaleClass  ), pointer                   :: darkMatterHaloScale_        => null()
     class           (darkMatterProfileDMOClass ), pointer                   :: darkMatterProfileDMO_       => null()
     class           (virialDensityContrastClass), pointer                   :: virialDensityContrast_      => null()
     double precision                            , allocatable, dimension(:) :: sigma_
     logical                                                                 :: growthIsWavenumberDependent          , nonConvergenceIsFatal , &
          &                                                                     requireCollapseBeforeHalo
     double precision                                                        :: alpha                                , beta                  , &
          &                                                                     kappa                                , C                     , &
          &                                                                     p                                    , coefficientScatter
     integer                                                                 :: promptCuspMassID                     , promptCuspAmplitudeID , &
          &                                                                     promptCuspNFWYID                     , promptCuspNFWScaleID  , &
          &                                                                     promptCuspNFWGrowthRateID            , promptCuspNFWDensityID, &
          &                                                                     labelUnformedID
   contains
     !![
     <methods docformat="rst">
       <method method="sigma" description="Evaluate :math:`\sigma_j^2 = \int_0^\infty \frac{\mathrm{d}k}{k} \mathcal{P}(k,t) k^{2j}` where :math:`\mathcal{P}(k) = k^3 P(k) / 2 \pi^2` is the dimensionless form of the power spectrum."/>
     </methods>
     !!]
     final     ::                                        darkMatterProfilePromptCuspsDestructor
     procedure :: nodeTreeInitialize                  => darkMatterProfilePromptCuspsNodeTreeInitialize
     procedure :: nodeInitialize                      => darkMatterProfilePromptCuspsNodeInitialize
     procedure :: nodePromote                         => darkMatterProfilePromptCuspsNodePromote
     procedure :: differentialEvolutionSolveAnalytics => darkMatterProfilePromptCuspsSolveAnalytics
     procedure :: sigma                               => darkMatterProfilePromptCuspsNodeSigma
  end type nodeOperatorDarkMatterProfilePromptCusps
  
  interface nodeOperatorDarkMatterProfilePromptCusps
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorDarkMatterProfilePromptCusps` node operator class.
     !!}
     module procedure darkMatterProfilePromptCuspsConstructorParameters
     module procedure darkMatterProfilePromptCuspsConstructorInternal
  end interface nodeOperatorDarkMatterProfilePromptCusps

  ! Submodule-scope variables used in root-finding.
  class           (nodeOperatorDarkMatterProfilePromptCusps), pointer   :: self_
  double precision                                                      :: sigma0Collapse          , time_                    , &
       &                                                                   expansionFactor_        , concentrationFactorTarget, &
       &                                                                   massHalo_               , amplitudeCusp_           , &
       &                                                                   concentration_          , radiusScale_
  integer                                                               :: j_ 
  !$omp threadprivate(self_,sigma0Collapse,time_,expansionFactor_,concentrationFactorTarget,massHalo_,concentration_,radiusScale_,amplitudeCusp_,j_)

  ! Maximum allowed value of the y-parameter in the cusp-NFW profile. Values of 1 or greater are not valid. We limit here to a
  ! value close to 1.
  double precision                                          , parameter :: yMaximum        =0.99999d0

  ! Minimum value of the y-parameter for which the cusp-NFW expressions are evaluated; below this the NFW limit is used instead,
  ! and y is stored as zero so that the cusp-NFW mass distribution class (which contains the same expressions, and so the same
  ! problems) likewise takes its NFW branch.
  !
  ! The asinh and atanh terms of the enclosed mass each diverge logarithmically as y→0. Their divergences cancel, so the enclosed
  ! mass remains finite - but the cancellation is increasingly ill-conditioned as y decreases, because the argument of the atanh
  ! approaches 1 and so loses precision. Two limits follow:
  !
  !  * Once y ≲ 10⁻⁹ that argument rounds to exactly 1 in double precision, giving atanh(1) = ∞ and hence a floating point
  !    exception.
  !  * Well before that, the evaluation is simply inaccurate: the relative error grows as 1/y², reaching ~10⁻⁴ by y = 10⁻⁶.
  !
  ! Meanwhile the cusp itself contributes to the enclosed mass only at order y², so the error incurred by taking the NFW limit
  ! *falls* as y². The two curves cross at y ≈ 10⁻⁴, where each is of order 10⁻⁸ - far below the tolerance to which the density
  ! normalization is solved, and so the point at which the switch is least visible. Note that this is also far below the y of any
  ! halo with a real cusp (values of order 10⁻² and above), so the threshold is never approached in practice.
  double precision                                          , parameter :: yMinimum        =1.00000d-4

contains
  
  function darkMatterProfilePromptCuspsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorDarkMatterProfilePromptCusps` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkMatterProfilePromptCusps)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (linearGrowthClass                       ), pointer       :: linearGrowth_
    class           (powerSpectrumClass                      ), pointer       :: powerSpectrum_
    class           (cosmologyParametersClass                ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass               ), pointer       :: darkMatterProfileDMO_
    class           (virialDensityContrastClass              ), pointer       :: virialDensityContrast_
    logical                                                                   :: nonConvergenceIsFatal , requireCollapseBeforeHalo
    double precision                                                          :: alpha                 , beta                     , &
         &                                                                       kappa                 , C                        , &
         &                                                                       p                     , coefficientScatter

    !![
    <inputParameter docformat="rst">
      <name>nonConvergenceIsFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, failure to converge on a solution for the scale radius, :math:`r_\mathrm{s}`, will result in a fatal error. Otherwise, only warnings are issued.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>requireCollapseBeforeHalo</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, a prompt cusp is assigned only to halos for which the cusp collapse epoch precedes the epoch of the halo itself, as
      in the reference implementation of :cite:t:`delos_cusp-halo_2025`. Since prompt cusp formation is, by definition, the
      formation event of a halo, halos failing this condition can not have formed: they are assigned no cusp and are labeled
      ``promptCuspUnformed`` so that they may be pruned from the merger tree (see the class description above). Halos for which no
      collapse epoch exists at all are treated identically, again matching the reference implementation.

      If false, the original behavior of this class is used instead: any collapse epoch found is accepted, even where it postdates
      the halo, and halos with no solution are assigned a collapse epoch at the limiting time of the search. No halo is labeled,
      and so none is pruned. This is provided only to allow existing models to reproduce their previous results; it is not
      physically preferable.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>alpha</name>
      <source>parameters</source>
      <defaultValue>24.0d0</defaultValue>
      <defaultSource>
      :cite:p:`delos_cusp-halo_2025`
      </defaultSource>
      <description>
      The coefficient, :math:`\alpha` of the cusp amplitude, :math:`A`, in the peak-cusp connection of the :cite:t:`delos_cusp-halo_2025` prompt cusp model.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>7.3d0</defaultValue>
      <defaultSource>
      :cite:p:`delos_cusp-halo_2025`
      </defaultSource>
      <description>
      The coefficient, :math:`\beta`, of the cusp mass, :math:`m`, in the peak-cusp connection of the :cite:t:`delos_cusp-halo_2025` prompt cusp model.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>C</name>
      <source>parameters</source>
      <defaultValue>0.8d0</defaultValue>
      <defaultSource>
      :cite:p:`delos_cusp-halo_2025`
      </defaultSource>
      <description>
      The coefficient, :math:`C`, of the cusp :math:`A`--:math:`m` relation in the :cite:t:`delos_cusp-halo_2025` prompt cusp model.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>p</name>
      <source>parameters</source>
      <defaultValue>1.9d0</defaultValue>
      <defaultSource>
      :cite:p:`delos_cusp-halo_2025`
      </defaultSource>
      <description>
      The exponent, :math:`p`, of the cusp :math:`A`--:math:`m` relation in the :cite:t:`delos_cusp-halo_2025` prompt cusp model.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>kappa</name>
      <source>parameters</source>
      <defaultValue>4.5d0</defaultValue>
      <defaultSource>
      :cite:p:`delos_cusp-halo_2025`
      </defaultSource>
      <description>
      The parameter, :math:`\kappa`, of the mass growth factor in the :cite:t:`delos_cusp-halo_2025` prompt cusp model.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>coefficientScatter</name>
      <source>parameters</source>
      <defaultValue>0.195d0</defaultValue>
      <defaultSource>
      Delos (private communication)
      </defaultSource>
      <description>
      The parameter, :math:`\mu`, in the expression for the scatter in cusp amplitude.
      </description>
    </inputParameter>
    <objectBuilder class="linearGrowth"          name="linearGrowth_"          source="parameters"/>
    <objectBuilder class="powerSpectrum"         name="powerSpectrum_"         source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfilePromptCusps(nonConvergenceIsFatal,requireCollapseBeforeHalo,alpha,beta,kappa,C,p,coefficientScatter,linearGrowth_,powerSpectrum_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="linearGrowth_"         />
    <objectDestructor name="powerSpectrum_"        />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="darkMatterProfileDMO_" />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function darkMatterProfilePromptCuspsConstructorParameters

  function darkMatterProfilePromptCuspsConstructorInternal(nonConvergenceIsFatal,requireCollapseBeforeHalo,alpha,beta,kappa,C,p,coefficientScatter,linearGrowth_,powerSpectrum_,cosmologyParameters_,cosmologyFunctions_,virialDensityContrast_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorDarkMatterProfilePromptCusps` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Nodes_Labels    , only : nodeLabelRegister
    implicit none
    type            (nodeOperatorDarkMatterProfilePromptCusps)                        :: self
    class           (linearGrowthClass                       ), intent(in   ), target :: linearGrowth_
    class           (powerSpectrumClass                      ), intent(in   ), target :: powerSpectrum_
    class           (cosmologyParametersClass                ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass               ), intent(in   ), target :: darkMatterProfileDMO_
    class           (virialDensityContrastClass              ), intent(in   ), target :: virialDensityContrast_
    logical                                                   , intent(in   )         :: nonConvergenceIsFatal , requireCollapseBeforeHalo
    double precision                                          , intent(in   )         :: alpha                 , beta                     , &
         &                                                                               kappa                 , C                        , &
         &                                                                               p                     , coefficientScatter
    !![
    <constructorAssign variables="nonConvergenceIsFatal, requireCollapseBeforeHalo, alpha, beta, kappa, C, p, coefficientScatter, *linearGrowth_, *powerSpectrum_, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *virialDensityContrast_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="promptCuspAmplitude"       id="self%promptCuspAmplitudeID"     isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspMass"            id="self%promptCuspMassID"          isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWY"            id="self%promptCuspNFWYID"          isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWScale"        id="self%promptCuspNFWScaleID"      isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWDensity"      id="self%promptCuspNFWDensityID"    isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWGrowthRateID" id="self%promptCuspNFWGrowthRateID" isEvolvable="no" isCreator="yes"/>
    !!]
    self%growthIsWavenumberDependent=self%linearGrowth_%isWavenumberDependent()
    ! Register the label used to identify halos in which no prompt cusp can form and which, therefore, can not themselves have
    ! formed. Such nodes are expected to be pruned from the tree prior to evolution - see
    ! `darkMatterProfilePromptCuspsNodeTreeInitialize()` for details.
    self%labelUnformedID=nodeLabelRegister(                                                                                                         &
         &                                 'promptCuspUnformed'                                                                                   , &
         &                                 'Halos for which the prompt cusp collapse epoch is later than the epoch of the halo itself, and which'// &
         &                                 ' therefore can not have formed.'                                                                        &
         &                                )
    return
  end function darkMatterProfilePromptCuspsConstructorInternal

  subroutine darkMatterProfilePromptCuspsDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodeOperatorDarkMatterProfilePromptCusps` node operator class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfilePromptCusps), intent(inout) :: self

    !![
    <objectDestructor name="self%linearGrowth_"         />
    <objectDestructor name="self%powerSpectrum_"        />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%darkMatterHaloScale_"  />
    <objectDestructor name="self%darkMatterProfileDMO_" />
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine darkMatterProfilePromptCuspsDestructor

  subroutine darkMatterProfilePromptCuspsNodeTreeInitialize(self,node)
    !!{RST
    Initialize dark matter profile prompt cusp properties.
    !!}
    use :: Calculations_Resets                 , only : Calculations_Reset
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Error                               , only : Error_Report                       , Warn                          , errorStatusSuccess
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , nodeComponentDarkMatterProfile
    use :: Lambert_Ws                          , only : Lambert_W0
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Comparison                , only : Values_Agree
    use :: Nodes_Labels                        , only : nodeLabelSet
    use :: Root_Finder                         , only : rootFinder                         , rangeExpandMultiplicative     , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    use :: ISO_Varying_String                  , only : var_str
    use :: String_Handling                     , only : operator(//)
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    type            (treeNode                                )               , pointer :: nodeChild                        , nodeTip
    class           (nodeComponentBasic                      )               , pointer :: basic
    class           (nodeComponentDarkMatterProfile          )               , pointer :: darkMatterProfile
    integer                                                   , parameter              :: iterationCountMaximum    =10000
    double precision                                          , parameter              :: toleranceRelative        =1.0d-6
    type            (rootFinder                              )               , save    :: finderCollapse                   , finderRadius
    logical                                                                  , save    :: finderCollapseInitialized=.false., finderRadiusInitialized=.false.
    !$omp threadprivate(finderCollapse,finderRadius,finderCollapseInitialized,finderRadiusInitialized)
    ! Factor by which the collapse time is allowed to exceed the present age of the Universe before the search is abandoned. σ₀(t)
    ! grows only until the Universe becomes Λ-dominated and is then constant to far better than double precision by this epoch, so
    ! no root can exist beyond it.
    double precision                                          , parameter              :: factorTimeCollapseMaximum=1.0d1
    double precision                                                                   :: timeCollapseMaximum
    double precision                                                         , save    :: errorFractionalMaximum   =0.0d+0
    logical                                                                            :: computeCusp                      , cuspCanForm
    integer                                                                            :: iterationCount                   , statusCollapse
    double precision                                                                   :: sigma0                           , sigma2                         , &
         &                                                                                densityMean                      , densityMeanCollapse            , &
         &                                                                                sigma2Collapse                   , timeCollapse                   , &
         &                                                                                amplitude                        , mass                           , &
         &                                                                                y                                , radiusScale                    , &
         &                                                                                concentration                    , densityScale                   , &
         &                                                                                radiusScalePrevious              , mass200Critical                , &
         &                                                                                gamma                            , zeta                           , &
         &                                                                                radiusMinus2                     , densityContrast                , &
         &                                                                                scatterRandom                    , errorFractional                , &
         &                                                                                densityScalePrevious             , massPrevious
    
    ! Compute cusp properties for branch tip nodes. Note that "branch tip" here means a node with no progenitors *which remain
    ! after pruning* - i.e. progenitors which have themselves been identified as unable to form are ignored. Since the tree is
    ! walked depth-first, all progenitors of this node have already been processed, so their status is known.
    nodeChild   =>      progenitorPrimary(self%labelUnformedID,node     )
    computeCusp = .not.associated        (                     nodeChild)
    if (.not.computeCusp) then
       ! For non-branch tip nodes, use the cusp properties from the branch tip node.
       nodeTip   => nodeChild
       nodeChild => progenitorPrimary(self%labelUnformedID,nodeTip)
       do while (associated(nodeChild))
          nodeTip   => nodeChild
          nodeChild => progenitorPrimary(self%labelUnformedID,nodeTip)
       end do
       darkMatterProfile =>  nodeTip          %darkMatterProfile        (                          )
       amplitude         =   darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspAmplitudeID)
       mass              =   darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspMassID     )
       y                 =   darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspNFWYID     )
    else
       ! For branch tip nodes, assume no cusp on the initial iteration.
       y                 =   0.0d0
    end if
    ! Initialize the root finder.
    if (.not.finderCollapseInitialized) then
       finderCollapse=rootFinder(                                                             &
            &                    rootFunction                 =timeCollapseRoot             , &
            &                    toleranceRelative            =1.0d-6                       , &
            &                    rangeExpandUpward            =2.0d+0                       , &
            &                    rangeExpandDownward          =0.5d+0                       , &
            &                    rangeExpandType              =rangeExpandMultiplicative    , &
            &                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                    rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
            &                   )
      finderCollapseInitialized=.true.
    end if
    ! Evaluate the required integrals over the power spectrum at the time of this node.
    self_             =>  self
    basic             =>  node                                  %basic               (                   )
    ! Bound the upward expansion of the collapse time search. For sufficiently low mass halos no solution exists at all: σ₀Collapse
    ! then exceeds σ₀(t) at every epoch, because σ₀(t) grows only until the Universe becomes Λ-dominated and is constant
    ! thereafter. With no limit the range expansion doubles the trial time without bound, the expansion factor table is rebuilt
    ! out to ever later times, and the cosmological functions eventually raise a floating point exception when a³ overflows -
    ! which is diagnosed far from its cause. The limit is placed well beyond any epoch at which σ₀ is still growing, so it does
    ! not alter any collapse time that is found; it only converts the runaway into a clean out-of-range status here.
    timeCollapseMaximum=+factorTimeCollapseMaximum                                  &
         &              *self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    call finderCollapse%rangeExpand(                                                             &
         &                          rangeExpandUpward            =2.0d+0                       , &
         &                          rangeExpandDownward          =0.5d+0                       , &
         &                          rangeExpandType              =rangeExpandMultiplicative    , &
         &                          rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                          rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
         &                          rangeUpwardLimit             =timeCollapseMaximum            &
         &                         )
    darkMatterProfile =>  node                                  %darkMatterProfile   (                   )
    sigma0            =   self                                  %sigma               (0,     basic%time())
    sigma2            =   self                                  %sigma               (2,     basic%time())
    densityMean       =   self             %cosmologyFunctions_ %matterDensityEpochal(  time=basic%time())
    ! We assume (following Delos 2025) that the "scale radius" that has been set is actually r₋₂. We must iteratively solve for
    ! the actual scale radius. We use rₛ=r₋₂ as our initial guess.
    radiusMinus2        =darkMatterProfile%scale()
    radiusScale         =radiusMinus2
    radiusScalePrevious =huge(0.0d0)
    densityScale        =     0.0d0
    densityScalePrevious=huge(0.0d0)
    massPrevious        =huge(0.0d0)
    iterationCount      =0
    cuspCanForm         =.true.
    ! Draw the deviate used to scatter the cusp amplitude. This is drawn here, prior to any test of whether a prompt cusp can form
    ! in this halo, so that the sequence of random numbers consumed from the tree's generator does not depend on the outcome of
    ! that test.
    if (computeCusp .and. self%coefficientScatter > 0.0d0) then
       scatterRandom=node%hostTree%randomNumberGenerator_%standardNormalSample()
    else
       scatterRandom=0.0d0
    end if
    ! Begin iteration.
    do while (                                                                                 &
         &     (                                                                               &
         &       .not.Values_Agree(radiusScale ,radiusScalePrevious ,relTol=toleranceRelative) &
         &      .or.                                                                           &
         &       .not.Values_Agree(densityScale,densityScalePrevious,relTol=toleranceRelative) &
         &      .or.                                                                           &
         &       .not.Values_Agree(basic%mass(),massPrevious        ,relTol=toleranceRelative) &
         &     )                                                                               &
         &    .and.                                                                            &
         &     iterationCount < iterationCountMaximum                                          &
         &   )
       iterationCount      =iterationCount+1
       radiusScalePrevious =radiusScale
       densityScalePrevious=densityScale
       if (computeCusp) then
          ! Compute the mass following the 200 times critical density definition as used by Delos (2025; note that Delos actually
          ! states 200 times mean density, but these are Einstein-de Sitter cosmologies where mean and critical are equivalent -
          ! when comparing with other simulations Sten Delos uses 200 times critical [private communication]).
          call Calculations_Reset(node)
          densityContrast=+200.0d0                                                                    & ! 200 ρ_crit/ρ_mean since the density contrast is
               &          /self%cosmologyFunctions_%OmegaMatterEpochal(time=basic%timeLastIsolated())   ! relative to mean density.
          mass200Critical=Dark_Matter_Profile_Mass_Definition(                                                    &
               &                                              node                  =node                       , &
               &                                              densityContrast       =densityContrast            , &
               &                                              cosmologyParameters_  =self%cosmologyParameters_  , &
               &                                              cosmologyFunctions_   =self%cosmologyFunctions_   , &
               &                                              virialDensityContrast_=self%virialDensityContrast_, &
               &                                              darkMatterProfileDMO_ =self%darkMatterProfileDMO_ , &
               &                                              useLastIsolatedTime   =.true.                       &
               &                                             )
          ! Find the collapse time of this node.
          sigma0Collapse     =+            (2.0d0*self%p-1.0d0)*self%kappa/3.0d0              &
               &              /Lambert_W0(                                                    &
               &                          +(2.0d0*self%p-1.0d0)*self%kappa/3.0d0              &
               &                          *(self%beta*self%C**2/self%alpha**2)**(1.0d0/3.0d0) &
               &                          *(                                                  &
               &                            +exp(self%kappa/sigma0)                           &
               &                            *mass200Critical                                  &
               &                            /densityMean                                      &
               &                            /(sigma0/sigma2)**1.5d0                           &
               &                           )**((2.0d0*self%p-1.0d0)/3.0d0)                    &
               &                         )
          timeCollapse       =finderCollapse            %find(                  rootGuess=basic%time        (),status=statusCollapse)
          ! Test whether a prompt cusp can form in this halo. Following Delos (2025), a cusp exists only if the collapse epoch
          ! precedes that of the halo itself - `CuspHalo.collapse_a()` in the reference implementation
          ! (https://github.com/galacticusorg/cusp-halo-relation) returns NaN otherwise. Since prompt cusp formation is, by
          ! definition, the formation event of a halo, a halo in which no cusp can form is one which can not itself exist. Such
          ! nodes are assigned no cusp and are labeled, so that they may subsequently be pruned from the tree (for example by a
          ! `mergerTreeOperatorPruneByFilter` operator using a `galacticFilterLabelled` filter). The case in which no solution
          ! exists at all (σ₀Collapse exceeds σ₀(t) at every epoch, which occurs for sufficiently low mass halos, since σ₀(t)
          ! grows only until the Universe becomes Λ-dominated and is constant thereafter) is treated identically - just as in the
          ! reference implementation, where such halos fall off the end of the interpolating table and are likewise assigned no
          ! cusp.
          if (statusCollapse /= errorStatusSuccess .or. timeCollapse > basic%time()) then
             if (self%requireCollapseBeforeHalo) then
                cuspCanForm=.false.
                exit
             else if (statusCollapse /= errorStatusSuccess) then
                ! Original behavior: no collapse time exists for this halo, so take the collapse to occur at the limiting time,
                ! which is the continuous limit of the solution as σ₀Collapse approaches the asymptotic value of σ₀. Where a
                ! solution does exist it is accepted unchanged, even though it postdates the halo.
                timeCollapse=timeCollapseMaximum
             end if
          end if
          sigma2Collapse     =self                      %sigma               (2,time=timeCollapse)
          densityMeanCollapse=self  %cosmologyFunctions_%matterDensityEpochal(  time=timeCollapse)
          ! Evaluate the prompt cusp parameters.
          amplitude=self%alpha*(self%alpha/self%C/self%beta**self%p)**(1.0d0/(2.0d0*self%p-1.0d0))*densityMeanCollapse*sigma0Collapse**((9.0d0-6.0d0*self%p)/(4.0d0-8.0d0*self%p))/sigma2Collapse**0.75d0
          mass     =self%beta *(self%alpha/self%C/self%beta**self%p)**(2.0d0/(2.0d0*self%p-1.0d0))*densityMeanCollapse*sigma0Collapse**((9.0d0-6.0d0*self%p)/(2.0d0-4.0d0*self%p))/sigma2Collapse**1.50d0
          ! Add scatter to the cusp amplitude.
          if (self%coefficientScatter > 0.0d0) then
             amplitude=+amplitude                         &
                  &    *10.0d0**(                         &
                  &              +self%coefficientScatter &
                  &              *exp(                    &
                  &                   -1.0d0              &
                  &                   /sigma0             &
                  &                  )                    &
                  &              *scatterRandom           &
                  &             )
          end if
       end if
       ! Adjust the initial concentration (and r₋₂ radius) if necessary. The equation in footnote 9 of Delos (2025;
       ! https://ui.adsabs.harvard.edu/abs/2025arXiv250603240D) gives a condition that must be satisfied by the concentration
       ! c=rᵥ/r₋₂. If that is not satisfied, adjust r₋₂ to satisfy that relation.
       concentration            =+self %darkMatterHaloScale_%radiusVirial(node)    &
            &                    /                           radiusMinus2
       concentrationFactorTarget=+self %darkMatterHaloScale_%densityMean (node)    &
            &                    *basic                     %mass        (    )    &
            &                    /                           amplitude         **2
       if (concentrationTargetRoot(concentration) < 0.0d0) then
          if (.not.finderRadiusInitialized) then
             finderRadius=rootFinder(                                                             &
                  &                  rootFunction                 =concentrationTargetRoot      , &
                  &                  toleranceRelative            =1.0d-6                       , &
                  &                  rangeExpandUpward            =2.0d+0                       , &
                  &                  rangeExpandDownward          =0.5d+0                       , &
                  &                  rangeExpandType              =rangeExpandMultiplicative    , &
                  &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
                  &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
                  &                 )
             finderRadiusInitialized=.true.
          end if
          ! Find a new concentration and use it to update r₋₂ (including some small tolerance to avoid being exactly at the
          ! limit).
          concentration= finderRadius                     %find         (concentration)
          radiusMinus2 =+self        %darkMatterHaloScale_%radiusVirial (node         ) &
               &        /                                  concentration                &
               &        *                                  yMaximum     **(2.0d0/3.0d0)
          y            =+yMaximum
          radiusScale  =+2.0d0        &
               &        *radiusMinus2
          densityScale =+amplitude          &
               &        /radiusScale**1.5d0
          ! Force convergence in this case.
          radiusScalePrevious =      radiusScale
          densityScalePrevious=      densityScale
          massPrevious        =basic%mass        ()
       else
          ! Compute the normalization of the cusp-NFW profile. The y=0 case (which is assumed on the first iteration) is handled
          ! by `factorMassCuspNFW()`, which takes the NFW limit for negligible y.
          concentration=+self%darkMatterHaloScale_%radiusVirial(node) &
               &        /                          radiusScale
          densityScale =+basic            %mass(               ) &
               &        /radiusScale**3                          &
               &        /4.0d0                                   &
               &        /Pi                                      &
               &        /factorMassCuspNFW     (concentration,y)
          ! Compute a new scale radius to obtain the target r₋₂, using equations (23) from Delos (2025).
          gamma      =+amplitude           &
               &      /densityScale        &
               &      /radiusMinus2**1.5d0
          zeta       =+(                         &
               &        +        1.00d0          &
               &        +       18.00d0*gamma**2 &
               &        +        0.75d0*gamma    &
               &        *sqrt(                   &
               &              + 72.00d0          &
               &              +564.00d0*gamma**2 &
               &              +  6.00d0*gamma**4 &
               &             )                   &
               &       )**(1.0d0/3.0d0)
          radiusScale=+(                                     &
               &        +(+1.0d0/3.0d0-gamma**2/2.0d0)/zeta  &
               &        +(+1.0d0      +zeta          )/3.0d0 &
               &       )                                     &
               &      *radiusMinus2
          ! Enforce a maximum scale radius (minimum concentration) as required to ensure that y<1 in the cusp-NFW profile.
          radiusScale=max(radiusScale,(amplitude/densityScale/yMaximum)**(2.0d0/3.0d0))
          ! Evaluate the prompt cusp y-parameter. Values below `yMinimum` are indistinguishable from a cuspless (i.e. NFW) profile
          ! and are set to zero, so that no negligible - but non-zero - y is stored for use by the cusp-NFW mass distribution.
          y          =+amplitude           &
               &      /densityScale        &
               &      /radiusScale **1.5d0
          if (y < yMinimum) y=0.0d0
          ! Evaluate the mass to use in convergence checks.
          massPrevious=+4.0d0                              &
               &       *Pi                                 &
               &       *radiusScale**3                     &
               &       *densityScale                       &
               &       *factorMassCuspNFW(concentration,y)
       end if
       ! Store prompt cusp parameters.
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspAmplitudeID ,amplitude   )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspMassID      ,mass        )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID      ,y           )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWScaleID  ,radiusScale )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWDensityID,densityScale)
    end do
    ! Handle halos in which no prompt cusp can form, and which therefore can not themselves have formed. Such nodes are given a
    ! cusp-free (i.e. NFW) profile - for which r₋₂ is simply the scale radius - and are labeled so that they may be pruned from
    ! the tree prior to evolution.
    if (.not.cuspCanForm) then
       amplitude    =+0.0d0
       mass         =+0.0d0
       y            =+0.0d0
       radiusScale  =+radiusMinus2
       concentration=+self %darkMatterHaloScale_%radiusVirial(node) &
            &        /                           radiusScale
       densityScale =+basic%mass       ()                           &
            &        /4.0d0                                         &
            &        /Pi                                            &
            &        /radiusScale**3                                &
            &        /(                                             &
            &          +              log(1.0d0+concentration)      &
            &          -concentration/   (1.0d0+concentration)      &
            &         )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspAmplitudeID ,amplitude   )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspMassID      ,mass        )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID      ,y           )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWScaleID  ,radiusScale )
       call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWDensityID,densityScale)
       call nodeLabelSet                              (self%labelUnformedID        ,node        )
       ! If this is the base node of the tree then no halo in the tree could have formed. The base node can not be pruned, so warn
       ! that this tree will be reduced to an inert base node.
       if (.not.associated(node%parent)) then
          block
            type(varying_string) :: message
            message=var_str('no prompt cusp can form in the base node of a tree {node index = ')//node%index()//'}'//char(10)// &
                 &  '  no halo in this tree can have formed - the tree will be reduced to an inert base node'
            call Warn(message)
          end block
       end if
       return
    end if
    ! Check for convergence.
    if     (                                                                               &
         &   .not.Values_Agree(radiusScale ,radiusScalePrevious ,relTol=toleranceRelative) &
         &  .or.                                                                           &
         &   .not.Values_Agree(densityScale,densityScalePrevious,relTol=toleranceRelative) &
         &  .or.                                                                           &
         &   .not.Values_Agree(basic%mass(),massPrevious        ,relTol=toleranceRelative) &
         & ) then
       if (self%nonConvergenceIsFatal) then
          block
            character(len=12        ) :: label
            type     (varying_string) :: message
            write (label,'(e12.6)') toleranceRelative
            message=var_str('failed to converge when seeking solution for rₛ {target fractional error ')//label//'; node uniqueID = '//node%uniqueID()//'}'
            errorFractional=+abs(+radiusScale-radiusScalePrevious) &
                 &          /   (+radiusScale+radiusScalePrevious) &
                 &          /0.5d0
            write (label,'(e12.6)') errorFractional
            message=message//char(10)//'  Δrₛ/rₛ = '//label
            errorFractional=+abs(+densityScale-densityScalePrevious) &
                 &          /   (+densityScale+densityScalePrevious) &
                 &          /0.5d0
            write (label,'(e12.6)') errorFractional
            message=message//char(10)//'  Δρₛ/ρₛ = '//label
            errorFractional=+abs(+basic%mass()-massPrevious) &
                 &          /   (+basic%mass()+massPrevious) &
                 &          /0.5d0
            write (label,'(e12.6)') errorFractional
            message=message//char(10)//'  ΔM /M  = '//label
            call Error_Report(message//{introspection:location})
         end block
       else
          block
            character(len=12        ) :: label
            type     (varying_string) :: message
            errorFractional=+abs(+radiusScale-radiusScalePrevious) &
                 &          /   (+radiusScale+radiusScalePrevious) &
                 &          /0.5d0
            if (errorFractional > errorFractionalMaximum) then
               write (label,'(e12.6)') errorFractional
               message='failed to converge when seeking solution for rₛ (relative error = '//label//')'
               call Warn(message)
               errorFractionalMaximum=errorFractional
            end if
          end block
       end if
    end if    
    return
  end subroutine darkMatterProfilePromptCuspsNodeTreeInitialize

  function progenitorPrimary(labelUnformedID,node) result(nodeProgenitor)
    !!{RST
    Return a pointer to the primary progenitor which ``node`` will have once all halos identified as unable to form have been
    pruned from the merger tree, or a ``null()`` pointer if ``node`` will have no progenitors. Pruning unlinks such nodes without
    preserving primary progenitor status (i.e. ``firstChild`` is replaced by ``sibling``), so the primary progenitor after pruning
    is simply the first progenitor which is not labeled as being unable to form.
    !!}
    use :: Nodes_Labels, only : nodeLabelIsPresent
    implicit none
    type   (treeNode), pointer                :: nodeProgenitor
    integer          , intent(in   )          :: labelUnformedID
    type   (treeNode), intent(inout), target  :: node

    nodeProgenitor => node%firstChild
    do while (associated(nodeProgenitor))
       if (.not.nodeLabelIsPresent(labelUnformedID,nodeProgenitor)) return
       nodeProgenitor => nodeProgenitor%sibling
    end do
    return
  end function progenitorPrimary

  double precision function factorMassCuspNFW(concentration,y)
    !!{RST
    Evaluate the dimensionless factor :math:`M/4 \pi \rho_\mathrm{s} r_\mathrm{s}^3` giving the mass enclosed within a radius
    :math:`c r_\mathrm{s}` of a cusp-NFW density profile of parameter :math:`y`. For :math:`y` below ``yMinimum`` the NFW limit
    of this expression is used instead---see the discussion accompanying that parameter.
    !!}
    implicit none
    double precision, intent(in   ) :: concentration, y

    if (y > yMinimum) then
       ! Cusp-NFW case.
       factorMassCuspNFW=+ 2.0d0                       *asinh(sqrt(concentration            )/               y     ) &
            &            -(2.0d0-y**2)/sqrt(1.0d0-y**2)*atanh(sqrt(concentration*(1.0d0-y**2)/(concentration+y**2))) &
            &            -sqrt(concentration*(concentration+y**2))/(1.0d0+concentration)
    else
       ! NFW limit.
       factorMassCuspNFW=+              log(1.0d0+concentration) &
            &            -concentration/   (1.0d0+concentration)
    end if
    return
  end function factorMassCuspNFW

  double precision function timeCollapseRoot(timeCollapse)
    !!{RST
    Root function used to find the time of collapse.
    !!}
    implicit none
    double precision, intent(in   ) :: timeCollapse

    timeCollapseRoot=+self_%sigma         (0,timeCollapse) &
         &           -      sigma0Collapse
    return
  end function timeCollapseRoot

  double precision function concentrationTargetRoot(concentration)
    !!{RST
    Implements the equation in footnote 9 of :cite:t:`delos_cusp-halo_2025`. Used in solving for the minimum allowed concentration in a cusp-NFW density profile.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: concentration

    concentrationTargetRoot=concentrationFactorTarget-384.0d0*Pi*(asinh(sqrt(concentration/2.0d0))-sqrt(concentration/(2.0d0+concentration)))**2/concentration**3
    return
  end function concentrationTargetRoot
  
  double precision function darkMatterProfilePromptCuspsNodeSigma(self,j,time) result(sigma)
    !!{RST
    Evaluate the integral

    .. math::

       \sigma_j^2(t) = \int_0^\infty \frac{\mathrm{d}k}{k} \mathcal{P}(k,t) k^{2j},

    where :math:`\mathcal{P}(k) = k^3 P(k) / 2 \pi^2` is the dimensionless form of the power spectrum.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Comparison , only : Values_Agree
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout), target       :: self
    integer                                                   , intent(in   )               :: j
    double precision                                          , intent(in   )               :: time
    double precision                                          , parameter                   :: wavenumberPhysicalMinimum           =1.0d-2 , wavenumberPhysicalMaximum           =1.0d+2, &
         &                                                                                     toleranceRelative                   =1.0d-3
    type            (integrator                              ), save                        :: integrator_
    logical                                                   , save                        :: integratorInitialized               =.false.
    !$omp threadprivate(integrator_,integratorInitialized)
    double precision                                          , allocatable  , dimension(:) :: sigmaTmp
    double precision                                                                        :: wavenumberPhysicalLogarithmicMinimum        , wavenumberPhysicalLogarithmicMaximum       , &
         &                                                                                     sigmaPrevious

    ! Initialize an integrator if necessary.
    if (.not.integratorInitialized) then       
       integrator_          =integrator(integrand,toleranceRelative=toleranceRelative)
       integratorInitialized=.true.
    end if
    ! Is growth is not wavenumber dependent, check if we have computed the value of sigma for this j.
    if (.not.self%growthIsWavenumberDependent) then
       if (.not.allocated(self%sigma_)) then
          allocate(self%sigma_(0:j))
          self%sigma_                    =-huge(0.0d0)
       else if (ubound(self%sigma_,dim=1) < j) then
          call move_alloc(self%sigma_,sigmaTmp)
          allocate(self%sigma_(0:j))
          self%sigma_                    =-huge(0.0d0)
          self%sigma_(0:size(sigmaTmp)-1)= sigmaTmp
       end if
    end if
    ! Determine if we need to perform the integral.
    if (self%growthIsWavenumberDependent .or. self%sigma_(j) < 0.0d0) then
       ! Make an initial guess for the range of wavenumbers over which to integrate the power spectrum.
       if (self%growthIsWavenumberDependent) then
          time_                             =  time
          expansionFactor_                  =  self%cosmologyFunctions_%expansionFactor(time            )
       else
          expansionFactor_                  =  1.0d0
          time_                             =  self%cosmologyFunctions_%cosmicTime     (expansionFactor_)
       end if
       self_                                => self
       j_                                   =  j
       wavenumberPhysicalLogarithmicMinimum =  log(wavenumberPhysicalMinimum)
       wavenumberPhysicalLogarithmicMaximum =  log(wavenumberPhysicalMaximum)
       sigmaPrevious                        =  -huge(0.0d0)
       sigma                                =  +     0.0d0
       ! Expand the range over wavenumbers integrated over until the integral is sufficiently well converged.
       do while (.not.Values_Agree(sigma,sigmaPrevious,relTol=toleranceRelative))
          sigmaPrevious                       =+sigma
          sigma                               =+sqrt(                                                            &
               &                                     integrator_%integrate(                                      &
               &                                                           wavenumberPhysicalLogarithmicMinimum, &
               &                                                           wavenumberPhysicalLogarithmicMaximum  &
               &                                                          )                                      &
               &                                    )
          wavenumberPhysicalLogarithmicMinimum=+wavenumberPhysicalLogarithmicMinimum-1.0d0
          wavenumberPhysicalLogarithmicMaximum=+wavenumberPhysicalLogarithmicMaximum+1.0d0
       end do
       if (.not.self%growthIsWavenumberDependent) self%sigma_(j)=sigma
    end if
    ! If growth is wavenumber independent simply scale sigma to the current time (accounting for both linear growth and the fact
    ! that these σ are defined in physical, not comoving coordinates).
    if (.not.self%growthIsWavenumberDependent)                           &
         & sigma=+self                    %sigma_         (     j   )    &
         &       *self%linearGrowth_      %value          (time=time)    &
         &       /self%cosmologyFunctions_%expansionFactor(time=time)**j
    return
  end function darkMatterProfilePromptCuspsNodeSigma
  
  double precision function integrand(wavenumberPhysicalLogarithmic)
    !!{RST
    Integrand used to compute the quantity :math:`\sigma^2_j(t)` in the prompt cusps model of :cite:t:`delos_cusp-halo_2025`.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumberPhysicalLogarithmic
    double precision                :: wavenumberComoving           , wavenumberPhysical
    
    wavenumberPhysical=+exp(wavenumberPhysicalLogarithmic)
    wavenumberComoving=+expansionFactor_   &
         &             *wavenumberPhysical
    integrand         =+self_%powerSpectrum_%powerDimensionless(wavenumberComoving        ,time_) &
         &             *                                        wavenumberPhysical**(2*j_)
    return
  end function integrand

  subroutine darkMatterProfilePromptCuspsNodeInitialize(self,node)
    !!{RST
    Compute the rate of growth of dark matter profile scale radius assuming a constant growth rate.
    !!}
    use :: Display         , only : displayBlue       , displayGreen                  , displayYellow, displayBold, &
         &                          displayReset
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    class           (nodeComponentBasic                      )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile          )               , pointer :: darkMatterProfile, darkMatterProfileParent
    double precision                                                                   :: timeInterval

    ! Set the growth rate for the scale radius.
    darkMatterProfile => node%darkMatterProfile()
    select type (darkMatterProfile)
    type is (nodeComponentDarkMatterProfile)
       call Error_Report(                                                                                                                                                                                            &
            &            displayBold()//'darkMatterProfile'//displayReset()//' component must be created prior to initialization cusp-NFW scale radius interpolation'                                   //char(10)// &
            &            '   For example, by using the following nodeOperator'                                                                                                                          //char(10)// &
            &            '    <'//displayBlue()//'nodeOperator'//displayReset()//' '//displayYellow()//'value'//displayReset()//'='//displayGreen()//'"darkMatterProfileScaleSet"'//displayReset()//'/>'          // &
            &            {introspection:location}                                                                                                                                                                    &
            &           )
    class default
       if (node%isPrimaryProgenitor()) then
          ! Node is the primary progenitor, so compute the scale radius growth rate.
          basic        =>  node              %basic()
          basicParent  =>  node       %parent%basic()
          timeInterval =  +basicParent       %time () &
               &          -basic             %time ()
          if (timeInterval > 0.0d0) then
             darkMatterProfileParent => node%parent%darkMatterProfile()
             call darkMatterProfile%floatRank0MetaPropertySet(                                                                                 &
                  &                                              self                   %promptCuspNFWGrowthRateID                           , &
                  &                                           +(                                                                               &
                  &                                             +darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWScaleID)  &
                  &                                             -darkMatterProfile      %floatRank0MetaPropertyGet(self%promptCuspNFWScaleID)  &
                  &                                            )                                                                               &
                  &                                           /                          timeInterval                                          &
                  &                                          )
          else
             ! Time interval is non-positive - assume zero growth rate.
             call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWGrowthRateID,0.0d0)
          end if
       else
          ! Node is a non-primary progenitor - assume zero growth rate.
          call    darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWGrowthRateID,0.0d0)
       end if
    end select
    return
  end subroutine darkMatterProfilePromptCuspsNodeInitialize

  subroutine darkMatterProfilePromptCuspsNodePromote(self,node)
    !!{RST
    Ensure that ``node`` is ready for promotion to its parent. In this case, we simply update the scale radius growth rate of ``node`` to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfilePromptCusps), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile, darkMatterProfileParent

    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    ! Note that the cusp amplitude and mass must be copied along with the remaining cusp properties. Along any branch which
    ! remains after halos unable to form have been pruned these are constant, so copying them is a no-op. They are not constant if
    ! such halos are left in the tree, however, in which case failing to copy them leaves the promoted node with the zero
    ! amplitude of an unformed halo alongside the non-zero y-parameter of its parent - a state which is not self-consistent.
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspAmplitudeID    ,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspAmplitudeID    ))
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspMassID         ,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspMassID         ))
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID         ,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWYID         ))
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWGrowthRateID,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWGrowthRateID))
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWScaleID     ,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWScaleID     ))
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWDensityID   ,darkMatterProfileParent%floatRank0MetaPropertyGet(self%promptCuspNFWDensityID   ))
    return
  end subroutine darkMatterProfilePromptCuspsNodePromote

  subroutine darkMatterProfilePromptCuspsSolveAnalytics(self,node,time)
    !!{RST
    Compute the value of the :math:`y`-parameter in the prompt cusp.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rootFinder        , rangeExpandMultiplicative     , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    implicit none
    class           (nodeOperatorDarkMatterProfilePromptCusps), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    double precision                                          , intent(in   ) :: time
    class           (nodeComponentBasic                      ), pointer       :: basic                    , basicParent
    class           (nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile        , darkMatterProfileParent
    type            (rootFinder                              ), save          :: finder
    logical                                                   , save          :: finderInitialized=.false.
    !$omp threadprivate(finder,finderInitialized)
    double precision                                                          :: densityScale             , y

    ! Primary progenitors do not evolve.
    if (.not.node%isPrimaryProgenitor()) return
    ! Initialize the root finder.
    if (.not.finderInitialized) then
       finder=rootFinder(                                                             &
            &            rootFunction                 =densityNormalizationRoot     , &
            &            toleranceRelative            =1.0d-6                       , &
            &            rangeExpandUpward            =2.0d+0                       , &
            &            rangeExpandDownward          =0.5d+0                       , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
            &           )
      finderInitialized=.true.
    end if
    ! Solve for the density normalization.
    basic                   =>  node                                        %basic                    (                              )
    basicParent             =>  node                   %parent              %basic                    (                              )
    darkMatterProfile       =>  node                                        %darkMatterProfile        (                              )
    darkMatterProfileParent =>  node                   %parent              %darkMatterProfile        (                              )
    massHalo_               =   basic                                       %mass                     (                              )
    amplitudeCusp_          =  +darkMatterProfile                           %floatRank0MetaPropertyGet(self%promptCuspAmplitudeID    )
    y                       =  +darkMatterProfile                           %floatRank0MetaPropertyGet(self%promptCuspNFWYID         )
    radiusScale_            =  +darkMatterProfileParent                     %floatRank0MetaPropertyGet(self%promptCuspNFWScaleID     ) &
         &                     +(                                                                                                      &
         &                       +                                           time                                                      &
         &                       -basicParent                               %time                     (                              ) &
         &                      )                                                                                                      &
         &                     *darkMatterProfile                           %floatRank0MetaPropertyGet(self%promptCuspNFWGrowthRateID)
    concentration_          =  +self                   %darkMatterHaloScale_%radiusVirial             (node                          ) &
         &                     /                                             radiusScale_
    if (amplitudeCusp_ > 0.0d0) then
       ! Cusp-NFW (y>0) case - the density normalization must be found by solving for a self-consistent y.
       densityScale         =  +massHalo_                           &
            &                  /radiusScale_**3                     &
            &                  /4.0d0                               &
            &                  /Pi                                  &
            &                  /factorMassCuspNFW(concentration_,y)
       densityScale         =   finder%find(rootGuess=densityScale)
       ! Compute the cusp y parameter, treating negligible values as zero - see the discussion accompanying `yMinimum`.
       y                    =min(                      &
            &                    +amplitudeCusp_       &
            &                    /densityScale         &
            &                    /radiusScale_**1.5d0, &
            &                    +yMaximum             &
            &                   )
       if (y < yMinimum) y=0.0d0
    else
       ! NFW (y=0) case - the profile has no cusp, so the density normalization follows directly and no root find is needed.
       densityScale         =  +massHalo_                               &
            &                  /4.0d0                                   &
            &                  /Pi                                      &
            &                  /radiusScale_**3                         &
            &                  /factorMassCuspNFW(concentration_,0.0d0)
       y                    =  +0.0d0
    end if
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWYID      ,y           )
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWScaleID  ,radiusScale_)    
    call darkMatterProfile%floatRank0MetaPropertySet(self%promptCuspNFWDensityID,densityScale)    
    return
  end subroutine darkMatterProfilePromptCuspsSolveAnalytics

  double precision function densityNormalizationRoot(densityScale)
    !!{RST
    Root function used in finding the density normalization for cusp-NFW density profiles.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: densityScale
    double precision                :: y

    y                 =min(                      &
         &                 +amplitudeCusp_       &
         &                 /densityScale         &
         &                 /radiusScale_**1.5d0, &
         &                 +yMaximum             &
         &                )
    densityNormalizationRoot=+massHalo_                           &
         &                   /radiusScale_**3                     &
         &                   /4.0d0                               &
         &                   /Pi                                  &
         &                   /factorMassCuspNFW(concentration_,y) &
         &                   -densityScale
    return
  end function densityNormalizationRoot

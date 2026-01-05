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
  An implementation of accretion from the \gls{igm} onto halos using simple truncation to
  mimic the effects of reionization.
  !!}

  use :: Accretion_Halo_Totals     , only : accretionHaloTotal                     , accretionHaloTotalClass
  use :: Chemical_States           , only : chemicalState                          , chemicalStateClass
  use :: Cosmology_Functions       , only : cosmologyFunctions                     , cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParameters                    , cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScale                    , darkMatterHaloScaleClass
  use :: Intergalactic_Medium_State, only : intergalacticMediumState               , intergalacticMediumStateClass
  use :: Radiation_Fields          , only : radiationFieldCosmicMicrowaveBackground

  !![
  <accretionHalo name="accretionHaloSimple">
   <description>
    Accretion onto halos using simple truncation to mimic the effects of reionization. The accretion rate of baryons into a
    halo is given by:
    \begin{equation}
     \dot{M}_\mathrm{accretion} = \left\{ \begin{array}{ll} (\Omega_\mathrm{b}/\Omega_\mathrm{M}) \dot{M}_\mathrm{halo} &amp;
     \hbox{ if } V_\mathrm{virial} &gt; V_\mathrm{reionization} \hbox{ or } z &gt; z_\mathrm{reionization} \\ 0 &amp; \hbox{
     otherwise,}\end{array} \right.
    \label{eq:accretionHalo:truncation}
    \end{equation}
    where $z_\mathrm{reionization}=${\normalfont \ttfamily [redshiftReionization]} is the redshift at which the Universe is
    reionized (alternatively, the optical depth to reionization can be specified via {\normalfont \ttfamily
    [opticalDepthReionization]} and the corresponding redshift will be computed) and $V_\mathrm{reionization}=${\normalfont
    \ttfamily [velocitySuppressionReionization]} is the virial velocity below which accretion is suppressed after
    reionization. Setting $V_\mathrm{reionization}$ to zero will effectively switch off the effects of reionization on the
    accretion of baryons. This algorithm attempts to offer a simple prescription for the effects of reionization and has been
    explored by multiple authors (e.g. \citealt{benson_effects_2002}). In particular, \cite{font_modelingmilky_2010} show that
    it produces results in good agreement with more elaborate treatments of reionization. For halos below the accretion
    threshold, any accretion rate that would have otherwise occurred is instead placed into the ``failed'' accretion rate. For
    halos which can accrete, and which have some mass in their ``failed'' reservoir, that mass will be added to the regular
    accretion rate at a rate equal to the mass of the ``failed'' reservoir times the specific growth rate of the halo. The gas
    accreted is assumed to be from a pristine \gls{igm} and so has zero abundances. Chemical abundances are computed from the
    chemical state functions (see \refPhysics{chemicalState}).
  
    Note that, if $\dot{M}_\mathrm{halo} &lt; 0$ then negative accretion rates of gas into the node can result. This can be
    prevented by setting {\normalfont \ttfamily [accretionNegativeAllowed]}$=${\normalfont \ttfamily false}.
  
    By default, gas is accreted whenever the halo is growing in total mass. However, setting {\normalfont \ttfamily
    [accretionNewGrowthOnly]}$=${\normalfont \ttfamily true} causes accretion to occur only if the node mass is growing and
    exceeds the previous maximum node mass achieved along this branch of the merger tree. This requires use of a basic
    component which tracks the maximum mass along the branch (i.e. the {\normalfont \ttfamily massMaximum} property).
   </description>
   <deepCopy>
    <functionClass variables="radiation"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="radiation"/>
   </stateStorable>
  </accretionHalo>
  !!]
  type, extends(accretionHaloClass) :: accretionHaloSimple
     !!{
     A halo accretion class using simple truncation to mimic the effects of reionization.
     !!}
     private
     class           (cosmologyParametersClass               ), pointer :: cosmologyParameters_      => null()
     class           (accretionHaloTotalClass                ), pointer :: accretionHaloTotal_       => null()
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_       => null()
     class           (darkMatterHaloScaleClass               ), pointer :: darkMatterHaloScale_      => null()
     class           (intergalacticMediumStateClass          ), pointer :: intergalacticMediumState_ => null()
     class           (chemicalStateClass                     ), pointer :: chemicalState_            => null()
     double precision                                                   :: timeReionization                   , velocitySuppressionReionization, &
          &                                                                opticalDepthReionization           , redshiftReionization
     logical                                                            :: accretionNegativeAllowed           , accretionNewGrowthOnly
     type            (radiationFieldCosmicMicrowaveBackground), pointer :: radiation                 => null()
     integer                                                            :: countChemicals                     , massProgenitorMaximumID        , &
          &                                                                abundancePatternSolar
   contains
     !![
     <methods>
       <method description="Returns the fraction of potential accretion onto a halo from the \gls{igm} which fails." method="failedFraction"/>
       <method description="Returns the velocity scale to use for {\normalfont \ttfamily node}."                     method="velocityScale" />
       <method description="Compute masses of chemical species given a total mass."                                  method="chemicalMasses"/>
     </methods>
     !!]
     final     ::                              simpleDestructor
     procedure :: branchHasBaryons          => simpleBranchHasBaryons
     procedure :: accretionRate             => simpleAccretionRate
     procedure :: accretedMass              => simpleAccretedMass
     procedure :: failedAccretionRate       => simpleFailedAccretionRate
     procedure :: failedAccretedMass        => simpleFailedAccretedMass
     procedure :: accretionRateMetals       => simpleAccretionRateMetals
     procedure :: accretedMassMetals        => simpleAccretedMassMetals
     procedure :: failedAccretionRateMetals => simpleFailedAccretionRateMetals
     procedure :: failedAccretedMassMetals  => simpleFailedAccretedMassMetals
     procedure :: accretionRateChemicals    => simpleAccretionRateChemicals
     procedure :: accretedMassChemicals     => simpleAccretedMassChemicals
     procedure :: velocityScale             => simpleVelocityScale
     procedure :: failedFraction            => simpleFailedFraction
     procedure :: chemicalMasses            => simpleChemicalMasses
  end type accretionHaloSimple

  interface accretionHaloSimple
     !!{
     Constructors for the \refClass{accretionHaloSimple} halo accretion class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface accretionHaloSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily simple} halo accretion class.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (accretionHaloSimple          )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (accretionHaloTotalClass      ), pointer       :: accretionHaloTotal_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    class           (intergalacticMediumStateClass), pointer       :: intergalacticMediumState_
    class           (chemicalStateClass           ), pointer       :: chemicalState_
    double precision                                               :: timeReionization         , velocitySuppressionReionization, &
         &                                                            opticalDepthReionization , redshiftReionization            
    logical                                                        :: accretionNegativeAllowed , accretionNewGrowthOnly

    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="intergalacticMediumState" name="intergalacticMediumState_" source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="accretionHaloTotal"       name="accretionHaloTotal_"       source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    <objectBuilder class="chemicalState"            name="chemicalState_"            source="parameters"/>
    !!]
    if (parameters%isPresent("redshiftReionization").or..not.parameters%isPresent("opticalDepthReionization")) then
       if (parameters%isPresent("opticalDepthReionization")) call Error_Report("only one of [opticalDepthReionization] and [redshiftReionization] should be specified"//{introspection:location})
       !![
       <inputParameter>
         <name>redshiftReionization</name>
         <defaultSource>(\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultSource>
         <defaultValue>9.97d0</defaultValue>
         <description>The redshift below which baryonic accretion is suppressed.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       timeReionization=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftReionization))
    else
       !![
       <inputParameter>
         <name>opticalDepthReionization</name>
         <description>The optical depth to electron scattering below which baryonic accretion is suppressed.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
       timeReionization=intergalacticMediumState_%electronScatteringTime(opticalDepthReionization,assumeFullyIonized=.true.)
    end if
    !![
    <inputParameter>
      <name>velocitySuppressionReionization</name>
      <defaultValue>35.0d0</defaultValue>
      <description>The velocity scale below which baryonic accretion is suppressed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>accretionNegativeAllowed</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether negative accretion (mass loss) is allowed in the simple halo accretion model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>accretionNewGrowthOnly</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether accretion from the \gls{igm} is allowed only when a halo is growing past its previous greatest mass.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=accretionHaloSimple(timeReionization,velocitySuppressionReionization,accretionNegativeAllowed,accretionNewGrowthOnly,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,accretionHaloTotal_,chemicalState_,intergalacticMediumState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="intergalacticMediumState_"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="accretionHaloTotal_"      />
    <objectDestructor name="darkMatterHaloScale_"     />
    <objectDestructor name="chemicalState_"           />
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(timeReionization,velocitySuppressionReionization,accretionNegativeAllowed,accretionNewGrowthOnly,cosmologyParameters_,cosmologyFunctions_,darkMatterHaloScale_,accretionHaloTotal_,chemicalState_,intergalacticMediumState_) result(self)
    !!{
    Internal constructor for the \refClass{accretionHaloSimple} halo accretion class.
    !!}
    use :: Atomic_Data                  , only : Abundance_Pattern_Lookup
    use :: Chemical_Abundances_Structure, only : Chemicals_Property_Count
    implicit none
    type            (accretionHaloSimple          ), target                :: self
    double precision                               , intent(in   )         :: timeReionization        , velocitySuppressionReionization
    logical                                        , intent(in   )         :: accretionNegativeAllowed, accretionNewGrowthOnly
    class           (cosmologyParametersClass     ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    class           (accretionHaloTotalClass      ), intent(in   ), target :: accretionHaloTotal_
    class           (darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    class           (chemicalStateClass           ), intent(in   ), target :: chemicalState_
    class           (intergalacticMediumStateClass), intent(in   ), target :: intergalacticMediumState_
    !![
    <constructorAssign variables="timeReionization, velocitySuppressionReionization, accretionNegativeAllowed, accretionNewGrowthOnly, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterHaloScale_, *accretionHaloTotal_, *chemicalState_, *intergalacticMediumState_"/>
    !!]

    allocate(self%radiation)
    !![
    <referenceConstruct isResult="yes" owner="self" object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
    !!]
    if (accretionNewGrowthOnly) then
       !![
       <addMetaProperty component="basic" name="massProgenitorMaximum" id="self%massProgenitorMaximumID" isEvolvable="no"/>
       !!]
    end if
    self%abundancePatternSolar   =Abundance_Pattern_Lookup(abundanceName='solar')
    self%countChemicals          =Chemicals_Property_Count()
    self%redshiftReionization    =self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeReionization))
    self%opticalDepthReionization=-huge(0.0d0)
    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Destructor for the \refClass{accretionHaloSimple} halo accretion class.
    !!}
    implicit none
    type(accretionHaloSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%darkMatterHaloScale_"     />
    <objectDestructor name="self%accretionHaloTotal_"      />
    <objectDestructor name="self%chemicalState_"           />
    <objectDestructor name="self%intergalacticMediumState_"/>
    <objectDestructor name="self%radiation"                />
    !!]
    return
  end subroutine simpleDestructor

  logical function simpleBranchHasBaryons(self,node)
    !!{
    Returns true if this branch can accrete any baryons.
    !!}
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodesBranch
    implicit none
    class(accretionHaloSimple                ), intent(inout)          :: self
    type (treeNode                           ), intent(inout), target  :: node
    type (treeNode                           )               , pointer :: branchNode
    type (mergerTreeWalkerIsolatedNodesBranch)                         :: treeWalker

    simpleBranchHasBaryons=.false.
    treeWalker             =mergerTreeWalkerIsolatedNodesBranch(node)
    do while (treeWalker%next(branchNode))
       if (self%failedFraction(branchNode) < 1.0d0) then
          simpleBranchHasBaryons=.true.
          exit
       end if
    end do
    return
  end function simpleBranchHasBaryons

  double precision function simpleAccretionRate(self,node,accretionMode)
    !!{
    Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentHotHalo, treeNode
    implicit none
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    class           (nodeComponentHotHalo        ), pointer       :: hotHalo
    double precision                                              :: growthRate    , unaccretedMass, &
         &                                                           failedFraction

    simpleAccretionRate=0.0d0
    if (accretionMode      == accretionModeCold) return
    if (node%isSatellite()                     ) return
    ! Get the failed accretion fraction.
    failedFraction=self%failedFraction(node)
    ! Get the default cosmology.
    basic                => node%basic         ()
    hotHalo              => node%hotHalo       ()
    simpleAccretionRate=(self%cosmologyParameters_%OmegaBaryon()/self%cosmologyParameters_%OmegaMatter())*self%accretionHaloTotal_%accretionRate(node)*(1.0d0-failedFraction)
    ! Test for negative accretion.
    if (.not.self%accretionNegativeAllowed.and.self%accretionHaloTotal_%accretionRate(node) < 0.0d0) then
       ! Accretion rate is negative, and not allowed. Return zero accretion rate.
       simpleAccretionRate=0.0d0
    else
       ! Return the standard accretion rate.
       unaccretedMass=hotHalo%unaccretedMass()
       growthRate=self%accretionHaloTotal_%accretionRate(node)/self%accretionHaloTotal_%accretedMass(node)
       simpleAccretionRate=simpleAccretionRate+unaccretedMass*growthRate*(1.0d0-failedFraction)
    end if
    ! If accretion is allowed only on new growth, check for new growth and shut off accretion if growth is not new.
    if (self%accretionNewGrowthOnly) then
       if (self%accretionHaloTotal_%accretedMass(node) < basic%floatRank0MetaPropertyGet(self%massProgenitorMaximumID)) simpleAccretionRate=0.0d0
    end if
    return
  end function simpleAccretionRate

  double precision function simpleAccretedMass(self,node,accretionMode)
    !!{
    Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    double precision                                              :: failedFraction

    simpleAccretedMass=0.0d0
    if (accretionMode      == accretionModeCold) return
    if (node%isSatellite()                     ) return
    ! Get the failed accretion fraction.
    failedFraction     =  self%failedFraction(node)
    basic              => node%basic         (    )
    simpleAccretedMass =  (self%cosmologyParameters_%OmegaBaryon()/self%cosmologyParameters_%OmegaMatter())*self%accretionHaloTotal_%accretedMass(node)*(1.0d0-failedFraction)
    return
  end function simpleAccretedMass

  double precision function simpleFailedAccretionRate(self,node,accretionMode)
    !!{
    Computes the baryonic accretion rate onto {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentHotHalo, treeNode
    implicit none
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    class           (nodeComponentHotHalo        ), pointer       :: hotHalo
    double precision                                              :: growthRate    , unaccretedMass, &
         &                                                           failedFraction

    simpleFailedAccretionRate=0.0d0
    if (accretionMode               == accretionModeCold) return
    if (node         %isSatellite()                     ) return
    ! Get the failed fraction.
    failedFraction=self%failedFraction(node)
    ! Test for negative accretion.
    if (.not.self%accretionNegativeAllowed.and.self%accretionHaloTotal_%accretionRate(node) < 0.0d0) then
       simpleFailedAccretionRate=(self%cosmologyParameters_%OmegaBaryon()/self%cosmologyParameters_%OmegaMatter())*self%accretionHaloTotal_%accretionRate(node)
    else
       hotHalo => node%hotHalo()
       simpleFailedAccretionRate=(self%cosmologyParameters_%OmegaBaryon()/self%cosmologyParameters_%OmegaMatter())*self%accretionHaloTotal_%accretionRate(node)*failedFraction
       unaccretedMass=hotHalo%unaccretedMass()
       growthRate=self%accretionHaloTotal_%accretionRate(node)/self%accretionHaloTotal_%accretedMass(node)
       simpleFailedAccretionRate=simpleFailedAccretionRate-unaccretedMass*growthRate*(1.0d0-failedFraction)
    end if
    ! If accretion is allowed only on new growth, check for new growth and shut off accretion if growth is not new.
    if (self%accretionNewGrowthOnly) then
       basic => node%basic()
       if (self%accretionHaloTotal_%accretedMass(node) < basic%floatRank0MetaPropertyGet(self%massProgenitorMaximumID)) simpleFailedAccretionRate=0.0d0
    end if
    return
  end function simpleFailedAccretionRate

  double precision function simpleFailedAccretedMass(self,node,accretionMode)
    !!{
    Computes the mass of baryons accreted into {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    double precision                                              :: failedFraction

    simpleFailedAccretedMass=0.0d0
    if (accretionMode      == accretionModeCold) return
    if (node%isSatellite()                     ) return
    ! Get the failed fraction.
    failedFraction=self%failedFraction(node)
    basic => node%basic()
    simpleFailedAccretedMass=(self%cosmologyParameters_%OmegaBaryon()/self%cosmologyParameters_%OmegaMatter())*self%accretionHaloTotal_%accretedMass(node)*failedFraction
    return
  end function simpleFailedAccretedMass

  function simpleAccretionRateMetals(self,node,accretionMode)
    !!{
    Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : metallicityTypeLinearByMass, adjustElementsReset
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    implicit none
    type            (abundances                  )                :: simpleAccretionRateMetals
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic           ), pointer      :: basic 
    double precision                                              :: metallicityIGM
    !$GLC attributes unused :: accretionMode

    basic          => node%basic                                (            )
    metallicityIGM =  self%intergalacticMediumState_%metallicity(basic%time())
    call simpleAccretionRateMetals%metallicitySet(metallicityIGM,metallicityTypeLinearByMass,adjustElementsReset,self%abundancePatternSolar)
    simpleAccretionRateMetals=+     simpleAccretionRateMetals                     &
         &                    *self%accretionRate            (node,accretionMode)
    return
  end function simpleAccretionRateMetals

  function simpleAccretedMassMetals(self,node,accretionMode)
    !!{
    Computes the mass of abundances accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : metallicityTypeLinearByMass, adjustElementsReset
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    implicit none
    type            (abundances                  )                :: simpleAccretedMassMetals
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    double precision                                              :: metallicityIGM
    !$GLC attributes unused :: accretionMode

    basic          => node%basic                                (            )
    metallicityIGM =  self%intergalacticMediumState_%metallicity(basic%time())
    call simpleAccretedMassMetals%metallicitySet(metallicityIGM,metallicityTypeLinearByMass,adjustElementsReset,self%abundancePatternSolar)
    simpleAccretedMassMetals=+     simpleAccretedMassMetals                     &
         &                   *self%accretionRate           (node,accretionMode)
    return
  end function simpleAccretedMassMetals

  function simpleFailedAccretionRateMetals(self,node,accretionMode)
    !!{
    Computes the rate of failed mass of abundance accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : metallicityTypeLinearByMass, adjustElementsReset
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    implicit none
    type            (abundances                  )                :: simpleFailedAccretionRateMetals
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    double precision                                              :: metallicityIGM
    !$GLC attributes unused :: accretionMode

    basic          => node%basic                                (            )
    metallicityIGM =  self%intergalacticMediumState_%metallicity(basic%time())
    call simpleFailedAccretionRateMetals%metallicitySet(metallicityIGM,metallicityTypeLinearByMass,adjustElementsReset,self%abundancePatternSolar)
    simpleFailedAccretionRateMetals=+     simpleFailedAccretionRateMetals                     &
         &                          *self%failedAccretionRate            (node,accretionMode)
    return
  end function simpleFailedAccretionRateMetals

  function simpleFailedAccretedMassMetals(self,node,accretionMode)
    !!{
    Computes the mass of abundances that failed to accrete (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure, only : metallicityTypeLinearByMass, adjustElementsReset
    use :: Galacticus_Nodes    , only : nodeComponentBasic
    implicit none
    type            (abundances                  )                :: simpleFailedAccretedMassMetals
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    double precision                                              :: metallicityIGM
    !$GLC attributes unused :: accretionMode

    basic          => node%basic                                (            )
    metallicityIGM =  self%intergalacticMediumState_%metallicity(basic%time())
    call simpleFailedAccretedMassMetals%metallicitySet(metallicityIGM,metallicityTypeLinearByMass,adjustElementsReset,self%abundancePatternSolar)
    simpleFailedAccretedMassMetals=+     simpleFailedAccretedMassMetals                     &
         &                         *self%accretionRate                 (node,accretionMode)
    return
  end function simpleFailedAccretedMassMetals
  
  function simpleAccretionRateChemicals(self,node,accretionMode)
    !!{
    Computes the rate of mass of chemicals accretion (in $M_\odot/$Gyr) onto {\normalfont \ttfamily node} from the intergalactic medium. Assumes a
    primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    temperature.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    implicit none
    type            (chemicalAbundances          )                :: simpleAccretionRateChemicals
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    double precision                                              :: massAccretionRate

    ! Ensure that chemicals are reset to zero.
    call simpleAccretionRateChemicals%reset()
    ! Return immediately if no chemicals are being tracked.
    if (self%countChemicals == 0) return
    ! Get the total mass accretion rate onto the halo.
    massAccretionRate=self%accretionRate(node,accretionMode)
    ! Get the mass accretion rates.
    simpleAccretionRateChemicals=self%chemicalMasses(node,massAccretionRate,accretionMode)
    return
  end function simpleAccretionRateChemicals

  function simpleAccretedMassChemicals(self,node,accretionMode)
    !!{
    Computes the mass of chemicals accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    implicit none
    type            (chemicalAbundances          )                :: simpleAccretedMassChemicals
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    double precision                                              :: massAccreted

    ! Ensure that chemicals are reset to zero.
    call simpleAccretedMassChemicals%reset()
    ! Return if no chemicals are being tracked.
    if (self%countChemicals == 0) return
    ! Total mass of material accreted.
    massAccreted=self%accretedMass(node,accretionMode)
    ! Get the masses of chemicals accreted.
    simpleAccretedMassChemicals=self%chemicalMasses(node,massAccreted,accretionMode)
    return
  end function simpleAccretedMassChemicals

  function simpleChemicalMasses(self,node,massAccreted,accretionMode)
    !!{
    Compute the masses of chemicals accreted (in $M_\odot$) onto {\normalfont \ttfamily node} from the intergalactic medium.
    !!}
    use :: Abundances_Structure             , only : zeroAbundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Astronomical , only : hydrogenByMassPrimordial
    use :: Numerical_Constants_Atomic       , only : atomicMassHydrogen
    implicit none
    class           (accretionHaloSimple         ), intent(inout) :: self
    type            (chemicalAbundances          )                :: simpleChemicalMasses
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: massAccreted
    type            (enumerationAccretionModeType), intent(in   ) :: accretionMode
    class           (nodeComponentBasic          ), pointer       :: basic
    type            (chemicalAbundances          ), save          :: chemicalDensities
    !$omp threadprivate(chemicalDensities)
    double precision                                              :: massToDensityConversion, numberDensityHydrogen, &
         &                                                           temperature

    ! Compute coefficient in conversion of mass to density for this node.
    massToDensityConversion=Chemicals_Mass_To_Density_Conversion(self%darkMatterHaloScale_%radiusVirial(node))/3.0d0
    ! Compute the temperature and density of accreting material, assuming
    ! accreted has is at the virial temperature and that the
    ! overdensity is one third of the mean overdensity of the halo.
    temperature          = self%darkMatterHaloScale_%temperatureVirial(node)
    numberDensityHydrogen= hydrogenByMassPrimordial*(self%cosmologyParameters_%OmegaBaryon()/self%cosmologyParameters_%OmegaMatter())*self%accretionHaloTotal_%accretedMass(node)*massToDensityConversion/atomicMassHydrogen
    ! Set the radiation field.
    basic => node%basic()
    call self%radiation%timeSet(basic%time())
    ! Get the chemical densities.
    call self%chemicalState_%chemicalDensities(chemicalDensities,numberDensityHydrogen,temperature,zeroAbundances,self%radiation)
    ! Convert from densities to masses.
    call chemicalDensities   %numberToMass(simpleChemicalMasses                                                          )
    call simpleChemicalMasses%scale       (massAccreted*hydrogenByMassPrimordial/numberDensityHydrogen/atomicMassHydrogen)
    return
  end function simpleChemicalMasses

  double precision function simpleVelocityScale(self,node)
    !!{
    Returns the velocity scale to use for {\normalfont \ttfamily node}. Use the virial velocity.
    !!}
    implicit none
    class(accretionHaloSimple), intent(inout) :: self
    type (treeNode           ), intent(inout) :: node

    simpleVelocityScale=self%darkMatterHaloScale_%velocityVirial(node)
    return
  end function simpleVelocityScale

  double precision function simpleFailedFraction(self,node)
    !!{
    Returns the fraction of potential accretion onto a halo from the \gls{igm} which fails.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(accretionHaloSimple), intent(inout) :: self
    type (treeNode           ), intent(inout) :: node
    class(nodeComponentBasic ), pointer       :: basic

    basic => node%basic()
    if     (                                                                 &
         &  basic%time         (    ) > self%timeReionization                &
         &   .and.                                                           &
         &  self %velocityScale(node) < self%velocitySuppressionReionization &
         & ) then
       simpleFailedFraction=1.0d0
    else
       simpleFailedFraction=0.0d0
    end if
    return
  end function simpleFailedFraction

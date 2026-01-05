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
  An implementation of calculations of chemical reaction rates for hydrogen using the fits from \cite{abel_modeling_1997} and
  \cite{tegmark_small_1997}.
  !!}

  use :: Atomic_Cross_Sections_Ionization_Photo, only : atomicCrossSectionIonizationPhotoClass
  use :: Atomic_Rates_Ionization_Collisional   , only : atomicIonizationRateCollisionalClass
  use :: Atomic_Rates_Recombination_Radiative  , only : atomicRecombinationRateRadiativeClass
  use :: Radiation_Fields                      , only : crossSectionFunctionTemplate

  !![
  <chemicalReactionRate name="chemicalReactionRateHydrogenNetwork">
   <description>
    A chemical reaction rate class that computes rates using the network of reactions and fitting functions from
    \cite{abel_modeling_1997} and \cite{tegmark_small_1997}. The parameter {\normalfont \ttfamily [fast]}
    controls the approximations made. If set {\normalfont \ttfamily true} then H$^-$ is assumed to be at equilibrium abundance,
    H$_2^+$ reactions are ignored and other slow reactions are ignored (see \citealt{abel_modeling_1997}).
   </description>
  </chemicalReactionRate>
  !!]
  type, extends(chemicalReactionRateClass) :: chemicalReactionRateHydrogenNetwork
     !!{
     A chemical reaction rates for hydrogen using the fits from \cite{abel_modeling_1997} and \cite{tegmark_small_1997}.
     !!}
     private
     class           (atomicIonizationRateCollisionalClass  ), pointer :: atomicIonizationRateCollisional_   => null()
     class           (atomicRecombinationRateRadiativeClass ), pointer :: atomicRecombinationRateRadiative_  => null()
     class           (atomicCrossSectionIonizationPhotoClass), pointer :: atomicCrossSectionIonizationPhoto_ => null()
     logical                                                           :: fast                                        , includeSelfShielding
     integer                                                           :: atomicHydrogenAnionIndex                    , atomicHydrogenCationIndex  , &
          &                                                               atomicHydrogenIndex                         , electronIndex
     double precision                                                  :: densityAtomicHydrogenAnion
   contains
     !![
     <methods>
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{e}^- \rightarrow \hbox{H}^+ + 2\hbox{e}^-$."    method="rateH_Electron_to_Hplus_2Electron"  />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^+ + \hbox{e}^- \rightarrow \hbox{H} + \gamma$."         method="rateHplus_Electron_to_H_Photon"     />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{e}^- \rightarrow \hbox{H}^- + \gamma$."         method="rateH_Electron_to_Hminus_Photon"    />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{H}^- \rightarrow \hbox{H}_2 + \hbox{e}^-$."     method="rateH_Hminus_to_H2_Electron"        />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \gamma$."       method="rateH_Hplus_to_H2plus_Photon"       />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$."   method="rateH2plus_H_to_H2_Hplus"           />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \hbox{H}$."   method="rateH2_Hplus_to_H2plus_H"           />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{e}^- \rightarrow 2\hbox{H} + \hbox{e}^-$."    method="rateH2_Electron_to_2H_Electron"     />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{H} \rightarrow 3\hbox{H}$."                   method="rateH2_H_to_3H"                     />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$."   method="rateHminus_Electron_to_H_2Electron" />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \hbox{H} \rightarrow 2 \hbox{H} + \hbox{e}^-$."     method="rateHminus_H_to_2H_Electron"        />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$."   method="rateHminus_Hplus_to_2H"             />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \hbox{e}^-$." method="rateHminus_Hplus_to_H2plus_Electron"/>
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{e}^- \rightarrow 2\hbox{H}$."               method="rateH2plus_Electron_to_2H"          />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H}^- \rightarrow \hbox{H}_2 + \hbox{H}$."   method="rateH2plus_Hminus_to_H2_H"          />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \gamma \rightarrow \hbox{H}^+ + \hbox{e}^-$."         method="rateH_Gamma_to_Hplus_Electron"      />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \gamma \rightarrow \hbox{H} + \hbox{e}^-$."         method="rateHminus_Gamma_to_H_Electron"     />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow \hbox{H}_2^+ + \hbox{e}^-$."     method="rateH2_Gamma_to_H2plus_Electron"    />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow \hbox{H} + \hbox{H}^+$."       method="rateH2plus_Gamma_to_H_Hplus"        />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow 2\hbox{H}^+ + \hbox{e}^-$."    method="rateH2plus_Gamma_to_2Hplus_Electron"/>
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow 2\hbox{H}^+ + \hbox{e}^-$."    method="rateH2_Gamma_to_H2star_to_2H"       />
       <method description="Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow 2\hbox{H}$."                     method="rateH2_Gamma_to_2H"                 />
     </methods>
     !!]
     final     ::                                        hydrogenNetworkDestructor
     procedure :: rates                               => hydrogenNetworkRates
     procedure :: rateH_Electron_to_Hplus_2Electron   => hydrogenNetworkRateH_Electron_to_Hplus_2Electron
     procedure :: rateHplus_Electron_to_H_Photon      => hydrogenNetworkRateHplus_Electron_to_H_Photon
     procedure :: rateH_Electron_to_Hminus_Photon     => hydrogenNetworkRateH_Electron_to_Hminus_Photon
     procedure :: rateH_Hminus_to_H2_Electron         => hydrogenNetworkRateH_Hminus_to_H2_Electron
     procedure :: rateH_Hplus_to_H2plus_Photon        => hydrogenNetworkRateH_Hplus_to_H2plus_Photon
     procedure :: rateH2plus_H_to_H2_Hplus            => hydrogenNetworkRateH2plus_H_to_H2_Hplus
     procedure :: rateH2_Hplus_to_H2plus_H            => hydrogenNetworkRateH2_Hplus_to_H2plus_H
     procedure :: rateH2_Electron_to_2H_Electron      => hydrogenNetworkRateH2_Electron_to_2H_Electron
     procedure :: rateH2_H_to_3H                      => hydrogenNetworkRateH2_H_to_3H
     procedure :: rateHminus_Electron_to_H_2Electron  => hydrogenNetworkRateHminus_Electron_to_H_2Electron
     procedure :: rateHminus_H_to_2H_Electron         => hydrogenNetworkRateHminus_H_to_2H_Electron
     procedure :: rateHminus_Hplus_to_2H              => hydrogenNetworkRateHminus_Hplus_to_2H
     procedure :: rateHminus_Hplus_to_H2plus_Electron => hydrogenNetworkRateHminus_Hplus_to_H2plus_Electron
     procedure :: rateH2plus_Electron_to_2H           => hydrogenNetworkRateH2plus_Electron_to_2H
     procedure :: rateH2plus_Hminus_to_H2_H           => hydrogenNetworkRateH2plus_Hminus_to_H2_H
     procedure :: rateH_Gamma_to_Hplus_Electron       => hydrogenNetworkRateH_Gamma_to_Hplus_Electron
     procedure :: rateHminus_Gamma_to_H_Electron      => hydrogenNetworkRateHminus_Gamma_to_H_Electron
     procedure :: rateH2_Gamma_to_H2plus_Electron     => hydrogenNetworkRateH2_Gamma_to_H2plus_Electron
     procedure :: rateH2plus_Gamma_to_H_Hplus         => hydrogenNetworkRateH2plus_Gamma_to_H_Hplus
     procedure :: rateH2plus_Gamma_to_2Hplus_Electron => hydrogenNetworkRateH2plus_Gamma_to_2Hplus_Electron
     procedure :: rateH2_Gamma_to_H2star_to_2H        => hydrogenNetworkRateH2_Gamma_to_H2star_to_2H
     procedure :: rateH2_Gamma_to_2H                  => hydrogenNetworkRateH2_Gamma_to_2H
  end type chemicalReactionRateHydrogenNetwork

  interface chemicalReactionRateHydrogenNetwork
     !!{
     Constructors for the \refClass{chemicalReactionRateHydrogenNetwork} chemical reaction rates class.
     !!}
     module procedure hydrogenNetworkConstructorParameters
     module procedure hydrogenNetworkConstructorInternal
  end interface chemicalReactionRateHydrogenNetwork

  ! Module-scope pointer to self used in integrations.
  class(chemicalReactionRateHydrogenNetwork), pointer   :: self_
  !$omp threadprivate(self_)

  ! Pointers to cross-section functions.
  procedure(crossSectionFunctionTemplate   ), pointer   :: crossSection_H_Gamma_to_Hplus_Electron_       => crossSection_H_Gamma_to_Hplus_Electron
  procedure(crossSectionFunctionTemplate   ), pointer   :: crossSection_H2_Gamma_to_2H_                  => crossSection_H2_Gamma_to_2H
  procedure(crossSectionFunctionTemplate   ), pointer   :: crossSection_H2plus_Gamma_to_2Hplus_Electron_ => crossSection_H2plus_Gamma_to_2Hplus_Electron
  procedure(crossSectionFunctionTemplate   ), pointer   :: crossSection_H2_Gamma_to_H2plus_Electron_     => crossSection_H2_Gamma_to_H2plus_Electron
  procedure(crossSectionFunctionTemplate   ), pointer   :: crossSection_H2plus_Gamma_to_H_Hplus_         => crossSection_H2plus_Gamma_to_H_Hplus
  procedure(crossSectionFunctionTemplate   ), pointer   :: crossSection_Hminus_Gamma_to_H_Electron_      => crossSection_Hminus_Gamma_to_H_Electron
  
contains

  function hydrogenNetworkConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{chemicalReactionRateHydrogenNetwork} chemical reaction rates class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (chemicalReactionRateHydrogenNetwork   )                :: self
    type   (inputParameters                       ), intent(inout) :: parameters
    class  (atomicIonizationRateCollisionalClass  ), pointer       :: atomicIonizationRateCollisional_
    class  (atomicRecombinationRateRadiativeClass ), pointer       :: atomicRecombinationRateRadiative_
    class  (atomicCrossSectionIonizationPhotoClass), pointer       :: atomicCrossSectionIonizationPhoto_
    logical                                                        :: fast                              , includeSelfShielding

    !![
    <inputParameter>
      <name>fast</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not to use simplifying assumptions to speed the hydrogen network calculation. If true, H$^-$
        is assumed to be at equilibrium abundance, H$_2^+$ reactions are ignored and other slow reactions are ignored (see
        \citealt{abel_modeling_1997}).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeSelfShielding</name>
      <defaultValue>.false.</defaultValue>
      <description>
	If true, include estimates of self-shielding when computing reaction rates involving the radiation field.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="atomicIonizationRateCollisional"   name="atomicIonizationRateCollisional_"   source="parameters"/>
    <objectBuilder class="atomicRecombinationRateRadiative"  name="atomicRecombinationRateRadiative_"  source="parameters"/>
    <objectBuilder class="atomicCrossSectionIonizationPhoto" name="atomicCrossSectionIonizationPhoto_" source="parameters"/>
    !!]
    self=chemicalReactionRateHydrogenNetwork(fast,includeSelfShielding,atomicIonizationRateCollisional_,atomicRecombinationRateRadiative_,atomicCrossSectionIonizationPhoto_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="atomicIonizationRateCollisional_"  />
    <objectDestructor name="atomicRecombinationRateRadiative_" />
    <objectDestructor name="atomicCrossSectionIonizationPhoto_"/>
    !!]
    return
  end function hydrogenNetworkConstructorParameters

  function hydrogenNetworkConstructorInternal(fast,includeSelfShielding,atomicIonizationRateCollisional_,atomicRecombinationRateRadiative_,atomicCrossSectionIonizationPhoto_) result(self)
    !!{
    Constructor for the \refClass{chemicalReactionRateHydrogenNetwork} chemical reaction rates class which takes a parameter set as
    input.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    implicit none
    type   (chemicalReactionRateHydrogenNetwork   )                        :: self
    logical                                        , intent(in   )         :: fast                              , includeSelfShielding
    class  (atomicIonizationRateCollisionalClass  ), intent(in   ), target :: atomicIonizationRateCollisional_
    class  (atomicRecombinationRateRadiativeClass ), intent(in   ), target :: atomicRecombinationRateRadiative_
    class  (atomicCrossSectionIonizationPhotoClass), intent(in   ), target :: atomicCrossSectionIonizationPhoto_
    !![
    <constructorAssign variables="fast, includeSelfShielding, *atomicIonizationRateCollisional_, *atomicRecombinationRateRadiative_, *atomicCrossSectionIonizationPhoto_"/>
    !!]

    if (self%fast) then
       ! Get indices needed for equilibrium calculation of H⁻.
       self%atomicHydrogenIndex      =Chemicals_Index("AtomicHydrogen"      )
       self%atomicHydrogenCationIndex=Chemicals_Index("AtomicHydrogenCation")
       self%electronIndex            =Chemicals_Index("Electron"            )
    else
       ! Get actual hydrogen anion index.
       self%atomicHydrogenAnionIndex=Chemicals_Index("AtomicHydrogenAnion"  )
    end if
    return
  end function hydrogenNetworkConstructorInternal

  subroutine hydrogenNetworkDestructor(self)
    !!{
    Destructor for the \refClass{chemicalReactionRateHydrogenNetwork} chemical reaction rates class.
    !!}
    implicit none
    type(chemicalReactionRateHydrogenNetwork), intent(inout) :: self

    !![
    <objectDestructor name="self%atomicIonizationRateCollisional_"  />
    <objectDestructor name="self%atomicRecombinationRateRadiative_" />
    <objectDestructor name="self%atomicCrossSectionIonizationPhoto_"/>
    !!]
    return
  end subroutine hydrogenNetworkDestructor

  subroutine hydrogenNetworkRates(self,lengthColumn,temperature,chemicalDensity,factorClumping,radiation,chemicalRates,node)
    !!{
    Compute rates of change of chemical abundances due to reactions involving chemical hydrogen species.
    !!}
    use :: Error           , only : Error_Report
    use :: Radiation_Fields, only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout), target :: self
    type            (chemicalAbundances                 ), intent(in   )         :: chemicalDensity
    double precision                                     , intent(in   )         :: temperature    , lengthColumn   , &
         &                                                                          factorClumping
    class           (radiationFieldClass                ), intent(inout)         :: radiation
    type            (chemicalAbundances                 ), intent(inout)         :: chemicalRates
    type            (treeNode                           ), intent(inout)         :: node
    double precision                                                             :: creationTerm   , destructionTerm

    ! Reset rates to zero initially.
    call chemicalRates%reset()
    ! Determine the atomic hydrogen anion density to use.
    if (self%fast) then
       ! For the fast network, assume the hydrogen anion is always at equilibrium density. (Following eqn. 24 of Abel et
       ! al. 1997.) No need to account for clumping factors here as all processes are proportional to density squared.
       creationTerm   =+hydrogenNetworkH_Electron_to_Hminus_Photon_RateCoefficient   (temperature)*chemicalDensity%abundance(self%atomicHydrogenIndex      ) &
            &                                                                                     *chemicalDensity%abundance(self%electronIndex            )
       destructionTerm=+hydrogenNetworkH_Hminus_to_H2_Electron_RateCoefficient       (temperature)*chemicalDensity%abundance(self%atomicHydrogenIndex      ) &
            &          +hydrogenNetworkHminus_Hplus_to_2H_RateCoefficient            (temperature)*chemicalDensity%abundance(self%atomicHydrogenCationIndex) &
            &          +hydrogenNetworkHminus_Electron_to_H_2Electron_RateCoefficient(temperature)*chemicalDensity%abundance(self%electronIndex            )
       if (destructionTerm /= 0.0d0) then
          self%densityAtomicHydrogenAnion=creationTerm/destructionTerm
       else
          if (creationTerm > 0.0d0) call Error_Report('hydrogen anion equilibrium density is infinite'//{introspection:location})
          self%densityAtomicHydrogenAnion=0.0d0
       end if
    else if (self%atomicHydrogenAnionIndex > 0) then
       ! For the slow network, if we have the hydrogen anion then use its density directly.
       self%densityAtomicHydrogenAnion=chemicalDensity%abundance(self%atomicHydrogenAnionIndex)
    else
       ! Otherwise, we have no way to compute the hydrogen anion density, so set it to zero.
       self%densityAtomicHydrogenAnion=0.0d0
    end if
    ! Compute rates. References after each call refer to the rate coefficient in Tegmark et al. (1997) and the equation number in
    ! Abel et al. (1997) respectively.
    call self%rateH_Electron_to_Hplus_2Electron  (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ;  1
    call self%rateHplus_Electron_to_H_Photon     (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) ! k_1;  2
    call self%rateH_Electron_to_Hminus_Photon    (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) ! k_2;  7
    call self%rateH_Hminus_to_H2_Electron        (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) ! k_3;  8
    call self%rateH_Hplus_to_H2plus_Photon       (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) ! k_5;  9
    call self%rateH2plus_H_to_H2_Hplus           (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) ! k_6; 10
    call self%rateH2_Hplus_to_H2plus_H           (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 11
    call self%rateH2_Electron_to_2H_Electron     (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 12
    call self%rateH2_H_to_3H                     (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 13
    call self%rateHminus_Electron_to_H_2Electron (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 14
    call self%rateHminus_H_to_2H_Electron        (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 15
    call self%rateHminus_Hplus_to_2H             (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 16
    call self%rateHminus_Hplus_to_H2plus_Electron(             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 17
    call self%rateH2plus_Electron_to_2H          (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 18
    call self%rateH2plus_Hminus_to_H2_H          (             temperature,radiation,chemicalDensity,factorClumping,chemicalRates     ) !    ; 19
    call self%rateH_Gamma_to_Hplus_Electron      (             temperature,radiation,chemicalDensity               ,chemicalRates,node) !    ; 20
    call self%rateHminus_Gamma_to_H_Electron     (             temperature,radiation,chemicalDensity               ,chemicalRates,node) ! k_4; 23
    call self%rateH2_Gamma_to_H2plus_Electron    (             temperature,radiation,chemicalDensity               ,chemicalRates,node) !    ; 24
    call self%rateH2plus_Gamma_to_H_Hplus        (             temperature,radiation,chemicalDensity               ,chemicalRates,node) ! k_7; 25
    call self%rateH2plus_Gamma_to_2Hplus_Electron(             temperature,radiation,chemicalDensity               ,chemicalRates,node) !    ; 26
    call self%rateH2_Gamma_to_H2star_to_2H       (lengthColumn,temperature,radiation,chemicalDensity               ,chemicalRates,node) !    ; 27
    call self%rateH2_Gamma_to_2H                 (             temperature,radiation,chemicalDensity               ,chemicalRates,node) !    ; 28
    return
  end subroutine hydrogenNetworkRates

  subroutine hydrogenNetworkRateH_Electron_to_Hplus_2Electron(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{e}^- \rightarrow \hbox{H}^+ + 2\hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                              , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                   =.false., reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenCationChemicalIndex        , atomicHydrogenChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                                     , rateCoefficient
    !$GLC attributes unused :: self, radiation

    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH_Electron_to_Hplus_2Electron_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex      =Chemicals_Index("AtomicHydrogen"      )
          atomicHydrogenCationChemicalIndex=Chemicals_Index("AtomicHydrogenCation")
          electronChemicalIndex            =Chemicals_Index("Electron"            )
          ! This reaction is active if all species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. atomicHydrogenCationChemicalIndex > 0 .and. electronChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH_Electron_to_Hplus_2Electron_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Get the rate coefficient.
       rateCoefficient=self%atomicIonizationRateCollisional_%rate(1,1,temperature)
       ! Compute the rate.
       rate           =rateCoefficient*chemicalDensity%abundance(electronChemicalIndex)*chemicalDensity%abundance(atomicHydrogenChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex) &
            & +rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex      , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex      ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex            , &
            & chemicalRates%abundance   (electronChemicalIndex            ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH_Electron_to_Hplus_2Electron

  subroutine hydrogenNetworkRateHplus_Electron_to_H_Photon(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^+ + \hbox{e}^- \rightarrow \hbox{H} + \gamma$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                              , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                   =.false., reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenCationChemicalIndex        , atomicHydrogenChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                                     , rateCoefficient
    !$GLC attributes unused :: self, radiation

    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateHplus_Electron_to_H_Photon_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex      =Chemicals_Index("AtomicHydrogen"      )
          atomicHydrogenCationChemicalIndex=Chemicals_Index("AtomicHydrogenCation")
          electronChemicalIndex            =Chemicals_Index("Electron"            )
          ! This reaction is active if all species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. atomicHydrogenCationChemicalIndex > 0 .and. electronChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateHplus_Electron_to_H_Photon_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Get the rate coefficient.
       rateCoefficient=self%atomicRecombinationRateRadiative_%rate(1,1,temperature)
       ! Compute the rate.
       rate=rateCoefficient*chemicalDensity%abundance(electronChemicalIndex)*chemicalDensity%abundance(atomicHydrogenCationChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex            , &
            & chemicalRates%abundance   (electronChemicalIndex            ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex      , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex      ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateHplus_Electron_to_H_Photon

  subroutine hydrogenNetworkRateH_Electron_to_Hminus_Photon(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{e}^- \rightarrow \hbox{H}^- + \gamma$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                             , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                  =.false., reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenAnionChemicalIndex        , atomicHydrogenChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                                    , rateCoefficient
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH_Electron_to_Hminus_Photon_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex     =Chemicals_Index("AtomicHydrogen"     )
          atomicHydrogenAnionChemicalIndex=Chemicals_Index("AtomicHydrogenAnion")
          electronChemicalIndex           =Chemicals_Index("Electron"           )
          ! This reaction is active if all species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. atomicHydrogenAnionChemicalIndex > 0 .and. electronChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH_Electron_to_Hminus_Photon_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Get the rate coefficient.
       rateCoefficient=hydrogenNetworkH_Electron_to_Hminus_Photon_RateCoefficient(temperature)
       ! Compute the rate.
       rate=rateCoefficient*chemicalDensity%abundance(electronChemicalIndex)*chemicalDensity%abundance(atomicHydrogenChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenAnionChemicalIndex, &
            & chemicalRates%abundance   (atomicHydrogenAnionChemicalIndex) &
            & +rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex     , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex     ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex           , &
            & chemicalRates%abundance   (electronChemicalIndex           ) &
            & -rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH_Electron_to_Hminus_Photon

  double precision function hydrogenNetworkH_Electron_to_Hminus_Photon_RateCoefficient(temperature)
    !!{
    Computes the rate coefficient (in units of cm$^3$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{e}^- \rightarrow \hbox{H}^- + \gamma$.
    !!}
    implicit none
    double precision, intent(in   ) :: temperature
    double precision, save          :: rateCoefficientStored, temperaturePrevious
    !$omp threadprivate(temperaturePrevious,rateCoefficientStored)
    double precision                :: log10Temperature

    ! Determine if we need to recompute the rate coefficient.
    if (temperature /= temperaturePrevious) then
       ! Store the new temperature.
       temperaturePrevious=temperature
       ! Compute base 10 logarithm of temperature.
       log10Temperature=log10(temperature)
       ! Compute rate coefficient.
       if      (temperature <=    1.0d0) then
          rateCoefficientStored=1.429d-18
       else if (temperature <= 6000.0d0) then
          rateCoefficientStored=1.429d-18*(temperature**0.7620d0)*(temperature**(0.1523d0*log10Temperature))*(temperature**(-3.274d-2*(log10Temperature**2)))
       else
          rateCoefficientStored=3.802d-17*(temperature**(0.1998d0*log10Temperature))*(10.0d0**((4.0415d-5*(log10Temperature**2)-5.447d-3)*(log10Temperature**4)))
       end if
    end if
    ! Return the store rate coefficient.
    hydrogenNetworkH_Electron_to_Hminus_Photon_RateCoefficient=rateCoefficientStored
    return
  end function hydrogenNetworkH_Electron_to_Hminus_Photon_RateCoefficient

  subroutine hydrogenNetworkRateH_Hminus_to_H2_Electron(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{H}^- \rightarrow \hbox{H}_2 + \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                             , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                  =.false., reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenAnionChemicalIndex        , atomicHydrogenChemicalIndex        , &
         &                                                                  chemicalHydrogenChemicalIndex           , electronChemicalIndex
    integer                                                              :: status
    double precision                                                     :: rate                                    , rateCoefficient
    !$GLC attributes unused :: radiation

    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH_Hminus_to_H2_Electron_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex     =Chemicals_Index("AtomicHydrogen"     ,status)
          atomicHydrogenAnionChemicalIndex=Chemicals_Index("AtomicHydrogenAnion",status)
          chemicalHydrogenChemicalIndex   =Chemicals_Index("MolecularHydrogen"  ,status)
          electronChemicalIndex           =Chemicals_Index("Electron"           ,status)
          ! This reaction is active if all species were found.
          reactionActive=       atomicHydrogenChemicalIndex      > 0                 &
               &         .and.  chemicalHydrogenChemicalIndex    > 0                 &
               &         .and.  electronChemicalIndex            > 0                 &
               &         .and. (atomicHydrogenAnionChemicalIndex > 0 .or. self%fast)
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH_Hminus_to_H2_Electron_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Get the rate coefficient.
       rateCoefficient=hydrogenNetworkH_Hminus_to_H2_Electron_RateCoefficient(temperature)
       ! Compute the rate.
       rate=rateCoefficient*chemicalDensity%abundance(atomicHydrogenChemicalIndex)*self%densityAtomicHydrogenAnion*factorClumping
       ! Record rate.
       if (.not.self%fast) &
            & call  chemicalRates%abundanceSet(atomicHydrogenAnionChemicalIndex, &
            &       chemicalRates%abundance   (atomicHydrogenAnionChemicalIndex) &
            &      -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex     , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex     ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex   , &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex   ) &
            & +rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex           , &
            & chemicalRates%abundance   (electronChemicalIndex           ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH_Hminus_to_H2_Electron

  double precision function hydrogenNetworkH_Hminus_to_H2_Electron_RateCoefficient(temperature)
    !!{
    Computes the rate coefficient (in units of c$^3$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{H}^- \rightarrow \hbox{H}_2 + \hbox{e}^-$.
    !!}
    use :: Numerical_Constants_Physical, only : boltzmannsConstant
    use :: Numerical_Constants_Units   , only : electronVolt
    implicit none
    double precision, intent(in   ) :: temperature
    double precision                :: logNaturalTemperatureElectronVolts, temperatureElectronVolts

    ! Compute the temperature in electron volts.
    temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
    ! Compute the rate coefficient.
    if (temperatureElectronVolts >= 0.1d0) then
       ! Get the natural logarithm of the temperature in electron volts.
       logNaturalTemperatureElectronVolts=log(temperatureElectronVolts)
       hydrogenNetworkH_Hminus_to_H2_Electron_RateCoefficient=exp(  &
            &                                      -20.069138970d0  &
            & +logNaturalTemperatureElectronVolts*(+ 0.228980000d0  &
            & +logNaturalTemperatureElectronVolts*(+ 3.599837700d-2 &
            & +logNaturalTemperatureElectronVolts*(- 4.555120000d-3 &
            & +logNaturalTemperatureElectronVolts*(- 3.105115440d-4 &
            & +logNaturalTemperatureElectronVolts*(+ 1.073294000d-4 &
            & +logNaturalTemperatureElectronVolts*(- 8.366719600d-6 &
            & +logNaturalTemperatureElectronVolts*(+ 2.238306230d-7 &
            &                                     )))))))           &
            &                                      )
    else
       hydrogenNetworkH_Hminus_to_H2_Electron_RateCoefficient=1.428d-9
    end if
    return
  end function hydrogenNetworkH_Hminus_to_H2_Electron_RateCoefficient

  subroutine hydrogenNetworkRateH_Hplus_to_H2plus_Photon(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \gamma$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : boltzmannsConstant
    use :: Numerical_Constants_Units    , only : electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                                , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                     =.false., reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenCationChemicalIndex          , atomicHydrogenChemicalIndex        , &
         &                                                                  chemicalHydrogenCationChemicalIndex
    double precision                                                     :: rate                                       , rateCoefficient                    , &
         &                                                                  temperatureElectronVolts
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH_Hplus_to_H2plus_Photon_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex        =Chemicals_Index("AtomicHydrogen"         )
          atomicHydrogenCationChemicalIndex  =Chemicals_Index("AtomicHydrogenCation"   )
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          ! This reaction is active if all species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. atomicHydrogenCationChemicalIndex > 0 .and. chemicalHydrogenCationChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH_Hplus_to_H2plus_Photon_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
       ! Compute the rate coefficient.
       if (temperatureElectronVolts < 0.577d0) then
          rateCoefficient=3.833d-16*(temperatureElectronVolts**1.8d0)
       else
          rateCoefficient=5.810d-16*((0.20651d0*temperatureElectronVolts)**(-0.2891d0*log(0.20651d0*temperatureElectronVolts)))
       end if
       ! Compute the rate.
       rate=rateCoefficient*chemicalDensity%abundance(atomicHydrogenChemicalIndex)*chemicalDensity%abundance(atomicHydrogenCationChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex        , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex        ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex  ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH_Hplus_to_H2plus_Photon

  subroutine hydrogenNetworkRateH2plus_H_to_H2_Hplus(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                                , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                                              :: reactionActive                     =.false., reactionInitialized          =.false.
    integer                                              , save          :: atomicHydrogenCationChemicalIndex          , atomicHydrogenChemicalIndex          , &
         &                                                                  chemicalHydrogenCationChemicalIndex        , chemicalHydrogenChemicalIndex
    double precision                                     , parameter     :: rateCoefficient                    =6.4d-10
    double precision                                                     :: rate
    !$GLC attributes unused :: radiation, temperature

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2plus_H_to_H2_Hplus_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex        =Chemicals_Index("AtomicHydrogen"         )
          atomicHydrogenCationChemicalIndex  =Chemicals_Index("AtomicHydrogenCation"   )
          chemicalHydrogenChemicalIndex      =Chemicals_Index("MolecularHydrogen"      )
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          ! This reaction is active if both species were found.
          reactionActive=atomicHydrogenChemicalIndex   > 0 .and. atomicHydrogenCationChemicalIndex   > 0 .and. &
               &         chemicalHydrogenChemicalIndex > 0 .and. chemicalHydrogenCationChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2plus_H_to_H2_Hplus_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute the rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenCationChemicalIndex)*chemicalDensity%abundance(atomicHydrogenChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex        , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex        ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex      , &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex      ) &
            & +rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex  ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH2plus_H_to_H2_Hplus

  subroutine hydrogenNetworkRateH2_Hplus_to_H2plus_H(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \hbox{H}$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : boltzmannsConstant
    use :: Numerical_Constants_Units    , only : electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                                , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                     =.false., reactionInitialized          =.false.
    integer                                              , save          :: atomicHydrogenCationChemicalIndex          , atomicHydrogenChemicalIndex          , &
         &                                                                  chemicalHydrogenCationChemicalIndex        , chemicalHydrogenChemicalIndex
    double precision                                                     :: logNaturalTemperatureElectronVolts         , rate                                 , &
         &                                                                  rateCoefficient                            , temperatureElectronVolts
    !$GLC attributes unused :: self, radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2_Hplus_to_H2plus_H_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex        =Chemicals_Index("AtomicHydrogen"         )
          atomicHydrogenCationChemicalIndex  =Chemicals_Index("AtomicHydrogenCation"   )
          chemicalHydrogenChemicalIndex      =Chemicals_Index("MolecularHydrogen"      )
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          ! This reaction is active if all species were found.
          reactionActive=atomicHydrogenChemicalIndex   > 0 .and. atomicHydrogenCationChemicalIndex   > 0 .and. &
               &         chemicalHydrogenChemicalIndex > 0 .and. chemicalHydrogenCationChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2_Hplus_to_H2plus_H_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
       logNaturalTemperatureElectronVolts=log(temperatureElectronVolts)
       ! Compute rate coefficient.
       if (temperature < 1000.0d0) then
          rateCoefficient=0.0d0
       else
          rateCoefficient=exp(                                           &
               & -24.24914687d0                                          &
               & +3.400824440d0 * logNaturalTemperatureElectronVolts     &
               & -3.898003960d0 *(logNaturalTemperatureElectronVolts**2) &
               & +2.045558782d0 *(logNaturalTemperatureElectronVolts**3) &
               & -0.541618285d0 *(logNaturalTemperatureElectronVolts**4) &
               & +8.410775030d-2*(logNaturalTemperatureElectronVolts**5) &
               & -7.879026150d-3*(logNaturalTemperatureElectronVolts**6) &
               & +4.138398420d-4*(logNaturalTemperatureElectronVolts**7) &
               & -9.363458880d-6*(logNaturalTemperatureElectronVolts**8) &
               &                     )
       end if
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenChemicalIndex)*chemicalDensity%abundance(atomicHydrogenCationChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex      , &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex      ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex  ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & +rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex        , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex        ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH2_Hplus_to_H2plus_H

  subroutine hydrogenNetworkRateH2_Electron_to_2H_Electron(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{e}^- \rightarrow 2\hbox{H} + \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                        , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive             =.false., reactionInitialized          =.false.
    integer                                              , save          :: atomicHydrogenChemicalIndex        , chemicalHydrogenChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                               , rateCoefficient
    !$GLC attributes unused :: self, radiation

    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2_Electron_to_2H_Electron_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex  =Chemicals_Index("AtomicHydrogen"   )
          chemicalHydrogenChemicalIndex=Chemicals_Index("MolecularHydrogen")
          electronChemicalIndex        =Chemicals_Index("Electron"         )
          ! This reaction is active if both species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. chemicalHydrogenChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2_Electron_to_2H_Electron_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       if (temperature < 1000.0d0) then
          rateCoefficient=0.0d0
       else
          rateCoefficient=5.6d-11*sqrt(temperature)*exp(-102124.0d0/temperature)
       end if
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenChemicalIndex)*chemicalDensity%abundance(electronChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex  ) &
            & +rate*2.0d0               )
    end if
    return
  end subroutine hydrogenNetworkRateH2_Electron_to_2H_Electron

  subroutine hydrogenNetworkRateH2_H_to_3H(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \hbox{H} \rightarrow 3\hbox{H}$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : boltzmannsConstant
    use :: Numerical_Constants_Units    , only : electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                        , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive             =.false., reactionInitialized          =.false.
    integer                                              , save          :: atomicHydrogenChemicalIndex        , chemicalHydrogenChemicalIndex
    double precision                                                     :: log10Temperature                   , rate                                 , &
         &                                                                  rateCoefficient                    , temperatureElectronVolts
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2_H_to_3H_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex  =Chemicals_Index("AtomicHydrogen"   )
          chemicalHydrogenChemicalIndex=Chemicals_Index("MolecularHydrogen")
          ! This reaction is active if both species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. chemicalHydrogenChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2_H_to_3H_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
       ! Compute base 10 logarithm of temperature.
       log10Temperature=log10(temperature)
       ! Compute the rate coefficient.
       if (log10Temperature < 3.0d0 .or. log10Temperature > 5.4d0) then
          rateCoefficient=0.0d0
       else
          rateCoefficient=1.067d-10*(temperatureElectronVolts**2.012d0)*exp(-(4.463d0/temperatureElectronVolts)*((1.0d0+0.2472d0*temperatureElectronVolts)**3.512d0))
       end if
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenChemicalIndex)*chemicalDensity%abundance(atomicHydrogenChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex  ) &
            & +rate*2.0d0               )
    end if
    return
  end subroutine hydrogenNetworkRateH2_H_to_3H

  subroutine hydrogenNetworkRateHminus_Electron_to_H_2Electron(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^- + \hbox{e}^- \rightarrow \hbox{H} + 2 \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                             , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                  =.false., reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenAnionChemicalIndex        , atomicHydrogenChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                                    , rateCoefficient
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH_Electron_to_Hminus_Photon_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex     =Chemicals_Index("AtomicHydrogen"     )
          atomicHydrogenAnionChemicalIndex=Chemicals_Index("AtomicHydrogenAnion")
          electronChemicalIndex           =Chemicals_Index("Electron"           )
          ! This reaction is active if both species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. atomicHydrogenAnionChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH_Electron_to_Hminus_Photon_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Get rate coefficient.
       rateCoefficient=hydrogenNetworkHminus_Electron_to_H_2Electron_RateCoefficient(temperature)
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(atomicHydrogenAnionChemicalIndex)*chemicalDensity%abundance(electronChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenAnionChemicalIndex, &
            & chemicalRates%abundance   (atomicHydrogenAnionChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex     , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex     ) &
            & +rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex           , &
            & chemicalRates%abundance   (electronChemicalIndex           ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateHminus_Electron_to_H_2Electron

  double precision function hydrogenNetworkHminus_Electron_to_H_2Electron_RateCoefficient(temperature)
    !!{
    Computes the rate coefficient (in units of cm$^3$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    !!}
    use :: Numerical_Constants_Physical, only : boltzmannsConstant
    use :: Numerical_Constants_Units   , only : electronVolt
    implicit none
    double precision, intent(in   ) :: temperature
    double precision                :: logNaturalTemperatureElectronVolts, temperatureElectronVolts

    ! Compute the temperature in electron volts.
    temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
    logNaturalTemperatureElectronVolts=log(temperatureElectronVolts)
    ! Compute rate coefficient.
    hydrogenNetworkHminus_Electron_to_H_2Electron_RateCoefficient=exp( &
         &                                      -18.01849334d0  &
         & +logNaturalTemperatureElectronVolts*(+ 2.36085220d0  &
         & +logNaturalTemperatureElectronVolts*(- 0.28274430d0  &
         & +logNaturalTemperatureElectronVolts*(+ 1.62331664d-2 &
         & +logNaturalTemperatureElectronVolts*(- 3.36501203d-2 &
         & +logNaturalTemperatureElectronVolts*(+ 1.17832978d-2 &
         & +logNaturalTemperatureElectronVolts*(- 1.65619470d-3 &
         & +logNaturalTemperatureElectronVolts*(+ 1.06827520d-4 &
         & +logNaturalTemperatureElectronVolts*(- 2.63128581d-6 &
         &                                     ))))))))         &
         &                                              )
    return
  end function hydrogenNetworkHminus_Electron_to_H_2Electron_RateCoefficient

  subroutine hydrogenNetworkRateHminus_H_to_2H_Electron(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \hbox{H} \rightarrow 2 \hbox{H} + \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : boltzmannsConstant
    use :: Numerical_Constants_Units    , only : electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                               , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                    =.false., reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenAnionChemicalIndex          , atomicHydrogenChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: logNaturalTemperatureElectronVolts        , rate                               , &
         &                                                                  rateCoefficient                           , temperatureElectronVolts
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH_Electron_to_Hminus_Photon_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex     =Chemicals_Index("AtomicHydrogen"     )
          atomicHydrogenAnionChemicalIndex=Chemicals_Index("AtomicHydrogenAnion")
          electronChemicalIndex           =Chemicals_Index("Electron"           )
          ! This reaction is active if both species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. atomicHydrogenAnionChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH_Electron_to_Hminus_Photon_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
       logNaturalTemperatureElectronVolts=log(temperatureElectronVolts)
       ! Compute rate coefficient.
       if (temperatureElectronVolts >= 0.1d0) then
          rateCoefficient=exp(                                           &
               & -20.37260896d0                                          &
               & + 1.13944330d0 * logNaturalTemperatureElectronVolts     &
               & - 0.14210136d0 *(logNaturalTemperatureElectronVolts**2) &
               & + 8.46445540d-3*(logNaturalTemperatureElectronVolts**3) &
               & - 1.43276410d-3*(logNaturalTemperatureElectronVolts**4) &
               & + 2.01225030d-4*(logNaturalTemperatureElectronVolts**5) &
               & + 8.66396320d-5*(logNaturalTemperatureElectronVolts**6) &
               & - 2.58500970d-5*(logNaturalTemperatureElectronVolts**7) &
               & + 2.45550120d-6*(logNaturalTemperatureElectronVolts**8) &
               & - 8.06838250d-8*(logNaturalTemperatureElectronVolts**9) &
               &                     )
       else
          rateCoefficient=2.5634d-9*(temperatureElectronVolts**1.78186d0)
       end if
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(atomicHydrogenAnionChemicalIndex)*chemicalDensity%abundance(atomicHydrogenChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenAnionChemicalIndex, &
            & chemicalRates%abundance   (atomicHydrogenAnionChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex     , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex     ) &
            & +rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex           , &
            & chemicalRates%abundance   (electronChemicalIndex           ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateHminus_H_to_2H_Electron

  subroutine hydrogenNetworkRateHminus_Hplus_to_2H(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                             , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                  =.false., reactionInitialized              =.false.
    integer                                              , save          :: atomicHydrogenAnionChemicalIndex        , atomicHydrogenCationChemicalIndex        , &
         &                                                                  atomicHydrogenChemicalIndex
    double precision                                                     :: rate                                    , rateCoefficient
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateHminus_Hplus_to_2H_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenCationChemicalIndex=Chemicals_Index("AtomicHydrogenCation")
          atomicHydrogenAnionChemicalIndex =Chemicals_Index("AtomicHydrogenAnion" )
          atomicHydrogenChemicalIndex      =Chemicals_Index("AtomicHydrogen"      )
          ! This reaction is active if all species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. atomicHydrogenAnionChemicalIndex > 0 .and. atomicHydrogenCationChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateHminus_Hplus_to_2H_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Get the rate coefficient.
       rateCoefficient=hydrogenNetworkHminus_Hplus_to_2H_RateCoefficient(temperature)
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(atomicHydrogenCationChemicalIndex)*chemicalDensity%abundance(atomicHydrogenAnionChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenAnionChemicalIndex , &
            & chemicalRates%abundance   (atomicHydrogenAnionChemicalIndex ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex      , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex      ) &
            & +rate*2.0d0               )
    end if
    return
  end subroutine hydrogenNetworkRateHminus_Hplus_to_2H

  double precision function hydrogenNetworkHminus_Hplus_to_2H_RateCoefficient(temperature)
    !!{
    Compute the rate coefficient (in units of c$^3$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H} \rightarrow \hbox{H}_2 + \hbox{H}^+$.
    !!}
    implicit none
    double precision, intent(in   ) :: temperature

    hydrogenNetworkHminus_Hplus_to_2H_RateCoefficient=7.0d-8/sqrt(temperature/100.0d0)
    return
  end function hydrogenNetworkHminus_Hplus_to_2H_RateCoefficient

  subroutine hydrogenNetworkRateHminus_Hplus_to_H2plus_Electron(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \hbox{H}^+ \rightarrow \hbox{H}_2^+ + \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : boltzmannsConstant
    use :: Numerical_Constants_Units    , only : electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                                , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                     =.false., reactionInitialized              =.false.
    integer                                              , save          :: atomicHydrogenAnionChemicalIndex           , atomicHydrogenCationChemicalIndex        , &
         &                                                                  chemicalHydrogenCationChemicalIndex        , electronChemicalIndex
    double precision                                                     :: rate                                       , rateCoefficient                          , &
         &                                                                  temperatureElectronVolts
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateHminus_Hplus_to_H2plus_Electron_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenCationChemicalIndex  =Chemicals_Index("AtomicHydrogenCation"   )
          atomicHydrogenAnionChemicalIndex   =Chemicals_Index("AtomicHydrogenAnion"    )
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          electronChemicalIndex              =Chemicals_Index("Electron"               )
          ! This reaction is active if both species were found.
          reactionActive= atomicHydrogenCationChemicalIndex   > 0 &
               &         .and.                                    &
               &          atomicHydrogenAnionChemicalIndex    > 0 &
               &         .and.                                    &
               &          chemicalHydrogenCationChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateHminus_Hplus_to_H2plus_Electron_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute the temperature in electron volts.
       temperatureElectronVolts=boltzmannsConstant*temperature/electronVolt
       ! Compute rate coefficient.
       if (temperatureElectronVolts < 1.719d0) then
          rateCoefficient=2.2910d-10/(temperatureElectronVolts**0.4d0)
       else
          rateCoefficient=8.4258d-10/(temperatureElectronVolts**1.4d0)*exp(-1.301d0/temperatureElectronVolts)
       end if
       ! Compute rate.
       rate   =+rateCoefficient                                              &
            &  *chemicalDensity%abundance(atomicHydrogenAnionChemicalIndex ) &
            &  *chemicalDensity%abundance(atomicHydrogenCationChemicalIndex) &
            &  *factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenAnionChemicalIndex   , &
            & chemicalRates%abundance   (atomicHydrogenAnionChemicalIndex   ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex  ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & +rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex              , &
            & chemicalRates%abundance   (electronChemicalIndex              ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateHminus_Hplus_to_H2plus_Electron

  subroutine hydrogenNetworkRateH2plus_Electron_to_2H(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{e}^- \rightarrow 2\hbox{H}$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                        , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive             =.false., reactionInitialized                =.false.
    integer                                              , save          :: atomicHydrogenChemicalIndex        , chemicalHydrogenCationChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                               , rateCoefficient
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2plus_Electron_to_2H_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex        =Chemicals_Index("AtomicHydrogen"         )
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          electronChemicalIndex              =Chemicals_Index("Electron"               )
          ! This reaction is active if both species were found.
          reactionActive=atomicHydrogenChemicalIndex > 0 .and. chemicalHydrogenCationChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2plus_Electron_to_2H_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       if (temperature < 617.0d0) then
          rateCoefficient=1.0d-8
       else
          rateCoefficient=1.32d-6/(temperature**0.76d0)
       end if
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenCationChemicalIndex)*chemicalDensity%abundance(electronChemicalIndex)*factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex        , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex        ) &
            & +rate*2.0d0               )
       call   chemicalRates%abundanceSet(electronChemicalIndex              , &
            & chemicalRates%abundance   (electronChemicalIndex              ) &
            & -rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH2plus_Electron_to_2H

  subroutine hydrogenNetworkRateH2plus_Hminus_to_H2_H(self,temperature,radiation,chemicalDensity,factorClumping,chemicalRates)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \hbox{H}^- \rightarrow \hbox{H}_2 + \hbox{H}$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature                                , factorClumping
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    logical                                              , save          :: reactionActive                     =.false., reactionInitialized          =.false.
    integer                                              , save          :: atomicHydrogenAnionChemicalIndex           , atomicHydrogenChemicalIndex          , &
         &                                                                  chemicalHydrogenCationChemicalIndex        , chemicalHydrogenChemicalIndex
    double precision                                                     :: rate                                       , rateCoefficient
    !$GLC attributes unused :: radiation

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2plus_Hminus_to_H2_H_Init)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          atomicHydrogenAnionChemicalIndex   =Chemicals_Index("AtomicHydrogenAnion"    )
          chemicalHydrogenChemicalIndex      =Chemicals_Index("MolecularHydrogen"      )
          atomicHydrogenChemicalIndex        =Chemicals_Index("AtomicHydrogen"         )
          ! This reaction is active if both species were found.
          reactionActive=chemicalHydrogenCationChemicalIndex > 0 .and. atomicHydrogenAnionChemicalIndex > 0 .and. &
               &         chemicalHydrogenChemicalIndex       > 0 .and. atomicHydrogenChemicalIndex      > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2plus_Hminus_to_H2_H_Init)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       rateCoefficient=5.0d-7*sqrt(100.0d0/temperature)
       ! Compute rate.
       rate   =+rateCoefficient                                                &
            &  *chemicalDensity%abundance(chemicalHydrogenCationChemicalIndex) &
            &  *chemicalDensity%abundance(atomicHydrogenAnionChemicalIndex   ) &
            &  *factorClumping
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenAnionChemicalIndex   , &
            & chemicalRates%abundance   (atomicHydrogenAnionChemicalIndex   ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex      , &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex      ) &
            & +rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex        , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex        ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH2plus_Hminus_to_H2_H

  subroutine hydrogenNetworkRateHminus_Gamma_to_H_Electron(self,temperature,radiation,chemicalDensity,chemicalRates,node)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}^- + \gamma \rightarrow \hbox{H} + \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : plancksConstant    , speedLight
    use :: Numerical_Constants_Units    , only : metersToAngstroms  , electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass, radiationFieldCosmicMicrowaveBackground
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    type            (treeNode                           ), intent(inout) :: node
    ! Energy range for the cross-section.
    double precision                                     , parameter     :: crossSectionEnergyLow           =0.755d0
    ! Wavelength range for the cross-section.
    double precision                                     , parameter     :: crossSectionWavelengthHigh      =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEnergyLow
    logical                                              , save          :: reactionActive                  =.false.                                                                        , reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenAnionChemicalIndex                                                                                , atomicHydrogenChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                                                                                                            , rateCoefficient
    integer                                                              :: status
    !$GLC attributes unused :: self, temperature

    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateHminus_Gamma_to_H_Electron)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenAnionChemicalIndex=Chemicals_Index("AtomicHydrogenAnion",status)
          atomicHydrogenChemicalIndex     =Chemicals_Index("AtomicHydrogen"            )
          electronChemicalIndex           =Chemicals_Index("Electron"                  )
          ! This reaction is active if all species were found.
          reactionActive=(atomicHydrogenAnionChemicalIndex > 0 .or. self%fast) .and. atomicHydrogenChemicalIndex > 0 .and. electronChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateHminus_Gamma_to_H_Electron)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       select type (radiation)
       class is (radiationFieldCosmicMicrowaveBackground)
          rateCoefficient=0.144d0*(radiation%temperature()**2.13d0)*exp(-8650.0d0/radiation%temperature())
       class default
          rateCoefficient=radiation%integrateOverCrossSection([0.0d0,crossSectionWavelengthHigh],crossSection_Hminus_Gamma_to_H_Electron_,node)
       end select
       ! Compute rate.
       rate=rateCoefficient*self%densityAtomicHydrogenAnion
       ! Record rate.
       if (.not.self%fast) call   chemicalRates%abundanceSet(atomicHydrogenAnionChemicalIndex, &
            &                     chemicalRates%abundance   (atomicHydrogenAnionChemicalIndex) &
            &                     -rate                     )
       call                       chemicalRates%abundanceSet(atomicHydrogenChemicalIndex     , &
            &                     chemicalRates%abundance   (atomicHydrogenChemicalIndex     ) &
            &                     +rate                     )
       call                       chemicalRates%abundanceSet(electronChemicalIndex           , &
            &                     chemicalRates%abundance   (electronChemicalIndex           ) &
            &                     +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateHminus_Gamma_to_H_Electron

  double precision function crossSection_Hminus_Gamma_to_H_Electron(wavelength)
    !!{
    Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}^- + \gamma \rightarrow \hbox{H} + \hbox{e}^-$
    using the fitting function given by \cite{shapiro_hydrogen_1987}, renormalized\footnote{It seems unclear what units were
    used in \protect\cite{shapiro_hydrogen_1987}, hence the recalibration.} to match the results of
    \cite{nascimento_photodetachment_1977}.
    !!}
    use :: Numerical_Constants_Physical, only : plancksConstant         , speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms       , electronVolt
    use :: Tables                      , only : table1DLogarithmicLinear
    use :: Table_Labels                , only : extrapolationTypeZero
    implicit none
    double precision                          , intent(in   ) :: wavelength
    double precision                          , parameter     :: energyThreshold=0.755d0
    double precision                          , parameter     :: energyMinimum  =energyThreshold, energyMaximum=1000.0d0
    integer                                   , parameter     :: energyCount    =100
    type            (table1DLogarithmicLinear), save          :: interpolator_
    logical                                   , save          :: initialized    =.false.
    !$omp threadprivate(interpolator_,initialized)
    double precision                                          :: energy                         , crossSection
    integer                                                   :: i

    if (.not.initialized) then
       call interpolator_%create(energyMinimum,energyMaximum,energyCount,extrapolationType=[extrapolationTypeZero,extrapolationTypeZero])
       do i=1,energyCount
          energy=interpolator_%x(i)
          ! Evaluate the fitting function for the cross-section.
          if (energy >=  energyThreshold) then
             crossSection=2.085d-16*(energy-energyThreshold)**1.5d0/energy**3
          else
             crossSection=0.000d+00
          end if
          call interpolator_%populate(crossSection,i)
       end do
       initialized=.true.
    end if
    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*metersToAngstroms/electronVolt/wavelength
    ! Evaluate the cross section.
    crossSection_Hminus_Gamma_to_H_Electron=interpolator_%interpolate(energy)
    return
  end function crossSection_Hminus_Gamma_to_H_Electron

  subroutine hydrogenNetworkRateH2plus_Gamma_to_H_Hplus(self,temperature,radiation,chemicalDensity,chemicalRates,node)
    !!{
    Computes the rate (in units of c$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow \hbox{H} + \hbox{H}^+$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : plancksConstant    , speedLight
    use :: Numerical_Constants_Units    , only : metersToAngstroms  , electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass, radiationFieldCosmicMicrowaveBackground
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    type            (treeNode                           ), intent(inout) :: node
    ! Energy range for the cross-section.
    double precision                                     , parameter     :: crossSectionEnergyLow              =2.65d0
    double precision                                     , parameter     :: crossSectionEnergyHigh             =21.00d0
    ! Wavelength range for the cross-section.
    double precision                                     , parameter     :: crossSectionWavelengthLow          =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEnergyHigh
    double precision                                     , parameter     :: crossSectionWavelengthHigh         =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEnergyLow
    logical                                              , save          :: reactionActive                     =.false.                                                                         , reactionInitialized        =.false.
    integer                                              , save          :: atomicHydrogenCationChemicalIndex                                                                                   , atomicHydrogenChemicalIndex        , &
         &                                                                  chemicalHydrogenCationChemicalIndex
    double precision                                                     :: rate                                                                                                                , rateCoefficient
    !$GLC attributes unused :: self, temperature

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2plus_Gamma_to_H_Hplus)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          atomicHydrogenChemicalIndex        =Chemicals_Index("AtomicHydrogen"         )
          atomicHydrogenCationChemicalIndex  =Chemicals_Index("AtomicHydrogenCation"   )
          ! This reaction is active if all species were found.
          reactionActive=chemicalHydrogenCationChemicalIndex > 0 .and. atomicHydrogenChemicalIndex > 0 .and. atomicHydrogenCationChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2plus_Gamma_to_H_Hplus)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       select type (radiation)
       class is (radiationFieldCosmicMicrowaveBackground)
          rateCoefficient=6.36d5*exp(-71600.0d0/radiation%temperature())
       class default
          rateCoefficient=radiation%integrateOverCrossSection([crossSectionWavelengthLow,crossSectionWavelengthHigh],crossSection_H2plus_Gamma_to_H_Hplus_,node)
       end select
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenCationChemicalIndex)
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex        , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex        ) &
            & +rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex  ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH2plus_Gamma_to_H_Hplus

  double precision function crossSection_H2plus_Gamma_to_H_Hplus(wavelength)
    !!{
    Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow \hbox{H} + \hbox{H}^+$
    as given by \cite{shapiro_hydrogen_1987}.
    !!}
    use :: Numerical_Constants_Physical, only : plancksConstant      , speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms    , electronVolt
    use :: Tables                      , only : table1DLinearLinear
    use :: Table_Labels                , only : extrapolationTypeZero
    implicit none
    double precision                     , intent(in   ) :: wavelength
    double precision                     , parameter     :: energyMinimum  =2.65d0, energyMaximum=21.00d0
    integer                              , parameter     :: energyCount    =100
    type            (table1DLinearLinear), save          :: interpolator_
    logical                              , save          :: initialized    =.false.
    !$omp threadprivate(interpolator_,initialized)
    double precision                                     :: energy                  , crossSection
    integer                                              :: i

    if (.not.initialized) then
       call interpolator_%create(energyMinimum,energyMaximum,energyCount,extrapolationType=[extrapolationTypeZero,extrapolationTypeZero])
       do i=1,energyCount
          energy=interpolator_%x(i)
          ! Evaluate the fitting function for the cross-section.
          if      (energy >=  2.65d0 .and. energy < 11.27d0) then
             crossSection=10.0d0**(-40.97d0+energy*(+6.030d+0 &
                  &                        +energy*(-0.504d+0 &
                  &                        +energy*(+1.387d-2 &
                  &                                )))        &
                  &               )
          else if (energy >= 11.27d0 .and. energy < 21.00d0) then
             crossSection=10.0d0**(-30.26d0+energy*(+2.790d+0 &
                  &                        +energy*(-0.184d+0 &
                  &                        +energy*(+3.535d-3 &
                  &                                )))        &
                  &               )
          else
             crossSection=0.0d0
          end if
          call interpolator_%populate(crossSection,i)
       end do
       initialized=.true.
    end if
    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*metersToAngstroms/electronVolt/wavelength
    ! Evaluate the cross-section.
    crossSection_H2plus_Gamma_to_H_Hplus=interpolator_%interpolate(energy)
    return
  end function crossSection_H2plus_Gamma_to_H_Hplus

  subroutine hydrogenNetworkRateH2_Gamma_to_H2star_to_2H(self,lengthColumn,temperature,radiation,chemicalDensity,chemicalRates,node)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow H_2^* \rightarrow
    2\hbox{H}$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Math     , only : Pi
    use :: Numerical_Constants_Physical , only : plancksConstant    , speedLight
    use :: Numerical_Constants_Units    , only : metersToAngstroms  , electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: lengthColumn                                                                                           , temperature
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    type            (treeNode                           ), intent(inout) :: node
    logical                                              , save          :: reactionActive             =.false.                                                                    , reactionInitialized                     =.false.
    integer                                              , save          :: atomicHydrogenChemicalIndex                                                                            , chemicalHydrogenChemicalIndex
    ! Median energy of the Lyman band in chemical hydrogen (in eV).
    double precision                                     , parameter     :: energyLymanBand            =12.87d0
    ! Corresponding median wavelength of the Lyman band in chemical hydrogen (in Å).
    double precision                                     , parameter     :: wavelengthLymanBand        =metersToAngstroms*plancksConstant*speedLight/(energyLymanBand*electronVolt)
    ! Exponent appearing in self-shielding equation.
    double precision                                     , parameter     :: alpha                      =1.1d0
    double precision                                                     :: rate                                                                                                   , rateCoefficient                                  , &
         &                                                                  columnDensityMolecularHydrogen                                                                         , columnDensityMolecularHydrogenNormalized         , &
         &                                                                  velocitySpread                                                                                         , velocitySpreadNormalized                         , &
         &                                                                  factorSelfShielding
    !$GLC attributes unused :: self, temperature

    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2_Gamma_to_H2star_to_2H)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          chemicalHydrogenChemicalIndex=Chemicals_Index("MolecularHydrogen")
          atomicHydrogenChemicalIndex  =Chemicals_Index("AtomicHydrogen"   )
          ! This reaction is active if all species were found.
          reactionActive=chemicalHydrogenChemicalIndex > 0 .and. atomicHydrogenChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2_Gamma_to_H2star_to_2H)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       rateCoefficient=1.1d8*(4.0d0*Pi*radiation%flux(wavelengthLymanBand,node))
       ! Compute the self-shielding factor.
       if (self%includeSelfShielding) then
          ! Apply the self-shielding model from Safranek-Shrader et al. (2012;
          ! https://ui.adsabs.harvard.edu/abs/2012MNRAS.426.1159S). Specifically their equation (11).
          columnDensityMolecularHydrogen          =+chemicalDensity%abundance(chemicalHydrogenChemicalIndex)          &
               &                                   *lengthColumn
          columnDensityMolecularHydrogenNormalized=+columnDensityMolecularHydrogen                                    &
               &                                   /5.0d14
          velocitySpread                          =+9.12d0                                                            &
               &                                   *sqrt(                                                             &
               &                                         +temperature                                                 &
               &                                         /1.0d4                                                       &
               &                                        )
          velocitySpreadNormalized                =+velocitySpread                                                    &
               &                                   /1.0d0
          if (columnDensityMolecularHydrogenNormalized > 0.0d0) then
             factorSelfShielding                     =+0.965d0                                                           &
                  &                                   /(                                                                 &
                  &                                     +1.0d0                                                           &
                  &                                     +columnDensityMolecularHydrogenNormalized                        &
                  &                                     /velocitySpreadNormalized                                        &
                  &                                    )**alpha                                                          &
                  &                                   +     3.5d-2/sqrt(1.0d0+columnDensityMolecularHydrogenNormalized)  &
                  &                                   *exp(-8.5d-4*sqrt(1.0d0+columnDensityMolecularHydrogenNormalized))
          else
             factorSelfShielding                  =+1.0d0
          end if
       else
          factorSelfShielding                     =+1.0d0
       end if
       ! Compute rate.
       rate   =+rateCoefficient                                          &
            &  *factorSelfShielding                                      &
            &  *chemicalDensity%abundance(chemicalHydrogenChemicalIndex)
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex) &
            & -      rate               )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex  ) &
            & +2.0d0*rate               )
    end if
    return
  end subroutine hydrogenNetworkRateH2_Gamma_to_H2star_to_2H

  subroutine hydrogenNetworkRateH2_Gamma_to_H2plus_Electron(self,temperature,radiation,chemicalDensity,chemicalRates,node)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow \hbox{H}_2^+ + \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : plancksConstant    , speedLight
    use :: Numerical_Constants_Units    , only : metersToAngstroms  , electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    type            (treeNode                           ), intent(inout) :: node
    logical                                              , save          :: reactionActive                     =.false.                                                                         , reactionInitialized          =.false.
    ! Energy of the edge in the cross-section.
    double precision                                     , parameter     :: crossSectionEdgeEnergy             =15.42d0
    ! Wavelength of the edge in the cross-section.
    double precision                                     , parameter     :: crossSectionEdgeWavelength         =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEdgeEnergy
    integer                                              , save          :: chemicalHydrogenCationChemicalIndex                                                                                 , chemicalHydrogenChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                                                                                                                , rateCoefficient
    !$GLC attributes unused :: self, temperature

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2_Gamma_to_H2plus_Electron)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          chemicalHydrogenChemicalIndex      =Chemicals_Index("MolecularHydrogen"      )
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          electronChemicalIndex              =Chemicals_Index("Electron"               )
          ! This reaction is active if all species were found.
          reactionActive=       chemicalHydrogenChemicalIndex       > 0 &
               &          .and. chemicalHydrogenCationChemicalIndex > 0 &
               &          .and. electronChemicalIndex               > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2_Gamma_to_H2plus_Electron)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       rateCoefficient=radiation%integrateOverCrossSection([0.0d0,crossSectionEdgeWavelength],crossSection_H2_Gamma_to_H2plus_Electron_,node)
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenChemicalIndex)
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex      , &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex      ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & +rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex              , &
            & chemicalRates%abundance   (electronChemicalIndex              ) &
            & +rate                     )

    end if
    return
  end subroutine hydrogenNetworkRateH2_Gamma_to_H2plus_Electron

  double precision function crossSection_H2_Gamma_to_H2plus_Electron(wavelength)
    !!{
    Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow \hbox{H}_2^+ +
    \hbox{e}^-$ as given by\footnote{\protect\cite{abel_modeling_1997} cite ``O'Neil \& Reinhardt (1978)'' as the source for
    this fit, but it is not listed in their bibliography, and I have not been able to locate by any other means.}
    \cite{abel_modeling_1997}.
    !!}
    use :: Numerical_Constants_Physical, only : plancksConstant         , speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms       , electronVolt
    use :: Tables                      , only : table1DLogarithmicLinear
    use :: Table_Labels                , only : extrapolationTypeZero
    implicit none
    double precision                          , intent(in   ) :: wavelength
    double precision                          , parameter     :: energyMinimum  =15.42d0, energyMaximum=1000.0d0
    integer                                   , parameter     :: energyCount    =100
    type            (table1DLogarithmicLinear), save          :: interpolator_
    logical                                   , save          :: initialized    =.false.
    !$omp threadprivate(interpolator_,initialized)
    double precision                                          :: energy                 , crossSection
    integer                                                   :: i

    if (.not.initialized) then
       call interpolator_%create(energyMinimum,energyMaximum,energyCount,extrapolationType=[extrapolationTypeZero,extrapolationTypeZero])
       do i=1,energyCount
          energy=interpolator_%x(i)
          ! Evaluate the fitting function for the cross-section.
          if      (energy < 15.42d0) then
             crossSection=0.0d0
          else if (energy < 16.50d0) then
             crossSection=6.2d-18*energy-9.40d-17
          else if (energy < 17.70d0) then
             crossSection=1.4d-18*energy-1.48d-17
          else
             crossSection=2.5d-14/energy**2.71d0
          end if
          call interpolator_%populate(crossSection,i)
       end do
       initialized=.true.
    end if
    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*metersToAngstroms/electronVolt/wavelength
    ! Evaluate the cross-section.
    crossSection_H2_Gamma_to_H2plus_Electron=interpolator_%interpolate(energy)
    return
  end function crossSection_H2_Gamma_to_H2plus_Electron

  subroutine hydrogenNetworkRateH2plus_Gamma_to_2Hplus_Electron(self,temperature,radiation,chemicalDensity,chemicalRates,node)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow 2\hbox{H}^+ +
    \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : plancksConstant    , speedLight
    use :: Numerical_Constants_Units    , only : metersToAngstroms  , electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    type            (treeNode                           ), intent(inout) :: node
    logical                                              , save          :: reactionActive                   =.false.                                                                         , reactionInitialized                =.false.
    ! Energy range for the cross-section.
    double precision                                     , parameter     :: crossSectionEnergyLow            =30.0d0
    double precision                                     , parameter     :: crossSectionEnergyHigh           =90.0d0
    ! Wavelength range for the cross-section.
    double precision                                     , parameter     :: crossSectionWavelengthLow        =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEnergyHigh
    double precision                                     , parameter     :: crossSectionWavelengthHigh       =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEnergyLow
    integer                                              , save          :: atomicHydrogenCationChemicalIndex                                                                                 , chemicalHydrogenCationChemicalIndex        , &
         &                                                                  electronChemicalIndex
    double precision                                                     :: rate                                                                                                              , rateCoefficient
    !$GLC attributes unused :: self, temperature

    ! If using the fast network, this reaction is ignored so simply return in such cases.
    if (self%fast) return
    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2plus_Gamma_to_2Hplus_Electron)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenCationChemicalIndex  =Chemicals_Index("AtomicHydrogenCation"   )
          chemicalHydrogenCationChemicalIndex=Chemicals_Index("MolecularHydrogenCation")
          electronChemicalIndex              =Chemicals_Index("Electron"               )
          ! This reaction is active if all species were found.
          reactionActive=       atomicHydrogenCationChemicalIndex   > 0 &
               &          .and. chemicalHydrogenCationChemicalIndex > 0 &
               &          .and. electronChemicalIndex               > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2plus_Gamma_to_2Hplus_Electron)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       rateCoefficient=radiation%integrateOverCrossSection([crossSectionWavelengthLow,crossSectionWavelengthHigh],crossSection_H2plus_Gamma_to_2Hplus_Electron_,node)
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenCationChemicalIndex)
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenCationChemicalIndex) &
            & -      rate               )
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex  ) &
            & +2.0d0*rate               )
       call   chemicalRates%abundanceSet(electronChemicalIndex              , &
            & chemicalRates%abundance   (electronChemicalIndex              ) &
            & +     rate                )
    end if
    return
  end subroutine hydrogenNetworkRateH2plus_Gamma_to_2Hplus_Electron

  double precision function crossSection_H2plus_Gamma_to_2Hplus_Electron(wavelength)
    !!{
    Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2^+ + \gamma \rightarrow 2\hbox{H}^+ +
    \hbox{e}^-$ as given by \cite{shapiro_hydrogen_1987}.
    !!}
    use :: Numerical_Constants_Physical, only : plancksConstant      , speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms    , electronVolt
    use :: Tables                      , only : table1DLinearLinear
    use :: Table_Labels                , only : extrapolationTypeZero
    implicit none
    double precision                     , intent(in   ) :: wavelength
    double precision                     , parameter     :: energyMinimum  =30.0   , energyMaximum=90.0d0
    integer                              , parameter     :: energyCount    =100
    type            (table1DLinearLinear), save          :: interpolator_
    logical                              , save          :: initialized    =.false.
    !$omp threadprivate(interpolator_,initialized)
    double precision                                     :: energy                 , crossSection
    integer                                              :: i

    if (.not.initialized) then
       call interpolator_%create(energyMinimum,energyMaximum,energyCount,extrapolationType=[extrapolationTypeZero,extrapolationTypeZero])
       do i=1,energyCount
          energy=interpolator_%x(i)
          ! Evaluate the fitting function for the cross-section.
          if (energy >= 30.0d0 .and. energy <= 90.0d0) then
             crossSection=10.0d0**(         -16.926d+0 &
                  &                +energy*(- 4.528d-2 &
                  &                +energy*(  2.238d-4 &
                  &                +energy*(  4.245d-7 &
                  &                        )))         &
                  &               )
          else
             crossSection=0.0d0
          end if
          call interpolator_%populate(crossSection,i)
       end do
       initialized=.true.
    end if
    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*metersToAngstroms/electronVolt/wavelength
    ! Evaluate the cross-section.
    crossSection_H2plus_Gamma_to_2Hplus_Electron=interpolator_%interpolate(energy)
    return
  end function crossSection_H2plus_Gamma_to_2Hplus_Electron

  subroutine hydrogenNetworkRateH2_Gamma_to_2H(self,temperature,radiation,chemicalDensity,chemicalRates,node)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow 2\hbox{H}$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : plancksConstant    , speedLight
    use :: Numerical_Constants_Units    , only : metersToAngstroms  , electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), intent(inout) :: self
    double precision                                     , intent(in   ) :: temperature
    class           (radiationFieldClass                ), intent(inout) :: radiation
    type            (chemicalAbundances                 ), intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 ), intent(inout) :: chemicalRates
    type            (treeNode                           ), intent(inout) :: node
    logical                                              , save          :: reactionActive             =.false.                                                                         , reactionInitialized          =.false.
    ! Energy range for the cross-section.
    double precision                                     , parameter     :: crossSectionEnergyLow      =14.159d0
    double precision                                     , parameter     :: crossSectionEnergyHigh     =17.700d0
    ! Wavelength range for the cross-section.
    double precision                                     , parameter     :: crossSectionWavelengthLow  =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEnergyHigh
    double precision                                     , parameter     :: crossSectionWavelengthHigh =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEnergyLow
    integer                                              , save          :: atomicHydrogenChemicalIndex                                                                                 , chemicalHydrogenChemicalIndex
    double precision                                                     :: rate                                                                                                        , rateCoefficient
    !$GLC attributes unused :: self, temperature

    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH2_Gamma_to_2H)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex  =Chemicals_Index("AtomicHydrogen"   )
          chemicalHydrogenChemicalIndex=Chemicals_Index("MolecularHydrogen")
          ! This reaction is active if all species were found.
          reactionActive=       atomicHydrogenChemicalIndex   > 0 &
               &          .and. chemicalHydrogenChemicalIndex > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH2_Gamma_to_2H)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       rateCoefficient=radiation%integrateOverCrossSection([crossSectionWavelengthLow,crossSectionWavelengthHigh],crossSection_H2_Gamma_to_2H_,node)
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(chemicalHydrogenChemicalIndex)
       ! Record rate.
       call   chemicalRates%abundanceSet(chemicalHydrogenChemicalIndex, &
            & chemicalRates%abundance   (chemicalHydrogenChemicalIndex) &
            & -      rate               )
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex  , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex  ) &
            & +2.0d0*rate               )
    end if
    return
  end subroutine hydrogenNetworkRateH2_Gamma_to_2H

  double precision function crossSection_H2_Gamma_to_2H(wavelength)
    !!{
    Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow 2\hbox{H}$ as given by
    \cite{abel_modeling_1997}.
    !!}
    use :: Numerical_Constants_Physical, only : plancksConstant      , speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms    , electronVolt
    use :: Tables                      , only : table1DLinearLinear
    use :: Table_Labels                , only : extrapolationTypeZero
    implicit none
    double precision                     , intent(in   ) :: wavelength
    double precision                     , parameter     :: ratioOrthoToPara      = 0.000d0  ! Assume all H₂ is in the para- configuration.
    double precision                     , parameter     :: energyMinimum         =14.159d0, energyMaximum          =17.700d0
    integer                              , parameter     :: energyCount           =100
    type            (table1DLinearLinear), save          :: interpolator_
    logical                              , save          :: initialized           =.false.
    !$omp threadprivate(interpolator_,initialized)
    double precision                                     :: energy                         , crossSection                    , &
         &                                                  crossSectionLymanPara          , crossSectionWernerPara          , &
         &                                                  crossSectionLymanOrtho         , crossSectionWernerOrtho
    integer                                              :: i

    if (.not.initialized) then
       call interpolator_%create(energyMinimum,energyMaximum,energyCount,extrapolationType=[extrapolationTypeZero,extrapolationTypeZero])
       do i=1,energyCount
          energy=interpolator_%x(i)
          ! Evaluate the Lyman and Werner band cross sections for para- and ortho- configurations.
          if         (energy > 14.675d0 .and. energy <= 16.820d0) then
             crossSectionLymanPara     =10.0d0**(-18.0d0+15.1289d0-1.0513900000d+0*energy                       )
          else if    (energy > 16.820d0 .and. energy <= 17.600d0) then
             crossSectionLymanPara     =10.0d0**(-18.0d0-31.4100d0+1.8042000000d-2*energy**3-4.2339d-5*energy**5)
          else
             crossSectionLymanPara     = 0.0d0
          end if
          if         (energy > 14.675d0 .and. energy <= 17.700d0) then
             crossSectionWernerPara    =10.0d0**(-18.0d0+13.5311d0-0.9182618000d0*energy                        )
          else
             crossSectionWernerPara    = 0.0d0
          end if
          if (ratioOrthoToPara > 0.0d0) then
             if      (energy > 14.159d0 .and. energy <= 15.302d0) then
                crossSectionLymanOrtho =10.0d0**(-18.0d0+12.0218406d0-0.8194290d0*energy                        )
             else if (energy > 15.302d0 .and. energy <= 17.200d0) then
                crossSectionLymanOrtho =10.0d0**(-18.0d0+16.0464400d0-1.0824380d0*energy                        )
             else
                crossSectionLymanOrtho = 0.0d0
             end if
             if      (energy > 14.159d0 .and. energy <= 17.200d0) then
                crossSectionWernerOrtho=10.0d0**(-18.0d0+12.8736700d0-0.85088597d0*energy                       )
             else
                crossSectionWernerOrtho= 0.0d0
             end if
          else
             crossSectionLymanOrtho =0.0d0
             crossSectionWernerOrtho=0.0d0
          end if
          ! Construct the combined cross-section weighted by the appropriate ortho- to para- ratio.
          crossSection=+(      1.0d0/(ratioOrthoToPara+1.0d0))*(crossSectionLymanPara +crossSectionWernerPara ) &
               &       +(1.0d0-1.0d0/(ratioOrthoToPara+1.0d0))*(crossSectionLymanOrtho+crossSectionWernerOrtho)
          call interpolator_%populate(crossSection,i)
       end do
       initialized=.true.
    end if
    ! Convert from wavelength (in Angstroms) to energy (in eV).
    energy=plancksConstant*speedLight*metersToAngstroms/electronVolt/wavelength
    ! Evaluate the cross-section.
    crossSection_H2_Gamma_to_2H=interpolator_%interpolate(energy)
    return
  end function crossSection_H2_Gamma_to_2H

  subroutine hydrogenNetworkRateH_Gamma_to_Hplus_Electron(self,temperature,radiation,chemicalDensity,chemicalRates,node)
    !!{
    Computes the rate (in units of cm$^{-3}$ s$^{-1}$) for the reaction $\hbox{H} + \gamma \rightarrow \hbox{H}^+ +
    \hbox{e}^-$.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Index
    use :: Numerical_Constants_Physical , only : plancksConstant    , speedLight
    use :: Numerical_Constants_Units    , only : metersToAngstroms  , electronVolt
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (chemicalReactionRateHydrogenNetwork), target   , intent(inout) :: self
    double precision                                                , intent(in   ) :: temperature
    class           (radiationFieldClass                )           , intent(inout) :: radiation
    type            (chemicalAbundances                 )           , intent(in   ) :: chemicalDensity
    type            (chemicalAbundances                 )           , intent(inout) :: chemicalRates
    type            (treeNode                           )           , intent(inout) :: node
    logical                                              , save                     :: reactionActive                   =.false.                                                                        , reactionInitialized        =.false.
    ! Energy range for the cross-section (in eV).
    double precision                                     , parameter                :: crossSectionEnergyLow            =13.60d0
    ! Wavelength range for the cross-section (in Angstroms).
    double precision                                     , parameter                :: crossSectionWavelengthHigh       =plancksConstant*speedLight*metersToAngstroms/electronVolt/crossSectionEnergyLow
    integer                                              , save                     :: atomicHydrogenCationChemicalIndex                                                                                , atomicHydrogenChemicalIndex        , &
         &                                                                             electronChemicalIndex
    double precision                                                                :: rate                                                                                                             , rateCoefficient
    !$GLC attributes unused :: self, temperature

    ! Check if this reaction needs initializing.
    if (.not.reactionInitialized) then
       !$omp critical(hydrogenNetworkRateH_Gamma_to_H_Electron)
       if (.not.reactionInitialized) then
          ! Find the chemicals in this reaction.
          atomicHydrogenChemicalIndex      =Chemicals_Index("AtomicHydrogen"      )
          atomicHydrogenCationChemicalIndex=Chemicals_Index("AtomicHydrogenCation")
          electronChemicalIndex            =Chemicals_Index("Electron"            )
          ! This reaction is active if all species were found.
          reactionActive=       atomicHydrogenChemicalIndex       > 0 &
               &          .and. atomicHydrogenCationChemicalIndex > 0 &
               &          .and. electronChemicalIndex             > 0
          ! Flag that the reaction is now initialized.
          reactionInitialized=.true.
       end if
       !$omp end critical(hydrogenNetworkRateH_Gamma_to_H_Electron)
    end if
    ! Do calculation if this reaction is active.
    if (reactionActive) then
       ! Compute rate coefficient.
       self_           => self
       rateCoefficient =  radiation%integrateOverCrossSection([0.0d0,crossSectionWavelengthHigh],crossSection_H_Gamma_to_Hplus_Electron_,node)
       ! Compute rate.
       rate=rateCoefficient*chemicalDensity%abundance(atomicHydrogenChemicalIndex)
       ! Record rate.
       call   chemicalRates%abundanceSet(atomicHydrogenChemicalIndex      , &
            & chemicalRates%abundance   (atomicHydrogenChemicalIndex      ) &
            & -rate                     )
       call   chemicalRates%abundanceSet(atomicHydrogenCationChemicalIndex, &
            & chemicalRates%abundance   (atomicHydrogenCationChemicalIndex) &
            & +rate                     )
       call   chemicalRates%abundanceSet(electronChemicalIndex            , &
            & chemicalRates%abundance   (electronChemicalIndex            ) &
            & +rate                     )
    end if
    return
  end subroutine hydrogenNetworkRateH_Gamma_to_Hplus_Electron

  double precision function crossSection_H_Gamma_to_Hplus_Electron(wavelength)
    !!{
    Compute the cross-section (in units of cm$^{2}$) for the reaction $\hbox{H}_2 + \gamma \rightarrow 2\hbox{H}$ as given by
    \cite{abel_modeling_1997}.
    !!}
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    use :: Tables                    , only : table1DLogarithmicLinear
    use :: Table_Labels              , only : extrapolationTypeZero
    implicit none
    double precision                          , intent(in   ) :: wavelength
    double precision                          , parameter     :: wavelengthFactor=1.0d-3
    integer                                   , parameter     :: wavelengthCount =100
    type            (table1DLogarithmicLinear), save          :: interpolator_
    logical                                   , save          :: initialized     =.false.
    !$omp threadprivate(interpolator_,initialized)
    double precision                                          :: crossSection
    integer                                                   :: i

    if (.not.initialized) then
       call interpolator_%create(wavelengthFactor*lymanSeriesLimitWavelengthHydrogen_atomic,lymanSeriesLimitWavelengthHydrogen_atomic,wavelengthCount,extrapolationType=[extrapolationTypeZero,extrapolationTypeZero])
       do i=1,wavelengthCount
          ! Use the hydrogen photoionization cross section method.
          crossSection=self_%atomicCrossSectionIonizationPhoto_%crossSection(1,1,1,interpolator_%x(i))
          call interpolator_%populate(crossSection,i)
       end do
       initialized=.true.
    end if
    ! Evaluate the cross-section.
    crossSection_H_Gamma_to_Hplus_Electron=interpolator_%interpolate(wavelength)
    return
  end function crossSection_H_Gamma_to_Hplus_Electron

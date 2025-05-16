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
  Implements a stellar population properties class based on the instantaneous recycling approximation.
  !!}

  use, intrinsic :: ISO_C_Binding                             , only : c_size_t
  use            :: Stellar_Luminosities_Structure            , only : stellarLuminosities
  use            :: Stellar_Population_Selectors              , only : stellarPopulationSelectorClass
  use            :: Stellar_Population_Broad_Band_Luminosities, only : stellarPopulationBroadBandLuminositiesClass

  !![
  <stellarPopulationProperties name="stellarPopulationPropertiesInstantaneous">
   <description>
    A stellar population properties class based on the instantaneous recycling approximation. Specifically, given a star
    formation rate $\phi$, this method assumes a rate of increase of stellar mass of $\dot{M}_\star=(1-R)\phi$, a corresponding
    rate of decrease in fuel mass. The rate of change of the metal content of stars follows from the fuel metallicity, while
    that of the fuel changes according to
    \begin{equation}
     \dot{M}_{fuel,Z} = - (1-R) Z_\mathrm{fuel} \phi + p \phi.
    \end{equation}
    In the above $R$ is the instantaneous recycled fraction and $p$ is the yield, both of which are supplied by the \gls{imf}
    subsystem. The rate of energy input from the stellar population is computed assuming that the canonical amount of energy
    from a single stellar population (as defined by the {\normalfont \ttfamily feedbackEnergyInputAtInfinityCanonical}) is
    input instantaneously.
   </description>
  </stellarPopulationProperties>
  !!]
  type, extends(stellarPopulationPropertiesClass) :: stellarPopulationPropertiesInstantaneous
     !!{
     A stellar population properties class based on the instantaneous recycling approximation.
     !!}
     private
     class           (stellarPopulationSelectorClass             ), pointer      :: stellarPopulationSelector_              => null()
     class           (stellarPopulationBroadBandLuminositiesClass), pointer      :: stellarPopulationBroadBandLuminosities_ => null()
     type            (stellarLuminosities                        ), dimension(2) :: rateLuminosityStellarPrevious
     double precision                                             , dimension(2) :: rateStarFormationPrevious                        , fuelMetallicityPrevious, &
          &                                                                         timePrevious
     integer         (c_size_t                                   ), dimension(2) :: populationIDPrevious
     integer                                                                     :: abundanceIndex
   contains
     final     ::                  instantaneousDestructor
     procedure :: rates         => instantaneousRates
     procedure :: scales        => instantaneousScales
     procedure :: historyCount  => instantaneousHistoryCount
     procedure :: historyCreate => instantaneousHistoryCreate
  end type stellarPopulationPropertiesInstantaneous

  interface stellarPopulationPropertiesInstantaneous
     !!{
     Constructors for the \refClass{stellarPopulationPropertiesInstantaneous} stellar population class.
     !!}
     module procedure instantaneousConstructorParameters
     module procedure instantaneousConstructorInternal
  end interface stellarPopulationPropertiesInstantaneous

contains

  function instantaneousConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationPropertiesInstantaneous} stellar population properties class which takes a parameter list
    as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (stellarPopulationPropertiesInstantaneous   )                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class(stellarPopulationSelectorClass             ), pointer       :: stellarPopulationSelector_
    class(stellarPopulationBroadBandLuminositiesClass), pointer       :: stellarPopulationBroadBandLuminosities_

    !![
    <objectBuilder class="stellarPopulationSelector"              name="stellarPopulationSelector_"              source="parameters"/>
    <objectBuilder class="stellarPopulationBroadBandLuminosities" name="stellarPopulationBroadBandLuminosities_" source="parameters"/>
    !!]
    self=stellarPopulationPropertiesInstantaneous(stellarPopulationSelector_,stellarPopulationBroadBandLuminosities_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarPopulationSelector_"             />
    <objectDestructor name="stellarPopulationBroadBandLuminosities_"/>
    !!]
    return
  end function instantaneousConstructorParameters

  function instantaneousConstructorInternal(stellarPopulationSelector_,stellarPopulationBroadBandLuminosities_) result(self)
    !!{
    Internal constructor for the \refClass{stellarPopulationPropertiesInstantaneous} stellar population properties class.
    !!}
    use :: Atomic_Data, only : Abundance_Pattern_Lookup
    implicit none
    type (stellarPopulationPropertiesInstantaneous   )                        :: self
    class(stellarPopulationSelectorClass             ), intent(in   ), target :: stellarPopulationSelector_
    class(stellarPopulationBroadBandLuminositiesClass), intent(in   ), target :: stellarPopulationBroadBandLuminosities_
    !![
    <constructorAssign variables="*stellarPopulationSelector_, *stellarPopulationBroadBandLuminosities_"/>
    !!]

    self%abundanceIndex           =Abundance_Pattern_Lookup(abundanceName="solar")
    self%rateStarFormationPrevious=-huge(0.0d0)
    self%fuelMetallicityPrevious  =-huge(0.0d0)
    self%timePrevious             =-huge(0.0d0)
    self%populationIDPrevious     =-1_c_size_t
    call self%rateLuminosityStellarPrevious(1)%reset()
    call self%rateLuminosityStellarPrevious(2)%reset()
    return
  end function instantaneousConstructorInternal

  subroutine instantaneousDestructor(self)
    !!{
    Destructor for the \refClass{stellarPopulationPropertiesInstantaneous} stellar population properties class.
    !!}
    implicit none
    type(stellarPopulationPropertiesInstantaneous), intent(inout) :: self

    !![
    <objectDestructor name="self%stellarPopulationSelector_"             />
    <objectDestructor name="self%stellarPopulationBroadBandLuminosities_"/>
    !!]
    return
  end subroutine instantaneousDestructor

  subroutine instantaneousRates(self,rateStarFormation,abundancesFuel,component,node,history_,rateMassStellar,rateMassFuel,rateEnergyInput&
       &,rateAbundancesFuel,rateAbundancesStellar,rateLuminosityStellar,computeRateLuminosityStellar)
    !!{
    Return stellar population property rates of change given a star formation rate and fuel abundances.
    !!}
    use :: Abundances_Structure      , only : Abundances_Get_Metallicity            , adjustElementsReset
    use :: Galactic_Structure_Options, only : componentTypeDisk                     , componentTypeSpheroid, componentTypeAll , enumerationComponentTypeType
    use :: Galacticus_Nodes          , only : nodeComponent                         , nodeComponentBasic   , nodeComponentDisk, nodeComponentSpheroid       , &
          &                                   treeNode
    use :: Stellar_Feedback          , only : feedbackEnergyInputAtInfinityCanonical
    use :: Stellar_Populations       , only : stellarPopulationClass
    implicit none
    class           (stellarPopulationPropertiesInstantaneous), intent(inout) :: self
    double precision                                          , intent(  out) :: rateEnergyInput              , rateMassFuel              , &
         &                                                                       rateMassStellar
    type            (abundances                              ), intent(inout) :: rateAbundancesFuel           , rateAbundancesStellar
    type            (stellarLuminosities                     ), intent(inout) :: rateLuminosityStellar
    double precision                                          , intent(in   ) :: rateStarFormation
    type            (abundances                              ), intent(in   ) :: abundancesFuel
    class           (nodeComponent                           ), intent(in   ) :: component
    type            (treeNode                                ), intent(inout) :: node
    type            (history                                 ), intent(inout) :: history_
    logical                                                   , intent(in   ) :: computeRateLuminosityStellar
    class           (stellarPopulationClass                  ), pointer       :: stellarPopulation_
    class           (nodeComponentBasic                      ), pointer       :: basic
    integer         (c_size_t                                )                :: populationID
    type            (enumerationComponentTypeType            )                :: componentIndex
    double precision                                                          :: fuelMetallicity               , fuelMetalsRateOfChange   , &
         &                                                                       recycledFractionInstantaneous , stellarMetalsRateOfChange, &
         &                                                                       time                          , yieldInstantaneous
    !$GLC attributes unused :: history_

    ! Get the instantaneous recycled fraction and yields for the selected stellar population.
    stellarPopulation_            => self%stellarPopulationSelector_%select                       (rateStarFormation,abundancesFuel,component)
    populationID                  =       stellarPopulation_        %uniqueID                     (                                          )
    recycledFractionInstantaneous =       stellarPopulation_        %recycledFractionInstantaneous(                                          )
    yieldInstantaneous            =       stellarPopulation_        %yieldInstantaneous           (                                          )
    ! Get the metallicity of the fuel supply.
    fuelMetallicity=Abundances_Get_Metallicity(abundancesFuel)
    ! Set the stellar and fuel mass rates of change.
    rateMassStellar          =+(1.0d0-recycledFractionInstantaneous) *rateStarFormation
    rateMassFuel             =-                                       rateMassStellar
    ! Set energy input rate to the canonical value assuming that all energy is injected instantaneously.
    rateEnergyInput          =+feedbackEnergyInputAtInfinityCanonical*rateStarFormation
    ! Set the rates of change of the stellar and fuel metallicities.
    stellarMetalsRateOfChange=+fuelMetallicity                       *rateMassStellar
    fuelMetalsRateOfChange   =-stellarMetalsRateOfChange                                &
         &                    +yieldInstantaneous                    *rateStarFormation
    call rateAbundancesStellar%metallicitySet(stellarMetalsRateOfChange,adjustElements=adjustElementsReset,abundanceIndex=self%abundanceIndex)
    call rateAbundancesFuel   %metallicitySet(fuelMetalsRateOfChange   ,adjustElements=adjustElementsReset,abundanceIndex=self%abundanceIndex)
    ! Get the current cosmological time for this node.
    basic => node %basic()
    time  =  basic%time ()
    ! Set luminosity rates of change.
    if (computeRateLuminosityStellar) then
       select type (component)
       class is (nodeComponentDisk    )
          componentIndex=componentTypeDisk
       class is (nodeComponentSpheroid)
          componentIndex=componentTypeSpheroid
       class default
          componentIndex=componentTypeAll
       end select
       select case (componentIndex%ID)
       case (componentTypeDisk%ID,componentTypeSpheroid%ID)
          if     (                                                                        &
               &   populationID      /= self%     populationIDPrevious(componentIndex%ID) &
               &  .or.                                                                    &
               &   rateStarFormation /= self%rateStarFormationPrevious(componentIndex%ID) &
               &  .or.                                                                    &
               &   time              /= self%             timePrevious(componentIndex%ID) &
               &  .or.                                                                    &
               &   fuelMetallicity   /= self%  fuelMetallicityPrevious(componentIndex%ID) &
               & ) then
             call self%rateLuminosityStellarPrevious(componentIndex%ID)%setLuminosities(rateStarFormation,stellarPopulation_,self%stellarPopulationBroadBandLuminosities_,time,abundancesFuel)
             self%populationIDPrevious     (componentIndex%ID)=populationID
             self%rateStarFormationPrevious(componentIndex%ID)=rateStarFormation
             self%timePrevious             (componentIndex%ID)=time
             self%fuelMetallicityPrevious  (componentIndex%ID)=fuelMetallicity
          end if
          rateLuminosityStellar=self%rateLuminosityStellarPrevious(componentIndex%ID)
       case default
          call rateLuminosityStellar%setLuminosities(rateStarFormation,stellarPopulation_,self%stellarPopulationBroadBandLuminosities_,time,abundancesFuel)
       end select
    end if
    return
  end subroutine instantaneousRates

  subroutine instantaneousScales(self,massStellar,abundancesStellar,history_)
    !!{
    Set the scalings for error control on the absolute values of stellar population properties. The instantaneous method
    requires none, so just return.
    !!}
    implicit none
    class           (stellarPopulationPropertiesInstantaneous), intent(inout) :: self
    double precision                                          , intent(in   ) :: massStellar
    type            (abundances                              ), intent(in   ) :: abundancesStellar
    type            (history                                 ), intent(inout) :: history_
    !$GLC attributes unused :: self, history_, massStellar, abundancesStellar

    return
  end subroutine instantaneousScales

  integer function instantaneousHistoryCount(self)
    !!{
    Returns the number of histories required by the instantaneous stellar populations properties class.
    !!}
    implicit none
    class(stellarPopulationPropertiesInstantaneous), intent(inout) :: self
    !$GLC attributes unused :: self

    instantaneousHistoryCount=0
    return
  end function instantaneousHistoryCount

  subroutine instantaneousHistoryCreate(self,node,history_)
    !!{
    Create any history required for storing stellar population properties. The instantaneous method requires none, so don't
    create one.
    !!}
    implicit none
    class(stellarPopulationPropertiesInstantaneous), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    type (history                                 ), intent(inout) :: history_
    !$GLC attributes unused :: self, node, history_

    return
  end subroutine instantaneousHistoryCreate

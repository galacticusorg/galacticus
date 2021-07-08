!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  Implements a node operator class that accumulates an estimate of the energy radiated from the hot halo due to cooling
  following the model of \cite{benson_galaxy_2010-1}.
  !!}

  use :: Hot_Halo_Temperature_Profiles, only : hotHaloTemperatureProfileClass
  use :: Radiation_Fields             , only : radiationFieldCosmicMicrowaveBackground
  use :: Cooling_Functions            , only : coolingFunctionClass
  use :: Cosmology_Functions          , only : cosmologyFunctionsClass
  use :: Chemical_States              , only : chemicalStateClass

  !![
  <nodeOperator name="nodeOperatorCoolingEnergyRadiated">
   <description>A node operator class that accumulates an estimate of the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.</description>
   <deepCopy>
    <functionClass variables="radiation"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="radiation"/>
   </stateStorable>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCoolingEnergyRadiated
     !!{
     A node operator class that accumulates an estimate of the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
     !!}
     private
     class(cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_        => null()
     class(coolingFunctionClass                   ), pointer :: coolingFunction_           => null()
     class(hotHaloTemperatureProfileClass         ), pointer :: hotHaloTemperatureProfile_ => null()
     class(chemicalStateClass                     ), pointer :: chemicalState_             => null()
     type (radiationFieldCosmicMicrowaveBackground), pointer :: radiation                  => null()
   contains
     final     ::                          coolingEnergyRadiatedDestructor
     procedure :: differentialEvolution => coolingEnergyRadiatedDifferentialEvolution
     procedure :: nodesMerge            => coolingEnergyRadiatedNodesMerge
  end type nodeOperatorCoolingEnergyRadiated
  
  interface nodeOperatorCoolingEnergyRadiated
     !!{
     Constructors for the {\normalfont \ttfamily coolingEnergyRadiated} node operator class.
     !!}
     module procedure coolingEnergyRadiatedConstructorParameters
     module procedure coolingEnergyRadiatedConstructorInternal
  end interface nodeOperatorCoolingEnergyRadiated
  
contains

  function coolingEnergyRadiatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily coolingEnergyRadiated} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorCoolingEnergyRadiated)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
        class(hotHaloTemperatureProfileClass     ), pointer       :: hotHaloTemperatureProfile_
    class(cosmologyFunctionsClass            ), pointer      :: cosmologyFunctions_
    class(coolingFunctionClass               ), pointer      :: coolingFunction_
    class(chemicalStateClass                 ), pointer      :: chemicalState_

    !![
    <objectBuilder class="hotHaloTemperatureProfile" name="hotHaloTemperatureProfile_" source="parameters"/>
    <objectBuilder class="coolingFunction"           name="coolingFunction_"           source="parameters"/>
    <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="chemicalState"             name="chemicalState_"            source="parameters"/>
    !!]
     self=nodeOperatorCoolingEnergyRadiated(cosmologyFunctions_,coolingFunction_,hotHaloTemperatureProfile_,chemicalState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloTemperatureProfile_"/>
    <objectDestructor name="coolingFunction_"          />
    <objectDestructor name="cosmologyFunctions_"       />
    <objectDestructor name="chemicalState_"            />
    !!]
    return
  end function coolingEnergyRadiatedConstructorParameters

  function coolingEnergyRadiatedConstructorInternal(cosmologyFunctions_,coolingFunction_,hotHaloTemperatureProfile_,chemicalState_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily coolingEnergyRadiated} node operator class.
    !!}
    implicit none
    type (nodeOperatorCoolingEnergyRadiated)                        :: self
     class(cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class(coolingFunctionClass               ), intent(in   ), target :: coolingFunction_
    class(hotHaloTemperatureProfileClass     ), intent(in   ), target :: hotHaloTemperatureProfile_
    class(chemicalStateClass                 ), intent(in   ), target :: chemicalState_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *coolingFunction_, *hotHaloTemperatureProfile_, *chemicalState_"/>
    !!]

    allocate(self%radiation)
    !![
    <referenceConstruct isResult="yes" owner="self" object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
    !!]
    return
  end function coolingEnergyRadiatedConstructorInternal

  subroutine coolingEnergyRadiatedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily coolingEnergyRadiated} node operator class.
    !!}
    implicit none
    type(nodeOperatorCoolingEnergyRadiated), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%coolingFunction_"          />
    <objectDestructor name="self%hotHaloTemperatureProfile_"/>
    <objectDestructor name="self%chemicalState_"            />
    <objectDestructor name="self%radiation"                 />
    !!]
    return
  end subroutine coolingEnergyRadiatedDestructor
  
  subroutine coolingEnergyRadiatedDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Accumulates an estimate of the energy radiated from the hot halo due to cooling following the model of \cite{benson_galaxy_2010-1}.
    !!}
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                  , nodeComponentHotHalo
    use :: Abundances_Structure             , only : abundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances                  , Chemicals_Property_Count
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options, only : radiusLarge,  massTypeGalactic
    use :: Numerical_Constants_Astronomical , only : gigaYear                            , massSolar               , megaParsec
    use :: Numerical_Constants_Atomic       , only : massHydrogenAtom
    use :: Numerical_Constants_Physical     , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes     , only : hecto                               , centi
    use :: Numerical_Constants_Units        , only : ergs
    use :: Numerical_Constants_Math         , only : Pi
    implicit none
    class           (nodeOperatorCoolingEnergyRadiated), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout)          :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentBasic               )               , pointer :: basic
    class           (nodeComponentHotHalo             )               , pointer :: hotHalo
    double precision                                                            :: density                , temperature          , &
         &                                                                         massToDensityConversion, numberDensityHydrogen, &
         &                                                                         numberDensityAllSpecies, coolingFunction      , &
         &                                                                         countParticles         , massNotional
    type            (abundances                       )                         :: abundances_
    type            (chemicalAbundances               )                         :: chemicalDensities_     , chemicalMasses_
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! Compute the mass in the notional hot halo.
    basic        =>  node   %basic        ()
    hotHalo      =>  node   %hotHalo      ()
    massNotional =  +hotHalo%mass         ()                                                      &
         &          +hotHalo%outflowedMass()                                                      &
         &          +Galactic_Structure_Enclosed_Mass(node,radiusLarge,massType=massTypeGalactic)
    if (massNotional <= 0.0d0) return
    ! Compute the mean density and temperature of the hot halo.
    density    =+massNotional    &
         &      *3.0d0                    &
         &      /4.0d0                    &
         &      /Pi                       &
         &      /hotHalo%outerRadius()**3
    temperature=self%hotHaloTemperatureProfile_%temperature(node,hotHalo%outerRadius())
    ! Get the abundances for this node.
    abundances_=hotHalo%abundances()
    call abundances_%massToMassFraction(hotHalo%mass())
    ! Get the chemicals for this node.
    if (Chemicals_Property_Count () > 0) then
       chemicalMasses_=hotHalo%chemicals()
       ! Scale all chemical masses by their mass in atomic mass units to get a number density.
       call chemicalMasses_%massToNumber(chemicalDensities_)
       ! Compute factor converting mass of chemicals in (M_Solar/M_Atomic) to number density in cm⁻³.
       if (hotHalo%outerRadius() > 0.0d0) then
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(hotHalo%outerRadius())
       else
          massToDensityConversion=0.0d0
       end if
       ! Convert to number density.
       chemicalDensities_=chemicalDensities_*massToDensityConversion
    end if
    ! Set epoch for radiation field.
    call self%radiation%timeSet(basic%time())
    ! Compute number density of hydrogen (in cm⁻³).
    numberDensityHydrogen  =  +density                                    &
         &                    *abundances_     %hydrogenMassFraction()    &
         &                    *massSolar                                  &
         &                    /massHydrogenAtom                           &
         &                    /hecto                                  **3 &
         &                    /megaParsec                             **3
    ! Get the number density of all species, including electrons.
    numberDensityAllSpecies=+                                                numberDensityHydrogen                                                            &
         &                  /     abundances_   %hydrogenNumberFraction(                                                                                    ) &
         &                  +self%chemicalState_%electronDensity       (     numberDensityHydrogen,temperature,abundances_                   ,self%radiation)
    ! Get the cooling function (in ergs cm³ s⁻¹).
    coolingFunction        =+self%coolingFunction_%coolingFunction     (node,numberDensityHydrogen,temperature,abundances_,chemicalDensities_,self%radiation)
    ! Compute the number of particles in the halo.
    countParticles         =+numberDensityAllSpecies &
         &                  *4.0d0                   &
         &                  *Pi                      &
         &                  /3.0d0                   &
         &                  *(                       &
         &                    +hotHalo%outerRadius() &
         &                    *megaParsec            &
         &                    /centi                 &
         &                   )**3
    ! Set the rate. Note that the "cooling function" here is λ(T,Z) = Λ(T,Z) nₕ², where Λ(T,Z) is the usual cooling function and
    ! nₕ is the number density of hydrogen atoms. We want to evaluate the integral ∫ dt Λ(T,Z) nₕ N where N is the total number of
    ! particles in the halo. This can be written as ∫ dt λ(T,Z) nₕ⁻¹ N.
    call hotHalo%energyRadiatedRate(                       &
         &                          +coolingFunction       &
         &                          /numberDensityHydrogen &
         &                          *countParticles        &
         &                         )
    return
  end subroutine coolingEnergyRadiatedDifferentialEvolution
  

  subroutine coolingEnergyRadiatedNodesMerge(self,node)
    !!{
    Zero the radiated energy of the hot halo component of nodes about to merge.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class(nodeOperatorCoolingEnergyRadiated), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    class(nodeComponentHotHalo             ), pointer       :: hotHalo

    ! We do not add the energy radiated from this node to that of its parent, as we assume that, on merging, the hot halo gas of
    ! this node is shock heated to the virial temperature of the parent, effectively negating the energy radiated.
    hotHalo => node%hotHalo()
    call hotHalo%energyRadiatedSet(       &
         &                         +0.0d0 &
         &                        )
    return
  end subroutine coolingEnergyRadiatedNodesMerge
  

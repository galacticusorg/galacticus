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
  Implementation of a time available for cooling class using the model of \cite{benson_galaxy_2010-1}.
  !!}
  
  use :: Radiation_Fields   , only : radiationFieldCosmicMicrowaveBackground
  use :: Cooling_Functions  , only : coolingFunctionClass
  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Chemical_States    , only : chemicalStateClass

  !![
  <coolingTimeAvailable name="coolingTimeAvailableBensonBower2010">
   <description>
    A time available for cooling class implementing the model of \cite{benson_galaxy_2010-1}.
   </description>
   <deepCopy>
    <functionClass variables="radiation"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="radiation"/>
   </stateStorable>
  </coolingTimeAvailable>
  !!]
  type, extends(coolingTimeAvailableClass) :: coolingTimeAvailableBensonBower2010
     !!{
     Implementation of a time available for cooling class using the model of \cite{benson_galaxy_2010-1}.
     !!}
     private
     class  (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_ => null()
     class  (coolingFunctionClass                   ), pointer :: coolingFunction_    => null()
     class  (chemicalStateClass                     ), pointer :: chemicalState_      => null()
     type   (radiationFieldCosmicMicrowaveBackground), pointer :: radiation           => null()
     integer                                                   :: energyRadiatedID

   contains
     final     ::                              bensonBower2010Destructor
     procedure :: timeAvailable             => bensonBower2010TimeAvailable
     procedure :: timeAvailableIncreaseRate => bensonBower2010TimeAvailableIncreaseRate
  end type coolingTimeAvailableBensonBower2010

  interface coolingTimeAvailableBensonBower2010
     !!{
     Constructors for the \cite{white_galaxy_1991} time available for cooling class.
     !!}
     module procedure bensonBower2010ConstructorParameters
     module procedure bensonBower2010ConstructorInternal
  end interface coolingTimeAvailableBensonBower2010

contains

  function bensonBower2010ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{benson_galaxy_2010-1} time available for cooling class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (coolingTimeAvailableBensonBower2010)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class(coolingFunctionClass               ), pointer       :: coolingFunction_
    class(chemicalStateClass                 ), pointer       :: chemicalState_

    !![
    <objectBuilder class="coolingFunction"    name="coolingFunction_"    source="parameters"/>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="chemicalState"      name="chemicalState_"      source="parameters"/>
    !!]
    self=coolingTimeAvailableBensonBower2010(cosmologyFunctions_,coolingFunction_,chemicalState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingFunction_"   />
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="chemicalState_"     />
    !!]
    return
  end function bensonBower2010ConstructorParameters

  function bensonBower2010ConstructorInternal(cosmologyFunctions_,coolingFunction_,chemicalState_) result(self)
    !!{
    Internal constructor for the \cite{benson_galaxy_2010-1} cooling rate class.
    !!}
    implicit none
    type (coolingTimeAvailableBensonBower2010)                        :: self
    class(cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class(coolingFunctionClass               ), intent(in   ), target :: coolingFunction_
    class(chemicalStateClass                 ), intent(in   ), target :: chemicalState_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *coolingFunction_, *chemicalState_"/>
    !!]

    allocate(self%radiation)
    !![
    <referenceConstruct isResult="yes" owner="self" object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
    <addMetaProperty component="hotHalo" name="energyRadiatedBensonBower2010" isEvolvable="yes" id="self%energyRadiatedID"/>
    !!]
    return
  end function bensonBower2010ConstructorInternal

  subroutine bensonBower2010Destructor(self)
    !!{
    Destructor for the simple cooling radius class.
    !!}
    implicit none
    type(coolingTimeAvailableBensonBower2010), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%coolingFunction_"   />
    <objectDestructor name="self%chemicalState_"     />
    <objectDestructor name="self%radiation"          />
    !!]
    return
  end subroutine bensonBower2010Destructor

  double precision function bensonBower2010TimeAvailable(self,node)
    !!{
    Returns the time available for cooling (in units of Gyr).
    !!}
    use :: Galacticus_Nodes                 , only : nodeComponentBasic                  , nodeComponentHotHalo
    use :: Abundances_Structure             , only : abundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances                  , Chemicals_Property_Count
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
    use :: Mass_Distributions               , only : massDistributionClass               , kinematicsDistributionClass, massDistributionZero
    use :: Coordinates                      , only : coordinateSpherical                 , assignment(=)
    use :: Galactic_Structure_Options       , only : componentTypeHotHalo                , massTypeGaseous            , massTypeGalactic
    use :: Numerical_Constants_Astronomical , only : gigaYear                            , massSolar                  , megaParsec
    use :: Numerical_Constants_Atomic       , only : massHydrogenAtom
    use :: Numerical_Constants_Physical     , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes     , only : hecto                               , centi
    use :: Numerical_Constants_Units        , only : ergs
    use :: Numerical_Constants_Math         , only : Pi
    implicit none
    class           (coolingTimeAvailableBensonBower2010), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentBasic                 ), pointer       :: basic
    class           (nodeComponentHotHalo               ), pointer       :: hotHalo
    class           (massDistributionClass              ), pointer       :: massDistribution_
    class           (kinematicsDistributionClass        ), pointer       :: kinematicsDistribution_
    type            (coordinateSpherical                )                :: coordinates
    double precision                                                     :: density                , temperature          , &
         &                                                                  massToDensityConversion, numberDensityHydrogen, &
         &                                                                  numberDensityAllSpecies, coolingFunction      , &
         &                                                                  countParticles         , massNotional
    type            (abundances                         )                :: abundances_
    type            (chemicalAbundances                 )                :: chemicalDensities_     , chemicalMasses_

    massDistribution_ =>  node             %massDistribution(massType=massTypeGalactic)
    basic             =>  node             %basic           (                         )
    hotHalo           =>  node             %hotHalo         (                         )
    massNotional      =  +hotHalo          %mass            (                         ) &
         &               +hotHalo          %outflowedMass   (                         ) &
         &               +massDistribution_%massTotal       (                         )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (massNotional <= 0.0d0) then
       bensonBower2010TimeAvailable=0.0d0
       return
    end if
    ! Get the mass distribution.
    massDistribution_       => node             %massDistribution      (componentType=componentTypeHotHalo,massType=massTypeGaseous)
    kinematicsDistribution_ => massDistribution_%kinematicsDistribution(                                                           )
    select type (massDistribution_)
    type is (massDistributionZero)
       ! No mass distribution exists for the hot halo (most likely it is in an unphysical state).
       bensonBower2010TimeAvailable=0.0d0
       !![
       <objectDestructor name="massDistribution_"      />
       <objectDestructor name="kinematicsDistribution_"/>
       !!]
       return
    end select
    ! Compute the mean density and temperature of the hot halo.
    density    =+massNotional             &
         &      *3.0d0                    &
         &      /4.0d0                    &
         &      /Pi                       &
         &      /hotHalo%outerRadius()**3
    coordinates=[hotHalo%outerRadius(),0.0d0,0.0d0]
    temperature=+kinematicsDistribution_%temperature(coordinates)
    !![
    <objectDestructor name="massDistribution_"      />
    <objectDestructor name="kinematicsDistribution_"/>
    !!]          
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
    ! Compute the time available for cooling (Benson & Bower 2010; eqn. 16).
    if (coolingFunction > 0.0d0) then
       ! Note that the "cooling function" here is λ(T,Z) = Λ(T,Z) nₕ², where Λ(T,Z) is the usual cooling function and nₕ is the
       ! number density of hydrogen atoms. We want to evaluate the integral ∫ dt Λ(T,Z) nₕ N where N is the total number of
       ! particles in the halo. This can be written as ∫ dt λ(T,Z) nₕ⁻¹ N.
       bensonBower2010TimeAvailable=+hotHalo%floatRank0MetaPropertyGet(self%energyRadiatedID) &
            &                       /coolingFunction                                          &
            &                       *numberDensityHydrogen                                    &
            &                       /countParticles
    else if (hotHalo%floatRank0MetaPropertyGet(self%energyRadiatedID) == 0.0d0) then
       bensonBower2010TimeAvailable=+0.0d0
    else
       bensonBower2010TimeAvailable=+huge(0.0d0)
    end if
    return
  end function bensonBower2010TimeAvailable

  double precision function bensonBower2010TimeAvailableIncreaseRate(self,node)
    !!{
    Compute the rate of increase of the time available for cooling. We return a rate of 1.
    !!}
    implicit none
    class(coolingTimeAvailableBensonBower2010), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    bensonBower2010TimeAvailableIncreaseRate=1.0d0
    return
  end function bensonBower2010TimeAvailableIncreaseRate

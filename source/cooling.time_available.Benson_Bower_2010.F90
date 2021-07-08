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
  Implementation of a time available for cooling class using the model of \cite{benson_galaxy_2010-1}.
  !!}
  
  use :: Hot_Halo_Temperature_Profiles, only : hotHaloTemperatureProfileClass
  use :: Radiation_Fields             , only : radiationFieldCosmicMicrowaveBackground
  use :: Cooling_Functions            , only : coolingFunctionClass
  use :: Cosmology_Functions          , only : cosmologyFunctionsClass
  use :: Chemical_States              , only : chemicalStateClass

  use Dark_Matter_Halo_Scales

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
     class(cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_        => null()
     class(coolingFunctionClass                   ), pointer :: coolingFunction_           => null()
     class(hotHaloTemperatureProfileClass         ), pointer :: hotHaloTemperatureProfile_ => null()
     class(chemicalStateClass                     ), pointer :: chemicalState_             => null()
     class(darkMatterHaloScaleClass                     ), pointer :: darkMatterHaloScale_             => null()
     type (radiationFieldCosmicMicrowaveBackground), pointer :: radiation                  => null()
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
    class(hotHaloTemperatureProfileClass     ), pointer       :: hotHaloTemperatureProfile_
    class(cosmologyFunctionsClass            ), pointer      :: cosmologyFunctions_
    class(coolingFunctionClass               ), pointer      :: coolingFunction_
    class(chemicalStateClass                 ), pointer      :: chemicalState_
    class(darkMatterHaloScaleClass                 ), pointer      :: darkMatterHaloScale_

    !![
    <objectBuilder class="hotHaloTemperatureProfile" name="hotHaloTemperatureProfile_" source="parameters"/>
    <objectBuilder class="coolingFunction"           name="coolingFunction_"           source="parameters"/>
    <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="chemicalState"             name="chemicalState_"            source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"             name="darkMatterHaloScale_"            source="parameters"/>
    !!]
    self=coolingTimeAvailableBensonBower2010(cosmologyFunctions_,coolingFunction_,hotHaloTemperatureProfile_,chemicalState_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloTemperatureProfile_"/>
    <objectDestructor name="coolingFunction_"          />
    <objectDestructor name="cosmologyFunctions_"       />
    <objectDestructor name="chemicalState_"            />
    <objectDestructor name="darkMatterHaloScale_"            />
    !!]
    return
  end function bensonBower2010ConstructorParameters

  function bensonBower2010ConstructorInternal(cosmologyFunctions_,coolingFunction_,hotHaloTemperatureProfile_,chemicalState_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \cite{benson_galaxy_2010-1} cooling rate class.
    !!}
    implicit none
    type (coolingTimeAvailableBensonBower2010)                        :: self
    class(cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class(coolingFunctionClass               ), intent(in   ), target :: coolingFunction_
    class(hotHaloTemperatureProfileClass     ), intent(in   ), target :: hotHaloTemperatureProfile_
    class(chemicalStateClass                 ), intent(in   ), target :: chemicalState_
    class(darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *coolingFunction_, *hotHaloTemperatureProfile_, *chemicalState_, *darkMatterHaloScale_"/>
    !!]

    allocate(self%radiation)
    !![
    <referenceConstruct isResult="yes" owner="self" object="radiation" constructor="radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)"/>
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
    <objectDestructor name="self%coolingFunction_"          />
    <objectDestructor name="self%hotHaloTemperatureProfile_"/>
    <objectDestructor name="self%chemicalState_"            />
    <objectDestructor name="self%darkMatterHaloScale_"            />
    <objectDestructor name="self%radiation"                 />
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
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options, only : radiusLarge,  massTypeGalactic
    use :: Numerical_Constants_Astronomical , only : gigaYear                            , massSolar               , megaParsec
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
    double precision                                                     :: density                , temperature          , &
         &                                                                  massToDensityConversion, numberDensityHydrogen, &
         &                                                                  numberDensityAllSpecies, coolingFunction      , &
         &                                                                  countParticles         , massNotional
    type            (abundances                         )                :: abundances_
    type            (chemicalAbundances                 )                :: chemicalDensities_     , chemicalMasses_

    basic        =>  node   %basic        ()
    hotHalo      =>  node   %hotHalo      ()
    massNotional =  +hotHalo%mass         ()                                                      &
         &          +hotHalo%outflowedMass()                                                      &
         &          +Galactic_Structure_Enclosed_Mass(node,radiusLarge,massType=massTypeGalactic)
    if (massNotional <= 0.0d0) then
       bensonBower2010TimeAvailable=0.0d0
       return
    end if
    ! Compute the mean density and temperature of the hot halo.
    density    =+massNotional             &
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
    ! Compute the time available for cooling (Benson & Bower 2010; eqn. 16).
    if (coolingFunction > 0.0d0) then
       ! Note that the "cooling function" here is λ(T,Z) = Λ(T,Z) nₕ², where Λ(T,Z) is the usual cooling function and nₕ is the
       ! number density of hydrogen atoms. We want to evaluate the integral ∫ dt Λ(T,Z) nₕ N where N is the total number of
       ! particles in the halo. This can be written as ∫ dt λ(T,Z) nₕ⁻¹ N.
       bensonBower2010TimeAvailable=+hotHalo%energyRadiated() &
            &                       /coolingFunction          &
            &                       *numberDensityHydrogen    &
            &                       /countParticles
    else if (hotHalo%energyRadiated() == 0.0d0) then
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

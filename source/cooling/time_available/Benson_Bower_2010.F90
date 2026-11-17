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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implementation of a time available for cooling class using the model of :cite:t:`benson_galaxy_2010-1`.
  !!}
  
  use :: Radiation_Fields   , only : radiationFieldCosmicMicrowaveBackground
  use :: Cooling_Functions  , only : coolingFunctionClass
  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Chemical_States    , only : chemicalStateClass

  !![
  <coolingTimeAvailable name="coolingTimeAvailableBensonBower2010" docformat="rst">
   <description>
   A time available for cooling class implementing the model of :cite:t:`benson_galaxy_2010-1`. The time available for cooling is the energy radiated by the hot halo, :math:`\mathcal{E}`, divided by its current rate of increase, :math:`\dot{\mathcal{E}}`. The radiated energy is accumulated by the :galacticus-class:`nodeOperatorCoolingEnergyRadiated` node operator, which must therefore be included, and which limits it such that the time available for cooling does not exceed the cooling time of the hot halo gas at the outer radius of the hot halo.
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
     !!{RST
     Implementation of a time available for cooling class using the model of :cite:t:`benson_galaxy_2010-1`.
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
     !!{RST
     Constructors for the :cite:t:`benson_galaxy_2010-1` time available for cooling class.
     !!}
     module procedure bensonBower2010ConstructorParameters
     module procedure bensonBower2010ConstructorInternal
  end interface coolingTimeAvailableBensonBower2010

contains

  function bensonBower2010ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :cite:t:`benson_galaxy_2010-1` time available for cooling class which builds the object from a parameter set.
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
    !!{RST
    Internal constructor for the :cite:t:`benson_galaxy_2010-1` time available for cooling class.
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
    !!{RST
    Destructor for the :cite:t:`benson_galaxy_2010-1` time available for cooling class.
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
    !!{RST
    Returns the time available for cooling (in units of Gyr).
    !!}
    use :: Galacticus_Nodes         , only : nodeComponentHotHalo
    use :: Cooling_Benson_Bower_2010, only : coolingBensonBower2010Properties
    implicit none
    class           (coolingTimeAvailableBensonBower2010), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentHotHalo               ), pointer       :: hotHalo
    double precision                                                     :: numberDensityHydrogen, coolingFunction, &
         &                                                                  countParticles

    if (.not.coolingBensonBower2010Properties(node,self%coolingFunction_,self%chemicalState_,self%radiation,coolingFunction,numberDensityHydrogen,countParticles)) then
       bensonBower2010TimeAvailable=0.0d0
       return
    end if
    ! Compute the time available for cooling (Benson & Bower 2010; eqn. 16).
    hotHalo => node%hotHalo()
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
    !!{RST
    Compute the rate of increase of the time available for cooling. We return a rate of 1.
    !!}
    implicit none
    class(coolingTimeAvailableBensonBower2010), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    bensonBower2010TimeAvailableIncreaseRate=1.0d0
    return
  end function bensonBower2010TimeAvailableIncreaseRate

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
  Implements a cooling function class which sums over other cooling functions.
  !!}

  !![
  <coolingFunction name="coolingFunctionSummation">
   <description>Class providing a cooling function which sums over other cooling functions.</description>
   <linkedList type="coolantList" variable="coolants" next="next" object="coolingFunction" objectType="coolingFunctionClass"/>
  </coolingFunction>
  !!]

  type, public :: coolantList
     class(coolingFunctionClass), pointer :: coolingFunction => null()
     type (coolantList         ), pointer :: next            => null()
  end type coolantList

  type, extends(coolingFunctionClass) :: coolingFunctionSummation
     !!{
     A cooling function class which sums over other cooling functions.
     !!}
     private
     type(coolantList), pointer :: coolants => null()
   contains
     final     ::                                       summationDestructor
     procedure :: coolingFunction                    => summationCoolingFunction
     procedure :: coolingFunctionFractionInBand      => summationCoolingFunctionFractionInBand
     procedure :: coolingFunctionTemperatureLogSlope => summationCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => summationCoolingFunctionDensityLogSlope
  end type coolingFunctionSummation

  interface coolingFunctionSummation
     !!{
     Constructors for the \refClass{coolingFunctionSummation} cooling function class.
     !!}
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface coolingFunctionSummation

contains

  function summationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{coolingFunctionSummation} cooling function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (coolingFunctionSummation)                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    type   (coolantList             ), pointer       :: coolant
    integer                                          :: i

    coolant => null()
    do i=1,parameters%copiesCount('coolingFunction',zeroIfNotPresent=.true.)
       if (associated(coolant)) then
          allocate(coolant%next)
          coolant => coolant%next
       else
          allocate(self%coolants)
          coolant => self%coolants
       end if
       !![
       <objectBuilder class="coolingFunction" name="coolant%coolingFunction" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="coolingFunction"/>
    !!]
    return
  end function summationConstructorParameters

  function summationConstructorInternal(coolants) result(self)
    !!{
    Internal constructor for the \refClass{coolingFunctionSummation} cooling function class.
    !!}
    implicit none
    type(coolingFunctionSummation)                        :: self
    type(coolantList             ), target, intent(in   ) :: coolants
    type(coolantList             ), pointer               :: coolant_

    self    %coolants => coolants
    coolant_          => coolants
    do while (associated(coolant_))
       !![
       <referenceCountIncrement owner="coolant_" object="coolingFunction"/>
       !!]
       coolant_ => coolant_%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !!{
    Destructor for the \refClass{coolingFunctionSummation} cooling function class.
    !!}
    implicit none
    type(coolingFunctionSummation), intent(inout) :: self
    type(coolantList             ), pointer       :: coolant_, coolantNext

    if (associated(self%coolants)) then
       coolant_ => self%coolants
       do while (associated(coolant_))
          coolantNext => coolant_%next
          !![
          <objectDestructor name="coolant_%coolingFunction"/>
          !!]
          deallocate(coolant_)
          coolant_ => coolantNext
       end do
    end if
    return
  end subroutine summationDestructor

  double precision function summationCoolingFunction(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the cooling function summed over other cooling functions.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionSummation), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances              ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances      ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass     ), intent(inout) :: radiation
    type            (coolantList             ), pointer       :: coolant

    summationCoolingFunction =  0.0d0
    coolant                  => self%coolants
    do while (associated(coolant))
       summationCoolingFunction=                                                                &
            &                   +summationCoolingFunction                                       &
            &                   +coolant%coolingFunction%coolingFunction(                       &
            &                                                            node                 , &
            &                                                            numberDensityHydrogen, &
            &                                                            temperature          , &
            &                                                            gasAbundances        , &
            &                                                            chemicalDensities    , &
            &                                                            radiation              &
            &                                                           )
       coolant => coolant%next
    end do
    return
  end function summationCoolingFunction

  double precision function summationCoolingFunctionFractionInBand(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation,energyLow,energyHigh)
    !!{
    Return the fraction of the cooling function summed over other cooling functions due to emission in the given band.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionSummation), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: numberDensityHydrogen, temperature          , &
         &                                                       energyLow            , energyHigh
    type            (abundances              ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances      ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass     ), intent(inout) :: radiation
    type            (coolantList             ), pointer       :: coolant
    double precision                                          :: coolingFunctionTotal , coolingFunctionInBand, &
         &                                                       coolingFunction      , fractionInBand

    coolingFunctionTotal  =  0.0d0
    coolingFunctionInBand =  0.0d0
    coolant               => self%coolants
    do while (associated(coolant))
       coolingFunction       =coolant%coolingFunction%coolingFunction             (                       &
            &                                                                      node                 , &
            &                                                                      numberDensityHydrogen, &
            &                                                                      temperature          , &
            &                                                                      gasAbundances        , &
            &                                                                      chemicalDensities    , &
            &                                                                      radiation              &
            &                                                                     )
       fractionInBand       =coolant%coolingFunction%coolingFunctionFractionInBand(                       &
            &                                                                      node                 , &
            &                                                                      numberDensityHydrogen, &
            &                                                                      temperature          , &
            &                                                                      gasAbundances        , &
            &                                                                      chemicalDensities    , &
            &                                                                      radiation            , &
            &                                                                      energyLow            , &
            &                                                                      energyHigh             &
            &                                                                     )
       coolingFunctionTotal =+coolingFunctionTotal                                                        &
            &                +coolingFunction
       coolingFunctionInBand=+coolingFunctionInBand                                                       &
            &                +coolingFunction                                                             &
            &                *fractionInBand
       coolant => coolant%next
    end do
    if (coolingFunctionTotal > 0.0d0) then
       summationCoolingFunctionFractionInBand=+coolingFunctionInBand &
            &                                 /coolingFunctionTotal
    else
       summationCoolingFunctionFractionInBand=+0.0d0
    end if
    return
  end function summationCoolingFunctionFractionInBand

  double precision function summationCoolingFunctionDensityLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the logarithmic gradient with respect to density of the cooling function due to Compton scattering off of \gls{cmb}
    photons.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Error                        , only : Error_Report
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionSummation), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: numberDensityHydrogen  , temperature
    type            (abundances              ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances      ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass     ), intent(inout) :: radiation
    type            (coolantList             ), pointer       :: coolant
    double precision                                          :: coolingFunction        , coolingFunctionCumulative, &
         &                                                       coolingFunctionGradient

    coolingFunctionCumulative =  0.0d0
    coolingFunctionGradient   =  0.0d0
    coolant                   => self%coolants
    do while (associated(coolant))
       coolingFunction          =+coolant%coolingFunction%coolingFunction               (                       &
            &                                                                            node                 , &
            &                                                                            numberDensityHydrogen, &
            &                                                                            temperature          , &
            &                                                                            gasAbundances        , &
            &                                                                            chemicalDensities    , &
            &                                                                            radiation              &
            &                                                                           )
       coolingFunctionCumulative=+coolingFunctionCumulative                                                     &
            &                    +coolingFunction
       coolingFunctionGradient  =+coolingFunctionGradient                                                       &
            &                    +coolingFunction                                                               &
            &                    *coolant%coolingFunction%coolingFunctionDensityLogSlope(                       &
            &                                                                            node                 , &
            &                                                                            numberDensityHydrogen, &
            &                                                                            temperature          , &
            &                                                                            gasAbundances        , &
            &                                                                            chemicalDensities    , &
            &                                                                            radiation              &
            &                                                                           )
       coolant => coolant%next
    end do
    if (coolingFunctionCumulative /= 0.0d0) then
       summationCoolingFunctionDensityLogSlope=+coolingFunctionGradient   &
            &                                  /coolingFunctionCumulative
    else
       summationCoolingFunctionDensityLogSlope=0.0d0
       if (coolingFunctionGradient /= 0.0d0) call Error_Report('cooling function is zero but has non-zero gradient with density - logarithmic slope is undefined'//{introspection:location})
    end if
    return
  end function summationCoolingFunctionDensityLogSlope

  double precision function summationCoolingFunctionTemperatureLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the logarithmic gradient with respect to temperature of the cooling function due to Compton scattering off of
    \gls{cmb} photons.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Error                        , only : Error_Report
    use :: Radiation_Fields             , only : radiationFieldClass
    implicit none
    class           (coolingFunctionSummation), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: numberDensityHydrogen                       , temperature
    type            (abundances              ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances      ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass     ), intent(inout) :: radiation
    type            (coolantList             ), pointer       :: coolant
    double precision                                          :: coolingFunction        , coolingFunctionCumulative, &
         &                                                       coolingFunctionGradient

    coolingFunctionCumulative =  0.0d0
    coolingFunctionGradient   =  0.0d0
    coolant                   => self%coolants
    do while (associated(coolant))
       coolingFunction          =+coolant%coolingFunction%coolingFunction                   (                       &
            &                                                                                node                 , &
            &                                                                                numberDensityHydrogen, &
            &                                                                                temperature          , &
            &                                                                                gasAbundances        , &
            &                                                                                chemicalDensities    , &
            &                                                                                radiation              &
            &                                                                               )
       coolingFunctionCumulative=+coolingFunctionCumulative                                                         &
            &                    +coolingFunction
       coolingFunctionGradient  =+coolingFunctionGradient                                                           &
            &                    +coolingFunction                                                                   &
            &                    /temperature                                                                       &
            &                    *coolant%coolingFunction%coolingFunctionTemperatureLogSlope(                       &
            &                                                                                node                 , &
            &                                                                                numberDensityHydrogen, &
            &                                                                                temperature          , &
            &                                                                                gasAbundances        , &
            &                                                                                chemicalDensities    , &
            &                                                                                radiation              &
            &                                                                               )
       coolant => coolant%next
    end do
    if (coolingFunctionCumulative /= 0.0d0) then
       summationCoolingFunctionTemperatureLogSlope=+coolingFunctionGradient   &
            &                                      *temperature               &
            &                                      /coolingFunctionCumulative
    else
       summationCoolingFunctionTemperatureLogSlope=0.0d0
       if (coolingFunctionGradient /= 0.0d0) call Error_Report('cooling function is zero but has non-zero gradient with temperature - logarithmic slope is undefined'//{introspection:location})
    end if
    return
  end function summationCoolingFunctionTemperatureLogSlope

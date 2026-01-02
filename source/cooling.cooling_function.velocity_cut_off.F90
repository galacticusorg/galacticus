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
  Implements a cooling function class which cuts off another cooling function below a certain velocity and before/after a certain epoch.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  ! Enumeration for whether cut off is before or after the given epoch.
  !![
  <enumeration>
   <name>cutOffWhen</name>
   <description>Specifies whether cooling is cut off before or after the given epoch.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <entry label="before"/>
   <entry label="after" />
  </enumeration>
  !!]

  !![
  <coolingFunction name="coolingFunctionVelocityCutOff">
   <description>A cooling function class which cuts off another cooling function below a certain velocity and before/after a certain epoch.</description>
  </coolingFunction>
  !!]
  type, extends(coolingFunctionClass) :: coolingFunctionVelocityCutOff
     !!{
     A cooling function class which cuts off another cooling function below a certain velocity and before/after a certain epoch.
     !!}
     private
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_  => null()
     class           (coolingFunctionClass     ), pointer :: coolingFunction_     => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_ => null()
     double precision                                     :: velocityCutOff                , timeCutOff, &
          &                                                  redshiftCutOff
     type            (enumerationCutOffWhenType)          :: whenCutOff
     logical                                              :: useFormationNode
   contains
     !![
     <methods>
      <method description="Returns true if the cooling function is cut off." method="isCutOff"/>
     </methods>
     !!]
     final     ::                                       velocityCutOffDestructor
     procedure :: coolingFunction                    => velocityCutOffCoolingFunction
     procedure :: coolingFunctionFractionInBand      => velocityCutOffCoolingFunctionFractionInBand
     procedure :: coolingFunctionTemperatureLogSlope => velocityCutOffCoolingFunctionTemperatureLogSlope
     procedure :: coolingFunctionDensityLogSlope     => velocityCutOffCoolingFunctionDensityLogSlope
     procedure :: isCutOff                           => velocityCutOffIsCutOff
  end type coolingFunctionVelocityCutOff

  interface coolingFunctionVelocityCutOff
     !!{
     Constructors for the \refClass{coolingFunctionVelocityCutOff} cooling function class.
     !!}
     module procedure velocityCutOffConstructorParameters
     module procedure velocityCutOffConstructorInternal
  end interface coolingFunctionVelocityCutOff

contains

  function velocityCutOffConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{coolingFunctionVelocityCutOff} cooling function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingFunctionVelocityCutOff)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (coolingFunctionClass         ), pointer       :: coolingFunction_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    double precision                                               :: velocityCutOff       , redshiftCutOff
    type            (varying_string               )                :: whenCutOff
    logical                                                        :: useFormationNode

    !![
    <inputParameter>
      <name>useFormationNode</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>Specifies whether to use the virial velocity of the formation node or current node in the cooling rate ``cut-off'' modifier.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityCutOff</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The velocity below which cooling is suppressed in the ``cut-off'' cooling rate modifier method.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftCutOff</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>The redshift below which cooling is suppressed in the ``cut-off'' cooling rate modifier method.</description>
    </inputParameter>
    <inputParameter>
      <name>whenCutOff</name>
      <defaultValue>var_str('after')</defaultValue>
      <source>parameters</source>
      <description>Specifies whether cooling is cut off before or after {\normalfont \ttfamily [redshiftCutOff]}.</description>
    </inputParameter>
    <objectBuilder class="coolingFunction"     name="coolingFunction_"     source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=coolingFunctionVelocityCutOff(                                                                          &
         &                                                                               velocityCutOff        , &
         &                             cosmologyFunctions_ %cosmicTime                 (                         &
         &                              cosmologyFunctions_%expansionFactorFromRedshift (                        &
         &                                                                               redshiftCutOff          &
         &                                                                              )                        &
         &                                                                             )                       , &
         &                             enumerationCutOffWhenEncode                     (                         &
         &                                                                               char(whenCutOff)      , &
         &                                                                               includesPrefix=.false.  &
         &                                                                             )                       , &
         &                                                                               useFormationNode      , &
         &                                                                               cosmologyFunctions_   , &
         &                                                                               darkMatterHaloScale_  , &
         &                                                                               coolingFunction_        &
         &                            )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="coolingFunction_"    />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function velocityCutOffConstructorParameters

  function velocityCutOffConstructorInternal(velocityCutOff,timeCutOff,whenCutOff,useFormationNode,cosmologyFunctions_,darkMatterHaloScale_,coolingFunction_) result(self)
    !!{
    Internal constructor for the \refClass{coolingFunctionVelocityCutOff} cooling function class.
    !!}
    implicit none
    type            (coolingFunctionVelocityCutOff)                        :: self
    double precision                               , intent(in   )         :: velocityCutOff      , timeCutOff
    type            (enumerationCutOffWhenType    ), intent(in   )         :: whenCutOff
    logical                                        , intent(in   )         :: useFormationNode
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    class           (coolingFunctionClass         ), intent(in   ), target :: coolingFunction_
    !![
    <constructorAssign variables="velocityCutOff, timeCutOff, whenCutOff, useFormationNode, *cosmologyFunctions_, *coolingFunction_, *darkMatterHaloScale_"/>
    !!]

    self%redshiftCutOff=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeCutOff))
    return
  end function velocityCutOffConstructorInternal

  subroutine velocityCutOffDestructor(self)
    !!{
    Destructor for the \refClass{coolingFunctionVelocityCutOff} cooling function class.
    !!}
    implicit none
    type(coolingFunctionVelocityCutOff), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%coolingFunction_"    />
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine velocityCutOffDestructor
  
  logical function velocityCutOffIsCutOff(self,node)
    !!{
    Return true if this cooling function is to be cut off.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (coolingFunctionVelocityCutOff), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    class           (nodeComponentBasic           ), pointer       :: basic
    double precision                                               :: velocityVirial

    select case (self%useFormationNode)
    case (.false.)
       velocityVirial=self%darkMatterHaloScale_%velocityVirial(node              )
    case (.true. )
       velocityVirial=self%darkMatterHaloScale_%velocityVirial(node%formationNode)
    end select
    basic                  =>  node%basic()
    velocityCutOffIsCutOff =  (                                                                               &
         &                      (basic%time()   >= self%timeCutOff .and. self%whenCutOff == cutOffWhenAfter ) &
         &                     .or.                                                                           &
         &                      (basic%time()   <= self%timeCutOff .and. self%whenCutOff == cutOffWhenBefore) &
         &                    )                                                                               &
         &                     .and.                                                                          &
         &                       velocityVirial <= self%velocityCutOff 
    return
  end function velocityCutOffIsCutOff

  double precision function velocityCutOffCoolingFunction(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the cooling function cut off above/below a given velocity.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (coolingFunctionVelocityCutOff), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                   ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances           ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass          ), intent(inout) :: radiation

    if (self%isCutOff(node)) then
       velocityCutOffCoolingFunction=0.0d0
    else
       velocityCutOffCoolingFunction=self%coolingFunction_%coolingFunction(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    end if
    return
  end function velocityCutOffCoolingFunction

  double precision function velocityCutOffCoolingFunctionFractionInBand(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation,energyLow,energyHigh)
    !!{
    Return the fraction of the cooling function due to emission in the given band.
    !!}
    implicit none
    class           (coolingFunctionVelocityCutOff), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: numberDensityHydrogen, temperature, &
         &                                                            energyLow            , energyHigh
    type            (abundances                   ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances           ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass          ), intent(inout) :: radiation

    if (self%isCutOff(node)) then
       velocityCutOffCoolingFunctionFractionInBand=0.0d0
    else
       velocityCutOffCoolingFunctionFractionInBand=self%coolingFunction_%coolingFunctionFractionInBand(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation,energyLow,energyHigh)
    end if
    return
  end function velocityCutOffCoolingFunctionFractionInBand

  double precision function velocityCutOffCoolingFunctionDensityLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the logarithmic gradient with respect to density of the cooling function.
    !!}
   implicit none
    class           (coolingFunctionVelocityCutOff), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                   ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances           ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass          ), intent(inout) :: radiation

    if (self%isCutOff(node)) then
       velocityCutOffCoolingFunctionDensityLogSlope=0.0d0
    else
       velocityCutOffCoolingFunctionDensityLogSlope=self%coolingFunction_%coolingFunctionDensityLogSlope(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    end if
    return
  end function velocityCutOffCoolingFunctionDensityLogSlope

  double precision function velocityCutOffCoolingFunctionTemperatureLogSlope(self,node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    !!{
    Return the logarithmic gradient with respect to temperature of the cooling function.
    !!}
    implicit none
    class           (coolingFunctionVelocityCutOff), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances                   ), intent(in   ) :: gasAbundances
    type            (chemicalAbundances           ), intent(in   ) :: chemicalDensities
    class           (radiationFieldClass          ), intent(inout) :: radiation

    if (self%isCutOff(node)) then
       velocityCutOffCoolingFunctionTemperatureLogSlope=0.0d0
    else
       velocityCutOffCoolingFunctionTemperatureLogSlope=self%coolingFunction_%coolingFunctionTemperatureLogSlope(node,numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
    end if
    return
  end function velocityCutOffCoolingFunctionTemperatureLogSlope

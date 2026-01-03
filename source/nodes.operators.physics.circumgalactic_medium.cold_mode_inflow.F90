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
  Implements a node operator class that implements cold mode inflow of gas from the \gls{cgm}.
  !!}

  use :: Cooling_Cold_Mode_Infall_Rates, only : coldModeInfallRateClass
  use :: Cooling_Infall_Torques        , only : coolingInfallTorqueClass

  !![
  <nodeOperator name="nodeOperatorCGMColdModeInflow">
   <description>
    A node operator class that implements cold mode inflow of gas from the \gls{cgm}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMColdModeInflow
     !!{
     A node operator class that implements cold mode inflow of gas from the \gls{cgm}.
     !!}
     private
     class(coldModeInfallRateClass     ), pointer :: coldModeInfallRate_  => null()
     class(coolingInfallTorqueClass    ), pointer :: coolingInfallTorque_ => null()
     type (enumerationComponentTypeType)          :: component
     type (enumerationCoolingFromType  )          :: coolingFrom
   contains
     final     ::                          cgmColdModeInflowDestructor
     procedure :: differentialEvolution => cgmColdModeInflowDifferentialEvolution
  end type nodeOperatorCGMColdModeInflow
  
  interface nodeOperatorCGMColdModeInflow
     !!{
     Constructors for the \refClass{nodeOperatorCGMColdModeInflow} node operator class.
     !!}
     module procedure cgmColdModeInflowConstructorParameters
     module procedure cgmColdModeInflowConstructorInternal
  end interface nodeOperatorCGMColdModeInflow
  
contains
  
  function cgmColdModeInflowConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMColdModeInflow} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    use :: Cooling_Options           , only : enumerationCoolingFromEncode
    implicit none
    type (nodeOperatorCGMColdModeInflow)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(coldModeInfallRateClass      ), pointer       :: coldModeInfallRate_
    class(coolingInfallTorqueClass     ), pointer       :: coolingInfallTorque_
    type (varying_string               )                :: component           , coolingFrom

    !![
    <inputParameter>
      <name>component</name>
      <source>parameters</source>
      <description>The component to which cooling gas should be directed.</description>
    </inputParameter>
    <inputParameter>
      <name>coolingFrom</name>
      <defaultValue>var_str('currentNode')</defaultValue>
      <description>Specifies whether the angular momentum of cooling gas should be computed from the {\normalfont \ttfamily currentNode} or the {\normalfont \ttfamily formationNode}.</description>
      <source>parameters</source>
      <variable>coolingFrom</variable>
    </inputParameter>
    <objectBuilder class="coldModeInfallRate"  name="coldModeInfallRate_"  source="parameters"/>
    <objectBuilder class="coolingInfallTorque" name="coolingInfallTorque_" source="parameters"/>
    !!]
    self=nodeOperatorCGMColdModeInflow(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),enumerationCoolingFromEncode(char(coolingFrom),includesPrefix=.false.),coldModeInfallRate_,coolingInfallTorque_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coldModeInfallRate_" />
    <objectDestructor name="coolingInfallTorque_"/>
    !!]
    return
  end function cgmColdModeInflowConstructorParameters

  function cgmColdModeInflowConstructorInternal(component,coolingFrom,coldModeInfallRate_,coolingInfallTorque_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMColdModeInflow} node operator class.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, componentTypeNone
    use :: Error                     , only : Error_Report
    implicit none
    type (nodeOperatorCGMColdModeInflow)                        :: self
    type (enumerationComponentTypeType ), intent(in   )         :: component
    type (enumerationCoolingFromType   ), intent(in   )         :: coolingFrom
    class(coldModeInfallRateClass      ), intent(in   ), target :: coldModeInfallRate_
    class(coolingInfallTorqueClass     ), intent(in   ), target :: coolingInfallTorque_
    !![
    <constructorAssign variables="component, coolingFrom, *coldModeInfallRate_, *coolingInfallTorque_"/>
    !!]

    if     (                                                                                                             &
         &   component /= componentTypeNone                                                                              &
         &  .and.                                                                                                        &
         &   component /= componentTypeDisk                                                                              &
         &  .and.                                                                                                        &
         &   component /= componentTypeSpheroid                                                                          &
         & ) call Error_Report("only 'disk', 'spheroid', and 'none' components are supported"//{introspection:location})
    return
  end function cgmColdModeInflowConstructorInternal

  subroutine cgmColdModeInflowDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorCGMColdModeInflow} node operator class.
    !!}
    implicit none
    type(nodeOperatorCGMColdModeInflow), intent(inout) :: self

    !![
    <objectDestructor name="self%coldModeInfallRate_" />
    <objectDestructor name="self%coolingInfallTorque_"/>
    !!]
    return
  end subroutine cgmColdModeInflowDestructor

  subroutine cgmColdModeInflowDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform inflow of \gls{cgm} cold mode gas.
    !!}
    use :: Abundances_Structure      , only : abundances            , operator(*)
    use :: Cooling_Options           , only : coolingFromCurrentNode, coolingFromFormationNode
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo  , nodeComponentDisk       , nodeComponentSpheroid
    use :: Galactic_Structure_Options, only : componentTypeDisk     , componentTypeSpheroid   , componentTypeNone
    use :: Error                     , only : Error_Report
    implicit none
    class           (nodeOperatorCGMColdModeInflow), intent(inout), target  :: self
    type            (treeNode                     ), intent(inout), target  :: node
    logical                                        , intent(inout)          :: interrupt
    procedure       (interruptTask                ), intent(inout), pointer :: functionInterrupt
    integer                                        , intent(in   )          :: propertyType
    type            (treeNode                     )               , pointer :: nodeCooling
    class           (nodeComponentHotHalo         )               , pointer :: hotHalo                    , hotHaloCooling
    class           (nodeComponentDisk            )               , pointer :: disk
    class           (nodeComponentSpheroid        )               , pointer :: spheroid
    type            (abundances                   )                         :: rateMassAbundancesInflow
    double precision                                                        :: rateMassInflow             , rateAngularMomentumInflow, &
         &                                                                     fractionLossAngularMomentum
    !$GLC attributes unused :: propertyType

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Compute the cold mode infall rate in this halo.
    rateMassInflow=self%coldModeInfallRate_%infallRate(node)
    ! Ignore zero rates.
    if     (                                       &
         &          rateMassInflow        == 0.0d0 &
         & .or.                                    &
         &  hotHalo%massCold           () <= 0.0d0 &
         & .or.                                    &
         &  hotHalo%angularMomentumCold() <= 0.0d0 &
         & ) return
    ! Find the node to use for cooling calculations.
    select case (self%coolingFrom%ID)
    case (coolingFromCurrentNode  %ID)
       nodeCooling => node
    case (coolingFromFormationNode%ID)
       nodeCooling => node%formationNode
    case default
       nodeCooling => null()
       call Error_Report('unknown `coolingFrom` - this should not happen'//{introspection:location})
    end select
    hotHaloCooling => nodeCooling%hotHalo()
    ! Compute the infall rate of angular momentum.
    rateAngularMomentumInflow  =rateMassInflow*hotHalo       %angularMomentumCold()/hotHalo       %massCold()
    ! Get the rate of change of abundances.
    rateMassAbundancesInflow   =rateMassInflow*hotHaloCooling%abundancesCold     ()/hotHaloCooling%massCold()
    ! Find the fraction of angular momentum lost during infall.
    fractionLossAngularMomentum=self%coolingInfallTorque_%fractionAngularMomentumLoss(node)
    ! Remove mass, angular momentum, and abundances from the hot component, and add to whichever component they should flow to.
    call    hotHalo %massColdRate           (-rateMassInflow                                       )
    call    hotHalo %angularMomentumColdRate(-rateAngularMomentumInflow                            )
    call    hotHalo %abundancesColdRate     (-rateMassAbundancesInflow                             )
    if      (self%component == componentTypeDisk    ) then
       disk     => node%disk    ()
       call disk    %massGasRate            (+rateMassInflow                                               ,interrupt,functionInterrupt)
       call disk    %angularMomentumRate    (sign(+rateAngularMomentumInflow*(1.0d0-fractionLossAngularMomentum),rateMassInflow),interrupt,functionInterrupt)
       call disk    %abundancesGasRate      (+rateMassAbundancesInflow                                     ,interrupt,functionInterrupt)
    else if (self%component == componentTypeSpheroid) then
       spheroid => node%spheroid()
       call spheroid%massGasRate            (+rateMassInflow                                               ,interrupt,functionInterrupt)
       call spheroid%angularMomentumRate    (sign(+rateAngularMomentumInflow*(1.0d0-fractionLossAngularMomentum),rateMassInflow),interrupt,functionInterrupt)
       call spheroid%abundancesGasRate      (+rateMassAbundancesInflow                                     ,interrupt,functionInterrupt)
    else if (self%component /= componentTypeNone    ) then
          call Error_Report('unexpected component - this should not happen'//{introspection:location})
    end if
    return
  end subroutine cgmColdModeInflowDifferentialEvolution
  

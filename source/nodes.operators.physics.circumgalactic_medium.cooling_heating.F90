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
  Implements a node operator class that implements inflow/outflow of gas from the \gls{cgm} due to cooling/heating.
  !!}

  use :: Cooling_Rates                   , only : coolingRateClass
  use :: Cooling_Infall_Radii            , only : coolingInfallRadiusClass
  use :: Cooling_Infall_Torques          , only : coolingInfallTorqueClass
  use :: Cooling_Specific_Angular_Momenta, only : coolingSpecificAngularMomentumClass
  use :: Circumgalactic_Medium_Heating   , only : circumgalacticMediumHeatingClass
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Galactic_Structure_Options      , only : enumerationComponentTypeType
  use :: Hot_Halo_Outflows_Stripping     , only : hotHaloOutflowStrippingClass
  use :: Cooling_Options                 , only : enumerationCoolingFromType

  !![
  <nodeOperator name="nodeOperatorCGMCoolingHeating">
   <description>
    A node operator class that implements inflow/outflow of gas from the \gls{cgm} due to cooling/heating.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMCoolingHeating
     !!{
     A node operator class that implements inflow/outflow of gas from the \gls{cgm} due to cooling/heating.
     !!}
     private
     class           (coolingRateClass                   ), pointer :: coolingRate_                    => null()
     class           (circumgalacticMediumHeatingClass   ), pointer :: circumgalacticMediumHeating_    => null()
     class           (darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_            => null()
     class           (hotHaloOutflowStrippingClass       ), pointer :: hotHaloOutflowStripping_        => null()
     class           (coolingInfallRadiusClass           ), pointer :: coolingInfallRadius_            => null()
     class           (coolingSpecificAngularMomentumClass), pointer :: coolingSpecificAngularMomentum_ => null()
     class           (coolingInfallTorqueClass           ), pointer :: coolingInfallTorque_            => null()
     type            (enumerationComponentTypeType       )          :: component
     type            (enumerationCoolingFromType         )          :: coolingFrom
     logical                                                        :: excessHeatDrivesOutflow
     double precision                                               :: rateMaximumExpulsion
   contains
     final     ::                          cgmCoolingHeatingDestructor
     procedure :: differentialEvolution => cgmCoolingHeatingDifferentialEvolution
  end type nodeOperatorCGMCoolingHeating
  
  interface nodeOperatorCGMCoolingHeating
     !!{
     Constructors for the \refClass{nodeOperatorCGMCoolingHeating} node operator class.
     !!}
     module procedure cgmCoolingHeatingConstructorParameters
     module procedure cgmCoolingHeatingConstructorInternal
  end interface nodeOperatorCGMCoolingHeating
  
contains
  
  function cgmCoolingHeatingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMCoolingHeating} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode
    use :: Cooling_Options           , only : enumerationCoolingFromEncode
    implicit none
    type            (nodeOperatorCGMCoolingHeating      )                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (coolingRateClass                   ), pointer       :: coolingRate_
    class           (circumgalacticMediumHeatingClass   ), pointer       :: circumgalacticMediumHeating_
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class           (hotHaloOutflowStrippingClass       ), pointer       :: hotHaloOutflowStripping_
    class           (coolingInfallRadiusClass           ), pointer       :: coolingInfallRadius_
    class           (coolingSpecificAngularMomentumClass), pointer       :: coolingSpecificAngularMomentum_
    class           (coolingInfallTorqueClass           ), pointer       :: coolingInfallTorque_
    type            (varying_string                     )                :: component                      , coolingFrom
    logical                                                              :: excessHeatDrivesOutflow
    double precision                                                     :: rateMaximumExpulsion

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
    <inputParameter>
      <name>excessHeatDrivesOutflow</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>Specifies whether heating of the halo in excess of its cooling rate will drive an outflow from the halo.</description>
    </inputParameter>
    <inputParameter>
      <name>rateMaximumExpulsion</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the maximum rate at which mass can be expelled from the hot halo in units of the inverse halo dynamical time.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="coolingRate"                    name="coolingRate_"                    source="parameters"/>
    <objectBuilder class="coolingInfallRadius"            name="coolingInfallRadius_"            source="parameters"/>
    <objectBuilder class="coolingSpecificAngularMomentum" name="coolingSpecificAngularMomentum_" source="parameters"/>
    <objectBuilder class="coolingInfallTorque"            name="coolingInfallTorque_"            source="parameters"/>
    <objectBuilder class="circumgalacticMediumHeating"    name="circumgalacticMediumHeating_"    source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
    <objectBuilder class="hotHaloOutflowStripping"        name="hotHaloOutflowStripping_"        source="parameters"/>
    !!]
    self=nodeOperatorCGMCoolingHeating(enumerationComponentTypeEncode(char(component),includesPrefix=.false.),enumerationCoolingFromEncode(char(coolingFrom),includesPrefix=.false.),excessHeatDrivesOutflow,rateMaximumExpulsion,coolingRate_,coolingInfallRadius_,coolingSpecificAngularMomentum_,coolingInfallTorque_,circumgalacticMediumHeating_,hotHaloOutflowStripping_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRate_"                   />
    <objectDestructor name="coolingInfallRadius_"           />
    <objectDestructor name="coolingSpecificAngularMomentum_"/>
    <objectDestructor name="coolingInfallTorque_"           />
    <objectDestructor name="circumgalacticMediumHeating_"   />
    <objectDestructor name="darkMatterHaloScale_"           />
    <objectDestructor name="hotHaloOutflowStripping_"       />
    !!]
    return
  end function cgmCoolingHeatingConstructorParameters

  function cgmCoolingHeatingConstructorInternal(component,coolingFrom,excessHeatDrivesOutflow,rateMaximumExpulsion,coolingRate_,coolingInfallRadius_,coolingSpecificAngularMomentum_,coolingInfallTorque_,circumgalacticMediumHeating_,hotHaloOutflowStripping_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMCoolingHeating} node operator class.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDisk, componentTypeSpheroid, componentTypeNone
    use :: Error                     , only : Error_Report
    implicit none
    type            (nodeOperatorCGMCoolingHeating      )                        :: self
    type            (enumerationComponentTypeType       ), intent(in   )         :: component
    type            (enumerationCoolingFromType         ), intent(in   )         :: coolingFrom
    logical                                              , intent(in   )         :: excessHeatDrivesOutflow
    double precision                                     , intent(in   )         :: rateMaximumExpulsion
    class           (coolingRateClass                   ), intent(in   ), target :: coolingRate_
    class           (coolingInfallRadiusClass           ), intent(in   ), target :: coolingInfallRadius_
    class           (coolingSpecificAngularMomentumClass), intent(in   ), target :: coolingSpecificAngularMomentum_
    class           (circumgalacticMediumHeatingClass   ), intent(in   ), target :: circumgalacticMediumHeating_
    class           (coolingInfallTorqueClass           ), intent(in   ), target :: coolingInfallTorque_
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class           (hotHaloOutflowStrippingClass       ), intent(in   ), target :: hotHaloOutflowStripping_
    !![
    <constructorAssign variables="component, coolingFrom, excessHeatDrivesOutflow, rateMaximumExpulsion, *coolingRate_, *coolingInfallRadius_, *coolingSpecificAngularMomentum_, *coolingInfallTorque_, *circumgalacticMediumHeating_, *hotHaloOutflowStripping_, *darkMatterHaloScale_"/>
    !!]

    if     (                                                                                                            &
         &   component /= componentTypeNone                                                                             &
         &  .and.                                                                                                       &
         &   component /= componentTypeDisk                                                                             &
         &  .and.                                                                                                       &
         &   component /= componentTypeSpheroid                                                                         &
         & ) call Error_Report("only 'disk' 'spheroid', and 'none' components are supported"//{introspection:location})
    return
  end function cgmCoolingHeatingConstructorInternal

  subroutine cgmCoolingHeatingDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorCGMCoolingHeating} node operator class.
    !!}
    implicit none
    type(nodeOperatorCGMCoolingHeating), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRate_"                   />
    <objectDestructor name="self%coolingInfallRadius_"           />
    <objectDestructor name="self%coolingSpecificAngularMomentum_"/>
    <objectDestructor name="self%coolingInfallTorque_"           />
    <objectDestructor name="self%circumgalacticMediumHeating_"   />
    <objectDestructor name="self%darkMatterHaloScale_"           />
    <objectDestructor name="self%hotHaloOutflowStripping_"       />
    !!]
    return
  end subroutine cgmCoolingHeatingDestructor

  subroutine cgmCoolingHeatingDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform inflow/outflow of \gls{cgm} gas due to cooling/heating.
    !!}
    use :: Abundances_Structure         , only : abundances            , operator(*)
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Cooling_Options              , only : coolingFromCurrentNode, coolingFromFormationNode
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo  , nodeComponentDisk       , nodeComponentSpheroid
    use :: Galactic_Structure_Options   , only : componentTypeDisk     , componentTypeSpheroid   , componentTypeNone
    use :: Error                        , only : Error_Report
    implicit none
    class           (nodeOperatorCGMCoolingHeating), intent(inout), target  :: self
    type            (treeNode                     ), intent(inout), target  :: node
    logical                                        , intent(inout)          :: interrupt
    procedure       (interruptTask                ), intent(inout), pointer :: functionInterrupt
    integer                                        , intent(in   )          :: propertyType
    type            (treeNode                     )               , pointer :: nodeCooling
    class           (nodeComponentHotHalo         )               , pointer :: hotHalo                    , hotHaloCooling
    class           (nodeComponentDisk            )               , pointer :: disk
    class           (nodeComponentSpheroid        )               , pointer :: spheroid
    type            (abundances                   )                         :: rateMassAbundancesCooling  , rateMassAbundancesOutflow
    type            (chemicalAbundances           )                         :: rateMassChemicalsCooling   , rateMassChemicalsOutflow
    logical                                                                 :: hasAngularMomentum         , hasOuterRadius            , &
         &                                                                     hasChemicals
    double precision                                                        :: rateMassCooling            , rateMassHeating           , &
         &                                                                     rateAngularMomentumCooling , rateAngularMomentumOutflow, &
         &                                                                     rateMassOutflow            , radiusInfall              , &
         &                                                                     fractionLossAngularMomentum
    !$GLC attributes unused :: propertyType

    ! Ignore cases with unphysical mass.
    hotHalo => node%hotHalo()
    if (hotHalo%mass() <= 0.0d0) return
    ! Find available properties.
    hasAngularMomentum=hotHalo%angularMomentumIsGettable()
    hasOuterRadius    =hotHalo%    outerRadiusIsGettable()
    hasChemicals      =hotHalo%      chemicalsIsGettable()
    ! Ignore cases with unphysical angular momentum or outer radius.
    if     (                                                               &
         &   (hasAngularMomentum .and. hotHalo%angularMomentum() <= 0.0d0) &
         &  .or.                                                           &
         &   (hasOuterRadius     .and. hotHalo%    outerRadius() <= 0.0d0) &
         & ) return
    ! Compute the rate of cooling.
    rateMassCooling=+self%coolingRate_                %rate          (node)
    ! Compute the heating rate.
    rateMassHeating=+self%circumgalacticMediumHeating_%heatingRate   (node)    &
         &          /self%darkMatterHaloScale_        %velocityVirial(node)**2
    ! Determine if we have net heating or cooling.
    if (rateMassHeating > rateMassCooling) then
       if (self%excessHeatDrivesOutflow) then
          rateMassOutflow           =min(                                                          &
               &                         +                             rateMassHeating             &
               &                         -                             rateMassCooling           , &
               &                         +self   %                     rateMaximumExpulsion        &
               &                         *hotHalo                     %mass                (    )  &
               &                         /self   %darkMatterHaloScale_%timescaleDynamical  (node)  &
               &                        )
          ! Get the rate of change of abundances, chemicals, and angular momentum.
          rateMassAbundancesOutflow        =hotHalo%abundances     ()*(rateMassOutflow/hotHalo%mass())
          if (hasAngularMomentum) &
               & rateAngularMomentumOutflow=hotHalo%angularMomentum()*(rateMassOutflow/hotHalo%mass())
          if (hasChemicals      ) &
               & rateMassChemicalsOutflow  =hotHalo%chemicals      ()*(rateMassOutflow/hotHalo%mass())
          call        hotHalo%           massRate(-rateMassOutflow           )
          call        hotHalo%     abundancesRate(-rateMassAbundancesOutflow )
          if (hasAngularMomentum) &
               & call hotHalo%angularMomentumRate(-rateAngularMomentumOutflow)
          if (hasChemicals      ) &
               & call hotHalo%      chemicalsRate(-rateMassChemicalsOutflow  )
          ! If this node is a satellite and stripped gas is being tracked, move mass and abundances to the stripped reservoir.
          if (.not.self%hotHaloOutflowStripping_%neverStripped(node)) then
             call hotHalo%      strippedMassRate(rateMassOutflow          )
             call hotHalo%strippedAbundancesRate(rateMassAbundancesOutflow)
             call hotHalo% strippedChemicalsRate(rateMassChemicalsOutflow )
          end if
          ! Trigger an event to allow other processes to respond to this outflow.
          !![
	  <eventHook name="hotHaloMassEjection">
	    <import>
	      <module name="Galacticus_Nodes" symbols="nodeComponentHotHalo"/>
	    </import>
	    <interface>
	      class           (nodeComponentHotHalo), intent(inout) :: hotHalo
	      double precision                      , intent(in   ) :: rateMassOutflow
	    </interface>
	    <callWith>hotHalo,rateMassOutflow</callWith>
	  </eventHook>
          !!]
       end if
    else if (rateMassCooling > rateMassHeating) then
       ! Reduce the cooling rate by the heating rate.
       rateMassCooling=max(0.0d0,rateMassCooling-rateMassHeating)
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
       !! Angular momentum.
       if (hasAngularMomentum) then
          ! Find the infall radius.
          radiusInfall              =+self          %coolingInfallRadius_           %radius                     (node                    )
          ! Find the fraction of angular momentum lost during infall.
          fractionLossAngularMomentum=self          %coolingInfallTorque_           %fractionAngularMomentumLoss(node                    )
          ! Calculate cooling rates of other quantities.
          rateAngularMomentumCooling=+                                               rateMassCooling                                       &
               &                     *self          %coolingSpecificAngularMomentum_%angularMomentumSpecific    (nodeCooling,radiusInfall)
       else
          rateAngularMomentumCooling=+0.0d0
       end if
       !! Abundances.
       rateMassAbundancesCooling    =+                                               rateMassCooling                                       &
            &                        *hotHaloCooling                                %abundances                 (                        ) &
            &                        /hotHaloCooling                                %mass                       (                        )
       !! Chemicals.
       if (hasChemicals      ) then
          rateMassChemicalsCooling  = hotHaloCooling                                %chemicals                  (                        )
          call rateMassChemicalsCooling%scale(                                  &
               &                              -               rateMassCooling   &
               &                              /hotHaloCooling%mass           () &
               &                             )
       end if
       ! Apply cooling.
       call        hotHalo %massRate           (-rateMassCooling                                                                           )
       if (hasAngularMomentum) &
            & call hotHalo %angularMomentumRate(-rateAngularMomentumCooling                                                                )
       call        hotHalo %abundancesRate     (-rateMassAbundancesCooling                                                                 )
       if (hasChemicals      ) &
            & call hotHalo %chemicalsRate      (+rateMassChemicalsCooling                                                                  )
       if      (self%component == componentTypeDisk    ) then
          disk     => node%disk    ()
          call     disk    %massGasRate        (+rateMassCooling                                               ,interrupt,functionInterrupt)
          call     disk    %abundancesGasRate  (+rateMassAbundancesCooling                                     ,interrupt,functionInterrupt)
          if (hasAngularMomentum) &
            & call disk    %angularMomentumRate(+rateAngularMomentumCooling*(1.0d0-fractionLossAngularMomentum),interrupt,functionInterrupt)
       else if (self%component == componentTypeSpheroid) then
          spheroid => node%spheroid()
          call     spheroid%massGasRate        (+rateMassCooling                                               ,interrupt,functionInterrupt)
          call     spheroid%abundancesGasRate  (+rateMassAbundancesCooling                                     ,interrupt,functionInterrupt)
          if (hasAngularMomentum) &
            & call spheroid%angularMomentumRate(+rateAngularMomentumCooling*(1.0d0-fractionLossAngularMomentum),interrupt,functionInterrupt)
       else if (self%component /= componentTypeNone    ) then
          call Error_Report('unexpected component - this should not happen'//{introspection:location})
       end if
       ! Trigger an event to allow other processes to respond to this inflow.
       !![
       <eventHook name="hotHaloMassInflow">
	 <import>
	   <module name="Galacticus_Nodes" symbols="nodeComponentHotHalo"/>
	 </import>
	 <interface>
	   class           (nodeComponentHotHalo), intent(inout) :: hotHalo
	   double precision                      , intent(in   ) :: rateMassCooling
	 </interface>
	 <callWith>hotHalo,rateMassCooling</callWith>
       </eventHook>
       !!]
    end if
    return
  end subroutine cgmCoolingHeatingDifferentialEvolution
  

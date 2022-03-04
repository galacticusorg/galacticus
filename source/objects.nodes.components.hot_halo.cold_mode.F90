!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements an extension to the standard hot halo node component which
supports a cold mode reservoir.
!!}

module Node_Component_Hot_Halo_Cold_Mode
  !!{
  Implements an extension to the standard hot halo node component which supports a cold mode
  reservoir.
  !!}
  use :: Accretion_Halos               , only : accretionHaloClass
  use :: Cooling_Cold_Mode_Infall_Rates, only : coldModeInfallRateClass
  use :: Cosmology_Parameters          , only : cosmologyParametersClass
  use :: Galactic_Structure            , only : galacticStructureClass
  implicit none
  private
  public :: Node_Component_Hot_Halo_Cold_Mode_Initialize       , Node_Component_Hot_Halo_Cold_Mode_Rate_Compute       , &
       &    Node_Component_Hot_Halo_Cold_Mode_Scale_Set        , Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize    , &
       &    Node_Component_Hot_Halo_Cold_Mode_Node_Merger      , Node_Component_Hot_Halo_Cold_Mode_Formation          , &
       &    Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize, Node_Component_Hot_Halo_Cold_Mode_Thread_Uninitialize, &
       &    Node_Component_Hot_Halo_Cold_Mode_State_Store      , Node_Component_Hot_Halo_Cold_Mode_State_Restore

  !![
  <component>
   <class>hotHalo</class>
   <name>coldMode</name>
   <extends>
     <class>hotHalo</class>
     <name>standard</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>massCold</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of cold-mode gas in the hot halo."/>
    </property>
    <property>
      <name>abundancesCold</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the cold-mode of the hot halo."/>
    </property>
    <property>
      <name>angularMomentumCold</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of cold-mode gas in the hot halo."/>
    </property>
    <property>
      <name>massTotal</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
      <getFunction>Node_Component_Hot_Halo_Cold_Mode_Mass_Total</getFunction>
    </property>
   </properties>
   <functions>objects.nodes.components.hot_halo.cold_mode.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(accretionHaloClass      ), pointer :: accretionHalo_
  class(coldModeInfallRateClass ), pointer :: coldModeInfallRate_
  class(cosmologyParametersClass), pointer :: cosmologyParameters_
  class(galacticStructureClass  ), pointer :: galacticStructure_
  !$omp threadprivate(accretionHalo_,coldModeInfallRate_,cosmologyParameters_,galacticStructure_)

  ! Options controlling the behavior of the cold mode gas.
  logical :: hotHaloOutflowToColdMode

  ! Internal count of abundances.
  integer :: abundancesCount

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Initialize(parameters_)
    !!{
    Initializes the tree node hot halo methods module.
    !!}
    use :: Abundances_Structure, only : Abundances_Property_Count
    use :: Galacticus_Nodes    , only : defaultHotHaloComponent  , nodeComponentHotHaloColdMode
    use :: Input_Parameters    , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters             ), intent(inout) :: parameters_
    type(nodeComponentHotHaloColdMode)                :: hotHalo

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Hot_Halo_Cold_Mode_Initialize)
    if (defaultHotHaloComponent%coldModeIsActive()) then
       ! Get numbers of abundance properties.
       abundancesCount=Abundances_Property_Count()
       ! Determine whether outflows go to the cold mode.
       !![
       <inputParameter>
         <name>hotHaloOutflowToColdMode</name>
         <defaultValue>.false.</defaultValue>
         <description>Specifies whether or not outflows from galaxies are returned to the cold or hot modes in the hot halo.</description>
         <source>parameters_</source>
       </inputParameter>
       !!]
       ! Bind the outflow return function if outflow returns to the cold mode. (If it does not, do
       ! not bind any function and let the parent class handle this behavior.)
       if (hotHaloOutflowToColdMode) call hotHalo%outflowReturnFunction(Node_Component_Hot_Halo_Cold_Mode_Outflow_Return)
    end if
    !$omp end critical (Node_Component_Hot_Halo_Cold_Mode_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node hot halo cold mode methods module.
    !!}
    use :: Events_Hooks                                     , only : nodePromotionEvent       , satelliteMergerEvent    , openMPThreadBindingAtLevel, dependencyRegEx, &
         &                                                           dependencyDirectionAfter
    use :: Galacticus_Nodes                                 , only : defaultHotHaloComponent
    use :: Hot_Halo_Cold_Mode_Density_Core_Radii            , only : hotHaloColdModeCoreRadii
    use :: Input_Parameters                                 , only : inputParameter           , inputParameters
    use :: Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks, only : darkMatterHaloScale_     , hotHaloColdModeCoreRadii_
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies

    if (defaultHotHaloComponent%coldModeIsActive()) then
       !![
       <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters_"/>
       <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters_"/>
       <objectBuilder class="accretionHalo"            name="accretionHalo_"            source="parameters_"/>
       <objectBuilder class="coldModeInfallRate"       name="coldModeInfallRate_"       source="parameters_"/>
       <objectBuilder class="hotHaloColdModeCoreRadii" name="hotHaloColdModeCoreRadii_" source="parameters_"/>
       <objectBuilder class="galacticStructure"        name="galacticStructure_"        source="parameters_"/>
       !!]
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call nodePromotionEvent  %attach(defaultHotHaloComponent,nodePromotion  ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloColdMode'                          )
       call satelliteMergerEvent%attach(defaultHotHaloComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentHotHaloColdMode',dependencies=dependencies)
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Uninitialize()
    !!{
    Uninitializes the tree node hot halo cold mode methods module.
    !!}
    use :: Events_Hooks                                     , only : nodePromotionEvent       , satelliteMergerEvent
    use :: Galacticus_Nodes                                 , only : defaultHotHaloComponent
    use :: Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks, only : darkMatterHaloScale_     , hotHaloColdModeCoreRadii_
    implicit none

    if (defaultHotHaloComponent%coldModeIsActive()) then
       !![
       <objectDestructor name="cosmologyParameters_"     />
       <objectDestructor name="darkMatterHaloScale_"     />
       <objectDestructor name="accretionHalo_"           />
       <objectDestructor name="coldModeInfallRate_"      />
       <objectDestructor name="hotHaloColdModeCoreRadii_"/>
       <objectDestructor name="galacticStructure_"       />
       !!]
       if (nodePromotionEvent  %isAttached(defaultHotHaloComponent,nodePromotion  )) call nodePromotionEvent  %detach(defaultHotHaloComponent,nodePromotion  )
       if (satelliteMergerEvent%isAttached(defaultHotHaloComponent,satelliteMerger)) call satelliteMergerEvent%detach(defaultHotHaloComponent,satelliteMerger)
   end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Uninitialize

  subroutine Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes(node,massRate,interrupt,interruptProcedure)
    !!{
    Push mass through the cooling pipes (along with appropriate amounts of metals and angular momentum) at the given rate.
    !!}
    use :: Abundances_Structure                 , only : abundances             , operator(*)
    use :: Galacticus_Error                     , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                     , only : interruptTask          , nodeComponentHotHalo, nodeComponentHotHaloColdMode      , treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : currentNode            , formationNode       , hotHaloAngularMomentumLossFraction, hotHaloCoolingFromNode
    use :: Numerical_Constants_Math             , only : Pi
    implicit none
    type            (treeNode                    ), intent(inout)          , target  :: node
    double precision                              , intent(in   )                    :: massRate
    logical                                       , intent(inout), optional          :: interrupt
    procedure       (interruptTask               ), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode                    )                         , pointer :: nodeCoolingFrom
    class           (nodeComponentHotHalo        )                         , pointer :: hotHaloCoolingFrom        , hotHalo
    type            (abundances                  ), save                             :: abundancesCoolingRate
    !$omp threadprivate(abundancesCoolingRate)
    double precision                                                                 :: angularMomentumCoolingRate

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Ignore zero rates.
       if (massRate /= 0.0d0 .and. hotHalo%massCold() > 0.0d0 .and. hotHalo%angularMomentumCold() > 0.0d0) then
          ! Remove mass from the hot component.
          call hotHalo%massColdRate(-massRate)
          ! Pipe the mass rate to whichever component claimed it.
          if (hotHalo%hotHaloCoolingMassRateIsAttached()) &
               & call hotHalo%hotHaloCoolingMassRate(+massRate,interrupt,interruptProcedure)
          ! Find the node to use for cooling calculations.
          select case (hotHaloCoolingFromNode)
          case (currentNode  )
             nodeCoolingFrom => node
          case (formationNode)
             nodeCoolingFrom => node%formationNode
          case default
             nodeCoolingFrom => null()
             call Galacticus_Error_Report('unknown cooling node'//{introspection:location})
          end select
          ! Compute the infall rate of angular momentum.
          angularMomentumCoolingRate=massRate*hotHalo%angularMomentumCold()/hotHalo%massCold()
          call hotHalo%angularMomentumColdRate(-angularMomentumCoolingRate)
          ! Pipe the cooling rate to which ever component claimed it.
          if (hotHalo%hotHaloCoolingAngularMomentumRateIsAttached()) &
               & call hotHalo%hotHaloCoolingAngularMomentumRate(sign(+angularMomentumCoolingRate*(1.0d0-hotHaloAngularMomentumLossFraction),massRate),interrupt,interruptProcedure)
          ! Get the rate of change of abundances.
          hotHaloCoolingFrom => nodeCoolingFrom   %hotHalo       ()
          abundancesCoolingRate=hotHaloCoolingFrom%abundancesCold()
          abundancesCoolingRate=massRate*abundancesCoolingRate/hotHaloCoolingFrom%massCold()
          call hotHalo%abundancesColdRate(-abundancesCoolingRate)
          ! Pipe the cooling rate to which ever component claimed it.
          if (hotHalo%hotHaloCoolingAbundancesRateIsAttached()) &
               & call hotHalo%hotHaloCoolingAbundancesRate(+abundancesCoolingRate,interrupt,interruptProcedure)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes

  !![
  <rateComputeTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the hot halo node mass rate of change.
    !!}
    use :: Abundances_Structure                             , only : abs
    use :: Accretion_Halos                                  , only : accretionModeCold
    use :: Galactic_Structure_Options                       , only : componentTypeColdHalo            , coordinateSystemSpherical         , massTypeGaseous
    use :: Galacticus_Nodes                                 , only : defaultHotHaloComponent          , interruptTask                     , nodeComponentBasic, nodeComponentHotHalo, &
          &                                                          nodeComponentHotHaloColdMode     , propertyInactive                  , treeNode          , nodeComponentSpin
    use :: Node_Component_Hot_Halo_Standard_Data            , only : hotHaloAngularMomentumAlwaysGrows, outerRadiusOverVirialRadiusMinimum
    use :: Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks, only : darkMatterHaloScale_
    use :: Numerical_Constants_Math                         , only : Pi
    implicit none
    type            (treeNode            ), intent(inout)          :: node
    logical                               , intent(inout)          :: interrupt
    procedure       (interruptTask       ), intent(inout), pointer :: interruptProcedure
    integer                               , intent(in   )          :: propertyType
    class           (nodeComponentSpin   )               , pointer :: spin
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                                               :: angularMomentumAccretionRate, densityAtOuterRadius , &
         &                                                            massAccretionRate           , massLossRate         , &
         &                                                            outerRadius                 , outerRadiusGrowthRate, &
         &                                                            gasMass                     , infallRate

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%coldModeIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Find the rate of gas mass accretion onto the halo.
    massAccretionRate=accretionHalo_%accretionRate(node,accretionModeCold)
    ! Get the basic component.
    basic => node%basic()
    ! Apply accretion rates.
    call hotHalo%massColdRate(massAccretionRate,interrupt,interruptProcedure)
    ! Next compute the cold mode infall rate in this halo.
    infallRate=coldModeInfallRate_%infallRate(node)
    ! Pipe the cooling rate to which ever component claimed it.
    call Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes(node,infallRate,interrupt,interruptProcedure)
    ! Get the rate at which abundances are accreted onto this halo.
    call hotHalo%abundancesColdRate(accretionHalo_%accretionRateMetals(node,accretionModeCold),interrupt,interruptProcedure)
    ! Next block of tasks occur only if the accretion rate is non-zero.
    if (massAccretionRate > 0.0d0) then
       ! Compute the rate of accretion of angular momentum.
       spin => node%spin()
       angularMomentumAccretionRate=+spin%angularMomentumGrowthRate() &
            &                       *     massAccretionRate           &
            &                       /basic%accretionRate           ()
           if (hotHaloAngularMomentumAlwaysGrows) angularMomentumAccretionRate=abs(angularMomentumAccretionRate)
       call hotHalo%angularMomentumColdRate(angularMomentumAccretionRate,interrupt,interruptProcedure)
    end if
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Test whether this halo is a satellite or not.
       if (node%isSatellite()) then
          ! For satellites, get the current ram pressure stripping radius for this hot halo.
          outerRadiusGrowthRate=hotHalo%outerRadiusGrowthRate()
          outerRadius          =hotHalo%outerRadius          ()
          gasMass              =hotHalo%massCold             ()
          if     (                                                                                                    &
               &   outerRadiusGrowthRate /= 0.0d0                                                                     &
               &  .and.                                                                                               &
               &   gasMass               >  0.0d0                                                                     &
               &  .and.                                                                                               &
               &   outerRadius           <=                                   darkMatterHaloScale_%radiusVirial(node) &
               &  .and.                                                                                               &
               &   outerRadius           > outerRadiusOverVirialRadiusMinimum*darkMatterHaloScale_%radiusVirial(node) &
               & ) then
             ! The ram pressure stripping radius is within the outer radius. Remove mass from the cold mode halo at the appropriate rate.
             densityAtOuterRadius=galacticStructure_%density(node,[outerRadius,0.0d0,0.0d0],coordinateSystemSpherical,componentTypeColdHalo,massTypeGaseous)
             ! Compute the mass loss rate.
             massLossRate=4.0d0*Pi*densityAtOuterRadius*outerRadius**2*outerRadiusGrowthRate
             ! Adjust the rates.
             ! Mass.
             call hotHalo%           massColdRate(+                              massLossRate        ,interrupt,interruptProcedure)
             ! Angular momentum.
             call hotHalo%angularMomentumColdRate(+hotHalo%angularMomentumCold()*massLossRate/gasMass,interrupt,interruptProcedure)
             ! Metal abundances.
             call hotHalo%     abundancesColdRate(+hotHalo%abundancesCold     ()*massLossRate/gasMass,interrupt,interruptProcedure)
             ! Mass.
             call hotHalo%       strippedMassRate(-                              massLossRate        ,interrupt,interruptProcedure)
             ! Metal abundances.
             call hotHalo% strippedAbundancesRate(-hotHalo%abundancesCold     ()*massLossRate/gasMass,interrupt,interruptProcedure)
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Rate_Compute

  subroutine Node_Component_Hot_Halo_Cold_Mode_Outflow_Return(self,interrupt,interruptProcedure)
    !!{
    Return outflowed gas to the cold mode reservoir.
    !!}
    use :: Abundances_Structure                             , only : abundances              , max                      , operator(*)
    use :: Galactic_Structure_Options                       , only : componentTypeColdHalo   , coordinateSystemSpherical, massTypeGaseous
    use :: Galacticus_Error                                 , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                                 , only : interruptTask           , nodeComponentBasic       , nodeComponentHotHaloColdMode, nodeComponentHotHaloStandard, &
          &                                                          treeNode
    use :: Node_Component_Hot_Halo_Standard_Data            , only : hotHaloOutflowReturnRate, starveSatellites
    use :: Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks, only : darkMatterHaloScale_
    use :: Numerical_Constants_Astronomical                 , only : gigaYear                , megaParsec
    use :: Numerical_Constants_Math                         , only : Pi
    use :: Numerical_Constants_Prefixes                     , only : kilo
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout)          :: self
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    type            (treeNode                    ), pointer                :: node
    class           (nodeComponentBasic          ), pointer                :: basic
    double precision                                                       :: outflowedMass            , massReturnRate, &
         &                                                                    angularMomentumReturnRate, radiusVirial  , &
         &                                                                    densityAtOuterRadius     , densityMinimum, &
         &                                                                    outerRadius
    type            (abundances                  ), save                   :: abundancesReturnRate
    !$omp threadprivate(abundancesReturnRate)

    select type (self)
    class is (nodeComponentHotHaloColdMode)
       ! Get the hosting node.
       node => self%hostNode
       ! Next tasks occur only for systems in which outflowed gas is being recycled.
       massReturnRate=0.0d0
       if (.not.starveSatellites.or..not.node%isSatellite()) then
          outflowedMass            =self%outflowedMass()
          massReturnRate           =hotHaloOutflowReturnRate*outflowedMass                  /darkMatterHaloScale_%timescaleDynamical(node)
          angularMomentumReturnRate=hotHaloOutflowReturnRate*self%outflowedAngularMomentum()/darkMatterHaloScale_%timescaleDynamical(node)
          abundancesReturnRate     =hotHaloOutflowReturnRate*self%outflowedAbundances     ()/darkMatterHaloScale_%timescaleDynamical(node)
          call self%           outflowedMassRate(-           massReturnRate,interrupt,interruptProcedure)
          call self%                massColdRate(+           massReturnRate,interrupt,interruptProcedure)
          call self%outflowedAngularMomentumRate(-angularMomentumReturnRate,interrupt,interruptProcedure)
          call self%     angularMomentumColdRate(+angularMomentumReturnRate,interrupt,interruptProcedure)
          call self%     outflowedAbundancesRate(-     abundancesReturnRate,interrupt,interruptProcedure)
          call self%          abundancesColdRate(+     abundancesReturnRate,interrupt,interruptProcedure)
       end if
       ! The outer radius must be increased as the halo fills up with gas.
       outerRadius =self%outerRadius()
       radiusVirial=darkMatterHaloScale_%radiusVirial(node)
       if (outerRadius < radiusVirial) then
          densityAtOuterRadius=galacticStructure_%density(node,[outerRadius,0.0d0,0.0d0],coordinateSystemSpherical,componentTypeColdHalo,massTypeGaseous)
          ! If the outer radius and density are non-zero we can expand the outer radius at a rate determined by the current
          ! density profile.
          if (outerRadius > 0.0d0 .and. densityAtOuterRadius > 0.0d0) then
             ! Limit the density at the outer radius to one third of the mean virial density (for baryons, assuming a
             ! universal baryon fraction) to prevent arbitrarily rapid growth of the outer radius in halos containing almost
             ! no gas.
             basic => node%basic()
             densityMinimum=(cosmologyParameters_%omegaBaryon()/cosmologyParameters_%omegaMatter())*basic%mass()/radiusVirial**3/4.0d0/Pi
             call self%outerRadiusRate(                           &
                  &                     massReturnRate            &
                  &                    /4.0d0                     &
                  &                    /Pi                        &
                  &                    /outerRadius**2            &
                  &                    /max(                      &
                  &                         densityAtOuterRadius, &
                  &                         densityMinimum        &
                  &                        )                      &
                  &                   )
          ! Otherwise, if we have a positive rate of mass return, simply grow the radius at the virial velocity.
          else if (massReturnRate > 0.0d0) then
             ! Force some growth here so the radius is not trapped at zero.
             call self%outerRadiusRate(darkMatterHaloScale_%velocityVirial(node)*kilo*gigaYear/megaParsec)
          end if
       end if
    class default
       call Galacticus_Error_Report('this function should not be called for non-coldMode class hot halo components'//{introspection:location})
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Outflow_Return

  !![
  <scaleSetTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure                             , only : unitAbundances
    use :: Galacticus_Nodes                                 , only : nodeComponentBasic     , nodeComponentHotHalo, nodeComponentHotHaloColdMode, treeNode, &
         &                                                           defaultHotHaloComponent
    use :: Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks, only : darkMatterHaloScale_
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                      , parameter              :: scaleMassRelative   =1.0d-3
    double precision                      , parameter              :: scaleRadiusRelative =1.0d+0
    double precision                                               :: massVirial                 , radiusVirial, &
         &                                                            velocityVirial

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%coldModeIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of the cold mode class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! The the basic component.
       basic => node%basic()
       ! Get virial properties.
       massVirial    =basic%mass()
       radiusVirial  =darkMatterHaloScale_%radiusVirial  (node)
       velocityVirial=darkMatterHaloScale_%velocityVirial(node)
       call    hotHalo%           massColdScale(               massVirial                            *scaleMassRelative)
       call    hotHalo%     abundancesColdScale(unitAbundances*massVirial                            *scaleMassRelative)
       call    hotHalo%angularMomentumColdScale(               massVirial*radiusVirial*velocityVirial*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Scale_Set

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize</unitName>
   <after>Node_Component_Hot_Halo_Standard_Tree_Initialize</after>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize(node)
    !!{
    Initialize the contents of the hot halo component for any sub-resolution accretion (i.e. the gas that would have been
    accreted if the merger tree had infinite resolution).
    !!}
    use :: Accretion_Halos , only : accretionModeCold
    use :: Galacticus_Nodes, only : defaultHotHaloComponent  , nodeComponentBasic, nodeComponentHotHalo, nodeEvent, &
          &                         nodeEventSubhaloPromotion, treeNode          , nodeComponentSpin
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    class           (nodeComponentSpin   )               , pointer :: spin
    class           (nodeEvent           )               , pointer :: event
    double precision                                               :: angularMomentum, coldModeMass

    ! If the node has a child or the standard hot halo is not active, then return immediately.
    if (associated(node%firstChild).or..not.defaultHotHaloComponent%coldModeIsActive()) return
    ! Search for a subhalo promotion events associated with this node.
    event => node%event
    do while (associated(event))
       ! Check if this event:
       !  a) is a subhalo promotion event;
       !  b) has no associated task (which means this is the node being promoted to, not the node being promoted itself).
       ! Do not assign any mass to such nodes, as they should receive gas from the node which is promoted to them.
       select type (event)
       type is (nodeEventSubhaloPromotion)
          if (.not.associated(event%task)) return
       end select
       event => event%next
    end do
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Get the mass of cold mode gas accreted.
    coldModeMass=accretionHalo_%accretedMass(node,accretionModeCold)
    ! If non-zero, then create a hot halo component and add to it.
    if (coldModeMass > 0.0d0) then
       ! Ensure that it is of unspecified class.
       hotHalo => node%hotHalo(autoCreate=.true.)
       basic   => node%basic  (                 )
       spin    => node%spin   (                 )
       call hotHalo%massColdSet(coldModeMass)
       ! Also add the appropriate angular momentum.
       angularMomentum=+      coldModeMass      &
            &          *spin %angularMomentum() &
            &          /basic%mass           ()
       call hotHalo%angularMomentumColdSet(angularMomentum)
       ! Add the appropriate abundances.
       call hotHalo%abundancesColdSet(accretionHalo_%accretedMassMetals(node,accretionModeCold))
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize

  !![
  <nodeMergerTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_Node_Merger</unitName>
   <before>Node_Component_Hot_Halo_Standard_Node_Merger</before>
  </nodeMergerTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Node_Merger(node)
    !!{
    Starve {\normalfont \ttfamily node} by transferring its hot halo to its parent.
    !!}
    use :: Abundances_Structure                 , only : abundances                          , operator(*)            , zeroAbundances
    use :: Accretion_Halos                      , only : accretionModeCold                   , accretionModeTotal
    use :: Galactic_Structure_Options           , only : componentTypeAll                    , massTypeBaryonic       , radiusLarge
    use :: Galacticus_Nodes                     , only : nodeComponentBasic                  , nodeComponentHotHalo   , nodeComponentHotHaloColdMode, nodeComponentSpin, &
          &                                              treeNode                            , defaultHotHaloComponent
    use :: Node_Component_Hot_Halo_Standard_Data, only : hotHaloNodeMergerLimitBaryonFraction, starveSatellites
    implicit none
    type            (treeNode            ), intent(inout) :: node
    type            (treeNode            ), pointer       :: nodeParent
    class           (nodeComponentHotHalo), pointer       :: hotHaloParent          , hotHalo
    class           (nodeComponentSpin   ), pointer       :: spinParent
    class           (nodeComponentBasic  ), pointer       :: basicParent            , basic
    double precision                                      :: baryonicMassCurrent    , baryonicMassMaximum   , &
         &                                                   fractionRemove         , massAccretedCold      , &
         &                                                   massAccreted           , massUnaccreted        , &
         &                                                   angularMomentumAccreted, massReaccreted        , &
         &                                                   fractionAccreted
    type            (abundances          ), save          :: massMetalsAccreted     , fractionMetalsAccreted, &
         &                                                   massMetalsReaccreted
    !$omp threadprivate(massMetalsAccreted,fractionMetalsAccreted,massMetalsReaccreted)

    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%coldModeIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of cold mode class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Find the parent node and its hot halo and angular momentum components.
       nodeParent    => node      %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       spinParent    => nodeParent%spin   (                 )
       basicParent   => nodeParent%basic  (                 )
       basic         => node      %basic  (                 )
       ! Since the parent node is undergoing mass growth through this merger we potentially return some of the unaccreted gas to
       ! the hot phase.
       !! First, find the masses of hot and failed mass the node would have if it formed instantaneously.
       massAccretedCold=accretionHalo_%      accretedMass(nodeParent,accretionModeCold )
       massAccreted    =accretionHalo_%      accretedMass(nodeParent,accretionModeTotal)
       massUnaccreted  =accretionHalo_%failedAccretedMass(nodeParent,accretionModeTotal)
       !! Find the fraction of mass that would be successfully accreted.
       fractionAccreted=+  massAccretedCold &
            &           /(                  &
            &             +massAccreted     &
            &             +massUnaccreted   &
            &            )
       !! Find the change in the unaccreted mass.
       massReaccreted=+hotHaloParent   %unaccretedMass() &
            &         *fractionAccreted                  &
            &         *basic           %          mass() &
            &         /basicParent     %          mass()
       !! Reaccrete the gas.
       call hotHaloParent%unaccretedMassSet(hotHaloParent%unaccretedMass()-massReaccreted)
       call hotHaloParent%      massColdSet(hotHaloParent%      massCold()+massReaccreted)
       ! Compute the reaccreted angular momentum.
       angularMomentumAccreted=+            massReaccreted    &
            &                  *spinParent %angularMomentum() &
            &                  /basicParent%mass           ()
       call hotHaloParent%angularMomentumColdSet(hotHaloParent%angularMomentumCold()+angularMomentumAccreted)
       ! Compute the reaccreted metals.
       !! First, find the metal mass the node would have if it formed instantaneously.
       massMetalsAccreted=accretionHalo_%accretedMassMetals(nodeParent,accretionModeCold)
       !! Find the mass fraction of metals that would be successfully accreted.
       fractionMetalsAccreted=+  massMetalsAccreted &
            &                 /(                    &
            &                   +massAccreted       &
            &                   +massUnaccreted     &
            &                  )
       !! Find the change in the unaccreted mass.
       massMetalsReaccreted=+hotHaloParent         %unaccretedMass() &
            &               *fractionMetalsAccreted                  &
            &               *basic                 %          mass() &
            &               /basicParent           %          mass()
       !! Reaccrete the metals.
       call hotHaloParent%abundancesColdSet(hotHaloParent%abundancesCold()+massMetalsReaccreted)
       ! Determine if starvation is to be applied.
       if (starveSatellites) then
          ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate to
          ! this hot halo (and will be moved to the parent at the end of the evolution timestep).
          call hotHaloParent%           massColdSet(                                      &
               &                                     hotHaloParent %massCold           () &
               &                                    +hotHalo       %massCold           () &
               &                                   )
          call hotHaloParent%angularMomentumColdSet(                                      &
               &                                     hotHaloParent %angularMomentumCold() &
               &                                    +hotHalo       %massCold           () &
               &                                    *spinParent    %angularMomentum    () &
               &                                    /basicParent   %mass               () &
               &                                   )
          call hotHalo      %           massColdSet(                                      &
               &                                     0.0d0                                &
               &                                   )
          call hotHalo      %angularMomentumColdSet(                                      &
               &                                     0.0d0                                &
               &                                   )
          call hotHaloParent%     abundancesColdSet(                                      &
               &                                     hotHaloParent %abundancesCold     () &
               &                                    +hotHalo       %abundancesCold     () &
               &                                   )
          call hotHalo      %     abundancesColdSet(                                      &
               &                                     zeroAbundances                       &
               &                                   )
          ! Check if the baryon fraction in the parent hot halo exceeds the universal value. If it does, mitigate this by moving
          ! some of the mass to the failed accretion reservoir.
          if (hotHaloNodeMergerLimitBaryonFraction) then
             baryonicMassMaximum=+basicParent         %mass        () &
                  &              *cosmologyParameters_%omegaBaryon () &
                  &              /cosmologyParameters_%omegaMatter ()
             baryonicMassCurrent= galacticStructure_  %massEnclosed(                                &
                  &                                                 nodeParent                    , &
                  &                                                 radiusLarge                   , &
                  &                                                 massType     =massTypeBaryonic, &
                  &                                                 componentType=componentTypeAll  &
                  &                                                )
             if (baryonicMassCurrent > baryonicMassMaximum .and. hotHaloParent%mass()+hotHaloParent%massCold() > 0.0d0) then
                fractionRemove=min((baryonicMassCurrent-baryonicMassMaximum)/hotHaloParent%massTotal(),1.0d0)
                call hotHaloParent%     unaccretedMassSet(                                                            &
                     &                                     hotHaloParent%unaccretedMass     ()                        &
                     &                                    +hotHaloParent%massCold           ()*       fractionRemove  &
                     &                                   )
                call hotHaloParent%           massColdSet( hotHaloParent%massCold           ()*(1.0d0-fractionRemove))
                call hotHaloParent%angularMomentumColdSet( hotHaloParent%angularMomentumCold()*(1.0d0-fractionRemove))
                call hotHaloParent%     abundancesColdSet( hotHaloParent%abundancesCold     ()*(1.0d0-fractionRemove))
             end if
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Node_Merger

  subroutine satelliteMerger(self,node)
    !!{
    Remove any cold mode gas associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Abundances_Structure                 , only : abundances          , zeroAbundances
    use :: Galacticus_Nodes                     , only : nodeComponentHotHalo, nodeComponentHotHaloColdMode, nodeComponentSpin, nodeComponentBasic, &
         &                                               treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : starveSatellites
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeHost
    class(nodeComponentBasic  ), pointer       :: basicHost
    class(nodeComponentHotHalo), pointer       :: hotHaloHost, hotHalo
    class(nodeComponentSpin   ), pointer       :: spinHost
    !$GLC attributes unused :: self

    ! Return immediately if satellites are starved, as in that case there is no hot halo to transfer.
    if (starveSatellites) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Find the node with which to merge.
       nodeHost    => node    %mergesWith(                 )
       hotHaloHost => nodeHost%hotHalo   (autoCreate=.true.)
       basicHost   => nodeHost%basic     (                 )
       spinHost    => nodeHost%spin      (                 )
       ! Move the cold mode to the host.
       call hotHaloHost%               massSet(                                 &
            &                                   hotHaloHost  %mass           () &
            &                                  +hotHalo      %massCold       () &
            &                                 )
       call hotHaloHost%    angularMomentumSet(                                 &
            &                                   hotHaloHost  %angularMomentum() &
            &                                  +hotHalo      %massCold       () &
            &                                  *spinHost     %angularMomentum() &
            &                                  /basicHost    %mass           () &
            &                                 )
       call hotHalo    %           massColdSet(                                 &
            &                                   0.0d0                           &
            &                                 )
       call hotHalo    %angularMomentumColdSet(                                 &
            &                                   0.0d0                           &
            &                                 )
       call hotHaloHost%         abundancesSet(                                 &
            &                                   hotHaloHost  %abundances     () &
            &                                  +hotHalo      %abundancesCold () &
            &                                 )
       call hotHalo    %     abundancesColdSet(                                 &
            &                                  zeroAbundances                   &
            &                                 )
    end select
    return
  end subroutine satelliteMerger

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply
    update the cold mode mass of {\normalfont \ttfamily node} to account for any cold mode gas already in the
    parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, nodeComponentHotHaloColdMode, treeNode
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeParent
    class(nodeComponentHotHalo), pointer       :: hotHaloParent, hotHalo
    !$GLC attributes unused :: self

    hotHalo       => node      %hotHalo(autoCreate=.true.)
    nodeParent    => node      %parent
    hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
    ! If the parent node has a hot halo component, then add its cold mode to that of this node,
    ! and perform other changes needed prior to promotion.
    select type (hotHaloParent)
    class is (nodeComponentHotHaloColdMode)
       call hotHalo%           massColdSet(                                      &
            &                                hotHalo      %massCold           () &
            &                               +hotHaloParent%massCold           () &
            &                              )
       call hotHalo%angularMomentumColdSet(                                      &
            &                                hotHalo      %angularMomentumCold() &
            &                               +hotHaloParent%angularMomentumCold() &
            &                              )
       call hotHalo%     abundancesColdSet(                                      &
            &                                hotHalo      %abundancesCold     () &
            &                               +hotHaloParent%abundancesCold     () &
            &                              )
    end select
    return
  end subroutine nodePromotion

  !![
  <haloFormationTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_Formation</unitName>
  </haloFormationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Formation(node)
    !!{
    Updates the hot halo gas distribution at a formation event, if requested.
    !!}
    use :: Abundances_Structure                 , only : abundances                     , zeroAbundances
    use :: Galacticus_Nodes                     , only : nodeComponentHotHalo           , nodeComponentHotHaloColdMode, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : hotHaloOutflowReturnOnFormation
    implicit none
    type (treeNode            ), intent(inout) :: node
    class(nodeComponentHotHalo), pointer       :: hotHalo

    ! Return immediately if return of outflowed gas on formation events is not requested.
    if (.not.hotHaloOutflowReturnOnFormation) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Transfer mass, angular momentum and abundances.
       call hotHalo%                    massSet(                                    &
            &                                    hotHalo%         mass           () &
            &                                   +hotHalo%outflowedMass           () &
            &                                  )
       call hotHalo%         angularMomentumSet(                                    &
            &                                    hotHalo%         angularMomentum() &
            &                                   +hotHalo%outflowedAngularMomentum() &
            &                                  )
       call hotHalo%              abundancesSet(                                    &
            &                                    hotHalo%         abundances     () &
            &                                   +hotHalo%outflowedAbundances     () &
            &                                  )
       call hotHalo%           outflowedMassSet(                                    &
            &                                    0.0d0                              &
            &                                  )
       call hotHalo%outflowedAngularMomentumSet(                                    &
            &                                    0.0d0                              &
            &                                  )
       call hotHalo%     outflowedAbundancesSet(                                    &
            &                                    zeroAbundances                     &
            &                                  )
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Formation

  !![
  <galacticusStateStoreTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_State_Store</unitName>
  </galacticusStateStoreTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentHotHalo -> coldMode',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="accretionHalo_ coldModeInfallRate_ cosmologyParameters_ galacticStructure_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_State_Store

  !![
  <galacticusStateRetrieveTask>
   <unitName>Node_Component_Hot_Halo_Cold_Mode_State_Restore</unitName>
  </galacticusStateRetrieveTask>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentHotHalo -> coldMode',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="accretionHalo_ coldModeInfallRate_ cosmologyParameters_ galacticStructure_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_State_Restore

end module Node_Component_Hot_Halo_Cold_Mode

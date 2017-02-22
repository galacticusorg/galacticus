!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements an extension to the standard hot halo node component which
!% supports a cold mode reservoir.

module Node_Component_Hot_Halo_Cold_Mode
  !% Implements an extension to the standard hot halo node component which supports a cold mode
  !% reservoir.
  use Galacticus_Nodes
  use Radiation_Structure
  use ISO_Varying_String
  implicit none
  private
  public :: Node_Component_Hot_Halo_Cold_Mode_Initialize       , Node_Component_Hot_Halo_Cold_Mode_Rate_Compute     , &
       &    Node_Component_Hot_Halo_Cold_Mode_Scale_Set        , Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize  , &
       &    Node_Component_Hot_Halo_Cold_Mode_Node_Merger      , Node_Component_Hot_Halo_Cold_Mode_Satellite_Merging, &
       &    Node_Component_Hot_Halo_Cold_Mode_Promote          , Node_Component_Hot_Halo_Cold_Mode_Formation        , &
       &    Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize

  !# <component>
  !#  <class>hotHalo</class>
  !#  <name>coldMode</name>
  !#  <extends>
  !#    <class>hotHalo</class>
  !#    <name>standard</name>
  !#  </extends>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>massCold</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of cold-mode gas in the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesCold</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the cold-mode of the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>angularMomentumCold</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of cold-mode gas in the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>massTotal</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <getFunction>Node_Component_Hot_Halo_Cold_Mode_Mass_Total</getFunction>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.hot_halo.cold_mode.bound_functions.inc</functions>
  !# </component>

  ! Options controlling the behavior of the cold mode gas.
  logical              :: hotHaloOutflowToColdMode   
  ! Internal count of abundances.
  integer              :: abundancesCount
  ! Name of the mass distribution to use for cold mode.
  type(varying_string) :: coldModeMassDistributionName
  ! Record of whether this module has been initialized.
  logical              :: moduleInitialized       =.false.

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Initialize()
    !% Initializes the standard hot halo component module.
    use Input_Parameters
    use Abundances_Structure
    implicit none
    type(nodeComponentHotHaloColdMode) :: hotHaloComponent

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Hot_Halo_Cold_Mode_Initialize)
    if (defaultHotHaloComponent%coldModeIsActive().and..not.moduleInitialized) then
       ! Get numbers of abundance properties.
       abundancesCount=Abundances_Property_Count()
       ! Determine whether outflows go to the cold mode.
       !@ <inputParameter>
       !@   <name>hotHaloOutflowToColdMode</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether or not outflows from galaxies are returned to the cold or hot modes in the hot halo.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloOutflowToColdMode',hotHaloOutflowToColdMode,defaultValue=.false.)
       ! Bind the outflow return function if outflow returns to the cold mode. (If it does not, do
       ! not bind any function and let the parent class handle this behavior.)
       if (hotHaloOutflowToColdMode) call hotHaloComponent%outflowReturnFunction(Node_Component_Hot_Halo_Cold_Mode_Outflow_Return)
       ! Create the cold mode mass distribution.
       !@ <inputParameter>
       !@   <name>coldModeMassDistribution</name>
       !@   <defaultValue>betaProfile</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The type of mass distribution to use for the cold mode component.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('coldModeMassDistribution',coldModeMassDistributionName,defaultValue="betaProfile")
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Hot_Halo_Cold_Mode_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Initialize

  !# <mergerTreeEvolveThreadInitialize>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize</unitName>
  !# </mergerTreeEvolveThreadInitialize>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize()
    !% Initializes the tree node hot halo methods module.
    use Node_Component_Hot_Halo_Cold_Mode_Structure_Tasks
    use Mass_Distributions
    implicit none

    ! Check if this implementation is selected. Define the radiation component to include both the CMB and the intergalactic background if it is.
    if (defaultHotHaloComponent%coldModeIsActive()) then
       call Node_Component_Hot_Halo_Cold_Mode_Initialize()
       coldModeMassDistribution => Mass_Distribution_Create(char(coldModeMassDistributionName))
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize

  subroutine Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes(node,massRate,interrupt,interruptProcedure)
    !% Push mass through the cooling pipes (along with appropriate amounts of metals and angular momentum) at the given rate.
    use Galacticus_Error
    use Cooling_Infall_Radii
    use Cooling_Specific_Angular_Momenta
    use Abundances_Structure
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode                    ), intent(inout)          , pointer :: node
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
          if (hotHalo%hotHaloCoolingMassRateIsAttached()) then
             call hotHalo%hotHaloCoolingMassRate(+massRate,interrupt,interruptProcedure)
             if (interrupt) return
          end if
          ! Find the node to use for cooling calculations.
          select case (hotHaloCoolingFromNode)
          case (currentNode  )
             nodeCoolingFrom => node
          case (formationNode)
             nodeCoolingFrom => node%formationNode
          case default
             nodeCoolingFrom => null()
             call Galacticus_Error_Report('Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes','unknown cooling node')
          end select
          ! Compute the infall rate of angular momentum.
          angularMomentumCoolingRate=massRate*hotHalo%angularMomentumCold()/hotHalo%massCold()
          call hotHalo%angularMomentumColdRate(-angularMomentumCoolingRate)
          ! Pipe the cooling rate to which ever component claimed it.
          if (hotHalo%hotHaloCoolingAngularMomentumRateIsAttached()) then
             call hotHalo%hotHaloCoolingAngularMomentumRate(sign(+angularMomentumCoolingRate*(1.0d0-hotHaloAngularMomentumLossFraction),massRate),interrupt,interruptProcedure)
             if (interrupt) return
          end if
          ! Get the rate of change of abundances.
          hotHaloCoolingFrom => nodeCoolingFrom   %hotHalo       ()
          abundancesCoolingRate=hotHaloCoolingFrom%abundancesCold()
          abundancesCoolingRate=massRate*abundancesCoolingRate/hotHaloCoolingFrom%massCold()
          call hotHalo%abundancesColdRate(-abundancesCoolingRate)
          ! Pipe the cooling rate to which ever component claimed it.
          if (hotHalo%hotHaloCoolingAbundancesRateIsAttached()) then
             call hotHalo%hotHaloCoolingAbundancesRate(+abundancesCoolingRate,interrupt,interruptProcedure)
             if (interrupt) return
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Rate_Compute(node,odeConverged,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Abundances_Structure
    use Accretion_Halos
    use Dark_Matter_Halo_Spins
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Hot_Halo_Ram_Pressure_Stripping
    use Node_Component_Hot_Halo_Standard_Data
    use Cooling_Cold_Mode_Infall_Rates
    use Galactic_Structure_Options
    use Galactic_Structure_Densities
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: node
    logical                                       , intent(in   )          :: odeConverged
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo        )               , pointer :: hotHalo
    class           (nodeComponentBasic          )               , pointer :: basic
    class           (darkMatterHaloScaleClass    )               , pointer :: darkMatterHaloScale_
    class           (accretionHaloClass          )               , pointer :: accretionHalo_
    double precision                                                       :: angularMomentumAccretionRate, densityAtOuterRadius , &
         &                                                                    massAccretionRate           , massLossRate         , &
         &                                                                    outerRadius                 , outerRadiusGrowthRate, &
         &                                                                    gasMass                     , infallRate
    !GCC$ attributes unused :: odeConverged
    
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that the standard hot halo implementation is active.
    if (defaultHotHaloComponent%coldModeIsActive()) then
       ! Get required objects.
       darkMatterHaloScale_ => darkMatterHaloScale()
       accretionHalo_       => accretionHalo      ()
       ! Find the rate of gas mass accretion onto the halo.
       massAccretionRate=accretionHalo_%accretionRate(node,accretionModeCold)
       ! Get the basic component.
       basic => node%basic()
       ! Apply accretion rates.
       if (massAccretionRate > 0.0d0 .or. hotHalo%massCold() > 0.0d0) call hotHalo%massColdRate(massAccretionRate,interrupt,interruptProcedure)
       ! Next compute the cold mode infall rate in this halo.
       infallRate=Cooling_Cold_Mode_Infall_Rate(node)
       ! Pipe the cooling rate to which ever component claimed it.
       call Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes(node,infallRate,interrupt,interruptProcedure)
       ! Get the rate at which abundances are accreted onto this halo.
       call hotHalo%abundancesColdRate(accretionHalo_%accretionRateMetals(node,accretionModeCold),interrupt,interruptProcedure)
       ! Next block of tasks occur only if the accretion rate is non-zero.
       if (massAccretionRate > 0.0d0) then
          ! Compute the rate of accretion of angular momentum.
          angularMomentumAccretionRate=Dark_Matter_Halo_Angular_Momentum_Growth_Rate(node)*(massAccretionRate &
               &/basic%accretionRate())
          if (hotHaloOutflowAngularMomentumAlwaysGrows) angularMomentumAccretionRate=abs(angularMomentumAccretionRate)
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
                  &   outerRadius           <=                                   darkMatterHaloScale_%virialRadius(node) &
                  &  .and.                                                                                               &
                  &   outerRadius           > outerRadiusOverVirialRadiusMinimum*darkMatterHaloScale_%virialRadius(node) &
                  & ) then
                ! The ram pressure stripping radius is within the outer radius. Cause the outer radius to shrink to the ram pressure
                ! stripping radius on the halo dynamical timescale.
                densityAtOuterRadius = Galactic_Structure_Density(node,[outerRadius,0.0d0,0.0d0],coordinateSystemSpherical,componentTypeColdHalo,massTypeGaseous,haloLoaded=.true.)
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
                call hotHalo% strippedAbundancesRate(-hotHalo%abundances         ()*massLossRate/gasMass,interrupt,interruptProcedure)
             end if
          end if
       end select
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Rate_Compute

  subroutine Node_Component_Hot_Halo_Cold_Mode_Outflow_Return(self,interrupt,interruptProcedure)
    !% Return outflowed gas to the cold mode reservoir.
    use Galacticus_Error
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Atomic
    use Numerical_Constants_Math
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Galactic_Structure_Densities
    use Galactic_Structure_Options
    use Node_Component_Hot_Halo_Standard_Data
    use Cosmology_Parameters
    implicit none
    class           (nodeComponentHotHaloStandard), intent(inout)          :: self
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    type            (treeNode                    ), pointer                :: selfNode
    class           (nodeComponentBasic          ), pointer                :: selfBasic
    class           (cosmologyParametersClass    ), pointer                :: cosmologyParameters_
    class           (darkMatterHaloScaleClass    ), pointer                :: darkMatterHaloScale_
    double precision                                                       :: outflowedMass            , massReturnRate, &
         &                                                                    angularMomentumReturnRate, radiusVirial  , &
         &                                                                    densityAtOuterRadius     , densityMinimum, &
         &                                                                    outerRadius
    type            (abundances                  ), save                   :: abundancesReturnRate
    !$omp threadprivate(abundancesReturnRate)

    select type (self)
    class is (nodeComponentHotHaloColdMode)
       ! Get required objects.
       darkMatterHaloScale_ => darkMatterHaloScale()
       ! Get the hosting node.
       selfNode => self%hostNode
       ! Next tasks occur only for systems in which outflowed gas is being recycled.
       massReturnRate=0.0d0
       if (.not.starveSatellites.or..not.selfNode%isSatellite()) then
          outflowedMass            =self%outflowedMass()
          massReturnRate           =hotHaloOutflowReturnRate*outflowedMass                  /darkMatterHaloScale_%dynamicalTimescale(selfNode)
          angularMomentumReturnRate=hotHaloOutflowReturnRate*self%outflowedAngularMomentum()/darkMatterHaloScale_%dynamicalTimescale(selfNode)
          abundancesReturnRate     =hotHaloOutflowReturnRate*self%outflowedAbundances     ()/darkMatterHaloScale_%dynamicalTimescale(selfNode)
          call self%           outflowedMassRate(-           massReturnRate,interrupt,interruptProcedure)
          call self%                massColdRate(+           massReturnRate,interrupt,interruptProcedure)
          call self%outflowedAngularMomentumRate(-angularMomentumReturnRate,interrupt,interruptProcedure)
          call self%     angularMomentumColdRate(+angularMomentumReturnRate,interrupt,interruptProcedure)
          call self%     outflowedAbundancesRate(-     abundancesReturnRate,interrupt,interruptProcedure)
          call self%          abundancesColdRate(+     abundancesReturnRate,interrupt,interruptProcedure)
       end if
       ! The outer radius must be increased as the halo fills up with gas.
       outerRadius =self%outerRadius()
       radiusVirial=darkMatterHaloScale_%virialRadius(selfNode)
       if (outerRadius < radiusVirial) then 
          densityAtOuterRadius=Galactic_Structure_Density(selfNode,[outerRadius,0.0d0,0.0d0],coordinateSystemSpherical,componentTypeColdHalo,massTypeGaseous,haloLoaded=.true.)
          ! If the outer radius and density are non-zero we can expand the outer radius at a rate determined by the current
          ! density profile.
          if (outerRadius > 0.0d0 .and. densityAtOuterRadius > 0.0d0) then
             ! Limit the density at the outer radius to one third of the mean virial density (for baryons, assuming a
             ! universal baryon fraction) to prevent arbitrarily rapid growth of the outer radius in halos containing almost
             ! no gas.
             cosmologyParameters_ => cosmologyParameters()
             selfBasic => selfNode%basic()
             densityMinimum=(cosmologyParameters_%omegaBaryon()/cosmologyParameters_%omegaMatter())*selfBasic%mass()/radiusVirial**3/4.0d0/Pi
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
             call self%outerRadiusRate(darkMatterHaloScale_%virialVelocity(selfNode)*kilo*gigaYear/megaParsec)
          end if
       end if
    class default
       call Galacticus_Error_Report('Node_Component_Hot_Halo_Cold_Mode_Outflow_Return','this function should not be called for non-coldMode class hot halo components')
    end select
    return   
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Outflow_Return

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo    )               , pointer :: hotHalo
    class           (nodeComponentBasic      )               , pointer :: basic
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                          , parameter              :: scaleMassRelative   =1.0d-3
    double precision                          , parameter              :: scaleRadiusRelative =1.0d+0
    double precision                                                   :: massVirial                 , radiusVirial, &
         &                                                                velocityVirial

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of the cold mode class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Get required objects.
       darkMatterHaloScale_ => darkMatterHaloScale()
       ! The the basic component.
       basic => node%basic()
       ! Get virial properties.
       massVirial    =basic%mass()
       radiusVirial  =darkMatterHaloScale_%virialRadius  (node)
       velocityVirial=darkMatterHaloScale_%virialVelocity(node)
       call    hotHalo%           massColdScale(               massVirial                            *scaleMassRelative  )
       call    hotHalo%     abundancesColdScale(unitAbundances*massVirial                            *scaleMassRelative  )
       call    hotHalo%angularMomentumColdScale(               massVirial*radiusVirial*velocityVirial*scaleMassRelative  )
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize</unitName>
  !#  <after>spin</after>
  !#  <after>darkMatterProfile</after>
  !#  <after>Node_Component_Hot_Halo_Standard_Tree_Initialize</after>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize(node)
    !% Initialize the contents of the hot halo component for any sub-resolution accretion (i.e. the gas that would have been
    !% accreted if the merger tree had infinite resolution).
    use Accretion_Halos
    use Dark_Matter_Halo_Spins
    use Abundances_Structure
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    class           (accretionHaloClass  )               , pointer :: accretionHalo_
    class           (nodeEvent           )               , pointer :: event
    double precision                                               :: angularMomentum   , coldModeMass

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
    ! Get required objects.
    accretionHalo_ => accretionHalo()
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Get the mass of cold mode gas accreted.
    coldModeMass=accretionHalo_%accretedMass(node,accretionModeCold)
    ! If non-zero, then create a hot halo component and add to it.
    if (coldModeMass > 0.0d0) then
       ! Ensure that it is of unspecified class.
       hotHalo => node%hotHalo(autoCreate=.true.)
       basic   => node%basic  (                 )
       call hotHalo%massColdSet(coldModeMass)
       ! Also add the appropriate angular momentum.
       angularMomentum=coldModeMass*Dark_Matter_Halo_Angular_Momentum(node)/basic%mass()
       call hotHalo%angularMomentumColdSet(angularMomentum)
       ! Add the appropriate abundances.
       call hotHalo%abundancesColdSet(accretionHalo_%accretedMassMetals(node,accretionModeCold))
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Node_Merger</unitName>
  !#  <before>Node_Component_Hot_Halo_Standard_Node_Merger</before>
  !# </nodeMergerTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Node_Merger(node)
    !% Starve {\normalfont \ttfamily node} by transferring its hot halo to its parent.
    use Accretion_Halos
    use Abundances_Structure
    use Dark_Matter_Halo_Spins
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Cosmology_Parameters
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode                ), intent(inout), pointer :: node
    type            (treeNode                )               , pointer :: nodeParent
    class           (nodeComponentHotHalo    )               , pointer :: hotHaloParent          , hotHalo
    class           (nodeComponentSpin       )               , pointer :: spinParent
    class           (nodeComponentBasic      )               , pointer :: basicParent            , basic
    class           (cosmologyParametersClass)               , pointer :: cosmologyParameters_
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    class           (accretionHaloClass      )               , pointer :: accretionHalo_
    double precision                                                   :: baryonicMassCurrent    , baryonicMassMaximum   , &
         &                                                                fractionRemove         , massAccretedCold      , &
         &                                                                massAccreted           , massUnaccreted        , &
         &                                                                angularMomentumAccreted, massReaccreted        , &
         &                                                                fractionAccreted
    type            (abundances              ), save                   :: massMetalsAccreted     , fractionMetalsAccreted, &
         &                                                                massMetalsReaccreted
    !$omp threadprivate(massMetalsAccreted,fractionMetalsAccreted,massMetalsReaccreted)

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of cold mode class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Get required objects.
       darkMatterHaloScale_ => darkMatterHaloScale()
       accretionHalo_       => accretionHalo      ()
       ! Find the parent node and its hot halo and spin components.
       nodeParent    => node  %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       spinParent    => nodeParent%spin   (                 )
       basicParent   => nodeParent%basic  (                 )
       basic         => node  %basic  (                 )
       ! Determine if starvation is to be applied and if .
       if (starveSatellites) then
          ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate to
          ! this hot halo (and will be moved to the parent at the end of the evolution timestep).
          call hotHaloParent%           massColdSet(                                                 &
               &                                     hotHaloParent%massCold             (          ) &
               &                                    +hotHalo      %massCold             (          ) &
               &                                   )
          call hotHaloParent%angularMomentumColdSet(                                                 &
               &                                     hotHaloParent%angularMomentumCold  (          ) &
               &                                    +hotHalo      %massCold             (          ) &
               &                                    *   spinParent%spin                 (          ) &
               &                                    *darkMatterHaloScale_%virialRadius  (nodeParent) &
               &                                    *darkMatterHaloScale_%virialVelocity(nodeParent) &
               &                                   )
          call hotHalo      %           massColdSet(                                                 &
               &                                     0.0d0                                           &
               &                                   )
          call hotHalo      %angularMomentumColdSet(                                                 &
               &                                     0.0d0                                           &
               &                                   )
          call hotHaloParent%     abundancesColdSet(                                                 &
               &                                     hotHaloParent%abundancesCold     (            ) &
               &                                    +hotHalo      %abundancesCold     (            ) &
               &                                   )
          call hotHalo      %     abundancesColdSet(                                                 &
               &                                     zeroAbundances                                  &
               &                                   )
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
          if (basicParent%accretionRate() /= 0.0d0) then
             angularMomentumAccreted=+Dark_Matter_Halo_Angular_Momentum_Growth_Rate(nodeParent) &
                  &                  *massReaccreted                                            &
                  &                  /basicParent%accretionRate()
             call hotHaloParent%angularMomentumColdSet(hotHaloParent%angularMomentumCold()+angularMomentumAccreted)
          end if
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
          ! Check if the baryon fraction in the parent hot halo exceeds the universal value. If it does, mitigate this by moving
          ! some of the mass to the failed accretion reservoir.
          if (hotHaloNodeMergerLimitBaryonFraction) then
             cosmologyParameters_ => cosmologyParameters()
             baryonicMassMaximum  =  basicParent         %mass       () &
                  &                 *cosmologyParameters_%omegaBaryon() &	
                  &                 /cosmologyParameters_%omegaMatter()
             baryonicMassCurrent  =  Galactic_Structure_Enclosed_Mass(                                &
                  &                                                   nodeParent                    , &
                  &                                                   radiusLarge                   , &
                  &                                                   massType     =massTypeBaryonic, &
                  &                                                   componentType=componentTypeAll  &
                  &                                                  )
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

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Satellite_Merging</unitName>
  !# </satelliteMergerTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Satellite_Merging(node)
    !% Remove any cold mode gas associated with {\normalfont \ttfamily node} before it merges with its host halo.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type (treeNode                ), intent(inout), pointer :: node
    type (treeNode                )               , pointer :: nodeHost
    class(nodeComponentHotHalo    )               , pointer :: hotHaloHost         , hotHalo
    class(nodeComponentSpin       )               , pointer :: spinHost
    class(darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_

    ! Return immediately if satellites are starved, as in that case there is no hot halo to transfer.
    if (starveSatellites) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Get required objects.
       darkMatterHaloScale_ => darkMatterHaloScale()
       ! Find the node with which to merge.
       nodeHost    => node    %mergesWith(                 )
       hotHaloHost => nodeHost%hotHalo   (autoCreate=.true.)
       spinHost    => nodeHost%spin      (                 )
       ! Move the cold mode to the host.
       call hotHaloHost%               massSet(                                               &
            &                                   hotHaloHost%mass                   (        ) &
            &                                  +hotHalo    %massCold               (        ) &
            &                                 )
       call hotHaloHost%    angularMomentumSet(                                               &
            &                                   hotHaloHost%angularMomentum        (        ) &
            &                                  +hotHalo    %massCold               (        ) &
            &                                  *   spinHost%spin                   (        ) &
            &                                  *darkMatterHaloScale_%virialRadius  (nodeHost) &
            &                                  *darkMatterHaloScale_%virialVelocity(nodeHost) &
            &                                 )
       call hotHalo    %           massColdSet(                                               &
            &                                   0.0d0                                         &
            &                                 )
       call hotHalo    %angularMomentumColdSet(                                               &
            &                                   0.0d0                                         &
            &                                 )
       call hotHaloHost%         abundancesSet(                                               &
            &                                   hotHaloHost%abundances             (        ) &
            &                                  +hotHalo    %abundancesCold         (        ) &
            &                                 )
       call hotHalo    %     abundancesColdSet(                                               &
            &                                  zeroAbundances                                 &
            &                                 )
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Satellite_Merging

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Promote(node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply
    !% update the cold mode mass of {\normalfont \ttfamily node} to account for any cold mode gas already in the
    !% parent.
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    type (treeNode            )               , pointer :: nodeParent
    class(nodeComponentHotHalo)               , pointer :: hotHaloParent, hotHalo

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of specified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Get the parent node of this node and its hot halo component.
       nodeParent    => node  %parent
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
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Promote

  !# <haloFormationTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Formation</unitName>
  !# </haloFormationTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Formation(node)
    !% Updates the hot halo gas distribution at a formation event, if requested.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo

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

end module Node_Component_Hot_Halo_Cold_Mode

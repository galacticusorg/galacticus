!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  public :: Node_Component_Hot_Halo_Cold_Mode_Initialize       , Node_Component_Hot_Halo_Cold_Mode_Rate_Compute       , &
       &    Node_Component_Hot_Halo_Cold_Mode_Scale_Set        , Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize    , &
       &    Node_Component_Hot_Halo_Cold_Mode_Node_Merger      , Node_Component_Hot_Halo_Cold_Mode_Satellite_Merger   , &
       &    Node_Component_Hot_Halo_Cold_Mode_Promote          , Node_Component_Hot_Halo_Cold_Mode_Formation          , &
       &    Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize

  !# <component>
  !#  <class>hotHalo</class>
  !#  <name>coldMode</name>
  !#  <extends>
  !#    <class>hotHalo</class>
  !#    <name>standard</name>
  !#  </extends>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>massCold</name>
  !#     <type>real</type>
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
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of cold-mode gas in the hot halo."/>
  !#   </property>
  !#  </properties>
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
       !@   <type>real</type>
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

  subroutine Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes(thisNode,massRate,interrupt,interruptProcedure)
    !% Push mass through the cooling pipes (along with appropriate amounts of metals and angular momentum) at the given rate.
    use Cooling_Infall_Radii
    use Cooling_Specific_Angular_Momenta
    use Abundances_Structure
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode                    ), intent(inout)          , pointer :: thisNode
    double precision                              , intent(in   )                    :: massRate
    logical                                       , intent(inout), optional          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), optional, pointer :: interruptProcedure
    type            (treeNode                    )                         , pointer :: coolingFromNode
    class           (nodeComponentHotHalo        )                         , pointer :: coolingFromHotHalo        , thisHotHalo
    type            (abundances                  ), save                             :: abundancesCoolingRate
    !$omp threadprivate(abundancesCoolingRate)
    double precision                                                                 :: angularMomentumCoolingRate

    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    select type (thisHotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Ignore zero rates.
       if (massRate /= 0.0d0 .and. thisHotHalo%massCold() > 0.0d0 .and. thisHotHalo%angularMomentumCold() > 0.0d0) then
          ! Remove mass from the hot component.
          call thisHotHalo%massColdRate(-massRate)
          ! Pipe the mass rate to whichever component claimed it.
          if (thisHotHalo%hotHaloCoolingMassRateIsAttached()) then
             call thisHotHalo%hotHaloCoolingMassRate(+massRate,interrupt,interruptProcedure)
             if (interrupt) return
          end if
          ! Find the node to use for cooling calculations.
          select case (hotHaloCoolingFromNode)
          case (currentNode  )
             coolingFromNode => thisNode
          case (formationNode)
             coolingFromNode => thisNode%formationNode
          end select
          ! Compute the infall rate of angular momentum.
          angularMomentumCoolingRate=massRate*thisHotHalo%angularMomentumCold()/thisHotHalo%massCold()
          call thisHotHalo%angularMomentumColdRate(-angularMomentumCoolingRate)
          ! Pipe the cooling rate to which ever component claimed it.
          if (thisHotHalo%hotHaloCoolingAngularMomentumRateIsAttached()) then
             call thisHotHalo%hotHaloCoolingAngularMomentumRate(sign(+angularMomentumCoolingRate*(1.0d0-hotHaloAngularMomentumLossFraction),massRate),interrupt,interruptProcedure)
             if (interrupt) return
          end if
          ! Get the rate of change of abundances.
          coolingFromHotHalo => coolingFromNode   %hotHalo       ()
          abundancesCoolingRate=coolingFromHotHalo%abundancesCold()
          abundancesCoolingRate=massRate*abundancesCoolingRate/coolingFromHotHalo%massCold()
          call thisHotHalo%abundancesColdRate(-abundancesCoolingRate)
          ! Pipe the cooling rate to which ever component claimed it.
          if (thisHotHalo%hotHaloCoolingAbundancesRateIsAttached()) then
             call thisHotHalo%hotHaloCoolingAbundancesRate(+abundancesCoolingRate,interrupt,interruptProcedure)
             if (interrupt) return
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the hot halo node mass rate of change.
    use Abundances_Structure
    use Accretion_Halos
    use Accretion_Halos_Options
    use Dark_Matter_Halo_Spins
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Hot_Halo_Ram_Pressure_Stripping
    use Node_Component_Hot_Halo_Standard_Data
    use Cooling_Cold_Mode_Infall_Rates
    use Galactic_Structure_Options
    use Galactic_Structure_Densities
    implicit none
    type            (treeNode                    )      , intent(inout), pointer :: thisNode
    logical                                             , intent(inout)          :: interrupt
    procedure       (Interrupt_Procedure_Template)      , intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo        )                     , pointer :: thisHotHalo
    class           (nodeComponentBasic          )                     , pointer :: thisBasic
    type            (abundances                  ), save                         :: accretionRateAbundances
    !$omp threadprivate(accretionRateAbundances)
    double precision                                                             :: angularMomentumAccretionRate, densityAtOuterRadius , &
         &                                                                          massAccretionRate           , massLossRate         , &
         &                                                                          outerRadius                 , outerRadiusGrowthRate, &
         &                                                                          gasMass                     , infallRate

    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    ! Ensure that the standard hot halo implementation is active.
    if (defaultHotHaloComponent%coldModeIsActive()) then
       ! Find the rate of gas mass accretion onto the halo.
       massAccretionRate=Halo_Baryonic_Accretion_Rate(thisNode,accretionModeCold)
       ! Get the basic component.
       thisBasic => thisNode%basic()
       ! Apply accretion rates.
       if (massAccretionRate > 0.0d0 .or. thisHotHalo%massCold() > 0.0d0) call thisHotHalo%massColdRate(massAccretionRate,interrupt,interruptProcedure)
       ! Next compute the cold mode infall rate in this halo.
       infallRate=Cooling_Cold_Mode_Infall_Rate(thisNode)
       ! Pipe the cooling rate to which ever component claimed it.
       call Node_Component_Hot_Halo_Cold_Mode_Push_To_Cooling_Pipes(thisNode,infallRate,interrupt,interruptProcedure)
       ! Get the rate at which abundances are accreted onto this halo.
       call Halo_Baryonic_Accretion_Rate_Abundances(thisNode,accretionRateAbundances,accretionModeCold)
       call thisHotHalo%abundancesColdRate(accretionRateAbundances,interrupt,interruptProcedure)
       ! Next block of tasks occur only if the accretion rate is non-zero.
       if (massAccretionRate > 0.0d0) then
          ! Compute the rate of accretion of angular momentum.
          angularMomentumAccretionRate=Dark_Matter_Halo_Angular_Momentum_Growth_Rate(thisNode)*(massAccretionRate &
               &/thisBasic%accretionRate())
          if (hotHaloOutflowAngularMomentumAlwaysGrows) angularMomentumAccretionRate=abs(angularMomentumAccretionRate)
          call thisHotHalo%angularMomentumColdRate(angularMomentumAccretionRate,interrupt,interruptProcedure)
       end if
       select type (thisHotHalo)
       class is (nodeComponentHotHaloColdMode)
          ! Test whether this halo is a satellite or not.
          if (thisNode%isSatellite()) then
             ! For satellites, get the current ram pressure stripping radius for this hot halo.
             outerRadiusGrowthRate=thisHotHalo%outerRadiusGrowthRate()
             outerRadius          =thisHotHalo%outerRadius          ()
             gasMass              =thisHotHalo%massCold             ()
             if     (                                                                                                     &
                  &   outerRadiusGrowthRate /= 0.0d0                                                                      &
                  &  .and.                                                                                                &
                  &   gasMass               >  0.0d0                                                                      &
                  &  .and.                                                                                                &
                  &   outerRadius           <=                                   Dark_Matter_Halo_Virial_Radius(thisNode) &
                  &  .and.                                                                                                &
                  &   outerRadius           > outerRadiusOverVirialRadiusMinimum*Dark_Matter_Halo_Virial_Radius(thisNode) &
                  & ) then
                ! The ram pressure stripping radius is within the outer radius. Cause the outer radius to shrink to the ram pressure
                ! stripping radius on the halo dynamical timescale.
                densityAtOuterRadius = Galactic_Structure_Density(thisNode,[outerRadius,0.0d0,0.0d0],coordinateSystemSpherical,componentTypeColdHalo,massTypeGaseous,haloLoaded=.true.)
                ! Compute the mass loss rate.
                massLossRate=4.0d0*Pi*densityAtOuterRadius*outerRadius**2*outerRadiusGrowthRate
                ! Adjust the rates.
                ! Mass.
                call thisHotHalo%           massColdRate(                                  massLossRate        ,interrupt,interruptProcedure)
                ! Angular momentum.
                call thisHotHalo%angularMomentumColdRate(thisHotHalo%angularMomentumCold()*massLossRate/gasMass,interrupt,interruptProcedure)
                ! Metal abundances.
                call thisHotHalo%     abundancesColdRate(thisHotHalo%abundancesCold     ()*massLossRate/gasMass,interrupt,interruptProcedure)
                ! Mass.
                call thisHotHalo%       strippedMassRate(                                  massLossRate        ,interrupt,interruptProcedure)
                ! Metal abundances.
                call thisHotHalo% strippedAbundancesRate(thisHotHalo%abundances         ()*massLossRate/gasMass,interrupt,interruptProcedure)
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
    procedure       (Interrupt_Procedure_Template), intent(inout), pointer :: interruptProcedure
    type            (treeNode                    ), pointer                :: selfNode
    class           (nodeComponentBasic          ), pointer                :: selfBasic
    class           (cosmologyParametersClass    ), pointer                :: thisCosmologyParameters
    double precision                                                       :: outflowedMass            , massReturnRate         , &
         &                                                                    angularMomentumReturnRate, radiusVirial           , &
         &                                                                    densityAtOuterRadius     , densityMinimum         , &
         &                                                                    outerRadius
    type            (abundances                  ), save                   :: abundancesReturnRate
    !$omp threadprivate(abundancesReturnRate)

    select type (self)
    class is (nodeComponentHotHaloColdMode)
       ! Get the hosting node.
       selfNode => self%hostNode
       ! Next tasks occur only for systems in which outflowed gas is being recycled.
       if (.not.starveSatellites.or..not.selfNode%isSatellite()) then
          outflowedMass            =self%outflowedMass()
          massReturnRate           =hotHaloOutflowReturnRate*outflowedMass                  /Dark_Matter_Halo_Dynamical_Timescale(selfNode)
          angularMomentumReturnRate=hotHaloOutflowReturnRate*self%outflowedAngularMomentum()/Dark_Matter_Halo_Dynamical_Timescale(selfNode)
          abundancesReturnRate     =hotHaloOutflowReturnRate*self%outflowedAbundances     ()/Dark_Matter_Halo_Dynamical_Timescale(selfNode)
          call self%           outflowedMassRate(-           massReturnRate,interrupt,interruptProcedure)
          call self%                massColdRate(+           massReturnRate,interrupt,interruptProcedure)
          call self%outflowedAngularMomentumRate(-angularMomentumReturnRate,interrupt,interruptProcedure)
          call self%     angularMomentumColdRate(+angularMomentumReturnRate,interrupt,interruptProcedure)
          call self%     outflowedAbundancesRate(-     abundancesReturnRate,interrupt,interruptProcedure)
          call self%          abundancesColdRate(+     abundancesReturnRate,interrupt,interruptProcedure)
       end if
       ! The outer radius must be increased as the halo fills up with gas.
       outerRadius =self%outerRadius()
       radiusVirial=Dark_Matter_Halo_Virial_Radius(selfNode)
       if (outerRadius < radiusVirial) then 
          densityAtOuterRadius=Galactic_Structure_Density(selfNode,[outerRadius,0.0d0,0.0d0],coordinateSystemSpherical,componentTypeColdHalo,massTypeGaseous,haloLoaded=.true.)
          ! If the outer radius and density are non-zero we can expand the outer radius at a rate determined by the current
          ! density profile.
          if (outerRadius > 0.0d0 .and. densityAtOuterRadius > 0.0d0) then
             ! Limit the density at the outer radius to one third of the mean virial density (for baryons, assuming a
             ! universal baryon fraction) to prevent arbitrarily rapid growth of the outer radius in halos containing almost
             ! no gas.
             thisCosmologyParameters => cosmologyParameters()
             selfBasic => selfNode%basic()
             densityMinimum=(thisCosmologyParameters%omegaBaryon()/thisCosmologyParameters%omegaMatter())*selfBasic%mass()/radiusVirial**3/4.0d0/Pi
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
             call self%outerRadiusRate(Dark_Matter_Halo_Virial_Velocity(selfNode)*kilo*gigaYear/megaParsec)
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
  subroutine Node_Component_Hot_Halo_Cold_Mode_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHalo
    class           (nodeComponentBasic  )               , pointer :: thisBasic
    double precision                      , parameter              :: scaleMassRelative   =1.0d-3
    double precision                      , parameter              :: scaleRadiusRelative =1.0d+0
    double precision                                               :: massVirial                 , radiusVirial, &
         &                                                            velocityVirial

    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    ! Ensure that it is of the cold mode class.
    select type (thisHotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! The the basic component.
       thisBasic => thisNode%basic()
       ! Get virial properties.
       massVirial    =thisBasic%mass()
       radiusVirial  =Dark_Matter_Halo_Virial_Radius  (thisNode)
       velocityVirial=Dark_Matter_Halo_Virial_Velocity(thisNode)
       call    thisHotHalo%           massColdScale(               massVirial                            *scaleMassRelative  )
       call    thisHotHalo%     abundancesColdScale(unitAbundances*massVirial                            *scaleMassRelative  )
       call    thisHotHalo%angularMomentumColdScale(               massVirial*radiusVirial*velocityVirial*scaleMassRelative  )
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize</unitName>
  !#  <after>spin</after>
  !#  <after>darkMatterProfile</after>
  !#  <after>Node_Component_Hot_Halo_Standard_Tree_Initialize</after>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize(thisNode)
    !% Initialize the contents of the hot halo component for any sub-resolution accretion (i.e. the gas that would have been
    !% accreted if the merger tree had infinite resolution).
    use Accretion_Halos
    use Accretion_Halos_Options
    use Dark_Matter_Halo_Spins
    use Abundances_Structure
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHalo
    class           (nodeComponentBasic  )               , pointer :: thisBasic
    type            (abundances          ), save                   :: accretedAbundances
    !$omp threadprivate(accretedAbundances)
    double precision                                               :: angularMomentum   , coldModeMass

    ! If the node has a child or the standard hot halo is not active, then return immediately.
    if (associated(thisNode%firstChild).or..not.defaultHotHaloComponent%coldModeIsActive()) return
    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    ! Get the mass of cold mode gas accreted.
    coldModeMass=Halo_Baryonic_Accreted_Mass(thisNode,accretionModeCold)
    ! If non-zero, then create a hot halo component and add to it.
    if (coldModeMass > 0.0d0) then
       ! Ensure that it is of unspecified class.
       thisHotHalo => thisNode%hotHalo(autoCreate=.true.)
       thisBasic   => thisNode%basic  (                 )
       call thisHotHalo%massColdSet(coldModeMass)
       ! Also add the appropriate angular momentum.
       angularMomentum=coldModeMass*Dark_Matter_Halo_Angular_Momentum(thisNode)/thisBasic%mass()
       call thisHotHalo%angularMomentumColdSet(angularMomentum)
       ! Add the appropriate abundances.
       call Halo_Baryonic_Accreted_Abundances(thisNode,accretedAbundances,accretionModeCold)
       call thisHotHalo%abundancesColdSet(accretedAbundances)
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Node_Merger</unitName>
  !#  <before>Node_Component_Hot_Halo_Standard_Node_Merger</before>
  !# </nodeMergerTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Node_Merger(thisNode)
    !% Starve {\tt thisNode} by transferring its hot halo to its parent.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Cosmology_Parameters
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    type            (treeNode            )               , pointer :: parentNode
    class           (nodeComponentHotHalo)               , pointer :: parentHotHalo      , thisHotHalo
    class           (nodeComponentSpin   )               , pointer :: parentSpin
    class           (nodeComponentBasic  )               , pointer :: parentBasic
    class           (cosmologyParametersClass)           , pointer :: thisCosmologyParameters
    double precision                                               :: baryonicMassCurrent, baryonicMassMaximum, fractionRemove

    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    ! Ensure that it is of cold mode class.
    select type (thisHotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Find the parent node and its hot halo and spin components.
       parentNode    => thisNode  %parent
       parentHotHalo => parentNode%hotHalo(autoCreate=.true.)
       parentSpin    => parentNode%spin   (                 )
       ! Determine if starvation is to be applied and if .
       if (starveSatellites) then
          ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate to
          ! this hot halo (and will be moved to the parent at the end of the evolution timestep).
          call parentHotHalo%           massColdSet(                                                &
               &                                      parentHotHalo%massCold           (          ) &
               &                                     +  thisHotHalo%massCold           (          ) &
               &                                    )
          call parentHotHalo%angularMomentumColdSet(                                                &
               &                                      parentHotHalo%angularMomentumCold(          ) &
               &                                     +  thisHotHalo%massCold           (          ) &
               &                                     *   parentSpin%spin               (          ) &
               &                                     *Dark_Matter_Halo_Virial_Radius   (parentNode) &
               &                                     *Dark_Matter_Halo_Virial_Velocity (parentNode) &
               &                                    )
          call   thisHotHalo%           massColdSet(                                                &
               &                                      0.0d0                                         &
               &                                    )
          call   thisHotHalo%angularMomentumColdSet(                                                &
               &                                      0.0d0                                         &
               &                                    )
          call parentHotHalo%     abundancesColdSet(                                                &
               &                                      parentHotHalo%abundancesCold     (          ) &
               &                                     +  thisHotHalo%abundancesCold     (          ) &
               &                                    )
          call   thisHotHalo%     abundancesColdSet(                                                &
               &                                      zeroAbundances                                &
               &                                    )
          ! Check if the baryon fraction in the parent hot halo exceeds the universal value. If it does, mitigate this by moving
          ! some of the mass to the failed accretion reservoir.
          if (hotHaloNodeMergerLimitBaryonFraction) then
             thisCosmologyParameters => cosmologyParameters()
             parentBasic             => parentNode%basic()
             baryonicMassMaximum     =  parentBasic%mass()&
                  &                    *thisCosmologyParameters%omegaBaryon() &	
                  &                    /thisCosmologyParameters%omegaMatter()
             baryonicMassCurrent     =  Galactic_Structure_Enclosed_Mass(                                &
                  &                                                      parentNode                    , &
                  &                                                      radiusLarge                   , &
                  &                                                      massType     =massTypeBaryonic, &
                  &                                                      componentType=componentTypeAll  &
                  &                                                     )
             if (baryonicMassCurrent > baryonicMassMaximum .and. parentHotHalo%mass()+parentHotHalo%massCold() > 0.0d0) then
                fractionRemove=min((baryonicMassCurrent-baryonicMassMaximum)/(parentHotHalo%mass()+parentHotHalo%massCold()),1.0d0)
                call parentHotHalo%     unaccretedMassSet(                                                       &
                     &                                     parentHotHalo%unaccretedMass ()                       &
                     &                                    +parentHotHalo%mass           ()*       fractionRemove &
                     &                                    +parentHotHalo%massCold       ()*       fractionRemove &
                     &                                   )
                call parentHotHalo%               massSet( parentHotHalo%mass               ()*(1.0d0-fractionRemove))
                call parentHotHalo%    angularMomentumSet( parentHotHalo%angularMomentum    ()*(1.0d0-fractionRemove))
                call parentHotHalo%         abundancesSet( parentHotHalo%abundances         ()*(1.0d0-fractionRemove))
                call parentHotHalo%           massColdSet( parentHotHalo%massCold           ()*(1.0d0-fractionRemove))
                call parentHotHalo%angularMomentumColdSet( parentHotHalo%angularMomentumCold()*(1.0d0-fractionRemove))
                call parentHotHalo%     abundancesColdSet( parentHotHalo%abundancesCold     ()*(1.0d0-fractionRemove))
             end if
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Node_Merger

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Satellite_Merger</unitName>
  !# </satelliteMergerTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Satellite_Merger(thisNode)
    !% Remove any cold mode gas associated with {\tt thisNode} before it merges with its host halo.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: hostNode
    class(nodeComponentHotHalo)               , pointer :: hostHotHalo, thisHotHalo
    class(nodeComponentSpin   )               , pointer :: hostSpin

    ! Return immediately if satellites are starved, as in that case there is no hot halo to transfer.
    if (starveSatellites) return
    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (thisHotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Find the node with which to merge.
       hostNode    => thisNode%mergesWith()
       hostHotHalo => hostNode%hotHalo   ()
       hostSpin    => hostNode%spin      ()
       ! Move the cold mode to the host.
       call hostHotHalo%               massSet(                                            &
            &                                   hostHotHalo%mass                (        ) &
            &                                  +thisHotHalo%massCold            (        ) &
            &                                 )
       call hostHotHalo%    angularMomentumSet(                                            &
            &                                   hostHotHalo%angularMomentum     (        ) &
            &                                  +thisHotHalo%massCold            (        ) &
            &                                  *   hostSpin%spin                (        ) &
            &                                  *Dark_Matter_Halo_Virial_Radius  (hostNode) &
            &                                  *Dark_Matter_Halo_Virial_Velocity(hostNode) &
            &                                 )
       call thisHotHalo%           massColdSet(                                            &
            &                                   0.0d0                                      &
            &                                 )
       call thisHotHalo%angularMomentumColdSet(                                            &
            &                                   0.0d0                                      &
            &                                 )
       call hostHotHalo%         abundancesSet(                                            &
            &                                   hostHotHalo%abundances          (        ) &
            &                                  +thisHotHalo%abundancesCold      (        ) &
            &                                 )
       call thisHotHalo%     abundancesColdSet(                                            &
            &                                  zeroAbundances                              &
            &                                 )
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Satellite_Merger

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply
    !% update the cold mode mass of {\tt thisNode} to account for any cold mode gas already in the
    !% parent.
    use Dark_Matter_Halo_Scales
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHalo, thisHotHalo

    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    ! Ensure that it is of specified class.
    select type (thisHotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Get the parent node of this node and its hot halo component.
       parentNode    => thisNode  %parent
       parentHotHalo => parentNode%hotHalo()
       ! If the parent node has a hot halo component, then add its cold mode to that of this node,
       ! and perform other changes needed prior to promotion.
       select type (parentHotHalo)
       class is (nodeComponentHotHaloColdMode)
          call thisHotHalo%           massColdSet(                                  &
               &                                  thisHotHalo%massCold           () &
               &                               +parentHotHalo%massCold           () &
               &                              )
          call thisHotHalo%angularMomentumColdSet(                                  &
               &                                  thisHotHalo%angularMomentumCold() &
               &                               +parentHotHalo%angularMomentumCold() &
               &                              )
          call thisHotHalo%     abundancesColdSet(                                  &
               &                                  thisHotHalo%abundancesCold     () &
               &                               +parentHotHalo%abundancesCold     () &
               &                              )
        end select
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Promote

  !# <haloFormationTask>
  !#  <unitName>Node_Component_Hot_Halo_Cold_Mode_Formation</unitName>
  !# </haloFormationTask>
  subroutine Node_Component_Hot_Halo_Cold_Mode_Formation(thisNode)
    !% Updates the hot halo gas distribution at a formation event, if requested.
    use Abundances_Structure
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Astronomical
    use Node_Component_Hot_Halo_Standard_Data
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHalo

    ! Return immediately if return of outflowed gas on formation events is not requested.
    if (.not.hotHaloOutflowReturnOnFormation) return
    ! Get the hot halo component.
    thisHotHalo => thisNode%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (thisHotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Transfer mass, angular momentum and abundances.
       call thisHotHalo%                    massSet(&
            &                                                 thisHotHalo%         mass           () &
            &                                                +thisHotHalo%outflowedMass           () &
            &)
       call thisHotHalo%         angularMomentumSet(                                                 &
            &                                                 thisHotHalo%         angularMomentum() &
            &                                                +thisHotHalo%outflowedAngularMomentum() &
            &                                               )
       call thisHotHalo%              abundancesSet(                                                 &
            &                                                 thisHotHalo%         abundances     () &
            &                                                +thisHotHalo%outflowedAbundances     () &
            &                                               )
       call thisHotHalo%           outflowedMassSet(                                                 &
            &                                                 0.0d0                                  &
            &                                               )
       call thisHotHalo%outflowedAngularMomentumSet(                                                 &
            &                                                 0.0d0                                  &
            &                                               )
       call thisHotHalo%     outflowedAbundancesSet(                                                 &
            &                                                 zeroAbundances                         &
            &                                               )
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Formation

end module Node_Component_Hot_Halo_Cold_Mode

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

!% Contains a module which implements the standard black hole node component.

module Node_Component_Black_Hole_Standard
  !% Implement black hole tree node methods.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Black_Hole_Standard_Rate_Compute     , Node_Component_Black_Hole_Standard_Scale_Set        , &
       &    Node_Component_Black_Hole_Standard_Satellite_Merging, Node_Component_Black_Hole_Standard_Output_Properties, &
       &    Node_Component_Black_Hole_Standard_Output_Names     , Node_Component_Black_Hole_Standard_Output_Count     , &
       &    Node_Component_Black_Hole_Standard_Output           , Node_Component_Black_Hole_Standard_Initialize

  !# <component>
  !#  <class>blackHole</class>
  !#  <name>standard</name>
  !#  <isDefault>yes</isDefault>
  !#  <output instances="first"/>
  !#  <methods>
  !#   <method>
  !#     <name>mass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <classDefault>defaultBlackHoleComponent%massSeed()</classDefault>
  !#     <output unitsInSI="massSolar" comment="Mass of the black hole."/>
  !#   </method>
  !#   <method>
  !#     <name>spin</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <getFunction>Node_Component_Black_Hole_Standard_Spin</getFunction>
  !#     <classDefault>self%spinSeed()</classDefault>
  !#     <output unitsInSI="0.0d0" comment="Spin of the black hole."/>
  !#   </method>
  !#   <method>
  !#     <name>radialPosition</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#   </method>
  !#   <method>
  !#     <name>tripleInteractionTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </method>
  !#   <method>
  !#     <name>massSeed</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <getFunction>Node_Component_Black_Hole_Standard_Seed_Mass</getFunction>
  !#   </method>
  !#   <method>
  !#     <name>spinSeed</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
  !#     <getFunction>Node_Component_Black_Hole_Standard_Seed_Spin</getFunction>
  !#   </method>
  !#  </methods>
  !#  <functions>objects.nodes.components.black_hole.standard.custom_methods.inc</functions>
  !# </component>
  
  ! Accretion model parameters.
  ! Enhancement factors for the accretion rate.
  double precision                            :: bondiHoyleAccretionEnhancementSpheroid,bondiHoyleAccretionEnhancementHotHalo
  ! Temperature of accreting gas.
  double precision                            :: bondiHoyleAccretionTemperatureSpheroid
  ! Control for hot mode only accretion.
  logical                                     :: bondiHoyleAccretionHotModeOnly

  ! Feedback parameters.
  double precision                            :: blackHoleWindEfficiency
  logical                                     :: blackHoleHeatsHotHalo

  ! Output options.
  logical                                     :: blackHoleOutputAccretion
  logical                                     :: blackHoleOutputData
  logical                                     :: blackHoleOutputMergers
  
  ! Option specifying whether the triple black hole interaction should be used.
  logical                                     :: tripleBlackHoleInteraction

  ! Index of black hole instance about to merge.
  integer                                     :: mergingInstance
  !$omp threadprivate(mergingInstance)

  ! Index of black hole involved in three-body interactions
  integer                                     :: binaryInstance,tripleInstance
  !$omp threadprivate(binaryInstance,tripleInstance)

  ! Record of whether this module has been initialized.
  logical                                     :: moduleInitialized   =.false.

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Black_Hole_Standard_Initialize()
    !% Initializes the standard black hole component module.
    use Input_Parameters
    implicit none
    type(nodeComponentBlackHoleStandard) :: blackHoleStandard

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Black_Hole_Standard_Initialize)
    if (.not.moduleInitialized) then
       ! Get accretion rate enhancement factors.
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionEnhancementSpheroid</name>
       !@   <defaultValue>420</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which the Bondi-Hoyle accretion rate of spheroid gas onto black holes in enhanced.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionEnhancementSpheroid",bondiHoyleAccretionEnhancementSpheroid,defaultValue&
            &=420.0d0)
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionEnhancementHotHalo</name>
       !@   <defaultValue>1.5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which the Bondi-Hoyle accretion rate of hot halo gas onto black holes in enhanced.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionEnhancementHotHalo",bondiHoyleAccretionEnhancementHotHalo,defaultValue=1.5d0)
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionHotModeOnly</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether accretion from the hot halo should only occur if the halo is in the hot accretion mode.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionHotModeOnly",bondiHoyleAccretionHotModeOnly,defaultValue=.true.)

       ! Get temperature of accreting gas.
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionTemperatureSpheroid</name>
       !@   <defaultValue>$10^2$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The assumed temperature (in Kelvin) of gas in the spheroid when computing Bondi-Hoyle accretion rates onto black holes.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionTemperatureSpheroid",bondiHoyleAccretionTemperatureSpheroid,defaultValue&
            &=1.0d2)

       ! Get temperature of accreting gas.
       !@ <inputParameter>
       !@   <name>blackHoleWindEfficiency</name>
       !@   <defaultValue>$2.4\times 10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The efficiency of the black hole-driven wind: $L_{\rm wind} = \epsilon_{\rm wind} \dot{M}_\bullet \clight^2$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleWindEfficiency",blackHoleWindEfficiency,defaultValue=2.4d-3)

       ! Options controlling AGN feedback.
       !@ <inputParameter>
       !@   <name>blackHoleHeatsHotHalo</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not the black hole launched jets should heat the hot halo.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleHeatsHotHalo",blackHoleHeatsHotHalo,defaultValue=.true.)

       ! Get options controlling output.
       !@ <inputParameter>
       !@   <name>blackHoleOutputAccretion</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not accretion rates and jet powers will be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleOutputAccretion",blackHoleOutputAccretion,defaultValue=.false.)

       ! Get options controlling output.
       !@ <inputParameter>
       !@   <name>blackHoleOutputData</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not properties for all black holes (rather than just the central black hole) will be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleOutputData",blackHoleOutputData,defaultValue=.false.)

       !@ <inputParameter>
       !@   <name>blackHoleOutputMergers</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not properties of black hole mergers will be output.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleOutputMergers",blackHoleOutputMergers,defaultValue=.false.)

       ! Get options controlling three body interactions.
       !@ <inputParameter>
       !@   <name>tripleBlackHoleInteraction</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not triple black hole interactions will be accounted for.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>blackHoles</group>
       !@ </inputParameter>
       call Get_Input_Parameter("tripleBlackHoleInteraction",tripleBlackHoleInteraction,defaultValue=.false.)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Black_Hole_Standard_Initialize)
    return
  end subroutine Node_Component_Black_Hole_Standard_Initialize
 
  !# <rateComputeTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Black_Hole_Standard_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the black hole node mass rate of change.
    use Accretion_Disks
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    use Black_Hole_Binary_Separations
    implicit none
    type     (treeNode              ), pointer, intent(inout) :: thisNode
    logical                          ,          intent(inout) :: interrupt
    procedure(                      ), pointer, intent(inout) :: interruptProcedure
    class    (nodeComponentBlackHole), pointer                :: thisBlackHoleComponent,centralBlackHoleComponent&
         &,binaryBlackHoleComponent
    class    (nodeComponentSpheroid ), pointer                :: thisSpheroidComponent
    class    (nodeComponentHotHalo  ), pointer                :: thisHotHaloComponent
    class    (nodeComponentBasic    ), pointer                :: thisBasicComponent
    double precision                 , parameter              :: windVelocity  =1.0d4 ! Velocity of disk wind.
    double precision                 , parameter              :: ismTemperature=1.0d4 ! Temperature of the ISM.
    double precision                 , parameter              :: criticalDensityNormalization=2.0d0*massHydrogenAtom*speedLight&
         &**2*megaParsec/3.0d0/Pi/boltzmannsConstant/gigaYear/ismTemperature/kilo/windVelocity
    integer                                                   :: iInstance,instanceCount
    double precision                                          :: restMassAccretionRate,massAccretionRate,radiativeEfficiency &
         &,energyInputRate ,spheroidDensityRadius2,spheroidGasMass,spheroidRadius,criticalDensityRadius2,windFraction &
         &,spheroidDensityOverCriticalDensity ,heatingRate,couplingEfficiency,jetEfficiency,accretionRateSpheroid &
         &,accretionRateHotHalo,binaryRadius,radialMigrationRate,radiusHardBinary
    logical                                                   :: binaryRadiusFound
    
    if (defaultBlackHoleComponent%standardIsActive()) then

       ! Get a count of the number of black holes associated with this node.
       instanceCount=thisNode%blackHoleCount()
       ! Get the central black hole.
       centralBlackHoleComponent => thisNode%blackHole(instance=1)
       ! Get the basic, spheroid, and hot halo components.
       thisBasicComponent    => thisNode%basic   ()
       thisSpheroidComponent => thisNode%spheroid()
       thisHotHaloComponent  => thisNode%hotHalo ()
       ! Iterate over instances.
       do iInstance=1,max(instanceCount,1)
          ! Get the black hole.
          thisBlackHoleComponent => thisNode%blackHole(instance=iInstance)
          ! Find the rate of rest mass accretion onto the black hole.
          call Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(thisBlackHoleComponent,accretionRateSpheroid&
               &,accretionRateHotHalo)
          restMassAccretionRate=accretionRateSpheroid+accretionRateHotHalo
          ! Finish if there is no accretion.
          if (restMassAccretionRate <= 0.0d0) cycle    
          ! Find the radiative efficiency of the accretion.
          radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisBlackHoleComponent,restMassAccretionRate)
          ! Find the jet efficiency.
          if (restMassAccretionRate > 0.0d0) then
             jetEfficiency=Accretion_Disk_Jet_Power(thisBlackHoleComponent,restMassAccretionRate)/restMassAccretionRate/(speedLight&
                  &/kilo)**2
          else
             jetEfficiency=0.0d0
          end if
          ! Find the rate of increase in mass of the black hole.
          massAccretionRate=restMassAccretionRate*(1.0d0-radiativeEfficiency-jetEfficiency)       
          ! If no black hole component currently exists and we have some accretion then interrupt and create a black hole.
          if (instanceCount == 0 .and. massAccretionRate /= 0.0d0) then
             interrupt=.true.
             interruptProcedure => Node_Component_Black_Hole_Standard_Create
             return
          end if
          ! Add the accretion to the black hole.
          call thisBlackHoleComponent%massRate(massAccretionRate)    
          ! Remove the accreted mass from the spheroid component.
          call thisSpheroidComponent%massGasSinkRate(-accretionRateSpheroid)
          ! Remove the accreted mass from the hot halo component.
          call thisHotHaloComponent %   massSinkSet (-accretionRateHotHalo )
          ! Set spin-up rate due to accretion.
          if (massAccretionRate > 0.0d0) call thisBlackHoleComponent%spinRate(Black_Hole_Spin_Up_Rate(thisBlackHoleComponent,massAccretionRate))
          ! Add heating to the hot halo component.
          if (blackHoleHeatsHotHalo) then
             ! Compute jet coupling efficiency based on whether halo is cooling quasistatically. Reduce this efficiency as the gas
             ! content in the halo drops below the cosmological mean.
             couplingEfficiency=Hot_Mode_Fraction(thisNode)*((Omega_Matter()/Omega_b())*thisHotHaloComponent%mass()/thisBasicComponent%mass())**2
             ! Get jet power.
             heatingRate=jetEfficiency*restMassAccretionRate*(speedLight/kilo)**2*couplingEfficiency
             ! Pipe this power to the hot halo.
             call thisHotHaloComponent%heatSourceRate(heatingRate,interrupt,interruptProcedure)
          end if
          ! Add energy to the spheroid component.
          if (blackHoleWindEfficiency > 0.0d0) then
             spheroidGasMass=thisSpheroidComponent%massGas()
             if (spheroidGasMass > 0.0d0) then
                spheroidRadius=thisSpheroidComponent%radius()
                if (spheroidRadius > 0.0d0) then
                   spheroidDensityRadius2=3.0d0*spheroidGasMass/4.0d0/Pi/spheroidRadius
                   criticalDensityRadius2=criticalDensityNormalization*blackHoleWindEfficiency*restMassAccretionRate
                   ! Construct an interpolating factor such that the energy input from the wind drops to zero below half of the
                   ! critical density.
                   spheroidDensityOverCriticalDensity=spheroidDensityRadius2/criticalDensityRadius2-0.5d0
                   if (spheroidDensityOverCriticalDensity <= 0.0d0) then
                      ! No energy input below half of critical density.
                      windFraction=0.0d0
                   else if (spheroidDensityOverCriticalDensity >= 1.0d0) then
                      ! Full energy input above 1.5 times critical density.
                      windFraction=1.0d0
                   else
                      ! Smooth polynomial interpolating function between these limits.
                      windFraction=3.0d0*spheroidDensityOverCriticalDensity**2-2.0d0*spheroidDensityOverCriticalDensity**3
                   end if
                   ! Compute the energy input and send it down the spheroid gas energy input pipe.
                   energyInputRate=windFraction*blackHoleWindEfficiency*restMassAccretionRate*(speedLight/kilo)**2
                   call thisSpheroidComponent%energyGasInputRate(energyInputRate)
                end if
             end if
          end if
          ! Do radial migration for non-central black holes.
          if (iInstance > 1) then
             ! Compute the hard binary radius.
             radiusHardBinary= (                                                 &
                  &              gravitationalConstantGalacticus                 &
                  &             *(                                               &
                  &                centralBlackHoleComponent%mass()              &
                  &               +   thisBlackHoleComponent%mass()              &
                  &              )                                               &
                  &            )                                                 &
                  &           /(                                                 &
                  &              4.0d0                                           &
                  &             *  Dark_Matter_Halo_Virial_Velocity(thisNode)**2 &
                  &            )
             ! Places a new black hole in the center of the galaxy in case there is no central one.
             if     (                                                                        &
                  &  centralBlackHoleComponent%mass          () == 0.0d0               .and. &
                  &     thisBlackHoleComponent%radialPosition() <= radiusHardBinary          &
                  & ) then
                mergingInstance=iInstance
                interrupt=.true.
                interruptProcedure => Node_Component_Black_Hole_Standard_Merge_Black_Holes
                return
             end if
             ! Check for a black hole that is about to merge.
             if (thisBlackHoleComponent%radialPosition() <= 0.0d0) then    
                ! Record which instance is merging, then trigger an interrupt.
                mergingInstance=iInstance
                interrupt=.true.
                interruptProcedure => Node_Component_Black_Hole_Standard_Merge_Black_Holes
                return
             end if
             ! Set the rate of radial migration.
             radialMigrationRate=Black_Hole_Binary_Separation_Growth_Rate(thisBlackHoleComponent)
             call thisBlackHoleComponent%radialPositionRate(radialMigrationRate)
          end if
       end do
       ! Loop over black holes, testing for triple black hole interactions. Find the three closest black holes then check if a
       ! three body interaction occurs using the radial condition derived in Hoffman and Loeb (2007).
       binaryRadiusFound=.false.
       if (tripleBlackHoleInteraction) then
          if (instanceCount >= 3 .and. centralBlackHoleComponent%mass() > 0.0d0) then
             do iInstance=2,instanceCount
                ! Get the black hole.
                thisBlackHoleComponent => thisNode%blackHole(instance=iInstance)
                if     (                                                                   &
                     &  (          thisBlackHoleComponent%radialPosition() <= binaryRadius &
                     &   .or. .not.binaryRadiusFound                                       &
                     &  )                                                                  &
                     &  .and. .not.thisBlackHoleComponent%mass          () <= 0.0d0        &
                     &  .and. .not.thisBlackHoleComponent%radialPosition() <= 0.0d0        &
                     & ) then
                   binaryRadius     =thisBlackHoleComponent%radialPosition()
                   binaryInstance   =iInstance
                   binaryRadiusFound=.true.
                end if
             end do
             if (binaryRadiusFound) then
                ! Get the binary black hole.
                binaryBlackHoleComponent => thisNode%blackHole(instance=binaryInstance)
                ! Compute the hard binary radius.
                radiusHardBinary= (                                               &
                     &              gravitationalConstantGalacticus               &
                     &             *(                                             &
                     &                centralBlackHoleComponent%mass()            &
                     &               + binaryBlackHoleComponent%mass()            &
                     &              )                                             &
                     &            )                                               &
                     &           /(                                               &
                     &              4.0d0                                         &
                     &             *Dark_Matter_Halo_Virial_Velocity(thisNode)**2 &
                     &            )
                ! Search for a third black hole.  
                do iInstance=2,instanceCount
                   ! Get the black hole.
                   thisBlackHoleComponent => thisNode%blackHole(instance=iInstance)
                   if     (      .not. iInstance                                      == binaryInstance   &
                        &  .and. .not. thisBlackHoleComponent%mass                 () <= 0.0d0            &
                        &  .and. .not. thisBlackHoleComponent%radialPosition       () <= 0.0d0            &
                        &  .and.       thisBlackHoleComponent%radialPosition       () <= radiusHardBinary &
                        &  .and.       thisBlackHoleComponent%tripleInteractionTime() == 0.0d0            &
                        & ) then
                      tripleInstance=iInstance
                      interrupt=.true.
                      interruptProcedure => Node_Component_Black_Hole_Standard_Triple_Interaction
                      return
                   end if
                end do
             end if
          end if
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Rate_Compute
  
  !# <scaleSetTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Black_Hole_Standard_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: scaleMassRelative=1.0d-4
    double precision, parameter              :: scaleSizeRelative=1.0d-4
    double precision, parameter              :: scaleSizeAbsolute=1.0d-6
    class(nodeComponentSpheroid ), pointer :: thisSpheroidComponent
    class(nodeComponentBlackHole), pointer :: thisBlackHoleComponent
    integer                                  :: instance
 
    ! Determine if the standard implementation is active and at least one black hole exists.
    if (defaultBlackHoleComponent%standardIsActive().and.thisNode%blackHoleCount() > 0) then
       ! Get the spheroid component.
       thisSpheroidComponent => thisNode%spheroid()
       ! Loop over instances.
       do instance=1,thisNode%blackHoleCount()
          ! Get the black hole.
          thisBlackHoleComponent => thisNode%blackHole(instance=instance)
          ! Set scale for mass.
          call thisBlackHoleComponent%massScale(                                                            &
               &                                max(                                                        &
               &                                    scaleMassRelative*thisSpheroidComponent %massStellar(), &
               &                                                      thisBlackHoleComponent%mass       ()  &
               &                                   )                                                        &
               &                               )
          ! Set scale for spin.
          call thisBlackHoleComponent%spinScale(1.0d0)
          
          ! Set scale for radius.
          call thisBlackHoleComponent%radialPositionScale(                                                                 &
               &                                          maxval(                                                          &
               &                                               [                                                           &
               &                                                scaleSizeAbsolute,                                         &
               &                                                scaleSizeRelative*thisSpheroidComponent %halfMassRadius(), &
               &                                                                  thisBlackHoleComponent%radialPosition()  &
               &                                               ]                                                           &
               &                                                )                                                          &
               &                                         )
       end do       
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Satellite_Merging</unitName>
  !# </satelliteMergerTask>
  subroutine Node_Component_Black_Hole_Standard_Satellite_Merging(thisNode)
    !% Merge any black hole associated with {\tt thisNode} before it merges with its host halo.
    use Black_Hole_Binary_Mergers
    use Black_Hole_Binary_Initial_Radii
    use Black_Hole_Binary_Recoil_Velocities
    use Galactic_Structure_Potentials
    use Galactic_Structure_Options
    implicit none
    type            (treeNode              ), pointer, intent(inout) :: thisNode
    type            (treeNode              ), pointer                :: hostNode
    class           (nodeComponentBlackHole), pointer                :: hostCentralBlackHoleComponent,thisBlackHoleComponent
    integer                                                          :: instance
    double precision                                                 :: radiusInitial,blackHoleMassNew,blackHoleSpinNew&
         &,recoilVelocity,massBlackHole1,massBlackHole2 ,spinBlackHole1,spinBlackHole2

    ! Check that the standard black hole implementation is active.
    if (defaultBlackHoleComponent%standardIsActive()) then       
       ! Find the node to merge with.
       hostNode => thisNode%mergesWith()       
       ! Find the initial radius of the satellite black hole in the remnant.
       radiusInitial=Black_Hole_Binary_Initial_Radius(thisNode,hostNode)       
       ! If the separation is non-positive, assume that the black holes merge instantaneously.
       if (radiusInitial <= 0.0d0) then
          ! Get the central black hole of the host galaxy.
          hostCentralBlackHoleComponent => hostNode%blackHole(instance=1)
          ! Loop over all black holes in the satellite galaxy.
          do instance=1,thisNode%blackHoleCount()
             ! Get the black hole.
             thisBlackHoleComponent => thisNode%blackHole(instance=instance)
             ! Compute the outcome of the merger,
             call Black_Hole_Binary_Merger(thisBlackHoleComponent       %mass(), &
                  &                        hostCentralBlackHoleComponent%mass(), &
                  &                        thisBlackHoleComponent       %spin(), &
                  &                        hostCentralBlackHoleComponent%spin(), &
                  &                        blackHoleMassNew                    , &
                  &                        blackHoleSpinNew                      &
                  &                       )     
             ! Merge the black holes instantaneously.
             ! Check which black hole is more massive in order to compute an appropriate recoil velocity
             if (hostCentralBlackHoleComponent%mass() >= thisBlackHoleComponent%mass()) then
                massBlackHole1=hostCentralBlackHoleComponent%mass()
                massBlackHole2=       thisBlackHoleComponent%mass()
                spinBlackHole1=hostCentralBlackHoleComponent%spin()
                spinBlackHole2=       thisBlackHoleComponent%spin()
             else
                massBlackHole2=hostCentralBlackHoleComponent%mass()
                massBlackHole1=       thisBlackHoleComponent%mass()
                spinBlackHole2=hostCentralBlackHoleComponent%spin()
                spinBlackHole1=       thisBlackHoleComponent%spin()
             end if
             ! Now calculate the recoil velocity of the binary black hole and check wether it escapes the galaxy. (Note that
             ! we subtract the black hole's own contribution to the potential here.)
             recoilVelocity=Black_Hole_Binary_Recoil_Velocity(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)
             if (recoilVelocity > 0.0d0) then
                if (0.5d0*recoilVelocity**2+Galactic_Structure_Potential(thisNode,0.0d0)-Galactic_Structure_Potential(thisNode&
                     &,0.0d0,componentType=componentTypeBlackHole) > 0.0d0) then
                   blackHoleMassNew=0.0d0
                   blackHoleSpinNew=0.0d0
                end if
             end if
             ! Move the black hole to the host.
             call Node_Component_Black_Hole_Standard_Output_Merger(thisNode,massBlackHole1,massBlackHole2)
             call hostCentralBlackHoleComponent%massSet(blackHoleMassNew)
             call hostCentralBlackHoleComponent%spinSet(blackHoleSpinNew)
             ! Reset the satellite black hole to zero mass.
             call thisBlackHoleComponent%massSet(0.0d0)
             call thisBlackHoleComponent%spinSet(0.0d0)
          end do
       else
          ! Adjust the radii of the black holes in the satellite galaxy.
          do instance=thisNode%blackHoleCount(),1,-1
             thisBlackHoleComponent => thisNode%blackHole(instance=instance)
             call thisBlackHoleComponent%radialPositionSet(radiusInitial)
             ! Declares them as not having interacted in a triple black hole interaction.
             call thisBlackHoleComponent%tripleInteractionTimeSet(0.0d0)
             ! Remove this black hole if it has no mass.
             if (thisBlackHoleComponent%mass() <= 0.0d0) call thisNode%blackHoleRemove(instance)
          end do
          ! Move black holes from the satellite to the host.
          call thisNode%blackHoleMove(hostNode)
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Satellite_Merging
  
  subroutine Node_Component_Black_Hole_Standard_Merge_Black_Holes(thisNode)
    !% Merge two black holes.
    use Black_Hole_Binary_Recoil_Velocities
    use Black_Hole_Binary_Mergers
    use Galactic_Structure_Options
    use Galactic_Structure_Potentials
    implicit none
    type (treeNode              ), pointer, intent(inout) :: thisNode
    class(nodeComponentBlackHole), pointer                :: thisBlackHoleComponent1,thisBlackHoleComponent2
    double precision                                      :: blackHoleMassNew,blackHoleSpinNew,recoilVelocity,massBlackHole1&
         &,massBlackHole2,spinBlackHole1,spinBlackHole2

    ! Get the black holes.
    thisBlackHoleComponent1 => thisNode%blackHole(instance=              1)
    thisBlackHoleComponent2 => thisNode%blackHole(instance=mergingInstance)
    ! Process the merger to get the mass and spin of the merged black hole.
    call Black_Hole_Binary_Merger(thisBlackHoleComponent2%mass(), &
         &                        thisBlackHoleComponent1%mass(), &
         &                        thisBlackHoleComponent2%spin(), &
         &                        thisBlackHoleComponent1%spin(), &
         &                        blackHoleMassNew              , &
         &                        blackHoleSpinNew                &
         &                       )    
    ! Check which black hole is more massive in order to compute an appropriate recoil velocity.
    if (thisBlackHoleComponent1%mass() >= thisBlackHoleComponent2%mass()) then
       massBlackHole1=thisBlackHoleComponent1%mass()
       massBlackHole2=thisBlackHoleComponent2%mass()
       spinBlackHole1=thisBlackHoleComponent1%spin()
       spinBlackHole2=thisBlackHoleComponent2%spin()
    else
       massBlackHole2=thisBlackHoleComponent1%mass()
       massBlackHole1=thisBlackHoleComponent2%mass()
       spinBlackHole2=thisBlackHoleComponent1%spin()
       spinBlackHole1=thisBlackHoleComponent2%spin()
    end if    
    ! Calculate the recoil velocity of the binary black hole and check wether it escapes the galaxy
    recoilVelocity=Black_Hole_Binary_Recoil_Velocity(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)    
    ! Compare the recoil velocity to the potential and determine wether the binary is ejected or stays in the galaxy.           
    if (0.5d0*recoilVelocity**2+Galactic_Structure_Potential(thisNode,0.0d0)-Galactic_Structure_Potential(thisNode&
         &,0.0d0,componentType=componentTypeBlackHole) > 0.0d0) then
       blackHoleMassNew=0.0d0
       blackHoleSpinNew=0.0d0
    end if
    ! Set the mass and spin of the central black hole.
    call Node_Component_Black_Hole_Standard_Output_Merger(thisNode,massBlackHole1,massBlackHole2)
    call thisBlackHoleComponent1%massSet(blackHoleMassNew)
    call thisBlackHoleComponent2%spinSet(blackHoleSpinNew)    
    ! Remove the merging black hole from the list.
    call thisNode%blackHoleRemove(mergingInstance)
    return
  end subroutine Node_Component_Black_Hole_Standard_Merge_Black_Holes

  subroutine Node_Component_Black_Hole_Standard_Triple_Interaction(thisNode)
    !% Handles triple black holes interactions, using conditions similar to those of \cite{volonteri_assembly_2003}.
    use Black_Hole_Binary_Mergers
    use Galacticus_Nodes
    use Memory_Management
    use Numerical_Constants_Physical
    use Galactic_Structure_Options
    use Galactic_Structure_Potentials
    use Merger_Tree_Active
    implicit none
    type            (treeNode              ), pointer, intent(inout) :: thisNode
    class           (nodeComponentBasic    ), pointer                :: thisBasicComponent
    class           (nodeComponentBlackHole), pointer                :: centralBlackHoleComponent,binaryBlackHoleComponent&
         &,tripleBlackHoleComponent,ejectedBlackHoleComponent,newBinaryBlackHoleComponent
    integer                                                          :: ejectedInstance,newBinaryInstance
    double precision                                                 :: massEjected,massBinary,massRatioIntruder,velocityEjected &
         &,kineticEnergyChange,bindingEnergy,newRadius,velocityBinary
    logical                                                          :: removeBinary,removeEjected

    ! Get the basic component.
    thisBasicComponent        => thisNode%basic    (                       )
    ! Get the black holes.
    centralBlackHoleComponent => thisNode%blackHole(instance=             1)
    binaryBlackHoleComponent  => thisNode%blackHole(instance=binaryInstance)
    tripleBlackHoleComponent  => thisNode%blackHole(instance=tripleInstance)
    ! We have to distinguish two cases, where a different black hole is ejected, the one with the lowest mass.
    massRatioIntruder=   tripleBlackHoleComponent %mass() &
         &            /( centralBlackHoleComponent%mass() &
         &              +binaryBlackHoleComponent %mass() &
         &             )
    ! Set the triple interaction time for the triple black hole.
    call tripleBlackHoleComponent%tripleInteractionTimeSet(thisBasicComponent%time())
    ! Branch on intruder mass ratio.
    if (massRatioIntruder <= 2.0d0) then
       if (tripleBlackHoleComponent%mass() <= binaryBlackHoleComponent%mass()) then
          newRadius           = binaryBlackHoleComponent%radialPosition()/(1.0d0+0.4d0*massRatioIntruder)
          call binaryBlackHoleComponent%radialPositionSet(newRadius)
          bindingEnergy       = gravitationalConstantGalacticus*( tripleBlackHoleComponent %mass          () &
               &                                                 *centralBlackHoleComponent%mass          () &
               &                                                )                                            &
               &               /                                  tripleBlackHoleComponent %radialPosition()
          kineticEnergyChange=0.4d0*massRatioIntruder*bindingEnergy
          ejectedInstance    =tripleInstance
          newBinaryInstance  =binaryInstance
          ejectedBlackHoleComponent   => tripleBlackHoleComponent
          newBinaryBlackHoleComponent => binaryBlackHoleComponent
       else  
          newRadius          = tripleBlackHoleComponent%radialPosition()/(1.0d0+0.4d0*massRatioIntruder )
          call tripleBlackHoleComponent%radialPositionSet(newRadius)
          bindingEnergy      = gravitationalConstantGalacticus*( binaryBlackHoleComponent %mass          () &
               &                                                *centralBlackHoleComponent%mass          () &
               &                                                )                                           &
               &               /                                 binaryBlackHoleComponent %radialPosition()
          kineticEnergyChange=0.4d0*massRatioIntruder*bindingEnergy
          ejectedInstance    =binaryInstance
          newBinaryInstance  =tripleInstance
          ejectedBlackHoleComponent   => binaryBlackHoleComponent
          newBinaryBlackHoleComponent => tripleBlackHoleComponent
       end if
    else
       ! This latter case can be referred to as head-on collision. 
       newRadius             =0.53d0*tripleBlackHoleComponent%radialPosition()
       call tripleBlackHoleComponent%radialPositionSet(newRadius)
       bindingEnergy         = gravitationalConstantGalacticus*(                                            &
            &                                                    binaryBlackHoleComponent %mass          () &
            &                                                   *centralBlackHoleComponent%mass          () &
            &                                                  )                                            &
            &                 /                                  binaryBlackHoleComponent %radialPosition()
       kineticEnergyChange=0.9d0*massRatioIntruder*bindingEnergy
       ejectedInstance    =binaryInstance
       newBinaryInstance  =tripleInstance
       ejectedBlackHoleComponent   => binaryBlackHoleComponent
       newBinaryBlackHoleComponent => tripleBlackHoleComponent
    end if          
    ! First we find the lightest black hole and tag it as being ejected.
    massEjected= ejectedBlackHoleComponent  %mass()
    massBinary = newBinaryBlackHoleComponent%mass() &
         &      +centralBlackHoleComponent  %mass()
    velocityEjected=dsqrt(kineticEnergyChange/(1.0d0+massEjected/massBinary )/massEjected*2.0d0)
    velocityBinary =dsqrt(kineticEnergyChange/(1.0d0+massBinary /massEjected)/massBinary *2.0d0)
    ! Determine whether the ejected black hole is actualy ejected.
    removeEjected=(0.5d0*velocityEjected**2+Galactic_Structure_Potential(thisNode,ejectedBlackHoleComponent%radialPosition()) >&
         & 0.0d0)
    ! Determine whether the binary black hole is ejected.
    removeBinary=(0.5d0*velocityBinary**2+Galactic_Structure_Potential(thisNode,newBinaryBlackHoleComponent%radialPosition())&
         &-Galactic_Structure_Potential(thisNode,newBinaryBlackHoleComponent%radialPosition(),componentType&
         &=componentTypeBlackHole) > 0.0d0)
    ! Remove the binary black hole from the list if required.
    if (removeBinary) then
       ! Set the central black hole as a zero mass component.
       call centralBlackHoleComponent%massSet(0.0d0)
       call centralBlackHoleComponent%spinSet(0.0d0)
       ! Remove the binary black hole.
       call thisNode%blackHoleRemove(newBinaryInstance)
       ! If this removal has changed the position of the ejected black hole in the list then update its index.
       if (ejectedInstance > newBinaryInstance) ejectedInstance=ejectedInstance-1
    end if
    ! Remove the ejected black hole from the list if required.
    if (removeEjected) call thisNode%blackHoleRemove(ejectedInstance)
    return
  end subroutine Node_Component_Black_Hole_Standard_Triple_Interaction
  
  subroutine Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(thisBlackHoleComponent,accretionRateSpheroid&
       &,accretionRateHotHalo)
    !% Returns the rate of mass accretion onto the black hole in {\tt thisNode}.
    use Cosmological_Parameters
    use Bondi_Hoyle_Lyttleton_Accretion
    use Galactic_Structure_Densities
    use Galactic_Structure_Options
    use Ideal_Gases_Thermodynamics
    use Black_Hole_Fundamentals
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use Accretion_Disks
    use Hot_Halo_Temperature_Profile
    use Memory_Management
    use Black_Hole_Binary_Separations
    implicit none
    class           (nodeComponentBlackHole), pointer, intent(inout) :: thisBlackHoleComponent
    double precision                        ,          intent(  out) :: accretionRateSpheroid,accretionRateHotHalo
    type            (treeNode              ), pointer                :: thisNode
    class           (nodeComponentSpheroid ), pointer                :: thisSpheroidComponent
    class           (nodeComponentHotHalo  ), pointer                :: thisHotHaloComponent
    double precision                        , parameter              :: gasDensityMinimum=1.0d0 ! Lowest gas density to consider when computing accretion rates onto black hole (in units of M_Solar/Mpc^3).
    double precision                                                 :: blackHoleMass,gasDensity,relativeVelocity,accretionRadius&
         &,jeansLength ,radiativeEfficiency,position(3),hotHaloTemperature,hotModeFraction,accretionRateMaximum

    ! Get the host node.
    thisNode => thisBlackHoleComponent%host()

    ! Get black hole mass.
    blackHoleMass=thisBlackHoleComponent%mass()
    ! Check black hole mass is positive.
    if (blackHoleMass > 0.0d0) then
       ! Compute the relative velocity of black hole and gas. We assume that relative motion arises only from the radial
       ! migration of the black hole.
       relativeVelocity=Black_Hole_Binary_Separation_Growth_Rate(thisBlackHoleComponent)*Mpc_per_km_per_s_To_Gyr       
       ! Contribution from spheroid:
       ! Get the accretion radius. We take this to be the larger of the Bondi-Hoyle radius and the current radius position of
       ! the black hole.
       accretionRadius=max(                                                                                              &
            &               Bondi_Hoyle_Lyttleton_Accretion_Radius(blackHoleMass,bondiHoyleAccretionTemperatureSpheroid) &
            &              ,thisBlackHoleComponent%radialPosition()                                                      &
            &             )
       ! Set the position.
       position=[accretionRadius,0.0d0,0.0d0]
       ! Get density of gas at the galactic center.
       gasDensity=Galactic_Structure_Density(thisNode,position,coordinateSystem=coordinateSystemSpherical,massType&
            &=massTypeGaseous,componentType=componentTypeSpheroid)
       ! Check if we have a non-negligible gas density.
       if (gasDensity > gasDensityMinimum) then
          ! Get the spheroid component.
          thisSpheroidComponent => thisNode%spheroid()
          ! Get the Jeans length scale.
          jeansLength=Ideal_Gas_Jeans_Length(bondiHoyleAccretionTemperatureSpheroid,gasDensity)
          ! Limit the smoothing scale to the scale of the spheroid.
          jeansLength=min(jeansLength,thisSpheroidComponent%radius())
          ! If the Jeans length exceeds the Bondi-Hoyle-Lyttleton accretion radius, then recompute gas density for a larger
          ! radius, as the gas should be smoothly distributed on scales below the Jeans length.
          if (jeansLength > accretionRadius) then
             ! Set the position.
             position=[jeansLength,0.0d0,0.0d0]
             ! Get density of gas at the galactic center.
             gasDensity=Galactic_Structure_Density(thisNode,position,coordinateSystem=coordinateSystemCylindrical,massType&
                  &=massTypeGaseous,componentType=componentTypeSpheroid)
          end if
          ! Compute the accretion rate.
          accretionRateSpheroid=bondiHoyleAccretionEnhancementSpheroid*Bondi_Hoyle_Lyttleton_Accretion_Rate(blackHoleMass&
               &,gasDensity ,relativeVelocity,bondiHoyleAccretionTemperatureSpheroid)             
          ! Get the radiative efficiency of the accretion.
          radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisBlackHoleComponent,accretionRateSpheroid)
          ! Limit the accretion rate to the Eddington limit.
          if (radiativeEfficiency > 0.0d0) accretionRateSpheroid=min(accretionRateSpheroid&
               &,Black_Hole_Eddington_Accretion_Rate(thisBlackHoleComponent) /radiativeEfficiency)
       else
          ! Gas density is negative - set zero accretion rate.
          accretionRateSpheroid=0.0d0
       end if
       ! Contribution from hot halo:
       ! Get the hot halo component.
       thisHotHaloComponent => thisNode%hotHalo()
       ! Get halo gas temperature.
       hotHaloTemperature=Hot_Halo_Temperature(thisNode,radius=0.0d0)
       ! Get the accretion radius.
       accretionRadius=Bondi_Hoyle_Lyttleton_Accretion_Radius(blackHoleMass,hotHaloTemperature)
       accretionRadius=min(accretionRadius,thisHotHaloComponent%outerRadius())          
       ! Set the position.
       position=[accretionRadius,0.0d0,0.0d0]
       ! Find the fraction of gas in the halo which is in the hot mode. Set this to unity if hot/cold mode is not to be
       ! considered.
       select case (bondiHoyleAccretionHotModeOnly)
       case (.true.)
          hotModeFraction=Hot_Mode_Fraction(thisNode)
       case (.false.)
          hotModeFraction=1.0d0
       end select             
       ! Get density of gas at the galactic center - scaled by the fraction in the hot accretion mode.
       gasDensity=hotModeFraction*Galactic_Structure_Density(thisNode,position,coordinateSystem=coordinateSystemSpherical&
            &,massType=massTypeGaseous,componentType=componentTypeHotHalo)
       ! Check if we have a non-zero gas density.
       if (gasDensity > gasDensityMinimum) then
          ! Compute the accretion rate.
          accretionRateHotHalo=bondiHoyleAccretionEnhancementHotHalo*Bondi_Hoyle_Lyttleton_Accretion_Rate(blackHoleMass&
               &,gasDensity,relativeVelocity,hotHaloTemperature,accretionRadius)
          ! Limit the accretion rate to the total mass of the hot halo, divided by the sound crossing time.
          accretionRateMaximum=thisHotHaloComponent%mass()/(thisHotHaloComponent%outerRadius()/(kilo*gigaYear/megaParsec)&
               &/Ideal_Gas_Sound_Speed(hotHaloTemperature))
          accretionRateHotHalo=min(accretionRateHotHalo,accretionRateMaximum)
          ! Get the radiative efficiency of the accretion.
          radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisBlackHoleComponent,accretionRateHotHalo)             
          ! Limit the accretion rate to the Eddington limit.
          if (radiativeEfficiency > 0.0d0) accretionRateHotHalo=min(accretionRateHotHalo&
               &,Black_Hole_Eddington_Accretion_Rate(thisBlackHoleComponent)/radiativeEfficiency)
       else
          ! No gas density, so zero accretion rate.
          accretionRateHotHalo=0.0d0
       end if
    else
       accretionRateSpheroid=0.0d0
       accretionRateHotHalo =0.0d0
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Mass_Accretion_Rate

  subroutine Node_Component_Black_Hole_Standard_Create(thisNode)
    !% Creates a black hole component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type (treeNode              ), pointer, intent(inout) :: thisNode
    class(nodeComponentBlackHole), pointer                :: thisBlackHoleComponent
  
    ! Create the black hole.
    thisBlackHoleComponent => thisNode%blackHole(autoCreate=.true.)
    ! Set to the seed mass.
    call thisBlackHoleComponent%          massSet(thisBlackHoleComponent%massSeed())
    call thisBlackHoleComponent%          spinSet(thisBlackHoleComponent%spinSeed())
    call thisBlackHoleComponent%radialPositionSet(                            0.0d0)
    return
  end subroutine Node_Component_Black_Hole_Standard_Create

  !# <mergerTreeOutputNames>
  !#  <unitName>Node_Component_Black_Hole_Standard_Output_Names</unitName>
  !#  <sortName>Node_Component_Black_Hole_Standard_Output</sortName>
  !# </mergerTreeOutputNames>
  subroutine Node_Component_Black_Hole_Standard_Output_Names(thisNode,integerProperty,integerPropertyNames&
       &,integerPropertyComments,integerPropertyUnitsSI ,doubleProperty,doublePropertyNames,doublePropertyComments&
       &,doublePropertyUnitsSI,time)
    !% Set names of black hole properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    use ISO_Varying_String
    implicit none
    type            (treeNode), intent(inout), pointer      :: thisNode
    double precision          , intent(in   )               :: time
    integer                   , intent(inout)               :: integerProperty,doubleProperty
    character       (len=*   ), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision          , intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI
    
    if (Node_Component_Black_Hole_Standard_Matches(thisNode)) then
       !@ <outputPropertyGroup>
       !@   <name>blackHole</name>
       !@   <description>Black hole properities</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>blackHoleCount</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Number of super-massive black holes in the galaxy.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>blackHole</group>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='blackHoleCount'
       integerPropertyComments(integerProperty)='Number of super-massive black holes in the galaxy.'
       integerPropertyUnitsSI (integerProperty)=0.0d0
       if (blackHoleOutputAccretion) then
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>blackHoleAccretionRate</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>Rest-mass accretion rate onto the black hole.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@   <group>blackHole</group>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='blackHoleAccretionRate'
          doublePropertyComments(doubleProperty)='Rest-mass accretion rate onto the black hole.'
          doublePropertyUnitsSI (doubleProperty)=massSolar/gigaYear
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>blackHoleJetPower</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>Power of the black hole-driven jet.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@   <group>blackHole</group>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='blackHoleJetPower'
          doublePropertyComments(doubleProperty)='Power of the black hole-driven jet.'
          doublePropertyUnitsSI (doubleProperty)=massSolar*kilo**2/gigaYear
          doubleProperty=doubleProperty+1
          !@ <outputProperty>
          !@   <name>blackHoleRadiativeEfficiency</name>
          !@   <datatype>real</datatype>
          !@   <cardinality>0..1</cardinality>
          !@   <description>The radiative efficiency of the black hole accretion system.</description>
          !@   <label>???</label>
          !@   <outputType>nodeData</outputType>
          !@   <group>blackHole</group>
          !@ </outputProperty>
          doublePropertyNames   (doubleProperty)='blackHoleRadiativeEfficiency'
          doublePropertyComments(doubleProperty)='The radiative efficiency of the black hole accretion system.'
          doublePropertyUnitsSI (doubleProperty)=0.0d0
       end if
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Node_Component_Black_Hole_Standard_Output_Count</unitName>
  !#  <sortName>Node_Component_Black_Hole_Standard_Output</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Node_Component_Black_Hole_Standard_Output_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of black hole properties to be written to the the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: integerPropertyCount,doublePropertyCount
    integer                   , parameter              :: extraPropertyCount=3

    if (Node_Component_Black_Hole_Standard_Matches(thisNode)) then
       integerPropertyCount=integerPropertyCount+1
       if (blackHoleOutputAccretion) doublePropertyCount=doublePropertyCount+extraPropertyCount
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Output</unitName>
  !#  <sortName>Node_Component_Black_Hole_Standard_Output</sortName>
  !# </mergerTreeOutputTask>
  subroutine Node_Component_Black_Hole_Standard_Output(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store black hole properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    use Accretion_Disks
    implicit none
    double precision                        , intent(in   )          :: time
    type            (treeNode              ), intent(inout), pointer :: thisNode  
    integer                                 , intent(inout)          :: integerProperty,integerBufferCount,doubleProperty&
         &,doubleBufferCount
    integer         (kind=kind_int8        ), intent(inout)          :: integerBuffer(:,:)
    double precision                        , intent(inout)          :: doubleBuffer (:,:)
    class           (nodeComponentBlackHole),                pointer :: thisBlackHoleComponent
    double precision                                                 :: restMassAccretionRate,accretionRateSpheroid&
         &,accretionRateHotHalo

    if (Node_Component_Black_Hole_Standard_Matches(thisNode)) then
       ! Store the properties.
       if (blackHoleOutputAccretion) then
          ! Get the black hole component.
          thisBlackHoleComponent => thisNode%blackHole(instance=1)
          ! Get the rest mass accretion rate.
          call Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(thisBlackHoleComponent,accretionRateSpheroid&
               &,accretionRateHotHalo)
          restMassAccretionRate=accretionRateSpheroid+accretionRateHotHalo
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=restMassAccretionRate
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Accretion_Disk_Jet_Power           (thisBlackHoleComponent,restMassAccretionRate)
          doubleProperty=doubleProperty+1
          doubleBuffer(doubleBufferCount,doubleProperty)=Accretion_Disk_Radiative_Efficiency(thisBlackHoleComponent,restMassAccretionRate)
       end if
       ! Count number of black holes associated with this galaxy.
       integerProperty=integerProperty+1
       integerBuffer(integerBufferCount,integerProperty)=thisNode%blackHoleCount()
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Output

  logical function Node_Component_Black_Hole_Standard_Matches(thisNode)
    !% Return true if the black hole component of {\tt thisNode} is a match to the standard implementation.
    implicit none
    type (treeNode              ), intent(inout), pointer :: thisNode
    class(nodeComponentBlackHole),                pointer :: thisBlackHoleComponent

    ! Get the black hole component.
    thisBlackHoleComponent => thisNode%blackHole()
    ! Ensure that it is of the standard class.
    Node_Component_Black_Hole_Standard_Matches=.false.
    select type (thisBlackHoleComponent)
    class is (nodeComponentBlackHoleStandard)
       Node_Component_Black_Hole_Standard_Matches=.true.
    type  is (nodeComponentBlackHole        )
       Node_Component_Black_Hole_Standard_Matches=defaultBlackHoleComponent%standardIsActive()
    end select
    return
  end function Node_Component_Black_Hole_Standard_Matches

  subroutine Node_Component_Black_Hole_Standard_Output_Merger(thisNode,massBlackHole1,massBlackHole2)
    !% Outputs properties of merging black holes.
    use Galacticus_Nodes
    use IO_HDF5
    use Galacticus_HDF5
    use Memory_Management
    use Kind_Numbers
    use ISO_Varying_String
    use String_Handling
    use Merger_Tree_Active
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode   
    double precision                    , intent(in   )          :: massBlackHole1,massBlackHole2
    class           (nodeComponentBasic),                pointer :: thisBasicComponent
    type            (hdf5Object        )                         :: mergersGroup

    ! Exit if merger data is not to be output.
    if (.not.blackHoleOutputMergers) return

    ! Ignore mergers with zero mass black holes.
    if (massBlackHole2 <= 0.0d0    ) return

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()

    ! Open the group to which black hole mergers should be written.
    !$omp critical (Node_Component_Black_Hole_Standard_Output_Merger)
    mergersGroup=galacticusOutputFile%openGroup("blackHoleMergers","Black hole mergers data.")
    ! Append to the datasets.
    call mergersGroup%writeDataset([massBlackHole1           ],"massBlackHole1","Mass of the first merging black hole." ,appendTo=.true.)
    call mergersGroup%writeDataset([massBlackHole2           ],"massBlackHole2","Mass of the second merging black hole.",appendTo=.true.)
    call mergersGroup%writeDataset([thisBasicComponent%time()],"timeOfMerger"  ,"The time of the black hole merger."    ,appendTo=.true.)
    call mergersGroup%writeDataset([activeTreeWeight         ],"volumeWeight"  ,"The weight for the black hole merger." ,appendTo=.true.)
    ! Close the group.
    call mergersGroup%close()  
    !$omp end critical (Node_Component_Black_Hole_Standard_Output_Merger)
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Merger

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Black_Hole_Standard_Output_Properties</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Black_Hole_Standard_Output_Properties(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Output properties for all black holes in {\tt thisNode}.
    use Galacticus_Nodes
    use IO_HDF5
    use Galacticus_HDF5
    use Memory_Management
    use Kind_Numbers
    use ISO_Varying_String
    use String_Handling
    use Dark_Matter_Profiles
    use Cosmology_Functions
    use Black_Hole_Binary_Separations
    use Black_Hole_Binary_Separations_Standard
    use Accretion_Disks
    use Cooling_Radii
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode              ), intent(inout), pointer      :: thisNode
    integer         (kind=kind_int8        ), intent(in   )               :: treeIndex
    integer                                 , intent(in   )               :: iOutput
    logical                                 , intent(in   )               :: nodePassesFilter
    class           (nodeComponentBlackHole),                pointer      :: thisBlackHoleComponent
    integer         (kind=kind_int8        ), allocatable,   dimension(:) :: nodeIndex,mergerTreeIndex
    double precision                        , allocatable,   dimension(:) :: massAccretionRate,radiativeEfficiency,mass,spin&
         &,timescale,radius
    double precision                                                      :: accretionRateSpheroid,accretionRateHotHalo
    integer                                                               :: instance,blackHoleCount
    type            (hdf5Object            )                              :: blackHolesGroup,outputGroup
    type            (varying_string        )                              :: groupName

    ! If black hole output was requested , output their properties.
    if (nodePassesFilter .and. blackHoleOutputData) then
       ! Get a count of the number of black holes present.
       blackHoleCount=thisNode%blackHoleCount()
       ! Open the output group.
       !$omp critical (HDF5_Access)
       blackHolesGroup=galacticusOutputFile%openGroup("blackHole","Black hole data.")
       groupName="Output"
       groupName=groupName//iOutput
       outputGroup=blackHolesGroup%openGroup(char(groupName),"Properties of black holes for all trees at each output.")  
       !$omp end critical (HDF5_Access)
       ! Allocate array to store profile.
       call Alloc_Array(radius             ,[blackHoleCount])
       call Alloc_Array(spin               ,[blackHoleCount])
       call Alloc_Array(mass               ,[blackHoleCount])
       call Alloc_Array(timescale          ,[blackHoleCount])
       call Alloc_Array(massAccretionRate  ,[blackHoleCount])
       call Alloc_Array(radiativeEfficiency,[blackHoleCount])
       call Alloc_Array(nodeIndex          ,[blackHoleCount])
       call Alloc_Array(mergerTreeIndex    ,[blackHoleCount])
       ! Construct arrays of black hole properties.
       do instance=1,blackHoleCount
          thisBlackHoleComponent => thisNode%blackHole(instance=instance)
          call  Node_Component_Black_Hole_Standard_Mass_Accretion_Rate(thisBlackHoleComponent,accretionRateSpheroid&
               &,accretionRateHotHalo)
          mass               (instance)=thisBlackHoleComponent%mass()
          spin               (instance)=thisBlackHoleComponent%spin()
          radius             (instance)=thisBlackHoleComponent%radialPosition() 
          massAccretionRate  (instance)=accretionRateSpheroid+accretionRateHotHalo
          radiativeEfficiency(instance)=Accretion_Disk_Radiative_Efficiency(thisBlackHoleComponent,massAccretionRate(instance))
          nodeIndex          (instance)=thisNode%index()
          mergerTreeIndex    (instance)=treeIndex     
          if (instance > 1) then
             if (Black_Hole_Binary_Separation_Growth_Rate(thisBlackHoleComponent) /= 0.0d0 )then
                timescale(instance)=-thisBlackHoleComponent%radialPosition()                 &
                     &              /Black_Hole_Binary_Separation_Growth_Rate(thisBlackHoleComponent) 
             else
                timescale(instance)=0.0d0
             end if
          else
             timescale   (instance)=0.0d0
          end if
       end do
       ! Write dataset to the group, first the arrays containing all data.
       !$omp critical (HDF5_Access)
       call outputGroup%writeDataset(mass,               "mass"               ,"The black hole masses.",                appendTo=.true.)
       call outputGroup%writeDataset(spin,               "spin"               ,"The black hole spins.",                 appendTo=.true.)
       call outputGroup%writeDataset(radius,             "radius"             ,"The black hole radial positions.",      appendTo=.true.)
       call outputGroup%writeDataset(timescale,          "timescale"          ,"The black hole timescales for merger.", appendTo=.true.)
       call outputGroup%writeDataset(radiativeEfficiency,"radiativeEfficiency","The black hole radiative efficiencies.",appendTo=.true.)
       call outputGroup%writeDataset(massAccretionRate,  "accretionRate"      ,"The black hole accretion rates.",       appendTo=.true.)
       call outputGroup%writeDataset(nodeIndex,          "nodeIndex"          ,"The black hole host galaxy inices.",    appendTo=.true.)
       call outputGroup%writeDataset(mergerTreeIndex,    "mergerTreeIndex"    ,"The black hole merger tree indices.",   appendTo=.true.)

       ! Deallocatate profile arrays.
       call Dealloc_Array(mass               )
       call Dealloc_Array(spin               )
       call Dealloc_Array(radius             )
       call Dealloc_Array(timescale          )
       call Dealloc_Array(radiativeEfficiency)
       call Dealloc_Array(massAccretionRate  )
       call Dealloc_Array(nodeIndex          )
       call Dealloc_Array(mergerTreeIndex    )
       call outputGroup    %close()
       call blackHolesGroup%close()   
       !$omp end critical (HDF5_Access)
    end if
    return
  end subroutine Node_Component_Black_Hole_Standard_Output_Properties

  double precision function Hot_Mode_Fraction(thisNode)
    !% A simple interpolating function which is used as a measure of the fraction of a halo which is in the hot accretion mode.
    use Cooling_Radii
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: coolingRadiusFractionalTransitionMinimum=0.9d0
    double precision, parameter              :: coolingRadiusFractionalTransitionMaximum=1.0d0
    double precision                         :: x,coolingRadiusFractional

    coolingRadiusFractional=Cooling_Radius(thisNode)/Dark_Matter_Halo_Virial_Radius(thisNode)
    if      (coolingRadiusFractional < coolingRadiusFractionalTransitionMinimum) then
       Hot_Mode_Fraction=1.0d0
    else if (coolingRadiusFractional > coolingRadiusFractionalTransitionMaximum) then
       Hot_Mode_Fraction=0.0d0
    else
       x=      (coolingRadiusFractional                 -coolingRadiusFractionalTransitionMinimum) &
            & /(coolingRadiusFractionalTransitionMaximum-coolingRadiusFractionalTransitionMinimum)
       Hot_Mode_Fraction=x**2*(2.0d0*x-3.0d0)+1.0d0
    end if
    return
  end function Hot_Mode_Fraction
  
end module Node_Component_Black_Hole_Standard

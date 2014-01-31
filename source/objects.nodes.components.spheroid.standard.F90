!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the standard spheroid component.

module Node_Component_Spheroid_Standard
  !% Implements the standard spheroid component.
  use Galacticus_Nodes
  use Histories
  use Stellar_Population_Properties
  use Node_Component_Spheroid_Standard_Data
  implicit none
  private
  public :: Node_Component_Spheroid_Standard_Rate_Compute     , Node_Component_Spheroid_Standard_Scale_Set                    , &
       &    Node_Component_Spheroid_Standard_Radius_Solver    , Node_Component_Spheroid_Standard_Star_Formation_History_Output, &
       &    Node_Component_Spheroid_Standard_Initialize       , Node_Component_Spheroid_Standard_Post_Evolve                  , &
       &    Node_Component_Spheroid_Standard_Pre_Evolve       , Node_Component_Spheroid_Standard_Radius_Solver_Plausibility   , &
       &    Node_Component_Spheroid_Standard_Satellite_Merging

  !# <component>
  !#  <class>spheroid</class>
  !#  <name>standard</name>
  !#  <isDefault>yes</isDefault>
  !#  <createFunction isDeferred="true" />
  !#  <properties>
  !#   <property>
  !#     <name>isInitialized</name>
  !#     <type>logical</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#   <property>
  !#     <name>massStellar</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of stars in the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesStellar</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the stellar phase of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>massGas</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>abundancesGas</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>angularMomentum</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>radius</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="megaParsec" comment="Scale length of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>halfMassRadius</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <getFunction>Node_Component_Spheroid_Standard_Half_Mass_Radius</getFunction>
  !#   </property>
  !#   <property>
  !#     <name>velocity</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="kilo" comment="Circular velocity at the scale length of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>starFormationRate</name>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" createIfNeeded="true" makeGeneric="true" />
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <output condition="[[spheroidOutputStarFormationRate]]" unitsInSI="massSolar/gigaYear" comment="Star formation rate of the standard spheroid."/>
  !#   </property>
  !#   <property>
  !#     <name>luminositiesStellar</name>
  !#     <type>stellarLuminosities</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#     <output unitsInSI="luminosityZeroPointAB" comment="Luminosity of spheroid stars."/>
  !#   </property>
  !#   <property>
  !#     <name>stellarPropertiesHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#   </property>
  !#   <property>
  !#     <name>starFormationHistory</name>
  !#     <type>history</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="true" />
  !#   </property>
  !#   <property>
  !#     <name>massGasSink</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" />
  !#   </property>
  !#   <property>
  !#     <name>energyGasInput</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <isVirtual>true</isVirtual>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" />
  !#   </property>
  !#  </properties>
  !#  <bindings>
  !#   <binding method="enclosedMass"          function="Node_Component_Spheroid_Standard_Enclosed_Mass"           bindsTo="component" />
  !#   <binding method="density"               function="Node_Component_Spheroid_Standard_Density"                 bindsTo="component" />
  !#   <binding method="rotationCurve"         function="Node_Component_Spheroid_Standard_Rotation_Curve"          bindsTo="component" />
  !#   <binding method="rotationCurveGradient" function="Node_Component_Spheroid_Standard_Rotation_Curve_Gradient" bindsTo="component" />
  !#   <binding method="potential"             function="Node_Component_Spheroid_Standard_Potential"               bindsTo="component" />
  !#  </bindings>
  !#  <functions>objects.nodes.components.spheroid.standard.bound_functions.inc</functions>
  !# </component>

  ! Internal count of abundances.
  integer                                     :: abundancesCount

  ! Storage for the star formation history time range, used whene extending this range.
  double precision, allocatable, dimension(:) :: starFormationHistoryTemplate
  !$omp threadprivate(starFormationHistoryTemplate)
  ! Parameters controlling the physical implementation.
  double precision                            :: spheroidEnergeticOutflowMassRate            , spheroidOutflowTimescaleMinimum
  double precision                            :: spheroidAngularMomentumAtScaleRadius        , spheroidMassToleranceAbsolute

  ! Record of whether this module has been initialized.
  logical                                     :: moduleInitialized                   =.false.

contains

  !# <mergerTreePreTreeConstructionTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Initialize</unitName>
  !# </mergerTreePreTreeConstructionTask>
  subroutine Node_Component_Spheroid_Standard_Initialize()
    !% Initializes the tree node standard spheroid methods module.
    use Input_Parameters
    use Stellar_Luminosities_Structure
    use Abundances_Structure
    use ISO_Varying_String
    use Galacticus_Error
    use Memory_Management
    implicit none
    type            (nodeComponentSpheroidStandard) :: spheroidStandardComponent
    double precision                                :: spheroidAngularMomentumAtScaleRadiusDefault, spheroidMassDistributionDensityMomentum2, &
         &                                             spheroidMassDistributionDensityMomentum3   , spheroidSersicIndex
    logical                                         :: densityMoment2IsInfinite                   , densityMoment3IsInfinite
    type            (varying_string               ) :: spheroidMassDistributionName

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Spheroid_Standard_Initialize)
    if (defaultSpheroidComponent%standardIsActive().and..not.moduleInitialized) then

       ! Get number of abundance properties.
       abundancesCount  =Abundances_Property_Count            ()

       ! Create the spheroid mass distribution.
       !@ <inputParameter>
       !@   <name>spheroidMassDistribution</name>
       !@   <defaultValue>hernquist</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The type of mass distribution to use for the standard spheroid component.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidMassDistribution',spheroidMassDistributionName,defaultValue="hernquist")
       spheroidMassDistribution => Mass_Distribution_Create(char(spheroidMassDistributionName))
       select type (spheroidMassDistribution)
       type is (massDistributionHernquist)
          call spheroidMassDistribution%initialize(isDimensionless=.true.)
       type is (massDistributionSersic   )
          !@ <inputParameter>
          !@   <name>spheroidSersicIndex</name>
          !@   <defaultValue>$4$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The S\'ersic index to use for the spheroid component mass distribution.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('spheroidSersicIndex',spheroidSersicIndex,defaultValue=4.0d0)
          call spheroidMassDistribution%initialize(index=spheroidSersicIndex,isDimensionless=.true.)
       class default
          call Galacticus_Error_Report('Node_Component_Spheroid_Standard_Initialize','unsupported mass distribution')
       end select

       ! Bind deferred functions.
       call spheroidStandardComponent% starFormationRateFunction(Node_Component_Spheroid_Standard_Star_Formation_Rate  )
       call spheroidStandardComponent%energyGasInputRateFunction(Node_Component_Spheroid_Standard_Energy_Gas_Input_Rate)
       call spheroidStandardComponent%   massGasSinkRateFunction(Node_Component_Spheroid_Standard_Mass_Gas_Sink_Rate   )
       call spheroidStandardComponent%         createFunctionSet(Node_Component_Spheroid_Standard_Initializor          )

       ! Determine the specific angular momentum at the scale radius in units of the mean specific angular
       ! momentum of the spheroid. This is equal to the ratio of the 2nd to 3rd radial moments of the density
       ! distribution (assuming a flat rotation curve).
       spheroidMassDistributionDensityMomentum2=spheroidMassDistribution%densityRadialMoment(2.0d0,isInfinite=densityMoment2IsInfinite)
       spheroidMassDistributionDensityMomentum3=spheroidMassDistribution%densityRadialMoment(3.0d0,isInfinite=densityMoment3IsInfinite)
       if (densityMoment2IsInfinite.or.densityMoment3IsInfinite) then
          ! One of the moments is infinte, so we can not compute the appropriate ratio. Simply assume a value
          ! of 0.5 as a default.
          spheroidAngularMomentumAtScaleRadiusDefault=0.5d0
       else
          ! Moments are well-defined, so compute their ratio.
          spheroidAngularMomentumAtScaleRadiusDefault=spheroidMassDistributionDensityMomentum2&
               &/spheroidMassDistributionDensityMomentum3
       end if

       ! Read parameters controlling the physical implementation.
       !@ <inputParameter>
       !@   <name>spheroidAngularMomentumAtScaleRadius</name>
       !@   <defaultValue>$I_2/I_3$ where $I_n=\int_0^\infty \rho(r) r^n {\rm d}r$, where $\rho(r)$ is the spheroid density profile, unless either $I_2$ or $I_3$ is infinite, in which case a default of $1/2$ is used instead</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The assumed ratio of the specific angular momentum at the scale radius to the mean specific angular momentum of the standard spheroid component.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidAngularMomentumAtScaleRadius',spheroidAngularMomentumAtScaleRadius,defaultValue=spheroidAngularMomentumAtScaleRadiusDefault)
       !@ <inputParameter>
       !@   <name>spheroidEnergeticOutflowMassRate</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The proportionallity factor relating mass outflow rate from the spheroid to the energy input rate divided by $V_{\rm spheroid}^2$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidEnergeticOutflowMassRate',spheroidEnergeticOutflowMassRate,defaultValue=1.0d-2)
       !@ <inputParameter>
       !@   <name>spheroidMassToleranceAbsolute</name>
       !@   <defaultValue>$10^{-6} M_\odot$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The mass tolerance used to judge whether the spheroid is physically plausible.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidMassToleranceAbsolute',spheroidMassToleranceAbsolute,defaultValue=1.0d-6)
       !@ <inputParameter>
       !@   <name>spheroidOutflowTimescaleMinimum</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The minimum timescale (in units of the spheroid dynamical time) on which outflows may deplete gas in the spheroid.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidOutflowTimescaleMinimum',spheroidOutflowTimescaleMinimum,defaultValue=1.0d-3)

       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Spheroid_Standard_Initialize)
    return
  end subroutine Node_Component_Spheroid_Standard_Initialize

  !# <preEvolveTask>
  !# <unitName>Node_Component_Spheroid_Standard_Pre_Evolve</unitName>
  !# </preEvolveTask>
  subroutine Node_Component_Spheroid_Standard_Pre_Evolve(thisNode)
    !% Ensure the spheroid has been initialized.
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentSpheroid)               , pointer :: thisSpheroidComponent

    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()
    ! Check if a standard spheroid component exists.
    select type (thisSpheroidComponent)
    class is (nodeComponentSpheroidStandard)
       ! Initialize the spheroid.
       call Node_Component_Spheroid_Standard_Initializor(thisSpheroidComponent)
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Pre_Evolve

  !# <postEvolveTask>
  !# <unitName>Node_Component_Spheroid_Standard_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Spheroid_Standard_Post_Evolve(thisNode)
    !% Trim histories attached to the spheroid.
    use Galacticus_Display
    use String_Handling
    use ISO_Varying_String
    use Abundances_Structure
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    class           (nodeComponentSpheroid )               , pointer :: thisSpheroidComponent
    class           (nodeComponentBasic    )               , pointer :: thisBasicComponent
    double precision                        , save                   :: fractionalErrorMaximum  =0.0d0
    double precision                                                 :: fractionalError               , specificAngularMomentum, &
         &                                                              spheroidMass
    character       (len=20                )                         :: valueString
    type            (varying_string        )                         :: message
    type            (history               )                         :: stellarPropertiesHistory

    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()
    ! Check if an exponential spheroid component exists.
    select type (thisSpheroidComponent)
       class is (nodeComponentSpheroidStandard)

          ! Trim the stellar populations properties future history.
       thisBasicComponent => thisNode%basic()
       stellarPropertiesHistory=thisSpheroidComponent%stellarPropertiesHistory()
       call stellarPropertiesHistory%trim(thisBasicComponent%time())
       call thisSpheroidComponent%stellarPropertiesHistorySet(stellarPropertiesHistory)

       ! Trap negative gas masses.
       if (thisSpheroidComponent%massGas() < 0.0d0) then

          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(thisSpheroidComponent%massGas    ()) &
               &          /(                                          &
               &             abs(thisSpheroidComponent%massGas    ()) &
               &            +abs(thisSpheroidComponent%massStellar()) &
               &           )
          !$omp critical (Standard_Spheroid_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: spheroid has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index            = '//thisNode%index() //char(10)
             write (valueString,'(e12.6)') thisSpheroidComponent%massGas    ()
             message=message//'  Spheroid gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') thisSpheroidComponent%massStellar()
             message=message//'  Spheroid stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure         = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Gas mass will be reset to zero (in future cases also).'//char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerically inaccuracy in the ODE solutions.'//char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call Galacticus_Display_Message(message,verbosityWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Spheroid_Post_Evolve_Check)

          ! Get the specific angular momentum of the spheroid material
          spheroidMass= thisSpheroidComponent%massGas    () &
               &       +thisSpheroidComponent%massStellar()
          if (spheroidMass == 0.0d0) then
             specificAngularMomentum=0.0d0
             call thisSpheroidComponent%        massStellarSet(                  0.0d0)
             call thisSpheroidComponent%  abundancesStellarSet(         zeroAbundances)
             call thisSpheroidComponent%luminositiesStellarSet(zeroStellarLuminosities)
          else
             specificAngularMomentum=thisSpheroidComponent%angularMomentum()/spheroidMass
          end if

          ! Reset the gas, abundances and angular momentum of the spheroid.
          call thisSpheroidComponent%        massGasSet(                                                      0.0d0)
          call thisSpheroidComponent%  abundancesGasSet(                                             zeroAbundances)
          call thisSpheroidComponent%angularMomentumSet(specificAngularMomentum*thisSpheroidComponent%massStellar())
       end if

    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Post_Evolve

  subroutine Node_Component_Spheroid_Standard_Mass_Gas_Sink_Rate(self,rate,interrupt,interruptProcedure)
    !% Account for a sink of gaseous material in the standard spheroid.
    use Galacticus_Error
    use Abundances_Structure
    implicit none
    class           (nodeComponentSpheroid       ), intent(inout)                    :: self
    logical                                       , intent(inout), optional          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), optional, pointer :: interruptProcedure
    double precision                              , intent(in   )                    :: rate
    double precision                                                                 :: gasMass           , stellarMass

    ! Trap cases where an attempt is made to add gas via this sink function.
    if (rate > 0.0d0) call Galacticus_Error_Report(                                                        &
         &                                         'Node_Component_Spheroid_Standard_Mass_Gas_Sink_Rate', &
         &                                         'attempt to add mass via sink in standard spheroid'    &
         &                                        )
    ! Return if no adjustment is being made.
    if (rate == 0.0d0) return
    ! Get the gas mass present.
    gasMass    =self%massGas    ()
    ! Get the stellar mass present.
    stellarMass=self%massStellar()
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0 .and. gasMass+stellarMass > 0.0d0) then
       call self%        massGasRate( rate                                              )
       call self%angularMomentumRate((rate/(gasMass+stellarMass))*self%angularMomentum())
       call self%  abundancesGasRate((rate/ gasMass             )*self%abundancesGas  ())
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Mass_Gas_Sink_Rate

  subroutine Node_Component_Spheroid_Standard_Energy_Gas_Input_Rate(self,rate,interrupt,interruptProcedure)
    !% Handles input of energy into the spheroid gas from other components (e.g. black holes). The energy input rate should be in
    !% units of $M_\odot$ km$^2$ s$^{-2}$ Gyr$^{-1}$.
    use Galacticus_Error
    use Abundances_Structure
    implicit none
    class           (nodeComponentSpheroid        ), intent(inout)                    :: self
    logical                                        , intent(inout), optional          :: interrupt
    procedure       (Interrupt_Procedure_Template ), intent(inout), optional, pointer :: interruptProcedure
    double precision                               , intent(in   )                    :: rate
    class           (nodeComponentHotHalo         )                         , pointer :: selfHotHaloComponent
    type            (treeNode                     )                         , pointer :: selfNode
    type            (abundances                   )                                   :: abundancesOutflowRate
    double precision                                                                  :: angularMomentumOutflowRate, gasMass         , &
         &                                                                               massOutflowRate           , spheroidVelocity, &
         &                                                                               stellarMass

    ! Trap cases where an attempt is made to remove energy via this input function.
    if (rate < 0.0d0) call Galacticus_Error_Report(                                                                &
         &                                         'Node_Component_Spheroid_Standard_Energy_Gas_Input_Rate'      , &
         &                                         'attempt to remove energy via input pipe in standard spheroid'  &
         &                                        )

    ! Return if no adjustment is being made.
    if (rate == 0.0d0) return
    ! Get the gas mass present.
    gasMass    =self%massGas    ()
    ! Get the stellar mass present.
    stellarMass=self%massStellar()
    ! If gas is present, adjust the rates.
    if (gasMass > 0.0d0 .and. gasMass+stellarMass > 0.0d0) then
       ! Compute outflow rates of quantities and adjust rates in the spheroid appropriately.
       spheroidVelocity=self%velocity()
       if (spheroidVelocity > 0.0d0) then
          massOutflowRate           =spheroidEnergeticOutflowMassRate*rate/spheroidVelocity**2
          angularMomentumOutflowRate=(massOutflowRate/(gasMass+stellarMass))*self%angularMomentum()
          abundancesOutflowRate     =(massOutflowRate/ gasMass             )*self%abundancesGas  ()
          call self%        massGasRate(-           massOutflowRate)
          call self%angularMomentumRate(-angularMomentumOutflowRate)
          call self%  abundancesGasRate(-     abundancesOutflowRate)
          ! Add outflowing rates to the hot halo component.
          selfNode             => self    %host   ()
          selfHotHaloComponent => selfNode%hotHalo()
          call selfHotHaloComponent%           outflowingMassRate(           massOutflowRate)
          call selfHotHaloComponent%outflowingAngularMomentumRate(angularMomentumOutflowRate)
          call selfHotHaloComponent%     outflowingAbundancesRate(     abundancesOutflowRate)
       end if
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Energy_Gas_Input_Rate

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Spheroid_Standard_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the standard spheroid node mass rate of change.
    use Star_Formation_Feedback_Spheroids
    use Star_Formation_Feedback_Expulsion_Spheroids
    use Abundances_Structure
    use Galactic_Structure_Options
    use Galacticus_Output_Star_Formation_Histories
    use Numerical_Constants_Astronomical
    use Dark_Matter_Halo_Scales
    use Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids
    use Tidal_Stripping_Mass_Loss_Rate_Spheroids
    use Satellites_Tidal_Fields
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    logical                                , intent(inout)          :: interrupt
    procedure       (                     ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroid
    class           (nodeComponentHotHalo )               , pointer :: thisHotHalo
    type            (abundances           ), save                   :: fuelAbundances            , fuelAbundancesRates     , &
         &                                                             stellarAbundancesRates
    !$omp threadprivate(fuelAbundances,stellarAbundancesRates,fuelAbundancesRates)
    double precision                                                :: angularMomentumOutflowRate, energyInputRate         , &
         &                                                             fractionGas               , fractionStellar         , &
         &                                                             fuelMass                  , fuelMassRate            , &
         &                                                             gasMass                   , massLossRate            , &
         &                                                             massOutflowRate           , massOutflowRateFromHalo , &
         &                                                             massOutflowRateToHotHalo  , outflowToHotHaloFraction, &
         &                                                             spheroidDynamicalTime     , spheroidMass            , &
         &                                                             starFormationRate         , stellarMassRate         , &
         &                                                             tidalField                , tidalTorque
    type            (history              )                         :: stellarHistoryRate
    type            (stellarLuminosities  )                         :: luminositiesStellarRates

    ! Get the disk and check that it is of our class.
    thisSpheroid => thisNode%spheroid()
    select type (thisSpheroid)
    class is (nodeComponentSpheroidStandard)

       ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
       if     (    thisSpheroid%angularMomentum() <  0.0d0 &
            & .or. thisSpheroid%radius         () <= 0.0d0 &
            & .or. thisSpheroid%massGas        () <  0.0d0 &
            & ) return

       ! Find the star formation timescale.
       starFormationRate=thisSpheroid%starFormationRate()

       ! Get the available fuel mass.
       fuelMass         =thisSpheroid%massGas          ()

       ! Find the metallicity of the fuel supply.
       fuelAbundances   =thisSpheroid%abundancesGas    ()
       call fuelAbundances%massToMassFraction(fuelMass)

       ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
       stellarHistoryRate=thisSpheroid%stellarPropertiesHistory()
       call Stellar_Population_Properties_Rates(starFormationRate,fuelAbundances,componentTypeSpheroid,thisNode &
            &,stellarHistoryRate,stellarMassRate,stellarAbundancesRates ,luminositiesStellarRates,fuelMassRate&
            &,fuelAbundancesRates,energyInputRate)
       if (stellarHistoryRate%exists()) call thisSpheroid%stellarPropertiesHistoryRate(stellarHistoryRate)

       ! Adjust rates.
       call thisSpheroid%        massStellarRate(         stellarMassRate)
       call thisSpheroid%            massGasRate(            fuelMassRate)
       call thisSpheroid%  abundancesStellarRate(  stellarAbundancesRates)
       call thisSpheroid%      abundancesGasRate(     fuelAbundancesRates)
       call thisSpheroid%luminositiesStellarRate(luminositiesStellarRates)

       ! Record the star formation history.
       stellarHistoryRate=thisSpheroid%starFormationHistory()
       call stellarHistoryRate%reset()
       call Star_Formation_History_Record(thisNode,stellarHistoryRate,fuelAbundances,starFormationRate)
       if (stellarHistoryRate%exists()) call thisSpheroid%starFormationHistoryRate(stellarHistoryRate)

       ! Find rate of outflow of material from the spheroid and pipe it to the outflowed reservoir.
       massOutflowRateToHotHalo=Star_Formation_Feedback_Spheroid_Outflow_Rate          (thisNode,starFormationRate,energyInputRate)
       massOutflowRateFromHalo =Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
       massOutflowRate         =massOutflowRateToHotHalo+massOutflowRateFromHalo
       if (massOutflowRate > 0.0d0) then
          ! Find the fraction of material which outflows to the hot halo.
          outflowToHotHaloFraction=massOutflowRateToHotHalo/massOutflowRate

          ! Get the masses of the spheroid.
          gasMass     =        thisSpheroid%massGas    ()
          spheroidMass=gasMass+thisSpheroid%massStellar()

          ! Limit the outflow rate timescale to a multiple of the dynamical time.
          spheroidDynamicalTime=Mpc_per_km_per_s_To_Gyr*thisSpheroid%radius()/thisSpheroid%velocity()

          ! Limit the mass outflow rate.
          massOutflowRate=min(massOutflowRate,gasMass/spheroidOutflowTimescaleMinimum/spheroidDynamicalTime)
          thisHotHalo => thisNode%hotHalo()
          call thisHotHalo %outflowingMassRate(+massOutflowRate*outflowToHotHaloFraction)
          call thisSpheroid%       massGasRate(-massOutflowRate                         )

          ! Compute the angular momentum outflow rate.
          if (spheroidMass > 0.0d0) then
             angularMomentumOutflowRate=(massOutflowRate/spheroidMass)*thisSpheroid%angularMomentum()
             call thisHotHalo %outflowingAngularMomentumRate(+angularMomentumOutflowRate*outflowToHotHaloFraction)
             call thisSpheroid%          angularMomentumRate(-angularMomentumOutflowRate                         )
          end if

          ! Compute the abundances outflow rate.
          fuelAbundancesRates=thisSpheroid%abundancesGas()
          call fuelAbundancesRates%massToMassFraction(gasMass)
          fuelAbundancesRates=fuelAbundancesRates*massOutflowRate
          call thisHotHalo %outflowingAbundancesRate(+fuelAbundancesRates*outflowToHotHaloFraction)
          call thisSpheroid%       abundancesGasRate(-fuelAbundancesRates                         )
       end if

       ! Apply mass loss rate due to ram pressure stripping.
       if (thisSpheroid%massGas() > 0.0d0) then
          massLossRate=Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid(thisNode)
          if (massLossRate > 0.0d0) then
             thisHotHalo => thisNode%hotHalo()
             call thisSpheroid%                  massGasRate(-massLossRate                                                                                   )
             call thisSpheroid%          angularMomentumRate(-massLossRate*thisSpheroid%angularMomentum()/(thisSpheroid%massGas()+thisSpheroid%massStellar()))
             call thisSpheroid%            abundancesGasRate(-massLossRate*thisSpheroid%abundancesGas  ()/ thisSpheroid%massGas()                            )
             call thisHotHalo %           outflowingMassRate(+massLossRate                                                                                   )
             call thisHotHalo %outflowingAngularMomentumRate(+massLossRate*thisSpheroid%angularMomentum()/(thisSpheroid%massGas()+thisSpheroid%massStellar()))
             call thisHotHalo %outflowingAbundancesRate     (+massLossRate*thisSpheroid%abundancesGas  ()/ thisSpheroid%massGas()                            )
          end if
       end if

       ! Apply mass loss rate due to tidal stripping.
       if (thisSpheroid%massGas()+thisSpheroid%massStellar() > 0.0d0) then
          massLossRate=Tidal_Stripping_Mass_Loss_Rate_Spheroid(thisNode)
          if (massLossRate > 0.0d0) then
             thisHotHalo    => thisNode%hotHalo()
             fractionGas    =  min(1.0d0,max(0.0d0,thisSpheroid%massGas()/(thisSpheroid%massGas()+thisSpheroid%massStellar())))
             fractionStellar=  1.0d0-fractionGas
             if (fractionGas    > 0.0d0 .and. thisSpheroid%massGas    () > 0.0d0) then
                call thisSpheroid%                  massGasRate(-fractionGas    *massLossRate                                                                                     )
                call thisSpheroid%            abundancesGasRate(-fractionGas    *massLossRate*thisSpheroid%abundancesGas    ()/ thisSpheroid%massGas()                            )
                call thisHotHalo %           outflowingMassRate(+fractionGas    *massLossRate                                                                                     )
                call thisHotHalo %outflowingAbundancesRate     (+fractionGas    *massLossRate*thisSpheroid%abundancesGas    ()/ thisSpheroid%massGas()                            )
                call thisHotHalo %outflowingAngularMomentumRate(+fractionGas    *massLossRate*thisSpheroid%angularMomentum  ()/(thisSpheroid%massGas()+thisSpheroid%massStellar()))
             end if
             if (fractionStellar > 0.0d0 .and. thisSpheroid%massStellar() > 0.0d0) then
                call thisSpheroid%              massStellarRate(-fractionStellar*massLossRate                                                                                     )
                call thisSpheroid%        abundancesStellarRate(-fractionStellar*massLossRate*thisSpheroid%abundancesStellar()/                        thisSpheroid%massStellar() )
             end if
             call    thisSpheroid%          angularMomentumRate(-                massLossRate*thisSpheroid%angularMomentum  ()/(thisSpheroid%massGas()+thisSpheroid%massStellar()))
          end if
       end if

       ! Apply tidal heating.
       if (thisNode%isSatellite() .and. thisSpheroid%angularMomentum() < (thisSpheroid%massGas()+thisSpheroid%massStellar())*Dark_Matter_Halo_Virial_Radius(thisNode)*Dark_Matter_Halo_Virial_Velocity(thisNode) .and. thisSpheroid%radius() < Dark_Matter_Halo_Virial_Radius(thisNode)) then
          tidalField =Satellite_Tidal_Field(thisNode)
          tidalTorque=abs(tidalField)*(thisSpheroid%massGas()+thisSpheroid%massStellar())*thisSpheroid%radius()**2
          call thisSpheroid%angularMomentumRate(+tidalTorque)
       end if

    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Rate_Compute

  subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Rate(self,rate,interrupt,interruptProcedure)
    !% Adjust the rates for the star formation history.
    implicit none
    class    (nodeComponentSpheroidStandard), intent(inout)                    :: self
    type     (history                      ), intent(in   )                    :: rate
    logical                                 , intent(inout), optional          :: interrupt
    procedure(Interrupt_Procedure_Template ), intent(inout), optional, pointer :: interruptProcedure
    type     (history                      )                                   :: starFormationHistory

    ! Get the star formation history in the spheroid.
    starFormationHistory=self%starFormationHistory()
    ! Ensure that the history already exists.
    if (.not.starFormationHistory%exists())                                                                &
         & call Galacticus_Error_Report(                                                                   &
         &                              'Tree_Node_Spheroid_Star_Formation_History_Rate_Adjust_Standard', &
         &                              'no star formation history has been created in spheroid'           &
         &                             )
    ! Check if the star formation history in the spheroid spans a sufficient range to accept the input rates.
    if     (                                                                                             &
         &       rate%time(              1) < starFormationHistory%time(                              1) &
         &  .or. rate%time(size(rate%time)) > starFormationHistory%time(size(starFormationHistory%time)) &
         & ) then
       ! It does not, so interrupt evolution and extend the history.
       if (allocated(starFormationHistoryTemplate)) call Dealloc_Array(starFormationHistoryTemplate)
       call Alloc_Array(starFormationHistoryTemplate,shape(rate%time))
       starFormationHistoryTemplate=rate%time
       interrupt=.true.
       interruptProcedure => Node_Component_Spheroid_Standard_Star_Formation_History_Extend
       return
    end if
    ! Adjust the rate.
    starFormationHistory=starFormationHistory+rate
    call self%starFormationHistorySet(starFormationHistory)
    return
  end subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Rate

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Spheroid_Standard_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}. Note that gas masses get an additional scaling down since they can approach
    !% zero and we'd like to prevent them from becoming negative.
    use Abundances_Structure
    use Galacticus_Output_Star_Formation_Histories
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroidComponent
    class           (nodeComponentDisk    )               , pointer :: thisDiskComponent
    double precision                       , parameter              :: massMinimum                   =1.0d0
    double precision                       , parameter              :: angularMomentumMinimum        =0.1d0
    double precision                       , parameter              :: gasMassScaling                =0.1d0
    double precision                       , parameter              :: luminosityMinimum             =1.0d0
    double precision                                                :: angularMomentum                     , mass
    type            (history              )                         :: stellarPopulationHistoryScales

    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()
    ! Check if a standard spheroid component exists.
    select type (thisSpheroidComponent)
    class is (nodeComponentSpheroidStandard)

       ! Get the disk component.
       thisDiskComponent => thisNode%disk()

       ! Set scale for angular momentum.
       angularMomentum=thisDiskComponent%angularMomentum()+thisSpheroidComponent%angularMomentum()
       call thisSpheroidComponent%angularMomentumScale(               max(angularMomentum,angularMomentumMinimum))

       ! Set scale for gas mass.
       mass           =thisDiskComponent%massGas        ()+thisSpheroidComponent%massGas        ()
       call thisSpheroidComponent%        massGasScale(gasMassScaling*max(           mass,           massMinimum))

       ! Set scale for stellar mass.
       mass           =thisDiskComponent%massStellar    ()+thisSpheroidComponent%massStellar    ()
       call thisSpheroidComponent%    massStellarScale(               max(           mass,           massMinimum))

       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for gas abundances.
          call thisSpheroidComponent%abundancesGasScale    (                                                  &
               &                                             gasMassScaling                                   &
               &                                            *max(                                             &
               &                                                      thisDiskComponent%abundancesGas     ()  &
               &                                                 +thisSpheroidComponent%abundancesGas     (), &
               &                                                  massMinimum*unitAbundances                  &
               &                                                )                                             &
               &                                           )

          ! Set scale for stellar abundances.
          call thisSpheroidComponent%abundancesStellarScale(                                                  &
               &                                            max(                                              &
               &                                                       thisDiskComponent%abundancesStellar()  &
               &                                                +  thisSpheroidComponent%abundancesStellar(), &
               &                                                massMinimum*unitAbundances                    &
               &                                               )                                              &
               &                                           )
       end if

       ! Set scales for stellar luminosities if necessary.
       call thisSpheroidComponent%luminositiesStellarScale(                                                   &
            &                                               max(                                              &
            &                                                        thisDiskComponent%luminositiesStellar()  &
            &                                                   +thisSpheroidComponent%luminositiesStellar(), &
            &                                                    unitStellarLuminosities                      &
            &                                                   *luminosityMinimum                            &
            &                                                  )                                              &
            &                                              )

       ! Set scales for stellar population properties history.
       stellarPopulationHistoryScales=thisSpheroidComponent%stellarPropertiesHistory()
       call Stellar_Population_Properties_Scales               (stellarPopulationHistoryScales,thisSpheroidComponent%massStellar(),thisSpheroidComponent%abundancesStellar())
       call thisSpheroidComponent%stellarPropertiesHistoryScale(stellarPopulationHistoryScales                                                                      )
       call stellarPopulationHistoryScales%destroy()
       stellarPopulationHistoryScales=thisSpheroidComponent%starFormationHistory()
       call Star_Formation_History_Scales                      (stellarPopulationHistoryScales,thisSpheroidComponent%massStellar(),thisSpheroidComponent%abundancesStellar())
       call thisSpheroidComponent%starFormationHistoryScale    (stellarPopulationHistoryScales                                                                      )
       call stellarPopulationHistoryScales%destroy()
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Scale_Set

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !#  <after>Satellite_Merging_Remnant_Size</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Spheroid_Standard_Satellite_Merging(thisNode)
    !% Transfer any standard spheroid associated with {\tt thisNode} to its host halo.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    use Satellite_Merging_Remnant_Sizes_Properties
    use Abundances_Structure
    use Stellar_Luminosities_Structure
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    type            (treeNode             )               , pointer :: hostNode
    class           (nodeComponentDisk    )               , pointer :: hostDiskComponent    , thisDiskComponent
    class           (nodeComponentSpheroid)               , pointer :: hostSpheroidComponent, thisSpheroidComponent
    type            (history              )                         :: historyDisk          , historySpheroid                , &
         &                                                             thisHistory
    double precision                                                :: angularMomentum      , diskSpecificAngularMomentum    , &
         &                                                             spheroidMass         , spheroidSpecificAngularMomentum

    ! Check that the standard spheroid is active.
    if (defaultSpheroidComponent%standardIsActive()) then

       ! Get the spheroid component, creating it if need be.
       thisSpheroidComponent => thisNode%spheroid(autoCreate=.true.)
       select type (thisSpheroidComponent)
       class is (nodeComponentSpheroidStandard)
          thisDiskComponent => thisNode%disk()

          ! Find the node to merge with.
          hostNode              => thisNode%mergesWith(                 )
          hostDiskComponent     => hostNode%disk      (autoCreate=.true.)
          hostSpheroidComponent => hostNode%spheroid  (autoCreate=.true.)

          ! Get specific angular momentum of the host spheroid and disk material.
          if (hostSpheroidComponent%massGas()+hostSpheroidComponent%massStellar() > 0.0d0) then
             spheroidSpecificAngularMomentum=   hostSpheroidComponent%angularMomentum() &
                  &                          /(                                         &
                  &                             hostSpheroidComponent%massGas        () &
                  &                            +hostSpheroidComponent%massStellar    () &
                  &                           )
          else
             spheroidSpecificAngularMomentum=0.0d0
          end if
          if (hostDiskComponent%massGas()+hostDiskComponent%massStellar() > 0.0d0) then
             diskSpecificAngularMomentum=   hostDiskComponent%angularMomentum() &
                  &                      /(                                     &
                  &                         hostDiskComponent%massGas        () &
                  &                        +hostDiskComponent%massStellar    () &
                  &                       )
          else
             diskSpecificAngularMomentum=0.0d0
          end if

          ! Move gas material within the host if necessary.
          select case (thisHostGasMovesTo)
          case (movesToDisk)
             call hostDiskComponent    %        massGasSet(                                         &
                  &                                         hostDiskComponent    %massGas        () &
                  &                                        +hostSpheroidComponent%massGas        () &
                  &                                       )
             call hostDiskComponent    %  abundancesGasSet(                                         &
                  &                                         hostDiskComponent    %abundancesGas  () &
                  &                                        +hostSpheroidComponent%abundancesGas  () &
                  &                                       )
             call hostDiskComponent    %angularMomentumSet(                                         &
                  &                                         hostDiskComponent    %angularMomentum() &
                  &                                        +hostSpheroidComponent%massGas        () &
                  &                                        *spheroidSpecificAngularMomentum         &
                  &                                       )
             call hostSpheroidComponent%angularMomentumSet(                                         &
                  &                                         hostSpheroidComponent%angularMomentum() &
                  &                                        -hostSpheroidComponent%massGas        () &
                  &                                        *spheroidSpecificAngularMomentum         &
                  &                                       )
             call hostSpheroidComponent%        massGasSet(                                         &
                  &                                         0.0d0                                   &
                  &    )
             call hostSpheroidComponent%  abundancesGasSet(                                         &
                  &                                         zeroAbundances                          &
                  &                                       )
          case (movesToSpheroid)
             call hostSpheroidComponent%        massGasSet(&
                  &                                         hostSpheroidComponent%massGas        () &
                  &                                        +hostDiskComponent    %massGas        () &
                  &                                       )
             call hostDiskComponent    %angularMomentumSet(                                         &
                  &                                         hostDiskComponent%angularMomentum    () &
                  &                                        -hostDiskComponent%massGas            () &
                  &                                        *diskSpecificAngularMomentum             &
                  &                                       )
             call hostSpheroidComponent%  abundancesGasSet(                                         &
                  &                                         hostSpheroidComponent%abundancesGas  () &
                  &                                        +hostDiskComponent    %abundancesGas  () &
                  &                                       )
             call hostDiskComponent    %        massGasSet(                                         &
                  &                                         0.0d0                                   &
                  &                                       )
             call hostDiskComponent    %  abundancesGasSet(                                         &
                  &                                         zeroAbundances                          &
                  &                                       )
          case (doesNotMove)
             ! Do nothing.
          case default
             call Galacticus_Error_Report('Node_Component_Spheroid_Standard_Satellite_Merging','unrecognized movesTo descriptor')
          end select

          ! Move stellar material within the host if necessary.
          select case (thisHostStarsMoveTo)
          case (movesToDisk)
             call hostDiskComponent    %        massStellarSet     (                                             &
                  &                                                  hostDiskComponent    %        massStellar() &
                  &                                                 +hostSpheroidComponent%        massStellar() &
                  &                                                )
             call hostDiskComponent    %  abundancesStellarSet     (                                             &
                  &                                                  hostDiskComponent    %  abundancesStellar() &
                  &                                                 +hostSpheroidComponent%  abundancesStellar() &
                  &                                                )
             call hostDiskComponent    %luminositiesStellarSet     (                                             &
                  &                                                  hostDiskComponent    %luminositiesStellar() &
                  &                                                 +hostSpheroidComponent%luminositiesStellar() &
                  &                                                )
             call hostDiskComponent    %    angularMomentumSet     (                                             &
                  &                                                  hostDiskComponent    %    angularMomentum() &
                  &                                                 +hostSpheroidComponent%        massStellar() &
                  &                                                 *spheroidSpecificAngularMomentum             &
                  &                                                )
             call hostSpheroidComponent%    angularMomentumSet     (                                             &
                  &                                                  hostSpheroidComponent%    angularMomentum() &
                  &                                                 -hostSpheroidComponent%        massStellar() &
                  &                                                 *spheroidSpecificAngularMomentum             &
                  &                                                )
             call hostSpheroidComponent%        massStellarSet     (                                             &
                  &                                                  0.0d0                                       &
                  &                                                )
             call hostSpheroidComponent%  abundancesStellarSet     (                                             &
                  &                                                  zeroAbundances                              &
                  &                                                )
             call hostSpheroidComponent%luminositiesStellarSet     (                                             &
                  &                                                  zeroStellarLuminosities                     &
                  &                                                )
             ! Also add stellar properties histories.
             historyDisk    =    hostDiskComponent%stellarPropertiesHistory()
             historySpheroid=hostSpheroidComponent%stellarPropertiesHistory()
             call historyDisk    %increment(historySpheroid    )
             call historySpheroid%reset    (                   )
             call hostDiskComponent    %stellarPropertiesHistorySet(historyDisk    )
             call hostSpheroidComponent%stellarPropertiesHistorySet(historySpheroid)
             ! Also add star formation histories.
             historyDisk    =    hostDiskComponent%starFormationHistory()
             historySpheroid=hostSpheroidComponent%starFormationHistory()
             call historyDisk    %combine(historySpheroid     )
             call historySpheroid%reset  (                    )
             call hostDiskComponent    %    starFormationHistorySet(historyDisk    )
             call hostSpheroidComponent%    starFormationHistorySet(historySpheroid)
             call historyDisk    %destroy(recordMemory=.false.)
             call historySpheroid%destroy(recordMemory=.false.)
          case (movesToSpheroid)
             call hostSpheroidComponent%             massStellarSet(                                             &
                  &                                                  hostSpheroidComponent%        massStellar() &
                  &                                                 +hostDiskComponent    %        massStellar() &
                  &                                                )
             call hostDiskComponent    %         angularMomentumSet( hostDiskComponent    %    angularMomentum() &
                  &                                                 -hostDiskComponent    %        massStellar() &
                  &                                                 *diskSpecificAngularMomentum                 &
                  &                                                )
             call hostSpheroidComponent%       abundancesStellarSet(                                             &
                  &                                                  hostSpheroidComponent%  abundancesStellar() &
                  &                                                 +hostDiskComponent    %  abundancesStellar() &
                  &                                                )
             call hostSpheroidComponent%     luminositiesStellarSet( hostSpheroidComponent%luminositiesStellar() &
                  &                                                 +hostDiskComponent    %luminositiesStellar() &
                  &                                                )
             call hostDiskComponent    %             massStellarSet(                                             &
                  &                                                  0.0d0                                       &
                  &                                                )
             call hostDiskComponent    %       abundancesStellarSet(                                             &
                  &                                                  zeroAbundances                              &
                  &                                                )
             call hostDiskComponent    %     luminositiesStellarSet(                                             &
                  &                                                  zeroStellarLuminosities                     &
                  &                                                )
             ! Also add stellar properties histories.
             historyDisk    =    hostDiskComponent%stellarPropertiesHistory()
             historySpheroid=hostSpheroidComponent%stellarPropertiesHistory()
             call historySpheroid%increment(historyDisk)
             call historyDisk    %reset    (           )
             call hostSpheroidComponent%stellarPropertiesHistorySet(historySpheroid)
             call hostDiskComponent    %stellarPropertiesHistorySet( historyDisk   )
             ! Also add star formation histories.
             historyDisk    =hostDiskComponent    %starFormationHistory()
             historySpheroid=hostSpheroidComponent%starFormationHistory()
             call historySpheroid%combine(historyDisk         )
             call historyDisk    %reset  (                    )
             call hostSpheroidComponent%starFormationHistorySet    (historySpheroid)
             call hostDiskComponent    %starFormationHistorySet    (historyDisk    )
             call historyDisk    %destroy(recordMemory=.false.)
             call historySpheroid%destroy(recordMemory=.false.)
             historyDisk    =hostDiskComponent    %starFormationHistory()
             historySpheroid=hostSpheroidComponent%starFormationHistory()
          case (doesNotMove)
             ! Do nothing.
          case default
             call Galacticus_Error_Report('Node_Component_Spheroid_Standard_Satellite_Merging','unrecognized movesTo descriptor')
          end select

          ! If the entire host disk/spheroid (gas plus stars) was moved to the spheroid/disk, ensure that the
          ! corresponding angular momentum is precisely zero.
          if (hostDiskComponent    %massStellar()+hostDiskComponent    %massGas() == 0.0d0) &
               & call hostDiskComponent%angularMomentumSet    (0.0d0)
          if (hostSpheroidComponent%massStellar()+hostSpheroidComponent%massGas() == 0.0d0) &
               & call hostSpheroidComponent%angularMomentumSet(0.0d0)

          ! Get specific angular momentum of the spheroid material.
          spheroidMass=thisSpheroidComponent%massGas()+thisSpheroidComponent%massStellar()
          if (spheroidMass > 0.0d0) then
             spheroidSpecificAngularMomentum=thisSpheroidComponent%angularMomentum()/spheroidMass

             ! Move the gas component of the standard spheroid to the host.
             select case (thisMergerGasMovesTo)
             case (movesToDisk)
                call hostDiskComponent    %            massGasSet(                                             &
                     &                                             hostDiskComponent    %            massGas() &
                     &                                            +thisSpheroidComponent%            massGas() &
                     &                                           )
                call hostDiskComponent    %      abundancesGasSet(                                             &
                     &                                             hostDiskComponent    %      abundancesGas() &
                     &                                            +thisSpheroidComponent%      abundancesGas() &
                     &                                           )
                call hostDiskComponent    %    angularMomentumSet( hostDiskComponent    %    angularMomentum() &
                     &                                            +thisSpheroidComponent%            massGas() &
                     &                                            *spheroidSpecificAngularMomentum             &
                     &                                           )
             case (movesToSpheroid)
                call hostSpheroidComponent%            massGasSet(&
                     &                                             hostSpheroidComponent%            massGas() &
                     &                                            +thisSpheroidComponent%            massGas() &
                     &                                           )
                call hostSpheroidComponent%      abundancesGasSet(                                             &
                     &                                             hostSpheroidComponent%      abundancesGas() &
                     &                                            +thisSpheroidComponent%      abundancesGas() &
                     &                                           )
             case default
                call Galacticus_Error_Report('Node_Component_Spheroid_Standard_Satellite_Merging','unrecognized movesTo descriptor')
             end select
             call thisSpheroidComponent%      massGasSet(0.0d0         )
             call thisSpheroidComponent%abundancesGasSet(zeroAbundances)

             ! Move the stellar component of the standard spheroid to the host.
             select case (thisMergerStarsMoveTo)
             case (movesToDisk)
                call hostDiskComponent    %        massStellarSet(                                             &
                     &                                             hostDiskComponent    %        massStellar() &
                     &                                            +thisSpheroidComponent%        massStellar() &
                     &                                           )
                call hostDiskComponent    %  abundancesStellarSet( hostDiskComponent    %  abundancesStellar() &
                     &                                            +thisSpheroidComponent%  abundancesStellar() &
                     &                                           )
                call hostDiskComponent    %luminositiesStellarSet( hostDiskComponent    %luminositiesStellar() &
                     &                                            +thisSpheroidComponent%luminositiesStellar() &
                     &                                           )
                call hostDiskComponent    %    angularMomentumSet( hostDiskComponent    %    angularMomentum() &
                     &                                            +thisSpheroidComponent%        massStellar() &
                     &                                            *spheroidSpecificAngularMomentum             &
                     &                                           )
                ! Also add stellar properties histories.
                historySpheroid=thisSpheroidComponent%stellarPropertiesHistory()
                thisHistory    =hostDiskComponent    %stellarPropertiesHistory()
                call thisHistory    %increment(historySpheroid    )
                call historySpheroid%reset    (                   )
                call hostDiskComponent    %stellarPropertiesHistorySet(thisHistory    )
                call thisSpheroidComponent%stellarPropertiesHistorySet(historySpheroid)
                ! Also add star formation histories.
                historySpheroid=thisSpheroidComponent%starFormationHistory()
                thisHistory    =hostDiskComponent    %starFormationHistory()
                call thisHistory    %combine (historySpheroid     )
                call historySpheroid%reset   (                    )
                call hostDiskComponent    %starFormationHistorySet(thisHistory    )
                call thisSpheroidComponent%starFormationHistorySet(historySpheroid)
                call thisHistory    %destroy (recordMemory=.false.)
                call historySpheroid%destroy (recordMemory=.false.)
             case (movesToSpheroid)
                call hostSpheroidComponent%        massStellarSet( hostSpheroidComponent%        massStellar() &
                     &                                            +thisSpheroidComponent%        massStellar() &
                     &                                           )
                call hostSpheroidComponent%  abundancesStellarSet( hostSpheroidComponent%  abundancesStellar() &
                     &                                            +thisSpheroidComponent%  abundancesStellar() &
                     &                                           )
                call hostSpheroidComponent%luminositiesStellarSet( hostSpheroidComponent%luminositiesStellar() &
                     &                                            +thisSpheroidComponent%luminositiesStellar() &
                     &                                           )
                ! Also add stellar properties histories.
                historySpheroid=thisSpheroidComponent%stellarPropertiesHistory()
                thisHistory    =hostSpheroidComponent%stellarPropertiesHistory()
                call thisHistory%increment(historySpheroid)
                call historySpheroid%reset()
                call hostSpheroidComponent%stellarPropertiesHistorySet(thisHistory    )
                call thisSpheroidComponent%stellarPropertiesHistorySet(historySpheroid)
                ! Also add star formation histories.
                historySpheroid=thisSpheroidComponent%starFormationHistory()
                thisHistory    =hostSpheroidComponent%starFormationHistory()
                call thisHistory%combine(historySpheroid)
                call historySpheroid%reset()
                call hostSpheroidComponent%starFormationHistorySet(thisHistory    )
                call thisSpheroidComponent%starFormationHistorySet(historySpheroid)
                call thisHistory    %destroy(recordMemory=.false.)
                call historySpheroid%destroy(recordMemory=.false.)
             case default
                call Galacticus_Error_Report('Node_Component_Spheroid_Standard_Satellite_Merging','unrecognized movesTo descriptor')
             end select
             call thisSpheroidComponent%        massStellarSet(0.0d0                  )
             call thisSpheroidComponent%  abundancesStellarSet(zeroAbundances         )
             call thisSpheroidComponent%luminositiesStellarSet(zeroStellarLuminosities)
             call thisSpheroidComponent%    angularMomentumSet(0.0d0                  )
          end if

          ! Set the angular momentum of the spheroid.
          if (remnantSpecificAngularMomentum /= remnantNoChangeValue) then
             ! Note that the remnant specific angular momentum computed by the merger remnant modules automatically gives the mean
             ! specific angular momentum of the component by virtue of the fact that it computes the ratio of the actual angular
             ! momentum to the contribution from the component's own rotation curve at its scale radius.
             angularMomentum=remnantSpecificAngularMomentum*(hostSpheroidComponent%massGas()+hostSpheroidComponent%massStellar())
             call hostSpheroidComponent%angularMomentumSet(angularMomentum)
          end if

       end select
    end if
    return
  end subroutine Node_Component_Spheroid_Standard_Satellite_Merging

  !# <radiusSolverPlausibility>
  !#  <unitName>Node_Component_Spheroid_Standard_Radius_Solver_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Node_Component_Spheroid_Standard_Radius_Solver_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
    !% Determines whether the spheroid is physically plausible for radius solving tasks. Require that it have non-zero mass and angular momentum.
    implicit none
    type   (treeNode             ), intent(inout), pointer :: thisNode
    logical                       , intent(inout)          :: galaxyIsPhysicallyPlausible
    class  (nodeComponentSpheroid)               , pointer :: thisSpheroidComponent

    ! Return immediately if our method is not selected.
    if (.not.defaultSpheroidComponent%standardIsActive()) return

    ! Determine the plausibility of the current spheroid.
    thisSpheroidComponent => thisNode%spheroid()
    select type (thisSpheroidComponent)
       class is (nodeComponentSpheroidStandard)
          ! Determine the plausibility of the current spheroid.
       if        (thisSpheroidComponent%massStellar          ()+thisSpheroidComponent%massGas() < -spheroidMassToleranceAbsolute) then
          galaxyIsPhysicallyPlausible=.false.
       else
          if     (      thisSpheroidComponent%massStellar    ()+thisSpheroidComponent%massGas() >  spheroidMassToleranceAbsolute &
               &  .and. thisSpheroidComponent%angularMomentum()                                 <                          0.0d0 &
               & ) galaxyIsPhysicallyPlausible=.false.
       end if
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Radius_Solver_Plausibility

  double precision function Node_Component_Spheroid_Standard_Radius_Solve(thisNode)
    !% Return the circular radius of the standard spheroid.
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentSpheroid)               , pointer :: thisSpheroidComponent

    thisSpheroidComponent => thisNode%spheroid()
    Node_Component_Spheroid_Standard_Radius_Solve=thisSpheroidComponent%radius()
    return
  end function Node_Component_Spheroid_Standard_Radius_Solve

  double precision function Node_Component_Spheroid_Standard_Velocity_Solve(thisNode)
    !% Return the circular velocity of the standard spheroid.
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentSpheroid)               , pointer :: thisSpheroidComponent

    thisSpheroidComponent => thisNode%spheroid()
    Node_Component_Spheroid_Standard_Velocity_Solve=thisSpheroidComponent%velocity()
    return
  end function Node_Component_Spheroid_Standard_Velocity_Solve

  subroutine Node_Component_Spheroid_Standard_Radius_Solve_Set(thisNode,radius)
    !% Set the scale radius of the standard spheroid.
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    double precision                       , intent(in   )          :: radius
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroidComponent

    thisSpheroidComponent => thisNode%spheroid()
    call thisSpheroidComponent%radiusSet(max(radius,0.0d0))
    return
  end subroutine Node_Component_Spheroid_Standard_Radius_Solve_Set

  subroutine Node_Component_Spheroid_Standard_Velocity_Solve_Set(thisNode,velocity)
    !% Set the scale velocity of the standard spheroid.
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    double precision                       , intent(in   )          :: velocity
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroidComponent

    thisSpheroidComponent => thisNode%spheroid()
    call thisSpheroidComponent%velocitySet(max(velocity,0.0d0))
    return
  end subroutine Node_Component_Spheroid_Standard_Velocity_Solve_Set

  !# <radiusSolverTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Radius_Solver</unitName>
  !# </radiusSolverTask>
  subroutine Node_Component_Spheroid_Standard_Radius_Solver(thisNode,componentActive,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get&
       &,Velocity_Set)
    !% Interface for the size solver algorithm.
    implicit none
    type            (treeNode                                                              ), intent(inout), pointer :: thisNode
    logical                                                                                 , intent(  out)          :: componentActive
    double precision                                                                        , intent(  out)          :: specificAngularMomentum
    procedure       (Node_Component_Spheroid_Standard_Radius_Solve_Set                     ), intent(  out), pointer :: Radius_Set             , Velocity_Set
    procedure       (Node_Component_Spheroid_Standard_Radius_Solve                         ), intent(  out), pointer :: Radius_Get             , Velocity_Get
    class           (nodeComponentSpheroid                                                 )               , pointer :: thisSpheroidComponent
    double precision                                                                                                 :: angularMomentum        , specificAngularMomentumMean, &
         &                                                                                                              spheroidMass

    ! Determine if thisNode has an active disk component supported by this module.
    componentActive=.false.
    thisSpheroidComponent => thisNode%spheroid()
    select type (thisSpheroidComponent)
       class is (nodeComponentSpheroidStandard)
       componentActive=.true.
       ! Get the angular momentum.
       angularMomentum=thisSpheroidComponent%angularMomentum()
       if (angularMomentum >= 0.0d0) then
          ! Compute the specific angular momentum at the scale radius, assuming a flat rotation curve.
          spheroidMass=thisSpheroidComponent%massGas()+thisSpheroidComponent%massStellar()
          if (spheroidMass > 0.0d0) then
             specificAngularMomentumMean=angularMomentum/spheroidMass
          else
             specificAngularMomentumMean=0.0d0
          end if
          specificAngularMomentum=spheroidAngularMomentumAtScaleRadius*specificAngularMomentumMean

          ! Associate the pointers with the appropriate property routines.
          Radius_Get   => Node_Component_Spheroid_Standard_Radius_Solve
          Radius_Set   => Node_Component_Spheroid_Standard_Radius_Solve_Set
          Velocity_Get => Node_Component_Spheroid_Standard_Velocity_Solve
          Velocity_Set => Node_Component_Spheroid_Standard_Velocity_Solve_Set
       else
          call Node_Component_Spheroid_Standard_Radius_Solve_Set  (thisNode,0.0d0)
          call Node_Component_Spheroid_Standard_Velocity_Solve_Set(thisNode,0.0d0)
          componentActive=.false.
       end if
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Radius_Solver

  subroutine Node_Component_Spheroid_Standard_Initializor(self)
    !% Initializes a standard spheroid component.
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type   (nodeComponentSpheroidStandard )          :: self
    type   (treeNode                      ), pointer :: selfNode
    type   (history                       )          :: starFormationHistory      , stellarPropertiesHistory
    logical                                          :: createStarFormationHistory, createStellarPropertiesHistory

    ! Return if already initialized.
    if (self%isInitialized()) return
    ! Get the associated node.
    selfNode => self%host()
    ! Determine which histories must be created.
    starFormationHistory          =self%starFormationHistory            ()
    createStarFormationHistory    =.not.starFormationHistory    %exists ()
    call                                starformationhistory    %destroy()
    stellarPropertiesHistory      =self%stellarPropertiesHistory        ()
    createStellarPropertiesHistory=.not.stellarPropertiesHistory%exists ()
    call                                stellarPropertiesHistory%destroy()
    ! Create the stellar properties history.
    if (createStellarPropertiesHistory) then
       call Stellar_Population_Properties_History_Create(selfNode,stellarPropertiesHistory)
       call self%stellarPropertiesHistorySet(                     stellarPropertiesHistory)
    end if
    ! Create the star formation history.
    if (createStarFormationHistory    ) then
       call Star_Formation_History_Create               (selfNode,    starFormationHistory)
       call self%    starFormationHistorySet(                         starFormationHistory)
    end if
    ! Record that the spheroid has been initialized.
    call self%isInitializedSet(.true.)
    return
  end subroutine Node_Component_Spheroid_Standard_Initializor

  double precision function Node_Component_Spheroid_Standard_Star_Formation_Rate(self)
    !% Return the star formation rate of the standard spheroid.
    use Star_Formation_Timescales_Spheroids
    implicit none
    class           (nodeComponentSpheroidStandard ), intent(inout) :: self
    type            (treeNode                      ), pointer       :: thisNode
    double precision                                                :: gasMass , starFormationTimescale

    ! Get the associated node.
    thisNode => self%host()

    ! Get the star formation timescale.
    starFormationTimescale=Star_Formation_Timescale_Spheroid(thisNode)

    ! Get the gas mass.
    gasMass=self%massGas()

    ! If timescale is finite and gas mass is positive, then compute star formation rate.
    if (starFormationTimescale > 0.0d0 .and. gasMass > 0.0d0) then
       Node_Component_Spheroid_Standard_Star_Formation_Rate=gasMass/starFormationTimescale
    else
       Node_Component_Spheroid_Standard_Star_Formation_Rate=0.0d0
    end if
    return
  end function Node_Component_Spheroid_Standard_Star_Formation_Rate

  subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Extend(thisNode)
    !% Extend the range of a star formation history in a standard spheroid component for {\tt thisNode}.
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentSpheroid)               , pointer :: thisSpheroidComponent
    type (history              )                         :: starFormationHistory

    ! Get the spheroid component.
    thisSpheroidComponent => thisNode%spheroid()

    ! Extend the range as necessary.
    starFormationHistory=thisSpheroidComponent%starFormationHistory()
    call starFormationHistory%extend(times=starFormationHistoryTemplate)
    call thisSpheroidComponent%starFormationHistorySet(starFormationHistory)
    return
  end subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Extend

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Spheroid_Standard_Star_Formation_History_Output</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Output(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Store the star formation history in the output file.
    use Kind_Numbers
    use Galacticus_Output_Star_Formation_Histories
    implicit none
    type   (treeNode             ), intent(inout), pointer :: thisNode
    integer                       , intent(in   )          :: iOutput
    integer(kind=kind_int8       ), intent(in   )          :: treeIndex
    logical                       , intent(in   )          :: nodePassesFilter
    class  (nodeComponentSpheroid)               , pointer :: thisSpheroidComponent
    type   (history              )                         :: starFormationHistory

    ! Output the star formation history if a spheroid exists for this component.
    thisSpheroidComponent => thisNode%spheroid()
    select type (thisSpheroidComponent)
    class is (nodeComponentSpheroidStandard)
       starFormationHistory=thisSpheroidComponent%starFormationHistory()
       call Star_Formation_History_Output(thisNode,nodePassesFilter,starFormationHistory,iOutput,treeIndex,'spheroid')
       call thisSpheroidComponent%starFormationHistorySet(starFormationHistory)
    end select
    return
  end subroutine Node_Component_Spheroid_Standard_Star_Formation_History_Output

end module Node_Component_Spheroid_Standard

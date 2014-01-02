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

!% Contains a module which implements calculations of properties of ADAFs based on the implementation of \cite{benson_maximum_2009}.

module Accretion_Disks_ADAF
  !% Implements calculations of properties of ADAFs based on the implementation of \cite{benson_maximum_2009}.
  use ISO_Varying_String
  use Tables
  implicit none
  private
  public :: Accretion_Disks_ADAF_Initialize, Accretion_Disk_Radiative_Efficiency_ADAF,&
       & Black_Hole_Spin_Up_Rate_ADAF, Accretion_Disk_Jet_Power_ADAF

  ! Flag indicating if the module has been initialized.
  logical                                               :: adafInitialized                    =.false.

  ! Option controlling type of radiative effiency to use.
  integer                                               :: adafRadiativeEfficiencyType
  integer                                   , parameter :: adafRadiativeEfficiencyTypeFixed   =0
  integer                                   , parameter :: adafRadiativeEfficiencyTypeThinDisk=1

  ! Radiative efficiency of the accretion flow.
  double precision                                      :: adafRadiativeEfficiency

  ! Adiabatic index of the accretion flow.
  double precision                                      :: adafAdiabaticIndex                         , adafThermalPressureFraction

  ! Limit to the jet efficiency.
  double precision                                      :: adafJetEfficiencyMaximum

  ! Options for the viscosity prescription.
  type            (varying_string          )            :: adafViscosityOption
  integer                                   , parameter :: adafViscosityFit                   =1      , adafViscosityFixed         =0
  integer                                               :: adafViscosity
  double precision                                      :: adafViscosityFixedAlpha

  ! Options for the field-enhancing shear.
  type            (varying_string          )            :: adafFieldEnhanceType
  integer                                   , parameter :: adafFieldEnhanceExponential        =0      , adafFieldEnhanceLinear     =1
  integer                                               :: adafFieldEnhance

  ! Variable determining whether ADAF energy is 1 or E_IS
  type            (varying_string          )            :: adafEnergyOption
  integer                                   , parameter :: adafEnergy1                        =1      , adafEnergyIsco             =0
  integer                                               :: adafEnergy

  ! Tables to store spin-up and jet power functions.
  logical                                               :: adafTableTabulated                 =.false.
  integer                                   , parameter :: adafTableCount                     =10000
  integer                                   , parameter :: jetPowerTable                      =1      , spinUpTable                =2
  type            (table1DLogarithmicLinear)            :: adafTable

contains

  !# <accretionDisksMethod>
  !#  <unitName>Accretion_Disks_ADAF_Initialize</unitName>
  !# </accretionDisksMethod>
  subroutine Accretion_Disks_ADAF_Initialize(accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get&
       &,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    implicit none
    type     (varying_string                          ), intent(in   )          :: accretionDisksMethod
    procedure(Accretion_Disk_Radiative_Efficiency_ADAF), intent(inout), pointer :: Accretion_Disk_Radiative_Efficiency_Get
    procedure(Black_Hole_Spin_Up_Rate_ADAF            ), intent(inout), pointer :: Black_Hole_Spin_Up_Rate_Get
    procedure(Accretion_Disk_Jet_Power_ADAF           ), intent(inout), pointer :: Accretion_Disk_Jet_Power_Get

    if (accretionDisksMethod == 'ADAF') then
       Accretion_Disk_Radiative_Efficiency_Get => Accretion_Disk_Radiative_Efficiency_ADAF
       Black_Hole_Spin_Up_Rate_Get             => Black_Hole_Spin_Up_Rate_ADAF
       Accretion_Disk_Jet_Power_Get            => Accretion_Disk_Jet_Power_ADAF
       call Accretion_Disks_ADAF_Get_Parameters
    end if
    return
  end subroutine Accretion_Disks_ADAF_Initialize

  subroutine Accretion_Disks_ADAF_Get_Parameters
    !% Initialize the module by reading in parameter values.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    double precision                 :: adafAdiabaticIndexDefault
    type            (varying_string) :: adafRadiativeEfficiencyTypeText

    if (.not.adafInitialized) then
       !$omp critical(adafInitalize)
       if (.not.adafInitialized) then
          !@ <inputParameter>
          !@   <name>adafRadiativeEfficiencyType</name>
          !@   <defaultValue>thinDisk</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies the specific energy of material at the inner edge of an ADAF. {\tt pureADAF} makes the specific energy equal
          !@     to 1 (i.e. all energy is advected with the flow); {\tt ISCO} makes the specific energy equal to that for the innermost
          !@     stable circular orbit.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("adafRadiativeEfficiencyType",adafRadiativeEfficiencyTypeText,defaultValue="thinDisk")
          select case (char(adafRadiativeEfficiencyTypeText))
          case ("fixed")
             adafRadiativeEfficiencyType=adafRadiativeEfficiencyTypeFixed
          case ("thinDisk")
             adafRadiativeEfficiencyType=adafRadiativeEfficiencyTypeThinDisk
          case default
             call Galacticus_Error_Report('Accretion_Disks_ADAF_Initialize','unknown adafRadiativeEfficiencyType')
          end select
          !@ <inputParameter>
          !@   <name>adafRadiativeEfficiency</name>
          !@   <defaultValue>0.01</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies the radiative efficiency of an ADAF (i.e. the fraction of $\dot{M}\clight^2$ that is emitted in radiation).
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("adafRadiativeEfficiency",adafRadiativeEfficiency,defaultValue=0.01d0)
          !@ <inputParameter>
          !@   <name>adafEnergyOption</name>
          !@   <defaultValue>pureADAF</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies the specific energy of material at the inner edge of an ADAF. {\tt pureADAF} makes the specific energy equal
          !@     to 1 (i.e. all energy is advected with the flow); {\tt ISCO} makes the specific energy equal to that for the innermost
          !@     stable circular orbit.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("adafEnergyOption",adafEnergyOption,defaultValue="pureADAF")
          select case (char(adafEnergyOption))
          case ("pureADAF")
             adafEnergy=adafEnergy1
          case ("ISCO")
             adafEnergy=adafEnergyIsco
          case default
             call Galacticus_Error_Report('Accretion_Disks_ADAF_Initialize','unknown adafEnergyType')
          end select
          !@ <inputParameter>
          !@   <name>adafFieldEnhanceType</name>
          !@   <defaultValue>exponential</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Controls how the field enhancing shear is determined. {\tt exponential} will cause the form $g=\exp(\omega t)$ \citep{benson_maximum_2009}
          !@    to be used, while {\tt linear} will cause $g=1+\omega t$ to be used instead. The functional form of $\alpha(j)$ (if used) will be adjusted
          !@    to achieve a sensible spin-up function in each case.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("adafFieldEnhanceType",adafFieldEnhanceType,defaultValue="exponential")
          select case (char(adafFieldEnhanceType))
          case ("exponential")
             adafFieldEnhance         =adafFieldEnhanceExponential
             adafAdiabaticIndexDefault=1.444d0
          case ("linear")
             adafFieldEnhance         =adafFieldEnhanceLinear
             adafAdiabaticIndexDefault=1.333d0
          case default
             call Galacticus_Error_Report('Accretion_Disks_ADAF_Initialize','unknown adafFieldEnhanceType')
          end select
          !@ <inputParameter>
          !@   <name>adafAdiabaticIndex</name>
          !@   <defaultValue>1.444 (for exponential form of field-enhancing shear) or 1.333 (for linear form)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies the effective adiabatic index of gas in an ADAF.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("adafAdiabaticIndex",adafAdiabaticIndex,defaultValue=adafAdiabaticIndexDefault)
          adafThermalPressureFraction=(8.0d0-6.0d0*adafAdiabaticIndex)/3.0d0/(1.0d0-adafAdiabaticIndex)
          !@ <inputParameter>
          !@   <name>adafViscosityOption</name>
          !@   <defaultValue>fit</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Controls how the viscosity parameter $\alpha$ in an ADAF is determined. {\tt fit} will cause $\alpha$ to be computed
          !@    using the fitting function of \cite{benson_maximum_2009}; {\tt fixed} will cause $\alpha=${\tt [adafViscosityFixedAlpha]}
          !@    to be used.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("adafViscosityOption",adafViscosityOption,defaultValue="fit")
          select case (char(adafViscosityOption))
          case ("fixed")
             adafViscosity=adafViscosityFixed
             !@ <inputParameter>
             !@   <name>adafViscosityFixedAlpha</name>
             !@   <defaultValue>0.1</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@    The value for the viscosity parameter $\alpha$ in an ADAF to be used if {\tt [adafViscosityOption]}$=${\tt fixed}.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter("adafViscosityFixedAlpha",adafViscosityFixedAlpha,defaultValue=0.1d0)
          case ("fit")
             adafViscosity=adafViscosityFit
          case default
             call Galacticus_Error_Report('Accretion_Disks_ADAF_Initialize','unknown adafViscosityOption')
          end select
          !@ <inputParameter>
          !@   <name>adafJetEfficiencyMaximum</name>
          !@   <defaultValue>2</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The maximum efficiency allowed for ADAF-driven jets (in units of the accretion power).
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("adafJetEfficiencyMaximum",adafJetEfficiencyMaximum,defaultValue=2.0d0)
          adafInitialized=.true.
       end if
       !$omp end critical(adafInitalize)
    end if
    return
  end subroutine Accretion_Disks_ADAF_Get_Parameters

  double precision function Accretion_Disk_Radiative_Efficiency_ADAF(thisBlackHole,massAccretionRate)
    !% Computes the radiative efficiency for an ADAF.
    use Accretion_Disks_Shakura_Sunyaev
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate

    ! Ensure that parameters have been read.
    call Accretion_Disks_ADAF_Get_Parameters

    select case (adafRadiativeEfficiencyType)
    case (adafRadiativeEfficiencyTypeFixed   )
       Accretion_Disk_Radiative_Efficiency_ADAF=adafRadiativeEfficiency
    case (adafRadiativeEfficiencyTypeThinDisk)
       Accretion_Disk_Radiative_Efficiency_ADAF=Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev(thisBlackHole,massAccretionRate)
    end select
    return
  end function Accretion_Disk_Radiative_Efficiency_ADAF

  subroutine Accretion_Disk_ADAF_Tabulate()
    !% Tabulate jet power and spin-up efficiency for an ADAF.
    use Numerical_Constants_Physical
    use Black_Hole_Fundamentals
    implicit none
    double precision, parameter :: blackHoleSpinParameterMaximum=1.0d0, blackHoleSpinParameterMinimum=1.0d-6
    integer                     :: iSpin
    double precision            :: adafEnergyValue                    , adafViscosityAlpha                  , &
         &                         blackHoleSpin                      , radiusIsco                          , &
         &                         radiusStatic

    if (.not.adafTableTabulated) then
       !$omp critical(ADAF_Interpolate)
       if (.not.adafTableTabulated) then
          call adafTable%destroy()
          call adafTable%create (                                                    &
               &                 blackHoleSpinParameterMinimum                     , &
               &                 blackHoleSpinParameterMaximum                     , &
               &                 adafTableCount                                    , &
               &                 tableCount                   =2                   , &
               &                 extrapolationType            =extrapolationTypeFix  &
               &                )
          do iSpin=1,adafTableCount
             ! Get the black hole spin. The "spin parameter" that we tabulate in is 1-j so that we can easily pack
             ! many points close to j=1.
             blackHoleSpin=1.0d0-adafTable%x(iSpin)
             ! Determine the ADAF viscosity.
             select case (adafViscosity)
             case (adafViscosityFixed)
                adafViscosityAlpha=adafViscosityFixedAlpha
             case (adafViscosityFit)
                adafViscosityAlpha=ADAF_alpha(blackHoleSpin)
             end select
             ! Determine the ADAF energy.
             select case (adafEnergy)
             case (adafEnergy1)
                adafEnergyValue=1.0d0
             case (adafEnergyIsco)
                adafEnergyValue=Black_Hole_ISCO_Specific_Energy(blackHoleSpin,orbitPrograde)
             end select
             ! Compute jet launch radii.
             radiusIsco  =Black_Hole_ISCO_Radius  (blackHoleSpin)
             radiusStatic=Black_Hole_Static_Radius(blackHoleSpin)
             ! Compute the jet power.
             call adafTable%populate(                                                                           &
                  &                   min(                                                                      &
                  &                       (                                                                     &
                  &                         ADAF_BH_Jet_Power  (radiusStatic,blackHoleSpin,adafViscosityAlpha)  &
                  &                        +ADAF_Disk_Jet_Power(radiusIsco  ,blackHoleSpin,adafViscosityAlpha)  &
                  &                       ),                                                                    &
                  &                       adafJetEfficiencyMaximum                                              &
                  &                      )                                                                      &
                  &                   *(speedLight/kilo)**2                                                   , &
                  &                  iSpin                                                                    , &
                  &                  table=jetPowerTable                                                        &
                  &                 )
             ! Compute the rate of spin up to mass rate of change ratio.
             call adafTable%populate(                                                                                        &
                  &                   ADAF_Angular_Momentum                 (radiusIsco  ,blackHoleSpin,adafViscosityAlpha)  &
                  &                  -2.0d0*blackHoleSpin*adafEnergyValue                                                    &
                  &                  -Black_Hole_Rotational_Energy_Spin_Down(             blackHoleSpin                   )  &
                  &                  *(                                                                                      &
                  &                     ADAF_BH_Jet_Power                   (radiusStatic,blackHoleSpin,adafViscosityAlpha)  &
                  &                    +ADAF_Disk_Jet_Power_From_Black_Hole (radiusIsco  ,blackHoleSpin,adafViscosityAlpha)  &
                  &                   )                                                                                    , &
                  &                  iSpin                                                                                 , &
                  &                  table=spinUpTable                                                                       &
                  &                 )
          end do
          adafTableTabulated=.true.
       end if
       !$omp end critical(ADAF_Interpolate)
    end if
    return
  end subroutine Accretion_Disk_ADAF_Tabulate

  double precision function Accretion_Disk_Jet_Power_ADAF(thisBlackHole,massAccretionRate)
    !% Computes the jet power for an ADAF in units of $M_\odot$ (km/s)$^2$ Gyr$^{-1}$.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                                        :: blackHoleSpin    , blackHoleSpinParameter

    ! Ensure that parameters have been read.
    call Accretion_Disks_ADAF_Get_Parameters()

    ! Ensure tables have been constructed.
    call Accretion_Disk_ADAF_Tabulate()

    ! Get the black hole spin.
    blackHoleSpin=thisBlackHole%spin()

    ! Get the "spin parameter".
    blackHoleSpinParameter=1.0d0-blackHoleSpin

    ! Compute the jet power.
    !$omp critical(ADAF_Interpolate)
    Accretion_Disk_Jet_Power_ADAF=massAccretionRate*adafTable%interpolate(blackHoleSpinParameter,jetPowerTable)
    !$omp end critical(ADAF_Interpolate)
    return
  end function Accretion_Disk_Jet_Power_ADAF

  double precision function Black_Hole_Spin_Up_Rate_ADAF(thisBlackHole,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisBlackHole} due to accretion from an ADAF.
    !% disk.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                                        :: blackHoleSpin              , blackHoleSpinParameter, &
         &                                                     spinToMassRateOfChangeRatio

    ! Ensure that parameters have been read.
    call Accretion_Disks_ADAF_Get_Parameters

    ! Ensure tables have been constructed.
    call Accretion_Disk_ADAF_Tabulate()

    ! Get the black hole spin.
    blackHoleSpin=thisBlackHole%spin()

    ! Get the "spin parameter".
    blackHoleSpinParameter=1.0d0-blackHoleSpin

    ! Compute the ratio of spin and mass rates of change.
    !$omp critical(ADAF_Interpolate)
    spinToMassRateOfChangeRatio=adafTable%interpolate(blackHoleSpinParameter,spinUpTable)
    !$omp end critical(ADAF_Interpolate)

    ! Scale to the mass rate of change.
    Black_Hole_Spin_Up_Rate_ADAF=spinToMassRateOfChangeRatio*massAccretionRate/thisBlackHole%mass()
    return
  end function Black_Hole_Spin_Up_Rate_ADAF

  double precision function ADAF_Disk_Jet_Power_From_Black_Hole(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the power extracted from the black hole by the disk-launched jet from an ADAF.
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha        , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious, blackHoleSpinPrevious=2.0d0, &
         &                             jetPowerPrevious          , radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,jetPowerPrevious)
    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       jetPowerPrevious=ADAF_Disk_Jet_Power(radius,blackHoleSpin,adafViscosityAlpha)*(1.0d0-1.0d0&
            &/ADAF_Field_Enhancement(radius,blackHoleSpin,adafViscosityAlpha)**2)
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_Disk_Jet_Power_From_Black_Hole=jetPowerPrevious
    return
  end function ADAF_Disk_Jet_Power_From_Black_Hole

  double precision function ADAF_Disk_Jet_Power(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the power of the disk-launched jet from an ADAF.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha        , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious, blackHoleSpinPrevious=2.0d0, &
         &                             diskPowerPrevious         , radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,diskPowerPrevious)
    double precision                :: betaPhi

    ! Check if arguments are the same as on the previous call.
    if (radius == radiusPrevious .and. blackHoleSpin == blackHoleSpinPrevious .and. adafViscosityAlpha == adafViscosityAlphaPrevious) then
       ! They are, so return the previously computed value.
       ADAF_Disk_Jet_Power=diskPowerPrevious
    else
       ! They are not, so compute (and store) a new value.
       betaPhi=sqrt(1.0d0-1.0d0/ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)**2)
       ADAF_Disk_Jet_Power=(3.0d0/80.0d0)*radius**2                                                                                                &
            & *(2.0d0*blackHoleSpin*betaPhi/radius**2+sqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius)))**2                                  &
            & *(1.0d0-adafThermalPressureFraction)                                                                                                 &
            & *(ADAF_Field_Enhancement(radius,blackHoleSpin,adafViscosityAlpha)*ADAF_gamma(radius,blackHoleSpin,adafViscosityAlpha))**2            &
            & *sqrt((1.0d0-ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)**2)/Black_Hole_Metric_D_Factor(blackHoleSpin,radius))                  &
            & *(ADAF_Fluid_Angular_Velocity(radius,blackHoleSpin,adafViscosityAlpha)+Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius))**2 &
            & *ADAF_Temperature          (radius,blackHoleSpin,adafViscosityAlpha)                                                                 &
            & /Black_Hole_Metric_A_Factor(blackHoleSpin,radius)                                                                                    &
            & /ADAF_V                    (radius,blackHoleSpin,adafViscosityAlpha)                                                                 &
            & /ADAF_Height               (radius,blackHoleSpin,adafViscosityAlpha)
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
       diskPowerPrevious         =ADAF_Disk_Jet_Power
    end if
    return
  end function ADAF_Disk_Jet_Power

  double precision function ADAF_BH_Jet_Power(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the power of the black hole-launched jet from an ADAF.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha               , blackHoleSpin              , &
         &                             radius
    double precision, parameter     :: blackHoleSpinMinimum      =5.0d-8
    double precision, save          :: adafViscosityAlphaPrevious       , blackHoleSpinPrevious=2.0d0, &
         &                             jetPowerPrevious                 , radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,jetPowerPrevious)
    double precision                :: betaPhi

    if (blackHoleSpin > blackHoleSpinMinimum) then
       ! Check if arguments are the same as on the previous call.
       if (radius == radiusPrevious .and. blackHoleSpin == blackHoleSpinPrevious .and. adafViscosityAlpha == adafViscosityAlphaPrevious) then
          ! They are, so return the previously computed value.
          ADAF_BH_Jet_Power=jetPowerPrevious
       else
          ! They are not, so compute (and store) a new value.
          betaPhi=sqrt(1.0d0-1.0d0/ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)**2)
          ADAF_BH_Jet_Power=(3.0d0/80.0d0)*radius**2                                                                                       &
               & *(2.0d0*blackHoleSpin*betaPhi/radius**2+sqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius)))**2                       &
               & *(1.0d0-adafThermalPressureFraction)                                                                                      &
               & *(ADAF_Field_Enhancement(radius,blackHoleSpin,adafViscosityAlpha)*ADAF_gamma(radius,blackHoleSpin,adafViscosityAlpha))**2 &
               & *sqrt((1.0d0-ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)**2)/Black_Hole_Metric_D_Factor(blackHoleSpin,radius))       &
               & *Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius)**2                                                             &
               & *ADAF_Temperature          (radius,blackHoleSpin,adafViscosityAlpha)                                                      &
               & /Black_Hole_Metric_A_Factor(blackHoleSpin,radius)                                                                         &
               & /ADAF_V                    (radius,blackHoleSpin,adafViscosityAlpha)                                                      &
               & /ADAF_Height               (radius,blackHoleSpin,adafViscosityAlpha)
          radiusPrevious            =radius
          blackHoleSpinPrevious     =blackHoleSpin
          adafViscosityAlphaPrevious=adafViscosityAlpha
          jetPowerPrevious          =ADAF_BH_Jet_Power
       end if
    else
       ADAF_BH_Jet_Power=0.0d0
    end if
    return
  end function ADAF_BH_Jet_Power

  double precision function ADAF_Field_Enhancement(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the field enhancement factor, $g$, in the ADAF.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha        , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious, blackHoleSpinPrevious=2.0d0, &
         &                             fieldEnhancementPrevious  , radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,fieldEnhancementPrevious)
    double precision                :: tau                       , tauPhi                     , &
         &                             tauR

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /= adafViscosityAlphaPrevious) then
       tauPhi=1.0d0/ADAF_Fluid_Angular_Velocity(radius,blackHoleSpin,adafViscosityAlpha)
       tauR  =radius*ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)/sqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius))&
            &/ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)
       tau=min(tauPhi,tauR)
       select case (adafFieldEnhance)
       case (adafFieldEnhanceExponential)
          fieldEnhancementPrevious= exp(Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius)*tau)
       case (adafFieldEnhanceLinear     )
          fieldEnhancementPrevious=1.0d0+Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius)*tau
       end select
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_Field_Enhancement=fieldEnhancementPrevious
    return
  end function ADAF_Field_Enhancement

  double precision function ADAF_Fluid_Angular_Velocity(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the angular velocity of the rotating fluid with respect to the local inertial observer (ZAMO).
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha              , blackHoleSpin          , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious      , angularVelocityPrevious, &
         &                             blackHoleSpinPrevious     =2.0d0, radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,angularVelocityPrevious)
    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       angularVelocityPrevious=                                          &
            & ADAF_Angular_Momentum(radius,blackHoleSpin,adafViscosityAlpha) &
            & *sqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius)        &
            & /Black_Hole_Metric_A_Factor(blackHoleSpin,radius)**3)          &
            & /radius**2                                                     &
            & /ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)       &
            & /ADAF_gamma_r(radius,blackHoleSpin,adafViscosityAlpha)
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_Fluid_Angular_Velocity=angularVelocityPrevious
    return
  end function ADAF_Fluid_Angular_Velocity

  double precision function ADAF_alpha(blackHoleSpin)
    !% Returns the effective value of $\alpha$ for an ADAF.
    implicit none
    double precision, intent(in   ) :: blackHoleSpin
    double precision, save          :: alphaPrevious, blackHoleSpinPrevious=2.0d0
    !$omp threadprivate(blackHoleSpinPrevious,alphaPrevious)
    ! Check if we are being called with the same arguments as the previous call.
    if (blackHoleSpin /= blackHoleSpinPrevious) then
       select case (adafEnergy)
       case (adafEnergyISCO)
          select case (adafFieldEnhance)
          case (adafFieldEnhanceExponential)
             alphaPrevious=0.015d0+0.02d0*blackHoleSpin**4
          case (adafFieldEnhanceLinear     )
             alphaPrevious=0.025d0+0.08d0*blackHoleSpin**4
          end select
       case (adafEnergy1   )
          select case (adafFieldEnhance)
          case (adafFieldEnhanceExponential)
             alphaPrevious=0.010d0
          case (adafFieldEnhanceLinear     )
             alphaPrevious=0.025d0+0.02d0*blackHoleSpin**4
          end select
       end select
       blackHoleSpinPrevious=blackHoleSpin
    end if
    ADAF_alpha=alphaPrevious
    return
  end function ADAF_alpha

  double precision function ADAF_gamma(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the net relativistic boost factor from the fluid frame of an ADAF to an observer at rest at infinity.
    !% The input quantities are in natural units.
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha        , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious, blackHoleSpinPrevious=2.0d0, &
         &                             gammaPrevious             , radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,gammaPrevious)
    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       gammaPrevious=ADAF_gamma_r(radius,blackHoleSpin,adafViscosityAlpha)*ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_gamma=gammaPrevious
    return
  end function ADAF_gamma

  double precision function ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the $\phi$ component relativistic boost factor from the fluid frame of an ADAF to an observer at rest at infinity.
    !% The input quantities are in natural units.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha        , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious, blackHoleSpinPrevious=2.0d0, &
         &                             gammaPrevious             , radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,gammaPrevious)
    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       gammaPrevious=sqrt(1.0d0+((ADAF_Angular_Momentum(radius,blackHoleSpin,adafViscosityAlpha)/radius/ADAF_gamma_r(radius&
            &,blackHoleSpin,adafViscosityAlpha))**2)/Black_Hole_Metric_A_Factor(blackHoleSpin,radius))
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_gamma_phi=gammaPrevious
    return
  end function ADAF_gamma_phi

  double precision function ADAF_gamma_r(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the $r$ component relativistic boost factor from the fluid frame of an ADAF to an observer at rest at infinity.
    !% The input quantities are in natural units.
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha        , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious, blackHoleSpinPrevious=2.0d0, &
         &                             gammaPrevious             , radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,gammaPrevious)
    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       gammaPrevious=sqrt(1.0d0/(1.0d0-ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)**2))
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_gamma_r=gammaPrevious
    return
  end function ADAF_gamma_r

  double precision function ADAF_Angular_Momentum(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the specific angular momentum of accreted material in the ADAF.
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha              , blackHoleSpin          , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious      , angularMomentumPrevious, &
         &                             blackHoleSpinPrevious     =2.0d0, radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,angularMomentumPrevious)
    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       angularMomentumPrevious=ADAF_Enthalpy_Angular_Momentum_Product(radius,blackHoleSpin,adafViscosityAlpha)/ADAF_Enthalpy(radius&
            &,blackHoleSpin,adafViscosityAlpha)
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_Angular_Momentum=angularMomentumPrevious
    return
  end function ADAF_Angular_Momentum

  double precision function ADAF_Enthalpy_Angular_Momentum_Product(radius,blackHoleSpin,adafViscosityAlpha)
    !% Return the product of enthalpy and angular momentum for the ADAF.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha             , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious     , blackHoleSpinPrevious=2.0d0, &
         &                             enthalpyAngularMomentumPrevious, radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,enthalpyAngularMomentumPrevious)
    double precision                :: etaLADAF1                      , etaLADAF2                  , &
         &                             etaLADAF3                      , etaLADAF4                  , &
         &                             etaLADAF5                      , etaLADAF6                  , &
         &                             logarithmAlpha                 , radiusISCO

    ! Check if we are being called with the same arguments as the previous call.
    if (radius == radiusPrevious .and. blackHoleSpin == blackHoleSpinPrevious .and. adafViscosityAlpha == adafViscosityAlphaPrevious) then
       ! We are, so just return the stored value.
       ADAF_Enthalpy_Angular_Momentum_Product=enthalpyAngularMomentumPrevious
    else
       ! We are not, so compute (and store) the value.
       logarithmAlpha=log10(adafViscosityAlpha)
       radiusISCO=Black_Hole_ISCO_Radius(blackHoleSpin,unitsGravitational)
       etaLADAF1=0.0871d0*radiusISCO-0.10282d0
       etaLADAF2=0.5d0-7.7983d0*(adafAdiabaticIndex-1.333d0)**1.26d0
       etaLADAF3=0.153d0*(radiusISCO-0.6d0)**0.30d0+0.105d0
       etaLADAF4=etaLADAF3*(0.9d0*adafAdiabaticIndex-0.2996d0)*(1.202d0-0.08d0*(logarithmAlpha+2.5d0)**2.6d0)
       etaLADAF5=-1.8d0*adafAdiabaticIndex+4.299d0-0.018d0+0.018d0*(logarithmAlpha+2.0d0)**3.571d0
       etaLADAF6=etaLADAF4*(((0.14d0*log10(radius)**etaLADAF5+0.23d0)/etaLADAF4)**10.0d0+1.0d0)**0.1d0
       ADAF_Enthalpy_Angular_Momentum_Product=etaLADAF2+(etaLADAF1+10.0d0**etaLADAF6)*(1.15d0-0.03d0*(3.0d0+logarithmAlpha)**2.37d0)
       radiusPrevious                 =radius
       blackHoleSpinPrevious          =blackHoleSpin
       adafViscosityAlphaPrevious     =adafViscosityAlpha
       enthalpyAngularMomentumPrevious=ADAF_Enthalpy_Angular_Momentum_Product
    end if
    return
  end function ADAF_Enthalpy_Angular_Momentum_Product

  double precision function ADAF_Enthalpy(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the relativistic enthalpy of the ADAF.
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha        , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious, blackHoleSpinPrevious=2.0d0, &
         &                             enthalpyPrevious          , radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,enthalpyPrevious)
    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       enthalpyPrevious=1.0d0+(adafAdiabaticIndex/(adafAdiabaticIndex-1.0d0))*ADAF_Temperature(radius,blackHoleSpin,adafViscosityAlpha)
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_Enthalpy=enthalpyPrevious
    return
  end function ADAF_Enthalpy

  double precision function ADAF_Temperature(radius,blackHoleSpin,adafViscosityAlpha)
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha        , blackHoleSpin              , &
         &                             radius
    double precision, save          :: adafViscosityAlphaPrevious, blackHoleSpinPrevious=2.0d0, &
         &                             radiusPrevious            , temperaturePrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,temperaturePrevious)
    double precision                :: logarithmAlpha            , radiusISCO                 , &
         &                             t1                        , t2                         , &
         &                             t3                        , t4                         , &
         &                             t5

    ! Check if we are being called with the same arguments as the previous call.
    if (radius == radiusPrevious .and. blackHoleSpin == blackHoleSpinPrevious .and. adafViscosityAlpha == adafViscosityAlphaPrevious) then
       ! We are, so just return the stored value.
       ADAF_Temperature=temperaturePrevious
    else
       ! We are not, so compute (and store) the value.
       logarithmAlpha=log10(adafViscosityAlpha)
       radiusISCO=Black_Hole_ISCO_Radius(blackHoleSpin,unitsGravitational)
       t1=-0.270278d0*adafAdiabaticIndex+1.36027d0
       t2=-0.94d0+4.4744d0*(adafAdiabaticIndex-1.444d0)-5.1402d0*(adafAdiabaticIndex-1.444d0)**2
       t3=0.84d0*logarithmAlpha+0.919d0-0.643d0*exp(-0.209d0/adafViscosityAlpha)
       t4=(0.6365d0*radiusISCO-0.4828d0)*(1.0d0+11.9d0*exp(-0.838d0*radiusISCO**4))
       t5=1.444d0*exp(-1.01d0*radiusISCO**0.86d0)+0.1d0
       ADAF_Temperature=0.31d0*((1.0d0+(t4/radius)**0.9d0)**(t2+t3))/((radius-t5)**t1)
       radiusPrevious                 =radius
       blackHoleSpinPrevious          =blackHoleSpin
       adafViscosityAlphaPrevious     =adafViscosityAlpha
       temperaturePrevious=ADAF_Temperature
    end if
    return
  end function ADAF_Temperature

  double precision function ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)
    !% Return the (dimensionless) velocity in an ADAF at given {\tt radius}, for a black hole of given {\tt blackHoleSpin} and a
    !% flow with viscosity parameter {\tt adafViscosityAlpha}.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha         , blackHoleSpin             , &
         &                             radius
    double precision, save          :: adafVelocityPrevious       , adafViscosityAlphaPrevious, &
         &                             blackHoleSpinPrevious=2.0d0, radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,adafVelocityPrevious)
    double precision                :: Phi                        , alpha_eff                 , &
         &                             rISCO                      , reff                      , &
         &                             rh                         , v1                        , &
         &                             v2                         , v3                        , &
         &                             v4                         , v5                        , &
         &                             z                          , zh

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /= adafViscosityAlphaPrevious) then
       rh=Black_Hole_Horizon_Radius(blackHoleSpin)
       rISCO=Black_Hole_ISCO_Radius(blackHoleSpin)
       z=radius/rISCO
       zh=rh/rISCO
       alpha_eff=adafViscosityAlpha*(1.0d0+6.450d0*(adafAdiabaticIndex-1.444d0)+1.355d0*(adafAdiabaticIndex-1.444d0)**2)
       v1=9.0d0*log(9.0d0*z)
       v2=exp(-0.66d0*(1.0d0-2.0d0*alpha_eff)*log(alpha_eff/0.1d0)*log(z/zh))
       v3=1.0d0-exp(-z*(0.16d0*(blackHoleSpin-1.0d0)+0.76d0))
       v4=1.4d0+0.29065d0*(blackHoleSpin-0.5d0)**4-0.8756d0*(blackHoleSpin-0.5d0)**2+(-0.33d0*blackHoleSpin+0.45035d0)*(1.0d0-exp(-(z-zh)))
       v5=2.3d0*exp(40.0d0*(blackHoleSpin-1.0d0))*exp(-15.0d0*rISCO*(z-zh))+1.0d0
       Phi=v1*v2*v3*v4*v5
       reff=rh+Phi*(radius-rh)
       adafVelocityPrevious      =sqrt(1.0d0-(1.0d0-2.0d0/reff+(blackHoleSpin/reff)**2))
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_V=adafVelocityPrevious
    return
  end function ADAF_V

  double precision function ADAF_Height(radius,blackHoleSpin,adafViscosityAlpha)
    !% Return the (dimensionless) height in an ADAF at given {\tt radius}, for a black hole of given {\tt blackHoleSpin} and a
    !% flow with viscosity parameter {\tt adafViscosityAlpha}.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in   ) :: adafViscosityAlpha         , blackHoleSpin             , &
         &                             radius
    double precision, save          :: adafHeightPrevious         , adafViscosityAlphaPrevious, &
         &                             blackHoleSpinPrevious=2.0d0, radiusPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,adafHeightPrevious)
    double precision                :: nuz2

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       nuz2=blackHoleSpin**2+(1.0d0-(blackHoleSpin*Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius))**2) &
            &*ADAF_Angular_Momentum(radius,blackHoleSpin,adafViscosityAlpha)**2-((blackHoleSpin*ADAF_gamma_phi(radius&
            &,blackHoleSpin,adafViscosityAlpha))**2/Black_Hole_Metric_A_Factor(blackHoleSpin,radius))*ADAF_gamma_r(radius&
            &,blackHoleSpin,adafViscosityAlpha)**2*Black_Hole_Metric_D_Factor(blackHoleSpin,radius)-ADAF_gamma_r(radius&
            &,blackHoleSpin,adafViscosityAlpha)*sqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius)&
            &/Black_Hole_Metric_A_Factor(blackHoleSpin,radius))*2.0d0*ADAF_Angular_Momentum(radius,blackHoleSpin&
            &,adafViscosityAlpha)*Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius)*ADAF_gamma_phi(radius,blackHoleSpin&
            &,adafViscosityAlpha)*blackHoleSpin**2
       nuz2=nuz2/radius**4
       adafHeightPrevious=sqrt(ADAF_Temperature(radius,blackHoleSpin,adafViscosityAlpha)/ADAF_Enthalpy(radius,blackHoleSpin &
            &,adafViscosityAlpha)/radius**2/nuz2)
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_Height=adafHeightPrevious
    return
  end function ADAF_Height

end module Accretion_Disks_ADAF

!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements calculations of properties of ADAFs based on the implementation of \cite{benson_maximum_2009}.

module Accretion_Disks_ADAF
  !% Implements calculations of properties of ADAFs based on the implementation of \cite{benson_maximum_2009}.
  use ISO_Varying_String
  use Black_Hole_Fundamentals
  private
  public :: Accretion_Disks_ADAF_Initialize, Accretion_Disk_Radiative_Efficiency_ADAF,&
       & Black_Hole_Spin_Up_Rate_ADAF, Accretion_Disk_Jet_Power_ADAF

  ! Flag indicating if the module has been initialized.
  logical                        :: adafInitialized=.false.

  ! Radiative efficiency of the accretion flow.
  double precision               :: adafRadiativeEfficiency

  ! Adiabatic index of the accretion flow.
  double precision               :: adafAdiabaticIndex,adafThermalPressureFraction

  ! Options for the viscosity prescription.
  type(varying_string)           :: adafViscosityOption
  integer,             parameter :: adafViscosityFixed=0, adafViscosityFit=1
  integer                        :: adafViscosity
  double precision               :: adafViscosityFixedAlpha

  ! Variable determining whether ADAF energy is 1 or E_ISCO.
  type(varying_string)           :: adafEnergyOption
  integer,             parameter :: adafEnergyIsco=0, adafEnergy1=1
  integer                        :: adafEnergy

contains

  !# <accretionDisksMethod>
  !#  <unitName>Accretion_Disks_ADAF_Initialize</unitName>
  !# </accretionDisksMethod>
  subroutine Accretion_Disks_ADAF_Initialize(accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get&
       &,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: accretionDisksMethod
    procedure(),          pointer, intent(inout) :: Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get&
         &,Accretion_Disk_Jet_Power_Get

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

    if (.not.adafInitialized) then
       !@ <inputParameter>
       !@   <name>adafRadiativeEfficiency</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies the radiative efficiency of an ADAF (i.e. the fraction of $\dot{M}\clight^2$ that is emitted in radiation).
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("adafRadiativeEfficiency",adafRadiativeEfficiency,defaultValue=0.01d0)
       !@ <inputParameter>
       !@   <name>adafEnergyOption</name>
       !@   <defaultValue>pure ADAF</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies the specific energy of material at the inner edge of an ADAF. {\tt pure ADAF} makes the specific energy equal
       !@     to 1 (i.e. all energy is advected with the flow); {\tt ISCO} makes the specific energy equal to that for the innermost
       !@     stable circular orbit.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("adafEnergyOption",adafEnergyOption,defaultValue="pure ADAF")
       select case (char(adafEnergyOption))
       case ("pure ADAF")
          adafEnergy=adafEnergy1
       case ("ISCO")
          adafEnergy=adafEnergyIsco
       case default
          call Galacticus_Error_Report('Accretion_Disks_ADAF_Initialize','unknown adafEnergyType')
       end select
       !@ <inputParameter>
       !@   <name>adafAdiabaticIndex</name>
       !@   <defaultValue>1.444</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies the effective adiabatic index of gas in an ADAF.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("adafAdiabaticIndex",adafAdiabaticIndex,defaultValue=1.444d0)
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
          !@ </inputParameter>
          call Get_Input_Parameter("adafViscosityFixedAlpha",adafViscosityFixedAlpha,defaultValue=0.1d0)
       case ("fit")
          adafViscosity=adafViscosityFit
       case default
          call Galacticus_Error_Report('Accretion_Disks_ADAF_Initialize','unknown adafEnergyType')
       end select
       adafInitialized=.true.
    end if
    return
  end subroutine Accretion_Disks_ADAF_Get_Parameters

  double precision function Accretion_Disk_Radiative_Efficiency_ADAF(thisNode,massAccretionRate)
    !% Computes the radiative efficiency for an ADAF.
    use Black_Hole_Fundamentals
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    ! Ensure that parameters have been read.
    call Accretion_Disks_ADAF_Get_Parameters

    Accretion_Disk_Radiative_Efficiency_ADAF=adafRadiativeEfficiency
    return
  end function Accretion_Disk_Radiative_Efficiency_ADAF

  double precision function Accretion_Disk_Jet_Power_ADAF(thisNode,massAccretionRate)
    !% Computes the jet power for an ADAF in units of $M_\odot$ (km/s)$^2$ Gyr$^{-1}$.
    use Tree_Nodes
    use Tree_Node_Methods
    use Black_Hole_Fundamentals
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision                         :: radiusIsco,radiusStatic,blackHoleSpin,adafViscosityAlpha

    ! Ensure that parameters have been read.
    call Accretion_Disks_ADAF_Get_Parameters

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    ! Determine the ADAF viscosity.
    select case (adafViscosity)
    case (adafViscosityFixed)
       adafViscosityAlpha=adafViscosityFixedAlpha
    case (adafViscosityFit)
       adafViscosityAlpha=ADAF_alpha(blackHoleSpin)
    end select

    ! Compute jet launch radii.
    radiusIsco=Black_Hole_ISCO_Radius(thisNode,unitsGravitational)
    radiusStatic=Black_Hole_Static_Radius(thisNode,units=unitsGravitational)

    ! Compute the jet power.
    Accretion_Disk_Jet_Power_ADAF=(ADAF_BH_Jet_Power(radiusStatic,blackHoleSpin,adafViscosityAlpha) &
         &+ADAF_Disk_Jet_Power(radiusIsco,blackHoleSpin,adafViscosityAlpha))*massAccretionRate*(speedLight/kilo)**2

    return
  end function Accretion_Disk_Jet_Power_ADAF

  double precision function Black_Hole_Spin_Up_Rate_ADAF(thisNode,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisNode} due to accretion from an ADAF.
    !% disk.
    use Tree_Nodes
    use Tree_Node_Methods
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision                         :: radiusIsco,radiusStatic,blackHoleSpin,adafEnergyValue,adafViscosityAlpha&
         &,spinToMassRateOfChangeRatio

    ! Ensure that parameters have been read.
    call Accretion_Disks_ADAF_Get_Parameters

    ! Determine the ADAF energy.
    select case (adafEnergy)
    case (adafEnergy1)
       adafEnergyValue=1.0d0
    case (adafEnergyIsco)
       adafEnergyValue=Black_Hole_ISCO_Specific_Energy(thisNode,unitsGravitational,orbitPrograde)
    end select

    ! Get the black hole spin.
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)

    ! Determine the ADAF viscosity.
    select case (adafViscosity)
    case (adafViscosityFixed)
       adafViscosityAlpha=adafViscosityFixedAlpha
    case (adafViscosityFit)
       adafViscosityAlpha=ADAF_alpha(blackHoleSpin)
    end select

    ! Compute the spin up rate including braking term due to jets.
    radiusIsco=Black_Hole_ISCO_Radius(thisNode,unitsGravitational)
    radiusStatic=Black_Hole_Static_Radius(thisNode,units=unitsGravitational)

    spinToMassRateOfChangeRatio=ADAF_Angular_Momentum(radiusIsco,blackHoleSpin,adafViscosityAlpha)-2.0d0*blackHoleSpin&
         &*adafEnergyValue-Black_Hole_Rotational_Energy_Spin_Down(blackHoleSpin)*(ADAF_BH_Jet_Power(radiusStatic,blackHoleSpin&
         &,adafViscosityAlpha)+ADAF_Disk_Jet_Power_From_Black_Hole(radiusIsco,blackHoleSpin,adafViscosityAlpha))

    Black_Hole_Spin_Up_Rate_ADAF=spinToMassRateOfChangeRatio*massAccretionRate/Tree_Node_Black_Hole_Mass(thisNode)
    return
  end function Black_Hole_Spin_Up_Rate_ADAF

  double precision function ADAF_Disk_Jet_Power_From_Black_Hole(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the power extracted from the black hole by the disk-launched jet from an ADAF.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,jetPowerPrevious
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,diskPowerPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,diskPowerPrevious)
    double precision             :: betaPhi

    ! Check if arguments are the same as on the previous call.
    if (radius == radiusPrevious .and. blackHoleSpin == blackHoleSpinPrevious .and. adafViscosityAlpha == adafViscosityAlphaPrevious) then
       ! They are, so return the previously computed value.
       ADAF_Disk_Jet_Power=diskPowerPrevious
    else
       ! They are not, so compute (and store) a new value.
       betaPhi=dsqrt(1.0d0-1.0d0/ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)**2)
       ADAF_Disk_Jet_Power=(3.0d0/80.0d0)*radius**2                                                                                                &
            & *(2.0d0*blackHoleSpin*betaPhi/radius**2+dsqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius)))**2                                  &
            & *(1.0d0-adafThermalPressureFraction)                                                                                                 &
            & *(ADAF_Field_Enhancement(radius,blackHoleSpin,adafViscosityAlpha)*ADAF_gamma(radius,blackHoleSpin,adafViscosityAlpha))**2            &
            & *dsqrt((1.0d0-ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)**2)/Black_Hole_Metric_D_Factor(blackHoleSpin,radius))                  &
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, parameter  :: blackHoleSpinMinimum=5.0d-8
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,jetPowerPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,jetPowerPrevious)
    double precision             :: betaPhi
    
    if (blackHoleSpin > blackHoleSpinMinimum) then
       ! Check if arguments are the same as on the previous call.
       if (radius == radiusPrevious .and. blackHoleSpin == blackHoleSpinPrevious .and. adafViscosityAlpha == adafViscosityAlphaPrevious) then
          ! They are, so return the previously computed value.
          ADAF_BH_Jet_Power=jetPowerPrevious
       else
          ! They are not, so compute (and store) a new value.
          betaPhi=dsqrt(1.0d0-1.0d0/ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)**2)
          ADAF_BH_Jet_Power=(3.0d0/80.0d0)*radius**2                                                                                       &
               & *(2.0d0*blackHoleSpin*betaPhi/radius**2+dsqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius)))**2                       &
               & *(1.0d0-adafThermalPressureFraction)                                                                                      &
               & *(ADAF_Field_Enhancement(radius,blackHoleSpin,adafViscosityAlpha)*ADAF_gamma(radius,blackHoleSpin,adafViscosityAlpha))**2 &
               & *dsqrt((1.0d0-ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)**2)/Black_Hole_Metric_D_Factor(blackHoleSpin,radius))       &
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,fieldEnhancementPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,fieldEnhancementPrevious)
    double precision             :: tauPhi,tauR,tau

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /= adafViscosityAlphaPrevious) then
       tauPhi=1.0d0/ADAF_Fluid_Angular_Velocity(radius,blackHoleSpin,adafViscosityAlpha)
       tauR  =radius*ADAF_gamma_phi(radius,blackHoleSpin,adafViscosityAlpha)/dsqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius))&
            &/ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)
       tau=min(tauPhi,tauR)
       fieldEnhancementPrevious  =dexp(Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius)*tau)
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,angularVelocityPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,angularVelocityPrevious)

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       angularVelocityPrevious=                                          &
            & ADAF_Angular_Momentum(radius,blackHoleSpin,adafViscosityAlpha) & 
            & *dsqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius)        &
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
    double precision, intent(in) :: blackHoleSpin
    double precision, save       :: blackHoleSpinPrevious=2.0d0,alphaPrevious
    !$omp threadprivate(blackHoleSpinPrevious,alphaPrevious)

    ! Check if we are being called with the same arguments as the previous call.
    if (blackHoleSpin /= blackHoleSpinPrevious) then
       select case (adafEnergy)
       case (adafEnergyISCO)
          alphaPrevious=0.025d0+0.400d0*blackHoleSpin**4
       case (adafEnergy1)
          alphaPrevious=0.025d0+0.055d0*blackHoleSpin**2
       end select
       blackHoleSpinPrevious     =blackHoleSpin
    end if
    ADAF_alpha=alphaPrevious
    return
  end function ADAF_alpha

  double precision function ADAF_gamma(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the net relativistic boost factor from the fluid frame of an ADAF to an observer at rest at infinity.
    !% The input quantities are in natural units.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,gammaPrevious
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,gammaPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,gammaPrevious)

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       gammaPrevious=dsqrt(1.0d0+((ADAF_Angular_Momentum(radius,blackHoleSpin,adafViscosityAlpha)/radius/ADAF_gamma_r(radius&
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,gammaPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,gammaPrevious)

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       gammaPrevious=dsqrt(1.0d0/(1.0d0-ADAF_V(radius,blackHoleSpin,adafViscosityAlpha)**2))
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_gamma_r=gammaPrevious
    return
  end function ADAF_gamma_r

  double precision function ADAF_Angular_Momentum(radius,blackHoleSpin,adafViscosityAlpha)
    !% Returns the specific angular momentum of accreted material in the ADAF.
    use Black_Hole_Fundamentals
    implicit none
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,angularMomentumPrevious
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,enthalpyAngularMomentumPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,enthalpyAngularMomentumPrevious)
    double precision             :: logarithmAlpha,radiusISCO,etaLADAF1,etaLADAF2,etaLADAF3,etaLADAF4,etaLADAF5,etaLADAF6

    ! Check if we are being called with the same arguments as the previous call.
    if (radius == radiusPrevious .and. blackHoleSpin == blackHoleSpinPrevious .and. adafViscosityAlpha == adafViscosityAlphaPrevious) then
       ! We are, so just return the stored value.
       ADAF_Enthalpy_Angular_Momentum_Product=enthalpyAngularMomentumPrevious
    else
       ! We are not, so compute (and store) the value.
       logarithmAlpha=dlog10(adafViscosityAlpha)
       radiusISCO=Black_Hole_ISCO_Radius(blackHoleSpin,unitsGravitational)
       etaLADAF1=0.0871d0*radiusISCO-0.10282d0
       etaLADAF2=0.5d0-7.7983d0*(adafAdiabaticIndex-1.333d0)**1.26d0
       etaLADAF3=0.153d0*(radiusISCO-0.6d0)**0.30d0+0.105d0
       etaLADAF4=etaLADAF3*(0.9d0*adafAdiabaticIndex-0.2996d0)*(1.202d0-0.08d0*(logarithmAlpha+2.5d0)**2.6d0)
       etaLADAF5=-1.8d0*adafAdiabaticIndex+4.299d0-0.018d0+0.018d0*(logarithmAlpha+2.0d0)**3.571d0
       etaLADAF6=etaLADAF4*(((0.14*dlog10(radius)**etaLADAF5+0.23)/etaLADAF4)**10.0+1.0)**0.1
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,enthalpyPrevious
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,temperaturePrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,enthalpyAngularMomentumPrevious)
    double precision             :: logarithmAlpha,radiusISCO,t1,t2,t3,t4,t5

    ! Check if we are being called with the same arguments as the previous call.
    if (radius == radiusPrevious .and. blackHoleSpin == blackHoleSpinPrevious .and. adafViscosityAlpha == adafViscosityAlphaPrevious) then
       ! We are, so just return the stored value.
       ADAF_Temperature=temperaturePrevious
    else
       ! We are not, so compute (and store) the value.
       logarithmAlpha=dlog10(adafViscosityAlpha)
       radiusISCO=Black_Hole_ISCO_Radius(blackHoleSpin,unitsGravitational)
       t1=-0.270278d0*adafAdiabaticIndex+1.36027d0
       t2=-0.94d0+4.4744d0*(adafAdiabaticIndex-1.444d0)-5.1402d0*(adafAdiabaticIndex-1.444d0)**2
       t3=0.84d0*logarithmAlpha+0.919d0-0.643d0*dexp(-0.209d0/adafViscosityAlpha)
       t4=(0.6365d0*radiusISCO-0.4828d0)*(1.0d0+11.9d0*dexp(-0.838d0*radiusISCO**4))
       t5=1.444d0*dexp(-1.01d0*radiusISCO**0.86d0)+0.1d0
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,adafVelocityPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,adafVelocityPrevious)
    double precision             :: rh,rISCO,z,zh,alpha_eff,v1,v2,v3,v4,v5,Phi,reff

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /= adafViscosityAlphaPrevious) then
       rh=Black_Hole_Horizon_Radius(blackHoleSpin)
       rISCO=Black_Hole_ISCO_Radius(blackHoleSpin)
       z=radius/rISCO
       zh=rh/rISCO
       alpha_eff=adafViscosityAlpha*(1.0d0+6.450d0*(adafAdiabaticIndex-1.444d0)+1.355d0*(adafAdiabaticIndex-1.444d0)**2)
       v1=9.0d0*dlog(9.0d0*z)
       v2=dexp(-0.66d0*(1.0d0-2.0d0*alpha_eff)*dlog(alpha_eff/0.1d0)*dlog(z/zh))
       v3=1.0d0-dexp(-z*(0.16d0*(blackHoleSpin-1.0d0)+0.76d0))
       v4=1.4d0+0.29065d0*(blackHoleSpin-0.5d0)**4-0.8756d0*(blackHoleSpin-0.5d0)**2+(-0.33d0*blackHoleSpin+0.45035d0)*(1.0d0-dexp(-(z-zh)))
       v5=2.3d0*dexp(40.0d0*(blackHoleSpin-1.0d0))*dexp(-15.0d0*rISCO*(z-zh))+1.0d0
       Phi=v1*v2*v3*v4*v5
       reff=rh+Phi*(radius-rh)
       adafVelocityPrevious      =dsqrt(1.0d0-(1.0d0-2.0d0/reff+(blackHoleSpin/reff)**2))
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
    double precision, intent(in) :: radius,blackHoleSpin,adafViscosityAlpha
    double precision, save       :: radiusPrevious,blackHoleSpinPrevious=2.0d0,adafViscosityAlphaPrevious,adafHeightPrevious
    !$omp threadprivate(radiusPrevious,blackHoleSpinPrevious,adafViscosityAlphaPrevious,adafHeightPrevious)
    double precision             :: nuz2

    ! Check if we are being called with the same arguments as the previous call.
    if (radius /= radiusPrevious .or. blackHoleSpin /= blackHoleSpinPrevious .or. adafViscosityAlpha /=&
         & adafViscosityAlphaPrevious) then
       nuz2=blackHoleSpin**2+(1.0d0-(blackHoleSpin*Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius))**2) &
            &*ADAF_Angular_Momentum(radius,blackHoleSpin,adafViscosityAlpha)**2-((blackHoleSpin*ADAF_gamma_phi(radius&
            &,blackHoleSpin,adafViscosityAlpha))**2/Black_Hole_Metric_A_Factor(blackHoleSpin,radius))*ADAF_gamma_r(radius&
            &,blackHoleSpin,adafViscosityAlpha)**2*Black_Hole_Metric_D_Factor(blackHoleSpin,radius)-ADAF_gamma_r(radius&
            &,blackHoleSpin,adafViscosityAlpha)*dsqrt(Black_Hole_Metric_D_Factor(blackHoleSpin,radius)&
            &/Black_Hole_Metric_A_Factor(blackHoleSpin,radius))*2.0d0*ADAF_Angular_Momentum(radius,blackHoleSpin&
            &,adafViscosityAlpha)*Black_Hole_Frame_Dragging_Frequency(blackHoleSpin,radius)*ADAF_gamma_phi(radius,blackHoleSpin&
            &,adafViscosityAlpha)*blackHoleSpin**2
       nuz2=nuz2/radius**4
       adafHeightPrevious=dsqrt(ADAF_Temperature(radius,blackHoleSpin,adafViscosityAlpha)/ADAF_Enthalpy(radius,blackHoleSpin &
            &,adafViscosityAlpha)/radius**2/nuz2)
       radiusPrevious            =radius
       blackHoleSpinPrevious     =blackHoleSpin
       adafViscosityAlphaPrevious=adafViscosityAlpha
    end if
    ADAF_Height=adafHeightPrevious
    return
  end function ADAF_Height

end module Accretion_Disks_ADAF

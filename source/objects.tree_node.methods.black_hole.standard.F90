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


!% Contains a module of black hole tree node methods.

module Tree_Node_Methods_Black_Hole
  !% Implement black hole tree node methods.
  use Tree_Nodes
  use Components
  private
  public :: Tree_Node_Methods_Black_Hole_Initialize, Galacticus_Output_Tree_Black_Hole_Standard,&
       & Galacticus_Output_Tree_Black_Hole_Standard_Property_Count, Galacticus_Output_Tree_Black_Hole_Standard_Names,&
       & Tree_Node_Black_Hole_Reset_Standard, Black_Hole_Satellite_Merging, Black_Hole_Hot_Halo_Heating,&
       & Tree_Node_Methods_Black_Hole_Standard_Dump
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=2, dataCount=0, historyCount=0
  integer, parameter :: massIndex=1
  integer, parameter :: spinIndex=2

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Black_Hole_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Black_Hole_Spin</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

  ! Seed mass for black holes.
  double precision :: blackHoleSeedMass

  ! Seed spin for black holes.
  double precision, parameter :: blackHoleSeedSpin=1.0d-3

  ! Accretion model parameters.
  ! Enhancement factors for the accretion rate.
  double precision :: bondiHoyleAccretionEnhancementSpheroid,bondiHoyleAccretionEnhancementHotHalo
  ! Temperature of accreting gas.
  double precision :: bondiHoyleAccretionTemperatureSpheroid

  ! Feedback parameters.
  double precision :: blackHoleWindEfficiency

  ! Quantities stored to avoid repeated computation.
  logical          :: gotAccretionRate=.false.
  double precision :: accretionRate,accretionRateHotHalo
  !$omp threadprivate(gotAccretionRate,accretionRate,accretionRateHotHalo)

contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Black_Hole_Initialize</unitName>
  !#  <optionName default="standard">treeNodeMethodBlackHole</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Black_Hole_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node black hole methods module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'standard') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Standard black hole method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Black_Hole_Mass              => Tree_Node_Black_Hole_Mass_Standard
       Tree_Node_Black_Hole_Mass_Set          => Tree_Node_Black_Hole_Mass_Set_Standard
       Tree_Node_Black_Hole_Mass_Rate_Adjust  => Tree_Node_Black_Hole_Mass_Rate_Adjust_Standard
       Tree_Node_Black_Hole_Mass_Rate_Compute => Tree_Node_Black_Hole_Mass_Rate_Compute_Standard
       Tree_Node_Black_Hole_Spin              => Tree_Node_Black_Hole_Spin_Standard
       Tree_Node_Black_Hole_Spin_Set          => Tree_Node_Black_Hole_Spin_Set_Standard
       Tree_Node_Black_Hole_Spin_Rate_Adjust  => Tree_Node_Black_Hole_Spin_Rate_Adjust_Standard
       Tree_Node_Black_Hole_Spin_Rate_Compute => Tree_Node_Black_Hole_Spin_Rate_Compute_Standard

       ! Get the black hole seed mass.
       !@ <inputParameter>
       !@   <name>blackHoleSeedMass</name>
       !@   <defaultValue>100</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass of the seed black hole placed at the center of each newly formed galaxy.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleSeedMass",blackHoleSeedMass,defaultValue=100.0d0)

       ! Get accretion rate enhancement factors.
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionEnhancementSpheroid</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which the Bondi-Hoyle accretion rate of spheroid gas onto black holes in enhanced.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionEnhancementSpheroid",bondiHoyleAccretionEnhancementSpheroid,defaultValue&
            &=1.0d0)
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionEnhancementHotHalo</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which the Bondi-Hoyle accretion rate of hot halo gas onto black holes in enhanced.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionEnhancementHotHalo",bondiHoyleAccretionEnhancementHotHalo,defaultValue=1.0d0)

       ! Get temperature of accreting gas.
       !@ <inputParameter>
       !@   <name>bondiHoyleAccretionTemperatureSpheroid</name>
       !@   <defaultValue>$10^2$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The assumed temperature (in Kelvin) of gas in the spheroid when computing Bondi-Hoyle accretion rates onto black holes.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("bondiHoyleAccretionTemperatureSpheroid",bondiHoyleAccretionTemperatureSpheroid,defaultValue&
            &=1.0d2)

       ! Get temperature of accreting gas.
       !@ <inputParameter>
       !@   <name>blackHoleWindEfficiency</name>
       !@   <defaultValue>$10^{-3}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The efficiency of the black hole-driven wind: $L_{\rm wind} = \epsilon_{\rm wind} \dot{M}_\bullet \clight^2$.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("blackHoleWindEfficiency",blackHoleWindEfficiency,defaultValue=1.0d-3)

    end if
    return
  end subroutine Tree_Node_Methods_Black_Hole_Initialize
  

  !# <calculationResetTask>
  !# <unitName>Tree_Node_Black_Hole_Reset_Standard</unitName>
  !# </calculationResetTask>
  subroutine Tree_Node_Black_Hole_Reset_Standard(thisNode)
    !% Remove memory of stored computed values as we're about to begin computing derivatives anew.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    gotAccretionRate=.false.
    return
  end subroutine Tree_Node_Black_Hole_Reset_Standard

  double precision function Tree_Node_Black_Hole_Mass_Standard(thisNode)
    !% Return the node black hole mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Black_Hole_Index(thisNode)
       Tree_Node_Black_Hole_Mass_Standard=thisNode%components(thisIndex)%properties(massIndex,propertyValue)
    else
       Tree_Node_Black_Hole_Mass_Standard=blackHoleSeedMass
    end if
    return
  end function Tree_Node_Black_Hole_Mass_Standard

  subroutine Tree_Node_Black_Hole_Mass_Set_Standard(thisNode,mass)
    !% Set the node black hole mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Black_Hole_Index(thisNode)
    thisNode%components(thisIndex)%properties(massIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Black_Hole_Mass_Set_Standard

  subroutine Tree_Node_Black_Hole_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node black hole mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Black_Hole_Index(thisNode)
    thisNode%components(thisIndex)%properties(massIndex,propertyDerivative)=thisNode%components(thisIndex)%properties(massIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Black_Hole_Mass_Rate_Adjust_Standard

  subroutine Tree_Node_Black_Hole_Mass_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the black hole node mass rate of change.
    use Accretion_Disks
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    procedure(),      pointer                :: interruptProcedurePassed
    double precision, parameter              :: windVelocity=1.0d4   ! Velocity of disk wind.
    double precision, parameter              :: ismTemperature=1.0d4 ! Temperature of the ISM.
    double precision, parameter              :: criticalDensityNormalization=2.0d0*massHydrogenAtom*speedLight**2*megaParsec/3.0&
         &/Pi/boltzmannsConstant/gigaYear/ismTemperature/kilo/windVelocity
    double precision                         :: restMassAccretionRate,massAccretionRate,radiativeEfficiency,energyInputRate &
         &,spheroidDensity,spheroidGasMass,spheroidRadius,criticalDensity,windFraction,spheroidDensityOverCriticalDensity

    ! Find the rate of rest mass accretion onto the black hole.
    restMassAccretionRate=Mass_Accretion_Rate(thisNode)

    ! Finish if there is no accretion.
    if (restMassAccretionRate <= 0.0d0) return

    ! Find the radiative efficiency of the accretion.
    radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisNode,restMassAccretionRate)

    ! Find the rate of increase in mass of the black hole.
    massAccretionRate=restMassAccretionRate*(1.0d0-radiativeEfficiency)

    ! If no black hole component currently exists and we have some accretion then interrupt and create a black hole.
    if (.not.thisNode%componentExists(componentIndex)) then    
       if (massAccretionRate /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Black_Hole_Create
       end if
       return
    end if
    call Tree_Node_Black_Hole_Mass_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure, massAccretionRate   )

    ! Remove the accreted mass from the spheroid component.
    call Tree_Node_Spheroid_Gas_Sink                   (thisNode,interrupt,interruptProcedure,-accretionRate       )

    ! Remove the accreted mass from the hot halo component.
    call Tree_Node_Hot_Halo_Hot_Gas_Sink               (thisNode,interrupt,interruptProcedure,-accretionRateHotHalo)

    ! Add energy to the spheroid component.
    if (blackHoleWindEfficiency > 0.0d0) then
       spheroidGasMass=Tree_Node_Spheroid_Gas_Mass(thisNode)
       if (spheroidGasMass > 0.0d0) then
          spheroidRadius=Tree_Node_Spheroid_Radius(thisNode)
          if (spheroidRadius > 0.0d0) then
             spheroidDensity=3.0d0*spheroidGasMass/4.0d0/Pi/spheroidRadius**3
             criticalDensity=criticalDensityNormalization*blackHoleWindEfficiency*restMassAccretionRate/spheroidRadius**2
             
             ! Construct an interpolating factor such that the energy input from the wind drops to zero below half of the critical density.
             spheroidDensityOverCriticalDensity=spheroidDensity/criticalDensity-0.5d0
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
             call Tree_Node_Spheroid_Gas_Energy_Input(thisNode,interrupt,interruptProcedure,energyInputRate)
          end if
       end if
    end if
    return
  end subroutine Tree_Node_Black_Hole_Mass_Rate_Compute_Standard

  double precision function Tree_Node_Black_Hole_Spin_Standard(thisNode)
    !% Return the node black hole angular momentum.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    ! Maximum allowed spin (useful to avoid infinities at maximal spin).
    double precision, parameter              :: blackHoleSpinMaximum=0.9999d0
    integer                                  :: thisIndex

    if (thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Black_Hole_Index(thisNode)
       Tree_Node_Black_Hole_Spin_Standard=thisNode%components(thisIndex)%properties(spinIndex,propertyValue)
       ! Ensure that the spin is within range.
       if (Tree_Node_Black_Hole_Spin_Standard < 0.0d0) then
          Tree_Node_Black_Hole_Spin_Standard=0.0d0
          call Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,0.0d0)
       end if
       if (Tree_Node_Black_Hole_Spin_Standard > blackHoleSpinMaximum) then
          Tree_Node_Black_Hole_Spin_Standard=blackHoleSpinMaximum
          call Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,blackHoleSpinMaximum)
       end if
    else
       Tree_Node_Black_Hole_Spin_Standard=blackHoleSeedSpin
    end if
    return
  end function Tree_Node_Black_Hole_Spin_Standard

  subroutine Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,spin)
    !% Set the node black hole angular momentum.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: spin
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Black_Hole_Index(thisNode)
    thisNode%components(thisIndex)%properties(spinIndex,propertyValue)=spin
    return
  end subroutine Tree_Node_Black_Hole_Spin_Set_Standard

  subroutine Tree_Node_Black_Hole_Spin_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node black hole mass rate of change.
    use Cosmological_Parameters
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Black_Hole_Index(thisNode)
    thisNode%components(thisIndex)%properties(spinIndex,propertyDerivative)&
         &=thisNode%components(thisIndex)%properties(spinIndex,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Black_Hole_Spin_Rate_Adjust_Standard

  subroutine Tree_Node_Black_Hole_Spin_Rate_Compute_Standard(thisNode,interrupt,interruptProcedure)
    !% Compute the black hole node mass rate of change.
    use Accretion_Disks
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    procedure(),    pointer                :: interruptProcedurePassed
    double precision                       :: massAccretionRate,spinUpRate

    ! Find the rate of gas mass accretion onto the halo.
    massAccretionRate=Mass_Accretion_Rate(thisNode)
    ! If no black hole component currently exists and we have some accretion then interrupt and create a black hole.
    if (.not.thisNode%componentExists(componentIndex)) then    
       if (massAccretionRate /= 0.0d0) then
          interrupt=.true.
          interruptProcedure => Black_Hole_Create
       end if
       return
    end if

    if (massAccretionRate > 0.0d0) then
       spinUpRate=Black_Hole_Spin_Up_Rate(thisNode,massAccretionRate)
       call Tree_Node_Black_Hole_Spin_Rate_Adjust_Standard(thisNode,interrupt,interruptProcedure,spinUpRate)
    end if

    return
  end subroutine Tree_Node_Black_Hole_Spin_Rate_Compute_Standard


  !# <hotHaloHeatingTask>
  !#  <unitName>Black_Hole_Hot_Halo_Heating</unitName>
  !# </hotHaloHeatingTask>
  subroutine Black_Hole_Hot_Halo_Heating(thisNode,heatingRate)
    !% Compute the heating rate of the hot halo due to jets launched from the black hole.
    use Tree_Nodes
    use Accretion_Disks
    use Cooling_Radii
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(out)   :: heatingRate
    double precision                         :: massAccretionRate

    ! Check if halo is cooling sufficiently slowly that jet can couple to it.
    if (.not.methodSelected .or. Cooling_Radius(thisNode) >= Dark_Matter_Halo_Virial_Radius(thisNode)) then

       ! Halo is cooling rapidly, assume no heating can occur as jet cannot couple to halo gas.
       heatingRate=0.0d0
       
    else
       
       ! Find the rate of gas mass accretion onto the halo.
       massAccretionRate=Mass_Accretion_Rate(thisNode)
       
       ! Get jet power.
       heatingRate=Accretion_Disk_Jet_Power(thisNode,massAccretionRate)

    end if

    return
  end subroutine Black_Hole_Hot_Halo_Heating

  !# <satelliteMergerTask>
  !#  <unitName>Black_Hole_Satellite_Merging</unitName>
  !# </satelliteMergerTask>
  subroutine Black_Hole_Satellite_Merging(thisNode)
    !% Merge (instantaneously) any black hole associated with {\tt thisNode} before it merges with its host halo.
    use Black_Hole_Binary_Mergers
    implicit none
    type(treeNode),   pointer, intent(inout)     :: thisNode
    type(treeNode),   pointer                    :: hostNode
    double precision                             :: blackHoleMassNew,blackHoleSpinNew

    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       
       ! Find the host node.
       hostNode => thisNode%parentNode

       call Black_Hole_Binary_Merger(Tree_Node_Black_Hole_Mass_Standard(thisNode), &
            &                        Tree_Node_Black_Hole_Mass_Standard(hostNode), &
            &                        Tree_Node_Black_Hole_Spin_Standard(thisNode), &
            &                        Tree_Node_Black_Hole_Spin_Standard(hostNode), &
            &                        blackHoleMassNew                            , &
            &                        blackHoleSpinNew                             )

       ! Move the black hole to the host.
       call Tree_Node_Black_Hole_Mass_Set_Standard(hostNode,blackHoleMassNew)
       call Tree_Node_Black_Hole_Spin_Set_Standard(hostNode,blackHoleSpinNew)
       call Tree_Node_Black_Hole_Mass_Set_Standard(thisNode,0.0d0           )
       call Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,0.0d0           )
    end if
    return
  end subroutine Black_Hole_Satellite_Merging
  

  double precision function Mass_Accretion_Rate(thisNode)
    !% Returns the rate of mass accretion onto the black hole in {\tt thisNode}.
    use Cosmological_Parameters
    use Bondi_Hoyle_Lyttleton_Accretion
    use Galactic_Structure_Densities
    use Galactic_Structure_Options
    use Ideal_Gases_Thermodynamics
    use Black_Hole_Fundamentals
    use Accretion_Disks
    use Hot_Halo_Temperature_Profile
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision                         :: blackHoleMass,gasDensity,relativeVelocity,accretionRadius,jeansLength&
         &,radiativeEfficiency, position(3),hotHaloTemperature

    if (.not.gotAccretionRate) then
       ! Get black hole mass.
       blackHoleMass=Tree_Node_Black_Hole_Mass_Standard(thisNode)

       ! Check black hole mass is positive.
       if (blackHoleMass > 0.0d0) then
          ! Assume black hole is stationary with respect to surrounding gas.
          relativeVelocity=0.0d0

          ! Contribution from spheroid:

          ! Get the accretion radius.
          accretionRadius=Bondi_Hoyle_Lyttleton_Accretion_Radius(blackHoleMass,bondiHoyleAccretionTemperatureSpheroid)
          
          ! Set the position.
          position=[accretionRadius,0.0d0,0.0d0]
          ! Get density of gas at the galactic center.
          gasDensity=Galactic_Structure_Density(thisNode,position,coordinateSystem=coordinateSystemSpherical,massType&
               &=massTypeGaseous,componentType=componentTypeSpheroid)
          
          ! Check if we have a non-zero gas density.
          if (gasDensity > 0.0d0) then
             ! Get the Jeans length scale.
             jeansLength=Ideal_Gas_Jeans_Length(bondiHoyleAccretionTemperatureSpheroid,gasDensity)

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
             accretionRate=bondiHoyleAccretionEnhancementSpheroid*Bondi_Hoyle_Lyttleton_Accretion_Rate(blackHoleMass,gasDensity&
                  &,relativeVelocity,bondiHoyleAccretionTemperatureSpheroid)
             
             ! Get the radiative efficiency of the accretion.
             radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisNode,accretionRate)
             
             ! Limit the accretion rate to the Eddington limit.
             if (radiativeEfficiency > 0.0d0) accretionRate=min(accretionRate,Black_Hole_Eddington_Accretion_Rate(thisNode)&
                  &/radiativeEfficiency)

          else
             ! Gas density is negative - set zero accretion rate.
             accretionRate=0.0d0
          end if

          ! Contribution from hot halo:

          ! Get halo gas temperature.
          hotHaloTemperature=Hot_Halo_Temperature(thisNode,radius=0.0d0)

          ! Get the accretion radius.
          accretionRadius=Bondi_Hoyle_Lyttleton_Accretion_Radius(blackHoleMass,hotHaloTemperature)
          
          ! Set the position.
          position=[accretionRadius,0.0d0,0.0d0]

          ! Get density of gas at the galactic center.
          gasDensity=Galactic_Structure_Density(thisNode,position,coordinateSystem=coordinateSystemSpherical,massType&
               &=massTypeGaseous,componentType=componentTypeHotHalo)

          ! Check if we have a non-zero gas density.
          if (gasDensity > 0.0d0) then
             
             ! Compute the accretion rate.
             accretionRateHotHalo=bondiHoyleAccretionEnhancementHotHalo*Bondi_Hoyle_Lyttleton_Accretion_Rate(blackHoleMass,gasDensity&
                  &,relativeVelocity,hotHaloTemperature)

             ! Get the radiative efficiency of the accretion.
             radiativeEfficiency=Accretion_Disk_Radiative_Efficiency(thisNode,accretionRateHotHalo)
             
             ! Limit the accretion rate to the Eddington limit.
             if (radiativeEfficiency > 0.0d0) accretionRateHotHalo=min(accretionRateHotHalo&
                  &,Black_Hole_Eddington_Accretion_Rate(thisNode)/radiativeEfficiency)

          else
             ! No gas density, so zero accretion rate.
             accretionRateHotHalo=0.0d0
          end if

       else
          accretionRate       =0.0d0
          accretionRateHotHalo=0.0d0
       end if
       gotAccretionRate=.true.
    end if
    Mass_Accretion_Rate=accretionRate+accretionRateHotHalo
    return
  end function Mass_Accretion_Rate

  integer function Tree_Node_Black_Hole_Index(thisNode)
    !% Ensure the black hole component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Black_Hole_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Black_Hole_Index

  subroutine Black_Hole_Create(thisNode)
    !% Creates a black hole component for {\tt thisNode}.
    use ISO_Varying_String
    use Galacticus_Display
    use String_Handling
    implicit none
    type(treeNode),      pointer, intent(inout) :: thisNode
    type(varying_string)                        :: message

    ! Display a message.
    message='Creating black hole component for node '
    message=message//thisNode%index()
    call Galacticus_Display_Message(message,verbosityInfo)
    ! Create the component.
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    ! Set to the seed mass.
    call Tree_Node_Black_Hole_Mass_Set_Standard(thisNode,blackHoleSeedMass)
    call Tree_Node_Black_Hole_Spin_Set_Standard(thisNode,blackHoleSeedSpin)
    return
  end subroutine Black_Hole_Create

  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Standard_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Standard</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Black_Hole_Standard_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,doubleProperty,doublePropertyNames,doublePropertyComments,time)
    !% Set names of black hole properties to be written to the \glc\ output file.
    use ISO_Varying_String
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='blackHoleMass'
       doublePropertyComments(doubleProperty)='Mass of the black hole.'
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='blackHoleSpin'
       doublePropertyComments(doubleProperty)='Spin of the black hole.'
     end if
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Standard_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Standard_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Standard</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Black_Hole_Standard_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of black hole properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+propertyCount
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Standard_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Black_Hole_Standard</unitName>
  !#  <sortName>Galacticus_Output_Tree_Black_Hole_Standard</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Black_Hole_Standard(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store black hole properties in the \glc\ output file buffers.
    use Tree_Nodes
    implicit none
    double precision, intent(in)                 :: time
    type(treeNode),   intent(inout), pointer     :: thisNode
    integer,          intent(inout)              :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer,          intent(inout)              :: integerBuffer(:,:)
    double precision, intent(inout)              :: doubleBuffer(:,:)

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Black_Hole_Mass_Standard(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Black_Hole_Spin_Standard(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Black_Hole_Standard

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Black_Hole_Standard_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Black_Hole_Standard_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'black hole component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'black hole mass:',Tree_Node_Black_Hole_Mass(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'black hole spin:',Tree_Node_Black_Hole_Spin(thisNode)
       else
          write (0,'(1x,a)'           ) 'black hole component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Black_Hole_Standard_Dump

end module Tree_Node_Methods_Black_Hole

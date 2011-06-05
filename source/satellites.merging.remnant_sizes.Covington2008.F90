!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements the \cite{covington_predicting_2008} algorithm for merger remnant sizes.

module Satellite_Merging_Remnant_Sizes_Covington2008
  !% Implements the \cite{covington_predicting_2008} algorithm for merger remnant sizes.
  private
  public :: Satellite_Merging_Remnant_Sizes_Covington2008_Initialize

  ! Parameter controlling the orbital energy used in the calculation.
  double precision :: mergerRemnantSizeOrbitalEnergy

  ! Parameter controlling the radiative efficiency used in the calculation.
  double precision :: mergerRemnantRadiativeEfficiency

contains

  !# <satelliteMergingRemnantSizeMethod>
  !#  <unitName>Satellite_Merging_Remnant_Sizes_Covington2008_Initialize</unitName>
  !# </satelliteMergingRemnantSizeMethod>
  subroutine Satellite_Merging_Remnant_Sizes_Covington2008_Initialize(satelliteMergingRemnantSizeMethod,Satellite_Merging_Remnant_Size_Do)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingRemnantSizeMethod
    procedure(),          pointer, intent(inout) :: Satellite_Merging_Remnant_Size_Do
    
    if (satelliteMergingRemnantSizeMethod == 'Covington2008') then
       Satellite_Merging_Remnant_Size_Do => Satellite_Merging_Remnant_Size_Covington2008
       !@ <inputParameter>
       !@   <name>mergerRemnantSizeOrbitalEnergy</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The orbital energy used in the ``Covington2008'' merger remnant sizes calculation in units of the characteristic orbital energy.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("mergerRemnantSizeOrbitalEnergy",mergerRemnantSizeOrbitalEnergy,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>mergerRemnantRadiativeEfficiency</name>
       !@   <defaultValue>2.75 \citep{covington_predicting_2008}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The coefficient, $C_{\rm rad}$ energy used in the \cite{covington_predicting_2008} merger remnant size algorithm.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("mergerRemnantRadiativeEfficiency",mergerRemnantRadiativeEfficiency,defaultValue=2.75d0)
    end if
    return
  end subroutine Satellite_Merging_Remnant_Sizes_Covington2008_Initialize

  subroutine Satellite_Merging_Remnant_Size_Covington2008(thisNode)
    !% Compute the size of the merger remnant for {\tt thisNode} using the \cite{covington_predicting_2008} algorithm.
    use Tree_Nodes
    use Numerical_Constants_Physical
    use Numerical_Comparison
    use Satellite_Merging_Remnant_Sizes_Properties
    use Galacticus_Error
    use String_Handling
    use ISO_Varying_String
    use Galacticus_Display
    use Satellite_Merging_Remnant_Sizes_Utilities
    implicit none
    type(treeNode),          intent(inout), pointer  :: thisNode
    type(treeNode),                         pointer  :: hostNode
    double precision,        parameter               :: bindingEnergyFormFactor=0.5d+0
    double precision,        parameter               :: absoluteMassTolerance  =1.0d-6
    double precision,        parameter               :: relativeMassTolerance  =1.0d-9
    double precision                                 :: satelliteMass,hostMass,satelliteRadius,hostRadius,satelliteSpheroidMass &
         &,hostSpheroidMass,progenitorsEnergy,hostSpheroidMassPreMerger,darkMatterFactor,remnantSpheroidGasMass&
         &,remnantSpheroidMass,gasFractionInitial,radiatedEnergy
    character(len= 2)                                :: joinString
    character(len=40)                                :: dataString
    type(varying_string)                             :: message
    logical                                          :: errorCondition

    ! Get the host node.
    call thisNode%mergesWith(hostNode)

    ! Get properties of the merging systems.
    call Satellite_Merging_Remnant_Progenitor_Properties(thisNode,hostNode,satelliteMass,hostMass,satelliteSpheroidMass &
         &,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,darkMatterFactor,remnantSpheroidMass&
         &,remnantSpheroidGasMass)

    if (satelliteSpheroidMass <= 0.0d0 .and. Values_Agree(hostSpheroidMass,hostSpheroidMassPreMerger,relTol=relativeMassTolerance)) then
       remnantRadius                 =remnantNoChangeValue
       remnantCircularVelocity       =remnantNoChangeValue
       remnantSpecificAngularMomentum=remnantNoChangeValue
    else       
       ! Check that the properties of the galaxies are physically reasonable.
       errorCondition=.false.
       if (satelliteRadius <= 0.0d0 .or. satelliteMass < -absoluteMassTolerance .or. satelliteSpheroidMass < -absoluteMassTolerance) then
          write (dataString,'(3(e12.6,":",e12.6,":",e12.6))') satelliteRadius,satelliteMass,satelliteSpheroidMass
          message='Satellite galaxy ['
          message=message//thisNode%index()//'] has '
          joinString=""
          if (satelliteRadius       <= 0.0d0         ) then
             message=message//trim(joinString)//'non-positive radius'
             joinString=", "
          end if
          if (satelliteMass         <  -absoluteMassTolerance) then
             message=message//trim(joinString)//'negative mass'
             joinString=", "
          end if
          if (satelliteSpheroidMass <  -absoluteMassTolerance) then
             message=message//trim(joinString)//'negative spheroid mass'
             joinString=", "
          end if
          message=message//' (radius:mass:spheroidMass='//trim(dataString)//')'
          call Galacticus_Display_Message(message)
          errorCondition=.true.
       end if
       if ((hostRadius <= 0.0d0 .and. hostMass > 0.0d0) .or. hostMass < -absoluteMassTolerance .or. hostSpheroidMass < -absoluteMassTolerance) then
          write (dataString,'(3(e12.6,":",e12.6,":",e12.6))') hostRadius,hostMass,hostSpheroidMass
          message='Host galaxy ['
          message=message//hostNode%index()//'] has '
          joinString=""
          if (hostRadius       <= 0.0d0         ) then
             message=message//trim(joinString)//'non-positive radius'
             joinString=", "
          end if
          if (hostMass         <  -absoluteMassTolerance) then
             message=message//trim(joinString)//'negative mass'
             joinString=", "
          end if
          if (hostSpheroidMass <  -absoluteMassTolerance) then
             message=message//trim(joinString)//'negative spheroid mass'
             joinString=", "
          end if
          message=message//' (radius:mass:spheroidMass='//trim(dataString)//')'
          call Galacticus_Display_Message(message)
          errorCondition=.true.
       end if
       if (errorCondition) call Galacticus_Error_Report('Satellite_Merging_Remnant_Size_Covington2008','error condition detected')
       ! Check if host has finite mass.
       if (hostMass > 0.0d0) then
          ! Apply the Covington et al. (2008) algorithm to compute the size of the new remnant.
          ! First calculate the energy of the progenitors.
          progenitorsEnergy= satelliteSpheroidMass*satelliteMass/satelliteRadius                                                &
               &            +hostSpheroidMass     *hostMass     /hostRadius                                                     &
               &            +mergerRemnantSizeOrbitalEnergy*satelliteSpheroidMass*hostSpheroidMass/(satelliteRadius+hostRadius) &
               &                                                                                  /bindingEnergyFormFactor
          ! Compute the gas fraction in the remnant.
          gasFractionInitial=remnantSpheroidGasMass/remnantSpheroidMass
          ! Compute the energy lost through radiation.
          radiatedEnergy=mergerRemnantRadiativeEfficiency*gasFractionInitial*progenitorsEnergy

          ! Compute the remnant radius.
          remnantRadius=(satelliteSpheroidMass+hostSpheroidMass)**2/(progenitorsEnergy+radiatedEnergy)
       else
          remnantRadius=satelliteRadius
       end if

       ! Also compute the specific angular momentum at the half-mass radius.
       remnantCircularVelocity=dsqrt(gravitationalConstantGalacticus*(satelliteSpheroidMass+hostSpheroidMass)/remnantRadius)
       remnantSpecificAngularMomentum=remnantRadius*remnantCircularVelocity*darkMatterFactor
    end if
    return
  end subroutine Satellite_Merging_Remnant_Size_Covington2008

end module Satellite_Merging_Remnant_Sizes_Covington2008

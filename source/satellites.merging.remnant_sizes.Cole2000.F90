!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the \cite{cole_hierarchical_2000} algorithm for merger remnant sizes.

module Satellite_Merging_Remnant_Sizes_Cole2000
  !% Implements the \cite{cole_hierarchical_2000} algorithm for merger remnant sizes.
  implicit none
  private
  public :: Satellite_Merging_Remnant_Sizes_Cole2000_Initialize

  ! Parameter controlling the orbital energy used in the calculation.
  double precision :: mergerRemnantSizeOrbitalEnergy

contains

  !# <satelliteMergingRemnantSizeMethod>
  !#  <unitName>Satellite_Merging_Remnant_Sizes_Cole2000_Initialize</unitName>
  !# </satelliteMergingRemnantSizeMethod>
  subroutine Satellite_Merging_Remnant_Sizes_Cole2000_Initialize(satelliteMergingRemnantSizeMethod,Satellite_Merging_Remnant_Size_Do)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingRemnantSizeMethod
    procedure(),          pointer, intent(inout) :: Satellite_Merging_Remnant_Size_Do
    
    if (satelliteMergingRemnantSizeMethod == 'Cole2000') then
       Satellite_Merging_Remnant_Size_Do => Satellite_Merging_Remnant_Size_Cole2000
       !@ <inputParameter>
       !@   <name>mergerRemnantSizeOrbitalEnergy</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The orbital energy used in the ``Cole2000'' merger remnant sizes calculation in units of the characteristic orbital energy.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("mergerRemnantSizeOrbitalEnergy",mergerRemnantSizeOrbitalEnergy,defaultValue=1.0d0)
    end if
    return
  end subroutine Satellite_Merging_Remnant_Sizes_Cole2000_Initialize

  subroutine Satellite_Merging_Remnant_Size_Cole2000(thisNode)
    !% Compute the size of the merger remnant for {\tt thisNode} using the \cite{cole_hierarchical_2000} algorithm.
    use Galacticus_Nodes
    use Numerical_Constants_Physical
    use Numerical_Comparison
    use Satellite_Merging_Remnant_Sizes_Properties
    use Galacticus_Error
    use String_Handling
    use ISO_Varying_String
    use Galacticus_Display
    use Satellite_Merging_Remnant_Sizes_Progenitors
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    type(treeNode),          intent(inout), pointer  :: thisNode
    type(treeNode),                         pointer  :: hostNode
    double precision,        parameter               :: bindingEnergyFormFactor=0.5d+0
    double precision,        parameter               :: absoluteMassTolerance  =1.0d-6
    double precision,        parameter               :: relativeMassTolerance  =1.0d-9
    double precision                                 :: satelliteMass,hostMass,satelliteRadius,hostRadius,satelliteSpheroidMass &
         &,hostSpheroidMass,progenitorsEnergy,hostSpheroidMassPreMerger,angularMomentumFactor,remnantSpheroidGasMass &
         &,remnantSpheroidMass,hostDarkMatterBoost,satelliteDarkMatterBoost,hostSpheroidMassTotal,satelliteSpheroidMassTotal
    character(len= 2)                                :: joinString
    character(len=40)                                :: dataString
    type(varying_string)                             :: message
    logical                                          :: errorCondition

    ! Find the node to merge with.
    hostNode => thisNode%mergesWith()

    ! Get properties of the merging systems.
    call Satellite_Merging_Remnant_Progenitor_Properties(thisNode,hostNode,satelliteMass,hostMass,satelliteSpheroidMass &
         &,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass&
         &,remnantSpheroidGasMass)
    if (satelliteSpheroidMass <= 0.0d0 .and. Values_Agree(hostSpheroidMass,hostSpheroidMassPreMerger,relTol=relativeMassTolerance)) then
       remnantRadius                 =remnantNoChangeValue
       remnantCircularVelocity       =remnantNoChangeValue
       remnantSpecificAngularMomentum=remnantNoChangeValue
    else       
       ! Check that the properties of the galaxies are physically reasonable.
       errorCondition=.false.
       if     (                                                                         &
            &      (satelliteRadius       <= 0.0d0 .and. satelliteSpheroidMass > 0.0d0) &
            &  .or. satelliteMass         < -absoluteMassTolerance                      &
            &  .or. satelliteSpheroidMass < -absoluteMassTolerance                      &
            & ) then
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
       if ((hostRadius <= 0.0d0 .and. hostSpheroidMass > 0.0d0) .or. hostMass < -absoluteMassTolerance .or. hostSpheroidMass < -absoluteMassTolerance) then
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
       if (errorCondition) call Galacticus_Error_Report('Satellite_Merging_Remnant_Size_Cole2000','error condition detected')
       ! Check if host has finite mass.
       if (satelliteSpheroidMass+hostSpheroidMass > 0.0d0) then
          ! Find the contribution of dark matter to the masses enclosed within the galaxy radii.
          if (hostSpheroidMass      > 0.0d0) then
             hostDarkMatterBoost     =1.0d0+Galactic_Structure_Enclosed_Mass(hostNode,hostRadius     ,massType=massTypeDark)/hostSpheroidMass
          else
             hostDarkMatterBoost     =1.0d0
          end if
          if (satelliteSpheroidMass > 0.0d0) then
             satelliteDarkMatterBoost=1.0d0+Galactic_Structure_Enclosed_Mass(thisNode,satelliteRadius,massType=massTypeDark)/satelliteSpheroidMass
          else
             satelliteDarkMatterBoost=1.0d0
          end if
          ! Scale masses to account for dark matter.
          hostSpheroidMassTotal     =hostSpheroidMass     *hostDarkMatterBoost
          satelliteSpheroidMassTotal=satelliteSpheroidMass*satelliteDarkMatterBoost
          ! Apply the Cole et al. (2000) algorithm to compute the size of the new remnant.
          progenitorsEnergy=0.0d0
          if (hostRadius                 > 0.0d0)                                                                                     &
               & progenitorsEnergy=progenitorsEnergy+                           hostSpheroidMassTotal**2/                 hostRadius
          if (           satelliteRadius > 0.0d0)                                                                                     & 
               & progenitorsEnergy=progenitorsEnergy+satelliteSpheroidMassTotal                      **2/ satelliteRadius  
          if (hostRadius+satelliteRadius > 0.0d0)                                                                                     &
               & progenitorsEnergy=progenitorsEnergy+satelliteSpheroidMassTotal*hostSpheroidMassTotal   /(satelliteRadius+hostRadius) &
               &                                    *mergerRemnantSizeOrbitalEnergy/bindingEnergyFormFactor
          remnantRadius=(satelliteSpheroidMassTotal+hostSpheroidMassTotal)**2/progenitorsEnergy
          
          ! Also compute the specific angular momentum at the half-mass radius.
          remnantCircularVelocity=dsqrt(gravitationalConstantGalacticus*(satelliteSpheroidMass+hostSpheroidMass)/remnantRadius)
          remnantSpecificAngularMomentum=remnantRadius*remnantCircularVelocity*angularMomentumFactor
       else
          ! Remnant has zero mass - don't do anything.
          remnantRadius                 =remnantNoChangeValue
          remnantCircularVelocity       =remnantNoChangeValue
          remnantSpecificAngularMomentum=remnantNoChangeValue
       end if       
    end if
    return
  end subroutine Satellite_Merging_Remnant_Size_Cole2000

end module Satellite_Merging_Remnant_Sizes_Cole2000

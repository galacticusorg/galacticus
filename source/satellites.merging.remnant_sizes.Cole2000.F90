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


!% Contains a module which implements the \cite{cole_hierarchical_2000} algorithm for merger remnant sizes.

module Satellite_Merging_Remnant_Sizes_Cole2000
  !% Implements the \cite{cole_hierarchical_2000} algorithm for merger remnant sizes.
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
       !@ </inputParameter>
       call Get_Input_Parameter("mergerRemnantSizeOrbitalEnergy",mergerRemnantSizeOrbitalEnergy,defaultValue=1.0d0)
    end if
    return
  end subroutine Satellite_Merging_Remnant_Sizes_Cole2000_Initialize

  subroutine Satellite_Merging_Remnant_Size_Cole2000(thisNode)
    !% Compute the size of the merger remnant for {\tt thisNode} using the \cite{cole_hierarchical_2000} algorithm.
    use Tree_Nodes
    use Tree_Node_Methods
    use Galactic_Structure_Radii
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Satellite_Merging_Mass_Movements_Descriptors
    use Numerical_Constants_Physical
    use Satellite_Merging_Remnant_Sizes_Properties
    use Galactic_Structure_Rotation_Curves
    use Numerical_Interpolation
    use FGSL
    use Galacticus_Error
    use String_Handling
    use ISO_Varying_String
    use Galacticus_Display
    implicit none
    type(treeNode),          intent(inout), pointer  :: thisNode
    type(treeNode),                         pointer  :: hostNode
    double precision,        parameter               :: bindingEnergyFormFactor=0.5d+0
    double precision,        parameter               :: massTolerance          =1.0d-6
    type(fgsl_interp),       save                    :: interpolationObject
    type(fgsl_interp_accel), save                    :: interpolationAccelerator
    !$omp threadprivate(interpolationObject,interpolationAccelerator)
    logical                                          :: interpolationReset
    double precision                                 :: satelliteMass,hostMass,satelliteRadius,hostRadius,satelliteSpheroidMass &
         &,hostSpheroidMass,progenitorsEnergy,hostSpheroidMassPreMerger,hostSpheroidDarkMatterFactor,hostDiskDarkMatterFactor&
         &,satelliteSpheroidDarkMatterFactor,satelliteDiskDarkMatterFactor,darkMatterFactor,componentMass
    character(len= 2)                                :: joinString
    character(len=40)                                :: dataString
    type(varying_string)                             :: message
    logical                                          :: errorCondition

    ! Get the host node.
    hostNode => thisNode%parentNode

    ! Solve for the radii of the host and satellite nodes, to ensure they are computed and up to date.
    call Galactic_Structure_Radii_Solve(hostNode)
    call Galactic_Structure_Radii_Solve(thisNode)

    ! Find the baryonic masses of the two galaxies.
    satelliteMass=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic)
    hostMass     =Galactic_Structure_Enclosed_Mass(hostNode,massType=massTypeGalactic)

    ! Find the half-mass radii of the two galaxies.
    satelliteRadius=Galactic_Structure_Radius_Enclosing_Mass(thisNode,fractionalMass=0.5d0,massType=massTypeGalactic)
    hostRadius     =Galactic_Structure_Radius_Enclosing_Mass(hostNode,fractionalMass=0.5d0,massType=massTypeGalactic)

    ! Compute dark matter factors. These are the specific angular momenta of components divided by sqrt(G M / r) where M is the
    ! component mass and r its half-mass radius. We use a weighted average of these factors to infer the specific angular momentum
    ! of the remnant from its mass and radius.
    componentMass=Tree_Node_Spheroid_Stellar_Mass(hostNode)+Tree_Node_Spheroid_Gas_Mass(hostNode)
    if (Tree_Node_Spheroid_Half_Mass_Radius(hostNode) > 0.0d0 .and. componentMass > 0.0d0) then
       hostSpheroidDarkMatterFactor=Tree_Node_Spheroid_Angular_Momentum(hostNode)/(componentMass**1.5d0) &
            &/dsqrt(gravitationalConstantGalacticus*Tree_Node_Spheroid_Half_Mass_Radius(hostNode))
    else
       hostSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass=Tree_Node_Disk_Stellar_Mass(hostNode) +Tree_Node_Disk_Gas_Mass(hostNode)
    if (Tree_Node_Disk_Half_Mass_Radius(hostNode) > 0.0d0 .and. componentMass > 0.0d0) then
       hostDiskDarkMatterFactor=Tree_Node_Disk_Angular_Momentum(hostNode)/(componentMass**1.5d0)&
            &/dsqrt(gravitationalConstantGalacticus *Tree_Node_Disk_Half_Mass_Radius(hostNode))
    else
       hostDiskDarkMatterFactor=0.0d0
    end if
    componentMass=Tree_Node_Spheroid_Stellar_Mass(thisNode) +Tree_Node_Spheroid_Gas_Mass(thisNode)
    if (Tree_Node_Spheroid_Half_Mass_Radius(thisNode) > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteSpheroidDarkMatterFactor=Tree_Node_Spheroid_Angular_Momentum(thisNode)/(componentMass**1.5d0)&
            &/dsqrt(gravitationalConstantGalacticus *Tree_Node_Spheroid_Half_Mass_Radius(thisNode))
    else
       satelliteSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass=Tree_Node_Disk_Stellar_Mass(thisNode) +Tree_Node_Disk_Gas_Mass(thisNode)
    if (Tree_Node_Disk_Half_Mass_Radius(thisNode) > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteDiskDarkMatterFactor=Tree_Node_Disk_Angular_Momentum(thisNode)/(componentMass**1.5d0)&
            &/dsqrt(gravitationalConstantGalacticus *Tree_Node_Disk_Half_Mass_Radius(thisNode))
    else
       satelliteDiskDarkMatterFactor=0.0d0
    end if

    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (thisHostGasMovesTo)
    case (movesToSpheroid)
       hostSpheroidMass=Tree_Node_Spheroid_Gas_Mass(hostNode)                             +Tree_Node_Disk_Gas_Mass(hostNode)
       darkMatterFactor=Tree_Node_Spheroid_Gas_Mass(hostNode)*hostSpheroidDarkMatterFactor+Tree_Node_Disk_Gas_Mass(hostNode)*hostDiskDarkMatterFactor
    case (movesToDisk)
       hostSpheroidMass=0.0d0
       darkMatterFactor=0.0d0
    case (doesNotMove)
       hostSpheroidMass=Tree_Node_Spheroid_Gas_Mass(hostNode)
       darkMatterFactor=Tree_Node_Spheroid_Gas_Mass(hostNode)*hostSpheroidDarkMatterFactor
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Size_Cole2000','unrecognized moveTo descriptor')
    end select
    select case (thisHostStarsMoveTo)
    case (movesToSpheroid)
       hostSpheroidMass=hostSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(hostNode)                             +Tree_Node_Disk_Stellar_Mass(hostNode)
       darkMatterFactor=darkMatterFactor+Tree_Node_Spheroid_Stellar_Mass(hostNode)*hostSpheroidDarkMatterFactor+Tree_Node_Disk_Stellar_Mass(hostNode)*hostDiskDarkMatterFactor
    case (movesToDisk)
       hostSpheroidMass=hostSpheroidMass
       darkMatterFactor=darkMatterFactor
    case (doesNotMove)
       hostSpheroidMass=hostSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(hostNode)
       darkMatterFactor=darkMatterFactor+Tree_Node_Spheroid_Stellar_Mass(hostNode)*hostSpheroidDarkMatterFactor
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Size_Cole2000','unrecognized moveTo descriptor')
    end select
    select case (thisMergerGasMovesTo)
    case (movesToSpheroid)
       satelliteSpheroidMass=                 Tree_Node_Spheroid_Gas_Mass(thisNode)                                  +Tree_Node_Disk_Gas_Mass(thisNode)
       darkMatterFactor     =darkMatterFactor+Tree_Node_Spheroid_Gas_Mass(thisNode)*satelliteSpheroidDarkMatterFactor+Tree_Node_Disk_Gas_Mass(thisNode)*satelliteDiskDarkMatterFactor
    case (movesToDisk)
       satelliteSpheroidMass=0.0d0
       darkMatterFactor     =darkMatterFactor
    case (doesNotMove)
       satelliteSpheroidMass=                 Tree_Node_Spheroid_Gas_Mass(thisNode)
       darkMatterFactor     =darkMatterFactor+Tree_Node_Spheroid_Gas_Mass(thisNode)*satelliteSpheroidDarkMatterFactor
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Size_Cole2000','unrecognized moveTo descriptor')
    end select
    select case (thisMergerStarsMoveTo)
    case (movesToSpheroid)
       satelliteSpheroidMass=satelliteSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(thisNode)                                  +Tree_Node_Disk_Stellar_Mass(thisNode)
       darkMatterFactor     =darkMatterFactor     +Tree_Node_Spheroid_Stellar_Mass(thisNode)*satelliteSpheroidDarkMatterFactor+Tree_Node_Disk_Stellar_Mass(thisNode)*satelliteDiskDarkMatterFactor
    case (movesToDisk)
       satelliteSpheroidMass=satelliteSpheroidMass
       darkMatterFactor     =darkMatterFactor
    case (doesNotMove)
       satelliteSpheroidMass=satelliteSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(thisNode)
       darkMatterFactor     =darkMatterFactor     +Tree_Node_Spheroid_Stellar_Mass(thisNode)*satelliteSpheroidDarkMatterFactor
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Size_Cole2000','unrecognized moveTo descriptor')
    end select
    if (satelliteSpheroidMass+hostSpheroidMass > 0.0d0) then
       darkMatterFactor=darkMatterFactor/(satelliteSpheroidMass+hostSpheroidMass)
    else
       darkMatterFactor=1.0d0
    end if

    hostSpheroidMassPreMerger=Tree_Node_Spheroid_Stellar_Mass(hostNode)+Tree_Node_Spheroid_Gas_Mass(hostNode)

    if (satelliteSpheroidMass <= 0.0d0 .and. hostSpheroidMass == hostSpheroidMassPreMerger) then
       remnantRadius                 =remnantNoChangeValue
       remnantCircularVelocity       =remnantNoChangeValue
       remnantSpecificAngularMomentum=remnantNoChangeValue
    else       
       ! Check that the properties of the galaxies are physically reasonable.
       errorCondition=.false.
       if (satelliteRadius <= 0.0d0 .or. satelliteMass < -massTolerance .or. satelliteSpheroidMass < -massTolerance) then
          write (dataString,'(3(e12.6,":",e12.6,":",e12.6))') satelliteRadius,satelliteMass,satelliteSpheroidMass
          message='Satellite galaxy ['
          message=message//thisNode%index()//'] has '
          joinString=""
          if (satelliteRadius       <= 0.0d0         ) then
             message=message//trim(joinString)//'non-positive radius'
             joinString=", "
          end if
          if (satelliteMass         <  -massTolerance) then
             message=message//trim(joinString)//'negative mass'
             joinString=", "
          end if
          if (satelliteSpheroidMass <  -massTolerance) then
             message=message//trim(joinString)//'negative spheroid mass'
             joinString=", "
          end if
          message=message//' (radius:mass:spheroidMass='//trim(dataString)//')'
          call Galacticus_Display_Message(message)
          errorCondition=.true.
       end if
       if ((hostRadius <= 0.0d0 .and. hostMass > 0.0d0) .or. hostMass < -massTolerance .or. hostSpheroidMass < -massTolerance) then
          write (dataString,'(3(e12.6,":",e12.6,":",e12.6))') hostRadius,hostMass,hostSpheroidMass
          message='Host galaxy ['
          message=message//hostNode%index()//'] has '
          joinString=""
          if (hostRadius       <= 0.0d0         ) then
             message=message//trim(joinString)//'non-positive radius'
             joinString=", "
          end if
          if (hostMass         <  -massTolerance) then
             message=message//trim(joinString)//'negative mass'
             joinString=", "
          end if
          if (hostSpheroidMass <  -massTolerance) then
             message=message//trim(joinString)//'negative spheroid mass'
             joinString=", "
          end if
          message=message//' (radius:mass:spheroidMass='//trim(dataString)//')'
          call Galacticus_Display_Message(message)
          errorCondition=.true.
       end if
       if (errorCondition) call Galacticus_Error_Report('Satellite_Merging_Remnant_Size_Cole2000','error condition detected')
       ! Check if host has finite mass.
       if (hostMass > 0.0d0) then
          ! Apply the Cole et al. (2000) algorithm to compute the size of the new remnant.
          progenitorsEnergy= satelliteSpheroidMass*satelliteMass/satelliteRadius &
               &            +hostSpheroidMass     *hostMass     /hostRadius      &
               &            +mergerRemnantSizeOrbitalEnergy*satelliteSpheroidMass*hostSpheroidMass/(satelliteRadius+hostRadius)&
               &                                                                                  /bindingEnergyFormFactor
          remnantRadius=(satelliteSpheroidMass+hostSpheroidMass)**2/progenitorsEnergy
       else
          remnantRadius=satelliteRadius
       end if

       ! Also compute the specific angular momentum at the half-mass radius.
       remnantCircularVelocity=dsqrt(gravitationalConstantGalacticus*(satelliteSpheroidMass+hostSpheroidMass)/remnantRadius)
       remnantSpecificAngularMomentum=remnantRadius*remnantCircularVelocity*darkMatterFactor
    end if
    return
  end subroutine Satellite_Merging_Remnant_Size_Cole2000

end module Satellite_Merging_Remnant_Sizes_Cole2000

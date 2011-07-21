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


!% Contains a module which calculates the specific angular momentum of cooling gas assuming a constant rotation velocity as a
!% function of radius.

module Cooling_Specific_Angular_Momenta_Constant_Rotation
  !% Calculates the specific angular momentum of cooling gas assuming a constant rotation velocity as a function of radius.
  use, intrinsic :: ISO_C_Binding
  use Kind_Numbers
  private
  public :: Cooling_Specific_AM_Constant_Rotation_Initialize, Cooling_Specific_AM_Constant_Rotation_Reset

  ! Record of unique ID of node which we last computed results for.
  integer(kind=kind_int8)    :: lastUniqueID=-1
  !$omp threadprivate(lastUniqueID)

  ! Record of whether or not specific angular momentum of cooling has has already been computed for this node.
  logical                    :: coolingSpecificAngularMomentumComputed=.false.
  !$omp threadprivate(coolingSpecificAngularMomentumComputed)

  ! Stored values of specific angular momentum of cooling gas.
  double precision           :: coolingSpecificAngularMomentumStored
  !$omp threadprivate(coolingSpecificAngularMomentumStored)

  ! Parameters controlling the calculation.
  integer,         parameter :: profileDarkMatter=0
  integer,         parameter :: profileHotGas    =1
  integer                    :: meanSpecificAngularMomentumFrom,rotationNormalizationFrom

contains

  !# <coolingSpecificAngularMomentumMethod>
  !#  <unitName>Cooling_Specific_AM_Constant_Rotation_Initialize</unitName>
  !# </coolingSpecificAngularMomentumMethod>
  subroutine Cooling_Specific_AM_Constant_Rotation_Initialize(coolingSpecificAngularMomentumMethod,Cooling_Specific_Angular_Momentum_Get)
    !% Initializes the ``constant rotation'' specific angular momentum of cooling gas module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: coolingSpecificAngularMomentumMethod
    procedure(double precision), pointer, intent(inout) :: Cooling_Specific_Angular_Momentum_Get
    type(varying_string)                                :: inputOption    

    if (coolingSpecificAngularMomentumMethod == 'constant rotation') then
       Cooling_Specific_Angular_Momentum_Get => Cooling_Specific_Angular_Momentum_Constant_Rotation

       !@ <inputParameter>
       !@   <name>coolingMeanAngularMomentumFrom</name>
       !@   <defaultValue>hot gas</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The component (``hot gas'' or ``dark matter'') from which the mean specific angular momentum should be computed for
       !@     calculations of cooling gas specific angular momentum.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingMeanAngularMomentumFrom',inputOption,defaultValue='hot gas')
       select case (char(inputOption))
       case ("dark matter")
          meanSpecificAngularMomentumFrom=profileDarkMatter
       case ("hot gas")
          meanSpecificAngularMomentumFrom=profileHotGas
       end select

       !@ <inputParameter>
       !@   <name>coolingRotationVelocityFrom</name>
       !@   <defaultValue>dark matter</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The component (``hot gas'' or ``dark matter'') from which the constant rotation speed should be computed for
       !@     calculations of cooling gas specific angular momentum.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingRotationVelocityFrom',inputOption,defaultValue='dark matter')
       select case (char(inputOption))
       case ("dark matter")
          rotationNormalizationFrom=profileDarkMatter
       case ("hot gas")
          rotationNormalizationFrom=profileHotGas
       end select

    end if
    return
  end subroutine Cooling_Specific_AM_Constant_Rotation_Initialize

  !# <calculationResetTask>
  !# <unitName>Cooling_Specific_AM_Constant_Rotation_Reset</unitName>
  !# </calculationResetTask>
  subroutine Cooling_Specific_AM_Constant_Rotation_Reset(thisNode)
    !% Reset the specific angular momentum of cooling gas calculation.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    coolingSpecificAngularMomentumComputed=.false.
    lastUniqueID                          =thisNode%uniqueID()
    return
  end subroutine Cooling_Specific_AM_Constant_Rotation_Reset

  double precision function Cooling_Specific_Angular_Momentum_Constant_Rotation(thisNode)
    !% Return the specific angular momentum of cooling gas in the constant rotation model.
    use Tree_Nodes
    use Cooling_Infall_Radii
    use Dark_Matter_Profiles
    use Hot_Halo_Density_Profile
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: meanSpecificAngularMomentum,rotationNormalization

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Cooling_Specific_AM_Constant_Rotation_Reset(thisNode)

    ! Check if specific angular momentum of cooling gas is already computed.
    if (.not.coolingSpecificAngularMomentumComputed) then
       ! Flag that cooling radius is now computed.
       coolingSpecificAngularMomentumComputed=.true.

       ! Compute the mean specific angular momentum.
       select case (meanSpecificAngularMomentumFrom)
       case (profileDarkMatter)
          ! Compute mean specific angular momentum of the dark matter halo from the spin parameter, mass and energy of the halo.
          meanSpecificAngularMomentum= gravitationalConstantGalacticus                          &
               &                      *           Tree_Node_Spin            (thisNode)          &
               &                      *           Tree_Node_Mass            (thisNode)  **1.5d0 &
               &                      /dsqrt(dabs(Dark_Matter_Profile_Energy(thisNode)))
       case (profileHotGas    )
          ! Compute mean specific angular momentum from the hot halo component.
          meanSpecificAngularMomentum= Tree_Node_Hot_Halo_Angular_Momentum(thisNode) &
               &                      /Tree_Node_Hot_Halo_Mass            (thisNode)
       end select

       ! Compute the rotation normalization.
       select case (rotationNormalizationFrom      )
       case (profileDarkMatter)
          rotationNormalization=Dark_Matter_Profile_Rotation_Normalization(thisNode)
       case (profileHotGas    )
          rotationNormalization=Hot_Halo_Profile_Rotation_Normalization   (thisNode)
       end select

       ! Compute the specific angular momentum of the cooling gas.
       coolingSpecificAngularMomentumStored= Infall_Radius(thisNode)     &
            &                               *rotationNormalization       &
            &                               *meanSpecificAngularMomentum
    end if

    ! Return the computed value.
    Cooling_Specific_Angular_Momentum_Constant_Rotation=coolingSpecificAngularMomentumStored
    return
  end function Cooling_Specific_Angular_Momentum_Constant_Rotation
  
end module Cooling_Specific_Angular_Momenta_Constant_Rotation

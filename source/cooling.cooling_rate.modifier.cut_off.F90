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

!% Contains a module which implements a cut off in the cooling rate at given redshift and virial velocity.

module Cooling_Rates_Modifier_Cut_Off
  !% Implements a cut off in the cooling rate at given redshift and virial velocity.
  implicit none
  private
  public :: Cooling_Rate_Modifier_Cut_Off

  ! Record of whether module has been initialized.
  logical          :: moduleInitialized=.false.

  ! Parameters controlling the time and velocity scale at which cooling is cut off.
  double precision :: coolingCutOffVelocity,coolingCutOffRedshift,coolingCutOffTime
  logical          :: coolingCutOffFormationNode
  integer          :: coolingCutOffWhen
  integer          :: coolingCutOffWhenBefore=0
  integer          :: coolingCutOffWhenAfter =1

contains

  !# <coolingRateModifierMethod>
  !#  <unitName>Cooling_Rate_Modifier_Cut_Off</unitName>
  !# </coolingRateModifierMethod>
  subroutine Cooling_Rate_Modifier_Cut_Off(thisNode,coolingRate)
    !% Modify cooling rates by truncating them to zero below a given redshift and virial velocity.
    use Input_Parameters
    use Cosmology_Functions
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(treeNode)      , intent(inout), pointer :: thisNode
    double precision    , intent(inout)          :: coolingRate
    double precision                             :: virialVelocity
    type(varying_string)                         :: coolingCutOffWhenText
    
    if (.not.moduleInitialized) then
       !$omp critical (Cooling_Rate_Modifier_Cut_Off_Initialize)
       if (.not.moduleInitialized) then
          !@ <inputParameter>
          !@   <name>coolingCutOffFormationNode</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies whether to use the virial velocity of the formation node or current node in the cooling rate ``cut-off'' modifier.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("coolingCutOffFormationNode",coolingCutOffFormationNode,defaultValue=.false.)
          !@ <inputParameter>
          !@   <name>coolingCutOffVelocity</name>
          !@   <defaultValue>0.0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The velocity below which cooling is suppressed in the ``cut-off'' cooling rate modifier method.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("coolingCutOffVelocity",coolingCutOffVelocity,defaultValue=0.0d0)
          !@ <inputParameter>
          !@   <name>coolingCutOffRedshift</name>
          !@   <defaultValue>0.0</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    The redshift below which cooling is suppressed in the ``cut-off'' cooling rate modifier method.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("coolingCutOffRedshift",coolingCutOffRedshift,defaultValue=0.0d0)
          coolingCutOffTime=Cosmology_Age(Expansion_Factor_from_Redshift(coolingCutOffRedshift))
          !@ <inputParameter>
          !@   <name>coolingCutOffWhen</name>
          !@   <defaultValue>after</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Specifies whether cooling is cut off before or after {\tt [coolingCutOffRedshift]}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("coolingCutOffWhen",coolingCutOffWhenText,defaultValue='after')
          select case (char(coolingCutOffWhenText))
          case ("before")
             coolingCutOffWhen=coolingCutOffWhenBefore
          case ("after" )
             coolingCutOffWhen=coolingCutOffWhenAfter
          case default
             call Galacticus_Error_Report('Cooling_Rate_Modifier_Cut_Off','[coolingCutOffWhen] must be either "before" or "after"')
          end select
          ! Record that the module is now initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical (Cooling_Rate_Modifier_Cut_Off_Initialize)
    end if
    
    ! Return immediately if cut-off is non-positive.
    if (coolingCutOffVelocity <= 0.0d0) return

    ! Test for halos where cooling should be cut off.
    select case (coolingCutOffFormationNode)
    case (.false.)
       virialVelocity=Dark_Matter_Halo_Virial_Velocity(thisNode              )
    case (.true. )
       virialVelocity=Dark_Matter_Halo_Virial_Velocity(thisNode%formationNode)
    end select
    if     (                                                                                                    &
         &  (                                                                                                   &
         &   (Tree_Node_Time(thisNode) >= coolingCutOffTime .and. coolingCutOffWhen == coolingCutOffWhenAfter ) &
         &    .or.                                                                                              &
         &   (Tree_Node_Time(thisNode) <= coolingCutOffTime .and. coolingCutOffWhen == coolingCutOffWhenBefore) &
         &  )                                                                                                   &
         &   .and.                                                                                              &
         &  virialVelocity             <= coolingCutOffVelocity                                                 &
         & ) coolingRate=0.0d0
    return
  end subroutine Cooling_Rate_Modifier_Cut_Off
  
end module Cooling_Rates_Modifier_Cut_Off

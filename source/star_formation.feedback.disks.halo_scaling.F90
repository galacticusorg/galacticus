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

!% Contains a module which implements an outflow rate due to star formation feedback in galactic disks that scales with halo
!% virial velocity and redshift.

module Star_Formation_Feedback_Disks_Halo_Scaling
  !% Implements an outflow rate due to star formation feedback in galactic disks that scales with halo
  !% virial velocity and redshift.
  use Galacticus_Nodes
  implicit none
  private
  public :: Star_Formation_Feedback_Disks_Halo_Scaling_Initialize

  ! Parameters of the feedback model.
  double precision :: diskOutflowFraction              , diskOutflowRedshiftExponent, &
       &              diskOutflowVirialVelocityExponent

  ! Normalization factor for the outflow rate.
  double precision :: outflowNormalization
  
contains

  !# <starFormationFeedbackDisksMethod>
  !#  <unitName>Star_Formation_Feedback_Disks_Halo_Scaling_Initialize</unitName>
  !# </starFormationFeedbackDisksMethod>
  subroutine Star_Formation_Feedback_Disks_Halo_Scaling_Initialize(starFormationFeedbackDisksMethod,Star_Formation_Feedback_Disk_Outflow_Rate_Get)
    !% Initializes the ``halo scaling'' disk star formation feedback module.
    use ISO_Varying_String
    use Input_Parameters
    use Stellar_Feedback
    implicit none
    type            (varying_string                                        ), intent(in   )          :: starFormationFeedbackDisksMethod              
    procedure       (Star_Formation_Feedback_Disk_Outflow_Rate_Halo_Scaling), intent(inout), pointer :: Star_Formation_Feedback_Disk_Outflow_Rate_Get 
    double precision                                                        , parameter              :: virialVelocityNormalization=200.0d0

    if (starFormationFeedbackDisksMethod == 'haloScaling') then
       Star_Formation_Feedback_Disk_Outflow_Rate_Get => Star_Formation_Feedback_Disk_Outflow_Rate_Halo_Scaling
       ! Get parameters of for the feedback calculation.
       !@ <inputParameter>
       !@   <name>diskOutflowFraction</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The ratio of outflow rate to star formation rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowFraction',diskOutflowFraction,defaultValue=0.01d0)
       !@ <inputParameter>
       !@   <name>diskOutflowVirialVelocityExponent</name>
       !@   <defaultValue>$-2.0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of virial velocity in the outflow rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowVirialVelocityExponent',diskOutflowVirialVelocityExponent,defaultValue=-2.0d0)
       !@ <inputParameter>
       !@   <name>diskOutflowRedshiftExponent</name>
       !@   <defaultValue>$0.0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of redshift in the outflow rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowRedshiftExponent',diskOutflowRedshiftExponent,defaultValue=0.0d0)
       ! Compute the normalization factor.
       outflowNormalization= diskOutflowFraction                                            &
            &               /feedbackEnergyInputAtInfinityCanonical                         &
            &               /virialVelocityNormalization**diskOutflowVirialVelocityExponent
    end if
    return
  end subroutine Star_Formation_Feedback_Disks_Halo_Scaling_Initialize

  double precision function Star_Formation_Feedback_Disk_Outflow_Rate_Halo_Scaling(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\tt thisNode}.
    use Cosmology_Functions
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode               ), intent(inout), pointer :: thisNode
    double precision                         , intent(in   )          :: energyInputRate                    , starFormationRate
    class           (nodeComponentBasic     )               , pointer :: thisBasicComponent
    double precision                         , save                   :: velocityPrevious           =-1.0d0, velocityFactorPrevious       =-1.0d0
    !$omp threadprivate(velocityPrevious,velocityFactorPrevious)
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctionsDefault
    double precision                         , save                   :: expansionFactorPrevious    =-1.0d0, expansionFactorFactorPrevious=-1.0d0
    !$omp threadprivate(expansionFactorPrevious,expansionFactorFactorPrevious)
    double precision                                                  :: expansionFactor                    , virialVelocity

    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Get the basic component.
    thisBasicComponent => thisNode%basic()

    ! Get virial velocity and expansion factor.
    virialVelocity =Dark_Matter_Halo_Virial_Velocity(thisNode                 )
    expansionFactor=cosmologyFunctionsDefault%expansionFactor                (thisBasicComponent%time())

    ! Compute the velocity factor.
    if (virialVelocity /= velocityPrevious) then
       velocityPrevious      =virialVelocity
       velocityFactorPrevious=virialVelocity**diskOutflowVirialVelocityExponent
    end if
    
    ! Compute the expansion-factor factor.
    if (expansionFactor /= expansionFactorPrevious) then
       expansionFactorPrevious      =      expansionFactor
       expansionFactorFactorPrevious=1.0d0/expansionFactor**diskOutflowRedshiftExponent
    end if

    ! Compute the outflow rate.
    Star_Formation_Feedback_Disk_Outflow_Rate_Halo_Scaling= &
         & outflowNormalization                             &
         & *energyInputRate                                 &
         & *velocityFactorPrevious                          &
         & *expansionFactorFactorPrevious
    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate_Halo_Scaling

end module Star_Formation_Feedback_Disks_Halo_Scaling

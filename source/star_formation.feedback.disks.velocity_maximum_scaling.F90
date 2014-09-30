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

!% Contains a module which implements an outflow rate due to star formation feedback in galactic disks that scales with halo
!% maximum velocity and redshift.

module Star_Formation_Feedback_Disks_VlctyMxSclng
  !% Implements an outflow rate due to star formation feedback in galactic disks that scales with halo
  !% maximum velocity and redshift.
  use Galacticus_Nodes
  implicit none
  private
  public :: Star_Formation_Feedback_Disks_VlctyMxSclng_Initialize

  ! Parameters of the feedback model.
  double precision :: diskOutflowFraction        , diskOutflowRedshiftExponent, &
       &              diskOutflowVelocityExponent

  ! Normalization factor for the outflow rate.
  double precision :: outflowNormalization
  
contains

  !# <starFormationFeedbackDisksMethod>
  !#  <unitName>Star_Formation_Feedback_Disks_VlctyMxSclng_Initialize</unitName>
  !# </starFormationFeedbackDisksMethod>
  subroutine Star_Formation_Feedback_Disks_VlctyMxSclng_Initialize(starFormationFeedbackDisksMethod,Star_Formation_Feedback_Disk_Outflow_Rate_Get)
    !% Initializes the ``halo scaling'' disk star formation feedback module.
    use ISO_Varying_String
    use Input_Parameters
    use Stellar_Feedback
    implicit none
    type            (varying_string                                        ), intent(in   )          :: starFormationFeedbackDisksMethod              
    procedure       (Star_Formation_Feedback_Disk_Outflow_Rate_VlctyMxSclng), intent(inout), pointer :: Star_Formation_Feedback_Disk_Outflow_Rate_Get 
    double precision                                                        , parameter              :: virialVelocityNormalization=200.0d0

    if (starFormationFeedbackDisksMethod == 'velocityMaximumScaling') then
       Star_Formation_Feedback_Disk_Outflow_Rate_Get => Star_Formation_Feedback_Disk_Outflow_Rate_VlctyMxSclng
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
       !@   <name>diskOutflowVelocityExponent</name>
       !@   <defaultValue>$-2.0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of virial velocity in the outflow rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowVelocityExponent',diskOutflowVelocityExponent,defaultValue=-2.0d0)
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
       outflowNormalization= diskOutflowFraction                                      &
            &               /feedbackEnergyInputAtInfinityCanonical                   &
            &               /virialVelocityNormalization**diskOutflowVelocityExponent
    end if
    return
  end subroutine Star_Formation_Feedback_Disks_VlctyMxSclng_Initialize

  double precision function Star_Formation_Feedback_Disk_Outflow_Rate_VlctyMxSclng(node,starFormationRate,energyInputRate)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\tt node}.
    use Cosmology_Functions
    use Dark_Matter_Profiles
    implicit none
    type            (treeNode               ), intent(inout), pointer :: node
    double precision                         , intent(in   )          :: energyInputRate                , starFormationRate
    class           (nodeComponentBasic     )               , pointer :: basic
    double precision                         , save                   :: velocityPrevious       =-1.0d0, velocityFactorPrevious       =-1.0d0
    !$omp threadprivate(velocityPrevious,velocityFactorPrevious)
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctions_
    class           (darkMatterProfileClass )               , pointer :: darkMatterProfile_
    double precision                         , save                   :: expansionFactorPrevious=-1.0d0, expansionFactorFactorPrevious=-1.0d0
    !$omp threadprivate(expansionFactorPrevious,expansionFactorFactorPrevious)
    double precision                                                  :: expansionFactor                , velocityMaximum

    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()
    darkMatterProfile_  => darkMatterProfile ()
    ! Get the basic component.
    basic               => node%basic        ()
    ! Get virial velocity and expansion factor.
    velocityMaximum=darkMatterProfile_ %circularVelocityMaximum(node        )
    expansionFactor=cosmologyFunctions_%expansionFactor        (basic%time())
    ! Compute the velocity factor.
    if (velocityMaximum /= velocityPrevious       ) then
       velocityPrevious             =      velocityMaximum
       velocityFactorPrevious       =      velocityMaximum**diskOutflowVelocityExponent
    end if
    ! Compute the expansion-factor factor.
    if (expansionFactor /= expansionFactorPrevious) then
       expansionFactorPrevious      =      expansionFactor
       expansionFactorFactorPrevious=1.0d0/expansionFactor**diskOutflowRedshiftExponent
    end if
    ! Compute the outflow rate.
    Star_Formation_Feedback_Disk_Outflow_Rate_VlctyMxSclng= &
         & outflowNormalization                             &
         & *energyInputRate                                 &
         & *velocityFactorPrevious                          &
         & *expansionFactorFactorPrevious
    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate_VlctyMxSclng

end module Star_Formation_Feedback_Disks_VlctyMxSclng

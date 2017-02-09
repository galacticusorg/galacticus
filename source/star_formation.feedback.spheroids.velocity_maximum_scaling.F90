!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module which implements an outflow rate due to star formation feedback in galactic spheroids that scales with halo
!% maximum velocity and redshift.

module Star_Formation_Feedback_Spheroids_VlctyMxSclng
  !% Implements an outflow rate due to star formation feedback in galactic spheroids that scales with halo
  !% maximum velocity and redshift.
  use Galacticus_Nodes
  implicit none
  private
  public :: Star_Formation_Feedback_Spheroids_VlctyMxSclng_Initialize

  ! Parameters of the feedback model.
  double precision :: spheroidOutflowFraction        , spheroidOutflowRedshiftExponent, &
       &              spheroidOutflowVelocityExponent

  ! Normalization factor for the outflow rate.
  double precision :: outflowNormalization
  
contains

  !# <starFormationFeedbackSpheroidsMethod>
  !#  <unitName>Star_Formation_Feedback_Spheroids_VlctyMxSclng_Initialize</unitName>
  !# </starFormationFeedbackSpheroidsMethod>
  subroutine Star_Formation_Feedback_Spheroids_VlctyMxSclng_Initialize(starFormationFeedbackSpheroidsMethod,Star_Formation_Feedback_Spheroid_Outflow_Rate_Get)
    !% Initializes the ``halo scaling'' spheroid star formation feedback module.
    use ISO_Varying_String
    use Input_Parameters
    use Stellar_Feedback
    implicit none
    type            (varying_string                                            ), intent(in   )          :: starFormationFeedbackSpheroidsMethod              
    procedure       (Star_Formation_Feedback_Spheroid_Outflow_Rate_VlctyMxSclng), intent(inout), pointer :: Star_Formation_Feedback_Spheroid_Outflow_Rate_Get 
    double precision                                                            , parameter              :: velocityNormalization=200.0d0

    if (starFormationFeedbackSpheroidsMethod == 'velocityMaximumScaling') then
       Star_Formation_Feedback_Spheroid_Outflow_Rate_Get => Star_Formation_Feedback_Spheroid_Outflow_Rate_VlctyMxSclng
       ! Get parameters of for the feedback calculation.
       !@ <inputParameter>
       !@   <name>spheroidOutflowFraction</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The ratio of outflow rate to star formation rate in spheroids.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidOutflowFraction',spheroidOutflowFraction,defaultValue=0.01d0)
       !@ <inputParameter>
       !@   <name>spheroidOutflowVelocityExponent</name>
       !@   <defaultValue>$-2.0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of maximum velocity in the outflow rate in spheroids.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidOutflowVelocityExponent',spheroidOutflowVelocityExponent,defaultValue=-2.0d0)
       !@ <inputParameter>
       !@   <name>spheroidOutflowRedshiftExponent</name>
       !@   <defaultValue>$0.0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of redshift in the outflow rate in spheroids.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidOutflowRedshiftExponent',spheroidOutflowRedshiftExponent,defaultValue=0.0d0)
       ! Compute the normalization factor.
       outflowNormalization=+spheroidOutflowFraction                                &
            &               /feedbackEnergyInputAtInfinityCanonical                 &
            &               /velocityNormalization**spheroidOutflowVelocityExponent
    end if
    return
  end subroutine Star_Formation_Feedback_Spheroids_VlctyMxSclng_Initialize

  double precision function Star_Formation_Feedback_Spheroid_Outflow_Rate_VlctyMxSclng(node,starFormationRate,energyInputRate)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic spheroid of {\normalfont \ttfamily node}.
    use Cosmology_Functions
    use Dark_Matter_Profiles
    implicit none
    type            (treeNode               ), intent(inout) :: node
    double precision                         , intent(in   ) :: energyInputRate                , starFormationRate
    class           (nodeComponentBasic     ), pointer       :: basic
    double precision                         , save          :: velocityPrevious       =-1.0d0, velocityFactorPrevious       =-1.0d0
    !$omp threadprivate(velocityPrevious,velocityFactorPrevious)
    class           (cosmologyFunctionsClass), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileClass ), pointer       :: darkMatterProfile_
    double precision                         , save          :: expansionFactorPrevious=-1.0d0, expansionFactorFactorPrevious=-1.0d0
    !$omp threadprivate(expansionFactorPrevious,expansionFactorFactorPrevious)
    double precision                                         :: expansionFactor                , velocityMaximum
    !GCC$ attributes unused :: starFormationRate
    
    ! Get the default cosmology functions object.
    cosmologyFunctions_ => cosmologyFunctions()
    darkMatterProfile_  => darkMatterProfile ()
    ! Get the basic component.
    basic               => node%basic        ()
    ! Get maximum velocity and expansion factor.
    velocityMaximum=darkMatterProfile_ %circularVelocityMaximum(node        )
    expansionFactor=cosmologyFunctions_%expansionFactor        (basic%time())
    ! Compute the velocity factor.
    if (velocityMaximum /= velocityPrevious       ) then
       velocityPrevious             =      velocityMaximum
       velocityFactorPrevious       =      velocityMaximum**spheroidOutflowVelocityExponent
    end if
    ! Compute the expansion-factor factor.
    if (expansionFactor /= expansionFactorPrevious) then
       expansionFactorPrevious      =      expansionFactor
       expansionFactorFactorPrevious=1.0d0/expansionFactor**spheroidOutflowRedshiftExponent
    end if
    ! Compute the outflow rate.
    Star_Formation_Feedback_Spheroid_Outflow_Rate_VlctyMxSclng= &
         & +outflowNormalization                                &
         & *energyInputRate                                     &
         & *velocityFactorPrevious                              &
         & *expansionFactorFactorPrevious
    return
  end function Star_Formation_Feedback_Spheroid_Outflow_Rate_VlctyMxSclng

end module Star_Formation_Feedback_Spheroids_VlctyMxSclng

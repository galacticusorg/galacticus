!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implementation of an outflow rate due to star formation feedback in galactic disks which scales with halo velocity.

  use Cosmology_Functions
  use Dark_Matter_Halo_Scales
  
  !# <starFormationFeedbackDisks name="starFormationFeedbackDisksHaloScaling" defaultThreadPrivate="yes">
  !#  <description>An outflow rate due to star formation feedback in galactic disks which scales with halo velocity.</description>
  !# </starFormationFeedbackDisks>
  type, extends(starFormationFeedbackDisksClass) :: starFormationFeedbackDisksHaloScaling
     !% Implementation of an outflow rate due to star formation feedback in galactic disks which scales with halo velocity.
     private
     double precision                                    :: fraction               , exponentRedshift     , &
          &                                                 exponentVelocity       , normalization        , &
          &                                                 velocityPrevious       , velocityFactor       , &
          &                                                 expansionFactorPrevious, expansionFactorFactor
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
   contains
     final     ::                haloScalingDestructor
     procedure :: outflowRate => haloScalingOutflowRate
  end type starFormationFeedbackDisksHaloScaling

  interface starFormationFeedbackDisksHaloScaling
     !% Constructors for the halo scaling fraction star formation feedback in disks class.
     module procedure haloScalingConstructorParameters
     module procedure haloScalingConstructorInternal
  end interface starFormationFeedbackDisksHaloScaling

contains

  function haloScalingConstructorParameters(parameters) result(self)
    !% Constructor for the halo scaling fraction star formation feedback in disks class which takes a parameter set as input.
    use Galacticus_Error
    implicit none
    type            (starFormationFeedbackDisksHaloScaling)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    double precision                                                       :: fraction            , exponentRedshift, &
         &                                                                    exponentVelocity
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_

    !# <inputParameter>
    !#   <name>fraction</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The ratio of outflow rate to star formation rate in disks.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentVelocity</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-2.0d0</defaultValue>
    !#   <description>The exponent of virial velocity in the outflow rate in disks.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentRedshift</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The exponent of redshift in the outflow rate in disks.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=starFormationFeedbackDisksHaloScaling(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function haloScalingConstructorParameters

  function haloScalingConstructorInternal(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !% Internal constructor for the halo scaling star formation feedback from disks class.
    use Stellar_Feedback
    implicit none
    type            (starFormationFeedbackDisksHaloScaling)                        :: self
    double precision                                       , intent(in   )         :: fraction                     , exponentRedshift, &
         &                                                                            exponentVelocity
    class           (cosmologyFunctionsClass              ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass             ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                       , parameter             :: velocityNormalization=200.0d0

    !# <constructorAssign variables="fraction, exponentRedshift, exponentVelocity, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    ! Initialize stored values.
    self%velocityPrevious       =-1.0d0
    self%expansionFactorPrevious=-1.0d0
    self%velocityFactor         =-1.0d0
    self%expansionFactorFactor  =-1.0d0
    ! Compute the normalization factor.
    self%normalization          =+self%fraction                                &
         &                       /feedbackEnergyInputAtInfinityCanonical       &
         &                       /velocityNormalization**self%exponentVelocity
    return
  end function haloScalingConstructorInternal

  subroutine haloScalingDestructor(self)
    !% Destructor for the halo scaling feedback from star formation in disks class.
    implicit none
    type(starFormationFeedbackDisksHaloScaling), intent(inout) :: self
  
    !# <objectDestructor name="self%cosmologyFunctions_"  />
    !# <objectDestructor name="self%darkMatterHaloScale_" />
    return
  end subroutine haloScalingDestructor
  
  double precision function haloScalingOutflowRate(self,node,rateEnergyInput,rateStarFormation)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\normalfont \ttfamily node}.
    use Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (starFormationFeedbackDisksHaloScaling), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    double precision                                       , intent(in   ) :: rateEnergyInput, rateStarFormation
    class           (nodeComponentBasic                   ), pointer       :: basic
    double precision                                                       :: expansionFactor, velocityVirial
    !GCC$ attributes unused :: rateStarFormation

    ! Get the basic component.
    basic => node%basic()
    ! Get virial velocity and expansion factor.
    velocityVirial =self%darkMatterHaloScale_%virialVelocity (node        )
    expansionFactor=self%cosmologyFunctions_ %expansionFactor(basic%time())
    ! Compute the velocity factor.
    if (velocityVirial /= self%velocityPrevious) then
       self%velocityPrevious=velocityVirial
       self%velocityFactor  =velocityVirial**self%exponentVelocity
    end if    
    ! Compute the expansion-factor factor.
    if (expansionFactor /= self%expansionFactorPrevious) then
       self%expansionFactorPrevious=      expansionFactor
       self%expansionFactorFactor  =1.0d0/expansionFactor**self%exponentRedshift
    end if
    ! Compute the outflow rate.
    haloScalingOutflowRate=+self%normalization         &
         &                 *rateEnergyInput            &
         &                 *self%velocityFactor        &
         &                 *self%expansionFactorFactor
    return
  end function haloScalingOutflowRate

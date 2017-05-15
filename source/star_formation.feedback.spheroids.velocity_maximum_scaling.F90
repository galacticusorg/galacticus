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

  !% Implementation of an outflow rate due to star formation feedback in galactic spheroids which scales with peak halo velocity.

  use Math_Exponentiation
  use Cosmology_Functions
  
  !# <starFormationFeedbackSpheroids name="starFormationFeedbackSpheroidsVlctyMxSclng" defaultThreadPrivate="yes">
  !#  <description>An outflow rate due to star formation feedback in galactic spheroids which scales with peak halo velocity.</description>
  !# </starFormationFeedbackSpheroids>
  type, extends(starFormationFeedbackSpheroidsClass) :: starFormationFeedbackSpheroidsVlctyMxSclng
     !% Implementation of an outflow rate due to star formation feedback in galactic spheroids which scales with peak halo velocity.
     private
     double precision                                    :: fraction               , exponentRedshift            , &
          &                                                 exponentVelocity       , normalization               , &
          &                                                 velocityPrevious       , velocityFactor              , &
          &                                                 expansionFactorPrevious, expansionFactorFactor
     type            (fastExponentiator       )          :: velocityExponentiator  , expansionFactorExponentiator
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
   contains
     final     ::                vlctyMxSclngDestructor
     procedure :: outflowRate => vlctyMxSclngOutflowRate
  end type starFormationFeedbackSpheroidsVlctyMxSclng

  interface starFormationFeedbackSpheroidsVlctyMxSclng
     !% Constructors for the velocity maximum scaling fraction star formation feedback in spheroids class.
     module procedure vlctyMxSclngConstructorParameters
     module procedure vlctyMxSclngConstructorInternal
  end interface starFormationFeedbackSpheroidsVlctyMxSclng

contains

  function vlctyMxSclngConstructorParameters(parameters) result(self)
    !% Constructor for the velocity maximum scaling fraction star formation feedback in spheroids class which takes a parameter set as
    !% input.
    use Galacticus_Error
    implicit none
    type            (starFormationFeedbackSpheroidsVlctyMxSclng)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    double precision                                                            :: fraction            , exponentRedshift, &
         &                                                                         exponentVelocity
    class           (cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>fraction</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The ratio of outflow rate to star formation rate in spheroids.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentVelocity</name>
    !#   <source>parameters</source>
    !#   <defaultValue>-2.0d0</defaultValue>
    !#   <description>The exponent of virial velocity in the outflow rate in spheroids.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponentRedshift</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The exponent of redshift in the outflow rate in spheroids.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    self=starFormationFeedbackSpheroidsVlctyMxSclng(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_)
    return
  end function vlctyMxSclngConstructorParameters

  function vlctyMxSclngConstructorInternal(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_) result(self)
    !% Internal constructor for the halo scaling star formation feedback from spheroids class.
    use Stellar_Feedback
    implicit none
    type            (starFormationFeedbackSpheroidsVlctyMxSclng)                        :: self
    double precision                                            , intent(in   )         :: fraction                     , exponentRedshift, &
         &                                                                                 exponentVelocity
    class           (cosmologyFunctionsClass                   ), intent(in   ), target :: cosmologyFunctions_
    double precision                                            , parameter             :: velocityNormalization=200.0d0

    !# <constructorAssign variables="fraction, exponentRedshift, exponentVelocity, *cosmologyFunctions_"/>
    ! Initialize stored values.
    self%velocityPrevious       =-1.0d0
    self%expansionFactorPrevious=-1.0d0
    self%velocityFactor         =-1.0d0
    self%expansionFactorFactor  =-1.0d0
    ! Compute the normalization factor.
    self%normalization          =+self%fraction                                &
         &                       /feedbackEnergyInputAtInfinityCanonical       &
         &                       /velocityNormalization**self%exponentVelocity
    ! Initialize fast exponentiators.
    self%velocityExponentiator       =fastExponentiator(1.0d+0,1.0d+3,self%exponentVelocity,1.0d+1,abortOutsideRange=.false.)
    self%expansionFactorExponentiator=fastExponentiator(1.0d-3,1.0d+0,self%exponentRedshift,1.0d+3,abortOutsideRange=.false.)
    return
  end function vlctyMxSclngConstructorInternal

  subroutine vlctyMxSclngDestructor(self)
    !% Destructor for the velocity maximum scaling feedback from star formation in spheroids class.
    implicit none
    type(starFormationFeedbackSpheroidsVlctyMxSclng), intent(inout) :: self
  
    !# <objectDestructor name="self%cosmologyFunctions_" />
    return
  end subroutine vlctyMxSclngDestructor
  
  double precision function vlctyMxSclngOutflowRate(self,node,rateEnergyInput,rateStarFormation)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic spheroid of {\normalfont \ttfamily node}.
    use Dark_Matter_Profiles
    implicit none
    class           (starFormationFeedbackSpheroidsVlctyMxSclng), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    double precision                                            , intent(in   ) :: rateEnergyInput   , rateStarFormation
    class           (nodeComponentBasic                        ), pointer       :: basic
    class           (darkMatterProfileClass                    ), pointer       :: darkMatterProfile_
    double precision                                                            :: expansionFactor   , velocityMaximum
    !GCC$ attributes unused :: rateStarFormation

    ! Get required objects.
    darkMatterProfile_ => darkMatterProfile      ()
    ! Get the basic component.
    basic              => node             %basic()
    ! Get virial velocity and expansion factor.
    velocityMaximum=     darkMatterProfile_ %circularVelocityMaximum(node        )
    expansionFactor=self%cosmologyFunctions_%expansionFactor        (basic%time())
    ! Compute the velocity factor.
    if (velocityMaximum /= self%velocityPrevious) then
       self%velocityPrevious       =                                                     velocityMaximum
       self%velocityFactor         =      self%velocityExponentiator       %exponentiate(velocityMaximum)
    end if    
    ! Compute the expansion-factor factor.
    if (expansionFactor /= self%expansionFactorPrevious) then
       self%expansionFactorPrevious=                                                     expansionFactor
       self%expansionFactorFactor  =1.0d0/self%expansionFactorExponentiator%exponentiate(expansionFactor)
    end if
    ! Compute the outflow rate.
    vlctyMxSclngOutflowRate=+self%normalization         &
         &                  *rateEnergyInput            &
         &                  *self%velocityFactor        &
         &                  *self%expansionFactorFactor
    return
  end function vlctyMxSclngOutflowRate

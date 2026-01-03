!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  Implementation of a timescale for star formation which scales with the circular velocity of the host halo.
  !!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Kind_Numbers            , only : kind_int8
  use :: Math_Exponentiation     , only : fastExponentiator

  !![
  <starFormationTimescale name="starFormationTimescaleVelocityMaxScaling">
   <description>A velocityMaxScaling timescale for star formation.</description>
  </starFormationTimescale>
  !!]
  type, extends(starFormationTimescaleClass) :: starFormationTimescaleVelocityMaxScaling
     !!{
     Implementation of a velocityMaxScaling timescale for star formation.
     !!}
     private
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_           => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_         => null()
     double precision                                     :: expansionFactorFactorPrevious          , exponentVelocity            , &
          &                                                  exponentRedshift                       , timescaleNormalization      , &
          &                                                  timescaleStored                        , velocityMaximumPrevious     , &
          &                                                  velocityFactorPrevious                 , expansionFactorPrevious     , &
          &                                                  timeScale_
     logical                                              :: timescaleComputed
     integer         (kind_int8                )          :: lastUniqueID
     type            (fastExponentiator        )          :: velocityExponentiator                  , expansionFactorExponentiator
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                     velocityMaxScalingDestructor
     procedure :: autoHook         => velocityMaxScalingAutoHook
     procedure :: timescale        => velocityMaxScalingTimescale
     procedure :: calculationReset => velocityMaxScalingCalculationReset
  end type starFormationTimescaleVelocityMaxScaling

  interface starFormationTimescaleVelocityMaxScaling
     !!{
     Constructors for the \refClass{starFormationTimescaleVelocityMaxScaling} timescale for star formation class.
     !!}
     module procedure velocityMaxScalingConstructorParameters
     module procedure velocityMaxScalingConstructorInternal
  end interface starFormationTimescaleVelocityMaxScaling

  double precision, parameter :: velocityNormalization=200.0d0

contains

  function velocityMaxScalingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationTimescaleVelocityMaxScaling} timescale for star formation class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationTimescaleVelocityMaxScaling)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass               ), pointer       :: darkMatterProfileDMO_
    double precision                                                          :: timescale            , exponentVelocity, &
         &                                                                       exponentRedshift

    ! Get parameters of for the timescale calculation.
    !![
    <inputParameter>
      <name>timescale</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The timescale for star formation in the velocity maximum scaling timescale model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponentVelocity</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of virial velocity in the timescale for star formation in the velocity maximum scaling timescale model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponentRedshift</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of redshift in the timescale for star formation in the velocity maximum scaling timescale model.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=starFormationTimescaleVelocityMaxScaling(timescale,exponentVelocity,exponentRedshift,cosmologyFunctions_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function velocityMaxScalingConstructorParameters

  function velocityMaxScalingConstructorInternal(timescale,exponentVelocity,exponentRedshift,cosmologyFunctions_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{starFormationTimescaleVelocityMaxScaling} timescale for star formation class.
    !!}
    implicit none
    type            (starFormationTimescaleVelocityMaxScaling)                        :: self
    double precision                                          , intent(in   )         :: timescale            , exponentVelocity, &
         &                                                                               exponentRedshift
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass               ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="exponentVelocity, exponentRedshift, *cosmologyFunctions_, *darkMatterProfileDMO_"/>
    !!]

    self%lastUniqueID                 =-1_kind_int8
    self%timescaleComputed            =.false.
    self%velocityMaximumPrevious      =-1.0d0
    self%velocityFactorPrevious       =-1.0d0
    self%expansionFactorPrevious      =-1.0d0
    self%expansionFactorFactorPrevious=-1.0d0
    ! Compute the normalization of the timescale.
    self%timeScale_            =+timescale
    self%timeScaleNormalization=+timescale                                    &
         &                      /velocityNormalization**self%exponentVelocity
    ! Initialize exponentiators.
    self%velocityExponentiator       =fastExponentiator(1.0d+0,1.0d+3,self%exponentVelocity,1.0d+1,abortOutsideRange=.false.)
    self%expansionFactorExponentiator=fastExponentiator(1.0d-3,1.0d+0,self%exponentRedshift,1.0d+3,abortOutsideRange=.false.)
    return
  end function velocityMaxScalingConstructorInternal

  subroutine velocityMaxScalingAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationTimescaleVelocityMaxScaling), intent(inout) :: self

    call calculationResetEvent%attach(self,velocityMaxScalingCalculationReset,openMPThreadBindingAllLevels,label='starFormationTimescaleVelocityMaxScaling')
    return
  end subroutine velocityMaxScalingAutoHook

  subroutine velocityMaxScalingDestructor(self)
    !!{
    Destructor for the \refClass{starFormationTimescaleVelocityMaxScaling} timescale for star formation class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(starFormationTimescaleVelocityMaxScaling), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    if (calculationResetEvent%isAttached(self,velocityMaxScalingCalculationReset)) call calculationResetEvent%detach(self,velocityMaxScalingCalculationReset)
    return
  end subroutine velocityMaxScalingDestructor

  subroutine velocityMaxScalingCalculationReset(self,node,uniqueID)
    !!{
    Reset the velocity maximum scaling star formation timescale calculation.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    class  (starFormationTimescaleVelocityMaxScaling), intent(inout) :: self
    type   (treeNode                                ), intent(inout) :: node
    integer(kind_int8                               ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%timescaleComputed=.false.
    self%lastUniqueID     =uniqueID
    return
  end subroutine velocityMaxScalingCalculationReset

  double precision function velocityMaxScalingTimescale(self,component)
    !!{
    Returns the timescale (in Gyr) for star formation in the {\normalfont \ttfamily component} in the velocity maximum scaling timescale model.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: Mass_Distributions, only : massDistributionClass
   implicit none
    class           (starFormationTimescaleVelocityMaxScaling), intent(inout) :: self
    class           (nodeComponent                           ), intent(inout) :: component
    class           (nodeComponentBasic                      ), pointer       :: basic
    class           (massDistributionClass                   ), pointer       :: massDistribution_
    double precision                                                          :: expansionFactor, velocityMaximum

    ! Check if node differs from previous one for which we performed calculations.
    if (component%hostNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(component%hostNode,component%hostNode%uniqueID())
    ! Compute the timescale if necessary.
    if (.not.self%timescaleComputed) then
       ! Get virial velocity and expansion factor.
       massDistribution_ => self              %darkMatterProfileDMO_%get                         (component%hostNode  )
       basic             => component%hostNode%basic                                             (                    )
       velocityMaximum   =  massDistribution_                       %velocityRotationCurveMaximum(                    )
       expansionFactor   =  self              %cosmologyFunctions_  %expansionFactor             (basic    %time    ())
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
       ! Compute the velocity factor.
       if (velocityMaximum /= self%velocityMaximumPrevious) then
           self%velocityMaximumPrevious=velocityMaximum
           self%velocityFactorPrevious =velocityMaximum**self%exponentVelocity
       end if
       ! Compute the expansion-factor factor.
       if (expansionFactor /= self%expansionFactorPrevious) then
          self%expansionFactorPrevious      =      expansionFactor
          self%expansionFactorFactorPrevious=1.0d0/expansionFactor**self%exponentRedshift
       end if
       ! Computed the timescale.
       self%timescaleStored=+self%timeScaleNormalization        &
            &               *self%velocityFactorPrevious        &
            &               *self%expansionFactorFactorPrevious
       ! Record that the timescale is now computed.
       self%timescaleComputed=.true.
    end if
    ! Return the stored timescale.
    velocityMaxScalingTimescale=self%timescaleStored
    return
  end function velocityMaxScalingTimescale

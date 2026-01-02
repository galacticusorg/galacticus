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
  Implementation of an stellar feedback model which scales with halo velocity.
  !!}

  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsHaloScaling">
   <description>An stellar feedback model which scales with halo velocity.</description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsHaloScaling
     !!{
     Implementation of an stellar feedback model which scales with halo velocity.
     !!}
     private
     double precision                                    :: fraction                         , exponentRedshift     , &
          &                                                 exponentVelocity                 , normalization        , &
          &                                                 velocityPrevious                 , velocityFactor       , &
          &                                                 expansionFactorPrevious          , expansionFactorFactor
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_     => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_    => null()
   contains
     !![
     <methods>
       <method method="node" description="Get the node from which to compute halo properties."/>
     </methods>
     !!]
     final     ::                haloScalingDestructor
     procedure :: outflowRate => haloScalingOutflowRate
     procedure :: node        => haloScalingNode
  end type stellarFeedbackOutflowsHaloScaling

  interface stellarFeedbackOutflowsHaloScaling
     !!{
     Constructors for the halo scaling fraction stellar feedback class.
     !!}
     module procedure haloScalingConstructorParameters
     module procedure haloScalingConstructorInternal
  end interface stellarFeedbackOutflowsHaloScaling

contains

  function haloScalingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the halo scaling fraction stellar feedback class which takes a parameter set as input.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (stellarFeedbackOutflowsHaloScaling)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    double precision                                                    :: fraction            , exponentRedshift, &
         &                                                                 exponentVelocity
    class           (cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_

    !![
    <inputParameter>
      <name>fraction</name>
      <source>parameters</source>
      <defaultValue>0.01d0</defaultValue>
      <description>The ratio of outflow rate to star formation rate in disks.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentVelocity</name>
      <source>parameters</source>
      <defaultValue>-2.0d0</defaultValue>
      <description>The exponent of virial velocity in the outflow rate in disks.</description>
    </inputParameter>
    <inputParameter>
      <name>exponentRedshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The exponent of redshift in the outflow rate in disks.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=stellarFeedbackOutflowsHaloScaling(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function haloScalingConstructorParameters

  function haloScalingConstructorInternal(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the halo scaling stellar feedback class.
    !!}
    use :: Stellar_Feedback, only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    type            (stellarFeedbackOutflowsHaloScaling)                        :: self
    double precision                                    , intent(in   )         :: fraction                     , exponentRedshift, &
         &                                                                         exponentVelocity
    class           (cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                    , parameter             :: velocityNormalization=200.0d0

    !![
    <constructorAssign variables="fraction, exponentRedshift, exponentVelocity, *cosmologyFunctions_, *darkMatterHaloScale_"/>
    !!]
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
    !!{
    Destructor for the halo scaling stellar feedback class.
    !!}
    implicit none
    type(stellarFeedbackOutflowsHaloScaling), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine haloScalingDestructor

  subroutine haloScalingOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the given {\normalfont \ttfamily component}.
    !!}
    use :: Galacticus_Nodes, only : treeNode, nodeComponentBasic
    implicit none
    class           (stellarFeedbackOutflowsHaloScaling), intent(inout) :: self
    class           (nodeComponent                     ), intent(inout) :: component
    double precision                                    , intent(in   ) :: rateEnergyInput    , rateStarFormation
    double precision                                    , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive
    type            (treeNode                          ), pointer       :: node
    class           (nodeComponentBasic                ), pointer       :: basic
    double precision                                                    :: expansionFactor    , velocityVirial
    !$GLC attributes unused :: rateStarFormation

    ! Get the basic component.
    node  => self%node (component)
    basic => node%basic(         )
    ! Get virial velocity and expansion factor.
    velocityVirial =self%darkMatterHaloScale_%velocityVirial (component%hostNode  )
    expansionFactor=self%cosmologyFunctions_ %expansionFactor(basic    %time    ())
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
    rateOutflowEjective =+self%normalization         &
         &               *rateEnergyInput            &
         &               *self%velocityFactor        &
         &               *self%expansionFactorFactor
    rateOutflowExpulsive=+0.0d0
    return
  end subroutine haloScalingOutflowRate

  function haloScalingNode(self,component) result(node)
    !!{
    Returns a pointer to the node from which to extract halo properties.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type (treeNode                          ), pointer       :: node
    class(stellarFeedbackOutflowsHaloScaling), intent(inout) :: self
    class(nodeComponent                     ), intent(in   ) :: component
    !$GLC attributes unused :: self

    node => component%hostNode
    return
  end function haloScalingNode

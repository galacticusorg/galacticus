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
  Implementation of an stellar feedback model which scales with peak halo velocity.
  !!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Math_Exponentiation     , only : fastExponentiator

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsVlctyMxSclng">
   <description>An stellar feedback model which scales with peak halo velocity.</description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsVlctyMxSclng
     !!{
     Implementation of an stellar feedback model which scales with peak halo velocity.
     !!}
     private
     double precision                                     :: fraction                         , exponentRedshift            , &
          &                                                  exponentVelocity                 , normalization               , &
          &                                                  velocityPrevious                 , velocityFactor              , &
          &                                                  expansionFactorPrevious          , expansionFactorFactor
     type            (fastExponentiator        )          :: velocityExponentiator            , expansionFactorExponentiator
     class           (cosmologyFunctionsClass  ), pointer :: cosmologyFunctions_     => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_   => null()
   contains
     final     ::                vlctyMxSclngDestructor
     procedure :: outflowRate => vlctyMxSclngOutflowRate
  end type stellarFeedbackOutflowsVlctyMxSclng

  interface stellarFeedbackOutflowsVlctyMxSclng
     !!{
     Constructors for the velocity maximum scaling fraction stellar feedback class.
     !!}
     module procedure vlctyMxSclngConstructorParameters
     module procedure vlctyMxSclngConstructorInternal
  end interface stellarFeedbackOutflowsVlctyMxSclng

contains

  function vlctyMxSclngConstructorParameters(parameters) result(self)
    !!{
    Constructor for the velocity maximum scaling fraction stellar feedback class which takes a parameter set as
    input.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (stellarFeedbackOutflowsVlctyMxSclng)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: fraction             , exponentRedshift, &
         &                                                                  exponentVelocity
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass          ), pointer       :: darkMatterProfileDMO_

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
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=stellarFeedbackOutflowsVlctyMxSclng(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"  />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function vlctyMxSclngConstructorParameters

  function vlctyMxSclngConstructorInternal(fraction,exponentRedshift,exponentVelocity,cosmologyFunctions_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the halo scaling stellar feedback class.
    !!}
    use :: Stellar_Feedback, only : feedbackEnergyInputAtInfinityCanonical
    implicit none
    type            (stellarFeedbackOutflowsVlctyMxSclng)                        :: self
    double precision                                     , intent(in   )         :: fraction                     , exponentRedshift, &
         &                                                                          exponentVelocity
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterProfileDMOClass          ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                     , parameter             :: velocityNormalization=200.0d0

    !![
    <constructorAssign variables="fraction, exponentRedshift, exponentVelocity, *cosmologyFunctions_, *darkMatterProfileDMO_"/>
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
    ! Initialize fast exponentiators.
    self%velocityExponentiator       =fastExponentiator(1.0d+0,1.0d+3,self%exponentVelocity,1.0d+1,abortOutsideRange=.false.)
    self%expansionFactorExponentiator=fastExponentiator(1.0d-3,1.0d+0,self%exponentRedshift,1.0d+3,abortOutsideRange=.false.)
    return
  end function vlctyMxSclngConstructorInternal

  subroutine vlctyMxSclngDestructor(self)
    !!{
    Destructor for the velocity maximum scaling stellar feedback class.
    !!}
    implicit none
    type(stellarFeedbackOutflowsVlctyMxSclng), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"  />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine vlctyMxSclngDestructor

  subroutine vlctyMxSclngOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the given {\normalfont \ttfamily component}.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (stellarFeedbackOutflowsVlctyMxSclng), intent(inout) :: self
    class           (nodeComponent                      ), intent(inout) :: component
    double precision                                     , intent(in   ) :: rateEnergyInput    , rateStarFormation
    double precision                                     , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive
    class           (nodeComponentBasic                 ), pointer       :: basic
    class           (massDistributionClass              ), pointer       :: massDistribution_
    double precision                                                     :: expansionFactor    , velocityMaximum
    !$GLC attributes unused :: rateStarFormation

    ! Get the basic component.
    basic => component%hostNode%basic()
    ! Get virial velocity and expansion factor.
    massDistribution_ => self             %darkMatterProfileDMO_%get                         (component%hostNode  )
    velocityMaximum   =  massDistribution_                      %velocityRotationCurveMaximum(                    )
    expansionFactor   =  self             %cosmologyFunctions_  %expansionFactor             (basic    %time    ())
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
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
    rateOutflowEjective =+self%normalization         &
         &               *rateEnergyInput            &
         &               *self%velocityFactor        &
         &               *self%expansionFactorFactor
    rateOutflowExpulsive=+0.0d0
    return
  end subroutine vlctyMxSclngOutflowRate

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a black hole CGM heating class where the coupling strength is determined by how quasistatic the \gls{cgm} is.
  !!}

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use :: Cooling_Radii             , only : coolingRadiusClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass

  !![
  <blackHoleCGMHeating name="blackHoleCGMHeatingQuasistatic">
   <description>
    A black hole CGM heating class where the coupling strength is determined by how quasistatic the \gls{cgm} is.
   </description>
  </blackHoleCGMHeating>
  !!]
  type, extends(blackHoleCGMHeatingClass) :: blackHoleCGMHeatingQuasistatic
     !!{
     A black hole CGM heating class where the coupling strength is determined by how quasistatic the \gls{cgm} is.
     !!}
     private
     class           (blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
     class           (darkMatterHaloScaleClass   ), pointer :: darkMatterHaloScale_    => null()
     class           (coolingRadiusClass         ), pointer :: coolingRadius_          => null()
     double precision                                       :: efficiencyHeating
   contains
     final     ::                quasistaticDestructor
     procedure :: heatingRate => quasistaticHeatingRate
  end type blackHoleCGMHeatingQuasistatic
  
  interface blackHoleCGMHeatingQuasistatic
     !!{
     Constructors for the {\normalfont \ttfamily quasistatic} black hole winds class.
     !!}
     module procedure quasistaticConstructorParameters
     module procedure quasistaticConstructorInternal
  end interface blackHoleCGMHeatingQuasistatic

contains

  function quasistaticConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily quasistatic} black hole winds class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleCGMHeatingQuasistatic)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (blackHoleAccretionRateClass   ), pointer       :: blackHoleAccretionRate_
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class           (coolingRadiusClass            ), pointer       :: coolingRadius_
    double precision                                                :: efficiencyHeating
    
    !![
    <inputParameter>
      <name>efficiencyHeating</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The efficiency with which accretion onto a black hole heats the \gls{cgm}.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"    name="darkMatterHaloScale_"    source="parameters"/>
    <objectBuilder class="coolingRadius"          name="coolingRadius_"          source="parameters"/>
    !!]
    self=blackHoleCGMHeatingQuasistatic(efficiencyHeating,blackHoleAccretionRate_,darkMatterHaloScale_,coolingRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="darkMatterHaloScale_"   />
    <objectDestructor name="coolingRadius_"         />
    !!]
    return
  end function quasistaticConstructorParameters

  function quasistaticConstructorInternal(efficiencyHeating,blackHoleAccretionRate_,darkMatterHaloScale_,coolingRadius_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily quasistatic} node operator class.
    !!}
    implicit none
    type            (blackHoleCGMHeatingQuasistatic)                        :: self
    class           (blackHoleAccretionRateClass   ), target, intent(in   ) :: blackHoleAccretionRate_
    class           (darkMatterHaloScaleClass      ), target, intent(in   ) :: darkMatterHaloScale_
    class           (coolingRadiusClass            ), target, intent(in   ) :: coolingRadius_
    double precision                                        , intent(in   ) :: efficiencyHeating
    !![
    <constructorAssign variables="efficiencyHeating, *blackHoleAccretionRate_, *darkMatterHaloScale_, *coolingRadius_"/>
    !!]

    return
  end function quasistaticConstructorInternal

  subroutine quasistaticDestructor(self)
    !!{
    Destructor for the quasistatic black hole CGM heating class.
    !!}
    implicit none
    type(blackHoleCGMHeatingQuasistatic), intent(inout) :: self
    
    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    <objectDestructor name="self%darkMatterHaloScale_"   />
    <objectDestructor name="self%coolingRadius_"         />
    !!]
    return
  end subroutine quasistaticDestructor
  
  double precision function quasistaticHeatingRate(self,blackHole) result(rateHeating)
    !!{
    Compute the heating rate of the CGM based on the accretion disk jet power.
    !!}
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class           (blackHoleCGMHeatingQuasistatic), intent(inout) :: self
    class           (nodeComponentBlackHole        ), intent(inout) :: blackHole
    double precision                                , parameter     :: radiusCoolingFractionalTransitionMinimum=0.9d0
    double precision                                , parameter     :: radiusCoolingFractionalTransitionMaximum=1.0d0
    double precision                                                :: rateAccretionSpheroid                         , rateAccretionHotHalo           , &
         &                                                             rateAccretion                                 , rateAccretionNuclearStarCluster, &
         &                                                             efficiencyCoupling                            , x                              , &
         &                                                             radiusCoolingFractional

    ! Compute accretion rate onto the black hole.
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
    rateAccretion=+rateAccretionSpheroid           &
         &        +rateAccretionHotHalo            &
         &        +rateAccretionNuclearStarCluster
    ! No heating for non-positive accretion rates.
    if (rateAccretion > 0.0d0) then
       ! Compute coupling efficiency based on whether halo is cooling quasistatically.
       radiusCoolingFractional=+self%coolingRadius_      %      radius(blackHole%hostNode) &
            &                  /self%darkMatterHaloScale_%radiusVirial(blackHole%hostNode)
       if      (radiusCoolingFractional < radiusCoolingFractionalTransitionMinimum) then
          efficiencyCoupling=1.0d0
       else if (radiusCoolingFractional > radiusCoolingFractionalTransitionMaximum) then
          efficiencyCoupling=0.0d0
       else
          x                 =     +(radiusCoolingFractional                 -radiusCoolingFractionalTransitionMinimum) &
               &                  /(radiusCoolingFractionalTransitionMaximum-radiusCoolingFractionalTransitionMinimum)
          efficiencyCoupling=x**2*(2.0d0*x-3.0d0)+1.0d0
       end if
       ! Compute the heating rate.
       rateHeating=+     efficiencyCoupling    &
            &      *self%efficiencyHeating     &
            &      *     rateAccretion         &
            &      *     speedLight**2/kilo**2
    else
       rateHeating=+0.0d0
    end if
    return
  end function quasistaticHeatingRate

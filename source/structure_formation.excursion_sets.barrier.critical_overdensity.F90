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
Implements a critical overdensity excursion set barrier class.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <excursionSetBarrier name="excursionSetBarrierCriticalOverdensity">
   <description>
    A excursion set barrier class that adopts a barrier equal to the critical linear theory overdensity for halo collapse.
   </description>
  </excursionSetBarrier>
  !!]
  type, extends(excursionSetBarrierClass) :: excursionSetBarrierCriticalOverdensity
     !!{
     A critical overdensity excursion set barrier class.
     !!}
     private
     class(criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class(cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
   contains
     final     ::                    criticalOverdensityDestructor
     procedure :: barrier         => criticalOverdensityBarrier
     procedure :: barrierGradient => criticalOverdensityBarrierGradient
  end type excursionSetBarrierCriticalOverdensity

  interface excursionSetBarrierCriticalOverdensity
     !!{
     Constructors for the critical overdensity excursion set barrier class.
     !!}
     module procedure criticalOverdensityConstructorParameters
     module procedure criticalOverdensityConstructorInternal
  end interface excursionSetBarrierCriticalOverdensity

contains

  function criticalOverdensityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the critical overdensity excursion set class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (excursionSetBarrierCriticalOverdensity)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(criticalOverdensityClass              ), pointer       :: criticalOverdensity_
    class(cosmologicalMassVarianceClass         ), pointer       :: cosmologicalMassVariance_

    ! Check and read parameters.
    !![
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=excursionSetBarrierCriticalOverdensity(criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function criticalOverdensityConstructorParameters

  function criticalOverdensityConstructorInternal(criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the critical overdensity excursion set class.
    !!}
    implicit none
    type (excursionSetBarrierCriticalOverdensity)                        :: self
    class(criticalOverdensityClass              ), intent(in   ), target :: criticalOverdensity_
    class(cosmologicalMassVarianceClass         ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="*criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    return
  end function criticalOverdensityConstructorInternal

  subroutine criticalOverdensityDestructor(self)
    !!{
    Destructor for the critical overdensity excursion set barrier class.
    !!}
    implicit none
    type(excursionSetBarrierCriticalOverdensity), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine criticalOverdensityDestructor

  double precision function criticalOverdensityBarrier(self,variance,time,node,rateCompute)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierCriticalOverdensity), intent(inout) :: self
    double precision                                        , intent(in   ) :: variance   , time
    type            (treeNode                              ), intent(inout) :: node
    logical                                                 , intent(in   ) :: rateCompute
    double precision                                                        :: mass
    !$GLC attributes unused :: rateCompute

    if (variance <= 0.0d0) then
       ! Return the critical overdensity at this time for infinite mass.
       criticalOverdensityBarrier=self%criticalOverdensity_     %value(time=time,mass       =huge(0.0d0   ),node=node)
    else
       ! Get the mass corresponding to this variance.
       mass                      =self%cosmologicalMassVariance_%mass(time=time,rootVariance=sqrt(variance)          )
       ! Return the critical overdensity at this time at the computed mass scale.
       criticalOverdensityBarrier=self%criticalOverdensity_     %value(time=time,mass       =     mass     ,node=node)
    end if
   return
  end function criticalOverdensityBarrier

  double precision function criticalOverdensityBarrierGradient(self,variance,time,node,rateCompute)
    !!{
    Return the gradient with respect to variance of the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierCriticalOverdensity), intent(inout) :: self
    double precision                                        , intent(in   ) :: variance   , time
    type            (treeNode                              ), intent(inout) :: node
    logical                                                 , intent(in   ) :: rateCompute
    double precision                                                        :: alpha      , mass
    !$GLC attributes unused :: rateCompute

    if (variance <= 0.0d0 .or. .not.self%criticalOverdensity_%isMassDependent()) then
       ! Return zero critical overdensity gradient at this time for infinite mass.
       criticalOverdensityBarrierGradient=0.0d0
    else
       ! Get the halo mass corresponding to this variance.
       mass =self%cosmologicalMassVariance_%mass                           (time=time,rootVariance=sqrt(variance))
       ! Get the logarithmic slope of σ(M).
       alpha=self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient(time=time,mass        =     mass     )
       ! Return the critical overdensity at this time at the computed mass scale.
       criticalOverdensityBarrierGradient=+0.5d0                                                                 &
            &                             *mass                                                                  &
            &                             /variance                                                              &
            &                             /alpha                                                                 &
            &                             *self%criticalOverdensity_%gradientMass(time=time,mass=mass,node=node)
    end if
    return
  end function criticalOverdensityBarrierGradient

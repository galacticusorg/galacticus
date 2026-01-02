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
  Implementation of a rate-limiting stellar feedback model.
  !!}

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsRateLimit">
   <description>A rate-limiting stellar feedback model.</description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsRateLimit
     !!{
     Implementation of a rate-limiting stellar feedback model.
     !!}
     private
     class           (stellarFeedbackOutflowsClass), pointer :: stellarFeedbackOutflows_          => null()
     double precision                                        :: timescaleOutflowFractionalMinimum
   contains
     final     ::                rateLimitDestructor
     procedure :: outflowRate => rateLimitOutflowRate
  end type stellarFeedbackOutflowsRateLimit

  interface stellarFeedbackOutflowsRateLimit
     !!{
     Constructors for the rate-limiting stellar feedback class.
     !!}
     module procedure rateLimitConstructorParameters
     module procedure rateLimitConstructorInternal
  end interface stellarFeedbackOutflowsRateLimit

contains

  function rateLimitConstructorParameters(parameters) result(self)
    !!{
    Constructor for the rate-limiting stellar feedback class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarFeedbackOutflowsRateLimit)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (stellarFeedbackOutflowsClass    ), pointer       :: stellarFeedbackOutflows_
    double precision                                                  :: timescaleOutflowFractionalMinimum

    !![
    <inputParameter>
      <name>timescaleOutflowFractionalMinimum</name>
      <source>parameters</source>
      <description>The minimum timescale (in units of the component dynamical time) for outflows due to stellar feedback.</description>
    </inputParameter>
    <objectBuilder class="stellarFeedbackOutflows" name="stellarFeedbackOutflows_" source="parameters"/>
    !!]
    self=stellarFeedbackOutflowsRateLimit(timescaleOutflowFractionalMinimum,stellarFeedbackOutflows_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stellarFeedbackOutflows_"/>
    !!]
    return
  end function rateLimitConstructorParameters

  function rateLimitConstructorInternal(timescaleOutflowFractionalMinimum,stellarFeedbackOutflows_) result(self)
    !!{
    Internal constructor for the rateLimit stellar feedback class.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsRateLimit)                        :: self
    class           (stellarFeedbackOutflowsClass    ), intent(in   ), target :: stellarFeedbackOutflows_
    double precision                                  , intent(in   )         :: timescaleOutflowFractionalMinimum
    !![
    <constructorAssign variables="timescaleOutflowFractionalMinimum, *stellarFeedbackOutflows_"/>
    !!]
    
    return
  end function rateLimitConstructorInternal

  subroutine rateLimitDestructor(self)
    !!{
    Internal constructor for the rateLimit stellar feedback class.
    !!}
    implicit none
    type(stellarFeedbackOutflowsRateLimit), intent(inout) :: self
    
    !![
    <objectDestructor name="self%stellarFeedbackOutflows_"/>
    !!]
    return
  end subroutine rateLimitDestructor

  subroutine rateLimitOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Limits the outflow rate from another stellar feedback class such that the outflow timescale never falls below a given
    fraction of the dynamical time.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk, nodeComponentSpheroid
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    implicit none
    class           (stellarFeedbackOutflowsRateLimit), intent(inout) :: self
    class           (nodeComponent                   ), intent(inout) :: component
    double precision                                  , intent(in   ) :: rateEnergyInput    , rateStarFormation
    double precision                                  , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive
    double precision                                                     radius             , velocity            , &
         &                                                               timescaleDynamical , massGas             , &
         &                                                               rateOutflowTotal   , rateOutflowMaximum

    call self%stellarFeedbackOutflows_%outflowRate(component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    ! Get properties.
    select type (component)
    class is (nodeComponentDisk    )
       radius  =component%radius  ()
       velocity=component%velocity()
       massGas =component%massGas ()
    class is (nodeComponentSpheroid)
       radius  =component%radius  ()
       velocity=component%velocity()
       massGas =component%massGas ()
    class default
       radius  =0.0d0
       velocity=0.0d0
       massGas =0.0d0
       call Error_Report('unsupported component'//{introspection:location})
    end select
    ! Compute the total outflow rate.
    rateOutflowTotal=+rateOutflowEjective  &
         &           +rateOutflowExpulsive
    ! Compute dynamical timescale.
    if (velocity <= 0.0d0 .or. radius <= 0.0d0) then
       ! Velocity and/or radius is zero, so dynamical timescale is undefined. This is acceptable only if the outflow rate is zero.
       timescaleDynamical=1.0d0
       if (rateOutflowTotal > 0.0d0) call Error_Report('outflow in unphysical component'//{introspection:location})
    else
       timescaleDynamical=+MpcPerKmPerSToGyr &
            &             *radius            &
            &             /velocity
    end if
    ! Compute the maximum outflow rate and current total outflow rate.
    rateOutflowMaximum=max(                                         & 
         &                 +     massGas                            &
         &                 /     timescaleDynamical                 &
         &                 /self%timescaleOutflowFractionalMinimum, &
         &                 +0.0d0                                   &
         &                )
    ! Limit the rate.
    if (rateOutflowTotal > rateOutflowMaximum) then
       rateOutflowEjective =+rateOutflowEjective *rateOutflowMaximum/rateOutflowTotal
       rateOutflowExpulsive=+rateOutflowExpulsive*rateOutflowMaximum/rateOutflowTotal
    end if
    return
  end subroutine rateLimitOutflowRate

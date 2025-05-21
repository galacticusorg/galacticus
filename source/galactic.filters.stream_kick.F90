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

!+    Contributions to this file made by: Paul Menker
  
!!{
Implements a filter for subhalos that could impact a stream during the timestep.
!!}

  use :: Output_Times, only : outputTimesClass

  !![
  <galacticFilter name="galacticFilterStreamKick">
   <description>
     A filter for the velocity kick imparted by subhalos to a stellar stream. We compute an upper limit on the velocity kick a
     subhalo can impart on the stream, for any orientation along a sphere of radius {\normalfont \ttfamily [radiusOrbitalStream]},
     and filter out any subhalos that do not create a total velocity kick greater than {\normalfont \ttfamily
     [cutoffVelocityKick]}.
  </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterStreamKick
     !!{
     A filter for subhalos that could impact a stream during the timestep.
     !!}
     private
     class           (outputTimesClass), pointer :: outputTimes_        => null()
     double precision                            :: radiusOrbitalStream          , cutoffVelocityKick, &
          &                                         speedOrbitalStream
   contains
     final     ::           streamKickDestructor
     procedure :: passes => streamKickPasses
  end type galacticFilterStreamKick

  interface galacticFilterStreamKick
     !!{
     Constructors for the \refClass{galacticFilterStreamKick} galactic filter class.
     !!}
     module procedure streamKickConstructorParameters
     module procedure streamKickConstructorInternal
  end interface galacticFilterStreamKick

contains

  function streamKickConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterStreamKick} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterStreamKick)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (outputTimesClass        ), pointer       :: outputTimes_
    double precision                                          :: radiusOrbitalStream, cutoffVelocityKick, &
         &                                                       speedOrbitalStream
    
    !![
    <inputParameter>
      <name>radiusOrbitalStream</name>
      <source>parameters</source>
      <description>The orbital radius of the stream (which is assumed to be on a circular orbit).</description>
    </inputParameter>
    <inputParameter>
      <name>speedOrbitalStream</name>
      <source>parameters</source>
      <description>The orbital speed of the stream (which is assumed to be on a circular orbit).</description>
    </inputParameter>
    <inputParameter>
      <name>cutoffVelocityKick</name>
      <source>parameters</source>
      <description>The minimum velocity kick to pass.</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=galacticFilterStreamKick(radiusOrbitalStream,speedOrbitalStream,cutoffVelocityKick,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function streamKickConstructorParameters

  function streamKickConstructorInternal(radiusOrbitalStream,speedOrbitalStream,cutoffVelocityKick,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterStreamKick} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterStreamKick)                        :: self
    double precision                          , intent(in   )         :: radiusOrbitalStream, speedOrbitalStream, &
         &                                                               cutoffVelocityKick
    class           (outputTimesClass        ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="radiusOrbitalStream, speedOrbitalStream, cutoffVelocityKick, *outputTimes_"/>
    !!]
    
    return
  end function streamKickConstructorInternal

  subroutine streamKickDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterStreamKick} galactic filter class.
    !!}
    implicit none
    type(galacticFilterStreamKick), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine streamKickDestructor

  logical function streamKickPasses(self,node) result(passes)
    !!{
    Filter based on whether a subhalo can impact a stream in the timestep.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentSatellite
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    class           (galacticFilterStreamKick), intent(inout)          :: self
    type            (treeNode                ), intent(inout), target  :: node
    class           (nodeComponentSatellite  )               , pointer :: satellite
    type            (treeNode                )               , pointer :: nodeWork
    double precision                          , dimension(3)           :: position              , velocity
    double precision                                                   :: speedRelativeMinimum  , velocityKick
    double precision                                                   :: massBound  
    double precision                                                   :: impactParameterMinimum, speed

    if (node%isSatellite()) then
       ! Find the position and velocity of the subhalo relative to its final host.
       nodeWork  => node
       satellite => node     %satellite()
       massBound =  satellite%boundMass()
       position  =  0.0d0
       velocity  =  0.0d0
       do while (nodeWork%isSatellite())
          satellite =>  nodeWork %satellite()
          position  =  +          position    &
               &       +satellite%position ()
          velocity  =  +          velocity    &
               &       +satellite%velocity ()
          nodeWork  =>  nodeWork %parent
       end do
       ! Find the velocity kick, Δv.
       speed                 =    Vector_Magnitude(         velocity                                        )
       impactParameterMinimum=    Vector_Magnitude(position-velocity*Dot_Product(velocity,position)/speed**2)-self%radiusOrbitalStream
       speedRelativeMinimum  =abs(                          speed                                            -self%speedOrbitalStream )
       velocityKick          =+2.0d0                               &
            &                 *    gravitationalConstant_internal  &
            &                 *    massBound                       &
            &                 /abs(speedRelativeMinimum          ) &
            &                 /    impactParameterMinimum
       ! Determine if the node passes.
       passes= velocityKick           >  self%cutoffVelocityKick &
            & .or.                                               &
            &  impactParameterMinimum <= 0.0d0
    else
       ! Non-subhalos are always passed.
       passes=.true.
    end if
    return
  end function streamKickPasses

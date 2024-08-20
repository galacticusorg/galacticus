!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
Implements a filter for subhalos that could impact a stream during the timestep.
!!}

  use :: Output_Times, only : outputTimesClass

  !![
  <galacticFilter name="galacticFilterStreamKick">
   <description>
   A filter for subhalos that could impact a stream during the timestep. We compute an upper limit on the velocity kick a subhalo can impart on the stream, for any orientation along a sphere at 13kpc, and filter
   out any subhalos that do not create a total velocity kick > 0.1km/s.
  </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterStreamKick
     !!{
     A filter for subhalos that could impact a stream during the timestep.
     !!}
     private
     class           (outputTimesClass), pointer :: outputTimes_        => null()
     double precision                            :: radiusOrbitalStream
     double precision                            :: cutoffVelocityKick
   contains
     final     ::           streamKickDestructor
     procedure :: passes => streamKickPasses
  end type galacticFilterStreamKick

  interface galacticFilterStreamKick
     !!{
     Constructors for the ``streamKick'' galactic filter class.
     !!}
     module procedure streamKickConstructorParameters
     module procedure streamKickConstructorInternal
  end interface galacticFilterStreamKick

contains

  function streamKickConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``streamKick'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterStreamKick)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (outputTimesClass          ), pointer       :: outputTimes_
    double precision                                            :: radiusOrbitalStream
    double precision                                            :: cutoffVelocityKick
    
    !![
    <inputParameter>
      <name>radiusOrbitalStream</name>
      <source>parameters</source>
      <description>The orbital radius of the stream (which is assumed to be on a circular orbit).</description>
      </inputParameter>
      <inputParameter>
      <name>cutoffVelocityKick</name>
      <source>parameters</source>
      <description>The minimum velocity kick we want to keep</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=galacticFilterStreamKick(radiusOrbitalStream, cutoffVelocityKick, outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function streamKickConstructorParameters

  function streamKickConstructorInternal(radiusOrbitalStream, cutoffVelocityKick, outputTimes_) result(self)
    !!{
    Internal constructor for the ``streamImpact'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterStreamKick)                        :: self
    double precision                            , intent(in   )         :: radiusOrbitalStream
    double precision                            , intent(in   )         :: cutoffVelocityKick
    class           (outputTimesClass          ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="radiusOrbitalStream, cutoffVelocityKick, *outputTimes_"/>
    !!]
    
    return
  end function streamKickConstructorInternal

  subroutine streamKickDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily streamImpact} galactic filter class.
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
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite
    use :: Numerical_Constants_Astronomical, only : Mpc_per_km_per_s_To_Gyr, gravitationalConstantGalacticus
    use :: Vectors, only : Vector_Magnitude
    implicit none
    class           (galacticFilterStreamKick), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), target  :: node
    class           (nodeComponentBasic        )               , pointer :: basic
    class           (nodeComponentSatellite    )               , pointer :: satellite
    type            (treeNode                  )               , pointer :: nodeWork
    double precision                            , dimension(3)           :: position         , velocity
    double precision                                                     :: w                , deltaV
    double precision                                                     :: massBound  
    double precision                                                     :: b, &
         &                                                                  speed

    if (node%isSatellite()) then
       ! Find the position and velocity of the subhalo relative to its final host.
       nodeWork => node
       satellite => node%satellite()
       position = 0.0d0
       velocity = 0.0d0
       massBound= satellite%boundMass()
       do while (nodeWork%isSatellite())
          satellite => nodeWork %satellite()
          position = position+satellite%position ()
          velocity = velocity+satellite%velocity ()
          nodeWork => nodeWork %parent
       end do


       ! Find deltav
       speed            =Vector_Magnitude(velocity)
       b=ABS(Vector_Magnitude(position - velocity * Dot_Product(velocity,position)/speed**2) - self%radiusOrbitalStream)
       w= ABS(speed - 220)
       deltaV= 2*gravitationalConstantGalacticus*massBound/(w*b)
       ! Determine if the node passes. Note that the impact times computed above are relative to the current time, so we must include that offset here.
       basic  => node%basic()
       passes =  deltaV > self%cutoffVelocityKick                                                                                           &
            &    .or.                                                                                                                       &
            &     Vector_Magnitude(position - velocity * Dot_Product(velocity,position)/speed**2) < self%radiusOrbitalStream
       ! FOR modification for elliptical streams
       ! speed            =Vector_Magnitude(velocity)
       ! b=ABS(Vector_Magnitude(position - velocity * Dot_Product(velocity,position)/speed**2) - maximum(self%radiusOrbitalStream))
       ! w= ABS(speed - maximum(vstream))
       ! deltaV= 2*gravitationalConstantGalacticus*massBound/(w*b)
       ! passes =  deltaV > self%cutoffVelocityKick                                                                                           &
       !      &    .or.                                                                                                                       &
       !      &     Vector_Magnitude(position - velocity * Dot_Product(velocity,position)/speed**2) < maximum(self%radiusOrbitalStream)
       ! write(0, *) "only subhalos with r>r_stream.  mass, w, b,  deltaV, passes, (r-. 013), node index are:", massBound, w, b, deltaV, passes, (Vector_Magnitude(position - velocity * Dot_Product(velocity,position)/speed**2) - .013), node%index()
    else
       ! Non-subhalos are always passed.
       passes=.true.
    end if
    return
  end function streamKickPasses

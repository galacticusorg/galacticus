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
Implements a filter for subhalos that could impact a stream during the timestep.
!!}

  use :: Output_Times, only : outputTimesClass

  !![
  <galacticFilter name="galacticFilterStreamImpact">
   <description>
     A filter for subhalos that could impact a stream during the timestep. We consider a subhalo, at time $t=0$ (defined as the
     current time), at position $\mathbf{r}$ and moving with velocity $\mathbf{v}$. The time of closest approach to a point
     $\mathbf{r}^\prime$ is given by
     \begin{equation}
      t_\mathrm{impact} = (\mathbf{v}\cdot \mathbf{r}^\prime - \mathbf{v}\cdot \mathbf{r})/v^2.
     \end{equation}     
     We want to keep subhalos which may impact upon a stream of radius $r^\prime$ during the timestep, allowing for arbitrary
     rotations of the subhalo system. Therefore, we must find the minimum and maximum possible values of $t_\mathrm{impact}$ when
     considering all points on the sphere of radius $r^\prime$. These extrema clearly occur when $\mathbf{r}^\prime$ is aligned,
     or anti-aligned with $\mathbf{v}$, i.e.:
     \begin{equation}
      t_\mathrm{impact, min/max} = (\pm v r^\prime - \mathbf{v}\cdot \mathbf{r})/v^2.
     \end{equation}
  </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterStreamImpact
     !!{
     A filter for subhalos that could impact a stream during the timestep.
     !!}
     private
     class           (outputTimesClass), pointer :: outputTimes_        => null()
     double precision                            :: radiusOrbitalStream
   contains
     final     ::           streamImpactDestructor
     procedure :: passes => streamImpactPasses
  end type galacticFilterStreamImpact

  interface galacticFilterStreamImpact
     !!{
     Constructors for the \refClass{galacticFilterStreamImpact} galactic filter class.
     !!}
     module procedure streamImpactConstructorParameters
     module procedure streamImpactConstructorInternal
  end interface galacticFilterStreamImpact

contains

  function streamImpactConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterStreamImpact} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterStreamImpact)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (outputTimesClass          ), pointer       :: outputTimes_
    double precision                                            :: radiusOrbitalStream

    !![
    <inputParameter>
      <name>radiusOrbitalStream</name>
      <source>parameters</source>
      <description>The orbital radius of the stream (which is assumed to be on a circular orbit).</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=galacticFilterStreamImpact(radiusOrbitalStream,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function streamImpactConstructorParameters

  function streamImpactConstructorInternal(radiusOrbitalStream,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterStreamImpact} galactic filter class.
    !!}
    implicit none
    type            (galacticFilterStreamImpact)                        :: self
    double precision                            , intent(in   )         :: radiusOrbitalStream
    class           (outputTimesClass          ), intent(in   ), target :: outputTimes_
    !![
    <constructorAssign variables="radiusOrbitalStream, *outputTimes_"/>
    !!]
    
    return
  end function streamImpactConstructorInternal

  subroutine streamImpactDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterStreamImpact} galactic filter class.
    !!}
    implicit none
    type(galacticFilterStreamImpact), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine streamImpactDestructor

  logical function streamImpactPasses(self,node) result(passes)
    !!{
    Filter based on whether a subhalo can impact a stream in the timestep.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic, nodeComponentSatellite
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    class           (galacticFilterStreamImpact), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), target  :: node
    class           (nodeComponentBasic        )               , pointer :: basic
    class           (nodeComponentSatellite    )               , pointer :: satellite
    type            (treeNode                  )               , pointer :: nodeWork
    double precision                            , dimension(3)           :: position         , velocity
    double precision                                                     :: timeImpactMinimum, timeImpactMaximum, &
         &                                                                  speed            , timeOutputNext

    if (node%isSatellite()) then
       ! Find the position and velocity of the subhalo relative to its final host.
       nodeWork => node
       position =  0.0d0
       velocity =  0.0d0
       do while (nodeWork%isSatellite())
          satellite =>          nodeWork %satellite()
          position  =  position+satellite%position ()
          velocity  =  velocity+satellite%velocity ()
          nodeWork  =>          nodeWork %parent
       end do
       ! Find the minimum and maximum times of possible impact on the stream.
       speed            =Vector_Magnitude(velocity)
       timeImpactMinimum=(-speed*self%radiusOrbitalStream-Dot_Product(velocity,position))/speed**2*MpcPerKmPerSToGyr
       timeImpactMaximum=(+speed*self%radiusOrbitalStream-Dot_Product(velocity,position))/speed**2*MpcPerKmPerSToGyr
       ! Determine if the node passes. Note that the impact times computed above are relative to the current time, so we must
       ! include that offset here.
       basic          => node             %basic   (            )
       timeOutputNext =  self%outputTimes_%timeNext(basic%time())
       passes         =  timeOutputNext    > 0.0d0                       & ! Negative next output time indicates no more outputs,
            &           .and.                                            & ! so no stream impact can occur.
            &            timeImpactMinimum < timeOutputNext-basic%time() &
            &           .and.                                            &
            &            timeImpactMaximum > 0.0d0
    else
       ! Non-subhalos are always passed.
       passes=.true.
    end if
    return
  end function streamImpactPasses

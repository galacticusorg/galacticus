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
  Implementation of a satellite dynamical friction class that truncates the acceleration to zero below a satellite/host mass
  ratio.
  !!}

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionMassRatioThreshold">
   <description>
    A satellite dynamical friction class that truncates the acceleration to zero below a satellite/host mass ratio.
   </description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionClass) :: satelliteDynamicalFrictionMassRatioThreshold
     !!{
     Implementation of a satellite dynamical friction class that truncates the acceleration to zero below a satellite/host mass ratio.
     !!}
     private
     class           (satelliteDynamicalFrictionClass), pointer :: satelliteDynamicalFriction_ => null()
     double precision                                           :: massRatioThreshold
   contains
     final     ::                 massRatioThresholdDestructor
     procedure :: acceleration => massRatioThresholdAcceleration
  end type satelliteDynamicalFrictionMassRatioThreshold

  interface satelliteDynamicalFrictionMassRatioThreshold
     !!{
     Constructors for the massRatioThreshold satellite dynamical friction class.
     !!}
     module procedure massRatioThresholdConstructorParameters
     module procedure massRatioThresholdConstructorInternal
  end interface satelliteDynamicalFrictionMassRatioThreshold

contains

  function massRatioThresholdConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteDynamicalFrictionMassRatioThreshold} satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteDynamicalFrictionMassRatioThreshold)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (satelliteDynamicalFrictionClass             ), pointer       :: satelliteDynamicalFriction_
    double precision                                                              :: massRatioThreshold

    !![
    <inputParameter>
      <name>massRatioThreshold</name>
      <description>The satellite-to-host mass ratio below which dynamical friction acceleration is truncated to zero.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="satelliteDynamicalFriction" name="satelliteDynamicalFriction_" source="parameters"/>
    !!]
    self=satelliteDynamicalFrictionMassRatioThreshold(massRatioThreshold,satelliteDynamicalFriction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteDynamicalFriction_"/>
    !!]
    return
  end function massRatioThresholdConstructorParameters

  function massRatioThresholdConstructorInternal(massRatioThreshold,satelliteDynamicalFriction_) result(self)
    !!{
    Internal constructor for the \refClass{satelliteDynamicalFrictionMassRatioThreshold} satellite dynamical friction class.
    !!}
    implicit none
    type            (satelliteDynamicalFrictionMassRatioThreshold)                        :: self
    class           (satelliteDynamicalFrictionClass             ), intent(in   ), target :: satelliteDynamicalFriction_
    double precision                                              , intent(in   )         :: massRatioThreshold
    !![
    <constructorAssign variables="massRatioThreshold, *satelliteDynamicalFriction_"/>
    !!]

    return
  end function massRatioThresholdConstructorInternal

  subroutine massRatioThresholdDestructor(self)
    !!{
    Destructor for the \refClass{satelliteDynamicalFrictionMassRatioThreshold} satellite dynamical friction class.
    !!}
    implicit none
    type(satelliteDynamicalFrictionMassRatioThreshold), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteDynamicalFriction_"/>
    !!]
    return
  end subroutine massRatioThresholdDestructor

  function massRatioThresholdAcceleration(self,node) result(acceleration)
    !!{
    Return an acceleration for satellites due to dynamical, truncating to zero below a given satellite-to-host mass ratio.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite
    implicit none
    double precision                                              , dimension(3)          :: acceleration
    class           (satelliteDynamicalFrictionMassRatioThreshold), intent(inout), target :: self
    type            (treeNode                                    ), intent(inout)         :: node
    class           (nodeComponentBasic                          ), pointer               :: basicHost
    class           (nodeComponentSatellite                      ), pointer               :: satellite
    type            (treeNode                                    ), pointer               :: nodeHost

    nodeHost                      =>  node    %mergesWith()
    basicHost                     =>  nodeHost%basic     ()
    satellite                     =>  node    %satellite ()
    if (satellite%boundMass() < self%massRatioThreshold*basicHost%mass()) then
       acceleration=0.0d0
    else
       acceleration=self%satelliteDynamicalFriction_%acceleration(node)
    end if
    return
  end function massRatioThresholdAcceleration

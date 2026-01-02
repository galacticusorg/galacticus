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
  An abstract implementation of the orphan satellite distribution which assumes an isotropic distribution with randomly
  assigned positions.
  !!}

  use :: Statistics_Distributions, only : distributionFunction1DNormal

  !![
  <satelliteOrphanDistribution name="satelliteOrphanDistributionRandomIsotropic" abstract="yes">
   <description>An abstract orphan satellite distribution which assumes an isotropic, random distribution of positions, and velocities drawn from an isotropic normal distribution. The radial distribution and velocity dispersion must be specified by the child class.</description>
  </satelliteOrphanDistribution>
  !!]
  type, abstract, extends(satelliteOrphanDistributionClass) :: satelliteOrphanDistributionRandomIsotropic
     !!{
     An abstract orphan satellite distribution which assumes an isotropic, random distribution.
     !!}
     private
     type(distributionFunction1DNormal) :: normalDistribution
   contains
     !![
     <methods>
       <method description="Return the radius enclosing the given fraction of the orphan satellite population." method="inverseCumulativeMassFunctionRadial" />
       <method description="Return the 1-D velocity dispersion of the orphan satellite population." method="velocityDispersion" />
       <method description="Initialize the class." method="initialize" />
     </methods>
     !!]
     procedure                                              :: position                            => randomIsotropicPosition
     procedure                                              :: velocity                            => randomIsotropicVelocity
     procedure                                              :: initialize                          => randomIsotropicInitialize
     procedure(randomIsotropicInverseCMFRadial  ), deferred :: inverseCumulativeMassFunctionRadial
     procedure(randomIsotropicVelocityDispersion), deferred :: velocityDispersion
  end type satelliteOrphanDistributionRandomIsotropic

  abstract interface
     !!{
     Abstract interface for the inverse cumulative mass function for the radial coordinate in the {\normalfont \ttfamily
     randomIsotropic} orphan satellite distribution class.
     !!}
     double precision function randomIsotropicInverseCMFRadial(self,node,fraction)
       import satelliteOrphanDistributionRandomIsotropic, treeNode
       class           (satelliteOrphanDistributionRandomIsotropic), intent(inout) :: self
       type            (treeNode                                  ), intent(inout) :: node
       double precision                                            , intent(in   ) :: fraction
     end function randomIsotropicInverseCMFRadial
  end interface

  abstract interface
     !!{
     Abstract interface for the velocity dispersion in the {\normalfont \ttfamily randomIsotropic} orphan satellite
     distribution class.
     !!}
     double precision function randomIsotropicVelocityDispersion(self,node)
       import satelliteOrphanDistributionRandomIsotropic, treeNode
       class(satelliteOrphanDistributionRandomIsotropic), intent(inout) :: self
       type (treeNode                                  ), intent(inout) :: node
     end function randomIsotropicVelocityDispersion
  end interface

contains

  subroutine randomIsotropicInitialize(self)
    !!{
    Perform initialization for the {\normalfont \ttfamily randomIsotropic} orphan satellite distribution class.
    !!}
    implicit none
    class           (satelliteOrphanDistributionRandomIsotropic), intent(inout) :: self
    double precision                                            , parameter     :: limitDistribution=5.0d0

    self%normalDistribution=distributionFunction1DNormal(0.0d0,1.0d0,-limitDistribution,+limitDistribution)
    return
  end subroutine randomIsotropicInitialize
  
  function randomIsotropicPosition(self,node)
    !!{
    Return the position of an orphan satellite in a random isotropic distribution.
    !!}
    use :: Coordinates             , only : assignment(=)        , coordinateCartesian, coordinateSpherical
    use :: Galacticus_Nodes        , only : nodeComponentPosition, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                                            , dimension(3)  :: randomIsotropicPosition
    class           (satelliteOrphanDistributionRandomIsotropic), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    type            (treeNode                                  ), pointer       :: nodeHost
    class           (nodeComponentPosition                     ), pointer       :: positionHost
    type            (coordinateSpherical                       )                :: positionSpherical
    type            (coordinateCartesian                       )                :: positionCartesian

    ! Select random spherical coordinates.
    positionSpherical=[                                                                                                                    &
         &             self%inverseCumulativeMassFunctionRadial(node,         node%hostTree%randomNumberGenerator_%uniformSample()      ), &
         &             acos                                    (     2.0d0   *node%hostTree%randomNumberGenerator_%uniformSample()-1.0d0), &
         &                                                           2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()         &
         &            ]
    ! Determine the position in Cartesian coordinates and offset to position of host halo.
    nodeHost               =>  node                   %parent
    positionHost           =>  nodeHost               %position()
    positionCartesian      =   positionSpherical
    randomIsotropicPosition=   positionCartesian
    randomIsotropicPosition=  +randomIsotropicPosition            &
         &                    +positionHost           %position()
    return
  end function randomIsotropicPosition

  function randomIsotropicVelocity(self,node)
    !!{
    Return the velocity of an orphan satellite in a random isotropic distribution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition, treeNode
    implicit none
    double precision                                            , dimension(3)  :: randomIsotropicVelocity
    class           (satelliteOrphanDistributionRandomIsotropic), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    type            (treeNode                                  ), pointer       :: nodeHost
    class           (nodeComponentPosition                     ), pointer       :: positionHost
    double precision                                                            :: velocityDispersion
    integer                                                                     :: i

    ! Get the 1-D velocity dispersion for the orphan.
    velocityDispersion=self%velocityDispersion(node)
    ! Find the velocity of the host.
    nodeHost               =>  node        %parent
    positionHost           =>  nodeHost    %position()
    randomIsotropicVelocity=  +positionHost%velocity()
    ! Add velocity of the orphan.
    do i=1,3
       randomIsotropicVelocity(i)=randomIsotropicVelocity(i)+velocityDispersion*self%normalDistribution%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
    end do
    return
  end function randomIsotropicVelocity

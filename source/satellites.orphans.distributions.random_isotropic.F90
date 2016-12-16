!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An abstract implementation of the orphan satellite distribution which assumes an isotropic distribution with randomly
  !% assigned positions.

  !# <satelliteOrphanDistribution name="satelliteOrphanDistributionRandomIsotropic" abstract="yes">
  !#  <description>An abstract orphan satellite distribution which assumes an isotropic, random distribution. The radial distribution must be specified by the child class.</description>
  !# </satelliteOrphanDistribution>
  type, abstract, extends(satelliteOrphanDistributionClass) :: satelliteOrphanDistributionRandomIsotropic
     !% An abstract orphan satellite distribution which assumes an isotropic, random distribution.
     private
   contains
     !@ <objectMethods>
     !@   <object>satelliteOrphanDistributionRandomIsotropic</object>
     !@   <objectMethod>
     !@     <method>inverseCumulativeMassFunctionRadial</method>
     !@     <description>Return the radius enclosing the given fraction of the orphan satellite population.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>textcolor{red}{\textless type(treeNode)\textgreater} node\arginout, \doublezero\ fraction\argin</arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                            :: position                            => randomIsotropicPosition
     procedure(randomIsotropicInverseCMFRadial), deferred :: inverseCumulativeMassFunctionRadial
  end type satelliteOrphanDistributionRandomIsotropic

  abstract interface
     !% Abstract interface for the inverse cumulative mass function for the radial coordinate in the {\normalfont \ttfamily
     !% randomIsotropic} orphan satellite distribution class.
     double precision function randomIsotropicInverseCMFRadial(self,node,fraction)
       import satelliteOrphanDistributionRandomIsotropic, treeNode
       class           (satelliteOrphanDistributionRandomIsotropic), intent(inout) :: self
       type            (treeNode                                  ), intent(inout) :: node
       double precision                                            , intent(in   ) :: fraction
     end function randomIsotropicInverseCMFRadial
  end interface
  
contains

  function randomIsotropicPosition(self,node)
    !% Return the position of an orphan satellite in a random isotropic distribution.
    use Pseudo_Random
    use Numerical_Constants_Math
    use Coordinates
    implicit none
    double precision                                            , dimension(3)  :: randomIsotropicPosition
    class           (satelliteOrphanDistributionRandomIsotropic), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    type            (treeNode                                  ), pointer       :: nodeHost
    class           (nodeComponentPosition                     ), pointer       :: positionHost
    type            (coordinateSpherical                       )                :: positionSpherical
    type            (coordinateCartesian                       )                :: positionCartesian
    
    ! Select random spherical coordinates.
    positionSpherical=[                                                                                                            &
         &             self%inverseCumulativeMassFunctionRadial(node,         node%hostTree%randomNumberGenerator%sample()      ), &
         &             acos                                    (     2.0d0   *node%hostTree%randomNumberGenerator%sample()-1.0d0), &
         &                                                           2.0d0*Pi*node%hostTree%randomNumberGenerator%sample()         &
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

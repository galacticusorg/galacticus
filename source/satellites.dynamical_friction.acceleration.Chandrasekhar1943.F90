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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!% Contains a module with a \cite{chandrasekhar_dynamical_1943} implementation of calculations of satellite acceleration due to dynamical friction.

module Dynamical_Friction_Acceleration_Chandrasekhar
  !% Implements \cite{chandrasekhar_dynamical_1943} value of calculations of satellite acceleration due to dynamical friction.
  implicit none
  private
  public :: Satellite_Dynamical_Friction_Chandrasekhar_Initialize

  ! Coulomb logarithm parameter.
  double precision :: satelliteDynamicalFrictionChandrasekharCoulombLogarithm

contains

  !# <satelliteDynamicalFrictionMethod>
  !#  <unitName>Satellite_Dynamical_Friction_Chandrasekhar_Initialize</unitName>
  !# </satelliteDynamicalFrictionMethod>
  subroutine Satellite_Dynamical_Friction_Chandrasekhar_Initialize(satelliteDynamicalFrictionMethod,Satellite_Dynamical_Friction_Acceleration)
    !% Determine if this method is to be used and set pointer appropriately.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type     (varying_string                                         ),          intent(in   ) :: satelliteDynamicalFrictionMethod
    procedure(Satellite_Dynamical_Friction_Acceleration_Chandrasekhar), pointer, intent(inout) :: Satellite_Dynamical_Friction_Acceleration

    if (satelliteDynamicalFrictionMethod == 'Chandrasekhar1943') then
       Satellite_Dynamical_Friction_Acceleration => Satellite_Dynamical_Friction_Acceleration_Chandrasekhar
       !@ <inputParameter>
       !@   <name>satelliteDynamicalFrictionChandrasekharCoulombLogarithm</name>
       !@   <defaultValue>2</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The Coulomb logarithm, $\ln \Lambda$, appearing in the \cite{chandrasekhar_dynamical_1943} formulation of the acceleration due to dynamical friction.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group></group>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteDynamicalFrictionChandrasekharCoulombLogarithm',satelliteDynamicalFrictionChandrasekharCoulombLogarithm,defaultValue=2.0d0)
   end if
    return
  end subroutine Satellite_Dynamical_Friction_Chandrasekhar_Initialize

  function Satellite_Dynamical_Friction_Acceleration_Chandrasekhar(thisNode)
    !% Return the \cite{chandrasekhar_dynamical_1943} acceleration for satellites due to dynamical friction.
    use Galacticus_Nodes
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    use Error_Functions
    use Galactic_Structure_Densities
    use Galactic_Structure_Options
    use Vectors
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode                ), pointer     , intent(inout) :: thisNode
    double precision                          , dimension(3)                :: Satellite_Dynamical_Friction_Acceleration_Chandrasekhar
    class           (nodeComponentSatellite  ), pointer                     :: thisSatellite
    type            (treeNode                ), pointer                     :: hostNode
    class           (darkMatterHaloScaleClass), pointer                     :: darkMatterHaloScale_
    double precision                          , dimension(3)                :: position,velocity
    double precision                                                        :: satelliteMass,parentDensity,velocityMagnitude,velocityDispersion,Xv

    darkMatterHaloScale_ => darkMatterHaloScale()
    hostNode             => thisNode     %mergesWith()
    thisSatellite        => thisNode     %satellite ()
    satelliteMass        =  thisSatellite%boundMass ()
    position             =  thisSatellite%position  ()
    velocity             =  thisSatellite%velocity  ()
    velocityMagnitude    =  Vector_Magnitude                   (velocity)
    velocityDispersion   =  darkMatterHaloScale_%virialVelocity(hostNode)
    Xv                   =  velocityMagnitude/velocityDispersion/sqrt(2.0d0)
    parentDensity        =  Galactic_Structure_Density(hostNode,position,coordinateSystemCartesian)
    Satellite_Dynamical_Friction_Acceleration_Chandrasekhar=        &
         & -4.0d0                                                   &
         & *Pi                                                      &
         & *satelliteDynamicalFrictionChandrasekharCoulombLogarithm &
         & *gravitationalConstantGalacticus**2                      &
         & *satelliteMass                                           &
         & *parentDensity                                           &
         & /velocityMagnitude**3                                    &
         & *(                                                       &
         &   Error_Function(Xv)                                     &
         &   -2.0d0                                                 &
         &   *Xv                                                    &
         &   *exp(-Xv**2)                                           &
         &   /sqrt(Pi)                                              &
         &  )                                                       &
         & *velocity                                                &
         & *kilo*gigaYear/megaParsec
    return
  end function Satellite_Dynamical_Friction_Acceleration_Chandrasekhar

end module Dynamical_Friction_Acceleration_Chandrasekhar

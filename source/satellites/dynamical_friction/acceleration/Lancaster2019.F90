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

  !!{RST
  Implementation of a satellite dynamical friction class which uses the model of :cite:t:`lancaster_dynamical_2020` for the Coulomb logarithm.
  !!}  

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionLancaster2020" docformat="rst">
   <description>
   A satellite dynamical friction class which computes the Coulomb logarithm following the prescription of :cite:t:`lancaster_dynamical_2020`.
   The Coulomb logarithm is evaluated using :math:`\Lambda = r_\mathrm{orbital}/l_{90}`, where :math:`l_{90}=\mathrm{G}M_\mathrm{sat}/v_\mathrm{orbital}^2`
   is the impact parameter for an encounter causing a 90-degree deflection.
   </description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionChandrasekhar1943) :: satelliteDynamicalFrictionLancaster2020
     !!{RST
     Implementation of a satellite dynamical friction class which uses the model of :cite:t:`lancaster_dynamical_2020` for the Coulomb logarithm.
     !!}
     private
   contains
     procedure :: coulombLogarithm => lancaster2020CoulombLogarithm
  end type satelliteDynamicalFrictionLancaster2020

  interface satelliteDynamicalFrictionLancaster2020
     !!{RST
     Constructors for the lancaster2020 satellite dynamical friction class.
     !!}
     module procedure lancaster2020ConstructorParameters
     module procedure lancaster2020ConstructorInternal
  end interface satelliteDynamicalFrictionLancaster2020

contains

  function lancaster2020ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteDynamicalFrictionLancaster2020` satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(satelliteDynamicalFrictionlancaster2020)                :: self
    type(inputParameters                        ), intent(inout) :: parameters
    
    self=satelliteDynamicalFrictionLancaster2020()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function lancaster2020ConstructorParameters

  function lancaster2020ConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteDynamicalFrictionLancaster2020` satellite dynamical friction class.
    !!}
    implicit none
    type(satelliteDynamicalFrictionLancaster2020) :: self

    ! Set the fixed (and unused by this class) `logarithmCoulomb` of the parent class to an
    ! unphysical value.
    self%logarithmCoulomb=-huge(0.0d0)
    return
  end function lancaster2020ConstructorInternal
  
  double precision function lancaster2020CoulombLogarithm(self,node) result(coulombLogarithm)
    !!{RST
    Evaluate the Coulomb logarithm for the :cite:t:`lancaster_dynamical_2020` dynamical friction model.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentSatellite, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    class           (satelliteDynamicalFrictionLancaster2020), intent   (inout) :: self
    type            (treeNode                               ), intent   (inout) :: node
    class           (nodeComponentSatellite                 ), pointer          :: satellite
    double precision                                         , dimension(3    ) :: position     , velocity    
    double precision                                                            :: speedOrbital , radiusOrbital          , &
         &                                                                         massSatellite, impactParameter90Degree

    satellite               =>  node     %satellite (        )
    massSatellite           =   satellite%boundMass (        )
    position                =   satellite%position  (        )
    velocity                =   satellite%velocity  (        )
    radiusOrbital           =   Vector_Magnitude    (position)
    speedOrbital            =   Vector_Magnitude    (velocity)
    impactParameter90Degree =  +gravitationalConstant_internal    &
         &                     *massSatellite                     &
         &                     /speedOrbital                  **2
    coulombLogarithm        =  +radiusOrbital                     &
         &                    /impactParameter90Degree
    return
  end function lancaster2020CoulombLogarithm

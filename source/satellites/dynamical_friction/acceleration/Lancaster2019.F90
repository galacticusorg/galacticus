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
  Implementation of a satellite dynamical friction class which uses the model of \cite{petts_semi-analytic_2015} for the
  Coulomb logarithm.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionLancaster2019">
   <description>
   </description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionChandrasekhar1943) :: satelliteDynamicalFrictionLancaster2019
     !!{
     Implementation of a satellite dynamical friction class which uses the model of \cite{petts_semi-analytic_2015} for the
     Coulomb logarithm.
     !!}
     private
   contains
     final     ::                     lancaster2019Destructor
     procedure :: coulombLogarithm => lancaster2019CoulombLogarithm
  end type satelliteDynamicalFrictionLancaster2019

  interface satelliteDynamicalFrictionLancaster2019
     !!{
     Constructors for the lancaster2019 satellite dynamical friction class.
     !!}
     module procedure lancaster2019ConstructorParameters
     module procedure lancaster2019ConstructorInternal
  end interface satelliteDynamicalFrictionLancaster2019

contains

  function lancaster2019ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteDynamicalFrictionLancaster2019} satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (satelliteDynamicalFrictionlancaster2019)                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    
    self=satelliteDynamicalFrictionLancaster2019()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function lancaster2019ConstructorParameters

  function lancaster2019ConstructorInternal( ) result(self)
    !!{
    Internal constructor for the \refClass{satelliteDynamicalFrictionlancaster2019} satellite dynamical friction class.
    !!}
    implicit none
    type   (satelliteDynamicalFrictionlancaster2019)                        :: self

    return
  end function lancaster2019ConstructorInternal

  subroutine lancaster2019Destructor(self)
    !!{
    Destructor for the \refClass{satelliteDynamicalFrictionLancaster2019} satellite dynamical friction class.
    !!}
    implicit none
    type(satelliteDynamicalFrictionLancaster2019), intent(inout) :: self

    return
  end subroutine lancaster2019Destructor
  
  double precision function lancaster2019CoulombLogarithm(self,node) result(coulombLogarithm)
    !!{
    Evaluate the Coulomb logarithm for the \cite{petts_semi-analytic_2015} dynamical friction model.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentSatellite, treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    class           (satelliteDynamicalFrictionLancaster2019), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentSatellite             ), pointer       :: satellite
    double precision                                     , dimension(3)  :: position               , velocity    
    double precision                                                     :: speedOrbital           , radiusOrbital          , &
         &                                                                  massSatellite          , l90

    satellite               =>  node     %satellite (        )
    massSatellite           =   satellite%boundMass (        )
    position                =   satellite%position  (        )
    velocity                =   satellite%velocity  (        )
    radiusOrbital           =   Vector_Magnitude    (position )
    speedOrbital            =   Vector_Magnitude    (velocity )

    l90   = gravitationalConstant_internal*massSatellite/(speedOrbital**2.0d0)
    coulombLogarithm=radiusOrbital/l90
    return
  end function lancaster2019CoulombLogarithm

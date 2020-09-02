!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !+ Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !% Implementation of a satellite dynamical friction class which uses the model of \cite{chandrasekhar_dynamical_1943}.

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !# <satelliteDynamicalFriction name="satelliteDynamicalFrictionChandrasekhar1943">
  !#  <description>A satellite dynamical friction class which uses the model of \cite{chandrasekhar_dynamical_1943}.</description>
  !# </satelliteDynamicalFriction>
  type, extends(satelliteDynamicalFrictionClass) :: satelliteDynamicalFrictionChandrasekhar1943
     !% Implementation of a satellite dynamical friction class which uses the model of \cite{chandrasekhar_dynamical_1943}.
     private
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     double precision                                     :: logarithmCoulomb
   contains
     final     ::                 chandrasekhar1943Destructor
     procedure :: acceleration => chandrasekhar1943Acceleration
  end type satelliteDynamicalFrictionChandrasekhar1943

  interface satelliteDynamicalFrictionChandrasekhar1943
     !% Constructors for the chandrasekhar1943 satellite dynamical friction class.
     module procedure chandrasekhar1943ConstructorParameters
     module procedure chandrasekhar1943ConstructorInternal
  end interface satelliteDynamicalFrictionChandrasekhar1943

contains

  function chandrasekhar1943ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily chandrasekhar1943} satellite dynamical friction class which builds the object from a parameter set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteDynamicalFrictionChandrasekhar1943)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                   ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                  ), pointer       :: darkMatterProfileDMO_
    double precision                                                             :: logarithmCoulomb

    !# <inputParameter>
    !#   <name>logarithmCoulomb</name>
    !#   <defaultValue>2.0d0</defaultValue>
    !#   <description>The Coulomb logarithm, $\ln \Lambda$, appearing in the \cite{chandrasekhar_dynamical_1943} formulation of the acceleration due to dynamical friction.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !# <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    self=satelliteDynamicalFrictionChandrasekhar1943(logarithmCoulomb,darkMatterHaloScale_,darkMatterProfileDMO_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_" />
    !# <objectDestructor name="darkMatterProfileDMO_"/>
    return
  end function chandrasekhar1943ConstructorParameters

  function chandrasekhar1943ConstructorInternal(logarithmCoulomb,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily chandrasekhar1943} satellite dynamical friction class.
    implicit none
    type            (satelliteDynamicalFrictionChandrasekhar1943)                        :: self
    class           (darkMatterHaloScaleClass                   ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                  ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                             , intent(in   )         :: logarithmCoulomb
    !# <constructorAssign variables="logarithmCoulomb, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>

    return
  end function chandrasekhar1943ConstructorInternal

  subroutine chandrasekhar1943Destructor(self)
    !% Destructor for the simple cooling radius class.
    implicit none
    type(satelliteDynamicalFrictionChandrasekhar1943), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_" />
    !# <objectDestructor name="self%darkMatterProfileDMO_"/>
    return
  end subroutine chandrasekhar1943Destructor

  function chandrasekhar1943Acceleration(self,node)
    !% Return an acceleration for satellites due to dynamical friction using the formulation of \cite{chandrasekhar_dynamical_1943}.
    use :: Error_Functions                           , only : Error_Function
    use :: Galactic_Structure_Densities              , only : Galactic_Structure_Density
    use :: Galactic_Structure_Options                , only : coordinateSystemCartesian
    use :: Galactic_Structure_Chandrasekhar_Integrals, only : Galactic_Structure_Chandrasekhar_Integral
    use :: Galacticus_Nodes                          , only : nodeComponentSatellite                   , treeNode
    use :: Numerical_Constants_Astronomical          , only : gigaYear                                 , megaParsec
    use :: Numerical_Constants_Math                  , only : Pi
    use :: Numerical_Constants_Astronomical          , only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes              , only : kilo
    use :: Vectors                                   , only : Vector_Magnitude
    implicit none
    double precision                                             , dimension(3)  :: chandrasekhar1943Acceleration
    class           (satelliteDynamicalFrictionChandrasekhar1943), intent(inout) :: self
    type            (treeNode                                   ), intent(inout) :: node
    class           (nodeComponentSatellite                     ), pointer       :: satellite
    type            (treeNode                                   ), pointer       :: nodeHost
    double precision                                             , dimension(3)  :: position                            , velocity
    double precision                                                             :: massSatellite                       , densityHost       , &
         &                                                                          velocityMagnitude                   , velocityDispersion, &
         &                                                                          Xv                                  , radius
    double precision                                             , parameter     :: XvMaximum                    =10.0d0

    nodeHost                     =>  node     %mergesWith                                    (               )
    satellite                    =>  node     %satellite                                     (               )
    massSatellite                =   satellite%boundMass                                     (               )
    position                     =   satellite%position                                      (               )
    velocity                     =   satellite%velocity                                      (               )
    radius                       =   Vector_Magnitude          (         position                            )
    velocityDispersion           =   self     %darkMatterProfileDMO_%radialVelocityDispersion(nodeHost,radius)
    velocityMagnitude            =   Vector_Magnitude          (         velocity                            )
    densityHost                  =   Galactic_Structure_Density(nodeHost,position,coordinateSystemCartesian  )
    if (velocityDispersion > 0.0d0) then
       Xv                        =  +velocityMagnitude                     &
            &                       /velocityDispersion                    &
            &                       /sqrt(2.0d0)
    else
       Xv                        =  +huge(0.0d0)
    end if
    if (Xv <= XvMaximum) then
       chandrasekhar1943Acceleration=  -4.0d0                              &
            &                          *Pi                                 &
            &                          *self%logarithmCoulomb              &
            &                          *gravitationalConstantGalacticus**2 &
            &                          *massSatellite                      &
            &                          *densityHost                        &
            &                          /velocityMagnitude**3               &
            &                          *(                                  &
            &                            +Error_Function(Xv)               &
            &                            -2.0d0                            &
            &                            *Xv                               &
            &                            *exp(-Xv**2)                      &
            &                            /sqrt(Pi)                         &
            &                           )                                  &
            &                          *velocity                           &
            &                          *kilo                               &
            &                          *gigaYear                           &
            &                          /megaParsec
    else
       chandrasekhar1943Acceleration=  -4.0d0                              &
            &                          *Pi                                 &
            &                          *self%logarithmCoulomb              &
            &                          *gravitationalConstantGalacticus**2 &
            &                          *massSatellite                      &
            &                          *densityHost                        &
            &                          /velocityMagnitude**3               &
            &                          *velocity                           &
            &                          *kilo                               &
            &                          *gigaYear                           &
            &                          /megaParsec
    end if
    return
  end function chandrasekhar1943Acceleration


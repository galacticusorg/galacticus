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

  !+ Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !!{
  Implementation of a satellite dynamical friction class which uses the model of \cite{chandrasekhar_dynamical_1943}.
  !!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionChandrasekhar1943">
   <description>
    A satellite dynamical friction class which uses the \cite{chandrasekhar_dynamical_1943} formula to compute the acceleration
    of a satellite at radius $r$ from the center of the host due to dynamical friction:
    \begin{equation}
    \mathbf{a}_{DF} = -\frac{4\pi \mathrm{G}^2M_\mathrm{sat}\rho_\mathrm{host}(r)}{v_\mathrm{sat}^3}\ln
    \Lambda\left[\mathrm{erf}(x)-\frac{2x}{\sqrt{\pi}}\exp(-x^2)\right]\mathbf{v}_\mathrm{sat},
    \end{equation}
    where $M_\mathrm{sat}$ and $\mathbf{v}_\mathrm{sat}$ are the satellite's mass and velocity, respectively,
    $v_\mathrm{sat}=|\mathbf{v}_\mathrm{sat}|$, $\rho_\mathrm{host}(r)$ is the host's density profile,
    $\ln\Lambda=${\normalfont \ttfamily [logarithmCoulomb]} is the Coulomb logarithm, and $x\equiv
    v_\mathrm{sat}/\sqrt{2}\sigma(r)$, where $\sigma(r)$ is the velocity dispersion of the host halo at radius $r$,
    approximated to be equal to the host virial velocity, $v_\mathrm{vir}$.
   </description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionClass) :: satelliteDynamicalFrictionChandrasekhar1943
     !!{
     Implementation of a satellite dynamical friction class which uses the model of \cite{chandrasekhar_dynamical_1943}.
     !!}
     private
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     double precision                                     :: logarithmCoulomb
   contains
     !![
     <methods>
       <method description="Compute the Coulomb logarithm, $\log \Lambda$, appearing in the \cite{chandrasekhar_dynamical_1943} dynamical friction equation." method="coulombLogarithm" />
     </methods>
     !!]
     final     ::                     chandrasekhar1943Destructor
     procedure :: acceleration     => chandrasekhar1943Acceleration
     procedure :: coulombLogarithm => chandrasekhar1943CoulombLogarithm
  end type satelliteDynamicalFrictionChandrasekhar1943

  interface satelliteDynamicalFrictionChandrasekhar1943
     !!{
     Constructors for the chandrasekhar1943 satellite dynamical friction class.
     !!}
     module procedure chandrasekhar1943ConstructorParameters
     module procedure chandrasekhar1943ConstructorInternal
  end interface satelliteDynamicalFrictionChandrasekhar1943

contains

  function chandrasekhar1943ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{satelliteDynamicalFrictionChandrasekhar1943} satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteDynamicalFrictionChandrasekhar1943)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                   ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                  ), pointer       :: darkMatterProfileDMO_
    double precision                                                             :: logarithmCoulomb

    !![
    <inputParameter>
      <name>logarithmCoulomb</name>
      <defaultValue>2.0d0</defaultValue>
      <description>The Coulomb logarithm, $\ln \Lambda$, appearing in the \cite{chandrasekhar_dynamical_1943} formulation of the acceleration due to dynamical friction.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=satelliteDynamicalFrictionChandrasekhar1943(logarithmCoulomb,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function chandrasekhar1943ConstructorParameters

  function chandrasekhar1943ConstructorInternal(logarithmCoulomb,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{satelliteDynamicalFrictionChandrasekhar1943} satellite dynamical friction class.
    !!}
    implicit none
    type            (satelliteDynamicalFrictionChandrasekhar1943)                        :: self
    class           (darkMatterHaloScaleClass                   ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                  ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                             , intent(in   )         :: logarithmCoulomb
    !![
    <constructorAssign variables="logarithmCoulomb, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function chandrasekhar1943ConstructorInternal

  subroutine chandrasekhar1943Destructor(self)
    !!{
    Destructor for the simple cooling radius class.
    !!}
    implicit none
    type(satelliteDynamicalFrictionChandrasekhar1943), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine chandrasekhar1943Destructor

  function chandrasekhar1943Acceleration(self,node)
    !!{
    Return an acceleration for satellites due to dynamical friction using the formulation of \cite{chandrasekhar_dynamical_1943}.
    !!}
    use :: Coordinates                     , only : coordinateCartesian      , assignment(=)
    use :: Error_Functions                 , only : Error_Function
    use :: Galactic_Structure_Options      , only : coordinateSystemCartesian, componentTypeDarkHalo         , massTypeDark
    use :: Galacticus_Nodes                , only : nodeComponentSatellite   , nodeComponentBasic            , treeNode
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gigaYear                 , gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    double precision                                             , dimension(3)          :: chandrasekhar1943Acceleration
    class           (satelliteDynamicalFrictionChandrasekhar1943), intent(inout), target :: self
    type            (treeNode                                   ), intent(inout)         :: node
    class           (nodeComponentBasic                         ), pointer               :: basic
    class           (nodeComponentSatellite                     ), pointer               :: satellite
    type            (treeNode                                   ), pointer               :: nodeHost
    class           (massDistributionClass                      ), pointer               :: massDistribution_            , massDistributionHost_
    double precision                                                                     :: massSatellite
    type            (coordinateCartesian                        )                        :: position                     , velocity

    nodeHost                      =>  node     %mergesWith      ()
    basic                         =>  node     %basic           ()
    satellite                     =>  node     %satellite       ()
    massSatellite                 =   satellite%boundMass       ()
    position                      =   satellite%position        ()
    velocity                      =   satellite%velocity        ()
    massDistribution_             =>  node     %massDistribution()
    massDistributionHost_         =>  nodeHost %massDistribution()
    chandrasekhar1943Acceleration =  +4.0d0                                                                                                                &
            &                        *Pi                                                                                                                   &
            &                        *massDistributionHost_%chandrasekharIntegral(massDistributionHost_,massDistribution_,massSatellite,position,velocity) &
            &                        *self                 %coulombLogarithm     (node                                                                   ) &
            &                        *gravitationalConstant_internal**2                                                                                    &
            &                        *massSatellite                                                                                                        &
            &                        *kilo                                                                                                                 &
            &                        *gigaYear                                                                                                             &
            &                        /megaParsec
    !![
    <objectDestructor name="massDistribution_"    />
    <objectDestructor name="massDistributionHost_"/>
    !!]
    return
  end function chandrasekhar1943Acceleration

  double precision function chandrasekhar1943CoulombLogarithm(self,node) result(coulombLogarithm)
    !!{
    Evaluate the Coulomb logarithm for the \cite{chandrasekhar_dynamical_1943} dynamical friction model.
    !!}
    implicit none
    class(satelliteDynamicalFrictionChandrasekhar1943), intent(inout) :: self
    type (treeNode                                   ), intent(inout) :: node

    coulombLogarithm=self%logarithmCoulomb
    return
  end function chandrasekhar1943CoulombLogarithm

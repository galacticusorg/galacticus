!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  A dark matter halo profile heating class which computes heating due to two-body relaxation.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingTwoBodyRelaxation">
   <description>A dark matter profile heating model which computes heating due to two-body relaxation.</description>
  </darkMatterProfileHeating>
  !!]

  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingTwoBodyRelaxation
     !!{
     A dark matter profile heating class which computes heating due to two-body relaxation.
     !!}
     private
     double precision :: massParticle, lengthSoftening, &
          &              timeStart   , efficiency
   contains
     procedure :: specificEnergy                 => twoBodyRelaxationSpecificEnergy
     procedure :: specificEnergyGradient         => twoBodyRelaxationSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero => twoBodyRelaxationSpecificEnergyIsEverywhereZero
  end type darkMatterProfileHeatingTwoBodyRelaxation

  interface darkMatterProfileHeatingTwoBodyRelaxation
     !!{
     Constructors for the {\normalfont \ttfamily twoBodyRelaxation} dark matter profile heating class.
     !!}
     module procedure twoBodyRelaxationConstructorParameters
     module procedure twoBodyRelaxationConstructorInternal
  end interface darkMatterProfileHeatingTwoBodyRelaxation

contains

  function twoBodyRelaxationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily twoBodyRelaxation} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileHeatingTwoBodyRelaxation), target        :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    double precision                                                           :: massParticle, lengthSoftening, &
         &                                                                        timeStart   , efficiency

    !![
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <description>The particle mass to use for two-body relaxation calculations.</description>
    </inputParameter>
    <inputParameter>
      <name>lengthSoftening</name>
      <source>parameters</source>
      <description>The softening length to use for two-body relaxation calculations.</description>
    </inputParameter>
    <inputParameter>
      <name>timeStart</name>
      <source>parameters</source>
      <description>The time at which two-body relaxation is assumed to have begun.</description>
    </inputParameter>
    <inputParameter>
      <name>efficiency</name>
      <source>parameters</source>
      <description>The fractional efficiency of two-body relaxation heating.</description>
    </inputParameter>
    !!]
    self=darkMatterProfileHeatingTwoBodyRelaxation(massParticle,lengthSoftening,timeStart,efficiency)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function twoBodyRelaxationConstructorParameters

  function twoBodyRelaxationConstructorInternal(massParticle,lengthSoftening,timeStart,efficiency) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily twoBodyRelaxation} dark matter profile heating scales class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingTwoBodyRelaxation)                :: self
    double precision                                           , intent(in   ) :: massParticle, lengthSoftening, &
         &                                                                        timeStart   , efficiency
    !![
    <constructorAssign variables="massParticle, lengthSoftening, timeStart, efficiency"/>
    !!]

    return
  end function twoBodyRelaxationConstructorInternal

  double precision function twoBodyRelaxationSpecificEnergy(self,node,darkMatterProfileDMO_,radius)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}. The assumption here is that the mean
    fractional change in energy for a particle per crossing time is $8 \log \Lambda / N$ where $N$ is the number of particles
    within radius $r=${\normalfont \ttfamily radius}. The crossing time is approximated by $r/V(r)$ where $V(r)$ is the
    circular velocity at $r$. The Coulomb logarithm is given by $\log\Lambda=\hbox{max}(\epsilon,b_{90})$ where $\epsilon$ is
    the softening length, $b_{90}=2\mathrm{G}m_\mathrm{p}/V^2(r)$, and $m_\mathrm{p}$ is the particle mass. Finally, the
    specific energy is assumed to be $\sigma^2(r)/2\approx V^2(r)/4$.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear          , megaParsec, gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (darkMatterProfileHeatingTwoBodyRelaxation), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    class           (darkMatterProfileDMOClass                ), intent(inout) :: darkMatterProfileDMO_
    double precision                                           , intent(in   ) :: radius
    class           (nodeComponentBasic                       ), pointer       :: basic
    double precision                                                           :: particleCount     , velocity               , &
         &                                                                        logarithmCoulomb  , impactParameterCritical

    basic => node%basic()
    if (basic%time() > self%timeStart) then
       velocity                       =+darkMatterProfileDMO_             %circularVelocity(node,radius)
       impactParameterCritical        =+2.0d0                                                                &
            &                          *gravitationalConstantGalacticus                                      &
            &                          *self                              %massParticle                      &
            &                          /velocity                                                        **2
       logarithmCoulomb               =+0.5d0                                                                &
            &                          *log(                                                                 &
            &                               +1.0d0                                                           &
            &                               +(                                                               &
            &                                 +radius                                                        &
            &                                 /max(                                                          &
            &                                      self                   %lengthSoftening                 , &
            &                                      impactParameterCritical                                   &
            &                                     )                                                          &
            &                                )                                                          **2  &
            &                              )
       particleCount                  =+darkMatterProfileDMO_             %enclosedMass    (node,radius)     &
            &                          /self                              %massParticle
       twoBodyRelaxationSpecificEnergy=+2.0d0                                                                &
            &                          *self                              %efficiency                        &
            &                          *logarithmCoulomb                                                     &
            &                          *(                                                                    &
            &                            +basic                           %time            (           )     &
            &                            -self                            %timeStart                         &
            &                           )                                                                    &
            &                          *velocity                                                        **3  &
            &                          /radius                                                               &
            &                          /particleCount                                                        &
            &                          *kilo                                                                 &
            &                          *gigaYear                                                             &
            &                          /megaParsec
    else
       twoBodyRelaxationSpecificEnergy=+0.0d0
    end if
    return
  end function twoBodyRelaxationSpecificEnergy

  double precision function twoBodyRelaxationSpecificEnergyGradient(self,node,darkMatterProfileDMO_,radius)
    !!{
    Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic             , treeNode
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileHeatingTwoBodyRelaxation), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    class           (darkMatterProfileDMOClass                ), intent(inout) :: darkMatterProfileDMO_
    double precision                                           , intent(in   ) :: radius
    class           (nodeComponentBasic                       ), pointer       :: basic
    double precision                                                           :: particleCount     , velocity               , &
         &                                                                        logarithmCoulomb  , impactParameterCritical, &
         &                                                                        gradientCoulomb

    basic => node%basic()
    if (basic%time() > self%timeStart) then
       velocity                       =+darkMatterProfileDMO_             %circularVelocity(node,radius)
       impactParameterCritical        =+2.0d0                                                                &
            &                          *gravitationalConstantGalacticus                                      &
            &                          *self                              %massParticle                      &
            &                          /velocity                                                       **2
       logarithmCoulomb               =+0.5d0                                                                &
            &                          *log(                                                                 &
            &                               +1.0d0                                                           &
            &                               +(                                                               &
            &                                 +radius                                                        &
            &                                 /max(                                                          &
            &                                      self                   %lengthSoftening                 , &
            &                                      impactParameterCritical                                   &
            &                                     )                                                          &
            &                                )                                                          **2  &
            &                              )
       particleCount                  =+darkMatterProfileDMO_             %enclosedMass    (node,radius)     &
            &                          /self                              %massParticle
       if (self%lengthSoftening > impactParameterCritical) then
          gradientCoulomb=+radius                                                            &
               &          /self                             %lengthSoftening
       else
          gradientCoulomb=+2.0d0                                                             &
               &          *radius                                                            &
               &          /impactParameterCritical                                           &
               &          *(                                                                 &
               &            -1.0d0                                                           &
               &            +8.0d0                                                           &
               &            *Pi                                                              &
               &            *gravitationalConstantGalacticus                                 &
               &            *radius                                                      **2 &
               &            *darkMatterProfileDMO_          %density        (node,radius)    &
               &            /velocity                                                    **2 &
               &           )
       end if
       twoBodyRelaxationSpecificEnergyGradient=+self                             %specificEnergy(node,darkMatterProfileDMO_,radius)    &
            &                                  /                                                                            radius     &
            &                                  *(                                                                                      &
            &                                    -2.5d0                                                                                &
            &                                    +6.0d0                                                                                &
            &                                    *Pi                                                                                   &
            &                                    *gravitationalConstantGalacticus                                                      &
            &                                    *darkMatterProfileDMO_          %density       (node,                      radius)    &
            &                                    *radius                                                                           **2 &
            &                                    /velocity                                                                         **2 &
            &                                    -                                                                          radius     &
            &                                    *                gradientCoulomb                                                      &
            &                                    /               logarithmCoulomb                                                      &
            &                                    *sqrt(exp(2.0d0*logarithmCoulomb)-1.0d0)                                              &
            &                                    /     exp(2.0d0*logarithmCoulomb)                                                     &
            &                                   )
    else
       twoBodyRelaxationSpecificEnergyGradient=0.0d0
    end if
    return
  end function twoBodyRelaxationSpecificEnergyGradient

  logical function twoBodyRelaxationSpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !!{
    Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(darkMatterProfileHeatingTwoBodyRelaxation), intent(inout) :: self
    type (treeNode                                 ), intent(inout) :: node
    class(darkMatterProfileDMOClass                ), intent(inout) :: darkMatterProfileDMO_
    class(nodeComponentBasic                       ), pointer       :: basic
    !$GLC attributes unused :: self, darkMatterProfileDMO_

    basic                                           => node %basic()
    twoBodyRelaxationSpecificEnergyIsEverywhereZero =  basic%time () <= self%timeStart
    return
  end function twoBodyRelaxationSpecificEnergyIsEverywhereZero

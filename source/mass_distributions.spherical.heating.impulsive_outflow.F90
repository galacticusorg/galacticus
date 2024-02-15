!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implements a mass distribution heating class that computes heating due to two-body relaxation.
  !!}

  !![
  <massDistributionHeating name="massDistributionHeatingImpulsiveOutflow">
    <description>A mass distribution heating class that computes heating due to impulsive outflows.</description>
  </massDistributionHeating>
  !!]
  type, extends(massDistributionHeatingClass) :: massDistributionHeatingImpulsiveOutflow
     !!{
     Implementation of a mass distribution heating class that computes heating due to impulsive outflows.
     !!}
     private
     double precision :: energyImpulsiveOutflowDisk, energyImpulsiveOutflowSpheroid, &
          &              impulsiveEnergyFactor
   contains
     procedure :: specificEnergy                 => impulsiveOutflowSpecificEnergy
     procedure :: specificEnergyGradient         => impulsiveOutflowSpecificEnergyGradient
     procedure :: specificEnergyIsEveryWhereZero => impulsiveOutflowSpecificEnergyIsEverywhereZero
  end type massDistributionHeatingImpulsiveOutflow

  interface massDistributionHeatingImpulsiveOutflow
     !!{
     Constructors for the {\normalfont \ttfamily impulsiveOutflow} mass distribution class.
     !!}
     module procedure impulsiveOutflowConstructorParameters
     module procedure impulsiveOutflowConstructorInternal
  end interface massDistributionHeatingImpulsiveOutflow

contains

  function impulsiveOutflowConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily impulsiveOutflow} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionHeatingImpulsiveOutflow)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    double precision                                                         :: energyImpulsiveOutflowDisk, energyImpulsiveOutflowSpheroid, &
          &                                                                     impulsiveEnergyFactor

    !![
    <inputParameter>
      <name>energyImpulsiveOutflowDisk</name>
      <source>parameters</source>
      <description>The impulsive energy of outflows from the disk.</description>
    </inputParameter>
    <inputParameter>
      <name>energyImpulsiveOutflowSpheroid</name>
      <source>parameters</source>
      <description>The impulsive energy of outflows from the spheroid.</description>
    </inputParameter>
    <inputParameter>
      <name>impulsiveEnergyFactor</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The parameter $\alpha$ appearing in the impulsive outflow heating rate.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=massDistributionHeatingImpulsiveOutflow(energyImpulsiveOutflowDisk,energyImpulsiveOutflowSpheroid,impulsiveEnergyFactor)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function impulsiveOutflowConstructorParameters
  
  function impulsiveOutflowConstructorInternal(energyImpulsiveOutflowDisk,energyImpulsiveOutflowSpheroid,impulsiveEnergyFactor) result(self)
    !!{
    Constructor for ``impulsiveOutflow'' dark matter profile heating class.
    !!}
    implicit none
    type             (massDistributionHeatingImpulsiveOutflow)                :: self
    double precision                                          , intent(in   ) :: energyImpulsiveOutflowDisk, energyImpulsiveOutflowSpheroid, &
         &                                                                       impulsiveEnergyFactor
    !![
    <constructorAssign variables="energyImpulsiveOutflowDisk, energyImpulsiveOutflowSpheroid, impulsiveEnergyFactor"/>
    !!]
 
    return
  end function impulsiveOutflowConstructorInternal

  double precision function impulsiveOutflowSpecificEnergy(self,radius,massDistribution_) result(energySpecific)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options      , only : componentTypeDisk              , componentTypeSpheroid
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionHeatingImpulsiveOutflow), intent(inout) :: self
    double precision                                         , intent(in   ) :: radius
    class           (massDistributionClass                  ), intent(inout) :: massDistribution_
    double precision                                                         :: massTotalDisk    , massTotalSpheroid   , &
         &                                                                      fractionMassDisk , fractionMassSpheroid
    
    massTotalDisk    =massDistribution_%massTotal(componentType=componentTypeDisk    )
    massTotalSpheroid=massDistribution_%massTotal(componentType=componentTypeSpheroid)
    if (massTotalDisk     > 0.0d0) then
       fractionMassDisk    =+massDistribution_%massEnclosedBySphere(radius,componentType=componentTypeDisk    ) &
            &               /                  massTotalDisk
    else
       fractionMassDisk    =+0.0d0
    end if
    if (massTotalSpheroid > 0.0d0) then
       fractionMassSpheroid=+massDistribution_%massEnclosedBySphere(radius,componentType=componentTypeSpheroid) &
            &               /                  massTotalSpheroid
    else
       fractionMassSpheroid=+0.0d0
    end if
    energySpecific=+  self%impulsiveEnergyFactor          &
         &         *gravitationalConstantGalacticus       &
         &         *(                                     &
         &           +self%energyImpulsiveOutflowDisk     &
         &           *fractionMassDisk                    &
         &           +self%energyImpulsiveOutflowSpheroid &
         &           *fractionMassSpheroid                &
         &          )                                     &
         &         /radius
    return
  end function impulsiveOutflowSpecificEnergy

  double precision function impulsiveOutflowSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating.
    !!}
    use :: Coordinates                     , only : coordinateSpherical            , assignment(=)
    use :: Galactic_Structure_Options      , only : componentTypeDisk              , componentTypeSpheroid
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionHeatingImpulsiveOutflow), intent(inout) :: self
    double precision                                         , intent(in   ) :: radius
    class           (massDistributionClass                  ), intent(inout) :: massDistribution_
    double precision                                                         :: massTotalDisk      , massTotalSpheroid      , &
         &                                                                      fractionMassDisk   , fractionMassSpheroid   , &
         &                                                                      fractionDensityDisk, fractionDensitySpheroid
    type            (coordinatespherical                    )                :: coordinates

    coordinates      =[radius,0.0d0,0.0d0]
    massTotalDisk    =massDistribution_%massTotal(componentType=componentTypeDisk    )
    massTotalSpheroid=massDistribution_%massTotal(componentType=componentTypeSpheroid)
    if (massTotalDisk     > 0.0d0) then
       fractionMassDisk       =+massDistribution_%massEnclosedBySphere(radius     ,componentType=componentTypeDisk    ) &
            &                  /                  massTotalDisk
       fractionDensityDisk    =+massDistribution_%density             (coordinates,componentType=componentTypeDisk    ) &
            &                  /                  massTotalDisk
    else
       fractionMassDisk       =+0.0d0
       fractionDensityDisk    =+0.0d0
    end if
    if (massTotalSpheroid > 0.0d0) then
       fractionMassSpheroid   =+massDistribution_%massEnclosedBySphere(radius     ,componentType=componentTypeSpheroid) &
            &                  /                  massTotalSpheroid
       fractionDensitySpheroid=+massDistribution_%density             (coordinates,componentType=componentTypeSpheroid) &
            &                  /                  massTotalSpheroid
    else
       fractionMassSpheroid   =+0.0d0
       fractionDensitySpheroid=+0.0d0
    end if
    energySpecificGradient=+self%impulsiveEnergyFactor              &
         &                 *gravitationalConstantGalacticus         &
         &                 *(                                       &
         &                   +(                                     &
         &                     +self%energyImpulsiveOutflowDisk     &
         &                     *fractionDensityDisk                 &
         &                     +self%energyImpulsiveOutflowSpheroid &
         &                     *fractionDensitySpheroid             &
         &                    )                                     &
         &                   *4.0d0                                 &
         &                   *Pi                                    &
         &                   *radius                                &
         &                   -(                                     &
         &                     +self%energyImpulsiveOutflowDisk     &
         &                     *fractionMassDisk                    &
         &                     +self%energyImpulsiveOutflowSpheroid &
         &                     *fractionMassSpheroid                &
         &                    )                                     &
         &                   /radius**2                             &
         &                  )
    return
  end function impulsiveOutflowSpecificEnergyGradient

  logical function impulsiveOutflowSpecificEnergyIsEverywhereZero(self) result(energySpecificIsEverywhereZero)
    !!{
    Returns true if the specific energy is everywhere zero.
    !!}
    implicit none
    class(massDistributionHeatingImpulsiveOutflow), intent(inout) :: self

    energySpecificIsEverywhereZero= self%energyImpulsiveOutflowDisk     <= 0.0d0 &
         &                         .and.                                         &
         &                          self%energyImpulsiveOutflowSpheroid <= 0.0d0
    return
  end function impulsiveOutflowSpecificEnergyIsEverywhereZero

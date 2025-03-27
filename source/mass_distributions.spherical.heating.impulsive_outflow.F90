!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
    <description>
      A mass distribution heating class that computes heating due to impulsive outflows---i.e. outflows occurring on
      timescales that are small relative to the dynamical time of the halo. The model assumed is that the energy injection is given by
      \begin{equation}
    \dot{\epsilon}(r) = \alpha \frac{\mathrm{G} \dot{M}_\mathrm{outflow}(r)}{r} f\left( \frac{t_\phi}{t_\mathrm{dyn}} \right),
    \end{equation}
      where $\alpha$ is a normalization factor, $t_\phi = M_\mathrm{gas}/\dot{M}_\mathrm{outflow}$ is the timescale for the
      outflow, and $t_\mathrm{dyn} = r_{1/2}/v_{1/2}$ is the dynamical time at the half-mass radius.
      
      The quantity
      \begin{equation}
    \dot{\epsilon}^\prime = \dot{M}_\mathrm{outflow} f\left( \frac{t_\phi}{t_\mathrm{dyn}} \right),
    \end{equation}
      if provided as an argument to the class constructor.
    </description>
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
    Constructor for {\normalfont \ttfamily impulsiveOutflow} dark matter profile heating class.
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
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , componentTypeSpheroid
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionHeatingImpulsiveOutflow), intent(inout) :: self
    double precision                                         , intent(in   ) :: radius
    class           (massDistributionClass                  ), intent(inout) :: massDistribution_
    class           (massDistributionClass                  ), pointer       :: massDistributionDisk, massDistributionSpheroid
    double precision                                                         :: massTotalDisk       , massTotalSpheroid       , &
         &                                                                      fractionMassDisk    , fractionMassSpheroid
    
    massDistributionDisk     => massDistribution_%subset(componentType=componentTypeDisk    )
    massDistributionSpheroid => massDistribution_%subset(componentType=componentTypeSpheroid)
    fractionMassDisk         =  0.0d0
    fractionMassSpheroid     =  0.0d0
    if (associated(massDistributionDisk    )) then
       massTotalDisk            =  massDistributionDisk    %massTotal()
       if (massTotalDisk     > 0.0d0) then
          fractionMassDisk    =+massDistributionDisk    %massEnclosedBySphere(radius) &
               &               /                         massTotalDisk
       end if
    end if
    if (associated(massDistributionSpheroid)) then
       massTotalSpheroid        =  massDistributionSpheroid%massTotal()
       if (massTotalSpheroid > 0.0d0) then
          fractionMassSpheroid=+massDistributionSpheroid%massEnclosedBySphere(radius) &
               &               /                         massTotalSpheroid
       end if
    end if
    energySpecific=+  self%impulsiveEnergyFactor          &
         &         *gravitationalConstant_internal        &
         &         *(                                     &
         &           +self%energyImpulsiveOutflowDisk     &
         &           *fractionMassDisk                    &
         &           +self%energyImpulsiveOutflowSpheroid &
         &           *fractionMassSpheroid                &
         &          )                                     &
         &         /radius
    !![
    <objectDestructor name="massDistributionDisk"    />
    <objectDestructor name="massDistributionSpheroid"/>
    !!]
    return
  end function impulsiveOutflowSpecificEnergy

  double precision function impulsiveOutflowSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Galactic_Structure_Options      , only : componentTypeDisk             , componentTypeSpheroid
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionHeatingImpulsiveOutflow), intent(inout) :: self
    double precision                                         , intent(in   ) :: radius
    class           (massDistributionClass                  ), intent(inout) :: massDistribution_
    class           (massDistributionClass                  ), pointer       :: massDistributionDisk, massDistributionSpheroid
    double precision                                                         :: massTotalDisk       , massTotalSpheroid       , &
         &                                                                      fractionMassDisk    , fractionMassSpheroid    , &
         &                                                                      fractionDensityDisk , fractionDensitySpheroid
    type            (coordinatespherical                    )                :: coordinates

    massDistributionDisk     => massDistribution_%subset(componentType=componentTypeDisk    )
    massDistributionSpheroid => massDistribution_%subset(componentType=componentTypeSpheroid)
    coordinates              =  [radius,0.0d0,0.0d0]
    fractionMassDisk         =  0.0d0
    fractionMassSpheroid     =  0.0d0
    fractionDensityDisk      =  0.0d0
    fractionDensitySpheroid  =  0.0d0
    if (associated(massDistributionDisk    )) then
       massTotalDisk            =  massDistributionDisk    %massTotal()
       if (massTotalDisk     > 0.0d0) then
          fractionMassDisk       =+massDistributionDisk    %massEnclosedBySphere(radius     ) &
               &                  /                         massTotalDisk
          fractionDensityDisk    =+massDistributionDisk    %density             (coordinates) &
               &                  /                         massTotalDisk
       end if
    end if
    if (associated(massDistributionSpheroid)) then
       massTotalSpheroid        =  massDistributionSpheroid%massTotal()
       if (massTotalSpheroid > 0.0d0) then
          fractionMassSpheroid   =+massDistributionSpheroid%massEnclosedBySphere(radius     ) &
               &                  /                         massTotalSpheroid
          fractionDensitySpheroid=+massDistributionSpheroid%density             (coordinates) &
               &                  /                         massTotalSpheroid
       end if
    end if
    energySpecificGradient=+self%impulsiveEnergyFactor              &
         &                 *gravitationalConstant_internal          &
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
    !![
    <objectDestructor name="massDistributionDisk"    />
    <objectDestructor name="massDistributionSpheroid"/>
    !!]
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

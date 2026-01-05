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
  A dark matter halo profile heating class which accounts for heating arising from impulsive outflows.
  !!}
  
  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingImpulsiveOutflow">
   <description>
    A dark matter profile heating model which accounts for heating due to impulsive outflows. The quantity
    \begin{equation}
     \dot{\epsilon}^\prime = \dot{M}_\mathrm{outflow} f\left( \frac{t_\phi}{t_\mathrm{dyn}} \right),
    \end{equation}    
    has been accumulated by the \refClass{nodeOperatorImpulsiveOutflowEnergy} object---radially-dependent factors are then applied
    in the \refClass{massDistributionHeatingImpulsiveOutflow} object returned from our factory.
   </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingImpulsiveOutflow
     !!{
     A dark matter profile heating class which accounts for heating arising from impulsive outflows.
     !!}
     private
     integer          :: energyImpulsiveOutflowDiskID, energyImpulsiveOutflowSpheroidID
     double precision :: impulsiveEnergyFactor
   contains
     procedure :: get  => impulsiveOutflowGet
  end type darkMatterProfileHeatingImpulsiveOutflow

  interface darkMatterProfileHeatingImpulsiveOutflow
     !!{
     Constructors for the \refClass{darkMatterProfileHeatingImpulsiveOutflow} dark matter profile heating class.
     !!}
     module procedure impulsiveOutflowConstructorParameters
     module procedure impulsiveOutflowConstructorInternal
  end interface darkMatterProfileHeatingImpulsiveOutflow

contains

  function impulsiveOutflowConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileHeatingImpulsiveOutflow} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileHeatingImpulsiveOutflow), target        :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    double precision                                                          :: impulsiveEnergyFactor

    !![
    <inputParameter>
      <name>impulsiveEnergyFactor</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The parameter $\alpha$ appearing in the impulsive outflow heating rate.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=darkMatterProfileHeatingImpulsiveOutflow(impulsiveEnergyFactor)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function impulsiveOutflowConstructorParameters

  function impulsiveOutflowConstructorInternal(impulsiveEnergyFactor) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileHeatingImpulsiveOutflow} dark matter profile heating scales class.
    !!}
    implicit none
    type            (darkMatterProfileHeatingImpulsiveOutflow)                :: self
    double precision                                          , intent(in   ) :: impulsiveEnergyFactor
    !![
    <constructorAssign variables="impulsiveEnergyFactor"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="energyImpulsiveOutflowDisk"     id="self%energyImpulsiveOutflowDiskID"     isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="energyImpulsiveOutflowSpheroid" id="self%energyImpulsiveOutflowSpheroidID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function impulsiveOutflowConstructorInternal

  function impulsiveOutflowGet(self,node) result(massDistributionHeating_)
    !!{
    Return the dark matter mass distribution heating for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentDarkMatterProfile
    use :: Mass_Distributions, only : massDistributionHeatingImpulsiveOutflow
    implicit none
    class(massDistributionHeatingClass            ), pointer       :: massDistributionHeating_
    class(darkMatterProfileHeatingImpulsiveOutflow), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile
 
    ! Create the mass distribution.
    allocate(massDistributionHeatingImpulsiveOutflow :: massDistributionHeating_)
    select type(massDistributionHeating_)
    type is (massDistributionHeatingImpulsiveOutflow)
       darkMatterProfile => node%darkMatterProfile()
       !![
       <referenceConstruct object="massDistributionHeating_">
	 <constructor>
           massDistributionHeatingImpulsiveOutflow(                                                                                                                   &amp;
           &amp;                                   energyImpulsiveOutflowDisk    =darkMatterProfile%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowDiskID    ), &amp;
           &amp;                                   energyImpulsiveOutflowSpheroid=darkMatterProfile%floatRank0MetaPropertyGet(self%energyImpulsiveOutflowSpheroidID), &amp;
           &amp;                                   impulsiveEnergyFactor         =self             %impulsiveEnergyFactor                                             &amp;
           &amp;                                  )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function impulsiveOutflowGet

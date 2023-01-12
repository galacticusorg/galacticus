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
  Implementation of an active mass for star formation class in which the mass of the ISM above a surface density threshold is active.
  !!}
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Galactic_Structure      , only : galacticStructureClass
  use :: Math_Exponentiation     , only : fastExponentiator

  !![
  <starFormationActiveMass name="starFormationActiveMassSurfaceDensityThreshold">
   <description>An active mass for star formation class in which the mass of the ISM above a surface density threshold is active.</description>
  </starFormationActiveMass>
  !!]
  type, extends(starFormationActiveMassClass) :: starFormationActiveMassSurfaceDensityThreshold
     !!{
     Implementation of an active mass for star formation class in which the mass of the ISM above a surface density threshold is active.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_   => null()
     class           (galacticStructureClass   ), pointer :: galacticStructure_      => null()
     double precision                                     :: surfaceDensityThreshold          , surfaceDensityNormalization, &
          &                                                  exponentVelocity
     type            (fastExponentiator        )          :: velocityExponentiator
   contains
     final :: surfaceDensityThresholdDestructor
     procedure :: massActive => surfaceDensityThresholdMassActive
  end type starFormationActiveMassSurfaceDensityThreshold

  interface starFormationActiveMassSurfaceDensityThreshold
     !!{
     Constructors for the {\normalfont \ttfamily surfaceDensityThreshold} active mass for star formation class.
     !!}
     module procedure surfaceDensityThresholdConstructorParameters
     module procedure surfaceDensityThresholdConstructorInternal
  end interface starFormationActiveMassSurfaceDensityThreshold

contains

  function surfaceDensityThresholdConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily surfaceDensityThreshold} active mass for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationActiveMassSurfaceDensityThreshold)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass                     ), pointer       :: darkMatterProfileDMO_
    class           (galacticStructureClass                        ), pointer       :: galacticStructure_
    double precision                                                                :: surfaceDensityThreshold, exponentVelocity

    !![
    <inputParameter>
      <name>surfaceDensityThreshold</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The surface density threshold above which ISM gas participates in star formation.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponentVelocity</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Tne exponent of velocity in the surface density threshold for star formation.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="galacticStructure"    name="galacticStructure_"    source="parameters"/>
    !!]
    self=starFormationActiveMassSurfaceDensityThreshold(surfaceDensityThreshold,exponentVelocity,darkMatterProfileDMO_,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="galacticStructure_"   />
    !!]
    return
  end function surfaceDensityThresholdConstructorParameters

  function surfaceDensityThresholdConstructorInternal(surfaceDensityThreshold,exponentVelocity,darkMatterProfileDMO_,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily surfaceDensityThreshold} active mass for star formation class.
    !!}
    implicit none
    type            (starFormationActiveMassSurfaceDensityThreshold)                        :: self
    double precision                                                , intent(in   )         :: surfaceDensityThreshold        , exponentVelocity
    class           (darkMatterProfileDMOClass                     ), intent(in   ), target :: darkMatterProfileDMO_
    class           (galacticStructureClass                        ), intent(in   ), target :: galacticStructure_
    double precision                                                , parameter             :: velocityNormalization  =200.0d0
    !![
    <constructorAssign variables="surfaceDensityThreshold, exponentVelocity, *darkMatterProfileDMO_, *galacticStructure_"/>
    !!]
    
    ! Initialize exponentiators.
    self%velocityExponentiator=fastExponentiator(1.0d+0,1.0d+3,exponentVelocity,1.0d+1,abortOutsideRange=.false.)
    ! Compute normalization factor for surface density.
    self%surfaceDensityNormalization=surfaceDensityThreshold/velocityNormalization**exponentVelocity
    return
  end function surfaceDensityThresholdConstructorInternal

  subroutine surfaceDensityThresholdDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily surfaceDensityThreshold} active mass for star formation class.
    !!}
    implicit none
    type(starFormationActiveMassSurfaceDensityThreshold), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%galacticStructure_"  />
    !!]
    return
  end subroutine surfaceDensityThresholdDestructor

  double precision function surfaceDensityThresholdMassActive(self,component)
    !!{
    Returns the mass (in $\mathrm{M}_\odot$) of gas actively undergoing star formation in the given {\normalfont \ttfamily
    component} as the mass of gas in the ISM above a given surface density threshold
    !!}
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : nodeComponentDisk
    use :: Galactic_Structure_Options, only : componentTypeDisk, coordinateSystemCartesian, coordinateSystemCylindrical, massTypeGaseous, &
          &                                   weightByMass     , weightIndexNull
    implicit none
    class           (starFormationActiveMassSurfaceDensityThreshold), intent(inout) :: self
    class           (nodeComponent                                 ), intent(inout) :: component
    double precision                                                                :: surfaceDensityThreshold, radiusBounding, &
         &                                                                             densitySurfaceCentral
    
    select type (component)
    class is (nodeComponentDisk)
       ! Compute the surface density threshold for this node.
       surfaceDensityThreshold=+self%surfaceDensityNormalization                                                                                       &
            &                  *self%velocityExponentiator      %exponentiate(self%darkMatterProfileDMO_ %circularVelocityMaximum(component%hostNode))
       ! We assume a monotonically decreasing surface density. So, if the central density is below threshold then the active mass
       ! is zero.
       densitySurfaceCentral               =self%galacticStructure_%surfaceDensity              (component%hostNode,[0.0d0,0.0d0,0.0d0]    ,coordinateSystemCylindrical,componentTypeDisk,massTypeGaseous,weightByMass,weightIndexNull)
       if (densitySurfaceCentral < surfaceDensityThreshold) then
          surfaceDensityThresholdMassActive=0.0d0
       else
          radiusBounding                   =self%galacticStructure_%radiusEnclosingSurfaceDensity(component%hostNode,surfaceDensityThreshold                            ,componentTypeDisk,massTypeGaseous,weightByMass,weightIndexNull)
          surfaceDensityThresholdMassActive=self%galacticStructure_%massEnclosed                 (component%hostNode,radiusBounding                                     ,componentTypeDisk,massTypeGaseous,weightByMass,weightIndexNull)
       end if
    class default
       surfaceDensityThresholdMassActive=0.0d0
       call Error_Report('unsupported class'//{introspection:location})
    end select
    return
  end function surfaceDensityThresholdMassActive

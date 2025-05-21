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
  An implementation of isothermal dark matter halo profiles.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOIsothermal">
   <description>
    A dark matter profile DMO class in which builds \refClass{massDistributionIsothermal} objects to implement isothermal density
    profiles, normalized such that the total mass of the \gls{node} is enclosed with the virial radius.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOIsothermal
     !!{
     A dark matter halo profile class implementing isothermal dark matter halos.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::        isothermalDestructor
     procedure :: get => isothermalGet
  end type darkMatterProfileDMOIsothermal

  interface darkMatterProfileDMOIsothermal
     !!{
     Constructors for the \refClass{darkMatterProfileDMOIsothermal} dark matter halo profile class.
     !!}
     module procedure isothermalConstructorParameters
     module procedure isothermalConstructorInternal
  end interface darkMatterProfileDMOIsothermal

contains

  function isothermalConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily isothermal} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMOIsothermal)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOIsothermal(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function isothermalConstructorParameters

  function isothermalConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOIsothermal} dark matter halo profile class.
    !!}
    implicit none
    type (darkMatterProfileDMOIsothermal)                        :: self
    class(darkMatterHaloScaleClass      ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    return
  end function isothermalConstructorInternal

  subroutine isothermalDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOIsothermal} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOIsothermal), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine isothermalDestructor

  function isothermalGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo     , massTypeDark                    , weightByMass
    use :: Mass_Distributions        , only : massDistributionIsothermal, kinematicsDistributionIsothermal
    implicit none
    class           (massDistributionClass           ), pointer                 :: massDistribution_
    type            (kinematicsDistributionIsothermal), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOIsothermal  ), intent(inout)           :: self
    type            (treeNode                        ), intent(inout)           :: node
    type            (enumerationWeightByType         ), intent(in   ), optional :: weightBy
    integer                                           , intent(in   ), optional :: weightIndex
    class           (nodeComponentBasic              ), pointer                 :: basic
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionIsothermal :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionIsothermal)
       basic => node%basic()
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionIsothermal(                                                                        &amp;
           &amp;                      mass           =basic                     %mass                 (    ), &amp;
           &amp;                      lengthReference=self %darkMatterHaloScale_%radiusVirial         (node), &amp;
           &amp;                      componentType  =                           componentTypeDarkHalo      , &amp;
           &amp;                      massType       =                           massTypeDark                 &amp;
           &amp;                     )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionIsothermal(                                                                                &amp;
        &amp;                            velocityDispersion_=self %darkMatterHaloScale_%velocityVirial(node)/sqrt(2.0d0) &amp;
        &amp;                           )
	 </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function isothermalGet

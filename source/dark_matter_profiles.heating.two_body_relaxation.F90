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
  A dark matter halo profile heating class which computes heating due to two-body relaxation.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingTwoBodyRelaxation">
    <description>
      A dark matter profile heating class which returns a \refClass{massDistributionHeatingTwoBodyRelaxation} objects to compute
      heating due to two-body relaxation.
    </description>
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
     procedure :: get => twoBodyRelaxationGet
  end type darkMatterProfileHeatingTwoBodyRelaxation

  interface darkMatterProfileHeatingTwoBodyRelaxation
     !!{
     Constructors for the \refClass{darkMatterProfileHeatingTwoBodyRelaxation} dark matter profile heating class.
     !!}
     module procedure twoBodyRelaxationConstructorParameters
     module procedure twoBodyRelaxationConstructorInternal
  end interface darkMatterProfileHeatingTwoBodyRelaxation

contains

  function twoBodyRelaxationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileHeatingTwoBodyRelaxation} dark matter profile heating scales class which takes a parameter set as input.
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
    Internal constructor for the \refClass{darkMatterProfileHeatingTwoBodyRelaxation} dark matter profile heating scales class.
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

  function twoBodyRelaxationGet(self,node) result(massDistributionHeating_)
    !!{
    Return the dark matter mass distribution heating for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: Mass_Distributions, only : massDistributionHeatingTwoBodyRelaxation
    implicit none
    class           (massDistributionHeatingClass             ), pointer       :: massDistributionHeating_
    class           (darkMatterProfileHeatingTwoBodyRelaxation), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    class           (nodeComponentBasic                       ), pointer       :: basic
 
    ! Create the mass distribution.
    allocate(massDistributionHeatingTwoBodyRelaxation :: massDistributionHeating_)
    select type(massDistributionHeating_)
    type is (massDistributionHeatingTwoBodyRelaxation)
       basic => node%basic()
       !![
       <referenceConstruct object="massDistributionHeating_">
	 <constructor>
           massDistributionHeatingTwoBodyRelaxation(                                          &amp;
            &amp;                                   massParticle   =+self %massParticle     , &amp;
            &amp;                                   lengthSoftening=+self %lengthSoftening  , &amp;
            &amp;                                   timeRelaxing   =+basic%time           ()  &amp;
	    &amp;                                                   -self %timeStart        , &amp;
            &amp;                                   efficiency     =+self %efficiency         &amp;
            &amp;                                  )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function twoBodyRelaxationGet

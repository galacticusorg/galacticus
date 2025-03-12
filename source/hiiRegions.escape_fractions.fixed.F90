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
!+ Contributions to this file made by Sachi Weerasooriya
  !!{
  Implementation of a power law luminosity function which scales with a given exponent.
  !!}
  !![
  <hiiRegionEscapeFraction name="hiiRegionEscapeFractionFixed">
   <description>
    Escape fraction is the fraction of ionizing photons escaping into IGM. Here we assume fixed escape fraction for all galaxies.
   </description>
  </hiiRegionEscapeFraction>
  !!]
  type, extends(hiiRegionEscapeFractionClass) :: hiiRegionEscapeFractionFixed
     !!{
     Implementation of a fixed escape fraction
     !!}
     private
     double precision :: escapeFraction, ageLimit     
   contains
     procedure :: escapeFractionMethod => escapeFractionFixed
  end type hiiRegionEscapeFractionFixed

  interface hiiRegionEscapeFractionFixed
     !!{
     Constructors for the hiiRegionEscapeFractionFixed for HII region luminosity function class.
     !!}
     module procedure escapeFractionFixedConstructorParameters
     module procedure escapeFractionFixedConstructorInternal
  end interface hiiRegionEscapeFractionFixed

contains

  function escapeFractionFixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily dynamicalTime} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hiiRegionEscapeFractionFixed)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: escapeFraction, ageLimit
    !![
    <inputParameter>
      <name>escapeFraction</name>
      <defaultValue>0.006d0</defaultValue>
      <description> Escape fraction of ionizing photons from galaxies into IGM. </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ageLimit</name>
      <defaultValue>0.03d0</defaultValue>
      <description> Escape fraction of ionizing photons from galaxies into IGM. </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=hiiRegionEscapeFractionFixed(escapeFraction,ageLimit)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function escapeFractionFixedConstructorParameters

  function escapeFractionFixedConstructorInternal(escapeFraction,ageLimit) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily fixed escape fraction} 
    !!}
    
    implicit none
    type            (hiiRegionEscapeFractionFixed)               :: self
    double precision                                     , intent(in   ) :: escapeFraction, ageLimit
    !![
    <constructorAssign variables="escapeFraction"/>
    !!]
    self%escapeFraction=escapeFraction
    self%ageLimit=ageLimit
    return
  end function escapeFractionFixedConstructorInternal

  double precision function escapeFractionFixed(self,age_pop) result(escapeFrac)
    !!{
    Returns the escape fraction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (hiiRegionEscapeFractionFixed), intent(inout) :: self
    double precision                                       , intent(in   ) :: age_pop
    !type            (treeNode                      ), intent(inout), optional :: node
    !$GLC attributes unused :: node
    if(age_pop .ge. self%ageLimit) then
      escapeFrac=1.0d0
    else
      escapeFrac=self%escapeFraction
    end if
    return
  end function escapeFractionFixed

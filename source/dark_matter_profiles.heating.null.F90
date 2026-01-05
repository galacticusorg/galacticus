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
  A null dark matter halo profile heating class.
  !!}

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingNull">
    <description>
      A dark matter profile heating model which constructs \refClass{massDistributionHeatingNull} objects to provide zero heating.
    </description>
  </darkMatterProfileHeating>
  !!]

  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingNull
     !!{
     A dark matter profile heating class with zero heating.
     !!}
     private
   contains
     procedure :: get => nullGet
  end type darkMatterProfileHeatingNull

  interface darkMatterProfileHeatingNull
     !!{
     Constructors for the \refClass{darkMatterProfileHeatingNull} dark matter profile heating class.
     !!}
     module procedure nullConstructorParameters
  end interface darkMatterProfileHeatingNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileHeatingNull} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterProfileHeatingNull), target        :: self
    type(inputParameters             ), intent(inout) :: parameters

    self=darkMatterProfileHeatingNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  function nullGet(self,node) result(massDistributionHeating_)
    !!{
    Return the dark matter mass distribution heating for the given {\normalfont \ttfamily node}.
    !!}
    use :: Mass_Distributions, only : massDistributionHeatingNull
    implicit none
    class(massDistributionHeatingClass), pointer       :: massDistributionHeating_
    class(darkMatterProfileHeatingNull), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
 
    ! Create the mass distribution.
    allocate(massDistributionHeatingNull :: massDistributionHeating_)
    select type(massDistributionHeating_)
    type is (massDistributionHeatingNull)       
       !![
       <referenceConstruct object="massDistributionHeating_">
	 <constructor>
           massDistributionHeatingNull()
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function nullGet
  

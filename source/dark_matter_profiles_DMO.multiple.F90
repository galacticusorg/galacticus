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

!+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

  !!{
  An implementation of multiple dark matter halo profiles which allow different profiles for the host and the satellite.
  !!}

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOMultiple">
   <description>
    A dark matter profile DMO class in which the density profiles of the host halo and the satellite halo can be set separately
    to any other {\normalfont \ttfamily darkMatterProfileDMO} available.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOMultiple
     !!{
     A dark matter halo profile class implementing multiple dark matter halos which allow different profiles for the host and the satellite.
     !!}
     private
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMOHost_ => null(), darkMatterProfileDMOSatellite_ => null()
   contains
     final     ::        multipleDestructor
     procedure :: get => multipleGet
  end type darkMatterProfileDMOMultiple

  interface darkMatterProfileDMOMultiple
     !!{
     Constructors for the \refClass{darkMatterProfileDMOMultiple} dark matter halo profile class.
     !!}
     module procedure multipleConstructorParameters
     module procedure multipleConstructorInternal
  end interface darkMatterProfileDMOMultiple

contains

  function multipleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOMultiple} dark matter halo profile class which takes a parameter set as input.
    !!}
    implicit none
    type   (darkMatterProfileDMOMultiple)                :: self
    type   (inputParameters             ), intent(inout) :: parameters
    class  (darkMatterProfileDMOClass   ), pointer       :: darkMatterProfileDMOHost_, darkMatterProfileDMOSatellite_

    !![
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMOHost_"      parameterName="darkMatterProfileDMOHost"      source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMOSatellite_" parameterName="darkMatterProfileDMOSatellite" source="parameters"/>
    !!]
    self=darkMatterProfileDMOMultiple(darkMatterProfileDMOHost_,darkMatterProfileDMOSatellite_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMOHost_"     />
    <objectDestructor name="darkMatterProfileDMOSatellite_"/>
    !!]
    return
  end function multipleConstructorParameters

  function multipleConstructorInternal(darkMatterProfileDMOHost_,darkMatterProfileDMOSatellite_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOMultiple} dark matter profile class.
    !!}
    implicit none
    type            (darkMatterProfileDMOMultiple)                        :: self
    class           (darkMatterProfileDMOClass   ), intent(in   ), target :: darkMatterProfileDMOHost_, darkMatterProfileDMOSatellite_
    !![
    <constructorAssign variables="*darkMatterProfileDMOHost_, *darkMatterProfileDMOSatellite_"/>
    !!]

    return
  end function multipleConstructorInternal

  subroutine multipleDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOMultiple} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOMultiple), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMOHost_"     />
    <objectDestructor name="self%darkMatterProfileDMOSatellite_"/>
    !!]
    return
  end subroutine multipleDestructor

  function multipleGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class  (massDistributionClass       ), pointer                 :: massDistribution_
    class  (darkMatterProfileDMOMultiple), intent(inout)           :: self
    type   (treeNode                    ), intent(inout)           :: node
    type   (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                              , intent(in   ), optional :: weightIndex

    if (node%isSatellite()) then
       massDistribution_ => self%darkMatterProfileDMOSatellite_%get(node,weightBy,weightIndex)
    else
       massDistribution_ => self%darkMatterProfileDMOHost_     %get(node,weightBy,weightIndex)
    end if
    return
  end function multipleGet

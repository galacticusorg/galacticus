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
  Implements a test harness for spherical mass distributions which forces the use of numerical solvers.
  !!}

  !![
  <massDistribution name="massDistributionSphericalTestHarness">
   <description>A test harness for spherical mass distributions which forces the use of numerical solvers.</description>
  </massDistribution>
  !!]
  type, extends(massDistributionSpherical) :: massDistributionSphericalTestHarness
     !!{
     Implementation of a test harness for spherical mass distributions which forces the use of numerical solvers.
     !!}
     private
     class(massDistributionSpherical), pointer :: massDistribution_ => null()
   contains
     final     ::            sphericalTestHarnessDestructor
     procedure :: density => sphericalTestHarnessDensity
  end type massDistributionSphericalTestHarness

  interface massDistributionSphericalTestHarness
     !!{
     Constructors for the {\normalfont \ttfamily sphericalTestHarness} mass distribution class.
     !!}
     module procedure sphericalTestHarnessConstructorParameters
     module procedure sphericalTestHarnessConstructorInternal
  end interface massDistributionSphericalTestHarness

contains

  function sphericalTestHarnessConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sphericalTestHarness} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (massDistributionSphericalTestHarness)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(massDistributionClass               ), pointer       :: massDistribution_

    !![
    <objectBuilder class="massDistribution" name="massDistribution_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       self=massDistributionSphericalTestHarness(massDistribution_)
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function sphericalTestHarnessConstructorParameters
  
  function sphericalTestHarnessConstructorInternal(massDistribution_) result(self)
    !!{
    Constructor for ``sphericalTestHarness'' convergence class.
    !!}
    implicit none
    type (massDistributionSphericalTestHarness)                        :: self
    class(massDistributionSpherical           ), intent(in   ), target :: massDistribution_
    !![
    <constructorAssign variables="*massDistribution_"/>
    !!]
 
    self%componentType=self%massDistribution_%componentType
    self%     massType=self%massDistribution_%     massType
    return
  end function sphericalTestHarnessConstructorInternal

  subroutine sphericalTestHarnessDestructor(self)
    !!{
    Destructor for the ``sphericalTestHarness'' mass distribution class.
    !!}
    implicit none
    type(massDistributionSphericalTestHarness), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"/>
    !!]
    return
  end subroutine sphericalTestHarnessDestructor

  double precision function sphericalTestHarnessDensity(self,coordinates,componentType,massType) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalTestHarness), intent(inout)              :: self
    class(coordinate                          ), intent(in   )              :: coordinates
    type (enumerationComponentTypeType        ), intent(in   ), optional    :: componentType
    type (enumerationMassTypeType             ), intent(in   ), optional    :: massType

    density=self%massDistribution_%density(coordinates,componentType,massType)
    return
  end function sphericalTestHarnessDensity

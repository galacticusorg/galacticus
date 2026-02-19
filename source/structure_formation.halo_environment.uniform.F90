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
Implements a uniform halo environment.
!!}

  !![
  <haloEnvironment name="haloEnvironmentUniform">
   <description>Implements a uniform halo environment.</description>
  </haloEnvironment>
  !!]
  type, extends(haloEnvironmentClass) :: haloEnvironmentUniform
     !!{
     A uniform halo environment class.
     !!}
     private
   contains
     procedure :: overdensityLinear             => uniformOverdensityLinear
     procedure :: overdensityLinearGradientTime => uniformOverdensityLinearGradientTime
     procedure :: overdensityNonLinear          => uniformOverdensityNonLinear
     procedure :: environmentRadius             => uniformEnvironmentRadius
     procedure :: environmentMass               => uniformEnvironmentMass
     procedure :: pdf                           => uniformPDF
     procedure :: cdf                           => uniformCDF
     procedure :: overdensityLinearSet          => uniformOverdensityLinearSet
     procedure :: overdensityIsSettable         => uniformOverdensityIsSettable
     procedure :: isNodeDependent               => uniformIsNodeDependent
     procedure :: isTreeDependent               => uniformIsTreeDependent
  end type haloEnvironmentUniform

  interface haloEnvironmentUniform
     !!{
     Constructors for the \refClass{haloEnvironmentUniform} halo environment class.
     !!}
     module procedure uniformConstructorParameters
  end interface haloEnvironmentUniform

contains

  function uniformConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloEnvironmentUniform} halo environment class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(haloEnvironmentUniform)                :: self
    type(inputParameters       ), intent(inout) :: parameters

    self=haloEnvironmentUniform()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function uniformConstructorParameters

  double precision function uniformOverdensityLinear(self,node,presentDay)
    !!{
    Return the environment of the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class  (haloEnvironmentUniform), intent(inout)           :: self
    type   (treeNode              ), intent(inout)           :: node
    logical                        , intent(in   ), optional :: presentDay
    !$GLC attributes unused :: self, node, presentDay

    uniformOverdensityLinear=0.0d0
    return
  end function uniformOverdensityLinear

  double precision function uniformOverdensityLinearGradientTime(self,node)
    !!{
    Return the time gradient of the environment of the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(haloEnvironmentUniform), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    uniformOverdensityLinearGradientTime=0.0d0
    return
  end function uniformOverdensityLinearGradientTime

  double precision function uniformOverdensityNonLinear(self,node)
    !!{
    Return the environment of the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(haloEnvironmentUniform), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    uniformOverdensityNonLinear=0.0d0
    return
  end function uniformOverdensityNonLinear

  double precision function uniformEnvironmentRadius(self)
    !!{
    Return the radius of the environment.
    !!}
    implicit none
    class(haloEnvironmentUniform), intent(inout) :: self
    !$GLC attributes unused :: self

    uniformEnvironmentRadius=huge(0.0d0)
    return
  end function uniformEnvironmentRadius

  double precision function uniformEnvironmentMass(self)
    !!{
    Return the mass of the environment.
    !!}
    implicit none
    class(haloEnvironmentUniform), intent(inout) :: self
    !$GLC attributes unused :: self

    uniformEnvironmentMass=huge(0.0d0)
    return
  end function uniformEnvironmentMass

  double precision function uniformPDF(self,overdensity)
    !!{
    Return the PDF of the environmental overdensity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (haloEnvironmentUniform), intent(inout) :: self
    double precision                        , intent(in   ) :: overdensity
    !$GLC attributes unused :: self, overdensity

    uniformPDF=0.0d0
    call Error_Report('PDF is a delta function'//{introspection:location})
    return
  end function uniformPDF

  double precision function uniformCDF(self,overdensity)
    !!{
    Return the CDF of the environmental overdensity.
    !!}
    implicit none
    class           (haloEnvironmentUniform), intent(inout) :: self
    double precision                        , intent(in   ) :: overdensity
    !$GLC attributes unused :: self

    if (overdensity >= 0.0d0) then
       uniformCDF=1.0d0
    else
       uniformCDF=0.0d0
    end if
    return
  end function uniformCDF

  subroutine uniformOverdensityLinearSet(self,node,overdensity)
    !!{
    Return the CDF of the environmental overdensity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (haloEnvironmentUniform), intent(inout) :: self
    type            (treeNode              ), intent(inout) :: node
    double precision                        , intent(in   ) :: overdensity
    !$GLC attributes unused :: self, node

    if (overdensity /= 0.0d0) call Error_Report('non-zero overdensity is inconsistent with uniform density field'//{introspection:location})
    return
  end subroutine uniformOverdensityLinearSet

  logical function uniformOverdensityIsSettable(self)
    !!{
    Return false as the overdensity is not settable.
    !!}
    implicit none
    class(haloEnvironmentUniform), intent(inout) :: self
    !$GLC attributes unused :: self

    uniformOverdensityIsSettable=.false.
    return
  end function uniformOverdensityIsSettable

  logical function uniformIsNodeDependent(self)
    !!{
    Return false as the environment is not dependent on the node.
    !!}
    implicit none
    class(haloEnvironmentUniform), intent(inout) :: self
    !$GLC attributes unused :: self

    uniformIsNodeDependent=.false.
    return
  end function uniformIsNodeDependent

  logical function uniformIsTreeDependent(self)
    !!{
    Return false as the environment is not dependent on the tree.
    !!}
    implicit none
    class(haloEnvironmentUniform), intent(inout) :: self
    !$GLC attributes unused :: self

    uniformIsTreeDependent=.false.
    return
  end function uniformIsTreeDependent

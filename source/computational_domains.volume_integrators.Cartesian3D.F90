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

  !![
  <computationalDomainVolumeIntegrator name="computationalDomainVolumeIntegratorCartesian3D">
   <description>A computational domain volume integrator for 3D Cartesian cells.</description>
  </computationalDomainVolumeIntegrator>
  !!]
  type, extends(computationalDomainVolumeIntegratorClass) :: computationalDomainVolumeIntegratorCartesian3D
     !!{
     Implementation of a computational domain for 3D Cartesian cells.
     !!}
     private
     double precision, dimension(  2) :: xBoundaries, yBoundaries, &
          &                              zBoundaries
     double precision, dimension(3,2) :: boundaries
     double precision                 :: volume_
   contains
     procedure :: volume    => cartesian3DVolume
     procedure :: integrate => cartesian3DIntegrate
  end type computationalDomainVolumeIntegratorCartesian3D

  interface computationalDomainVolumeIntegratorCartesian3D
     !!{
     Constructors for the \refClass{computationalDomainVolumeIntegratorCartesian3D} computational domain.
     !!}
     module procedure cartesian3DConstructorParameters
     module procedure cartesian3DConstructorInternal
  end interface computationalDomainVolumeIntegratorCartesian3D

contains

  function cartesian3DConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{computationalDomainVolumeIntegratorCartesian3D} computational domain volume integrator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (computationalDomainVolumeIntegratorCartesian3D)                 :: self
    type            (inputParameters                               ), intent(inout)  :: parameters
    double precision                                                , dimension(  2) :: xBoundaries             , yBoundaries, &
         &                                                                              zBoundaries
    double precision                                                , dimension(3,2) :: boundaries

    !![
    <inputParameter>
      <name>xBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>The $x$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>yBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>The $y$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>zBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>The $z$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    boundaries(1,:)=xBoundaries
    boundaries(2,:)=yBoundaries
    boundaries(3,:)=zBoundaries
    self=computationalDomainVolumeIntegratorCartesian3D(boundaries)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cartesian3DConstructorParameters

  function cartesian3DConstructorInternal(boundaries) result(self)
    !!{
    Internal constructor for the \refClass{computationalDomainVolumeIntegratorCartesian3D} computational domain volume integrator class.
    !!}
    implicit none
    type            (computationalDomainVolumeIntegratorCartesian3D)                                :: self
    double precision                                                , dimension(3,2), intent(in   ) :: boundaries
    !![
    <constructorAssign variables="boundaries"/>
    !!]

    self%xBoundaries=boundaries(1,:)
    self%yBoundaries=boundaries(2,:)
    self%zBoundaries=boundaries(3,:)
    self%volume_    =+(boundaries(1,2)-boundaries(1,1)) &
         &           *(boundaries(2,2)-boundaries(2,1)) &
         &           *(boundaries(3,2)-boundaries(3,1)) 
    return
  end function cartesian3DConstructorInternal

  double precision function cartesian3DVolume(self)
    !!{
    Return the volume of the computational domain cell.
    !!}
    implicit none
    class(computationalDomainVolumeIntegratorCartesian3D), intent(inout) :: self

    cartesian3DVolume=self%volume_
    return
  end function cartesian3DVolume

  double precision function cartesian3DIntegrate(self,integrand)
    !!{
    Integrate over the computational domain cell.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Coordinates          , only : coordinateCartesian
    implicit none
    class    (computationalDomainVolumeIntegratorCartesian3D), intent(inout), target :: self
    procedure(computationalDomainVolumeIntegrand            )                        :: integrand
    type     (integrator                                    )                        :: integrator_
    type     (coordinateCartesian                           )                        :: coordinates

    integrator_         = integrator           (                                         &
         &                                                        cartesian3DIntegrandX, &
         &                                      toleranceRelative=1.0d-2                 &
         &                                     )
    cartesian3DIntegrate=+integrator_%integrate(                                         &
         &                                                         self%boundaries(1,1), &
         &                                                         self%boundaries(1,2)  &
         &                                     )  
    return

  contains

    double precision function cartesian3DIntegrandX(x)
      !!{
      $x$-integrand over Cartesian 3D computational domain cells.
      !!}
      implicit none
      double precision            , intent(in   ) :: x
      type            (integrator)                :: integrator_

      call coordinates%xSet(x)
      integrator_         = integrator            (                                         &
           &                                                         cartesian3DIntegrandY, &
           &                                       toleranceRelative=1.0d-2                 &
           &                                      )
      cartesian3DIntegrandX=+integrator_%integrate(                                         &
           &                                                         self%boundaries(2,1) , &
           &                                                         self%boundaries(2,2)   &
           &                                      )  
      return
    end function cartesian3DIntegrandX

    double precision function cartesian3DIntegrandY(y)
      !!{
      $y$-integrand over Cartesian 3D computational domain cells.
      !!}
      implicit none
      double precision            , intent(in   ) :: y
      type            (integrator)                :: integrator_

      call coordinates%ySet(y)
      integrator_         = integrator            (                                         &
           &                                                         cartesian3DIntegrandZ, &
           &                                       toleranceRelative=1.0d-2                 &
           &                                      )
      cartesian3DIntegrandY=+integrator_%integrate(                                         &
           &                                                         self%boundaries(3,1) , &
           &                                                         self%boundaries(3,2)   &
           &                                      )  
      return
    end function cartesian3DIntegrandY

    double precision function cartesian3DIntegrandZ(z)
      !!{
      $z$-integrand over Cartesian 3D computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: z

      call coordinates%zSet(z)
      cartesian3DIntegrandZ=integrand(coordinates)
      return
    end function cartesian3DIntegrandZ

  end function cartesian3DIntegrate

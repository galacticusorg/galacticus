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
  <computationalDomainVolumeIntegrator name="computationalDomainVolumeIntegratorCylindrical">
   <description>A computational domain volume integrator for cylindrical cells.</description>
  </computationalDomainVolumeIntegrator>
  !!]
  type, extends(computationalDomainVolumeIntegratorClass) :: computationalDomainVolumeIntegratorCylindrical
     !!{
     Implementation of a computational domain for cylindrical cells.
     !!}
     private
     double precision, dimension(  2) :: rBoundaries, zBoundaries
     double precision, dimension(2,2) :: boundaries
     double precision                 :: volume_
   contains
     procedure :: volume    => cylindricalVolume
     procedure :: integrate => cylindricalIntegrate
  end type computationalDomainVolumeIntegratorCylindrical

  interface computationalDomainVolumeIntegratorCylindrical
     !!{
     Constructors for the \refClass{computationalDomainVolumeIntegratorCylindrical} computational domain.
     !!}
     module procedure cylindricalConstructorParameters
     module procedure cylindricalConstructorInternal
  end interface computationalDomainVolumeIntegratorCylindrical

contains

  function cylindricalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{computationalDomainVolumeIntegratorCylindrical} computational domain volume integrator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (computationalDomainVolumeIntegratorCylindrical)                 :: self
    type            (inputParameters                               ), intent(inout)  :: parameters
    double precision                                                , dimension(  2) :: rBoundaries, zBoundaries
    double precision                                                , dimension(2,2) :: boundaries

    !![
    <inputParameter>
      <name>rBoundaries</name>
      <defaultValue>[0.0d0,1.0d0]</defaultValue>
      <description>The $r$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>zBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>The $z$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    boundaries(1,:)=rBoundaries
    boundaries(2,:)=zBoundaries
    self=computationalDomainVolumeIntegratorCylindrical(boundaries)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cylindricalConstructorParameters

  function cylindricalConstructorInternal(boundaries) result(self)
    !!{
    Internal constructor for the \refClass{computationalDomainVolumeIntegratorCylindrical} computational domain volume integrator class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (computationalDomainVolumeIntegratorCylindrical)                              :: self
    double precision                                              , dimension(2,2), intent(in   ) :: boundaries
    !![
    <constructorAssign variables="boundaries"/>
    !!]

    self%rBoundaries=boundaries(1,:)
    self%zBoundaries=boundaries(2,:)
    self%volume_    =+2.0d0                                   &
         &           *Pi                                      &
         &           *(boundaries(1,2)**2-boundaries(1,1)**2) &
         &           *(boundaries(2,2)   -boundaries(2,1)   )
    return
  end function cylindricalConstructorInternal

  double precision function cylindricalVolume(self)
    !!{
    Return the volume of the computational domain cell.
    !!}
    implicit none
    class(computationalDomainVolumeIntegratorCylindrical), intent(inout) :: self

    cylindricalVolume=self%volume_
    return
  end function cylindricalVolume

  double precision function cylindricalIntegrate(self,integrand)
    !!{
    Integrate over the computational domain cell.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Coordinates          , only : coordinateCylindrical
    implicit none
    class    (computationalDomainVolumeIntegratorCylindrical), intent(inout), target :: self
    procedure(computationalDomainVolumeIntegrand            )                        :: integrand
    type     (integrator                                    )                        :: integrator_
    type     (coordinateCylindrical                         )                        :: coordinates

    integrator_         =integrator           (                                         &
         &                                                       cylindricalIntegrandR, &
         &                                     toleranceRelative=1.0d-2                 &
         &                                    )
    cylindricalIntegrate=integrator_%integrate(                                         &
         &                                                       self%boundaries(1,1) , &
         &                                                       self%boundaries(1,2)   &
         &)                                    
    return

  contains
    
    double precision function cylindricalIntegrandR(r)
      !!{
      $r$-integrand over cylindrical computational domain cells.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision            , intent(in   ) :: r
      type            (integrator)                :: integrator_

      call coordinates%rSet(r)
      integrator_         =  integrator           (cylindricalIntegrandPhi,toleranceRelative=1.0d-2   )
      cylindricalIntegrandR=+integrator_%integrate(0.0d0                  ,                  2.0d+0*Pi) &
           &                *r
      return
    end function cylindricalIntegrandR

    double precision function cylindricalIntegrandPhi(phi)
      !!{
      $\phi$-integrand over cylindrical computational domain cells.
      !!}
      implicit none
      double precision            , intent(in   ) :: phi
      type            (integrator)                :: integrator_

      call coordinates%phiSet(phi)
      integrator_            = integrator           (                                         &
           &                                                           cylindricalIntegrandZ, &
           &                                         toleranceRelative=1.0d-2                 &
           &                                        )
      cylindricalIntegrandPhi=+integrator_%integrate(                                         &
           &                                                           self%boundaries(2,1) , &
           &                                                           self%boundaries(2,2)   &
           &                                        )
      return
    end function cylindricalIntegrandPhi

    double precision function cylindricalIntegrandZ(z)
      !!{
      $z$-integrand over cylindrical computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: z

      call coordinates%zSet(z)
      cylindricalIntegrandZ=integrand(coordinates)
      return
    end function cylindricalIntegrandZ

  end function cylindricalIntegrate

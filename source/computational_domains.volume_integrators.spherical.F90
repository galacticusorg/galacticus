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
  <computationalDomainVolumeIntegrator name="computationalDomainVolumeIntegratorSpherical">
   <description>A computational domain volume integrator for spherical cells.</description>
  </computationalDomainVolumeIntegrator>
  !!]
  type, extends(computationalDomainVolumeIntegratorClass) :: computationalDomainVolumeIntegratorSpherical
     !!{
     Implementation of a computational domain for spherical cells.
     !!}
     private
     double precision, dimension(2) :: boundaries
     double precision               :: volume_
   contains
     procedure :: volume    => sphericalVolume
     procedure :: integrate => sphericalIntegrate
  end type computationalDomainVolumeIntegratorSpherical

  interface computationalDomainVolumeIntegratorSpherical
     !!{
     Constructors for the \refClass{computationalDomainVolumeIntegratorSpherical} computational domain.
     !!}
     module procedure sphericalConstructorParameters
     module procedure sphericalConstructorInternal
  end interface computationalDomainVolumeIntegratorSpherical

contains

  function sphericalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{computationalDomainVolumeIntegratorSpherical} computational domain volume integrator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (computationalDomainVolumeIntegratorSpherical)                 :: self
    type            (inputParameters                             ), intent(inout)  :: parameters
    double precision                                              , dimension(  2) :: boundaries

    !![
    <inputParameter>
      <name>boundaries</name>
      <defaultValue>[0.0d0,1.0d0]</defaultValue>
      <description>The $r$-interval spanned by the computational domain.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
     self=computationalDomainVolumeIntegratorSpherical(boundaries)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sphericalConstructorParameters
  
  function sphericalConstructorInternal(boundaries) result(self)
    !!{
    Internal constructor for the \refClass{computationalDomainVolumeIntegratorSpherical} computational domain volume integrator class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (computationalDomainVolumeIntegratorSpherical)                              :: self
    double precision                                              , dimension(2), intent(in   ) :: boundaries
    !![
    <constructorAssign variables="boundaries"/>
    !!]

    self%volume_=+4.0d0              &
         &       *Pi                 &
         &       /3.0d0              &
         &       *(                  &
         &         +boundaries(2)**3 &
         &         -boundaries(1)**3 &
         &        )
    return
  end function sphericalConstructorInternal

  double precision function sphericalVolume(self)
    !!{
    Return the volume of the computational domain cell.
    !!}
    implicit none
    class(computationalDomainVolumeIntegratorSpherical), intent(inout) :: self

    sphericalVolume=self%volume_
    return
  end function sphericalVolume

  double precision function sphericalIntegrate(self,integrand)
    !!{
    Integrate over the computational domain cell.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Coordinates          , only : coordinateSpherical
    implicit none
    class    (computationalDomainVolumeIntegratorSpherical), intent(inout), target :: self
    procedure(computationalDomainVolumeIntegrand          )                        :: integrand
    type     (integrator                                  )                        :: integrator_
    type     (coordinateSpherical                         )                        :: coordinates

    integrator_       = integrator           (                                       &
         &                                                      sphericalIntegrandR, &
         &                                    toleranceRelative=1.0d-2               &
         &                                   )
    sphericalIntegrate=+integrator_%integrate(                                       &
         &                                                      self%boundaries(1) , &
         &                                                      self%boundaries(2)   &
         &                                   )  
    return

  contains
    
    double precision function sphericalIntegrandR(r)
      !!{
      $r$-integrand over spherical computational domain cells.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision            , intent(in   ) :: r
      type            (integrator)                :: integrator_

      call coordinates%rSet(r)
      integrator_        = integrator           (                                           &
         &                                                         sphericalIntegrandTheta, &
         &                                       toleranceRelative=1.0d-2                   &
         &                                      )
      sphericalIntegrandR=+integrator_%integrate(                                           &
           &                                                       0.0d+0                 , &
           &                                                       Pi                       &
           &                                    )                                           &
           &              *r**2
      return
    end function sphericalIntegrandR

    double precision function sphericalIntegrandTheta(theta)
      !!{
      $\theta$-integrand over spherical computational domain cells.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision            , intent(in   ) :: theta
      type            (integrator)                :: integrator_

      call coordinates%thetaSet(theta)
      integrator_            = integrator           (                                         &
         &                                                             sphericalIntegrandPhi, &
         &                                           toleranceRelative=1.0d-2                 &
         &                                          )
      sphericalIntegrandTheta=+integrator_%integrate(                                         &
           &                                                           0.0d+0               , &
           &                                                           2.0d+0*Pi              &
           &                                        )                                         &
           &                  *sin(theta)
      return
    end function sphericalIntegrandTheta

    double precision function sphericalIntegrandPhi(phi)
      !!{
      $\phi$-integrand over spherical computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: phi

      call coordinates%phiSet(phi)
      sphericalIntegrandPhi=integrand(coordinates)
      return
    end function sphericalIntegrandPhi

  end function sphericalIntegrate

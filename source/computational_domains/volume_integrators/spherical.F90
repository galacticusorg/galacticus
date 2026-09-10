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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !![
  <computationalDomainVolumeIntegrator name="computationalDomainVolumeIntegratorSpherical" docformat="rst">
   <description>
   Computes volume integrals over spherical shell cells, with the radial boundaries of the domain specified by ``[boundaries]``. Cell volumes are computed from the difference of the enclosed spherical volumes at the inner and outer radial boundaries of each shell.
   </description>
  </computationalDomainVolumeIntegrator>
  !!]
  type, extends(computationalDomainVolumeIntegratorClass) :: computationalDomainVolumeIntegratorSpherical
     !!{RST
     Implementation of a computational domain for spherical cells.
     !!}
     private
     double precision, dimension(2) :: boundaries
     double precision               :: volume_          , toleranceRelative, &
          &                            toleranceAbsolute
   contains
     procedure :: volume       => sphericalVolume
     procedure :: toleranceSet => sphericalToleranceSet
     procedure :: integrate    => sphericalIntegrate
  end type computationalDomainVolumeIntegratorSpherical

  interface computationalDomainVolumeIntegratorSpherical
     !!{RST
     Constructors for the :galacticus-class:`computationalDomainVolumeIntegratorSpherical` computational domain.
     !!}
     module procedure sphericalConstructorParameters
     module procedure sphericalConstructorInternal
  end interface computationalDomainVolumeIntegratorSpherical

contains

  function sphericalConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`computationalDomainVolumeIntegratorSpherical` computational domain volume integrator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (computationalDomainVolumeIntegratorSpherical)                 :: self
    type            (inputParameters                             ), intent(inout)  :: parameters
    double precision                                              , dimension(  2) :: boundaries

    !![
    <inputParameter docformat="rst">
      <name>boundaries</name>
      <defaultValue>[0.0d0,1.0d0]</defaultValue>
      <description>
      A two-element array :math:`[r_\mathrm{min}, r_\mathrm{max}]` specifying the radial extent of the spherically symmetric integration domain.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
     self=computationalDomainVolumeIntegratorSpherical(boundaries)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function sphericalConstructorParameters
  
  function sphericalConstructorInternal(boundaries,toleranceAbsolute,toleranceRelative) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`computationalDomainVolumeIntegratorSpherical` computational domain volume integrator class. The optional toleranceRelative and toleranceAbsolute arguments specify the tolerances to which the volume integral is to be evaluated---the absolute tolerance is apportioned between the nested one-dimensional integrals in proportion to the measure which each is integrated against by the level above it.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (computationalDomainVolumeIntegratorSpherical)                              :: self
    double precision                                              , dimension(2), intent(in   ) :: boundaries
    double precision                                              , optional    , intent(in   ) :: toleranceAbsolute, toleranceRelative
    !![
    <constructorAssign variables="boundaries"/>
    <optionalArgument name="toleranceAbsolute" defaultsTo="0.0d0"                          />
    <optionalArgument name="toleranceRelative" defaultsTo="toleranceRelativeVolumeIntegral"/>
    !!]

    self%toleranceAbsolute=toleranceAbsolute_
    self%toleranceRelative=toleranceRelative_
    self%volume_          =+4.0d0              &
         &                 *Pi                 &
         &                 /3.0d0              &
         &                 *(                  &
         &                   +boundaries(2)**3 &
         &                   -boundaries(1)**3 &
         &                  )
    return
  end function sphericalConstructorInternal

  double precision function sphericalVolume(self)
    !!{RST
    Return the volume of the computational domain cell.
    !!}
    implicit none
    class(computationalDomainVolumeIntegratorSpherical), intent(inout) :: self

    sphericalVolume=self%volume_
    return
  end function sphericalVolume

  subroutine sphericalToleranceSet(self,toleranceAbsolute,toleranceRelative)
    !!{RST
    Set the tolerances to which the volume integral over the computational domain cell is to be evaluated.
    !!}
    implicit none
    class           (computationalDomainVolumeIntegratorSpherical), intent(inout) :: self
    double precision                                              , intent(in   ) :: toleranceAbsolute, toleranceRelative

    self%toleranceAbsolute=toleranceAbsolute
    self%toleranceRelative=toleranceRelative
    return
  end subroutine sphericalToleranceSet

  double precision function sphericalIntegrate(self,integrand)
    !!{RST
    Integrate over the computational domain cell.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator         , GSL_Integ_Gauss15
    use :: Coordinates             , only : coordinateSpherical
    implicit none
    class           (computationalDomainVolumeIntegratorSpherical), intent(inout), target :: self
    procedure       (computationalDomainVolumeIntegrand          )                        :: integrand
    type            (integrator                                  )                        :: integratorR         , integratorTheta      , &
         &                                                                                   integratorPhi
    type            (coordinateSpherical                         )                        :: coordinates
    double precision                                                                      :: toleranceAbsoluteR  , toleranceAbsoluteTheta, &
         &                                                                                   toleranceAbsolutePhi

    ! Apportion the absolute tolerance between the nested integrals. The radial integral is required to meet the requested absolute
    ! tolerance directly. Each inner integral supplies the integrand of the integral one level out, and so need only be accurate to
    ! the tolerance of that outer integral divided by the measure, ∫r²dr and ∫sin(θ)dθ respectively, against which it is there
    ! integrated.
    toleranceAbsoluteR    =+self%toleranceAbsolute
    toleranceAbsoluteTheta=+toleranceAbsoluteR         &
         &                 /(                          &
         &                   +self%boundaries(2)**3    &
         &                   -self%boundaries(1)**3    &
         &                  )                          &
         &                 *3.0d0
    toleranceAbsolutePhi  =+toleranceAbsoluteTheta     &
         &                 /2.0d0
    ! Construct the integrators. These are built once here, rather than in the integrands where they would be rebuilt (along with
    ! their GSL workspaces) on every evaluation.
    integratorR       =integrator(sphericalIntegrandR    ,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsoluteR    ,integrationRule=GSL_Integ_Gauss15)
    integratorTheta   =integrator(sphericalIntegrandTheta,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsoluteTheta,integrationRule=GSL_Integ_Gauss15)
    integratorPhi     =integrator(sphericalIntegrandPhi  ,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsolutePhi  ,integrationRule=GSL_Integ_Gauss15)
    sphericalIntegrate=+integratorR%integrate(                    &
         &                                    self%boundaries(1), &
         &                                    self%boundaries(2)  &
         &                                   )
    return

  contains
    
    double precision function sphericalIntegrandR(r)
      !!{RST
      :math:`r`-integrand over spherical computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: r

      call coordinates%rSet(r)
      sphericalIntegrandR=+integratorTheta%integrate(        &
           &                                         0.0d+0, &
           &                                         Pi      &
           &                                        )        &
           &              *r**2
      return
    end function sphericalIntegrandR

    double precision function sphericalIntegrandTheta(theta)
      !!{RST
      :math:`\theta`-integrand over spherical computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: theta

      call coordinates%thetaSet(theta)
      sphericalIntegrandTheta=+integratorPhi%integrate(           &
           &                                           0.0d+0   , &
           &                                           2.0d+0*Pi  &
           &                                          )           &
           &                  *sin(theta)
      return
    end function sphericalIntegrandTheta

    double precision function sphericalIntegrandPhi(phi)
      !!{RST
      :math:`\phi`-integrand over spherical computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: phi

      call coordinates%phiSet(phi)
      sphericalIntegrandPhi=integrand(coordinates)
      return
    end function sphericalIntegrandPhi

  end function sphericalIntegrate

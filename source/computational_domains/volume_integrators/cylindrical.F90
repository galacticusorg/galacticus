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
  <computationalDomainVolumeIntegrator name="computationalDomainVolumeIntegratorCylindrical" docformat="rst">
   <description>
   Computes volume integrals over cylindrical grid cells, accounting for the annular geometry of each cell. The radial and vertical extents of the domain are specified by ``[rBoundaries]`` and ``[zBoundaries]``, with cell volumes computed from the cylindrical shell geometry.
   </description>
  </computationalDomainVolumeIntegrator>
  !!]
  type, extends(computationalDomainVolumeIntegratorClass) :: computationalDomainVolumeIntegratorCylindrical
     !!{RST
     Implementation of a computational domain for cylindrical cells.
     !!}
     private
     double precision, dimension(  2) :: rBoundaries      , zBoundaries
     double precision, dimension(2,2) :: boundaries
     double precision                 :: volume_          , toleranceRelative, &
          &                              toleranceAbsolute
   contains
     procedure :: volume       => cylindricalVolume
     procedure :: toleranceSet => cylindricalToleranceSet
     procedure :: integrate    => cylindricalIntegrate
  end type computationalDomainVolumeIntegratorCylindrical

  interface computationalDomainVolumeIntegratorCylindrical
     !!{RST
     Constructors for the :galacticus-class:`computationalDomainVolumeIntegratorCylindrical` computational domain.
     !!}
     module procedure cylindricalConstructorParameters
     module procedure cylindricalConstructorInternal
  end interface computationalDomainVolumeIntegratorCylindrical

contains

  function cylindricalConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`computationalDomainVolumeIntegratorCylindrical` computational domain volume integrator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (computationalDomainVolumeIntegratorCylindrical)                 :: self
    type            (inputParameters                               ), intent(inout)  :: parameters
    double precision                                                , dimension(  2) :: rBoundaries, zBoundaries
    double precision                                                , dimension(2,2) :: boundaries

    !![
    <inputParameter docformat="rst">
      <name>rBoundaries</name>
      <defaultValue>[0.0d0,1.0d0]</defaultValue>
      <description>
      A two-element array :math:`[r_\mathrm{min}, r_\mathrm{max}]` specifying the radial extent of the cylindrical integration domain.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>zBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>
      A two-element array :math:`[z_\mathrm{min}, z_\mathrm{max}]` specifying the vertical extent of the cylindrical integration domain along the symmetry axis.
      </description>
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

  function cylindricalConstructorInternal(boundaries,toleranceAbsolute,toleranceRelative) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`computationalDomainVolumeIntegratorCylindrical` computational domain volume integrator class. The optional  toleranceAbsolute and  toleranceRelative arguments specify the tolerances to which the volume integral is to be evaluated---the absolute tolerance is apportioned between the nested one-dimensional integrals in proportion to the measure which each is integrated against by the level above it.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (computationalDomainVolumeIntegratorCylindrical)                                :: self
    double precision                                                , dimension(2,2), intent(in   ) :: boundaries
    double precision                                                , optional      , intent(in   ) :: toleranceAbsolute, toleranceRelative
    !![
    <constructorAssign variables="boundaries"/>
    <optionalArgument name="toleranceAbsolute" defaultsTo="0.0d0"                          />
    <optionalArgument name="toleranceRelative" defaultsTo="toleranceRelativeVolumeIntegral"/>
    !!]

    self%toleranceAbsolute=toleranceAbsolute_
    self%toleranceRelative=toleranceRelative_
    self%rBoundaries      =boundaries(1,:)
    self%zBoundaries      =boundaries(2,:)
    ! The cell is an annulus, spanning the full 2π in azimuth. Its volume is therefore ∫r dr ∫dφ ∫dz = π (r₂²-r₁²) (z₂-z₁).
    self%volume_    =+Pi                                      &
         &           *(boundaries(1,2)**2-boundaries(1,1)**2) &
         &           *(boundaries(2,2)   -boundaries(2,1)   )
    return
  end function cylindricalConstructorInternal

  double precision function cylindricalVolume(self)
    !!{RST
    Return the volume of the computational domain cell.
    !!}
    implicit none
    class(computationalDomainVolumeIntegratorCylindrical), intent(inout) :: self

    cylindricalVolume=self%volume_
    return
  end function cylindricalVolume

  subroutine cylindricalToleranceSet(self,toleranceAbsolute,toleranceRelative)
    !!{RST
    Set the tolerances to which the volume integral over the computational domain cell is to be evaluated.
    !!}
    implicit none
    class           (computationalDomainVolumeIntegratorCylindrical), intent(inout) :: self
    double precision                                                , intent(in   ) :: toleranceAbsolute, toleranceRelative

    self%toleranceAbsolute=toleranceAbsolute
    self%toleranceRelative=toleranceRelative
    return
  end subroutine cylindricalToleranceSet

  double precision function cylindricalIntegrate(self,integrand)
    !!{RST
    Integrate over the computational domain cell.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator           , GSL_Integ_Gauss15
    use :: Coordinates             , only : coordinateCylindrical
    implicit none
    class           (computationalDomainVolumeIntegratorCylindrical), intent(inout), target :: self
    procedure       (computationalDomainVolumeIntegrand            )                        :: integrand
    type            (integrator                                    )                        :: integratorR       , integratorPhi       , &
         &                                                                                     integratorZ
    type            (coordinateCylindrical                         )                        :: coordinates
    double precision                                                                        :: toleranceAbsoluteR, toleranceAbsolutePhi, &
         &                                                                                     toleranceAbsoluteZ

    ! Apportion the absolute tolerance between the nested integrals. The radial integral is required to meet the requested absolute
    ! tolerance directly. Each inner integral supplies the integrand of the integral one level out, and so need only be accurate to
    ! the tolerance of that outer integral divided by the measure, ∫rdr and ∫dφ respectively, against which it is there integrated.
    toleranceAbsoluteR   =+self%toleranceAbsolute
    toleranceAbsolutePhi =+toleranceAbsoluteR           &
         &                /(                            &
         &                  +self%boundaries(1,2)**2    &
         &                  -self%boundaries(1,1)**2    &
         &                 )                            &
         &                *2.0d0
    toleranceAbsoluteZ   =+toleranceAbsolutePhi         &
         &                /2.0d0                        &
         &                /Pi
    ! Construct the integrators. These are built once here, rather than in the integrands where they would be rebuilt (along with
    ! their GSL workspaces) on every evaluation.
    integratorR          =integrator(cylindricalIntegrandR  ,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsoluteR  ,integrationRule=GSL_Integ_Gauss15)
    integratorPhi        =integrator(cylindricalIntegrandPhi,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsolutePhi,integrationRule=GSL_Integ_Gauss15)
    integratorZ          =integrator(cylindricalIntegrandZ  ,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsoluteZ  ,integrationRule=GSL_Integ_Gauss15)
    cylindricalIntegrate =+integratorR%integrate(                         &
         &                                       self%boundaries(1,1)   , &
         &                                       self%boundaries(1,2)     &
         &                                      )
    return

  contains
    
    double precision function cylindricalIntegrandR(r)
      !!{RST
      :math:`r`-integrand over cylindrical computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: r

      call coordinates%rSet(r)
      cylindricalIntegrandR=+integratorPhi%integrate(          &
           &                                         0.0d0   , &
           &                                         2.0d0*Pi  &
           &                                        )          &
           &                *r
      return
    end function cylindricalIntegrandR

    double precision function cylindricalIntegrandPhi(phi)
      !!{RST
      :math:`\phi`-integrand over cylindrical computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: phi

      call coordinates%phiSet(phi)
      cylindricalIntegrandPhi=+integratorZ%integrate(                      &
           &                                         self%boundaries(2,1), &
           &                                         self%boundaries(2,2)  &
           &                                        )
      return
    end function cylindricalIntegrandPhi

    double precision function cylindricalIntegrandZ(z)
      !!{RST
      :math:`z`-integrand over cylindrical computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: z

      call coordinates%zSet(z)
      cylindricalIntegrandZ=integrand(coordinates)
      return
    end function cylindricalIntegrandZ

  end function cylindricalIntegrate

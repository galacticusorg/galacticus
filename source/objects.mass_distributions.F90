!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the mass distribution class.

module Mass_Distributions
  !% Implements the mass distribution class.
  use FGSL
  implicit none
  private
  public :: Mass_Distribution_Create

  type, public                            :: massDistribution
     !% The basic mass distribution class. Has no symmetry and will abort on inqueries.
     logical :: dimensionless
   contains
     !@ <objectMethods>
     !@   <object>massDistribution</object>
     !@   <objectMethod>
     !@     <method>symmetry</method>
     !@     <description>Returns a label specifying the symmetry of the mass distribution (see \S\ref{sec:MassDistributions}).</description>
     !@     <type>\enumMassDistributionSymmetry</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isDimensionless</method>
     !@     <description>Returns {\tt true} is this is a dimensionless mass distribution, {\tt false} otherwise.</description>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>density</method>
     !@     <description>Returns the density of the mass distribution at the supplied {\tt coordinates}.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(coordinate)\textgreater} coordinates\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>densityRadialMoment</method>
     !@     <description>Returns the $n^{\rm th}$ moment of the integral of the density over radius, $\int_0^\infty \rho({\bf x}) |x|^n {\rm d} {\bf x}$.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ moment\argin, \logicalzero\ [isInfinite]\argout</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>massEnclosedBySphere</method>
     !@     <description>Returns the mass enclosed by a sphere of given {\tt radius} centered on the origin.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ radius\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>potential</method>
     !@     <description>Returns the gravitational potential at the specified {\tt coordinates}.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(coordinate)\textgreater} coordinates\argin</arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure, nopass :: symmetry            =>Mass_Distribution_Symmetry_None
     procedure         :: isDimensionless     =>Mass_Distribution_Is_Dimensionless
     procedure         :: density             =>Mass_Distribution_Density_Null
     procedure         :: densityRadialMoment =>Mass_Distribution_Density_Radial_Moment_Null
     procedure         :: massEnclosedBySphere=>Mass_Distribution_Mass_Enc_By_Sphere_Null
     procedure         :: potential           =>Mass_Distribution_Potential_Null
  end type massDistribution

  type, public, extends(massDistribution) :: massDistributionSpherical
     !% A spherical mass distribution class.
   contains
     !@ <objectMethods>
     !@   <object>massDistributionSpherical</object>
     !@   <objectMethod>
     !@     <method>halfMassRadius</method>
     !@     <description>Returns the radius enclosing half of the mass of the mass distribution.</description>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure, nopass :: symmetry            =>Mass_Distribution_Symmetry_Spherical
     procedure         :: massEnclosedBySphere=>Mass_Distribution_Mass_Enc_By_Sphere_Spherical
     procedure         :: halfMassRadius      =>Mass_Distribution_Half_Mass_Radius_Spherical
  end type massDistributionSpherical
  class(massDistributionSpherical), pointer :: massDistributionSphericalActive
  !$omp threadprivate(massDistributionSphericalActive)
  type, public, extends(massDistribution) :: massDistributionCylindrical
     !% A cylindrical mass distribution class.
   contains
     procedure, nopass :: symmetry=>Mass_Distribution_Symmetry_Cylindrical
  end type massDistributionCylindrical

  ! Include type definitions.
  include 'objects.mass_distributions.NFW.type.inc'
  include 'objects.mass_distributions.beta_profile.type.inc'
  include 'objects.mass_distributions.Hernquist.type.inc'
  include 'objects.mass_distributions.Sersic.type.inc'

  ! Labels for mass distribution symmetries.
  !@ <enumeration>
  !@  <name>massDistributionSymmetry</name>
  !@  <description>Used to specify the symmetry of {\tt massDistribution} objects.</description>
  !@  <entry label="massDistributionSymmetryNone"        />
  !@  <entry label="massDistributionSymmetryCylindrical" />
  !@  <entry label="massDistributionSymmetrySpherical"   />
  !@ </enumeration>
  integer                             , parameter, public :: massDistributionSymmetryNone       =0
  integer                             , parameter, public :: massDistributionSymmetryCylindrical=1
  integer                             , parameter, public :: massDistributionSymmetrySpherical  =2

  ! Template distributions.
  type   (massDistributionNFW        )                    :: massDistributionTemplateNFW
  type   (massDistributionBetaProfile)                    :: massDistributionTemplateBetaProfile
  type   (massDistributionHernquist  )                    :: massDistributionTemplateHernquist
  type   (massDistributionSersic     )                    :: massDistributionTemplateSersic

contains

  function Mass_Distribution_Create(type) result (newMassDistribution)
    !% Create a mass distribution given the name.
    use Galacticus_Error
    implicit none
    class    (massDistribution), pointer       :: newMassDistribution
    character(len=*           ), intent(in   ) :: type

    if      (trim(type) == "NFW"        ) then
       allocate(newMassDistribution,source=massDistributionTemplateNFW        )
    else if (trim(type) == "betaProfile") then
       allocate(newMassDistribution,source=massDistributionTemplateBetaProfile)
    else if (trim(type) == "hernquist") then
       allocate(newMassDistribution,source=massDistributionTemplateHernquist  )
    else if (trim(type) == "sersic") then
       allocate(newMassDistribution,source=massDistributionTemplateSersic     )
    else
       call Galacticus_Error_Report('Mass_Distribution_Create','unrecognized mass distribution')
    end if
    return
  end function Mass_Distribution_Create

  integer function Mass_Distribution_Symmetry_None()
    !% Returns symmetry label for mass dsitributions with no symmetry.
    implicit none

    Mass_Distribution_Symmetry_None=massDistributionSymmetryNone
    return
  end function Mass_Distribution_Symmetry_None

  integer function Mass_Distribution_Symmetry_Cylindrical()
    !% Returns symmetry label for mass dsitributions with cylindrical symmetry.
    implicit none

    Mass_Distribution_Symmetry_Cylindrical=massDistributionSymmetryCylindrical
    return
  end function Mass_Distribution_Symmetry_Cylindrical

  integer function Mass_Distribution_Symmetry_Spherical()
    !% Returns symmetry label for mass dsitributions with spherical symmetry.
    implicit none

    Mass_Distribution_Symmetry_Spherical=massDistributionSymmetrySpherical
    return
  end function Mass_Distribution_Symmetry_Spherical

  logical function Mass_Distribution_Is_Dimensionless(self)
    !% Return true if {\tt self} is a dimensionless mass distribution.
    implicit none
    class(massDistribution), intent(in   ) :: self

    Mass_Distribution_Is_Dimensionless=self%dimensionless
    return
  end function Mass_Distribution_Is_Dimensionless

  double precision function Mass_Distribution_Density_Null(self,coordinates)
    !% Aborts on attempts to get density of mass distributions with no density defined.
    use Coordinates
    use Galacticus_Error
    implicit none
    class(massDistribution), intent(in   ) :: self
    class(coordinate      ), intent(in   ) :: coordinates

    call Galacticus_Error_Report('Mass_Distribution_Density_Null','this mass distribution has no density method defined')
    return
  end function Mass_Distribution_Density_Null

  double precision function Mass_Distribution_Density_Radial_Moment_Null(self,moment,isInfinite)
    !% Aborts on attempts to get radial density moment of mass distributions with no density defined.
    use Galacticus_Error
    implicit none
    class           (massDistribution), intent(in   )           :: self
    double precision                  , intent(in   )           :: moment
    logical                           , intent(  out), optional :: isInfinite

    call Galacticus_Error_Report('Mass_Distribution_Density_Radial_Moment_Null','this mass distribution has no radial density moment method defined')
    return
  end function Mass_Distribution_Density_Radial_Moment_Null

  double precision function Mass_Distribution_Mass_Enc_By_Sphere_Null(self,radius)
    !% Aborts on attempts to get the mass enclosed by a sphere for mass distributions with no density defined.
    use Galacticus_Error
    implicit none
    class           (massDistribution), intent(in   ), target :: self
    double precision                  , intent(in   )         :: radius

    call Galacticus_Error_Report('Mass_Distribution_Mass_Enc_By_Sphere_Null','this mass distribution has no massEnclosedBySphere defined')
    return
  end function Mass_Distribution_Mass_Enc_By_Sphere_Null

  double precision function Mass_Distribution_Mass_Enc_By_Sphere_Spherical(self,radius)
    !% Computes the mass enclosed within a sphere of given {\tt radius} for spherically-symmetric mass distributions using
    !% numerical integration.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Numerical_Constants_Math
    implicit none
    class           (massDistributionSpherical ), intent(in   ), target :: self
    double precision                            , intent(in   )         :: radius
    type            (c_ptr                     )                        :: parameterPointer
    type            (fgsl_function             )                        :: integrandFunction
    type            (fgsl_integration_workspace)                        :: integrationWorkspace
    logical                                                             :: integrationReset

    massDistributionSphericalActive => self
    integrationReset=.true.
    Mass_Distribution_Mass_Enc_By_Sphere_Spherical=4.0d0*Pi*Integrate(0.0d0,radius&
         &,Mass_Distribution_Mass_Enc_By_Sphere_Spherical_Integrand,parameterPointer,integrandFunction,integrationWorkspace,reset&
         &=integrationReset,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function Mass_Distribution_Mass_Enc_By_Sphere_Spherical

  function Mass_Distribution_Mass_Enc_By_Sphere_Spherical_Integrand(radius,parameterPointer) bind(c)
    !% Enclosed mass integrand for spherical mass distributions.
    use, intrinsic :: ISO_C_Binding
    use Coordinates
    implicit none
    real(kind=c_double      )        :: Mass_Distribution_Mass_Enc_By_Sphere_Spherical_Integrand
    real(kind=c_double      ), value :: radius
    type(c_ptr              ), value :: parameterPointer
    type(coordinateSpherical)        :: position

    position=[radius,0.0d0,0.0d0]
    Mass_Distribution_Mass_Enc_By_Sphere_Spherical_Integrand=radius**2*massDistributionSphericalActive%density(position)
    return
  end function Mass_Distribution_Mass_Enc_By_Sphere_Spherical_Integrand

  double precision function Mass_Distribution_Potential_Null(self,coordinates)
    !% Aborts on attempts to get potential of mass distributions with no density defined.
    use Coordinates
    use Galacticus_Error
    implicit none
    class(massDistribution), intent(in   ) :: self
    class(coordinate      ), intent(in   ) :: coordinates

    call Galacticus_Error_Report('Mass_Distribution_Potential_Null','this mass distribution has no potential method defined')
    return
  end function Mass_Distribution_Potential_Null

  double precision function Mass_Distribution_Half_Mass_Radius_Spherical(self)
    !% Aborts on attempts to get half-mass radius in spherical mass distributions.
    use Galacticus_Error
    implicit none
    class(massDistributionSpherical), intent(in   ) :: self

    call Galacticus_Error_Report('Mass_Distribution_Half_Mass_Radius_Spherical','this mass distribution has no halfMassRadius method defined')
    return
  end function Mass_Distribution_Half_Mass_Radius_Spherical

  ! Include mass distribution methods.
  include 'objects.mass_distributions.NFW.methods.inc'
  include 'objects.mass_distributions.beta_profile.methods.inc'
  include 'objects.mass_distributions.Hernquist.methods.inc'
  include 'objects.mass_distributions.Sersic.methods.inc'

end module Mass_Distributions

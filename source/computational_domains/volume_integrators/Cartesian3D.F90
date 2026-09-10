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
  <computationalDomainVolumeIntegrator name="computationalDomainVolumeIntegratorCartesian3D" docformat="rst">
   <description>
   Computes volume integrals over three-dimensional Cartesian grid cells, with the spatial extent of the domain defined by ``[xBoundaries]``, ``[yBoundaries]``, and ``[zBoundaries]``. Each cell volume is the product of its axis-aligned extents, enabling accurate integration of physical quantities over the full Cartesian domain.
   </description>
  </computationalDomainVolumeIntegrator>
  !!]
  type, extends(computationalDomainVolumeIntegratorClass) :: computationalDomainVolumeIntegratorCartesian3D
     !!{RST
     Implementation of a computational domain for 3D Cartesian cells.
     !!}
     private
     double precision, dimension(  2) :: xBoundaries      , yBoundaries      , &
          &                              zBoundaries
     double precision, dimension(3,2) :: boundaries
     double precision                 :: volume_          , toleranceRelative, &
          &                              toleranceAbsolute
   contains
     procedure :: volume       => cartesian3DVolume
     procedure :: toleranceSet => cartesian3DToleranceSet
     procedure :: integrate    => cartesian3DIntegrate
  end type computationalDomainVolumeIntegratorCartesian3D

  interface computationalDomainVolumeIntegratorCartesian3D
     !!{RST
     Constructors for the :galacticus-class:`computationalDomainVolumeIntegratorCartesian3D` computational domain.
     !!}
     module procedure cartesian3DConstructorParameters
     module procedure cartesian3DConstructorInternal
  end interface computationalDomainVolumeIntegratorCartesian3D

contains

  function cartesian3DConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`computationalDomainVolumeIntegratorCartesian3D` computational domain volume integrator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (computationalDomainVolumeIntegratorCartesian3D)                 :: self
    type            (inputParameters                               ), intent(inout)  :: parameters
    double precision                                                , dimension(  2) :: xBoundaries             , yBoundaries, &
         &                                                                              zBoundaries
    double precision                                                , dimension(3,2) :: boundaries

    !![
    <inputParameter docformat="rst">
      <name>xBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>
      A two-element array :math:`[x_\mathrm{min}, x_\mathrm{max}]` specifying the extent of the 3D Cartesian integration domain along the :math:`x`-axis.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>yBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>
      A two-element array :math:`[y_\mathrm{min}, y_\mathrm{max}]` specifying the extent of the 3D Cartesian integration domain along the :math:`y`-axis.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>zBoundaries</name>
      <defaultValue>[-1.0d0,+1.0d0]</defaultValue>
      <description>
      A two-element array :math:`[z_\mathrm{min}, z_\mathrm{max}]` specifying the extent of the 3D Cartesian integration domain along the :math:`z`-axis.
      </description>
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

  function cartesian3DConstructorInternal(boundaries,toleranceAbsolute,toleranceRelative) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`computationalDomainVolumeIntegratorCartesian3D` computational domain volume integrator class. The optional  toleranceAbsolute and  toleranceRelative arguments specify the tolerances to which the volume integral is to be evaluated---the absolute tolerance is apportioned between the nested one-dimensional integrals in proportion to the measure which each is integrated against by the level above it.
    !!}
    implicit none
    type            (computationalDomainVolumeIntegratorCartesian3D)                                :: self
    double precision                                                , dimension(3,2), intent(in   ) :: boundaries
    double precision                                                , optional      , intent(in   ) :: toleranceAbsolute, toleranceRelative
    !![
    <constructorAssign variables="boundaries"/>
    <optionalArgument name="toleranceAbsolute" defaultsTo="0.0d0"                          />
    <optionalArgument name="toleranceRelative" defaultsTo="toleranceRelativeVolumeIntegral"/>
    !!]

    self%toleranceAbsolute=toleranceAbsolute_
    self%toleranceRelative=toleranceRelative_
    self%xBoundaries      =boundaries(1,:)
    self%yBoundaries      =boundaries(2,:)
    self%zBoundaries      =boundaries(3,:)
    self%volume_          =+(boundaries(1,2)-boundaries(1,1)) &
         &                 *(boundaries(2,2)-boundaries(2,1)) &
         &                 *(boundaries(3,2)-boundaries(3,1)) 
    return
  end function cartesian3DConstructorInternal

  double precision function cartesian3DVolume(self)
    !!{RST
    Return the volume of the computational domain cell.
    !!}
    implicit none
    class(computationalDomainVolumeIntegratorCartesian3D), intent(inout) :: self

    cartesian3DVolume=self%volume_
    return
  end function cartesian3DVolume

  subroutine cartesian3DToleranceSet(self,toleranceAbsolute,toleranceRelative)
    !!{RST
    Set the tolerances to which the volume integral over the computational domain cell is to be evaluated.
    !!}
    implicit none
    class           (computationalDomainVolumeIntegratorCartesian3D), intent(inout) :: self
    double precision                                                , intent(in   ) :: toleranceAbsolute, toleranceRelative

    self%toleranceAbsolute=toleranceAbsolute
    self%toleranceRelative=toleranceRelative
    return
  end subroutine cartesian3DToleranceSet

  double precision function cartesian3DIntegrate(self,integrand)
    !!{RST
    Integrate over the computational domain cell.
    !!}
    use :: Numerical_Integration, only : integrator         , GSL_Integ_Gauss15
    use :: Coordinates          , only : coordinateCartesian
    implicit none
    class           (computationalDomainVolumeIntegratorCartesian3D), intent(inout), target :: self
    procedure       (computationalDomainVolumeIntegrand            )                        :: integrand
    type            (integrator                                    )                        :: integratorX       , integratorY       , &
         &                                                                                     integratorZ
    type            (coordinateCartesian                           )                        :: coordinates
    double precision                                                                        :: toleranceAbsoluteX, toleranceAbsoluteY, &
         &                                                                                     toleranceAbsoluteZ

    ! Apportion the absolute tolerance between the nested integrals. The :math:`x`-integral is required to meet the requested
    ! absolute tolerance directly. Each inner integral supplies the integrand of the integral one level out, and so need only be
    ! accurate to the tolerance of that outer integral divided by the measure, ∫dx and ∫dy respectively, against which it is there
    ! integrated.
    toleranceAbsoluteX  =+self%toleranceAbsolute
    toleranceAbsoluteY  =+toleranceAbsoluteX          &
         &               /(                           &
         &                 +self%boundaries(1,2)      &
         &                 -self%boundaries(1,1)      &
         &                )
    toleranceAbsoluteZ  =+toleranceAbsoluteY          &
         &               /(                           &
         &                 +self%boundaries(2,2)      &
         &                 -self%boundaries(2,1)      &
         &                )
    ! Construct the integrators. These are built once here, rather than in the integrands where they would be rebuilt (along with
    ! their GSL workspaces) on every evaluation.
    integratorX         =integrator(cartesian3DIntegrandX,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsoluteX,integrationRule=GSL_Integ_Gauss15)
    integratorY         =integrator(cartesian3DIntegrandY,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsoluteY,integrationRule=GSL_Integ_Gauss15)
    integratorZ         =integrator(cartesian3DIntegrandZ,toleranceRelative=self%toleranceRelative,toleranceAbsolute=toleranceAbsoluteZ,integrationRule=GSL_Integ_Gauss15)
    cartesian3DIntegrate=+integratorX%integrate(                      &
         &                                      self%boundaries(1,1), &
         &                                      self%boundaries(1,2)  &
         &                                     )
    return

  contains

    double precision function cartesian3DIntegrandX(x)
      !!{RST
      :math:`x`-integrand over Cartesian 3D computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: x

      call coordinates%xSet(x)
      cartesian3DIntegrandX=+integratorY%integrate(                      &
           &                                       self%boundaries(2,1), &
           &                                       self%boundaries(2,2)  &
           &                                      )
      return
    end function cartesian3DIntegrandX

    double precision function cartesian3DIntegrandY(y)
      !!{RST
      :math:`y`-integrand over Cartesian 3D computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: y

      call coordinates%ySet(y)
      cartesian3DIntegrandY=+integratorZ%integrate(                      &
           &                                       self%boundaries(3,1), &
           &                                       self%boundaries(3,2)  &
           &                                      )
      return
    end function cartesian3DIntegrandY

    double precision function cartesian3DIntegrandZ(z)
      !!{RST
      :math:`z`-integrand over Cartesian 3D computational domain cells.
      !!}
      implicit none
      double precision, intent(in   ) :: z

      call coordinates%zSet(z)
      cartesian3DIntegrandZ=integrand(coordinates)
      return
    end function cartesian3DIntegrandZ

  end function cartesian3DIntegrate

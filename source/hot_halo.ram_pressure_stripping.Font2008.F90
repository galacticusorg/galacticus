!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a model of ram pressure stripping of hot halos based on the methods of
!% \cite{font_colours_2008}.

module Hot_Halo_Ram_Pressure_Stripping_Font2008
  !% Implements a module which implements a model of ram pressure stripping of hot halos based on the methods of
  !% \cite{font_colours_2008}.
  use Galacticus_Nodes
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Stripping_Font2008_Initialize

  ! Pointers to the host and satellite nodes.
  type(treeNode),   pointer :: hostNode,satelliteNode
  !$omp threadprivate(hostNode,satelliteNode)

  ! The ram pressure force (per unit area) used in root finding.
  double precision          :: ramPressureForce
  !$omp threadprivate(ramPressureForce)

  ! Parameters of the ram pressure stripping model.
  double precision          :: ramPressureStrippingFormFactor

contains

  !# <hotHaloRamPressureStrippingMethod>
  !#  <unitName>Hot_Halo_Ram_Pressure_Stripping_Font2008_Initialize</unitName>
  !# </hotHaloRamPressureStrippingMethod>
  subroutine Hot_Halo_Ram_Pressure_Stripping_Font2008_Initialize(hotHaloRamPressureStrippingMethod,Hot_Halo_Ram_Pressure_Stripping_Get)
    !% Initializes the ``Font2008'' hot halo ram pressure stripping module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in   ) :: hotHaloRamPressureStrippingMethod
    procedure(double precision), pointer, intent(inout) :: Hot_Halo_Ram_Pressure_Stripping_Get
    
    if (hotHaloRamPressureStrippingMethod == 'Font2008') then
       Hot_Halo_Ram_Pressure_Stripping_Get => Hot_Halo_Ram_Pressure_Stripping_Font2008_Get
       !@ <inputParameter>
       !@   <name>ramPressureStrippingFormFactor</name>
       !@   <defaultValue>2</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The form factor appearing in the gravitational binding force (per unit area) in the ram pressure stripping model
       !@     of \citeauthor{font_colours_2008}~(\citeyear{font_colours_2008}; their eqn.~4).
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('ramPressureStrippingFormFactor',ramPressureStrippingFormFactor,defaultValue=2.0d0)
    end if
    return
  end subroutine Hot_Halo_Ram_Pressure_Stripping_Font2008_Initialize

  double precision function Hot_Halo_Ram_Pressure_Stripping_Font2008_Get(thisNode)
    !% Computes the hot halo ram pressure stripping radius, assuming a null calculation in which that radius always equals the
    !% virial radius.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Kepler_Orbits
    use Satellite_Orbits
    use Hot_Halo_Density_Profile
    use Root_Finder
    use FGSL
    use, intrinsic :: ISO_C_Binding
    implicit none
    type (treeNode              ), intent(inout), pointer :: thisNode
    class(nodeComponentSatellite),                pointer :: thisSatelliteComponent
    type (fgsl_function         ), save                   :: rootFunction
    type (fgsl_root_fsolver     ), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision             , parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3
    double precision             , parameter              :: radiusSmallestOverRadiusVirial=1.0d-6
    type (c_ptr                 )                         :: parameterPointer
    type (keplerOrbit           )                         :: thisOrbit
    double precision                                      :: virialRadius,orbitalRadius,orbitalVelocity,densityHotHaloHost&
         &,radiusMinimum,radiusMaximum

    ! Get the virial radius of the satellite.
    virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
    ! Test whether thisNode is a satellite.
    if (thisNode%isSatellite()) then
       ! Find the host node.
       hostNode      => thisNode%parent
       ! Set a pointer to the satellite node.
       satelliteNode => thisNode
       ! Get the satellite component.
       thisSatelliteComponent => thisNode%satellite() 
       ! Get the orbit for this node.
       thisOrbit=thisSatelliteComponent%virialOrbit()
       ! Get the orbital radius and velocity at pericenter.
       call Satellite_Orbit_Pericenter_Phase_Space_Coordinates(hostNode,thisOrbit,orbitalRadius,orbitalVelocity)
       ! Find the density of the host node hot halo at the pericentric radius.
       densityHotHaloHost=Hot_Halo_Density(hostNode,orbitalRadius)
       ! Find the ram pressure force at pericenter.
       ramPressureForce=densityHotHaloHost*orbitalVelocity**2
       ! Find the radial range within which the pericenter must lie. The pericenter must be smaller than (or equal to) the
       ! current radius of the orbit.
       if      (Ram_Pressure_Stripping_Radius_Solver(                               virialRadius,parameterPointer) >= 0.0d0) then
          ! The ram pressure force is not sufficiently strong to strip even at the satellite virial radius - simply return the
          ! virial radius as the stripping radius in this case.
          Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=virialRadius
       else if (Ram_Pressure_Stripping_Radius_Solver(radiusSmallestOverRadiusVirial*virialRadius,parameterPointer) <= 0.0d0) then
          ! The ram pressure force can strip to (essentially) arbitrarily small radii.
          Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=0.0d0
       else
          radiusMinimum=virialRadius
          radiusMaximum=virialRadius
          do while (Ram_Pressure_Stripping_Radius_Solver(radiusMinimum,parameterPointer) <= 0.0d0 .and. radiusMinimum > radiusSmallestOverRadiusVirial*virialRadius)
             radiusMinimum=0.5d0*radiusMinimum
          end do
          if (radiusMinimum > radiusSmallestOverRadiusVirial*virialRadius) then
             ! Solve for the pericentric radius.
             Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=Root_Find(radiusMinimum,radiusMaximum,Ram_Pressure_Stripping_Radius_Solver&
                  &,parameterPointer,rootFunction,rootFunctionSolver ,toleranceAbsolute,toleranceRelative)
          else
             Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=0.0d0
          end if
       end if
    else
       ! If thisNode is not a satellite, return a ram pressure stripping radius equal to the virial radius.
       Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=virialRadius
    end if
    return
  end function Hot_Halo_Ram_Pressure_Stripping_Font2008_Get
  
  function Ram_Pressure_Stripping_Radius_Solver(radius,parameterPointer) bind(c)
    !% Root function used in finding the ram pressure stripping radius.
    use, intrinsic :: ISO_C_Binding
    use Dark_Matter_Profiles
    use Hot_Halo_Density_Profile
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    real(c_double), value :: radius
    type(c_ptr),    value :: parameterPointer
    real(c_double)        :: Ram_Pressure_Stripping_Radius_Solver
    double precision      :: enclosedMass,hotHaloDensity,gravitationalBindingForce

    enclosedMass             =Galactic_Structure_Enclosed_Mass(satelliteNode,radius,massType=massTypeAll,componentType=componentTypeAll)
    hotHaloDensity           =Hot_Halo_Density(satelliteNode,radius)
    gravitationalBindingForce=ramPressureStrippingFormFactor*gravitationalConstantGalacticus*enclosedMass*hotHaloDensity/radius
    if (gravitationalBindingForce >= 0.0d0) then
       Ram_Pressure_Stripping_Radius_Solver=gravitationalBindingForce-ramPressureForce
    else
       Ram_Pressure_Stripping_Radius_Solver=                         -ramPressureForce
    end if
    return
  end function Ram_Pressure_Stripping_Radius_Solver

end module Hot_Halo_Ram_Pressure_Stripping_Font2008

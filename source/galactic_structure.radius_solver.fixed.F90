!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements a ``fixed'' galactic radii solver in which sizes are always equal to
!% the halo virial or turnaround radius multiplied by its spin parameter and a multiplicative constant.

module Galactic_Structure_Radii_Fixed
  !% Implements a ``fixed'' galactic radii solver in which sizes are always equal to
  !% the halo virial or turnaround radius multiplied by its spin parameter and a multiplicative constant.
  use Galactic_Structure_Radius_Solver_Procedures
  implicit none
  private
  public :: Galactic_Structure_Radii_Fixed_Initialize

  ! The ratio of galaxy size to the product of spin parameter and virial radius.
  double precision            :: galacticStructureRadiiFixedFactor

  ! Radius option.
  integer                     :: galacticStructureRadiiFixedRadius
  integer         , parameter :: galacticStructureRadiiFixedRadiusVirial    =0
  integer         , parameter :: galacticStructureRadiiFixedRadiusTurnaround=1

contains

  !# <galacticStructureRadiusSolverMethod>
  !#  <unitName>Galactic_Structure_Radii_Fixed_Initialize</unitName>
  !# </galacticStructureRadiusSolverMethod>
  subroutine Galactic_Structure_Radii_Fixed_Initialize(galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do,Galactic_Structure_Radii_Revert_Do)
    !% Initializes the ``fixed'' galactic radii solver module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Nodes, only : defaultBasicComponent
    use Galacticus_Error
    implicit none
    type     (varying_string                       ), intent(in   )          :: galacticStructureRadiusSolverMethod
    procedure(Galactic_Structure_Radii_Solve_Fixed ), intent(inout), pointer :: Galactic_Structure_Radii_Solve_Do
    procedure(Galactic_Structure_Radii_Revert_Fixed), intent(inout), pointer :: Galactic_Structure_Radii_Revert_Do
    type     (varying_string                       )                         :: galacticStructureRadiiFixedRadiusText

    if (galacticStructureRadiusSolverMethod == 'fixed') then
       Galactic_Structure_Radii_Solve_Do  => Galactic_Structure_Radii_Solve_Fixed
       Galactic_Structure_Radii_Revert_Do => Galactic_Structure_Radii_Revert_Fixed
       !# <inputParameter>
       !#   <name>galacticStructureRadiiFixedFactor</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultSource>\citep{mo_formation_1998}</defaultSource>
       !#   <defaultValue>sqrt(0.5d0)</defaultValue>
       !#   <description>The ratio of galaxy radius to $\lambda r_\mathrm{vir}$ in the ``fixed'' galactic structure radius solver algorithm.</description>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>galacticStructureRadiiFixedRadius</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>var_str('virial')</defaultValue>
       !#   <description>The radius to use in the ``fixed'' galactic structure radius solver algorithm. Allowed options are ``virial'' and ``turnaround''.</description>
       !#   <source>globalParameters</source>
       !#   <type>string</type>
       !#   <variable>galacticStructureRadiiFixedRadiusText</variable>
       !# </inputParameter>
       select case (char(galacticStructureRadiiFixedRadiusText))
       case ('virial'    )
          galacticStructureRadiiFixedRadius=galacticStructureRadiiFixedRadiusVirial
       case ('turnaround')
          galacticStructureRadiiFixedRadius=galacticStructureRadiiFixedRadiusTurnaround
          if (.not.defaultBasicComponent%radiusTurnaroundIsGettable())                                                         &
               & call Galacticus_Error_Report                                                                                  &
               &   (                                                                                                           &
               &    'the "radiusTurnaround" property of the basic component must be gettable.'                              // &
               &    Galacticus_Component_List(                                                                                 &
               &                              'basic'                                                                        , &
               &                                defaultBasicComponent%radiusTurnaroundAttributeMatch(requireGettable=.true.)   &
               &                             )                                                                              // &
               &    {introspection:location}                                                                                   &
               &   )
        case default
          call Galacticus_Error_Report('[galacticStructureRadiiFixedRadius] must be either "virial" or "turnaround"'//{introspection:location})
       end select
    end if
    return
  end subroutine Galactic_Structure_Radii_Fixed_Initialize

  subroutine Galactic_Structure_Radii_Solve_Fixed(node)
    !% Find the radii of galactic components in {\normalfont \ttfamily node} using the ``fixed'' method.
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    implicit none
    type            (treeNode                  ), intent(inout), target :: node
    procedure       (Radius_Solver_Get_Template), pointer               :: Radius_Get                     => null(), Velocity_Get => null()
    procedure       (Radius_Solver_Set_Template), pointer               :: Radius_Set                     => null(), Velocity_Set => null()
    !$omp threadprivate(Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    logical                                     , parameter             :: specificAngularMomentumRequired=  .false.
    logical                                                             :: componentActive
    double precision                                                    :: specificAngularMomentum

    ! Determine if the node is physically plausible.
    node%isPhysicallyPlausible=.true.
    node%isSolvable           =.true.
    include 'galactic_structure.radius_solver.plausible.inc'
    if (.not.node%isPhysicallyPlausible) return
    
    ! Solve for radii.
    include 'galactic_structure.radius_solver.tasks.inc'

    return
  end subroutine Galactic_Structure_Radii_Solve_Fixed

  subroutine Solve_For_Radius(node,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    !% Solve for the equilibrium radius of the given component.
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profiles
    use Galacticus_Nodes       , only : treeNode, nodeComponentSpin, nodeComponentBasic
    implicit none
    type            (treeNode                  ), intent(inout)          :: node
    double precision                            , intent(in   )          :: specificAngularMomentum
    procedure       (Radius_Solver_Get_Template), intent(in   ), pointer :: Radius_Get             , Velocity_Get
    procedure       (Radius_Solver_Set_Template), intent(in   ), pointer :: Radius_Set             , Velocity_Set
    class           (nodeComponentSpin         )               , pointer :: thisSpinComponent
    class           (nodeComponentBasic        )               , pointer :: basic
    class           (darkMatterHaloScaleClass  )               , pointer :: darkMatterHaloScale_
    class           (darkMatterProfileClass    )               , pointer :: darkMatterProfile_
    double precision                                                     :: radius                 , velocity
    !GCC$ attributes unused :: Radius_Get, Velocity_Get, specificAngularMomentum
    
    ! Find the radius of the component, assuming radius is a fixed fraction of radius times spin parameter.
    thisSpinComponent    => node%spin      ()
    select case (galacticStructureRadiiFixedRadius)
    case (galacticStructureRadiiFixedRadiusVirial    )
       darkMatterHaloScale_ => darkMatterHaloScale                         (    )
       velocity             =  darkMatterHaloScale_%virialVelocity         (node)
       radius               =  darkMatterHaloScale_%virialRadius           (node)*thisSpinComponent%spin()*galacticStructureRadiiFixedFactor
    case (galacticStructureRadiiFixedRadiusTurnaround)
       darkMatterProfile_   => darkMatterProfile                           (    )
       basic                => node                %basic                  (    )
       velocity             =  darkMatterProfile_  %circularVelocityMaximum(node)
       radius               =  basic               %radiusTurnaround       (    )*thisSpinComponent%spin()*galacticStructureRadiiFixedFactor
    end select
    ! Set the component size to new radius and velocity.
    call Radius_Set  (node,radius  )
    call Velocity_Set(node,velocity)
    return
  end subroutine Solve_For_Radius

  subroutine Galactic_Structure_Radii_Revert_Fixed(node)
    !% Revert radii for the fixed galactic structure solve. Not necessary for this algorithm.
    implicit none
    type(treeNode), intent(inout), target :: node
    !GCC$ attributes unused :: node
    
    return
  end subroutine Galactic_Structure_Radii_Revert_Fixed

end module Galactic_Structure_Radii_Fixed

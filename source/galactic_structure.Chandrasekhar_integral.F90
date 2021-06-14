!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!+    Contributions to this file made by: Andrew Benson.

!% Contains a module which implements calculations of the integral appearing in the \cite{chandrasekhar_dynamical_1943} dynamical
!% friction model.

module Galactic_Structure_Chandrasekhar_Integrals
  !% Implements calculations of the integral appearing in the \cite{chandrasekhar_dynamical_1943} dynamical friction model:
  !% \begin{equation}
  !%  \rho(\boldsymbol{x}_\mathrm{s}) \int \mathrm{d}\boldsymbol{v} f(\boldsymbol{v}) {\boldsymbol{v}-\boldsymbol{v}_\mathrm{s} \over |\boldsymbol{v}-\boldsymbol{v}_\mathrm{s}|^3},
  !% \end{equation}  
  !% where $\rho(\boldsymbol{x}_\mathrm{s})$ is the density at the position of the perturber, $\boldsymbol{x}_\mathrm{s}$,
  !% $f(\boldsymbol{v})$ is the velocity distribution function at velocity $\boldsymbol{v}$, and $\boldsymbol{v}_\mathrm{s}$ is
  !% the velocity of the perturber.
  use :: Kind_Numbers, only : kind_int8
  implicit none
  private
  public :: Galactic_Structure_Chandrasekhar_Integral

  ! Module scope variables used in mapping over components.
  integer                        :: componentTypeShared    , massTypeShared
  double precision, dimension(3) :: positionCartesianShared, velocityCartesianShared
  !$omp threadprivate(massTypeShared,componentTypeShared,positionCartesianShared,velocityCartesianShared)

contains

  function Galactic_Structure_Chandrasekhar_Integral(node,positionCartesian,velocityCartesian,componentType,massType)
    !% Compute the \cite{chandrasekhar_dynamical_1943} integral at a given position and velocity.
    use :: Galactic_Structure_Options, only : componentTypeAll                         , massTypeAll
    use :: Galacticus_Nodes          , only : optimizeForChandrasekharIntegralSummation, reductionSummation, treeNode
    !# <include directive="chandrasekharIntegralTask" type="moduleUse">
    include 'galactic_structure.chandrasekharIntegral.tasks.modules.inc'
    !# </include>
    implicit none
    double precision                                                 , dimension(3) :: Galactic_Structure_Chandrasekhar_Integral
    type            (treeNode                        ), intent(inout)               :: node
    double precision                                  , intent(in   ), dimension(3) :: positionCartesian                          , velocityCartesian
    integer                                           , intent(in   ), optional     :: componentType                              , massType
    integer                                           , parameter                   :: chandrasekharIntegralSize                =3
    procedure       (Component_Chandrasekhar_Integral), pointer                     :: componentChandrasekharIntegralFunction
    double precision                                                 , dimension(3) :: componentChandrasekharIntegral

    ! Copy the position and velocity to module-scope.
    positionCartesianShared=positionCartesian
    velocityCartesianShared=velocityCartesian
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared=massType
    else
       massTypeShared=massTypeAll
    end if
    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Initialize pointer to function that supplies the Chandrasekhar integral for all components.
    componentChandrasekharIntegralFunction => Component_Chandrasekhar_Integral
    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared     =massType
    else
       massTypeShared     =massTypeAll
    end if
    ! Compute the Chandrasekhar integral.
    Galactic_Structure_Chandrasekhar_Integral=node%mapDouble1(componentChandrasekharIntegralFunction,chandrasekharIntegralSize,reductionSummation,optimizeFor=optimizeForChandrasekharIntegralSummation)
    !# <include directive="chandrasekharIntegralTask" type="functionCall" functionType="function" returnParameter="componentChandrasekharIntegral">
    !#  <functionArgs>node,positionCartesianShared,velocityCartesianShared,componentTypeShared,massTypeShared</functionArgs>
    !#  <onReturn>Galactic_Structure_Chandrasekhar_Integral=Galactic_Structure_Chandrasekhar_Integral+componentChandrasekharIntegral</onReturn>
    include 'galactic_structure.chandrasekharIntegral.tasks.inc'
    !# </include>
    return
  end function Galactic_Structure_Chandrasekhar_Integral

  function Component_Chandrasekhar_Integral(component,resultSize)
    !% Function returning the Chandrasekhar integral in a component. Suitable for mapping over components.
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    integer                        , intent(in   )         :: resultSize
    class           (nodeComponent), intent(inout)         :: component
    double precision               , dimension(resultSize) :: Component_Chandrasekhar_Integral

    Component_Chandrasekhar_Integral=component%chandrasekharIntegral(positionCartesianShared,velocityCartesianShared,componentTypeShared,massTypeShared)
    return
  end function Component_Chandrasekhar_Integral

end module Galactic_Structure_Chandrasekhar_Integrals

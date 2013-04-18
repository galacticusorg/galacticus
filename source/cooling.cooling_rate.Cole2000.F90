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

!% Contains a module which implements a \cite{cole_hierarchical_2000} cooling rate calculation.

module Cooling_Rates_Cole2000
  !% Implements a \cite{cole_hierarchical_2000} cooling rate calculation.
  use, intrinsic :: ISO_C_Binding
  use Galacticus_Nodes
  implicit none
  private
  public :: Cooling_Rate_Cole2000_Initialize

contains

  !# <coolingRateMethod>
  !#  <unitName>Cooling_Rate_Cole2000_Initialize</unitName>
  !# </coolingRateMethod>
  subroutine Cooling_Rate_Cole2000_Initialize(coolingRateMethod,Cooling_Rate_Get)
    !% Initializes the ``Cole et al. (2000)'' cooling rate module.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: coolingRateMethod
    procedure(Cooling_Rate_Cole2000), pointer, intent(inout) :: Cooling_Rate_Get
    
    ! Return a pointer to our implementation.
    if (coolingRateMethod == 'Cole2000') Cooling_Rate_Get => Cooling_Rate_Cole2000
    return
  end subroutine Cooling_Rate_Cole2000_Initialize

  double precision function Cooling_Rate_Cole2000(thisNode)
    !% Computes the mass cooling rate in a hot gas halo utilizing the \cite{cole_hierarchical_2000} method. This is based on the
    !% properties of the halo at formation time, and gives a zero cooling rate when the cooling radius exceeds the virial radius.
    use Galacticus_Nodes
    use Cooling_Infall_Radii
    use Numerical_Constants_Math
    use Hot_Halo_Density_Profile
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic  ),                pointer :: formationBasicComponent
    class(nodeComponentHotHalo),                pointer :: formationHotHaloComponent
    double precision                                    :: infallRadius,infallDensity,outerRadius,infallRadiusGrowthRate
    
    ! Get node components.
    formationBasicComponent   => thisNode%formationNode%basic  ()
    formationHotHaloComponent => thisNode%formationNode%hotHalo()

    ! Check for empty halos.
    if (formationHotHaloComponent%mass() <= 0.0d0) then
       Cooling_Rate_Cole2000=0.0d0
       return
    end if

    ! Get the outer radius of the hot halo.
    outerRadius=formationHotHaloComponent%outerRadius()

    ! Get the infall radius.
    infallRadius=Infall_Radius                 (thisNode%formationNode)

    if (infallRadius >= outerRadius) then
       ! Cooling radius exceeds the outer radius - zero infall rate.
       Cooling_Rate_Cole2000=0.0d0
    else
       ! Find the density at the infall radius.
       infallDensity=Hot_Halo_Density                  (thisNode%formationNode,infallRadius)
       ! Find infall radius growth rate.
       infallRadiusGrowthRate=Infall_Radius_Growth_Rate(thisNode%formationNode             )
       ! Compute the infall rate.
       Cooling_Rate_Cole2000=4.0d0*Pi*(infallRadius**2)*infallDensity*infallRadiusGrowthRate
    end if
    return
  end function Cooling_Rate_Cole2000

end module Cooling_Rates_Cole2000

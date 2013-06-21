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

!% Contains a module which implements a \cite{white_galaxy_1991} cooling rate calculation.

module Cooling_Rates_White_Frenk
  !% Implements a \cite{white_galaxy_1991} cooling rate calculation.
  use, intrinsic :: ISO_C_Binding
  use Galacticus_Nodes
  implicit none
  private
  public :: Cooling_Rate_White_Frenk_Initialize

  ! Velocity (in km/s) above which cooling rates in gas are forced to zero.
  double precision :: zeroCoolingRateAboveVelocity

contains

  !# <coolingRateMethod>
  !#  <unitName>Cooling_Rate_White_Frenk_Initialize</unitName>
  !# </coolingRateMethod>
  subroutine Cooling_Rate_White_Frenk_Initialize(coolingRateMethod,Cooling_Rate_Get)
    !% Initializes the ``White-Frenk1991'' cooling rate module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type     (varying_string              ), intent(in   )          :: coolingRateMethod
    procedure(Cooling_Rate_White_Frenk    ), intent(inout), pointer :: Cooling_Rate_Get

    if (coolingRateMethod == 'White-Frenk1991') then
       Cooling_Rate_Get => Cooling_Rate_White_Frenk

       ! Get cooling rate parameters.
       !@ <inputParameter>
       !@   <name>zeroCoolingRateAboveVelocity</name>
       !@   <defaultValue>10000</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The halo virial velocity (in km/s) above which cooling rates are forced to zero in the {\tt White-Frenk1991} cooling rate model.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('zeroCoolingRateAboveVelocity',zeroCoolingRateAboveVelocity,defaultValue=1.0d4)

       ! Check that the properties we need are gettable.
       if (.not.defaultHotHaloComponent%massIsGettable())                                                  &
            & call Galacticus_Error_Report(                                                                &
            &                              'Cooling_Rate_White_Frenk_Initialize'                         , &
            &                              'mass property of hot halo component must be gettable'          &
            &                             )
       if (.not.defaultHotHaloComponent%outerRadiusIsGettable())                                           &
            & call Galacticus_Error_Report(                                                                &
            &                              'Cooling_Rate_White_Frenk_Initialize'                         , &
            &                              'outer radius property of hot halo component must be gettable'  &
            &                             )
    end if
    return
  end subroutine Cooling_Rate_White_Frenk_Initialize

  double precision function Cooling_Rate_White_Frenk(thisNode)
    !% Computes the mass cooling rate in a hot gas halo utilizing the \cite{white_galaxy_1991} method.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Cooling_Times_Available
    use Cooling_Infall_Radii
    use Numerical_Constants_Math
    use Hot_Halo_Density_Profile
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    double precision                                               :: coolingDensity      , infallRadius  , infallRadiusGrowthRate, &
         &                                                            outerRadius         , virialVelocity

    ! Get the virial velocity.
    virialVelocity=Dark_Matter_Halo_Virial_Velocity(thisNode)

    ! Return zero cooling rate if virial velocity exceeds critical value.
    if (virialVelocity > zeroCoolingRateAboveVelocity) then
       Cooling_Rate_White_Frenk=0.0d0
       return
    end if

    ! Get the outer radius of the hot halo.
    thisHotHaloComponent => thisNode            %hotHalo    ()
    outerRadius          =  thisHotHaloComponent%outerRadius()

    ! Get the cooling radius.
    infallRadius=Infall_Radius(thisNode)

    if (infallRadius >= outerRadius) then
       ! Cooling radius exceeds the outer radius. Limit infall to the dynamical timescale.
       Cooling_Rate_White_Frenk=thisHotHaloComponent%mass()/Dark_Matter_Halo_Dynamical_Timescale(thisNode)
    else
       ! Find the density at the cooling radius.
       coolingDensity=Hot_Halo_Density(thisNode,infallRadius)
       ! Find cooling radius growth rate.
       infallRadiusGrowthRate=Infall_Radius_Growth_Rate(thisNode)
       ! Compute the cooling rate.
       Cooling_Rate_White_Frenk=4.0d0*Pi*(infallRadius**2)*coolingDensity*infallRadiusGrowthRate
    end if
    return
  end function Cooling_Rate_White_Frenk

end module Cooling_Rates_White_Frenk

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

!% Contains a module which implements simple mass loss rates from spheroids due to ram pressure stripping.

module Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids_Simple
  !% Implements simple mass loss rates from spheroids due to ram pressure stripping.
  use Galacticus_Nodes
  implicit none
  private
  public :: Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids_Simple_Init

  ! Parameter controlling the maximum mass loss fraction per dynamical time.
  double precision :: ramPressureStrippingMassLossRateSpheroidSimpleFractionalRateMax

contains

  !# <ramPressureStrippingMassLossRateSpheroidsMethod>
  !#  <unitName>Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids_Simple_Init</unitName>
  !# </ramPressureStrippingMassLossRateSpheroidsMethod>
  subroutine Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids_Simple_Init(ramPressureStrippingMassLossRateSpheroidsMethod,Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Get)
    !% Initializes the ``simple'' ram pressure stripping mass loss rate from spheroids module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string  ),          intent(in   ) :: ramPressureStrippingMassLossRateSpheroidsMethod
    procedure(Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Simple), pointer, intent(inout) :: Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Get
    
    if (ramPressureStrippingMassLossRateSpheroidsMethod == 'simple') then
       Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Get => Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Simple
       !@ <inputParameter>
       !@   <name>ramPressureStrippingMassLossRateSpheroidSimpleFractionalRateMaximum</name>
       !@   <defaultValue>$10$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum fractional mass loss rate per dynamical time in the simple model of mass loss from spheroids due to ram pressure stripping.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('ramPressureStrippingMassLossRateSpheroidSimpleFractionalRateMaximum',ramPressureStrippingMassLossRateSpheroidSimpleFractionalRateMax,defaultValue=10.0d0)
    end if
    return
  end subroutine Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids_Simple_Init

  double precision function Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Simple(thisNode)
    !% Computes the mass loss rate from spheroids due to ram pressure stripping assuming a simple model. Specifically, the mass loss
    !% rate is
    !% \begin{equation}
    !% \dot{M}_{\rm gas} = -\alpha M_{\rm gas}/\tau_{\rm spheroid},
    !% \end{equation}
    !% where
    !% \begin{equation}
    !% \alpha = F_{\rm ram}/F_{\rm gravity},
    !% \end{equation}
    !% $F_{\rm ram}$ is the ram pressure force from the hot halo (see \S\ref{sec:HotHaloRamPressureForce}), and
    !% \begin{equation}
    !% F_{\rm gravity} = {4\over 3} \rho_{\rm gas}(r_{1/2}) {{\rm G} M_{\rm total}(r_{1/2})\over r_{1/2}}
    !% \end{equation}
    !% is the gravitational restoring force in the spheroid at the half-mass radius, $r_{\rm 1/2}$ \citep{takeda_ram_1984}.
    use Galacticus_Nodes
    use Hot_Halo_Ram_Pressure_Forces
    use Galactic_Structure_Options
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode             ), intent(inout), pointer :: thisNode
    class           (nodeComponentSpheroid),                pointer :: thisSpheroid
    double precision                                                :: forceRamPressure,radiusHalfMass,densityGas&
         &,massHalf,massLossRateFractional,forceGravitational,timeDynamical

    ! Assume no mass loss rate due to ram pressure by default.
    Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Simple=0.0d0
    ! Return immediately if this is not a satellite.
    if (.not.thisNode%isSatellite()) return
    ! Get the spheroid component.
    thisSpheroid       => thisNode%spheroid()
    ! Get the ram pressure force due to the hot halo.
    forceRamPressure   =  Hot_Halo_Ram_Pressure_Force(thisNode)
    ! Get the spheroid half-mass radius.
    radiusHalfMass     =  thisSpheroid%halfMassRadius()
    ! Compute the spheroid densities at the half mass radius.
    densityGas=Galactic_Structure_Density      (                                            &
         &                                      thisNode                                  , &
         &                                      [radiusHalfMass,0.0d0,0.0d0]              , &
         &                                      coordinateSystem=coordinateSystemSpherical, &
         &                                      massType        =massTypeGaseous          , &
         &                                      componentType   =componentTypeSpheroid      &
         &                                     )
    massHalf  =Galactic_Structure_Enclosed_Mass(                                            &
         &                                      thisNode                                  , &
         &                                      radiusHalfMass                            , &
         &                                      massType        =massTypeAll              , &
         &                                      componentType   =componentTypeSpheroid      &
         &                                     )
    ! Compute the gravitational restoring force in the spheroid midplane.
    forceGravitational =4.0*densityGas*gravitationalConstantGalacticus*massHalf/3.0d0/radiusHalfMass
    ! Return zero rate if the gravitational force is zero.
    if (forceGravitational <= 0.0d0) return
    ! Compute the mass loss fraction per dynamical time.
    if (forceRamPressure < ramPressureStrippingMassLossRateSpheroidSimpleFractionalRateMax*forceGravitational) then
       massLossRateFractional=forceRamPressure/forceGravitational
    else
       massLossRateFractional=ramPressureStrippingMassLossRateSpheroidSimpleFractionalRateMax
    end if
    ! Compute the dynamical time.
    timeDynamical         =(megaParsec/kilo/gigaYear)*thisSpheroid%radius()/thisSpheroid%velocity()
    ! Compute the mass loss rate.
    Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Simple=massLossRateFractional*thisSpheroid%massGas()/timeDynamical
    return
  end function Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroid_Simple
  
end module Ram_Pressure_Stripping_Mass_Loss_Rate_Spheroids_Simple

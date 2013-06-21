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

!% Contains a module which implements simple mass loss rates from disks due to ram pressure stripping.

module Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Simple
  !% Implements simple mass loss rates from disks due to ram pressure stripping.
  use Galacticus_Nodes
  implicit none
  private
  public :: Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Simple_Init

  ! Parameter controlling the maximum mass loss fraction per dynamical time.
  double precision :: ramPressureStrippingMassLossRateDiskSimpleFractionalRateMaximum 
  
contains

  !# <ramPressureStrippingMassLossRateDisksMethod>
  !#  <unitName>Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Simple_Init</unitName>
  !# </ramPressureStrippingMassLossRateDisksMethod>
  subroutine Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Simple_Init(ramPressureStrippingMassLossRateDisksMethod,Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Get)
    !% Initializes the ``simple'' ram pressure stripping mass loss rate from disks module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                                   ), intent(in   )          :: ramPressureStrippingMassLossRateDisksMethod    
    procedure(Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Simple), intent(inout), pointer :: Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Get 
    
    if (ramPressureStrippingMassLossRateDisksMethod == 'simple') then
       Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Get => Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Simple
       !@ <inputParameter>
       !@   <name>ramPressureStrippingMassLossRateDiskSimpleFractionalRateMaximum</name>
       !@   <defaultValue>$10$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum fractional mass loss rate per dynamical time in the simple model of mass loss from disks due to ram pressure stripping.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('ramPressureStrippingMassLossRateDiskSimpleFractionalRateMaximum',ramPressureStrippingMassLossRateDiskSimpleFractionalRateMaximum,defaultValue=10.0d0)
    end if
    return
  end subroutine Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Simple_Init

  double precision function Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Simple(thisNode)
    !% Computes the mass loss rate from disks due to ram pressure stripping assuming a simple model. Specifically, the mass loss
    !% rate is
    !% \begin{equation}
    !% \dot{M}_{\rm gas} = -\alpha M_{\rm gas}/\tau_{\rm disk},
    !% \end{equation}
    !% where
    !% \begin{equation}
    !% \alpha = F_{\rm ram}/F_{\rm gravity},
    !% \end{equation}
    !% $F_{\rm ram}$ is the ram pressure force from the hot halo (see \S\ref{sec:HotHaloRamPressureForce}), and
    !% \begin{equation}
    !% F_{\rm gravity} = 2 \pi {\rm G} \Sigma_{\rm gas}(r_{1/2}) \Sigma_{\rm total}(r_{1/2})
    !% \end{equation}
    !% is the gravitational restoring force in the disk at the half-mass radius, $r_{\rm 1/2}$.
    use Galacticus_Nodes
    use Hot_Halo_Ram_Pressure_Forces
    use Galactic_Structure_Options
    use Galactic_Structure_Surface_Densities
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode         ), intent(inout), pointer :: thisNode                                                         
    class           (nodeComponentDisk)               , pointer :: thisDisk                                                         
    double precision                                            :: forceGravitational, forceRamPressure , massLossRateFractional, & 
         &                                                         radiusHalfMass    , surfaceDensityGas, surfaceDensityTotal   , & 
         &                                                         timeDynamical                                                    
    
    ! Assume no mass loss rate due to ram pressure by default.
    Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Simple=0.0d0
    ! Return immediately if this is not a satellite.
    if (.not.thisNode%isSatellite()) return
    ! Get the disk component.
    thisDisk           => thisNode%disk()
    ! Get the ram pressure force due to the hot halo.
    forceRamPressure   =  Hot_Halo_Ram_Pressure_Force(thisNode)
    ! Get the disk half-mass radius.
    radiusHalfMass     =  thisDisk%halfMassRadius()
    ! Compute the disk densities at the half mass radius.
    surfaceDensityGas  =Galactic_Structure_Surface_Density(                                              &
         &                                                 thisNode                                    , &
         &                                                 [radiusHalfMass,0.0d0,0.0d0]                , &
         &                                                 coordinateSystem=coordinateSystemCylindrical, &
         &                                                 massType        =massTypeGaseous            , &
         &                                                 componentType   =componentTypeDisk            &
         &                                              )
    surfaceDensityTotal=Galactic_Structure_Surface_Density(                                              &
         &                                                 thisNode                                    , &
         &                                                 [radiusHalfMass,0.0d0,0.0d0]                , &
         &                                                 coordinateSystem=coordinateSystemCylindrical, &
         &                                                 massType        =massTypeAll                , &
         &                                                 componentType   =componentTypeDisk            &
         &                                                )
    ! Compute the gravitational restoring force in the disk midplane.
    forceGravitational =2.0d0*Pi*gravitationalConstantGalacticus*surfaceDensityGas*surfaceDensityTotal
    ! Return zero rate if the gravitational force is zero.
    if (forceGravitational <= 0.0d0) return
    ! Compute the mass loss fraction per dynamical time.
    massLossRateFractional=min(forceRamPressure/forceGravitational,ramPressureStrippingMassLossRateDiskSimpleFractionalRateMaximum)
    ! Compute the dynamical time.
    timeDynamical         =(megaParsec/kilo/gigaYear)*thisDisk%radius()/thisDisk%velocity()
    ! Compute the mass loss rate.
    Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Simple=massLossRateFractional*thisDisk%massGas()/timeDynamical
    return
  end function Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Simple
  
end module Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Simple

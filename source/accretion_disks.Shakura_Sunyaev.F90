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

!% Contains a module which implements calculations of properties of thin Shakura-Sunyaev accretion disks.

module Accretion_Disks_Shakura_Sunyaev
  !% Implements calculations of properties of thin Shakura-Sunyaev accretion disks.
  implicit none
  private
  public :: Accretion_Disks_Shakura_Sunyaev_Initialize, Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev,&
       & Black_Hole_Spin_Up_Rate_Shakura_Sunyaev, Accretion_Disk_Jet_Power_Shakura_Sunyaev

contains

  !# <accretionDisksMethod>
  !#  <unitName>Accretion_Disks_Shakura_Sunyaev_Initialize</unitName>
  !# </accretionDisksMethod>
  subroutine Accretion_Disks_Shakura_Sunyaev_Initialize(accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get&
       &,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: accretionDisksMethod
    procedure(double precision), pointer, intent(inout) :: Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get&
         &,Accretion_Disk_Jet_Power_Get
    
    if (accretionDisksMethod == 'Shakura-Sunyaev') then
       Accretion_Disk_Radiative_Efficiency_Get => Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev
       Black_Hole_Spin_Up_Rate_Get             => Black_Hole_Spin_Up_Rate_Shakura_Sunyaev
       Accretion_Disk_Jet_Power_Get            => Accretion_Disk_Jet_Power_Shakura_Sunyaev
    end if
    return
  end subroutine Accretion_Disks_Shakura_Sunyaev_Initialize

  double precision function Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev(thisBlackHole,massAccretionRate)
    !% Computes the radiative efficiency for a Shakura-Sunyaev (thin) accretion disk.
    use Black_Hole_Fundamentals
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate

    Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev=1.0d0-Black_Hole_ISCO_Specific_Energy(thisBlackHole,units=unitsGravitational,orbit=orbitPrograde)
    return
  end function Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev

  double precision function Accretion_Disk_Jet_Power_Shakura_Sunyaev(thisBlackHole,massAccretionRate)
    !% Computes the jet power for a Shakura-Sunyaev (thin) accretion disk, using the expressions from
    !% \citeauthor{meier_association_2001}~(\citeyear{meier_association_2001}; his equations 4 and 5).
    use Galacticus_Nodes
    use Black_Hole_Fundamentals
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                        , parameter     :: alphaViscosity                =0.01d0
    double precision                        , parameter     :: alphaViscosityNormalized      =alphaViscosity/0.01d0
    double precision                        , parameter     :: powerNormalizationKerr        =(10.0d0**42.7)*ergs*gigaYear/massSolar/kilo**2
    double precision                        , parameter     :: powerNormalizationSchwarzchild=(10.0d0**41.7)*ergs*gigaYear/massSolar/kilo**2
    double precision                        , parameter     :: meierMassNormalization        =1.0d9
    double precision                                        :: blackHoleSpin,accretionRateDimensionless,blackHoleMassDimensionless

    ! Return immediately for non-positive accretion rates.
    if (massAccretionRate <= 0.0d0) then
       Accretion_Disk_Jet_Power_Shakura_Sunyaev=0.0d0
       return
    end if

    ! Get the black hole spin and dimensionless accretion rate and mass as defined by Meier (2001).
    blackHoleSpin=thisBlackHole%spin()
    accretionRateDimensionless=massAccretionRate/Black_Hole_Eddington_Accretion_Rate(thisBlackHole)
    blackHoleMassDimensionless=thisBlackHole%mass()/meierMassNormalization
    if (blackHoleMassDimensionless > 0.0d0 .and. accretionRateDimensionless > 0.0d0) then
       if (blackHoleSpin > 0.8d0) then
          ! Use Meier's rapidly rotating (Kerr) solution for high spin black holes.
          Accretion_Disk_Jet_Power_Shakura_Sunyaev=powerNormalizationKerr*(blackHoleMassDimensionless**0.9d0)&
               &*(accretionRateDimensionless**1.2d0)*(1.0d0+1.1d0*blackHoleSpin+0.29d0*blackHoleSpin**2)/(alphaViscosityNormalized&
               &**0.1d0)
       else
          ! Use Meier's Schwarzchild solution for low spin black holes. We, somewhat arbitrarily, interpolate
          ! from the Schwarzchild solution assuming power grows exponentially with spin and matched to the Kerr solution at the
          ! transition spin.
          Accretion_Disk_Jet_Power_Shakura_Sunyaev=powerNormalizationSchwarzchild*(blackHoleMassDimensionless**0.9d0) &
               &*(accretionRateDimensionless**1.2d0)*dexp(3.785d0*blackHoleSpin)/(alphaViscosityNormalized**0.1d0)
       end if
    end if
    return
  end function Accretion_Disk_Jet_Power_Shakura_Sunyaev

  double precision function Black_Hole_Spin_Up_Rate_Shakura_Sunyaev(thisBlackHole,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisBlackHole} due to accretion from a Shakura-Sunyaev (thin) accretion
    !% disk.
    use Galacticus_Nodes
    use Black_Hole_Fundamentals
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                                        :: spinToMassRateOfChangeRatio

    spinToMassRateOfChangeRatio=Black_Hole_ISCO_Specific_Angular_Momentum(thisBlackHole,units=unitsGravitational,orbit=orbitPrograde)&
         &-2.0d0*thisBlackHole%spin()*Black_Hole_ISCO_Specific_Energy(thisBlackHole,units=unitsGravitational,orbit&
         &=orbitPrograde)
    Black_Hole_Spin_Up_Rate_Shakura_Sunyaev=spinToMassRateOfChangeRatio*massAccretionRate/thisBlackHole%mass()
    return
  end function Black_Hole_Spin_Up_Rate_Shakura_Sunyaev

end module Accretion_Disks_Shakura_Sunyaev

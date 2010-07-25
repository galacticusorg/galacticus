!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements calculations of properties of thin Shakura-Sunyaev accretion disks.

module Accretion_Disks_Shakura_Sunyaev
  !% Implements calculations of properties of thin Shakura-Sunyaev accretion disks.
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
    procedure(),          pointer, intent(inout) :: Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get&
         &,Accretion_Disk_Jet_Power_Get
    
    if (accretionDisksMethod == 'Shakura-Sunyaev') then
       Accretion_Disk_Radiative_Efficiency_Get => Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev
       Black_Hole_Spin_Up_Rate_Get             => Black_Hole_Spin_Up_Rate_Shakura_Sunyaev
       Accretion_Disk_Jet_Power_Get            => Accretion_Disk_Jet_Power_Shakura_Sunyaev
    end if
    return
  end subroutine Accretion_Disks_Shakura_Sunyaev_Initialize

  double precision function Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev(thisNode,massAccretionRate)
    !% Computes the radiative efficiency for a Shakura-Sunyaev (thin) accretion disk.
    use Black_Hole_Fundamentals
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate

    Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev=1.0d0-Black_Hole_ISCO_Specific_Energy(thisNode,units=unitsGravitational,orbit=orbitPrograde)
    return
  end function Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev

  double precision function Accretion_Disk_Jet_Power_Shakura_Sunyaev(thisNode,massAccretionRate)
    !% Computes the jet power for a Shakura-Sunyaev (thin) accretion disk, using the expressions from
    !% \citeauthor{meier_association_2001}~(\citeyear{meier_association_2001}; his equations 4 and 5).
    use Tree_Nodes
    use Black_Hole_Fundamentals
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Units
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision, parameter              :: alphaViscosity                =0.01d0
    double precision, parameter              :: alphaViscosityNormalized      =alphaViscosity/0.01d0
    double precision, parameter              :: powerNormalizationKerr        =(10.0d0**42.7)*ergs*gigaYear/massSolar/kilo**2
    double precision, parameter              :: powerNormalizationSchwarzchild=(10.0d0**41.7)*ergs*gigaYear/massSolar/kilo**2
    double precision, parameter              :: meierMassNormalization        =1.0d9
    double precision                         :: blackHoleSpin,accretionRateDimensionless,blackHoleMassDimensionless

    ! Get the black hole spin and dimensionless accretion rate and mass as defined by Meier (2001).
    blackHoleSpin=Tree_Node_Black_Hole_Spin(thisNode)
    accretionRateDimensionless=massAccretionRate/Black_Hole_Eddington_Accretion_Rate(thisNode)
    blackHoleMassDimensionless=Tree_Node_Black_Hole_Mass(thisNode)/meierMassNormalization
    if (blackHoleMassDimensionless > 0.0d0 .and. accretionRateDimensionless > 0.0d0) then
       if (blackHoleSpin > 0.8d0) then
          ! Use Meier's rapidly rotating (Kerr) solution for high spin black holes.
          Accretion_Disk_Jet_Power_Shakura_Sunyaev=powerNormalizationKerr*(blackHoleMassDimensionless**0.9d0)&
               &*(accretionRateDimensionless**1.2d0)*(1.0d0+1.1d0*blackHoleSpin+0.29d0*blackHoleSpin**2)/(alphaViscosityNormalized&
               &**0.1d0)
       else
          ! Use Meier's rapidly rotating (Schwarzchild) solution for low spin black holes. We, somewhat arbitrarily, interpolate
          ! from the Schwarzchild solution assuming power grows exponentially with spin and matched to the Kerr solution at the
          ! transition spin.
         Accretion_Disk_Jet_Power_Shakura_Sunyaev=powerNormalizationSchwarzchild*(blackHoleMassDimensionless**0.9d0) &
              &*(accretionRateDimensionless**1.2d0)*dexp(3.785d0*blackHoleSpin)/(alphaViscosityNormalized**0.1d0)
      end if
    end if
    return
  end function Accretion_Disk_Jet_Power_Shakura_Sunyaev

  double precision function Black_Hole_Spin_Up_Rate_Shakura_Sunyaev(thisNode,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisNode} due to accretion from a Shakura-Sunyaev (thin) accretion
    !% disk.
    use Tree_Nodes
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision                         :: accretionEfficiency,spinToMassRateOfChangeRatio

    spinToMassRateOfChangeRatio=Black_Hole_ISCO_Specific_Angular_Momentum(thisNode,units=unitsGravitational,orbit=orbitPrograde)&
         &-2.0d0*Tree_Node_Black_Hole_Spin(thisNode)*Black_Hole_ISCO_Specific_Energy(thisNode,units=unitsGravitational,orbit&
         &=orbitPrograde)
    Black_Hole_Spin_Up_Rate_Shakura_Sunyaev=spinToMassRateOfChangeRatio*massAccretionRate/Tree_Node_Black_Hole_Mass(thisNode)
    return
  end function Black_Hole_Spin_Up_Rate_Shakura_Sunyaev

end module Accretion_Disks_Shakura_Sunyaev

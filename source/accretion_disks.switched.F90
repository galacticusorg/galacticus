!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations of properties of accretion disks which switch between thin and ADAF depening on
!% the accretion rate.

module Accretion_Disks_Switched
  !% Implements calculations of properties of accretion disks which switch between thin and ADAF depening on
  !% the accretion rate.
  use Accretion_Disks_ADAF
  use Accretion_Disks_Shakura_Sunyaev
  implicit none
  private
  public :: Accretion_Disks_Switched_Initialize

  ! Values used to indicate which type of accretion disk is being used.
  integer, parameter :: accretionDiskThin=0
  integer, parameter :: accretionDiskADAF=1

  ! Parameters controlling the range of accretion rates over which the accretion disk will be an ADAF.
  double precision :: accretionRateThinDiskMinimum,accretionRateThinDiskMaximum,accretionRateTransitionWidth
  double precision :: accretionRateThinDiskMinimumLogarithmic,accretionRateThinDiskMaximumLogarithmic
  logical          :: accretionRateThinDiskMinimumExists,accretionRateThinDiskMaximumExists

contains

  !# <accretionDisksMethod>
  !#  <unitName>Accretion_Disks_Switched_Initialize</unitName>
  !# </accretionDisksMethod>
  subroutine Accretion_Disks_Switched_Initialize(accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get&
       &,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: accretionDisksMethod
    procedure(double precision), pointer, intent(inout) :: Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get&
         &,Accretion_Disk_Jet_Power_Get
    character(len=30)                            :: accretionRateThin
    
    if (accretionDisksMethod == 'switched') then
       Accretion_Disk_Radiative_Efficiency_Get => Accretion_Disk_Radiative_Efficiency_Switched
       Black_Hole_Spin_Up_Rate_Get             => Black_Hole_Spin_Up_Rate_Switched
       Accretion_Disk_Jet_Power_Get            => Accretion_Disk_Jet_Power_Switched
       !@ <inputParameter>
       !@   <name>accretionRateThinDiskMinimum</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The accretion rate (in Eddington units) below which a switched accretion disk becomes an ADAF.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionRateThinDiskMinimum",accretionRateThin,defaultValue='0.01d0')
       if (trim(accretionRateThin) == "none") then
          accretionRateThinDiskMinimumExists=.false.
       else
          accretionRateThinDiskMinimumExists=.true.
          read (accretionRateThin,*) accretionRateThinDiskMinimum
          accretionRateThinDiskMinimumLogarithmic=dlog(accretionRateThinDiskMinimum)
       end if
       !@ <inputParameter>
       !@   <name>accretionRateThinDiskMaximum</name>
       !@   <defaultValue>0.3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The accretion rate (in Eddington units) above which a switched accretion disk becomes an ADAF.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionRateThinDiskMaximum",accretionRateThin,defaultValue="0.30d0")
       if (trim(accretionRateThin) == "none") then
          accretionRateThinDiskMaximumExists=.false.
       else
          accretionRateThinDiskMaximumExists=.true.
          read (accretionRateThin,*) accretionRateThinDiskMaximum
          accretionRateThinDiskMaximumLogarithmic=dlog(accretionRateThinDiskMaximum)
       end if
       !@ <inputParameter>
       !@   <name>accretionRateTransitionWidth</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The width (in $\ln[\dot{M}/\dot{M}_{\rm Eddington}]$) over which transitions between accretion disk states occur.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionRateTransitionWidth",accretionRateTransitionWidth,defaultValue=0.1d0)
    end if
    return
  end subroutine Accretion_Disks_Switched_Initialize

  double precision function Accretion_Disk_Radiative_Efficiency_Switched(thisNode,massAccretionRate)
    !% Computes the radiative efficiency for a switching accretion disk.
    use Black_Hole_Fundamentals
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision                         :: adafFraction

    adafFraction=Accretion_Disk_Switched_ADAF_Fraction(thisNode,massAccretionRate)
    Accretion_Disk_Radiative_Efficiency_Switched= (1.0d0-adafFraction)*Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev(thisNode,massAccretionRate) &
         &                                       +       adafFraction *Accretion_Disk_Radiative_Efficiency_ADAF           (thisNode,massAccretionRate)
    return
  end function Accretion_Disk_Radiative_Efficiency_Switched

  double precision function Accretion_Disk_Jet_Power_Switched(thisNode,massAccretionRate)
    !% Computes the jet power for a switching accretion disk.
    use Tree_Nodes
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision                         :: adafFraction

    adafFraction=Accretion_Disk_Switched_ADAF_Fraction(thisNode,massAccretionRate)
    Accretion_Disk_Jet_Power_Switched= (1.0d0-adafFraction)*Accretion_Disk_Jet_Power_Shakura_Sunyaev(thisNode,massAccretionRate) &
         &                            +       adafFraction *Accretion_Disk_Jet_Power_ADAF(thisNode,massAccretionRate)
    return
  end function Accretion_Disk_Jet_Power_Switched

  double precision function Black_Hole_Spin_Up_Rate_Switched(thisNode,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisNode} due to accretion from a switching accretion disk.
    !% disk.
    use Tree_Nodes
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision                         :: adafFraction

    adafFraction=Accretion_Disk_Switched_ADAF_Fraction(thisNode,massAccretionRate)

    Black_Hole_Spin_Up_Rate_Switched= (1.0d0-adafFraction)*Black_Hole_Spin_Up_Rate_Shakura_Sunyaev(thisNode,massAccretionRate) &
         &                           +       adafFraction* Black_Hole_Spin_Up_Rate_ADAF(thisNode,massAccretionRate)

    return
  end function Black_Hole_Spin_Up_Rate_Switched

  double precision function Accretion_Disk_Switched_ADAF_Fraction(thisNode,massAccretionRate)
    !% Decide which type of accretion disk to use.
    use Tree_Nodes
    use Black_Hole_Fundamentals
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: massAccretionRate
    double precision, parameter              :: exponentialArgumentMaximum=60.0d0
    double precision                         :: eddingtonAccretionRate,massAccretionRateDimensionless,adafFraction,accretionRateLogarithmic,argument

    ! Get the Eddington accretion rate.
    eddingtonAccretionRate=Black_Hole_Eddington_Accretion_Rate(thisNode)

    ! Check that a black hole is present.
    if (eddingtonAccretionRate > 0.0d0 .and. massAccretionRate > 0.0d0) then

       ! Compute the accretion rate in Eddington units.
       massAccretionRateDimensionless=massAccretionRate/eddingtonAccretionRate

       ! Compute the ADAF fraction.
       accretionRateLogarithmic=dlog(massAccretionRateDimensionless)
       adafFraction=0.0d0
       if (accretionRateThinDiskMinimumExists) then
          argument=min( (accretionRateLogarithmic-accretionRateThinDiskMinimumLogarithmic)/accretionRateTransitionWidth,exponentialArgumentMaximum)
          adafFraction=adafFraction+1.0d0/(1.0d0+dexp(argument))
       end if
       if (accretionRateThinDiskMaximumExists) then
          argument=min(-(accretionRateLogarithmic-accretionRateThinDiskMaximumLogarithmic)/accretionRateTransitionWidth,exponentialArgumentMaximum)
          adafFraction=adafFraction+1.0d0/(1.0d0+dexp(argument))
       end if
       Accretion_Disk_Switched_ADAF_Fraction=adafFraction

    else
       
       ! No black hole present: assume a thin disk.
       Accretion_Disk_Switched_ADAF_Fraction=0.0d0

    end if

    return
  end function Accretion_Disk_Switched_ADAF_Fraction

end module Accretion_Disks_Switched

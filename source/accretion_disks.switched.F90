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

  ! Option controlling ADAF radiative efficiency.
  logical          :: accretionDiskSwitchedScaleAdafRadiativeEfficiency

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
    type     (varying_string  ),          intent(in   ) :: accretionDisksMethod
    procedure(double precision), pointer, intent(inout) :: Accretion_Disk_Radiative_Efficiency_Get,Black_Hole_Spin_Up_Rate_Get &
         &,Accretion_Disk_Jet_Power_Get
    character(len=30          )                         :: accretionRateThin
    
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
        !@ <inputParameter>
       !@   <name>accretionDiskSwitchedScaleAdafRadiativeEfficiency</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Specifies whether the radiative efficiency of the ADAF component in a switched accretion disk scales with accretion rate.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionDiskSwitchedScaleAdafRadiativeEfficiency",accretionDiskSwitchedScaleAdafRadiativeEfficiency,defaultValue=.true.)
    end if
    return
  end subroutine Accretion_Disks_Switched_Initialize

  double precision function Accretion_Disk_Radiative_Efficiency_Switched(thisBlackHole,massAccretionRate)
    !% Computes the radiative efficiency for a switching accretion disk.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                                        :: adafFraction,thinDiskRadiativeEfficiency,adafRadiativeEfficiency

    adafFraction               =Accretion_Disk_Switched_ADAF_Fraction              (thisBlackHole,massAccretionRate)
    thinDiskRadiativeEfficiency=Accretion_Disk_Radiative_Efficiency_Shakura_Sunyaev(thisBlackHole,massAccretionRate)
    adafRadiativeEfficiency    =Accretion_Disk_Radiative_Efficiency_ADAF           (thisBlackHole,massAccretionRate)
    if (accretionDiskSwitchedScaleAdafRadiativeEfficiency) adafRadiativeEfficiency=adafRadiativeEfficiency&
         &*Accretion_Disk_Switched_ADAF_Radiative_Efficiency_Scaling(thisBlackHole,massAccretionRate)
    Accretion_Disk_Radiative_Efficiency_Switched= (1.0d0-adafFraction)*thinDiskRadiativeEfficiency &
         &                                       +       adafFraction *    adafRadiativeEfficiency
    return
  end function Accretion_Disk_Radiative_Efficiency_Switched

  double precision function Accretion_Disk_Jet_Power_Switched(thisBlackHole,massAccretionRate)
    !% Computes the jet power for a switching accretion disk.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                                        :: adafFraction

    adafFraction=Accretion_Disk_Switched_ADAF_Fraction(thisBlackHole,massAccretionRate)
    Accretion_Disk_Jet_Power_Switched= (1.0d0-adafFraction)*Accretion_Disk_Jet_Power_Shakura_Sunyaev(thisBlackHole,massAccretionRate) &
         &                            +       adafFraction *Accretion_Disk_Jet_Power_ADAF(thisBlackHole,massAccretionRate)
    return
  end function Accretion_Disk_Jet_Power_Switched

  double precision function Black_Hole_Spin_Up_Rate_Switched(thisBlackHole,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\tt thisBlackHole} due to accretion from a switching accretion disk.
    !% disk.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                                        :: adafFraction

    adafFraction=Accretion_Disk_Switched_ADAF_Fraction(thisBlackHole,massAccretionRate)

    Black_Hole_Spin_Up_Rate_Switched= (1.0d0-adafFraction)*Black_Hole_Spin_Up_Rate_Shakura_Sunyaev(thisBlackHole,massAccretionRate) &
         &                           +       adafFraction* Black_Hole_Spin_Up_Rate_ADAF(thisBlackHole,massAccretionRate)

    return
  end function Black_Hole_Spin_Up_Rate_Switched

  double precision function Accretion_Disk_Switched_ADAF_Fraction(thisBlackHole,massAccretionRate)
    !% Decide which type of accretion disk to use.
    use Galacticus_Nodes
    use Black_Hole_Fundamentals
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                        , parameter     :: exponentialArgumentMaximum=60.0d0
    double precision                                        :: eddingtonAccretionRate,massAccretionRateDimensionless,adafFraction,accretionRateLogarithmic,argument

    ! Get the Eddington accretion rate.
    eddingtonAccretionRate=Black_Hole_Eddington_Accretion_Rate(thisBlackHole)

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

  double precision function Accretion_Disk_Switched_ADAF_Radiative_Efficiency_Scaling(thisBlackHole,massAccretionRate)
    !% Determine the scaling of radiative efficiency of the ADAF component in a switched accretion disk.
    use Galacticus_Nodes
    use Black_Hole_Fundamentals
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                                        :: eddingtonAccretionRate,massAccretionRateDimensionless

    ! Get the Eddington accretion rate.
    eddingtonAccretionRate=Black_Hole_Eddington_Accretion_Rate(thisBlackHole)

    ! Check that a black hole is present.
    if (eddingtonAccretionRate > 0.0d0 .and. massAccretionRate > 0.0d0) then

       ! Compute the accretion rate in Eddington units.
       massAccretionRateDimensionless=massAccretionRate/eddingtonAccretionRate

       ! If below the critical accretion rate for transition to a thin disk, reduce the radiative efficiency by a factor
       ! proportional to the accretion rate.
       if (accretionRateThinDiskMinimumExists.and.massAccretionRateDimensionless < accretionRateThinDiskMinimum) then
          Accretion_Disk_Switched_ADAF_Radiative_Efficiency_Scaling=massAccretionRateDimensionless/accretionRateThinDiskMinimum
       else
          Accretion_Disk_Switched_ADAF_Radiative_Efficiency_Scaling=1.0d0
       end if

    else
       
       ! No black hole present: return unit scaling.
       Accretion_Disk_Switched_ADAF_Radiative_Efficiency_Scaling=1.0d0

    end if

    return
  end function Accretion_Disk_Switched_ADAF_Radiative_Efficiency_Scaling
  
end module Accretion_Disks_Switched

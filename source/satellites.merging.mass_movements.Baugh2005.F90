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

!% Contains a module which implements the \cite{baugh_can_2005} model of mass movements during satellite mergers.

module Satellite_Merging_Mass_Movements_Baugh2005
  !% Implements the \cite{baugh_can_2005} model of mass movements during satellite mergers.
  use Satellite_Merging_Mass_Movements_Descriptors
  implicit none
  private
  public :: Satellite_Merging_Mass_Movements_Baugh2005_Initialize

  ! Mass ratio above which a merger is considered to be "major".
  double precision :: majorMergerMassRatio,burstMassRatio,burstCriticalGasFraction

  ! Location to which gas from satellite galaxy in minor merger is moved.
  integer          :: minorMergerGasMovesToValue

contains

  !# <satelliteMergingMassMovementsMethod>
  !#  <unitName>Satellite_Merging_Mass_Movements_Baugh2005_Initialize</unitName>
  !# </satelliteMergingMassMovementsMethod>
  subroutine Satellite_Merging_Mass_Movements_Baugh2005_Initialize(satelliteMergingMassMovementsMethod,Satellite_Merging_Mass_Movement_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingMassMovementsMethod
    procedure(),          pointer, intent(inout) :: Satellite_Merging_Mass_Movement_Get
    character(len=10)                            :: minorMergerGasMovesTo

    if (satelliteMergingMassMovementsMethod == 'Baugh2005') then
       Satellite_Merging_Mass_Movement_Get => Satellite_Merging_Mass_Movement_Baugh2005
       !@ <inputParameter>
       !@   <name>majorMergerMassRatio</name>
       !@   <defaultValue>0.3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass ratio above which mergers are considered to be ``major'' in the \cite{baugh_can_2005} merger mass movements method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("majorMergerMassRatio"    ,majorMergerMassRatio    ,defaultValue=0.30d0    )
       !@ <inputParameter>
       !@   <name>burstMassRatio</name>
       !@   <defaultValue>0.05</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass ratio above which mergers are considered to trigger a burst in the \cite{baugh_can_2005} merger mass movements method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("burstMassRatio"          ,burstMassRatio          ,defaultValue=0.05d0    )
       !@ <inputParameter>
       !@   <name>burstCriticalGasFraction</name>
       !@   <defaultValue>0.75</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The host gas fraction above which mergers are considered to trigger a burst in the \cite{baugh_can_2005} merger mass movements method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("burstCriticalGasFraction",burstCriticalGasFraction,defaultValue=0.75d0    )
       !@ <inputParameter>
       !@   <name>minorMergerGasMovesTo</name>
       !@   <defaultValue>spheroid</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The component to which satellite galaxy gas moves to as a result of a minor merger.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("minorMergerGasMovesTo"   ,minorMergerGasMovesTo   ,defaultValue="spheroid")
       select case (trim(minorMergerGasMovesTo))
       case ("disk")
          minorMergerGasMovesToValue=movesToDisk
       case ("spheroid")
          minorMergerGasMovesToValue=movesToSpheroid
       case default
          call Galacticus_Error_Report('Satellite_Merging_Mass_Movements_Baugh2005_Initialize','unrecognized location for minor merger satellite gas')
       end select
    end if
    return
  end subroutine Satellite_Merging_Mass_Movements_Baugh2005_Initialize

  subroutine Satellite_Merging_Mass_Movement_Baugh2005(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo,mergerIsMajor)
    !% Determine how different mass components should be redistributed as the result of a merger according to the model of \cite{baugh_can_2005}.
    use Galacticus_Nodes
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(out)             :: gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo
    logical,        intent(out)             :: mergerIsMajor
    type(treeNode),                pointer  :: hostNode
    logical                                 :: triggersBurst
    double precision                        :: satelliteMass,hostMass,hostGasMass,hostSpheroidMass

    ! Find the node to merge with.
   hostNode => thisNode%mergesWith()

    ! Find the baryonic masses of the two galaxies.
    satelliteMass=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic)
    hostMass        =Galactic_Structure_Enclosed_Mass(hostNode,massType=massTypeGalactic                                    )
    hostGasMass     =Galactic_Structure_Enclosed_Mass(hostNode,massType=massTypeGaseous                                     )
    hostSpheroidMass=Galactic_Structure_Enclosed_Mass(hostNode,componentType=componentTypeSpheroid,massType=massTypeGalactic)
    ! Decide if the mass ratio is large enough to trigger a major merger.
    mergerIsMajor=satelliteMass >= majorMergerMassRatio*hostMass

    ! Determine if the merger will trigger a burst, in which case gas will be moved to the spheroid.
    triggersBurst=mergerIsMajor                                          &
         &         .or.                                                  &
         &        (                                                      &
         &         hostSpheroidMass <  burstMassRatio          *hostMass &
         &          .and.                                                &
         &         hostGasMass      >= burstCriticalGasFraction*hostMass &
         &        )

    if (mergerIsMajor) then
       gasMovesTo     =movesToSpheroid
       starsMoveTo    =movesToSpheroid
       hostGasMovesTo =movesToSpheroid
       hostStarsMoveTo=movesToSpheroid
    else
       if (triggersBurst) then
          gasMovesTo     =movesToSpheroid
          starsMoveTo    =movesToSpheroid
          hostGasMovesTo =movesToSpheroid
          hostStarsMoveTo=doesNotMove
       else
          gasMovesTo     =minorMergerGasMovesToValue
          starsMoveTo    =movesToSpheroid
          hostGasMovesTo =doesNotMove
          hostStarsMoveTo=doesNotMove
       end if
    end if

    return
  end subroutine Satellite_Merging_Mass_Movement_Baugh2005

end module Satellite_Merging_Mass_Movements_Baugh2005

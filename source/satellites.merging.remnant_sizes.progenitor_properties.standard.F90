!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements standard calculations of progenitor properties for merger remnant calculations.

module Satellite_Merging_Remnant_Progenitors_Properties_Standard
  !% Implements standard calculations of progenitor properties for merger remnant calculations.
  implicit none
  private
  public :: Satellite_Merging_Remnant_Progenitor_Properties_Standard_Init

contains

  !# <satelliteMergingRemnantProgenitorPropertiesMethod>
  !#  <unitName>Satellite_Merging_Remnant_Progenitor_Properties_Standard_Init</unitName>
  !# </satelliteMergingRemnantProgenitorPropertiesMethod>
  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Standard_Init(satelliteMergingRemnantProgenitorPropertiesMethod&
       &,Satellite_Merging_Remnant_Progenitor_Properties_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingRemnantProgenitorPropertiesMethod
    procedure(),          pointer, intent(inout) :: Satellite_Merging_Remnant_Progenitor_Properties_Get
    
    if (satelliteMergingRemnantProgenitorPropertiesMethod == 'standard') Satellite_Merging_Remnant_Progenitor_Properties_Get => Satellite_Merging_Remnant_Progenitor_Properties_Standard
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Standard_Init

  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Standard(satelliteNode,hostNode,satelliteMass,hostMass,satelliteSpheroidMass&
       &,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass&
         &,remnantSpheroidGasMass)
    !% Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    use Tree_Nodes
    use Galactic_Structure_Radii
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Satellite_Merging_Mass_Movements_Descriptors
    use Numerical_Constants_Physical
    use Galacticus_Error
    implicit none
    type(treeNode),   intent(inout), pointer :: satelliteNode,hostNode
    double precision, intent(out)            :: satelliteMass,hostMass,satelliteSpheroidMass,hostSpheroidMass,hostSpheroidMassPreMerger&
         &,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass,remnantSpheroidGasMass
    double precision                         :: componentMass,hostSpheroidDarkMatterFactor,hostDiskDarkMatterFactor&
         &,satelliteSpheroidDarkMatterFactor,satelliteDiskDarkMatterFactor,hostSpheroidHalfMassRadius,hostDiskHalfMassRadius&
         &,satelliteSpheroidHalfMassRadius,satelliteDiskHalfMassRadius

    ! Solve for the radii of the host and satellite nodes, to ensure they are computed and up to date.
    call Galactic_Structure_Radii_Solve(hostNode     )
    call Galactic_Structure_Radii_Solve(satelliteNode)

    ! Find the baryonic masses of the two galaxies.
    satelliteMass=Galactic_Structure_Enclosed_Mass(satelliteNode,massType=massTypeGalactic)
    hostMass     =Galactic_Structure_Enclosed_Mass(hostNode     ,massType=massTypeGalactic)

    ! Compute dark matter factors. These are the specific angular momenta of components divided by sqrt(G M r) where M is the
    ! component mass and r its half-mass radius. We use a weighted average of these factors to infer the specific angular momentum
    ! of the remnant from its mass and radius.
    componentMass                 =Tree_Node_Spheroid_Stellar_Mass     (hostNode     )+Tree_Node_Spheroid_Gas_Mass(hostNode     )
    hostSpheroidHalfMassRadius    =Tree_Node_Spheroid_Half_Mass_Radius (hostNode     )
    if (hostSpheroidHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       hostSpheroidDarkMatterFactor     = Tree_Node_Spheroid_Angular_Momentum(hostNode     )                     &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /dsqrt(gravitationalConstantGalacticus*hostSpheroidHalfMassRadius     )
    else
       hostSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass                  =Tree_Node_Disk_Stellar_Mass        (hostNode     )+Tree_Node_Disk_Gas_Mass    (hostNode     )
    hostDiskHalfMassRadius         =Tree_Node_Disk_Half_Mass_Radius    (hostNode     )
    if (hostDiskHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       hostDiskDarkMatterFactor         = Tree_Node_Disk_Angular_Momentum    (hostNode     )                     &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /dsqrt(gravitationalConstantGalacticus*hostDiskHalfMassRadius         )
    else
       hostDiskDarkMatterFactor=0.0d0
    end if
    componentMass                  =Tree_Node_Spheroid_Stellar_Mass    (satelliteNode)+Tree_Node_Spheroid_Gas_Mass(satelliteNode)
    satelliteSpheroidHalfMassRadius=Tree_Node_Spheroid_Half_Mass_Radius(satelliteNode)
    if (satelliteSpheroidHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteSpheroidDarkMatterFactor= Tree_Node_Spheroid_Angular_Momentum(satelliteNode)                     &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /dsqrt(gravitationalConstantGalacticus*satelliteSpheroidHalfMassRadius)
    else
       satelliteSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass                  =Tree_Node_Disk_Stellar_Mass        (satelliteNode)+Tree_Node_Disk_Gas_Mass    (satelliteNode)
    satelliteDiskHalfMassRadius    =Tree_Node_Disk_Half_Mass_Radius    (satelliteNode)
    if (satelliteDiskHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteDiskDarkMatterFactor    = Tree_Node_Disk_Angular_Momentum    (satelliteNode)                     &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /dsqrt(gravitationalConstantGalacticus*satelliteDiskHalfMassRadius    )
    else
       satelliteDiskDarkMatterFactor=0.0d0
    end if

    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (thisHostGasMovesTo)
    case (movesToSpheroid)
       hostSpheroidMass      =Tree_Node_Spheroid_Gas_Mass(hostNode)                             +Tree_Node_Disk_Gas_Mass(hostNode)
       hostRadius            =Tree_Node_Spheroid_Gas_Mass(hostNode)*hostSpheroidHalfMassRadius  +Tree_Node_Disk_Gas_Mass(hostNode)*hostDiskHalfMassRadius
       angularMomentumFactor =Tree_Node_Spheroid_Gas_Mass(hostNode)*hostSpheroidDarkMatterFactor+Tree_Node_Disk_Gas_Mass(hostNode)*hostDiskDarkMatterFactor
       remnantSpheroidGasMass=Tree_Node_Spheroid_Gas_Mass(hostNode)                             +Tree_Node_Disk_Gas_Mass(hostNode)
       remnantSpheroidMass   =Tree_Node_Spheroid_Gas_Mass(hostNode)                             +Tree_Node_Disk_Gas_Mass(hostNode)
    case (movesToDisk)
       hostSpheroidMass      =0.0d0
       hostRadius            =0.0d0
       angularMomentumFactor =0.0d0
       remnantSpheroidGasMass=0.0d0
       remnantSpheroidMass   =0.0d0
    case (doesNotMove)
       hostSpheroidMass      =Tree_Node_Spheroid_Gas_Mass(hostNode)
       hostRadius            =Tree_Node_Spheroid_Gas_Mass(hostNode)*hostSpheroidHalfMassRadius
       angularMomentumFactor      =Tree_Node_Spheroid_Gas_Mass(hostNode)*hostSpheroidDarkMatterFactor
       remnantSpheroidGasMass=Tree_Node_Spheroid_Gas_Mass(hostNode)
       remnantSpheroidMass   =Tree_Node_Spheroid_Gas_Mass(hostNode)
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisHostStarsMoveTo)
    case (movesToSpheroid)
       hostSpheroidMass     =hostSpheroidMass   +Tree_Node_Spheroid_Stellar_Mass(hostNode)                             +Tree_Node_Disk_Stellar_Mass(hostNode)
       hostRadius           =hostRadius         +Tree_Node_Spheroid_Stellar_Mass(hostNode)*hostSpheroidHalfMassRadius  +Tree_Node_Disk_Stellar_Mass(hostNode)*hostDiskHalfMassRadius
       angularMomentumFactor=angularMomentumFactor   +Tree_Node_Spheroid_Stellar_Mass(hostNode)*hostSpheroidDarkMatterFactor+Tree_Node_Disk_Stellar_Mass(hostNode)*hostDiskDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(hostNode)                             +Tree_Node_Disk_Stellar_Mass(hostNode)
    case (movesToDisk)
       hostSpheroidMass     =hostSpheroidMass
       hostRadius           =hostRadius
       angularMomentumFactor=angularMomentumFactor
    case (doesNotMove)
       hostSpheroidMass   =hostSpheroidMass   +Tree_Node_Spheroid_Stellar_Mass(hostNode)
       hostRadius         =hostRadius         +Tree_Node_Spheroid_Stellar_Mass(hostNode)*hostSpheroidHalfMassRadius
       angularMomentumFactor   =angularMomentumFactor   +Tree_Node_Spheroid_Stellar_Mass(hostNode)*hostSpheroidDarkMatterFactor
       remnantSpheroidMass=remnantSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(hostNode)
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisMergerGasMovesTo)
    case (movesToSpheroid)
       satelliteSpheroidMass =                       Tree_Node_Spheroid_Gas_Mass(satelliteNode)                                  +Tree_Node_Disk_Gas_Mass(satelliteNode)
       satelliteRadius       =                       Tree_Node_Spheroid_Gas_Mass(satelliteNode)*satelliteSpheroidHalfMassRadius  +Tree_Node_Disk_Gas_Mass(satelliteNode)*satelliteDiskHalfMassRadius
       angularMomentumFactor      =angularMomentumFactor      +Tree_Node_Spheroid_Gas_Mass(satelliteNode)*satelliteSpheroidDarkMatterFactor+Tree_Node_Disk_Gas_Mass(satelliteNode)*satelliteDiskDarkMatterFactor
       remnantSpheroidGasMass=remnantSpheroidGasMass+Tree_Node_Spheroid_Gas_Mass(satelliteNode)                                  +Tree_Node_Disk_Gas_Mass(satelliteNode)
       remnantSpheroidMass   =remnantSpheroidMass   +Tree_Node_Spheroid_Gas_Mass(satelliteNode)                                  +Tree_Node_Disk_Gas_Mass(satelliteNode)
    case (movesToDisk)
       satelliteSpheroidMass =0.0d0
       satelliteRadius       =0.0d0
       angularMomentumFactor      =angularMomentumFactor
    case (doesNotMove)
       satelliteSpheroidMass=                        Tree_Node_Spheroid_Gas_Mass(satelliteNode)
       satelliteRadius       =                       Tree_Node_Spheroid_Gas_Mass(satelliteNode)*satelliteSpheroidHalfMassRadius
       angularMomentumFactor      =angularMomentumFactor      +Tree_Node_Spheroid_Gas_Mass(satelliteNode)*satelliteSpheroidDarkMatterFactor
       remnantSpheroidGasMass=remnantSpheroidGasMass+Tree_Node_Spheroid_Gas_Mass(satelliteNode)
       remnantSpheroidMass   =remnantSpheroidMass   +Tree_Node_Spheroid_Gas_Mass(satelliteNode)
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisMergerStarsMoveTo)
    case (movesToSpheroid)
       satelliteSpheroidMass=satelliteSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(satelliteNode)                                  +Tree_Node_Disk_Stellar_Mass(satelliteNode)
       satelliteRadius      =satelliteRadius      +Tree_Node_Spheroid_Stellar_Mass(satelliteNode)*satelliteSpheroidHalfMassRadius  +Tree_Node_Disk_Stellar_Mass(satelliteNode)*satelliteDiskHalfMassRadius
       angularMomentumFactor     =angularMomentumFactor     +Tree_Node_Spheroid_Stellar_Mass(satelliteNode)*satelliteSpheroidDarkMatterFactor+Tree_Node_Disk_Stellar_Mass(satelliteNode)*satelliteDiskDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +Tree_Node_Spheroid_Stellar_Mass(satelliteNode)                                  +Tree_Node_Disk_Stellar_Mass(satelliteNode)
    case (movesToDisk)
       satelliteSpheroidMass=satelliteSpheroidMass
       satelliteRadius      =satelliteRadius
       angularMomentumFactor     =angularMomentumFactor
    case (doesNotMove)
       satelliteSpheroidMass=satelliteSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(satelliteNode)
       satelliteRadius      =satelliteRadius      +Tree_Node_Spheroid_Stellar_Mass(satelliteNode)*satelliteSpheroidHalfMassRadius
       angularMomentumFactor     =angularMomentumFactor     +Tree_Node_Spheroid_Stellar_Mass(satelliteNode)*satelliteSpheroidDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +Tree_Node_Spheroid_Stellar_Mass(satelliteNode)
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select

    ! Compute the angular momentum factor.
    if (satelliteSpheroidMass+hostSpheroidMass > 0.0d0) then
       angularMomentumFactor=angularMomentumFactor/(satelliteSpheroidMass+hostSpheroidMass)
    else
       angularMomentumFactor=1.0d0
    end if

    ! Trap cases where radius is zero, but mass is finite (due to numerical inaccuracies).
    if (hostRadius      <= 0.0d0) hostSpheroidMass     =0.0d0
    if (satelliteRadius <= 0.0d0) satelliteSpheroidMass=0.0d0
    
    ! Compute the radii of the spheroid components.
    if (hostSpheroidMass      > 0.0d0) hostRadius     =hostRadius     /hostSpheroidMass
    if (satelliteSpheroidMass > 0.0d0) satelliteRadius=satelliteRadius/satelliteSpheroidMass

    ! Compute the mass of the host spheroid before the merger.
    hostSpheroidMassPreMerger=Tree_Node_Spheroid_Stellar_Mass(hostNode)+Tree_Node_Spheroid_Gas_Mass(hostNode)
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Standard

end module Satellite_Merging_Remnant_Progenitors_Properties_Standard

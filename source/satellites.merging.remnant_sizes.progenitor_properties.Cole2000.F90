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


!% Contains a module which implements calculations of progenitor properties for merger remnant calculations using the algorithm of
!% \cite{cole_hierarchical_2000}.

module Satellite_Merging_Remnant_Progenitors_Properties_Cole2000
  !% Implements calculations of progenitor properties for merger remnant calculations using the algorithm of
  !% \cite{cole_hierarchical_2000}.
  use Tree_Nodes
  implicit none
  private
  public :: Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init

  ! Module global variables used in root finding.
  type(treeNode), pointer :: activeNode
  integer                 :: activeGasMovesTo,activeStarsMoveTo
  double precision        :: activeHalfMass
  !$omp threadprivate(activeNode,activeGasMovesTo,activeStarsMoveTo,activeHalfMass)

contains

  !# <satelliteMergingRemnantProgenitorPropertiesMethod>
  !#  <unitName>Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init</unitName>
  !# </satelliteMergingRemnantProgenitorPropertiesMethod>
  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init(satelliteMergingRemnantProgenitorPropertiesMethod&
       &,Satellite_Merging_Remnant_Progenitor_Properties_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingRemnantProgenitorPropertiesMethod
    procedure(),          pointer, intent(inout) :: Satellite_Merging_Remnant_Progenitor_Properties_Get
    
    if (satelliteMergingRemnantProgenitorPropertiesMethod == 'Cole2000') Satellite_Merging_Remnant_Progenitor_Properties_Get =>&
         & Satellite_Merging_Remnant_Progenitor_Properties_Cole2000
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init

  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Cole2000(satelliteNode,hostNode,satelliteMass,hostMass,satelliteSpheroidMass&
       &,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass&
         &,remnantSpheroidGasMass)
    !% Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    use Galactic_Structure_Radii
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Satellite_Merging_Mass_Movements_Descriptors
    use Numerical_Constants_Physical
    use Galacticus_Error
    use Root_Finder
    use FGSL
    implicit none
    type(treeNode),          intent(inout), pointer  :: satelliteNode,hostNode
    double precision,        intent(out)             :: satelliteMass,hostMass,satelliteSpheroidMass,hostSpheroidMass,hostSpheroidMassPreMerger&
         &,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass,remnantSpheroidGasMass
    type(fgsl_function),     save                    :: rootFunction
    type(fgsl_root_fsolver), save                    :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision                                 :: radiusMinimum,radiusMaximum,satelliteSpheroidHalfMassRadius&
         &,hostSpheroidHalfMassRadius,satelliteDiskHalfMassRadius,hostDiskHalfMassRadius,satelliteSpheroidDarkMatterFactor&
         &,hostSpheroidDarkMatterFactor,satelliteDiskDarkMatterFactor,hostDiskDarkMatterFactor,componentMass
    type(c_ptr)                                      :: parameterPointer

    ! Solve for the radii of the host and satellite nodes, to ensure they are computed and up to date.
    call Galactic_Structure_Radii_Solve(hostNode     )
    call Galactic_Structure_Radii_Solve(satelliteNode)

    ! Find the baryonic masses of the two galaxies.
    satelliteMass=Galactic_Structure_Enclosed_Mass(satelliteNode,massType=massTypeGalactic)
    hostMass     =Galactic_Structure_Enclosed_Mass(hostNode     ,massType=massTypeGalactic)

    ! Compute dark matter factors. These are the specific angular momenta of components divided by sqrt(G M r) where M is the
    ! component mass and r its half-mass radius. We use a weighted average of these factors to infer the specific angular momentum
    ! of the remnant from its mass and radius.
    componentMass=Tree_Node_Spheroid_Stellar_Mass(hostNode)+Tree_Node_Spheroid_Gas_Mass(hostNode)
    hostSpheroidHalfMassRadius=Tree_Node_Spheroid_Half_Mass_Radius(hostNode)
    if (hostSpheroidHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       hostSpheroidDarkMatterFactor=Tree_Node_Spheroid_Angular_Momentum(hostNode)/(componentMass**1.5d0) &
            &/dsqrt(gravitationalConstantGalacticus*hostSpheroidHalfMassRadius)
    else
       hostSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass=Tree_Node_Disk_Stellar_Mass(hostNode) +Tree_Node_Disk_Gas_Mass(hostNode)
    hostDiskHalfMassRadius=Tree_Node_Disk_Half_Mass_Radius(hostNode)
    if (hostDiskHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       hostDiskDarkMatterFactor=Tree_Node_Disk_Angular_Momentum(hostNode)/(componentMass**1.5d0)&
            &/dsqrt(gravitationalConstantGalacticus*hostDiskHalfMassRadius)
    else
       hostDiskDarkMatterFactor=0.0d0
    end if
    componentMass=Tree_Node_Spheroid_Stellar_Mass(satelliteNode) +Tree_Node_Spheroid_Gas_Mass(satelliteNode)
    satelliteSpheroidHalfMassRadius=Tree_Node_Spheroid_Half_Mass_Radius(satelliteNode)
    if (satelliteSpheroidHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteSpheroidDarkMatterFactor=Tree_Node_Spheroid_Angular_Momentum(satelliteNode)/(componentMass**1.5d0)&
            &/dsqrt(gravitationalConstantGalacticus*satelliteSpheroidHalfMassRadius)
    else
       satelliteSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass=Tree_Node_Disk_Stellar_Mass(satelliteNode) +Tree_Node_Disk_Gas_Mass(satelliteNode)
    satelliteDiskHalfMassRadius=Tree_Node_Disk_Half_Mass_Radius(satelliteNode)
    if (satelliteDiskHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteDiskDarkMatterFactor=Tree_Node_Disk_Angular_Momentum(satelliteNode)/(componentMass**1.5d0)&
            &/dsqrt(gravitationalConstantGalacticus*satelliteDiskHalfMassRadius)
    else
       satelliteDiskDarkMatterFactor=0.0d0
    end if

    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (thisHostGasMovesTo)
    case (movesToSpheroid)
       hostSpheroidMass      =Tree_Node_Spheroid_Gas_Mass(hostNode)+Tree_Node_Disk_Gas_Mass(hostNode)
       angularMomentumFactor =Tree_Node_Spheroid_Gas_Mass(hostNode)*hostSpheroidDarkMatterFactor+Tree_Node_Disk_Gas_Mass(hostNode)*hostDiskDarkMatterFactor
       remnantSpheroidGasMass=Tree_Node_Spheroid_Gas_Mass(hostNode)+Tree_Node_Disk_Gas_Mass(hostNode)
       remnantSpheroidMass   =Tree_Node_Spheroid_Gas_Mass(hostNode)+Tree_Node_Disk_Gas_Mass(hostNode)
    case (movesToDisk)
       hostSpheroidMass      =0.0d0
       angularMomentumFactor      =0.0d0
       remnantSpheroidGasMass=0.0d0
       remnantSpheroidMass   =0.0d0
    case (doesNotMove)
       hostSpheroidMass      =Tree_Node_Spheroid_Gas_Mass(hostNode)
       angularMomentumFactor =Tree_Node_Spheroid_Gas_Mass(hostNode)*hostSpheroidDarkMatterFactor
       remnantSpheroidGasMass=Tree_Node_Spheroid_Gas_Mass(hostNode)
       remnantSpheroidMass   =Tree_Node_Spheroid_Gas_Mass(hostNode)
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisHostStarsMoveTo)
    case (movesToSpheroid)
       hostSpheroidMass     =hostSpheroidMass   +Tree_Node_Spheroid_Stellar_Mass(hostNode)+Tree_Node_Disk_Stellar_Mass(hostNode)
       angularMomentumFactor=angularMomentumFactor   +Tree_Node_Spheroid_Stellar_Mass(hostNode)*hostSpheroidDarkMatterFactor+Tree_Node_Disk_Stellar_Mass(hostNode)*hostDiskDarkMatterFactor
      remnantSpheroidMass   =remnantSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(hostNode)+Tree_Node_Disk_Stellar_Mass(hostNode)
    case (movesToDisk)
       hostSpheroidMass     =hostSpheroidMass
       angularMomentumFactor=angularMomentumFactor
    case (doesNotMove)
       hostSpheroidMass     =hostSpheroidMass   +Tree_Node_Spheroid_Stellar_Mass(hostNode)
       angularMomentumFactor=angularMomentumFactor   +Tree_Node_Spheroid_Stellar_Mass(hostNode)*hostSpheroidDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(hostNode)
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisMergerGasMovesTo)
    case (movesToSpheroid)
       satelliteSpheroidMass =                       Tree_Node_Spheroid_Gas_Mass(satelliteNode)+Tree_Node_Disk_Gas_Mass(satelliteNode)
       angularMomentumFactor =angularMomentumFactor +Tree_Node_Spheroid_Gas_Mass(satelliteNode)*satelliteSpheroidDarkMatterFactor+Tree_Node_Disk_Gas_Mass(satelliteNode)*satelliteDiskDarkMatterFactor
       remnantSpheroidGasMass=remnantSpheroidGasMass+Tree_Node_Spheroid_Gas_Mass(satelliteNode)+Tree_Node_Disk_Gas_Mass(satelliteNode)
       remnantSpheroidMass   =remnantSpheroidMass   +Tree_Node_Spheroid_Gas_Mass(satelliteNode)+Tree_Node_Disk_Gas_Mass(satelliteNode)
    case (movesToDisk)
       satelliteSpheroidMass =0.0d0
       angularMomentumFactor =angularMomentumFactor
    case (doesNotMove)
       satelliteSpheroidMass=                        Tree_Node_Spheroid_Gas_Mass(satelliteNode)
       angularMomentumFactor =angularMomentumFactor +Tree_Node_Spheroid_Gas_Mass(satelliteNode)*satelliteSpheroidDarkMatterFactor
       remnantSpheroidGasMass=remnantSpheroidGasMass+Tree_Node_Spheroid_Gas_Mass(satelliteNode)
       remnantSpheroidMass   =remnantSpheroidMass   +Tree_Node_Spheroid_Gas_Mass(satelliteNode)
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisMergerStarsMoveTo)
    case (movesToSpheroid)
       satelliteSpheroidMass=satelliteSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(satelliteNode)+Tree_Node_Disk_Stellar_Mass(satelliteNode)
       angularMomentumFactor=angularMomentumFactor+Tree_Node_Spheroid_Stellar_Mass(satelliteNode)*satelliteSpheroidDarkMatterFactor+Tree_Node_Disk_Stellar_Mass(satelliteNode)*satelliteDiskDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +Tree_Node_Spheroid_Stellar_Mass(satelliteNode)+Tree_Node_Disk_Stellar_Mass(satelliteNode)
    case (movesToDisk)
       satelliteSpheroidMass=satelliteSpheroidMass
       angularMomentumFactor=angularMomentumFactor
    case (doesNotMove)
       satelliteSpheroidMass=satelliteSpheroidMass+Tree_Node_Spheroid_Stellar_Mass(satelliteNode)
       remnantSpheroidMass  =remnantSpheroidMass  +Tree_Node_Spheroid_Stellar_Mass(satelliteNode)
       angularMomentumFactor=angularMomentumFactor+Tree_Node_Spheroid_Stellar_Mass(satelliteNode)*satelliteSpheroidDarkMatterFactor
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    
    ! Compute the half-mass radii of the material that will end up in the remnant spheroid.
    ! Host node.
    if (hostSpheroidMass > 0.0d0) then
       activeNode       => hostNode
       activeGasMovesTo =  thisHostGasMovesTo
       activeStarsMoveTo=  thisHostStarsMoveTo
       activeHalfMass   =  0.0d0
       activeHalfMass   =  0.5d0*Half_Mass_Radius_Root_Cole2000(radiusLarge,parameterPointer)
       radiusMinimum    =  Galactic_Structure_Radius_Enclosing_Mass(activeNode,fractionalMass=0.50d0,massType=massTypeGalactic)
       radiusMaximum    =  radiusMinimum
       do while (Half_Mass_Radius_Root_Cole2000(radiusMinimum,parameterPointer) >= 0.0d0)
          radiusMinimum=0.5d0*radiusMinimum
       end do
       do while (Half_Mass_Radius_Root_Cole2000(radiusMaximum,parameterPointer) <= 0.0d0)
          radiusMaximum=2.0d0*radiusMaximum
       end do
       hostRadius        =  Root_Find(radiusMinimum,radiusMaximum,Half_Mass_Radius_Root_Cole2000,parameterPointer,rootFunction&
            &,rootFunctionSolver,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    else
       hostRadius=0.0d0
    end if
    if (satelliteSpheroidMass > 0.0d0) then
       activeNode       => satelliteNode
       activeGasMovesTo =  thisMergerGasMovesTo
       activeStarsMoveTo=  thisMergerStarsMoveTo
       activeHalfMass   =  0.0d0
       activeHalfMass   =  0.50d0*Half_Mass_Radius_Root_Cole2000(radiusLarge,parameterPointer)
       radiusMinimum    =  Galactic_Structure_Radius_Enclosing_Mass(activeNode,fractionalMass=0.50d0,massType=massTypeGalactic)
       radiusMaximum    =  radiusMinimum
       do while (Half_Mass_Radius_Root_Cole2000(radiusMinimum,parameterPointer) >= 0.0d0)
          radiusMinimum=0.5d0*radiusMinimum
       end do
       do while (Half_Mass_Radius_Root_Cole2000(radiusMaximum,parameterPointer) <= 0.0d0)
          radiusMaximum=2.0d0*radiusMaximum
       end do
       satelliteRadius   =  Root_Find(radiusMinimum,radiusMaximum,Half_Mass_Radius_Root_Cole2000,parameterPointer,rootFunction&
            &,rootFunctionSolver,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
    else
       satelliteRadius=0.0d0
    end if

    ! Compute the angular momentum factor.
    ! <v0.9.1> This is actually not the way that Cole et al. (2000) do this. Instead they directly compute the mass of dark matter
    ! inside the half mass radius of each galaxy. The angularMomentumFactor is then equal to a mass-weighted mean of (1+f_DM) for
    ! the host and satellite. This will be implemented once Galacticus has a generic function for getting the enclosed mass of
    ! dark matter in an adiabatically contracted halo.
    if (satelliteSpheroidMass+hostSpheroidMass > 0.0d0) then
       angularMomentumFactor=angularMomentumFactor/(satelliteSpheroidMass+hostSpheroidMass)
    else
       angularMomentumFactor=1.0d0
    end if

    ! Compute the mass of the host spheroid before the merger.
    hostSpheroidMassPreMerger=Tree_Node_Spheroid_Stellar_Mass(hostNode)+Tree_Node_Spheroid_Gas_Mass(hostNode)
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Cole2000

  function Half_Mass_Radius_Root_Cole2000(radius,parameterPointer) bind(c)
    !% Function used in root finding for progenitor galaxy half-mass radii.
    use, intrinsic :: ISO_C_Binding
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: Half_Mass_Radius_Root_Cole2000

    ! Initialize enclosed mass to negative of the half mass.
    Half_Mass_Radius_Root_Cole2000=-activeHalfMass

    ! Account for gas mass.
    select case (activeGasMovesTo)
    case (movesToSpheroid)
       Half_Mass_Radius_Root_Cole2000= Half_Mass_Radius_Root_Cole2000                                                                          &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeGaseous,componentType=componentTypeSpheroid) &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeGaseous,componentType=componentTypeDisk    )
    case (doesNotMove    )
       Half_Mass_Radius_Root_Cole2000= Half_Mass_Radius_Root_Cole2000                                                                          &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeGaseous,componentType=componentTypeSpheroid)
    end select
  
    ! Account for gas mass.
    select case (activeStarsMoveTo)
    case (movesToSpheroid)
       Half_Mass_Radius_Root_Cole2000= Half_Mass_Radius_Root_Cole2000                                                                          &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeStellar,componentType=componentTypeSpheroid) &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeStellar,componentType=componentTypeDisk    )
    case (doesNotMove    )
       Half_Mass_Radius_Root_Cole2000= Half_Mass_Radius_Root_Cole2000                                                                          &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,massType=massTypeStellar,componentType=componentTypeSpheroid)
    end select
    return
  end function Half_Mass_Radius_Root_Cole2000

end module Satellite_Merging_Remnant_Progenitors_Properties_Cole2000

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
    use Galacticus_Error
    use Galacticus_Nodes
    implicit none
    type(varying_string),          intent(in   ) :: satelliteMergingRemnantProgenitorPropertiesMethod
    procedure(Satellite_Merging_Remnant_Progenitor_Properties_Standard),          pointer, intent(inout) :: Satellite_Merging_Remnant_Progenitor_Properties_Get
    
    if (satelliteMergingRemnantProgenitorPropertiesMethod == 'standard') then
       Satellite_Merging_Remnant_Progenitor_Properties_Get => Satellite_Merging_Remnant_Progenitor_Properties_Standard
       ! Ensure that required methods are supported.
       if     (                                                                &
            &  .not.                                                           &
            &       (                                                          &
            &        defaultDiskComponent    %    massStellarIsGettable().and. &
            &        defaultDiskComponent    %        massGasIsGettable().and. &
            &        defaultDiskComponent    % halfMassRadiusIsGettable().and. &
            &        defaultDiskComponent    %angularMomentumIsGettable().and. &
            &        defaultSpheroidComponent%    massStellarIsGettable().and. &
            &        defaultSpheroidComponent%        massGasIsGettable().and. &
            &        defaultSpheroidComponent% halfMassRadiusIsGettable().and. &
            &        defaultSpheroidComponent%angularMomentumIsGettable()      &
            &  )                                                               &
            & ) call Galacticus_Error_Report('Satellite_Merging_Remnant_Progenitor_Properties_Standard_Init','this method requires that massStellar, massGas, halfMassRadius, and angularMomentum properties must all be gettable for both disk and spheroid')
    end if
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Standard_Init

  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Standard(satelliteNode,hostNode,satelliteMass,hostMass,satelliteSpheroidMass&
       &,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass &
       &,remnantSpheroidGasMass)
    !% Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    use Galacticus_Nodes
    use Galactic_Structure_Radii
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Satellite_Merging_Mass_Movements_Descriptors
    use Numerical_Constants_Physical
    use Galacticus_Error
    implicit none
    type (treeNode             ), intent(inout), pointer :: satelliteNode,hostNode
    double precision            , intent(  out)          :: satelliteMass,hostMass,satelliteSpheroidMass,hostSpheroidMass,hostSpheroidMassPreMerger&
         &,satelliteRadius,hostRadius,angularMomentumFactor,remnantSpheroidMass,remnantSpheroidGasMass
    class(nodeComponentDisk    ),                pointer :: hostDiskComponent    ,satelliteDiskComponent
    class(nodeComponentSpheroid),                pointer :: hostSpheroidComponent,satelliteSpheroidComponent
    double precision                                     :: componentMass,hostSpheroidDarkMatterFactor,hostDiskDarkMatterFactor &
         &,satelliteSpheroidDarkMatterFactor,satelliteDiskDarkMatterFactor,hostSpheroidHalfMassRadius,hostDiskHalfMassRadius &
         &,satelliteSpheroidHalfMassRadius,satelliteDiskHalfMassRadius

    ! Get the disk and spheroid components of host and satellite.
    hostDiskComponent          =>      hostNode%disk    ()
    hostSpheroidComponent      =>      hostNode%spheroid()
    satelliteDiskComponent     => satelliteNode%disk    ()
    satelliteSpheroidComponent => satelliteNode%spheroid()
    
    ! Solve for the radii of the host and satellite nodes, to ensure they are computed and up to date.
    call Galactic_Structure_Radii_Solve(hostNode     )
    call Galactic_Structure_Radii_Solve(satelliteNode)

    ! Find the baryonic masses of the two galaxies.
    satelliteMass=Galactic_Structure_Enclosed_Mass(satelliteNode,massType=massTypeGalactic)
    hostMass     =Galactic_Structure_Enclosed_Mass(hostNode     ,massType=massTypeGalactic)

    ! Compute dark matter factors. These are the specific angular momenta of components divided by sqrt(G M r) where M is the
    ! component mass and r its half-mass radius. We use a weighted average of these factors to infer the specific angular momentum
    ! of the remnant from its mass and radius.
    componentMass                 =hostSpheroidComponent%massStellar()+hostSpheroidComponent%massGas()
    hostSpheroidHalfMassRadius    =hostSpheroidComponent%halfMassRadius()
    if (hostSpheroidHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       hostSpheroidDarkMatterFactor     = hostSpheroidComponent%angularMomentum()                     &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*hostSpheroidHalfMassRadius     )
    else
       hostSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass                  =hostDiskComponent%massStellar()+hostDiskComponent%massGas()
    hostDiskHalfMassRadius         =hostDiskComponent%halfMassRadius()
    if (hostDiskHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       hostDiskDarkMatterFactor         = hostDiskComponent%angularMomentum()                     &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*hostDiskHalfMassRadius         )
    else
       hostDiskDarkMatterFactor=0.0d0
    end if
    componentMass                  =satelliteSpheroidComponent%massStellar()+satelliteSpheroidComponent%massGas()
    satelliteSpheroidHalfMassRadius=satelliteSpheroidComponent%halfMassRadius()
    if (satelliteSpheroidHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteSpheroidDarkMatterFactor= satelliteSpheroidComponent%angularMomentum()                     &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*satelliteSpheroidHalfMassRadius)
    else
       satelliteSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass                  =satelliteDiskComponent%massStellar()+satelliteDiskComponent%massGas()
    satelliteDiskHalfMassRadius    =satelliteDiskComponent%halfMassRadius()
    if (satelliteDiskHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteDiskDarkMatterFactor    = satelliteDiskComponent%angularMomentum()                     &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*satelliteDiskHalfMassRadius    )
    else
       satelliteDiskDarkMatterFactor=0.0d0
    end if

    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (thisHostGasMovesTo)
    case (movesToSpheroid)
       hostSpheroidMass      =hostSpheroidComponent%massGas()                             +hostDiskComponent%massGas()
       hostRadius            =hostSpheroidComponent%massGas()*hostSpheroidHalfMassRadius  +hostDiskComponent%massGas()*hostDiskHalfMassRadius
       angularMomentumFactor =hostSpheroidComponent%massGas()*hostSpheroidDarkMatterFactor+hostDiskComponent%massGas()*hostDiskDarkMatterFactor
       remnantSpheroidGasMass=hostSpheroidComponent%massGas()                             +hostDiskComponent%massGas()
       remnantSpheroidMass   =hostSpheroidComponent%massGas()                             +hostDiskComponent%massGas()
    case (movesToDisk)
       hostSpheroidMass      =0.0d0
       hostRadius            =0.0d0
       angularMomentumFactor =0.0d0
       remnantSpheroidGasMass=0.0d0
       remnantSpheroidMass   =0.0d0
    case (doesNotMove)
       hostSpheroidMass      =hostSpheroidComponent%massGas()
       hostRadius            =hostSpheroidComponent%massGas()*hostSpheroidHalfMassRadius
       angularMomentumFactor =hostSpheroidComponent%massGas()*hostSpheroidDarkMatterFactor
       remnantSpheroidGasMass=hostSpheroidComponent%massGas()
       remnantSpheroidMass   =hostSpheroidComponent%massGas()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisHostStarsMoveTo)
    case (movesToSpheroid)
       hostSpheroidMass     =hostSpheroidMass     +hostSpheroidComponent%massStellar()                             +hostDiskComponent%massStellar()
       hostRadius           =hostRadius           +hostSpheroidComponent%massStellar()*hostSpheroidHalfMassRadius  +hostDiskComponent%massStellar()*hostDiskHalfMassRadius
       angularMomentumFactor=angularMomentumFactor+hostSpheroidComponent%massStellar()*hostSpheroidDarkMatterFactor+hostDiskComponent%massStellar()*hostDiskDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +hostSpheroidComponent%massStellar()                             +hostDiskComponent%massStellar()
    case (movesToDisk)
       hostSpheroidMass     =hostSpheroidMass
       hostRadius           =hostRadius
       angularMomentumFactor=angularMomentumFactor
    case (doesNotMove)
       hostSpheroidMass     =hostSpheroidMass     +hostSpheroidComponent%massStellar()
       hostRadius           =hostRadius           +hostSpheroidComponent%massStellar()*hostSpheroidHalfMassRadius
       angularMomentumFactor=angularMomentumFactor+hostSpheroidComponent%massStellar()*hostSpheroidDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +hostSpheroidComponent%massStellar()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisMergerGasMovesTo)
    case (movesToSpheroid)
       satelliteSpheroidMass =                       satelliteSpheroidComponent%massGas()                                  +satelliteDiskComponent%massGas()
       satelliteRadius       =                       satelliteSpheroidComponent%massGas()*satelliteSpheroidHalfMassRadius  +satelliteDiskComponent%massGas()*satelliteDiskHalfMassRadius
       angularMomentumFactor =angularMomentumFactor +satelliteSpheroidComponent%massGas()*satelliteSpheroidDarkMatterFactor+satelliteDiskComponent%massGas()*satelliteDiskDarkMatterFactor
       remnantSpheroidGasMass=remnantSpheroidGasMass+satelliteSpheroidComponent%massGas()                                  +satelliteDiskComponent%massGas()
       remnantSpheroidMass   =remnantSpheroidMass   +satelliteSpheroidComponent%massGas()                                  +satelliteDiskComponent%massGas()
    case (movesToDisk)
       satelliteSpheroidMass =0.0d0
       satelliteRadius       =0.0d0
       angularMomentumFactor =angularMomentumFactor
    case (doesNotMove)
       satelliteSpheroidMass =                       satelliteSpheroidComponent%massGas()
       satelliteRadius       =                       satelliteSpheroidComponent%massGas()*satelliteSpheroidHalfMassRadius
       angularMomentumFactor =angularMomentumFactor +satelliteSpheroidComponent%massGas()*satelliteSpheroidDarkMatterFactor
       remnantSpheroidGasMass=remnantSpheroidGasMass+satelliteSpheroidComponent%massGas()
       remnantSpheroidMass   =remnantSpheroidMass   +satelliteSpheroidComponent%massGas()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisMergerStarsMoveTo)
    case (movesToSpheroid)
       satelliteSpheroidMass=satelliteSpheroidMass+satelliteSpheroidComponent%massStellar()                                  +satelliteDiskComponent%massStellar()
       satelliteRadius      =satelliteRadius      +satelliteSpheroidComponent%massStellar()*satelliteSpheroidHalfMassRadius  +satelliteDiskComponent%massStellar()*satelliteDiskHalfMassRadius
       angularMomentumFactor     =angularMomentumFactor     +satelliteSpheroidComponent%massStellar()*satelliteSpheroidDarkMatterFactor+satelliteDiskComponent%massStellar()*satelliteDiskDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +satelliteSpheroidComponent%massStellar()                                  +satelliteDiskComponent%massStellar()
    case (movesToDisk)
       satelliteSpheroidMass=satelliteSpheroidMass
       satelliteRadius      =satelliteRadius
       angularMomentumFactor     =angularMomentumFactor
    case (doesNotMove)
       satelliteSpheroidMass=satelliteSpheroidMass+satelliteSpheroidComponent%massStellar()
       satelliteRadius      =satelliteRadius      +satelliteSpheroidComponent%massStellar()*satelliteSpheroidHalfMassRadius
       angularMomentumFactor     =angularMomentumFactor     +satelliteSpheroidComponent%massStellar()*satelliteSpheroidDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +satelliteSpheroidComponent%massStellar()
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
    hostSpheroidMassPreMerger=hostSpheroidComponent%massStellar()+hostSpheroidComponent%massGas()
    if (hostRadius <= 0.0d0) hostSpheroidMassPreMerger=0.0d0
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Standard

end module Satellite_Merging_Remnant_Progenitors_Properties_Standard

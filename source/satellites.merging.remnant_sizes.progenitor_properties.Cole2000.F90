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

!% Contains a module which implements calculations of progenitor properties for merger remnant calculations using the algorithm of
!% \cite{cole_hierarchical_2000}.

module Satellite_Merging_Remnant_Progenitors_Properties_Cole2000
  !% Implements calculations of progenitor properties for merger remnant calculations using the algorithm of
  !% \cite{cole_hierarchical_2000}.
  use Galacticus_Nodes
  use Root_Finder
  implicit none
  private
  public :: Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init

  ! Module global variables used in root finding.
  type            (treeNode  ), pointer :: activeNode
  integer                               :: activeGasMovesTo, activeStarsMoveTo
  double precision                      :: activeHalfMass
  !$omp threadprivate(activeNode,activeGasMovesTo,activeStarsMoveTo,activeHalfMass)
contains

  !# <satelliteMergingRemnantProgenitorPropertiesMethod>
  !#  <unitName>Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init</unitName>
  !# </satelliteMergingRemnantProgenitorPropertiesMethod>
  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init(satelliteMergingRemnantProgenitorPropertiesMethod&
       &,Satellite_Merging_Remnant_Progenitor_Properties_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type     (varying_string                                          ), intent(in   )          :: satelliteMergingRemnantProgenitorPropertiesMethod
    procedure(Satellite_Merging_Remnant_Progenitor_Properties_Cole2000), intent(inout), pointer :: Satellite_Merging_Remnant_Progenitor_Properties_Get

    if (satelliteMergingRemnantProgenitorPropertiesMethod == 'Cole2000') then
       Satellite_Merging_Remnant_Progenitor_Properties_Get =>&
         & Satellite_Merging_Remnant_Progenitor_Properties_Cole2000
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
            & ) call Galacticus_Error_Report('Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init','this method requires that massStellar, massGas, halfMassRadius, and angularMomentum properties must all be gettable for both disk and spheroid')
    end if
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Cole2000_Init

  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Cole2000(satelliteNode,hostNode,satelliteMass,hostMass&
       &,satelliteSpheroidMass ,hostSpheroidMass,hostSpheroidMassPreMerger,satelliteRadius,hostRadius,angularMomentumFactor&
       &,remnantSpheroidMass ,remnantSpheroidGasMass)
    !% Computes various properties of the progenitor galaxies useful for calculations of merger remnant sizes.
    use Galactic_Structure_Radii
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Satellite_Merging_Mass_Movements_Descriptors
    use Numerical_Constants_Physical
    use Galacticus_Error
    implicit none
    type            (treeNode             ), intent(inout), pointer :: hostNode                       , satelliteNode
    double precision                       , intent(  out)          :: angularMomentumFactor          , hostMass                         , &
         &                                                             hostRadius                     , hostSpheroidMass                 , &
         &                                                             hostSpheroidMassPreMerger      , remnantSpheroidGasMass           , &
         &                                                             remnantSpheroidMass            , satelliteMass                    , &
         &                                                             satelliteRadius                , satelliteSpheroidMass
    class           (nodeComponentDisk    )               , pointer :: hostDiskComponent              , satelliteDiskComponent
    class           (nodeComponentSpheroid)               , pointer :: hostSpheroidComponent          , satelliteSpheroidComponent
    type            (rootFinder           ), save                   :: finder
    !$omp threadprivate(finder)
    double precision                                                :: componentMass                  , hostDiskDarkMatterFactor         , &
         &                                                             hostDiskHalfMassRadius         , hostSpheroidDarkMatterFactor     , &
         &                                                             hostSpheroidHalfMassRadius     , satelliteDiskDarkMatterFactor    , &
         &                                                             satelliteDiskHalfMassRadius    , satelliteSpheroidDarkMatterFactor, &
         &                                                             satelliteSpheroidHalfMassRadius

    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(Half_Mass_Radius_Root_Cole2000                  )
       call finder%tolerance   (toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6)
       call finder%rangeExpand (                                                             &
            &                   rangeExpandUpward            =2.0d0                        , &
            &                   rangeExpandDownward          =0.5d0                        , &
            &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
            &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
            &                   rangeExpandType              =rangeExpandMultiplicative      &
            &                  )
    end if

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
    componentMass                  =          hostSpheroidComponent%massStellar   () &
         &                          +         hostSpheroidComponent%massGas       ()
    hostSpheroidHalfMassRadius     =          hostSpheroidComponent%halfMassRadius()
    if (hostSpheroidHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       hostSpheroidDarkMatterFactor=           hostSpheroidComponent%angularMomentum()                           &
            &                       /(componentMass**1.5d0)                                                      &
            &                       /sqrt(gravitationalConstantGalacticus*          hostSpheroidHalfMassRadius)
    else
       hostSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass                  =          hostDiskComponent%massStellar   () &
         &                          +         hostDiskComponent%massGas       ()
    hostDiskHalfMassRadius         =          hostDiskComponent%halfMassRadius()
    if (hostDiskHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       hostDiskDarkMatterFactor         =          hostDiskComponent%angularMomentum()                           &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*         hostDiskHalfMassRadius)
    else
       hostDiskDarkMatterFactor=0.0d0
    end if
    componentMass                  = satelliteSpheroidComponent%massStellar   () &
         &                          +satelliteSpheroidComponent%massGas       ()
    satelliteSpheroidHalfMassRadius= satelliteSpheroidComponent%halfMassRadius()
    if (satelliteSpheroidHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteSpheroidDarkMatterFactor= satelliteSpheroidComponent%angularMomentum()                           &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*satelliteSpheroidHalfMassRadius)
    else
       satelliteSpheroidDarkMatterFactor=0.0d0
    end if
    componentMass                  =     satelliteDiskComponent%massStellar   () &
        &                           +    satelliteDiskComponent%massGas       ()
    satelliteDiskHalfMassRadius    =     satelliteDiskComponent%halfMassRadius()
    if (satelliteDiskHalfMassRadius > 0.0d0 .and. componentMass > 0.0d0) then
       satelliteDiskDarkMatterFactor    =     satelliteDiskComponent%angularMomentum()                           &
            &                            /(componentMass**1.5d0)                                                 &
            &                            /sqrt(gravitationalConstantGalacticus*    satelliteDiskHalfMassRadius)
    else
       satelliteDiskDarkMatterFactor=0.0d0
    end if

    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (thisHostGasMovesTo)
    case (movesToSpheroid)
       hostSpheroidMass      =hostSpheroidComponent%massGas()                             +hostDiskComponent%massGas()
       angularMomentumFactor =hostSpheroidComponent%massGas()*hostSpheroidDarkMatterFactor+hostDiskComponent%massGas()*hostDiskDarkMatterFactor
       remnantSpheroidGasMass=hostSpheroidComponent%massGas()                             +hostDiskComponent%massGas()
       remnantSpheroidMass   =hostSpheroidComponent%massGas()                             +hostDiskComponent%massGas()
    case (movesToDisk)
       hostSpheroidMass      =0.0d0
       angularMomentumFactor =0.0d0
       remnantSpheroidGasMass=0.0d0
       remnantSpheroidMass   =0.0d0
    case (doesNotMove)
       hostSpheroidMass      =hostSpheroidComponent%massGas()
       angularMomentumFactor =hostSpheroidComponent%massGas()*hostSpheroidDarkMatterFactor
       remnantSpheroidGasMass=hostSpheroidComponent%massGas()
       remnantSpheroidMass   =hostSpheroidComponent%massGas()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisHostStarsMoveTo)
    case (movesToSpheroid)
       hostSpheroidMass     =hostSpheroidMass     +hostSpheroidComponent%massStellar()                             +hostDiskComponent%massStellar()
       angularMomentumFactor=angularMomentumFactor+hostSpheroidComponent%massStellar()*hostSpheroidDarkMatterFactor+hostDiskComponent%massStellar()*hostDiskDarkMatterFactor
      remnantSpheroidMass   =remnantSpheroidMass  +hostSpheroidComponent%massStellar()                             +hostDiskComponent%massStellar()
    case (movesToDisk)
       hostSpheroidMass     =hostSpheroidMass
       angularMomentumFactor=angularMomentumFactor
    case (doesNotMove)
       hostSpheroidMass     =hostSpheroidMass     +hostSpheroidComponent%massStellar()
       angularMomentumFactor=angularMomentumFactor+hostSpheroidComponent%massStellar()*hostSpheroidDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +hostSpheroidComponent%massStellar()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisMergerGasMovesTo)
    case (movesToSpheroid)
       satelliteSpheroidMass =                       satelliteSpheroidComponent%massGas()                                  +satelliteDiskComponent%massGas()
       angularMomentumFactor =angularMomentumFactor +satelliteSpheroidComponent%massGas()*satelliteSpheroidDarkMatterFactor+satelliteDiskComponent%massGas()*satelliteDiskDarkMatterFactor
       remnantSpheroidGasMass=remnantSpheroidGasMass+satelliteSpheroidComponent%massGas()                                  +satelliteDiskComponent%massGas()
       remnantSpheroidMass   =remnantSpheroidMass   +satelliteSpheroidComponent%massGas()                                  +satelliteDiskComponent%massGas()
    case (movesToDisk)
       satelliteSpheroidMass =0.0d0
       angularMomentumFactor =angularMomentumFactor
    case (doesNotMove)
       satelliteSpheroidMass=                        satelliteSpheroidComponent%massGas()
       angularMomentumFactor =angularMomentumFactor +satelliteSpheroidComponent%massGas()*satelliteSpheroidDarkMatterFactor
       remnantSpheroidGasMass=remnantSpheroidGasMass+satelliteSpheroidComponent%massGas()
       remnantSpheroidMass   =remnantSpheroidMass   +satelliteSpheroidComponent%massGas()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Sizes_Utilities','unrecognized moveTo descriptor')
    end select
    select case (thisMergerStarsMoveTo)
    case (movesToSpheroid)
       satelliteSpheroidMass=satelliteSpheroidMass+satelliteSpheroidComponent%massStellar()                                  +satelliteDiskComponent%massStellar()
       angularMomentumFactor=angularMomentumFactor+satelliteSpheroidComponent%massStellar()*satelliteSpheroidDarkMatterFactor+satelliteDiskComponent%massStellar()*satelliteDiskDarkMatterFactor
       remnantSpheroidMass  =remnantSpheroidMass  +satelliteSpheroidComponent%massStellar()                                  +satelliteDiskComponent%massStellar()
    case (movesToDisk)
       satelliteSpheroidMass=satelliteSpheroidMass
       angularMomentumFactor=angularMomentumFactor
    case (doesNotMove)
       satelliteSpheroidMass=satelliteSpheroidMass+satelliteSpheroidComponent%massStellar()
       remnantSpheroidMass  =remnantSpheroidMass  +satelliteSpheroidComponent%massStellar()
       angularMomentumFactor=angularMomentumFactor+satelliteSpheroidComponent%massStellar()*satelliteSpheroidDarkMatterFactor
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
       activeHalfMass   =  0.5d0*Half_Mass_Radius_Root_Cole2000(radiusLarge)
       hostRadius=finder%find(rootGuess=Galactic_Structure_Radius_Enclosing_Mass(                                 &
            &                                                                    activeNode                     , &
            &                                                                    fractionalMass=0.50d0          , &
            &                                                                    massType      =massTypeGalactic  &
            &                                                                   )                                 &
            &                )
    else
       hostRadius=0.0d0
    end if
    if (satelliteSpheroidMass > 0.0d0) then
       activeNode       => satelliteNode
       activeGasMovesTo =  thisMergerGasMovesTo
       activeStarsMoveTo=  thisMergerStarsMoveTo
       activeHalfMass   =  0.0d0
       activeHalfMass   =  0.50d0*Half_Mass_Radius_Root_Cole2000(radiusLarge)
       satelliteRadius=finder%find(rootGuess=Galactic_Structure_Radius_Enclosing_Mass(                                 &
            &                                                                         activeNode                     , &
            &                                                                         fractionalMass=0.50d0          , &
            &                                                                         massType      =massTypeGalactic  &
            &                                                                        )                                 &
            &                     )
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
    hostSpheroidMassPreMerger=hostSpheroidComponent%massStellar()+hostSpheroidComponent%massGas()
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Cole2000

  double precision function Half_Mass_Radius_Root_Cole2000(radius)
    !% Function used in root finding for progenitor galaxy half-mass radii.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    double precision, intent(in   ) :: radius

    ! Initialize enclosed mass to negative of the half mass.
    Half_Mass_Radius_Root_Cole2000=-activeHalfMass

    ! Account for gas mass.
    select case (activeGasMovesTo)
    case (movesToSpheroid)
       Half_Mass_Radius_Root_Cole2000= Half_Mass_Radius_Root_Cole2000                                                                          &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,componentType=componentTypeSpheroid                         ) &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,componentType=componentTypeDisk    ,massType=massTypeGaseous)
    case (doesNotMove    )
       Half_Mass_Radius_Root_Cole2000= Half_Mass_Radius_Root_Cole2000                                                                          &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,componentType=componentTypeSpheroid,massType=massTypeGaseous)
    end select

    ! Account for gas mass.
    select case (activeStarsMoveTo)
    case (movesToSpheroid)
       Half_Mass_Radius_Root_Cole2000= Half_Mass_Radius_Root_Cole2000                                                                          &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,componentType=componentTypeSpheroid                         ) &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,componentType=componentTypeDisk    ,massType=massTypeStellar)
    case (doesNotMove    )
       Half_Mass_Radius_Root_Cole2000= Half_Mass_Radius_Root_Cole2000                                                                          &
            &                +Galactic_Structure_Enclosed_Mass(activeNode,radius,componentType=componentTypeSpheroid,massType=massTypeStellar)
    end select
    return
  end function Half_Mass_Radius_Root_Cole2000

end module Satellite_Merging_Remnant_Progenitors_Properties_Cole2000

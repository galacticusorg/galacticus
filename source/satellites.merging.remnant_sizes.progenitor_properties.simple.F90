!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of progenitor properties for merger remnant
!% calculations using a simplified algorithm.

module Satellite_Merging_Remnant_Progenitors_Properties_Simple
  !% Implements calculations of progenitor properties for merger remnant calculations using a
  !% simplified algorithm.
  use Galacticus_Nodes
  implicit none
  private
  public :: Satellite_Merging_Remnant_Progenitor_Properties_Simple_Init

contains

  !# <satelliteMergingRemnantProgenitorPropertiesMethod>
  !#  <unitName>Satellite_Merging_Remnant_Progenitor_Properties_Simple_Init</unitName>
  !# </satelliteMergingRemnantProgenitorPropertiesMethod>
  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Simple_Init(satelliteMergingRemnantProgenitorPropertiesMethod&
       &,Satellite_Merging_Remnant_Progenitor_Properties_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Galacticus_Error
    use Array_Utilities
    implicit none
    type     (varying_string                                        ), intent(in   )          :: satelliteMergingRemnantProgenitorPropertiesMethod
    procedure(Satellite_Merging_Remnant_Progenitor_Properties_Simple), intent(inout), pointer :: Satellite_Merging_Remnant_Progenitor_Properties_Get

    if (satelliteMergingRemnantProgenitorPropertiesMethod == 'simple') then
       Satellite_Merging_Remnant_Progenitor_Properties_Get => Satellite_Merging_Remnant_Progenitor_Properties_Simple
       ! Ensure that required methods are supported.
       if     (                                                                                                                                          &
            &  .not.                                                                                                                                     &
            &       (                                                                                                                                    &
            &        defaultDiskComponent    %    massStellarIsGettable().and.                                                                           &
            &        defaultDiskComponent    %        massGasIsGettable().and.                                                                           &
            &        defaultDiskComponent    % halfMassRadiusIsGettable()                                                                                &
            &  )                                                                                                                                         &
            & ) call Galacticus_Error_Report                                                                                                             &
            &        (                                                                                                                                   &
            &         'Satellite_Merging_Remnant_Progenitor_Properties_Simple_Init'                                                                   ,  &
            &         'this method requires that massStellar, massGas, and halfMassRadius properties must all be gettable for the disk component.'    // &
            &         Galacticus_Component_List(                                                                                                         &
            &                                   'disk'                                                                                                ,  &
            &                                   defaultDiskComponent    %    massStellarAttributeMatch(requireGettable=.true.).intersection.             &
            &                                   defaultDiskComponent    %        massGasAttributeMatch(requireGettable=.true.).intersection.             &
            &                                   defaultDiskComponent    % halfMassRadiusAttributeMatch(requireGettable=.true.)                           &
            &                                  )                                                                                                         &
            &        )
       if     (                                                                                                                                          &
            &  .not.                                                                                                                                     &
            &       (                                                                                                                                    &
            &        defaultSpheroidComponent%    massStellarIsGettable().and.                                                                           &
            &        defaultSpheroidComponent%        massGasIsGettable().and.                                                                           &
            &        defaultSpheroidComponent% halfMassRadiusIsGettable()                                                                                &
            &  )                                                                                                                                         &
            & ) call Galacticus_Error_Report                                                                                                             &
            &        (                                                                                                                                   &
            &         'Satellite_Merging_Remnant_Progenitor_Properties_Simple_Init'                                                                   ,  &
            &         'this method requires that massStellar, massGas, and halfMassRadius properties must all be gettable for the spheroid component.'// &
            &         Galacticus_Component_List(                                                                                                         &
            &                                   'spheroid'                                                                                            ,  &
            &                                   defaultSpheroidComponent%    massStellarAttributeMatch(requireGettable=.true.).intersection.             &
            &                                   defaultSpheroidComponent%        massGasAttributeMatch(requireGettable=.true.).intersection.             &
            &                                   defaultSpheroidComponent% halfMassRadiusAttributeMatch(requireGettable=.true.)                           &
            &                                  )                                                                                                         &
            &        )
    end if
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Simple_Init

  subroutine Satellite_Merging_Remnant_Progenitor_Properties_Simple(satelliteNode,hostNode,satelliteMass,hostMass&
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
    type            (treeNode             ), intent(inout), target :: hostNode                       , satelliteNode
    double precision                       , intent(  out)         :: angularMomentumFactor          , hostMass                       , &
         &                                                            hostRadius                     , hostSpheroidMass               , &
         &                                                            hostSpheroidMassPreMerger      , remnantSpheroidGasMass         , &
         &                                                            remnantSpheroidMass            , satelliteMass                  , &
         &                                                            satelliteRadius                , satelliteSpheroidMass
    class           (nodeComponentDisk    ), pointer               :: hostDiskComponent              , satelliteDiskComponent
    class           (nodeComponentSpheroid), pointer               :: hostSpheroidComponent          , satelliteSpheroidComponent

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

    ! Find the masses of material that will end up in the spheroid component of the remnant.
    select case (thisHostGasMovesTo)
    case (movesToSpheroid)
       hostSpheroidMass      =+hostSpheroidComponent%massGas()+hostDiskComponent    %massGas       ()
       remnantSpheroidGasMass=+hostSpheroidComponent%massGas()+hostDiskComponent    %massGas       ()
       remnantSpheroidMass   =+hostSpheroidComponent%massGas()+hostDiskComponent    %massGas       ()
       hostRadius            =+hostSpheroidComponent%massGas()*hostSpheroidComponent%halfMassRadius() &
            &                 +hostDiskComponent    %massGas()*hostDiskComponent    %halfMassRadius()
    case (movesToDisk)
       hostSpheroidMass      =0.0d0
       remnantSpheroidGasMass=0.0d0
       remnantSpheroidMass   =0.0d0
    case (doesNotMove)
       hostSpheroidMass      =hostSpheroidComponent%massGas()
       remnantSpheroidGasMass=hostSpheroidComponent%massGas()
       remnantSpheroidMass   =hostSpheroidComponent%massGas()
       hostRadius            =hostSpheroidComponent%massGas()*hostSpheroidComponent%halfMassRadius()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Progenitor_Properties_Simple','unrecognized moveTo descriptor')
    end select
    select case (thisHostStarsMoveTo)
    case (movesToSpheroid)
       hostSpheroidMass     =hostSpheroidMass       +hostSpheroidComponent%massStellar()+hostDiskComponent    %massStellar   ()
       remnantSpheroidMass  =remnantSpheroidMass    +hostSpheroidComponent%massStellar()+hostDiskComponent    %massStellar   ()
       hostRadius           =hostRadius             +hostSpheroidComponent%massStellar()*hostSpheroidComponent%halfMassRadius() &
            &                                       +hostDiskComponent    %massStellar()*hostDiskComponent    %halfMassRadius()
    case (movesToDisk)
       hostSpheroidMass     =hostSpheroidMass
    case (doesNotMove)
       hostSpheroidMass     =hostSpheroidMass       +hostSpheroidComponent%massStellar()
       remnantSpheroidMass  =remnantSpheroidMass    +hostSpheroidComponent%massStellar()
       hostRadius           =hostRadius             +hostSpheroidComponent%massStellar()*hostSpheroidComponent%halfMassRadius()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Progenitor_Properties_Simple','unrecognized moveTo descriptor')
    end select
    select case (thisMergerGasMovesTo)
    case (movesToSpheroid)
       satelliteSpheroidMass =                       satelliteSpheroidComponent%massGas()+satelliteDiskComponent    %massGas       ()
       remnantSpheroidGasMass=remnantSpheroidGasMass+satelliteSpheroidComponent%massGas()+satelliteDiskComponent    %massGas       ()
       remnantSpheroidMass   =remnantSpheroidMass   +satelliteSpheroidComponent%massGas()+satelliteDiskComponent    %massGas       ()
       satelliteRadius       =                      +satelliteSpheroidComponent%massGas()*satelliteSpheroidComponent%halfMassRadius() &
            &                                       +satelliteDiskComponent    %massGas()*satelliteDiskComponent    %halfMassRadius()
    case (movesToDisk)
       satelliteSpheroidMass =0.0d0
    case (doesNotMove)
       satelliteSpheroidMass=                        satelliteSpheroidComponent%massGas()
       remnantSpheroidGasMass=remnantSpheroidGasMass+satelliteSpheroidComponent%massGas()
       remnantSpheroidMass   =remnantSpheroidMass   +satelliteSpheroidComponent%massGas()
       satelliteRadius       =                      +satelliteSpheroidComponent%massGas()*satelliteSpheroidComponent%halfMassRadius()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Progenitor_Properties_Simple','unrecognized moveTo descriptor')
    end select
    select case (thisMergerStarsMoveTo)
    case (movesToSpheroid)
       satelliteSpheroidMass=satelliteSpheroidMass+satelliteSpheroidComponent%massStellar()+satelliteDiskComponent    %massStellar   ()
       remnantSpheroidMass  =remnantSpheroidMass  +satelliteSpheroidComponent%massStellar()+satelliteDiskComponent    %massStellar   ()
       satelliteRadius      =satelliteRadius      +satelliteSpheroidComponent%massStellar()*satelliteSpheroidComponent%halfMassRadius() &
            &                                     +satelliteDiskComponent    %massStellar()*satelliteDiskComponent    %halfMassRadius()
    case (movesToDisk)
       satelliteSpheroidMass=satelliteSpheroidMass
    case (doesNotMove)
       satelliteSpheroidMass=satelliteSpheroidMass+satelliteSpheroidComponent%massStellar()
       remnantSpheroidMass  =remnantSpheroidMass  +satelliteSpheroidComponent%massStellar()
       satelliteRadius      =satelliteRadius      +satelliteSpheroidComponent%massStellar()*satelliteSpheroidComponent%halfMassRadius()
    case default
       call Galacticus_Error_Report('Satellite_Merging_Remnant_Progenitor_Properties_Simple','unrecognized moveTo descriptor')
    end select

    ! Compute the half-mass radii of the material that will end up in the remnant spheroid.
    if (hostSpheroidMass      > 0.0d0) hostRadius     =hostRadius     /hostSpheroidMass
    if (satelliteSpheroidMass > 0.0d0) satelliteRadius=satelliteRadius/satelliteSpheroidMass

    ! Compute the angular momentum factor.
    angularMomentumFactor=1.0d0

    ! Compute the mass of the host spheroid before the merger.
    hostSpheroidMassPreMerger=hostSpheroidComponent%massStellar()+hostSpheroidComponent%massGas()
    return
  end subroutine Satellite_Merging_Remnant_Progenitor_Properties_Simple

end module Satellite_Merging_Remnant_Progenitors_Properties_Simple

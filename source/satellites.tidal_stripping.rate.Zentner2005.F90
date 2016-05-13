!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

!% Contains a module with a \cite{king_structure_1962} implementation of calculations of satellite mass loss due to tidal
!% stripping.

module Tidal_Stripping_Rate_Zentner2005
  !% Implements the \cite{king_structure_1962} calculation of satellite mass loss due to tidal stripping.
  use Galacticus_Nodes
  implicit none
  private
  public :: Satellite_Tidal_Stripping_Rate_Zentner2005_Initialize

  ! Module scope variables used in root finding.
  double precision                     :: tidalPullGlobal
  type            (treeNode),  pointer :: activeNode
  !$omp threadprivate(tidalPullGlobal,activeNode)

  ! Dimensionless mass loss rate.
  double precision                     :: satelliteTidalStrippingZentner2005Rate

contains

  !# <satelliteTidalStrippingMethod>
  !#  <unitName>Satellite_Tidal_Stripping_Rate_Zentner2005_Initialize</unitName>
  !# </satelliteTidalStrippingMethod>
  subroutine Satellite_Tidal_Stripping_Rate_Zentner2005_Initialize(satelliteTidalStrippingMethod,Satellite_Tidal_Stripping_Rate)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                     ),          intent(in   ) :: satelliteTidalStrippingMethod
    procedure(Satellite_Tidal_Stripping_Rate_Zentner2005), pointer, intent(inout) :: Satellite_Tidal_Stripping_Rate

    if (satelliteTidalStrippingMethod == 'Zentner2005') then
       Satellite_Tidal_Stripping_Rate => Satellite_Tidal_Stripping_Rate_Zentner2005
       !@ <inputParameter>
       !@   <name>satelliteTidalStrippingZentner2005Rate</name>
       !@   <defaultValue>2.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The dimensionless rate coefficient apeparing in the \cite{zentner_physics_2005} expression for the tidal mass loss rate from subhalos.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group></group>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteTidalStrippingZentner2005Rate',satelliteTidalStrippingZentner2005Rate,defaultValue=2.5d0)
   end if
    return
  end subroutine Satellite_Tidal_Stripping_Rate_Zentner2005_Initialize

  double precision function Satellite_Tidal_Stripping_Rate_Zentner2005(thisNode)
    !% Return the Zentner2005 mass loss rate for satellites due to tidal stripping.
    use Galacticus_Nodes
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    use Error_Functions
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Vectors
    use Dark_Matter_Halo_Scales
    use Root_Finder
    implicit none
    class           (nodeComponentSatellite), pointer                     :: thisSatellite
    type            (treeNode              ), pointer     , intent(inout) :: thisNode
    type            (treeNode              ), pointer                     :: hostNode
    double precision                        , dimension(3)                :: position,velocity
    double precision                        , parameter                   :: toleranceAbsolute      =0.0d0 , toleranceRelative=1.0d-3
    double precision                        , parameter                   :: radiusZero             =0.0d0
    double precision                        , parameter                   :: tidalRadiusTinyFraction=1.0d-6
    type            (rootFinder            ), save                        :: finder
    !$omp threadprivate(finder)
    double precision                                                      :: satelliteMass          , parentDensity           , &
         &                                                                   parentEnclosedMass     , angularVelocity         , &
         &                                                                   orbitalPeriod          , radius                  , &
         &                                                                   tidalTensor            , tidalRadius             , &
         &                                                                   outerSatelliteMass     , radialTimeScale

    ! Get required quantities from satellite and host nodes.
    hostNode          => thisNode     %mergesWith()
    thisSatellite     => thisNode     %satellite ()
    satelliteMass     =  thisSatellite%boundMass ()
    position          =  thisSatellite%position  ()
    velocity          =  thisSatellite%velocity  ()
    radius            =  Vector_Magnitude                (         position                          )
    parentDensity     =  Galactic_Structure_Density      (hostNode,position,coordinateSystemCartesian)
    parentEnclosedMass=  Galactic_Structure_Enclosed_Mass(hostNode,radius                            )
    ! Compute the orbital period.
    angularVelocity   =  Vector_Magnitude(Vector_Product(position,velocity)) &
         &              /radius**2                                           &
         &              *kilo                                                &
         &              *gigaYear                                            &
         &              /megaParsec
    radialTimescale   =  abs             (   Dot_Product(position,velocity)) &
         &              /radius**2                                           &
         &              *kilo                                                &
         &              *gigaYear                                            &
         &              /megaParsec
    orbitalPeriod     =  1.0d0/max(angularVelocity/2.0d0/Pi,radialTimescale)
    ! Find the tidal tensor.
    tidalTensor       = -gravitationalConstantGalacticus &
         &              *(                               &
         &                 2.0d0                         &
         &                *parentEnclosedMass            &
         &                /radius**3                     &
         &                -4.0d0                         &
         &                *Pi                            &
         &                *parentDensity                 &
         &               )                               &
         &              *(kilo*gigaYear/megaParsec)**2
    ! If the tidal force is stretching (not compressing), compute the tidal radius.
    if     (                                                                     &
         &   angularVelocity**2                                    > tidalTensor &
         &  .and.                                                                &
         &   satelliteMass                                         >  0.0d0      &
         &  .and.                                                                &
         &   Galactic_Structure_Enclosed_Mass(thisNode,radiusZero) >= 0.0d0      &
         & ) then
       ! Initial estimate of the tidal radius.
       tidalPullGlobal =  angularVelocity**2-tidalTensor
       tidalRadius     =                                  &
            &           (                                 &
            &            +gravitationalConstantGalacticus &
            &            *satelliteMass                   &
            &            /tidalPullGlobal                 &
            &            *(kilo*gigaYear/megaParsec)**2   &
            &           )**(1.0d0/3.0d0)       
       ! Find the tidal radius in the dark matter profile.
       if (.not.finder%isInitialized()) then
          call finder%rootFunction(Tidal_Radius_Solver)
          call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          call finder%rangeExpand (                                                             &
               &                   rangeExpandUpward            =2.0d0                        , &
               &                   rangeExpandDownward          =0.5d0                        , &
               &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                   rangeExpandType              =rangeExpandMultiplicative      &
               &                  )
       end if
       activeNode      => thisNode
       ! Check for complete stripping.
       if (Tidal_Radius_Solver(tidalRadiusTinyFraction*tidalRadius) > 0.0d0) then
          tidalRadius=0.0d0
       else
          tidalRadius=finder%find(rootGuess=tidalRadius)
       end if
       outerSatelliteMass=max(                                                                                                   &
            &                 Galactic_Structure_Enclosed_Mass(thisNode)-Galactic_Structure_Enclosed_Mass(thisNode,tidalRadius), &
            &                 0.0d0                                                                                              &
            &                )
    else
       outerSatelliteMass=0.0d0
    end if
    ! Compute the rate of mass loss.
    Satellite_Tidal_Stripping_Rate_Zentner2005=-satelliteTidalStrippingZentner2005Rate*outerSatelliteMass/orbitalPeriod
    return
  end function Satellite_Tidal_Stripping_Rate_Zentner2005

  double precision function Tidal_Radius_Solver(radius)
    !% Root function used to find the tidal radius within a subhalo.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    double precision                        , intent(in   ) :: radius
    class           (nodeComponentSatellite), pointer       :: satelliteComponent
    double precision                                        :: enclosedMass

    ! Get the satellite component.
    satelliteComponent => activeNode%satellite()
    enclosedMass       =  Galactic_Structure_Enclosed_Mass(activeNode,radius)
    Tidal_Radius_Solver=+tidalPullGlobal                    &
         &              -gravitationalConstantGalacticus    &
         &              *enclosedMass                       &
         &              /radius                         **3 &
         &              *(                                  &
         &                +kilo                             &
         &                *gigaYear                         &
         &                /megaParsec                       &
         &               )                              **2
    return
  end function Tidal_Radius_Solver

end module Tidal_Stripping_Rate_Zentner2005

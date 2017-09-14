!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a model of ram pressure stripping of hot halos based on the methods of
!% \cite{font_colours_2008}.

module Hot_Halo_Ram_Pressure_Stripping_Font2008
  !% Implements a module which implements a model of ram pressure stripping of hot halos based on the methods of
  !% \cite{font_colours_2008}.
  use Kind_Numbers
  use Galacticus_Nodes
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Stripping_Font2008_Initialize

  ! Pointers to the satellite node.
  type            (treeNode      ), pointer :: satelliteNode
  !$omp threadprivate(satelliteNode)
  ! The ram pressure force (per unit area) used in root finding.
  double precision                          :: ramPressureForce
  !$omp threadprivate(ramPressureForce)
  ! Parameters of the ram pressure stripping model.
  double precision                          :: ramPressureStrippingFormFactor
  ! Optimization.
  integer         (kind=kind_int8)          :: lastUniqueID                  =-1
  double precision                          :: radiusLast                    =-1.0d0
  !$omp threadprivate(lastUniqueID,radiusLast)

  
contains
 
  !# <hotHaloRamPressureStrippingMethod>
  !#  <unitName>Hot_Halo_Ram_Pressure_Stripping_Font2008_Initialize</unitName>
  !# </hotHaloRamPressureStrippingMethod>
  subroutine Hot_Halo_Ram_Pressure_Stripping_Font2008_Initialize(hotHaloRamPressureStrippingMethod,Hot_Halo_Ram_Pressure_Stripping_Get)
    !% Initializes the ``Font2008'' hot halo ram pressure stripping module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                              ), intent(in   )          :: hotHaloRamPressureStrippingMethod
    procedure(Hot_Halo_Ram_Pressure_Stripping_Font2008_Get), intent(inout), pointer :: Hot_Halo_Ram_Pressure_Stripping_Get

    if (hotHaloRamPressureStrippingMethod == 'Font2008') then
       Hot_Halo_Ram_Pressure_Stripping_Get => Hot_Halo_Ram_Pressure_Stripping_Font2008_Get
       !# <inputParameter>
       !#   <name>ramPressureStrippingFormFactor</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>2.0d0</defaultValue>
       !#   <description>The form factor appearing in the gravitational binding force (per unit area) in the ram pressure stripping model
       !#      of \citeauthor{font_colours_2008}~(\citeyear{font_colours_2008}; their eqn.~4).</description>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
    end if
    return
  end subroutine Hot_Halo_Ram_Pressure_Stripping_Font2008_Initialize

  double precision function Hot_Halo_Ram_Pressure_Stripping_Font2008_Get(node)
    !% Computes the hot halo ram pressure stripping radius, assuming a null calculation in which that radius always equals the
    !% virial radius.
    use Dark_Matter_Halo_Scales
    use Hot_Halo_Ram_Pressure_Forces
    use Root_Finder
    implicit none
    type            (treeNode                ), intent(inout), target :: node
    double precision                          , parameter             :: toleranceAbsolute             =0.0d0 , toleranceRelative=1.0d-3
    double precision                          , parameter             :: radiusSmallestOverRadiusVirial=1.0d-6
    class           (darkMatterHaloScaleClass), pointer               :: darkMatterHaloScale_
    type            (rootFinder              ), save                  :: finder
    !$omp threadprivate(finder)
    double precision                                                  :: virialRadius
 
    ! Get the virial radius of the satellite.
    darkMatterHaloScale_ => darkMatterHaloScale              (        )
    virialRadius         =  darkMatterHaloScale_%virialRadius(node)
    ! Test whether node is a satellite.
    if (node%isSatellite()) then
       ! Set a pointer to the satellite node.
       satelliteNode => node
       ! Get the ram pressure force due to the hot halo.
       ramPressureForce=Hot_Halo_Ram_Pressure_Force(node)
       ! Find the radial range within which the ram pressure radius must lie.
       if      (Hot_Halo_Ram_Pressure_Stripping_Radius_Solver(                               virialRadius) >= 0.0d0) then
          ! The ram pressure force is not sufficiently strong to strip even at the satellite virial radius - simply return the
          ! virial radius as the stripping radius in this case.
          Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=virialRadius
       else
          if (Hot_Halo_Ram_Pressure_Stripping_Radius_Solver(radiusSmallestOverRadiusVirial*virialRadius) <= 0.0d0) then
             ! The ram pressure force can strip to (essentially) arbitrarily small radii.
             Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=0.0d0
          else
             ! Solver for the ram pressure stripping radius.
             if (.not.finder%isInitialized()) then
                call finder%rootFunction(Hot_Halo_Ram_Pressure_Stripping_Radius_Solver)
                call finder%tolerance   (toleranceAbsolute,toleranceRelative          )
             end if
             ! If we have a previously found radius, and if the node is the same as the previous node for which this function was
             ! called, then use that previous radius as a guess for the new solution. Note that we do not reset the previous
             ! radius between ODE steps (i.e. we do not make use of "calculationReset" events to reset the radius) as we
             ! specifically want to retain knowledge from the previous step.
             if (radiusLast > 0.0d0 .and. node%uniqueID() == lastUniqueID) then
                call finder%rangeExpand(                                                                           &
                     &                  rangeExpandDownward          =0.9d0                                      , &
                     &                  rangeExpandUpward            =1.1d0                                      , &
                     &                  rangeExpandType              =rangeExpandMultiplicative                  , &
                     &                  rangeDownwardLimit           =radiusSmallestOverRadiusVirial*virialRadius, &
                     &                  rangeUpwardLimit             =                               virialRadius, &
                     &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive              , &
                     &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative                &
                     &                 )
             else
                call finder%rangeExpand(                                                                           &
                     &                  rangeExpandDownward          =0.5d0                                      , &
                     &                  rangeExpandType              =rangeExpandMultiplicative                  , &
                     &                  rangeDownwardLimit           =radiusSmallestOverRadiusVirial*virialRadius, &
                     &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive                &
                     &                 )
                radiusLast  =virialRadius
                lastUniqueID=node%uniqueID()
             end if
             Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=finder%find(rootGuess=radiusLast)
             radiusLast=Hot_Halo_Ram_Pressure_Stripping_Font2008_Get
          end if
       end if
    else
       ! If node is not a satellite, return a ram pressure stripping radius equal to the virial radius.
       Hot_Halo_Ram_Pressure_Stripping_Font2008_Get=virialRadius
    end if
    return
  end function Hot_Halo_Ram_Pressure_Stripping_Font2008_Get

  double precision function Hot_Halo_Ram_Pressure_Stripping_Radius_Solver(radius)
    !% Root function used in finding the ram pressure stripping radius.
    use Hot_Halo_Mass_Distributions
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    double precision                              , intent(in   ) :: radius
    class           (hotHaloMassDistributionClass), pointer       :: hotHaloMassDistribution_
    double precision                                              :: enclosedMass            , gravitationalBindingForce, &
         &                                                           hotHaloDensity

    ! Get the hot halo mass distribution.
    hotHaloMassDistribution_ => hotHaloMassDistribution()
    enclosedMass             =Galactic_Structure_Enclosed_Mass(satelliteNode,radius,massType=massTypeAll,componentType=componentTypeAll)
    hotHaloDensity           =hotHaloMassDistribution_%density(satelliteNode,radius)
    gravitationalBindingForce=ramPressureStrippingFormFactor*gravitationalConstantGalacticus*enclosedMass*hotHaloDensity/radius
    if (gravitationalBindingForce >= 0.0d0) then
       Hot_Halo_Ram_Pressure_Stripping_Radius_Solver=gravitationalBindingForce-ramPressureForce
    else
       Hot_Halo_Ram_Pressure_Stripping_Radius_Solver=                         -ramPressureForce
    end if
    return
  end function Hot_Halo_Ram_Pressure_Stripping_Radius_Solver

end module Hot_Halo_Ram_Pressure_Stripping_Font2008

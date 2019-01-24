!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which provides a hot halo mass distribution class.

module Hot_Halo_Mass_Distributions
  !% Provides an object which provides a hot halo mass distribution class.
  use Galacticus_Nodes, only : treeNode
  private
  public :: hotHaloMassDistributionDensity     , hotHaloMassDistributionRotationCurve        , &
       &    hotHaloMassDistributionEnclosedMass, hotHaloMassDistributionRotationCurveGradient

  !# <functionClass>
  !#  <name>hotHaloMassDistribution</name>
  !#  <descriptiveName>Hot Halo Mass Distributions</descriptiveName>
  !#  <description>Object implementing hot halo mass distributions.</description>
  !#  <default>betaProfile</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="density" >
  !#   <description>Return the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: radius</argument>
  !#  </method>
  !#  <method name="densityLogSlope" >
  !#   <description>Return the logarithmic slope of the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   )          :: radius</argument>
  !#  </method>
  !#  <method name="enclosedMass" >
  !#   <description>Return the mass enclosed in the hot halo at the given {\normalfont \ttfamily radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout), target :: node</argument>
  !#   <argument>double precision          , intent(in   )         :: radius</argument>
  !#  </method>
  !#  <method name="radialMoment" >
  !#   <description>Return the density of the hot halo at the given {\normalfont \ttfamily radius}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: moment, radius</argument>
  !#  </method>
  !#  <method name="rotationNormalization" >
  !#   <description>Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant in radius) for {\normalfont \ttfamily node}. Specifically, the normalization, $A$, returned is such that $V_\mathrm{rot} = A J/M$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#  </method>
  !# </functionClass>

contains
  
  !# <enclosedMassTask>
  !#  <unitName>hotHaloMassDistributionEnclosedMass</unitName>
  !# </enclosedMassTask>
  double precision function hotHaloMassDistributionEnclosedMass(node,radius,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the mass within a given radius for a dark matter profile.
    use Galactic_Structure_Options
    use Galacticus_Nodes          , only : nodeComponentHotHalo
    implicit none
    type            (treeNode                    ), intent(inout)           :: node
    integer                                       , intent(in   )           :: componentType           , massType   , &
         &                                                                     weightBy                , weightIndex
    double precision                              , intent(in   )           :: radius
    logical                                       , intent(in   ), optional :: haloLoaded
    class           (hotHaloMassDistributionClass)               , pointer  :: hotHaloMassDistribution_
    class           (nodeComponentHotHalo        )               , pointer  :: hotHalo
    !GCC$ attributes unused :: haloLoaded, weightIndex

    ! Return zero mass if the requested mass type or component is not matched.
    hotHaloMassDistributionEnclosedMass=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    ! Return the enclosed mass.
    if (radius >= radiusLarge) then
       hotHalo                             => node   %hotHalo()
       hotHaloMassDistributionEnclosedMass =  hotHalo%mass   ()
    else
       hotHaloMassDistribution_            =>     hotHaloMassDistribution              (           )
       hotHaloMassDistributionEnclosedMass =  max(hotHaloMassDistribution_%enclosedMass(node,radius),0.0d0)
    end if
    return
  end function hotHaloMassDistributionEnclosedMass

  !# <rotationCurveTask>
  !#  <unitName>hotHaloMassDistributionRotationCurve</unitName>
  !# </rotationCurveTask>
  double precision function hotHaloMassDistributionRotationCurve(node,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve at a given radius for the hot halo density profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    type            (treeNode), intent(inout)           :: node
    integer                   , intent(in   )           :: componentType, massType
    double precision          , intent(in   )           :: radius
    logical                   , intent(in   ), optional :: haloLoaded
    double precision                                    :: componentMass

    ! Compute rotation curve if radius is non-zero.
    hotHaloMassDistributionRotationCurve=0.0d0
    if (radius > 0.0d0) then
       componentMass=hotHaloMassDistributionEnclosedMass(node,radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
       if (componentMass > 0.0d0)                                                         &
            & hotHaloMassDistributionRotationCurve=+sqrt(                                 &
            &                                            +gravitationalConstantGalacticus &
            &                                            *componentMass                   &
            &                                            /radius                          &
            &                                           )
    end if
    return
  end function hotHaloMassDistributionRotationCurve

  !# <rotationCurveGradientTask>
  !#  <unitName>hotHaloMassDistributionRotationCurveGradient</unitName>
  !# </rotationCurveGradientTask>
  double precision function hotHaloMassDistributionRotationCurveGradient(node,radius,componentType,massType,haloLoaded)
    !% Computes the rotation curve gradient at a given radius for the hot halo density profile.
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    implicit none
    type            (treeNode                    ), intent(inout)           :: node
    integer                                       , intent(in   )           :: componentType   , massType
    double precision                              , intent(in   )           :: radius
    logical                                       , intent(in   ), optional :: haloLoaded
    class           (hotHaloMassDistributionClass)               , pointer  :: hotHalo
    double precision                                                        :: componentDensity, componentMass

    ! Set to zero by default.
    hotHaloMassDistributionRotationCurveGradient=0.0d0
    ! Compute if a spheroid is present.
    if (radius > 0.0d0) then
       componentMass=hotHaloMassDistributionEnclosedMass(node,radius,componentType,massType,weightByMass,weightIndexNull,haloLoaded)
       if (componentMass > 0.0d0) then
          hotHalo                                      =>  hotHaloMassDistribution        (           )
          componentDensity                             =   hotHalo                %density(node,radius)
          hotHaloMassDistributionRotationCurveGradient =  +gravitationalConstantGalacticus    &
               &                                          *(                                  &
               &                                            -componentMass                    &
               &                                            /radius                       **2 &
               &                                            +4.0d0                            &
               &                                            *Pi                               &
               &                                            *radius                           &
               &                                            *componentDensity                 &
               &                                           )
       end if
    end if
    return
  end function hotHaloMassDistributionRotationCurveGradient

  !# <densityTask>
  !#  <unitName>hotHaloMassDistributionDensity</unitName>
  !# </densityTask>
  double precision function hotHaloMassDistributionDensity(node,positionSpherical,componentType,massType,weightBy,weightIndex,haloLoaded)
    !% Computes the density at a given position for a dark matter profile.
    use Galactic_Structure_Options
    implicit none
    type            (treeNode                    ), intent(inout)           :: node
    integer                                       , intent(in   )           :: componentType       , massType, weightBy, &
         &                                                                     weightIndex
    double precision                              , intent(in   )           :: positionSpherical(3)
    class           (hotHaloMassDistributionClass)               , pointer  :: hotHalo
    logical                                       , intent(in   ), optional :: haloLoaded
    !GCC$ attributes unused :: haloLoaded, weightIndex
    
    hotHaloMassDistributionDensity=0.0d0
    if (.not.(componentType == componentTypeAll .or. componentType == componentTypeHotHalo                                 )) return
    if (.not.(massType      == massTypeAll      .or. massType      == massTypeBaryonic     .or. massType == massTypeGaseous)) return
    if (.not.(weightBy      == weightByMass                                                                                )) return
    
    hotHalo                        =>     hotHaloMassDistribution        (                         )
    hotHaloMassDistributionDensity =  max(hotHalo                %density(node,positionSpherical(1)),0.0d0)
    return
  end function hotHaloMassDistributionDensity

end module Hot_Halo_Mass_Distributions

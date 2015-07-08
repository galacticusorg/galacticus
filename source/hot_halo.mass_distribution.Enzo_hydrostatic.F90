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

!% An implementation of the hot halo mass distribution class which uses the ``hydrostatic'' profile used by the Enzo simulation code.

  !# <hotHaloMassDistribution name="hotHaloMassDistributionEnzoHydrostatic">
  !#  <description>Provides an implementation of the hot halo mass distribution class which uses the ``hydrostatic'' profile used by the Enzo simulation code.</description>
  !# </hotHaloMassDistribution>
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionEnzoHydrostatic
     !% An implementation of the hot halo mass distribution class which uses the ``hydrostatic'' profile used by the Enzo simulation code.
     private
   contains
     procedure :: densityNormalization  => enzoHydrostaticDensityNormalization
     procedure :: density               => enzoHydrostaticDensity
     procedure :: densityLogSlope       => enzoHydrostaticDensityLogSlope
     procedure :: enclosedMass          => enzoHydrostaticEnclosedMass
     procedure :: radialMoment          => enzoHydrostaticRadialMoment
     procedure :: rotationNormalization => enzoHydrostaticRotationNormalization
  end type hotHaloMassDistributionEnzoHydrostatic

  interface hotHaloMassDistributionEnzoHydrostatic
     !% Constructors for the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution class.
     module procedure enzoHydrostaticDefaultConstructor
  end interface hotHaloMassDistributionEnzoHydrostatic

  type            (treeNode                      ), pointer :: enzoHydrostaticNode
  double precision                                          :: enzoHydrostaticRadiusScale
  class           (hotHaloTemperatureProfileClass), pointer :: enzoHydrostaticNodeHotHaloTemperatureProfile
  !$omp threadprivate(enzoHydrostaticNode,enzoHydrostaticRadiusScale,enzoHydrostaticNodeHotHaloTemperatureProfile)
  
contains

  function enzoHydrostaticDefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution class.
    use Input_Parameters
    implicit none
    type(hotHaloMassDistributionEnzoHydrostatic) :: enzoHydrostaticDefaultConstructor

    return
  end function enzoHydrostaticDefaultConstructor

  double precision function enzoHydrostaticDensityNormalization(self,node)
    !% Return the density normalization in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    use Hot_Halo_Mass_Distributions_Core_Radii
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo
    class           (hotHaloMassDistributionCoreRadiusClass)               , pointer :: hotHaloMassDistributionCoreRadius_
    double precision                                        , parameter              :: toleranceRelative   =1.0d-3
    type            (c_ptr                                 )                         :: parameterPointer
    type            (fgsl_function                         )                         :: integrandFunction
    type            (fgsl_integration_workspace            )                         :: integrationWorkspace
    logical                                                                          :: integrationReset
    double precision                                                                 :: radiusInner                , radiusOuter

    hotHaloMassDistributionCoreRadius_           => hotHaloMassDistributionCoreRadius             (    )
    enzoHydrostaticRadiusScale                   =  hotHaloMassDistributionCoreRadius_%radius     (node)
    enzoHydrostaticNodeHotHaloTemperatureProfile => hotHaloTemperatureProfile                     (    )
    enzoHydrostaticNode                          => node
    hotHalo                                      => node                              %hotHalo    (    )
    if     (                                &
         &   hotHalo%mass       () <= 0.0d0 &
         &  .or.                            &
         &   hotHalo%outerRadius() <= 0.0d0 &
         & ) then
       enzoHydrostaticDensityNormalization=0.0d0
    else
       radiusInner                        =+0.0d0
       radiusOuter                        =+hotHalo%outerRadius()
       integrationReset                   =.true.
       enzoHydrostaticDensityNormalization=+hotHalo%mass       ()                                              &
            &                              /Integrate(                                                         &
            &                                         radiusInner                                            , &
            &                                         radiusOuter                                            , &                  
            &                                         enzoHydrostaticEnclosedMassIntegrand                   , &
            &                                         parameterPointer                                       , &
            &                                         integrandFunction                                      , &
            &                                         integrationWorkspace                                   , &
            &                                         reset                               =integrationReset  , &
            &                                         toleranceAbsolute                   =0.0d+0            , &
            &                                         toleranceRelative                   =toleranceRelative   &
            &                                        )
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end if
    return
  end function enzoHydrostaticDensityNormalization
    
  function enzoHydrostaticEnclosedMassIntegrand(radius,parameterPointer) bind(c)
    !% Integrand used in finding the normalization of the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    implicit none
    real(kind=c_double)        :: enzoHydrostaticEnclosedMassIntegrand
    real(kind=c_double), value :: radius
    type(c_ptr        ), value :: parameterPointer
    real(kind=c_double)        :: temperature                         , radiusEffective

    if (radius <= 0.0d0) then
       enzoHydrostaticEnclosedMassIntegrand=0.0d0
    else
       radiusEffective                     =max(radius,enzoHydrostaticRadiusScale)
       temperature                         =enzoHydrostaticNodeHotHaloTemperatureProfile%temperature(enzoHydrostaticNode,radiusEffective)
       enzoHydrostaticEnclosedMassIntegrand=+4.0d0              &
            &                               *Pi                 &
            &                               /temperature        &
            &                               *radius         **2 &
            &                               /radiusEffective**3
    end if
    return
  end function enzoHydrostaticEnclosedMassIntegrand
  
  double precision function enzoHydrostaticDensity(self,node,radius)
    !% Return the density in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use Hot_Halo_Mass_Distributions_Core_Radii
    use Hot_Halo_Temperature_Profiles
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    double precision                                        , intent(in   )          :: radius
    class           (hotHaloTemperatureProfileClass        )               , pointer :: hotHaloTemperatureProfile_
    class           (hotHaloMassDistributionCoreRadiusClass)               , pointer :: hotHaloMassDistributionCoreRadius_
    double precision                                                                 :: temperature                       , radiusScale, &
         &                                                                              radiusEffective
    
    hotHaloMassDistributionCoreRadius_ => hotHaloMassDistributionCoreRadius             (                    )
    radiusScale                        =  hotHaloMassDistributionCoreRadius_%radius     (node                )
    radiusEffective                    =  max(radius,radiusScale)
    hotHaloTemperatureProfile_         => hotHaloTemperatureProfile                     (                    )
    temperature                        = +hotHaloTemperatureProfile_        %temperature(node,radiusEffective)
    enzoHydrostaticDensity             = +self%densityNormalization                     (node                )     &
         &                               /temperature                                                              &
         &                               /radiusEffective                                                     **3
    return
  end function enzoHydrostaticDensity
  
  double precision function enzoHydrostaticDensityLogSlope(self,node,radius)
    !% Return the logarithmic slope of the density profile in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass
    !% distribution.
    use Hot_Halo_Mass_Distributions_Core_Radii
    use Hot_Halo_Temperature_Profiles
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    double precision                                        , intent(in   )          :: radius
    class           (hotHaloTemperatureProfileClass        )               , pointer :: hotHaloTemperatureProfile_
    class           (hotHaloMassDistributionCoreRadiusClass)               , pointer :: hotHaloMassDistributionCoreRadius_
    double precision                                                                 :: radiusScale

    hotHaloMassDistributionCoreRadius_ => hotHaloMassDistributionCoreRadius             (           )
    radiusScale                        =  hotHaloMassDistributionCoreRadius_%radius     (node       )
    if (radius > radiusScale) then
       hotHaloTemperatureProfile_      => hotHaloTemperatureProfile                     (           )
       enzoHydrostaticDensityLogSlope  = -hotHaloTemperatureProfile_%temperatureLogSlope(node,radius) &
            &                            -3.0d0
    else
       enzoHydrostaticDensityLogSlope  = +0.0d0
    end if
    return
  end function enzoHydrostaticDensityLogSlope
  
  double precision function enzoHydrostaticEnclosedMass(self,node,radius)
    !% Return the enclosed mass in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    use Hot_Halo_Mass_Distributions_Core_Radii
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    double precision                                        , intent(in   )          :: radius 
    class           (hotHaloMassDistributionCoreRadiusClass)               , pointer :: hotHaloMassDistributionCoreRadius_
    double precision                                        , parameter              :: toleranceRelative   =1.0d-3
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo
    type            (c_ptr                                 )                         :: parameterPointer
    type            (fgsl_function                         )                         :: integrandFunction
    type            (fgsl_integration_workspace            )                         :: integrationWorkspace
    logical                                                                          :: integrationReset
    double precision                                                                 :: radiusInner                , radiusOuter
    

    hotHalo => node%hotHalo()
    if (radius > hotHalo%outerRadius()) then
       enzoHydrostaticEnclosedMass=hotHalo%mass()
    else
       hotHaloMassDistributionCoreRadius_           => hotHaloMassDistributionCoreRadius        (    )
       enzoHydrostaticRadiusScale                   =  hotHaloMassDistributionCoreRadius_%radius(node)
       enzoHydrostaticNodeHotHaloTemperatureProfile => hotHaloTemperatureProfile                (    )
       enzoHydrostaticNode                          => node
       radiusInner                                  =  0.0d0
       radiusOuter                                  =  radius
       integrationReset                             =  .true.
       enzoHydrostaticEnclosedMass=+self%densityNormalization(node)                                    &
            &                      *Integrate(                                                         &
            &                                 radiusInner                                            , &
            &                                 radiusOuter                                            , &                  
            &                                 enzoHydrostaticEnclosedMassIntegrand                   , &
            &                                 parameterPointer                                       , &
            &                                 integrandFunction                                      , &
            &                                 integrationWorkspace                                   , &
            &                                 reset                               =integrationReset  , &
            &                                 toleranceAbsolute                   =0.0d+0            , &
            &                                 toleranceRelative                   =toleranceRelative   &
            &                                )
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end if
    return
  end function enzoHydrostaticEnclosedMass
  
  double precision function enzoHydrostaticRadialMoment(self,node,moment,radius)
    !% Return a radial moment of an {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Integration
    use Hot_Halo_Mass_Distributions_Core_Radii
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), pointer :: node
    double precision                                        , intent(in   )          :: moment                                   , radius
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo
    class           (hotHaloTemperatureProfileClass        )               , pointer :: hotHaloTemperatureProfile_
    class           (hotHaloMassDistributionCoreRadiusClass)               , pointer :: hotHaloMassDistributionCoreRadius_
    double precision                                        , parameter              :: toleranceRelative                 =1.0d-3
    type            (c_ptr                                 )                         :: parameterPointer
    type            (fgsl_function                         )                         :: integrandFunction
    type            (fgsl_integration_workspace            )                         :: integrationWorkspace
    logical                                                                          :: integrationReset
    double precision                                                                 :: radiusInner                              , radiusOuter, &
         &                                                                              radiusScale

    hotHaloMassDistributionCoreRadius_ => hotHaloMassDistributionCoreRadius          (    )
    radiusScale                        =  hotHaloMassDistributionCoreRadius_%radius  (node)
    hotHaloTemperatureProfile_         => hotHaloTemperatureProfile                  (    )
    hotHalo                            => node                               %hotHalo(    )
    radiusInner                        =  0.0d0
    radiusOuter                        =  min(                       &
         &                                    radius               , &
         &                                    hotHalo%outerRadius()  &
         &                                   )
    integrationReset                   =  .true.
    enzoHydrostaticRadialMoment=+self%densityNormalization(node)                                    &
         &                      *Integrate(                                                         &
         &                                 radiusInner                                            , &
         &                                 radiusOuter                                            , &
         &                                 enzoHydrostaticRadialMomentIntegrand                   , &
         &                                 parameterPointer                                       , &
         &                                 integrandFunction                                      , &
         &                                 integrationWorkspace                                   , &
         &                                 reset                               =integrationReset  , &
         &                                 toleranceAbsolute                   =0.0d+0            , &
         &                                 toleranceRelative                   =toleranceRelative   &
         &                                )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
    
  contains
    
    function enzoHydrostaticRadialMomentIntegrand(radius,parameterPointer) bind(c)
      !% Integrand used in finding the normalization of the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
      use Hot_Halo_Temperature_Profiles
      implicit none
      real(kind=c_double)        :: enzoHydrostaticRadialMomentIntegrand
      real(kind=c_double), value :: radius
      type(c_ptr        ), value :: parameterPointer
      real(kind=c_double)        :: temperature                         , radiusEffective
      
      if (radius <= 0.0d0) then
         enzoHydrostaticRadialMomentIntegrand=0.0d0
      else
         radiusEffective=max(radius,radiusScale)
         temperature=hotHaloTemperatureProfile_%temperature(enzoHydrostaticNode,radiusEffective)
         enzoHydrostaticRadialMomentIntegrand=+radius         **moment &
              &                               /radiusEffective**3      &
              &                               /temperature
      end if
      return
    end function enzoHydrostaticRadialMomentIntegrand

  end function enzoHydrostaticRadialMoment

  double precision function enzoHydrostaticRotationNormalization(self,node)
    !% Return the relation between specific angular momentum and rotation velocity (assuming a
    !% rotation velocity that is constant in radius) for {\normalfont \ttfamily node}. Specifically, the
    !% normalization, $A$, returned is such that $V_{\mathrm rot} = A J/M$.
    implicit none
    class(hotHaloMassDistributionEnzoHydrostatic), intent(inout)          :: self
    type (treeNode                              ), intent(inout), pointer :: node
    class(nodeComponentHotHalo                  )               , pointer :: hotHalo

    hotHalo                             => node%hotHalo()
    enzoHydrostaticRotationNormalization=                       &
         &  self%radialMoment(node,2.0d0,hotHalo%outerRadius()) &
         & /self%radialMoment(node,3.0d0,hotHalo%outerRadius())
    return
  end function enzoHydrostaticRotationNormalization

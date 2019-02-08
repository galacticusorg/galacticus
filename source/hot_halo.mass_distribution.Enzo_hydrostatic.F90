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

!% An implementation of the hot halo mass distribution class which uses the ``hydrostatic'' profile used by the Enzo simulation code.

  use Hot_Halo_Temperature_Profiles         , only : hotHaloTemperatureProfileClass, hotHaloTemperatureProfile
  use Hot_Halo_Mass_Distributions_Core_Radii

  !# <hotHaloMassDistribution name="hotHaloMassDistributionEnzoHydrostatic">
  !#  <description>Provides an implementation of the hot halo mass distribution class which uses the ``hydrostatic'' profile used by the Enzo simulation code.</description>
  !# </hotHaloMassDistribution>
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionEnzoHydrostatic
     !% An implementation of the hot halo mass distribution class which uses the ``hydrostatic'' profile used by the Enzo simulation code.
     private
     class(hotHaloTemperatureProfileClass        ), pointer :: hotHaloTemperatureProfile_ => null()
     class(hotHaloMassDistributionCoreRadiusClass), pointer :: hotHaloMassDistributionCoreRadius_ => null()
   contains
     !@ <objectMethods>
     !@   <object>hotHaloMassDistributionEnzoHydrostatic</object>
     !@   <objectMethod>
     !@     <method>densityNormalization</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater} node\arginout</arguments>
     !@     <description>Return the normalization of the density profile.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final                                 enzoHydrostaticDestructor
     procedure :: densityNormalization  => enzoHydrostaticDensityNormalization
     procedure :: density               => enzoHydrostaticDensity
     procedure :: densityLogSlope       => enzoHydrostaticDensityLogSlope
     procedure :: enclosedMass          => enzoHydrostaticEnclosedMass
     procedure :: radialMoment          => enzoHydrostaticRadialMoment
     procedure :: rotationNormalization => enzoHydrostaticRotationNormalization
  end type hotHaloMassDistributionEnzoHydrostatic

  interface hotHaloMassDistributionEnzoHydrostatic
     !% Constructors for the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution class.
     module procedure enzoHydrostaticConstructorParameters
     module procedure enzoHydrostaticConstructorInternal
  end interface hotHaloMassDistributionEnzoHydrostatic
  
  type            (treeNode                      ), pointer :: enzoHydrostaticNode
  double precision                                          :: enzoHydrostaticRadiusScale
  class           (hotHaloTemperatureProfileClass), pointer :: enzoHydrostaticNodeHotHaloTemperatureProfile
  !$omp threadprivate(enzoHydrostaticNode,enzoHydrostaticRadiusScale,enzoHydrostaticNodeHotHaloTemperatureProfile)
  
contains

  function enzoHydrostaticConstructorParameters(parameters) result(self)
    !% Default constructor for the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution class.
    use Input_Parameters
    implicit none
    type (hotHaloMassDistributionEnzoHydrostatic)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(hotHaloTemperatureProfileClass        ), pointer       :: hotHaloTemperatureProfile_
    class(hotHaloMassDistributionCoreRadiusClass), pointer       :: hotHaloMassDistributionCoreRadius_

    !# <objectBuilder class="hotHaloTemperatureProfile"         name="hotHaloTemperatureProfile_"         source="parameters"/>
    !# <objectBuilder class="hotHaloMassDistributionCoreRadius" name="hotHaloMassDistributionCoreRadius_" source="parameters"/>
    self=hotHaloMassDistributionEnzoHydrostatic(hotHaloTemperatureProfile_,hotHaloMassDistributionCoreRadius_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="hotHaloTemperatureProfile_"        />
    !# <objectDestructor name="hotHaloMassDistributionCoreRadius_"/>
    return
  end function enzoHydrostaticConstructorParameters

  function enzoHydrostaticConstructorInternal(hotHaloTemperatureProfile_,hotHaloMassDistributionCoreRadius_) result(self)
    !% Generic constructor for the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution class.
    use Galacticus_Error
    use Array_Utilities
    implicit none
    type (hotHaloMassDistributionEnzoHydrostatic)                        :: self
    class(hotHaloTemperatureProfileClass        ), intent(in   ), target :: hotHaloTemperatureProfile_
    class(hotHaloMassDistributionCoreRadiusClass), intent(in   ), target :: hotHaloMassDistributionCoreRadius_
    !# <constructorAssign variables="*hotHaloTemperatureProfile_, *hotHaloMassDistributionCoreRadius_"/>

    return
  end function enzoHydrostaticConstructorInternal

  subroutine enzoHydrostaticDestructor(self)
    !% Destructor for the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution class.
    implicit none
    type(hotHaloMassDistributionEnzoHydrostatic), intent(inout) :: self

    !# <objectDestructor name="self%hotHaloTemperatureProfile_"        />
    !# <objectDestructor name="self%hotHaloMassDistributionCoreRadius_"/>
    return
  end subroutine enzoHydrostaticDestructor

  double precision function enzoHydrostaticDensityNormalization(self,node)
    !% Return the density normalization in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use FGSL                 , only : fgsl_function       , fgsl_integration_workspace
    use Galacticus_Nodes     , only : nodeComponentHotHalo
    use Numerical_Integration
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), target  :: node
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo
    double precision                                        , parameter              :: toleranceRelative    =1.0d-3
    type            (fgsl_function                         )                         :: integrandFunction
    type            (fgsl_integration_workspace            )                         :: integrationWorkspace
    logical                                                                          :: integrationReset
    double precision                                                                 :: radiusInner                 , radiusOuter
    
    enzoHydrostaticRadiusScale                   =  self%hotHaloMassDistributionCoreRadius_%radius                    (node)
    enzoHydrostaticNodeHotHaloTemperatureProfile => self                                   %hotHaloTemperatureProfile_
    enzoHydrostaticNode                          => node
    hotHalo                                      => node                                    %hotHalo                   (    )
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
       enzoHydrostaticDensityNormalization=+hotHalo%mass       ()                                             &
            &                              /Integrate(                                                        &
            &                                                           radiusInner                         , &
            &                                                           radiusOuter                         , &                  
            &                                                           enzoHydrostaticEnclosedMassIntegrand, &
            &                                                           integrandFunction                   , &
            &                                                           integrationWorkspace                , &
            &                                         reset            =integrationReset                    , &
            &                                         toleranceAbsolute=0.0d+0                              , &
            &                                         toleranceRelative=toleranceRelative                     &
            &                                        )
       call Integrate_Done(integrandFunction,integrationWorkspace)
    end if
    return
  end function enzoHydrostaticDensityNormalization
    
  double precision function enzoHydrostaticEnclosedMassIntegrand(radius)
    !% Integrand used in finding the normalization of the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: temperature, radiusEffective
    
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
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    double precision                                        , intent(in   ) :: radius
    double precision                                                        :: temperature    , radiusScale, &
         &                                                                     radiusEffective
    
    radiusScale                        =  self%hotHaloMassDistributionCoreRadius_%radius     (node                )
    radiusEffective                    =  max(radius,radiusScale)
    temperature                        = +self%hotHaloTemperatureProfile_        %temperature(node,radiusEffective)
    enzoHydrostaticDensity             = +self%densityNormalization                          (node                )     &
         &                               /temperature                                                                   &
         &                               /radiusEffective                                                          **3
    return
  end function enzoHydrostaticDensity
  
  double precision function enzoHydrostaticDensityLogSlope(self,node,radius)
    !% Return the logarithmic slope of the density profile in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass
    !% distribution.
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    double precision                                        , intent(in   ) :: radius
    double precision                                                        :: radiusScale

    radiusScale                        =  self%hotHaloMassDistributionCoreRadius_%radius             (node       )
    if (radius > radiusScale) then
       enzoHydrostaticDensityLogSlope  = -self%hotHaloTemperatureProfile_        %temperatureLogSlope(node,radius) &
            &                            -3.0d0
    else
       enzoHydrostaticDensityLogSlope  = +0.0d0
    end if
    return
  end function enzoHydrostaticDensityLogSlope
  
  double precision function enzoHydrostaticEnclosedMass(self,node,radius)
    !% Return the enclosed mass in a {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
    use FGSL                 , only : fgsl_function       , fgsl_integration_workspace
    use Galacticus_Nodes     , only : nodeComponentHotHalo
    use Numerical_Integration
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout)          :: self
    type            (treeNode                              ), intent(inout), target  :: node
    double precision                                        , intent(in   )          :: radius 
    double precision                                        , parameter              :: toleranceRelative   =1.0d-3
    class           (nodeComponentHotHalo                  )               , pointer :: hotHalo
    type            (fgsl_function                         )                         :: integrandFunction
    type            (fgsl_integration_workspace            )                         :: integrationWorkspace
    logical                                                                          :: integrationReset
    double precision                                                                 :: radiusInner                , radiusOuter

    hotHalo => node%hotHalo()
    if (radius > hotHalo%outerRadius()) then
       enzoHydrostaticEnclosedMass=hotHalo%mass()
    else
       enzoHydrostaticRadiusScale                   =  self%hotHaloMassDistributionCoreRadius_%radius(node)
       enzoHydrostaticNodeHotHaloTemperatureProfile => self%hotHaloTemperatureProfile_
       enzoHydrostaticNode                          => node
       radiusInner                                  =  0.0d0
       radiusOuter                                  =  radius
       integrationReset                             =  .true.
       enzoHydrostaticEnclosedMass=+self%densityNormalization(node)                                    &
            &                      *Integrate(                                                         &
            &                                 radiusInner                                            , &
            &                                 radiusOuter                                            , &                  
            &                                 enzoHydrostaticEnclosedMassIntegrand                   , &
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
    use               FGSL                 , only : fgsl_function       , fgsl_integration_workspace
    use               Galacticus_Nodes     , only : nodeComponentHotHalo
    use               Numerical_Integration
    implicit none
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    double precision                                        , intent(in   ) :: moment                                   , radius
    class           (nodeComponentHotHalo                  ), pointer       :: hotHalo
    double precision                                        , parameter     :: toleranceRelative                 =1.0d-3
    type            (fgsl_function                         )                :: integrandFunction
    type            (fgsl_integration_workspace            )                :: integrationWorkspace
    logical                                                                 :: integrationReset
    double precision                                                        :: radiusInner                              , radiusOuter, &
         &                                                                     radiusScale

    radiusScale                        =  self%hotHaloMassDistributionCoreRadius_%radius (node)
    hotHalo                            => node                                   %hotHalo(    )
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
         &                                 integrandFunction                                      , &
         &                                 integrationWorkspace                                   , &
         &                                 reset                               =integrationReset  , &
         &                                 toleranceAbsolute                   =0.0d+0            , &
         &                                 toleranceRelative                   =toleranceRelative   &
         &                                )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
    
  contains
    
    double precision function enzoHydrostaticRadialMomentIntegrand(radius)
      !% Integrand used in finding the normalization of the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution.
      use Hot_Halo_Temperature_Profiles
      implicit none
      double precision, intent(in   ) :: radius
      double precision                :: temperature, radiusEffective
      
      if (radius <= 0.0d0) then
         enzoHydrostaticRadialMomentIntegrand=0.0d0
      else
         radiusEffective                     =max(radius,radiusScale)
         temperature                         =self%hotHaloTemperatureProfile_%temperature(enzoHydrostaticNode,radiusEffective)
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
    !% normalization, $A$, returned is such that $V_\mathrm{rot} = A J/M$.
    use Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class(hotHaloMassDistributionEnzoHydrostatic), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node
    class(nodeComponentHotHalo                  ), pointer       :: hotHalo

    hotHalo                             => node%hotHalo()
    enzoHydrostaticRotationNormalization=                       &
         &  self%radialMoment(node,2.0d0,hotHalo%outerRadius()) &
         & /self%radialMoment(node,3.0d0,hotHalo%outerRadius())
    return
  end function enzoHydrostaticRotationNormalization

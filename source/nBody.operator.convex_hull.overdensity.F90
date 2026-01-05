!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Implements an N-body data operator which computes the overdensity within the convex hull of the particles.
!!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass 
  use :: Cosmology_Functions , only : cosmologyFunctionsClass 
  use :: Linear_Growth       , only : linearGrowthClass

  !![
  <nbodyOperator name="nbodyOperatorConvexHullOverdensity">
   <description>An N-body data operator which computes the overdensity within the convex hull of the particles.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorConvexHullOverdensity
     !!{
     An N-body data operator which computes the overdensity within the convex hull of the particles.
     !!}
     private
     class(cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     class(cosmologyFunctionsClass ), pointer :: cosmologyFunctions_  => null()
     class(linearGrowthClass       ), pointer :: linearGrowth_        => null()
   contains
     final     ::            convexHullOverdensityDestructor
     procedure :: operate => convexHullOverdensityOperate
  end type nbodyOperatorConvexHullOverdensity

  interface nbodyOperatorConvexHullOverdensity
     !!{
     Constructors for the \refClass{nbodyOperatorConvexHullOverdensity} N-body operator class.
     !!}
     module procedure convexHullOverdensityConstructorParameters
     module procedure convexHullOverdensityConstructorInternal
  end interface nbodyOperatorConvexHullOverdensity

contains

  function convexHullOverdensityConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorConvexHullOverdensity} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nbodyOperatorConvexHullOverdensity)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    class(cosmologyFunctionsClass           ), pointer       :: cosmologyFunctions_
    class(linearGrowthClass                 ), pointer       :: linearGrowth_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="linearGrowth"        name="linearGrowth_"        source="parameters"/>
    !!]
    self=nbodyOperatorConvexHullOverdensity(cosmologyParameters_,cosmologyFunctions_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="linearGrowth_"      />
    !!]
    return
  end function convexHullOverdensityConstructorParameters

  function convexHullOverdensityConstructorInternal(cosmologyParameters_,cosmologyFunctions_,linearGrowth_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorConvexHullOverdensity} N-body operator class.
    !!}
    implicit none
    type (nbodyOperatorConvexHullOverdensity)                        :: self
    class(cosmologyParametersClass          ), intent(in   ), target :: cosmologyParameters_
    class(cosmologyFunctionsClass           ), intent(in   ), target :: cosmologyFunctions_
    class(linearGrowthClass                 ), intent(in   ), target :: linearGrowth_
    !![
    <constructorAssign variables="*cosmologyParameters_, *cosmologyFunctions_, *linearGrowth_"/>
    !!]

    return
  end function convexHullOverdensityConstructorInternal

  subroutine convexHullOverdensityDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorConvexHullOverdensity} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorConvexHullOverdensity), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%linearGrowth_"       />
    !!]
    return
  end subroutine convexHullOverdensityDestructor

  subroutine convexHullOverdensityOperate(self,simulations)
    !!{
    Compute the overdensity of the points.
    !!}
    use :: Display, only : displayIndent, displayUnindent, verbosityLevelStandard
    implicit none
    class           (nbodyOperatorConvexHullOverdensity), intent(inout)               :: self
    type            (nBodyData                         ), intent(inout), dimension(:) :: simulations
    integer                                                                           :: i
    double precision                                                                  :: overdensity

    call displayIndent('compute convex hull overdensity',verbosityLevelStandard)
    do i=1,size(simulations)
       overdensity=+(                                                                                                                                             &
            &        +simulations(i)%attributesReal      %value          ('massTotal'       )                                                                     &
            &        /simulations(i)%attributesReal      %value          ('convexHullVolume')                                                                     &
            &        /self          %cosmologyParameters_%densityCritical(                  )                                                                     &
            &        /self          %cosmologyParameters_%OmegaMatter    (                  )                                                                     &
            &        -1.0d0                                                                                                                                       &
            &       )                                                                                                                                             &
            &      /self%linearGrowth_%value(                                                                                                                     &
            &                                expansionFactor=self%cosmologyFunctions_%expansionFactorFromRedshift(                                                &
            &                                                                                                     simulations(i)%attributesReal%value('redshift') &
            &                                                                                                    )                                                &
            &                               )
       call simulations(i)%attributesReal%set           (keyCH        ='convexHullOverdensity',value         =overdensity)
       call simulations(i)%analysis      %writeAttribute(attributeName='convexHullOverdensity',attributeValue=overdensity)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine convexHullOverdensityOperate

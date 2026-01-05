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
An implementation of the hot halo mass distribution class which uses the ``hydrostatic'' profile used by the Enzo simulation code.
!!}

  use :: Hot_Halo_Mass_Distributions_Core_Radii, only : hotHaloMassDistributionCoreRadiusClass
  use :: Hot_Halo_Temperature_Profiles         , only : hotHaloTemperatureProfile             , hotHaloTemperatureProfileClass

  !![
  <hotHaloMassDistribution name="hotHaloMassDistributionEnzoHydrostatic">
   <description>
    A hot halo mass distribution class which adopts a spherically symmetric density profile for the hot halo motivated by the
    ``hydrostatic'' profile available in the \gls{enzo} code. Specifically,
    \begin{equation}
     \rho_\mathrm{hot halo}(r) \propto \left\{ \begin{array}{ll} T^{-1} r^{-1} &amp; \hbox{ if } r &gt; r_\mathrm{core} \\ T^{-1}
     r_\mathrm{core}^{-1} &amp; \hbox{ if } r \le r_\mathrm{core}, \end{array} \right.
    \end{equation}
    where the core radius, $r_\mathrm{core}$, is set using the selected cored profile core radius method (see
    \refPhysics{hotHaloMassDistributionCoreRadius}). The profile is normalized such that the current mass in the
    hot gas profile is contained within the outer radius of the hot halo, $r_\mathrm{hot, outer}$. Note that the \gls{enzo}
    hydrostatic profile does not include this core, but without introducing this the profile mass can be divergent at small
    radii.
   </description>
  </hotHaloMassDistribution>
  !!]
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionEnzoHydrostatic
     !!{
     An implementation of the hot halo mass distribution class which uses the ``hydrostatic'' profile used by the Enzo simulation code.
     !!}
     private
     class(hotHaloTemperatureProfileClass        ), pointer :: hotHaloTemperatureProfile_         => null()
     class(hotHaloMassDistributionCoreRadiusClass), pointer :: hotHaloMassDistributionCoreRadius_ => null()
   contains
     final     ::        enzoHydrostaticDestructor
     procedure :: get => enzoHydrostaticGet
  end type hotHaloMassDistributionEnzoHydrostatic

  interface hotHaloMassDistributionEnzoHydrostatic
     !!{
     Constructors for the \refClass{hotHaloMassDistributionEnzoHydrostatic} hot halo mass distribution class.
     !!}
     module procedure enzoHydrostaticConstructorParameters
     module procedure enzoHydrostaticConstructorInternal
  end interface hotHaloMassDistributionEnzoHydrostatic

contains

  function enzoHydrostaticConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily enzoHydrostatic} hot halo mass distribution class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloMassDistributionEnzoHydrostatic)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(hotHaloTemperatureProfileClass        ), pointer       :: hotHaloTemperatureProfile_
    class(hotHaloMassDistributionCoreRadiusClass), pointer       :: hotHaloMassDistributionCoreRadius_

    !![
    <objectBuilder class="hotHaloTemperatureProfile"         name="hotHaloTemperatureProfile_"         source="parameters"/>
    <objectBuilder class="hotHaloMassDistributionCoreRadius" name="hotHaloMassDistributionCoreRadius_" source="parameters"/>
    !!]
    self=hotHaloMassDistributionEnzoHydrostatic(hotHaloTemperatureProfile_,hotHaloMassDistributionCoreRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloTemperatureProfile_"        />
    <objectDestructor name="hotHaloMassDistributionCoreRadius_"/>
    !!]
    return
  end function enzoHydrostaticConstructorParameters

  function enzoHydrostaticConstructorInternal(hotHaloTemperatureProfile_,hotHaloMassDistributionCoreRadius_) result(self)
    !!{
    Internal constructor for the \refClass{hotHaloMassDistributionEnzoHydrostatic} hot halo mass distribution class.
    !!}
    implicit none
    type (hotHaloMassDistributionEnzoHydrostatic)                        :: self
    class(hotHaloTemperatureProfileClass        ), intent(in   ), target :: hotHaloTemperatureProfile_
    class(hotHaloMassDistributionCoreRadiusClass), intent(in   ), target :: hotHaloMassDistributionCoreRadius_
    !![
    <constructorAssign variables="*hotHaloTemperatureProfile_, *hotHaloMassDistributionCoreRadius_"/>
    !!]

    return
  end function enzoHydrostaticConstructorInternal

  subroutine enzoHydrostaticDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloMassDistributionEnzoHydrostatic} hot halo mass distribution class.
    !!}
    implicit none
    type(hotHaloMassDistributionEnzoHydrostatic), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloTemperatureProfile_"        />
    <objectDestructor name="self%hotHaloMassDistributionCoreRadius_"/>
    !!]
    return
  end subroutine enzoHydrostaticDestructor

  function enzoHydrostaticGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the Enzo hydrostatic hot halo mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo           , treeNode
    use :: Galactic_Structure_Options, only : componentTypeHotHalo           , massTypeGaseous, weightByMass
    use :: Mass_Distributions        , only : massDistributionEnzoHydrostatic
    implicit none
    class           (massDistributionClass                 ), pointer                 :: massDistribution_
    class           (hotHaloMassDistributionEnzoHydrostatic), intent(inout)           :: self
    type            (treeNode                              ), intent(inout)           :: node
    type            (enumerationWeightByType               ), intent(in   ), optional :: weightBy
    integer                                                 , intent(in   ), optional :: weightIndex
    class           (nodeComponentHotHalo                  ), pointer                 :: hotHalo
    double precision                                                                  :: radiusScale      , radiusOuter, &
         &                                                                               mass
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Get properties of the hot halo.
    radiusScale =  self   %hotHaloMassDistributionCoreRadius_%radius     (node)
    hotHalo     => node                                      %hotHalo    (    )
    radiusOuter =  hotHalo                                   %outerRadius(    )
    mass        =  hotHalo                                   %mass       (    )
    ! If outer radius is non-positive return a null profile.
    if (radiusOuter <= 0.0d0 .or. mass <= 0.0d0) return
    ! Create the mass distribution.
    allocate(massDistributionEnzoHydrostatic :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionEnzoHydrostatic)
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionEnzoHydrostatic(                                             &amp;
             &amp;                     mass                 =     mass                , &amp;
             &amp;                     radiusOuter          =     radiusOuter         , &amp;
             &amp;                     radiusScale          =     radiusScale         , &amp;
             &amp;                     truncateAtOuterRadius=     .true.              , &amp;
             &amp;                     componentType        =     componentTypeHotHalo, &amp;
             &amp;                     massType             =     massTypeGaseous       &amp;
             &amp;                    )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function enzoHydrostaticGet

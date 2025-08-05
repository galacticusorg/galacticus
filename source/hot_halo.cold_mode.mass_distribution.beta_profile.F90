!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
An implementation of the hot halo mass distribution class for $\beta$-profile distributions.
!!}

  use :: Hot_Halo_Cold_Mode_Density_Core_Radii, only : hotHaloColdModeCoreRadiiClass
  use :: Mass_Distributions                   , only : massDistributionBetaProfile

  !![
  <hotHaloColdModeMassDistribution name="hotHaloColdModeMassDistributionBetaProfile">
   <description>
    A hot halo cold mode mass distribution class which adopts a spherically symmetric $\beta$-profile density profile for the hot
    halo. Specifically,
    \begin{equation}
     \rho_\mathrm{hot halo}(r) \propto \left[ r^2 + r_\mathrm{core}^2 \right]^{3\beta/2},
    \end{equation}
    where the core radius, $r_\mathrm{core}$, is set using the selected cored profile core radius method (see
    \refPhysics{hotHaloColdModeCoreRadii}). The value of $\beta$ is specified by the {\normalfont
    \ttfamily [beta]} parameter. The profile is normalized such that the current mass in the hot gas profile is contained
    within the outer radius of the hot halo, $r_\mathrm{hot, outer}$.
   </description>
  </hotHaloColdModeMassDistribution>
  !!]
  type, extends(hotHaloColdModeMassDistributionClass) :: hotHaloColdModeMassDistributionBetaProfile
     !!{
     A $\beta$-profile implementation of the hot halo mass distribution class.
     !!}
     private
     double precision                                         :: beta
     class           (hotHaloColdModeCoreRadiiClass), pointer :: hotHaloColdModeCoreRadii_ => null()
   contains
     final     ::        betaProfileDestructor
     procedure :: get => betaProfileGet
  end type hotHaloColdModeMassDistributionBetaProfile

  interface hotHaloColdModeMassDistributionBetaProfile
     !!{
     Constructors for the $\beta$-profile hot halo cold mode mass distribution class.
     !!}
     module procedure betaProfileConstructorParameters
     module procedure betaProfileConstructorInternal
  end interface hotHaloColdModeMassDistributionBetaProfile

contains

  function betaProfileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hotHaloColdModeMassDistributionBetaProfile} hot halo cold mode mass distributionclass which builds the object from a
    parameter set.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List          , Error_Report
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type            (hotHaloColdModeMassDistributionBetaProfile)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    logical                                                     , save          :: initialized              =.false.
    class           (hotHaloColdModeCoreRadiiClass             ), pointer       :: hotHaloColdModeCoreRadii_
    double precision                                                            :: beta

    if (.not.initialized) then
       !$omp critical(betaProfileColdModeInitialized)
       if (.not.initialized) then
          ! Check that required property is gettable.
          if     (                                                                                            &
               &  .not.(                                                                                      &
               &         defaultHotHaloComponent%   massColdIsGettable()                                      &
               &        .and.                                                                                 &
               &         defaultHotHaloComponent%outerRadiusIsGettable()                                      &
               &       )                                                                                      &
               & ) call Error_Report                                                                          &
               & (                                                                                            &
               &  'This method requires that the "massCold" property of the hot halo is gettable.'//          &
               &  Component_List(                                                                             &
               &                 'hotHalo'                                                                 ,  &
               &                  defaultHotHaloComponent%   massColdAttributeMatch(requireGettable=.true.)   &
               &                 .intersection.                                                               &
               &                  defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.)   &
               &                )                                                                          // &
               &  {introspection:location}                                                                    &
               & )
          initialized=.true.
       end if
       !$omp end critical(betaProfileColdModeInitialized)
    end if

    !![
    <inputParameter>
      <name>beta</name>
      <defaultValue>2.0d0/3.0d0</defaultValue>
      <description>The value of $\beta$ in $\beta$-profile hot halo cold mode mass distributions.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="hotHaloColdModeCoreRadii" name="hotHaloColdModeCoreRadii_" source="parameters"/>
    !!]
    self=hotHaloColdModeMassDistributionBetaProfile(beta,hotHaloColdModeCoreRadii_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloColdModeCoreRadii_"/>
    !!]
    return
  end function betaProfileConstructorParameters

  function betaProfileConstructorInternal(beta,hotHaloColdModeCoreRadii_) result(self)
    !!{
    Internal constructor for the \refClass{hotHaloColdModeMassDistributionBetaProfile} hot halo mass distribution class.
    !!}
    implicit none
    type            (hotHaloColdModeMassDistributionBetaProfile)                        :: self
    double precision                                            , intent(in   )         :: beta
    class           (hotHaloColdModeCoreRadiiClass             ), intent(in   ), target :: hotHaloColdModeCoreRadii_
    !![
    <constructorAssign variables="beta, *hotHaloColdModeCoreRadii_"/>
    !!]

    return
  end function betaProfileConstructorInternal

  subroutine betaProfileDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloColdModeMassDistributionBetaProfile} hot halo mass distribution class.
    !!}
    implicit none
    type(hotHaloColdModeMassDistributionBetaProfile), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloColdModeCoreRadii_"/>
    !!]
    return
  end subroutine betaProfileDestructor

  function betaProfileGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the $\beta$-profile hot halo mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo , treeNode
    use :: Galactic_Structure_Options, only : componentTypeColdHalo, massTypeGaseous, weightByMass
    implicit none
    class           (massDistributionClass                     ), pointer                 :: massDistribution_
    class           (hotHaloColdModeMassDistributionBetaProfile), intent(inout)           :: self
    type            (treeNode                                  ), intent(inout)           :: node
    type            (enumerationWeightByType                   ), intent(in   ), optional :: weightBy
    integer                                                     , intent(in   ), optional :: weightIndex
    class           (nodeComponentHotHalo                      ), pointer                 :: hotHalo
    double precision                                                                      :: radiusScale      , radiusOuter, &
         &                                                                                   mass
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Get properties of the hot halo.
    radiusScale =  self   %hotHaloColdModeCoreRadii_%radius     (node)
    hotHalo     => node                             %hotHalo    (    )
    radiusOuter =  hotHalo                          %outerRadius(    )
    mass        =  hotHalo                          %massCold   (    )
    ! If outer radius is non-positive return a null profile.
    if (radiusOuter <= 0.0d0 .or. mass <= 0.0d0) return
    ! Create the beta-profile distribution.
    allocate(massDistributionBetaProfile :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionBetaProfile)
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionBetaProfile(                                                  &amp;
             &amp;                     beta                 =self%beta                 , &amp;
             &amp;                     coreRadius           =     radiusScale          , &amp;
             &amp;                     mass                 =     mass                 , &amp;
             &amp;                     outerRadius          =     radiusOuter          , &amp;
             &amp;                     truncateAtOuterRadius=     .true.               , &amp;
             &amp;                     componentType        =     componentTypeColdHalo, &amp;
             &amp;                     massType             =     massTypeGaseous        &amp;
             &amp;                    )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function betaProfileGet

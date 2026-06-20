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

!!{RST
An implementation of the hot halo mass distribution class for :math:`\beta`-profile distributions.
!!}

  use :: Hot_Halo_Mass_Distributions_Core_Radii, only : hotHaloMassDistributionCoreRadiusClass
  use :: Mass_Distributions                    , only : massDistributionBetaProfile
  use :: Object_Pools                          , only : objectPool

  !![
  <hotHaloMassDistribution name="hotHaloMassDistributionBetaProfile" docformat="rst">
   <description>
   A hot halo mass distribution class which adopts a spherically symmetric :math:`\beta`-profile density profile for the hot halo. Specifically,

   .. math::

      \rho_\mathrm{hot halo}(r) \propto \left[ r^2 + r_\mathrm{core}^2 \right]^{3\beta/2},

   where the core radius, :math:`r_\mathrm{core}`, is set using the selected cored profile core radius method (see :galacticus-class:`hotHaloMassDistributionCoreRadius`). The value of :math:`\beta` is specified by the ``[beta]`` parameter. The profile is normalized such that the current mass in the hot gas profile is contained within the outer radius of the hot halo, :math:`r_\mathrm{hot, outer}`.
   </description>
   <deepCopy>
     <deallocate variables="pool"/>
   </deepCopy>
  </hotHaloMassDistribution>
  !!]
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionBetaProfile
     !!{RST
     A :math:`\beta`-profile implementation of the hot halo mass distribution class.
     !!}
     private
     double precision                                                      :: beta
     class           (hotHaloMassDistributionCoreRadiusClass), pointer     :: hotHaloMassDistributionCoreRadius_ => null()
     type            (objectPool                            ), allocatable :: pool
   contains
     final     ::        betaProfileDestructor
     procedure :: get => betaProfileGet
  end type hotHaloMassDistributionBetaProfile

  interface hotHaloMassDistributionBetaProfile
     !!{RST
     Constructors for the :math:`\beta`-profile hot halo mass distribution class.
     !!}
     module procedure betaProfileConstructorParameters
     module procedure betaProfileConstructorInternal
  end interface hotHaloMassDistributionBetaProfile

contains

  function betaProfileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`hotHaloMassDistributionBetaProfile` hot halo mass distribution class which takes a parameter set as input.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List          , Error_Report
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    use :: Input_Parameters, only : inputParameter          , inputParameters
    implicit none
    type            (hotHaloMassDistributionBetaProfile    )                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    logical                                                 , save          :: initialized                       =.false.
    class           (hotHaloMassDistributionCoreRadiusClass), pointer       :: hotHaloMassDistributionCoreRadius_
    double precision                                                        :: beta

    if (.not.initialized) then
       !$omp critical(betaProfileInitialized)
       if (.not.initialized) then
          ! Check that required property is gettable.
          if     (                                                                                            &
               &  .not.(                                                                                      &
               &         defaultHotHaloComponent%       massIsGettable()                                      &
               &        .and.                                                                                 &
               &         defaultHotHaloComponent%outerRadiusIsGettable()                                      &
               &       )                                                                                      &
               & ) call Error_Report                                                                          &
               & (                                                                                            &
               &  'This method requires that the "mass" property of the hot halo is gettable.'//              &
               &  Component_List(                                                                             &
               &                 'hotHalo'                                                                 ,  &
               &                  defaultHotHaloComponent%       massAttributeMatch(requireGettable=.true.)   &
               &                 .intersection.                                                               &
               &                  defaultHotHaloComponent%outerRadiusAttributeMatch(requireGettable=.true.)   &
               &                )                                                                          // &
               &  {introspection:location}                                                                    &
               & )
          initialized=.true.
       end if
       !$omp end critical(betaProfileInitialized)
    end if

    !![
    <inputParameter docformat="rst">
      <name>beta</name>
      <defaultValue>2.0d0/3.0d0</defaultValue>
      <description>
      The value of :math:`\beta` in :math:`\beta`-profile hot halo mass distributions.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="hotHaloMassDistributionCoreRadius" name="hotHaloMassDistributionCoreRadius_" source="parameters"/>
    !!]
    self=hotHaloMassDistributionBetaProfile(beta,hotHaloMassDistributionCoreRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloMassDistributionCoreRadius_"/>
    !!]
    return
  end function betaProfileConstructorParameters

  function betaProfileConstructorInternal(beta,hotHaloMassDistributionCoreRadius_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`hotHaloMassDistributionBetaProfile` hot halo mass distribution class.
    !!}
    implicit none
    type            (hotHaloMassDistributionBetaProfile    )                        :: self
    double precision                                        , intent(in   )         :: beta
    class           (hotHaloMassDistributionCoreRadiusClass), intent(in   ), target :: hotHaloMassDistributionCoreRadius_
    !![
    <constructorAssign variables="beta, *hotHaloMassDistributionCoreRadius_"/>
    !!]

    return
  end function betaProfileConstructorInternal

  subroutine betaProfileDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`hotHaloMassDistributionBetaProfile` hot halo mass distribution class.
    !!}
    implicit none
    type(hotHaloMassDistributionBetaProfile), intent(inout) :: self

    ! Release any pooled mass distributions so that they are not leaked when this object is destroyed.
    if (allocated(self%pool)) call self%pool%destroy()
    !![
    <objectDestructor name="self%hotHaloMassDistributionCoreRadius_"/>
    !!]
    return
  end subroutine betaProfileDestructor

  function betaProfileGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{RST
    Return the :math:`\beta`-profile hot halo mass distribution for the given ``node``.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo, treeNode
    use :: Galactic_Structure_Options, only : componentTypeHotHalo, massTypeGaseous, weightByMass
    implicit none
    class           (massDistributionClass             ), pointer                 :: massDistribution_
    class           (hotHaloMassDistributionBetaProfile), intent(inout)           :: self
    type            (treeNode                          ), intent(inout)           :: node
    type            (enumerationWeightByType           ), intent(in   ), optional :: weightBy
    integer                                             , intent(in   ), optional :: weightIndex
    class           (nodeComponentHotHalo              ), pointer                 :: hotHalo
    double precision                                                              :: radiusScale      , radiusOuter, &
         &                                                                           mass
    logical                                                                       :: reused
    integer                                                                       :: i
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
    ! Acquire a pool slot, creating the pool itself on first use. If an existing object is available
    ! for re-use "reused" is returned true, otherwise we must create a new object.
    if (.not.allocated(self%pool)) allocate(self%pool)
    call self%pool%acquire(i,reused)
    if (.not.reused) allocate(massDistributionBetaProfile :: self%pool%slots(i)%object_)
    select type (massDistribution__ => self%pool%slots(i)%object_)
    type is (massDistributionBetaProfile)
       if (reused) then
          ! An existing distribution in the pool is available for re-use - update its properties for
          ! this node.
          call massDistribution__%initialize(                                   &
               &                             beta                 =self%beta  , &
               &                             coreRadius           =radiusScale, &
               &                             mass                 =mass       , &
               &                             outerRadius          =radiusOuter, &
               &                             truncateAtOuterRadius=.true.       &
               &                            )
       else
          ! No pool object was available - construct a new distribution.
          !![
          <referenceConstruct object="massDistribution__">
	    <constructor>
              massDistributionBetaProfile(                                                 &amp;
                &amp;                     beta                 =self%beta                , &amp;
                &amp;                     coreRadius           =     radiusScale         , &amp;
                &amp;                     mass                 =     mass                , &amp;
                &amp;                     outerRadius          =     radiusOuter         , &amp;
                &amp;                     truncateAtOuterRadius=     .true.              , &amp;
                &amp;                     componentType        =     componentTypeHotHalo, &amp;
                &amp;                     massType             =     massTypeGaseous       &amp;
                &amp;                    )
	    </constructor>
          </referenceConstruct>
          !!]
       end if
       ! Increment the reference count since we are returning the object to a calling function which
       ! will therefore hold a reference to it.
       !![
       <referenceCountIncrement object="massDistribution__"/>
       !!]
       massDistribution_ => massDistribution__
    end select
    return
  end function betaProfileGet

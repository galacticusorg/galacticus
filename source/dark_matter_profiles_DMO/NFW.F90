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
  An implementation of :cite:t:`navarro_universal_1997` dark matter halo profiles.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Mass_Distributions     , only : massDistributionNFW
  use :: Object_Pools           , only : objectPool

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMONFW" docformat="rst">
   <description>
   A dark matter profile DMO class which builds :galacticus-class:`massDistributionNFW` objects to implement the :term:`NFW` density profile :cite:p:`navarro_universal_1997`, normalized such that the total mass of the :term:`node` is enclosed with the virial radius and with the scale length :math:`r_\mathrm{s} = r_\mathrm{virial}/c` where :math:`c` is the halo concentration (see :galacticus-class:`darkMatterProfileConcentration`).
   </description>
   <deepCopy>
     <deallocate variables="pool"/>
   </deepCopy>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMONFW
     !!{RST
     A dark matter halo profile class implementing :cite:t:`navarro_universal_1997` dark matter halos.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer     :: darkMatterHaloScale_                 => null()
     logical                                        :: velocityDispersionUseSeriesExpansion
     type   (objectPool              ), allocatable :: pool
   contains
     final     ::        nfwDestructor
     procedure :: get => nfwGet
  end type darkMatterProfileDMONFW

  interface darkMatterProfileDMONFW
     !!{RST
     Constructors for the :galacticus-class:`darkMatterProfileDMONFW` dark matter halo profile class.
     !!}
     module procedure nfwConstructorParameters
     module procedure nfwConstructorInternal
  end interface darkMatterProfileDMONFW

contains

  function nfwConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterProfileDMONFW` dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (darkMatterProfileDMONFW )                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    class  (darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    logical                                          :: velocityDispersionUseSeriesExpansion

    !![
    <inputParameter docformat="rst">
      <name>velocityDispersionUseSeriesExpansion</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>
      If ``true``, radial velocity dispersion is computed using series expansion.
      </description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMONFW(velocityDispersionUseSeriesExpansion,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function nfwConstructorParameters

  function nfwConstructorInternal(velocityDispersionUseSeriesExpansion,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterProfileDMONFW` dark matter halo profile class.
    !!}
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none
    type   (darkMatterProfileDMONFW )                        :: self
    class  (darkMatterHaloScaleClass), intent(in   ), target :: darkMatterHaloScale_
    logical                          , intent(in   )         :: velocityDispersionUseSeriesExpansion
    !![
    <constructorAssign variables="velocityDispersionUseSeriesExpansion,*darkMatterHaloScale_"/>
    !!]

    ! Ensure that the dark matter profile component supports a "scale" property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                        &
         & call Error_Report                                                                                             &
         &      (                                                                                                        &
         &       'NFW dark matter profile requires a dark matter profile component with a gettable "scale" property.'//  &
         &       Component_List(                                                                                         &
         &                      'darkMatterProfile'                                                                   ,  &
         &                      defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)            &
         &                     )                                                                                      // &
         &      {introspection:location}                                                                                 &
         &      )
    return
  end function nfwConstructorInternal

  subroutine nfwDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`darkMatterProfileDMONFW` dark matter halo profile class.
    !!}
    implicit none
    type   (darkMatterProfileDMONFW), intent(inout) :: self

    ! Release any pooled mass distributions (and their attached kinematics
    ! distributions) so that they are not leaked when this profile is destroyed.
    if (allocated(self%pool)) call self%pool%destroy()
    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine nfwDestructor

  function nfwGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{RST
    Return the dark matter mass distribution for the given ``node``.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic   , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo, massTypeDark                   , weightByMass
    use :: Mass_Distributions        , only : massDistributionNFW  , kinematicsDistributionNFW
    implicit none
    class  (massDistributionClass         ), pointer                     :: massDistribution_
    type   (kinematicsDistributionNFW     ), pointer                     :: kinematicsDistribution_
    class  (darkMatterProfileDMONFW       ), intent(inout)               :: self
    type   (treeNode                      ), intent(inout)               :: node
    type   (enumerationWeightByType       ), intent(in   ), optional     :: weightBy
    integer                                , intent(in   ), optional     :: weightIndex
    class  (nodeComponentBasic            ), pointer                     :: basic
    class  (nodeComponentDarkMatterProfile), pointer                     :: darkMatterProfile
    logical                                                              :: reused
    integer                                                              :: i
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Get the components needed to define the NFW profile.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    ! Acquire a pool slot, creating the pool itself on first use. If an existing object is
    ! available for re-use "reused" is returned true, otherwise we must create a new object.
    if (.not.allocated(self%pool)) allocate(self%pool)
    call self%pool%acquire(i,reused)
    if (.not.reused) allocate(massDistributionNFW :: self%pool%slots(i)%object_)
    select type (massDistribution__ => self%pool%slots(i)%object_)
    type is (massDistributionNFW)
       if (reused) then
          ! An existing mass distribution in the pool is available for re-use - update its
          ! properties for this node.
          call massDistribution__%initialize(                                                                          &
               &                              mass         =basic            %mass                             (    ), &
               &                              virialRadius =self             %darkMatterHaloScale_%radiusVirial(node), &
               &                              scaleLength  =darkMatterProfile%scale                            (    )  &
               &                             )
       else
          ! No pool object was available - construct a new mass distribution and attach its
          ! kinematics distribution.
          !![
          <referenceConstruct object="massDistribution__">
	    <constructor>
              massDistributionNFW(                                                                                  &amp;
              &amp;               mass         =basic            %mass                                      (    ), &amp;
              &amp;               virialRadius =self             %darkMatterHaloScale_%radiusVirial         (node), &amp;
              &amp;               scaleLength  =darkMatterProfile%scale                                     (    ), &amp;
              &amp;               componentType=                                       componentTypeDarkHalo      , &amp;
              &amp;               massType     =                                       massTypeDark                 &amp;
              &amp;              )
	    </constructor>
          </referenceConstruct>
          !!]
          allocate(kinematicsDistribution_)
          !![
          <referenceConstruct object="kinematicsDistribution_">
	    <constructor>
              kinematicsDistributionNFW(                                                                 &amp;
	       &amp;                    useSeriesApproximation=self%velocityDispersionUseSeriesExpansion &amp;
	       &amp;                   )
	    </constructor>
          </referenceConstruct>
          !!]
          call massDistribution__%setKinematicsDistribution(kinematicsDistribution_)
          !![
          <objectDestructor name="kinematicsDistribution_"/>
          !!]
       end if
       ! Increment the reference count for the object since we are returning it to a calling
       ! function which will therefore hold a reference to it.
       !![
       <referenceCountIncrement object="massDistribution__"/>
       !!]
       massDistribution_ => massDistribution__
    end select
    return
  end function nfwGet

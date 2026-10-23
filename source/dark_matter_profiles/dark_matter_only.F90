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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  An implementation of non-dark-matter-only dark matter halo profiles which are unchanged from their dark-matter-only counterpart.
  !!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Object_Pools            , only : objectPool

  !![
  <darkMatterProfile name="darkMatterProfileDarkMatterOnly" docformat="rst">
   <description>
   An implementation of non-dark-matter-only dark matter halo profiles that returns the unmodified dark-matter-only profile, i.e. baryonic effects on the dark matter distribution are ignored.
   </description>
   <deepCopy>
     <deallocate variables="pool"/>
   </deepCopy>
  </darkMatterProfile>
  !!]
  type, extends(darkMatterProfileClass) :: darkMatterProfileDarkMatterOnly
     !!{RST
     A class implementing non-dark-matter-only dark matter halo profiles which are unchanged from their dark-matter-only counterpart.
     !!}
     private
     class           (cosmologyParametersClass ), pointer     :: cosmologyParameters_                      => null()
     class           (darkMatterProfileDMOClass), pointer     :: darkMatterProfileDMO_                     => null()
     double precision                                         :: darkMatterFraction
     logical                                                  :: chandrasekharIntegralSuppressExtendedMass =  .true.
     type            (objectPool               ), allocatable :: pool
   contains
     final     ::        darkMatterOnlyDestructor
     procedure :: get => darkMatterOnlyGet
  end type darkMatterProfileDarkMatterOnly

  interface darkMatterProfileDarkMatterOnly
     !!{RST
     Constructors for the :galacticus-class:`darkMatterProfileDarkMatterOnly` non-dark-matter-only dark matter halo profile class.
     !!}
     module procedure darkMatterOnlyConstructorParameters
     module procedure darkMatterOnlyConstructorInternal
  end interface darkMatterProfileDarkMatterOnly

contains

  function darkMatterOnlyConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterProfileDarkMatterOnly` non-dark-matter-only dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (darkMatterProfileDarkMatterOnly)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class(darkMatterProfileDMOClass      ), pointer       :: darkMatterProfileDMO_
    logical                                               :: chandrasekharIntegralSuppressExtendedMass

    !![
    <inputParameter docformat="rst">
      <name>chandrasekharIntegralSuppressExtendedMass</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>
      If true, the Chandrasekhar integral (used to compute dynamical friction) is suppressed by a factor accounting for the finite extent of the perturbing subhalo. If false, no such suppression is applied (restoring the behavior prior to the introduction of this factor).
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=darkMatterProfileDarkMatterOnly(cosmologyParameters_,darkMatterProfileDMO_,chandrasekharIntegralSuppressExtendedMass)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function darkMatterOnlyConstructorParameters

  function darkMatterOnlyConstructorInternal(cosmologyParameters_,darkMatterProfileDMO_,chandrasekharIntegralSuppressExtendedMass) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterProfileDarkMatterOnly` non-dark-matter-only dark matter halo profile class.
    !!}
    implicit none
    type   (darkMatterProfileDarkMatterOnly)                          :: self
    class  (cosmologyParametersClass       ), intent(in   ), target   :: cosmologyParameters_
    class  (darkMatterProfileDMOClass      ), intent(in   ), target   :: darkMatterProfileDMO_
    logical                                 , intent(in   ), optional :: chandrasekharIntegralSuppressExtendedMass
    !![
    <constructorAssign variables="chandrasekharIntegralSuppressExtendedMass, *cosmologyParameters_, *darkMatterProfileDMO_"/>
    !!]

    ! Evaluate the dark matter fraction.
    self%darkMatterFraction=+1.0d0                                   &
         &                  -self%cosmologyParameters_%OmegaBaryon() &
         &                  /self%cosmologyParameters_%OmegaMatter()
    return
  end function darkMatterOnlyConstructorInternal

  subroutine darkMatterOnlyDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`darkMatterProfileDarkMatterOnly` non-dark-matter-only dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDarkMatterOnly), intent(inout) :: self

    ! Release any pooled mass distributions (and their attached kinematics distributions). This
    ! must precede destruction of the dark matter-only profile, as the pooled distributions hold
    ! references to distributions drawn from that profile's own pool.
    if (allocated(self%pool)) call self%pool%destroy()
    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine darkMatterOnlyDestructor

  function darkMatterOnlyGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{RST
    Return the dark matter mass distribution for the given ``node``.
    !!}
    use :: Galactic_Structure_Options, only : weightByMass
    use :: Mass_Distributions        , only : massDistributionSpherical, massDistributionSphericalScaler, kinematicsDistributionSphericalScaler, kinematicsDistributionClass
    use :: Error                     , only : Error_Report
    implicit none
    class  (massDistributionClass                ), pointer                 :: massDistribution_
    type   (kinematicsDistributionSphericalScaler), pointer                 :: kinematicsDistribution_
    class  (kinematicsDistributionClass          ), pointer                 :: kinematicsDistributionDMO, kinematicsDistributionScaler
    class  (darkMatterProfileDarkMatterOnly      ), intent(inout), target   :: self
    type   (treeNode                             ), intent(inout), target   :: node
    type   (enumerationWeightByType              ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    class  (massDistributionClass                ), pointer                 :: massDistributionDMO
    logical                                                                 :: reused
    integer                                                                 :: i
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Get the dark matter-only mass distribution.
    massDistributionDMO => self%darkMatterProfileDMO_%get(node,weightBy,weightIndex)
    if (.not.associated(massDistributionDMO)) return
    kinematicsDistributionDMO => massDistributionDMO%kinematicsDistribution()
    select type (massDistributionDMO)
    class is (massDistributionSpherical)
       ! Acquire a pool slot, creating the pool itself on first use. If an existing object is
       ! available for re-use "reused" is returned true, otherwise we must create a new object.
       if (.not.allocated(self%pool)) allocate(self%pool)
       call self%pool%acquire(i,reused)
       if (.not.reused) allocate(massDistributionSphericalScaler :: self%pool%slots(i)%object_)
       select type (massDistribution__ => self%pool%slots(i)%object_)
       type is (massDistributionSphericalScaler)
          if (reused) then
             ! Re-use an existing scaled mass distribution: re-point it, and its attached kinematics
             ! distribution, to the dark matter-only distribution for this node (releasing the
             ! previous ones - in particular returning the previously-scaled distribution to its own
             ! pool) and clear its stale memoized state. The scaling factors are unchanged.
             call massDistribution__%reinitialize(massDistributionDMO)
             kinematicsDistributionScaler => massDistribution__%kinematicsDistribution()
             select type (kinematicsDistributionScaler)
             type is (kinematicsDistributionSphericalScaler)
                call kinematicsDistributionScaler%reinitialize(kinematicsDistributionDMO)
             class default
                call Error_Report('expected a spherical scaler kinematics distribution'//{introspection:location})
             end select
             !![
             <objectDestructor name="kinematicsDistributionScaler"/>
             !!]
          else
             ! No pool object was available - construct a new scaled mass distribution and attach
             ! a scaled kinematics distribution to it.
             !![
             <referenceConstruct object="massDistribution__">
               <constructor>
                 massDistributionSphericalScaler(                                                                                          &amp;
                   &amp;                         factorScalingLength                      =     1.0d0                                    , &amp;
                   &amp;                         factorScalingMass                        =self%darkMatterFraction                       , &amp;
                   &amp;                         massDistribution_                        =     massDistributionDMO                      , &amp;
                   &amp;                         chandrasekharIntegralSuppressExtendedMass=self%chandrasekharIntegralSuppressExtendedMass  &amp;
                   &amp;                        )
               </constructor>
             </referenceConstruct>
             !!]
             allocate(kinematicsDistribution_)
             !![
             <referenceConstruct object="kinematicsDistribution_">
               <constructor>
                 kinematicsDistributionSphericalScaler(                                                        &amp;
                   &amp;                               factorScalingLength    =     1.0d0                    , &amp;
                   &amp;                               factorScalingMass      =self%darkMatterFraction       , &amp;
                   &amp;                               kinematicsDistribution_=     kinematicsDistributionDMO  &amp;
                   &amp;                              )
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
    class default
       call Error_Report('a spherical mass distribution is required'//{introspection:location})
    end select
    !![
    <objectDestructor name="massDistributionDMO"      />
    <objectDestructor name="kinematicsDistributionDMO"/>
    !!]
    return
  end function darkMatterOnlyGet

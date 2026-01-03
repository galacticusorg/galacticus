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
  An implementation of ``Einasto'' dark matter halo profiles.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOEinasto">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionEinasto} objects to compute the Einasto density profile
    (e.g. \citealt{cardone_spherical_2005}), normalized such that the total mass of the \gls{node} is enclosed with the virial
    radius and with the characteristic length $r_{-2} = r_\mathrm{virial}/c$ where $c$ is the halo concentration (see
    \refPhysics{darkMatterProfileConcentration}). The shape parameter, $\alpha$, is set using the density profile shape method
    (see \refPhysics{darkMatterProfileShape}).
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOEinasto
     !!{
     A dark matter halo profile class implementing ``Einasto'' dark matter halos.
     !!}
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::        einastoDestructor
     procedure :: get => einastoGet
  end type darkMatterProfileDMOEinasto

  interface darkMatterProfileDMOEinasto
     !!{
     Constructors for the \refClass{darkMatterProfileDMOEinasto} dark matter halo profile class.
     !!}
     module procedure einastoConstructorParameters
     module procedure einastoConstructorInternal
  end interface darkMatterProfileDMOEinasto

contains

  function einastoConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOEinasto} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMOEinasto )                :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOEinasto(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function einastoConstructorParameters

  function einastoConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOEinasto} dark matter halo profile class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none
    type (darkMatterProfileDMOEinasto)                        :: self
    class(darkMatterHaloScaleClass   ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]
    
    ! Ensure that the dark matter profile component supports both "scale" and "shape" properties. Since we've been called with
    ! a treeNode to process, it should have been initialized by now.
    if     (                                                                                                            &
         &  .not.(                                                                                                      &
         &         defaultDarkMatterProfileComponent%scaleIsGettable()                                                  &
         &        .and.                                                                                                 &
         &         defaultDarkMatterProfileComponent%shapeIsGettable()                                                  &
         &       )                                                                                                      &
         & ) then
       call Error_Report                                                                                                &
            &        (                                                                                                  &
            &         'Einasto dark matter profile requires a dark matter profile component that supports gettable '//  &
            &         '"scale" and "shape" properties.'                                                             //  &
            &         Component_List(                                                                                   &
            &                        'darkMatterProfile'                                                              , &
            &                         defaultDarkMatterProfileComponent%shapeAttributeMatch(requireGettable=.true.)     &
            &                        .intersection.                                                                     &
            &                         defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)     &
            &                       )                                                                               //  &
            &         {introspection:location}                                                                          &
            &        )
    end if
    return
  end function einastoConstructorInternal

  subroutine einastoDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOEinasto} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOEinasto), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine einastoDestructor

  function einastoGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic     , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo  , massTypeDark                       , weightByMass
    use :: Mass_Distributions        , only : massDistributionEinasto, kinematicsDistributionCollisionless
    implicit none
    class           (massDistributionClass              ), pointer                 :: massDistribution_
    type            (kinematicsDistributionCollisionless), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOEinasto        ), intent(inout)           :: self
    type            (treeNode                           ), intent(inout)           :: node
    type            (enumerationWeightByType            ), intent(in   ), optional :: weightBy
    integer                                              , intent(in   ), optional :: weightIndex
    class           (nodeComponentBasic                 ), pointer                 :: basic
    class           (nodeComponentDarkMatterProfile     ), pointer                 :: darkMatterProfile
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionEinasto :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionEinasto)
       basic             => node%basic            ()
       darkMatterProfile => node%darkMatterProfile()
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionEinasto(                                                                                  &amp;
           &amp;                   mass          =basic            %mass                                     (    ), &amp;
           &amp;                   virialRadius  =self             %darkMatterHaloScale_%radiusVirial        (node), &amp;
           &amp;                   scaleLength   =darkMatterProfile%scale                                    (    ), &amp;
           &amp;                   shapeParameter=darkMatterProfile%shape                                    (    ), &amp;
           &amp;                   componentType=                                       componentTypeDarkHalo      , &amp;
           &amp;                   massType     =                                       massTypeDark                 &amp;
           &amp;                  )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionCollisionless()
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function einastoGet

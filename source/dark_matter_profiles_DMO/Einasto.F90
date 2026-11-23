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
  An implementation of "Einasto" dark matter halo profiles.
  !!}

  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_Shape, only : darkMatterProfileShapeClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOEinasto" docformat="rst">
   <description>
   A dark matter profile DMO class which builds :galacticus-class:`massDistributionEinasto` objects to compute the Einasto density
   profile (e.g. :cite:author:`cardone_spherical_2005` :cite:year:`cardone_spherical_2005`), normalized such that the total mass
   of the :term:`node` is enclosed with the virial radius and with the characteristic length :math:`r_{-2} = r_\mathrm{virial}/c`
   where :math:`c` is the halo concentration (see :galacticus-class:`darkMatterProfileConcentration`). The shape parameter,
   :math:`\alpha`, is that of the node's dark matter profile component, set using the density profile shape method (see
   :galacticus-class:`darkMatterProfileShape`). Alternatively, when this class is constructed internally with a
   :galacticus-class:`darkMatterProfileShape` object, :math:`\alpha` is taken from that object, and the dark matter profile
   component need not have a shape.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOEinasto
     !!{RST
     A dark matter halo profile class implementing "Einasto" dark matter halos.
     !!}
     class(darkMatterHaloScaleClass   ), pointer :: darkMatterHaloScale_    => null()
     class(darkMatterProfileShapeClass), pointer :: darkMatterProfileShape_ => null()
   contains
     final     ::        einastoDestructor
     procedure :: get => einastoGet
  end type darkMatterProfileDMOEinasto

  interface darkMatterProfileDMOEinasto
     !!{RST
     Constructors for the :galacticus-class:`darkMatterProfileDMOEinasto` dark matter halo profile class.
     !!}
     module procedure einastoConstructorParameters
     module procedure einastoConstructorInternal
  end interface darkMatterProfileDMOEinasto

contains

  function einastoConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterProfileDMOEinasto` dark matter halo profile class which takes a parameter set as input.
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

  function einastoConstructorInternal(darkMatterHaloScale_,darkMatterProfileShape_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterProfileDMOEinasto` dark matter halo profile class. If
    ``darkMatterProfileShape_`` is given, the shape parameter is taken from it, rather than from the node's dark matter profile
    component.
    !!}
    implicit none
    type (darkMatterProfileDMOEinasto)                                  :: self
    class(darkMatterHaloScaleClass   ), intent(in   ), target           :: darkMatterHaloScale_
    class(darkMatterProfileShapeClass), intent(in   ), target, optional :: darkMatterProfileShape_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *darkMatterProfileShape_"/>
    !!]

    ! Ensure that the dark matter profile component supports the "scale" property, and the "shape" property unless the shape
    ! is to be taken from a shape object. Since we've been called with a treeNode to process, it should have been initialized
    ! by now.
    if (associated(self%darkMatterProfileShape_)) then
       !![
       <componentPropertyAssert class="darkMatterProfile" properties="scale"       require="gettable"/>
       !!]
    else
       !![
       <componentPropertyAssert class="darkMatterProfile" properties="scale shape" require="gettable"/>
       !!]
    end if
    return
  end function einastoConstructorInternal

  subroutine einastoDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`darkMatterProfileDMOEinasto` dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOEinasto), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"   />
    <objectDestructor name="self%darkMatterProfileShape_"/>
    !!]
    return
  end subroutine einastoDestructor

  function einastoGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{RST
    Return the dark matter mass distribution for the given ``node``.
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
    double precision                                                               :: shapeParameter
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
       if (associated(self%darkMatterProfileShape_)) then
          shapeParameter=self%darkMatterProfileShape_%shape(node)
       else
          shapeParameter=     darkMatterProfile      %shape(    )
       end if
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionEinasto(                                                                                  &amp;
           &amp;                   mass          =basic            %mass                                     (    ), &amp;
           &amp;                   radiusVirial  =self             %darkMatterHaloScale_%radiusVirial        (node), &amp;
           &amp;                   scaleLength   =self                                  %scaleRadiusValidated(node), &amp;
           &amp;                   shapeParameter=                                       shapeParameter            , &amp;
           &amp;                   componentType =                                       componentTypeDarkHalo     , &amp;
           &amp;                   massType      =                                       massTypeDark                &amp;
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

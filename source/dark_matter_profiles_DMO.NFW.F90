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
  An implementation of \cite{navarro_universal_1997} dark matter halo profiles.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMONFW">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionNFW} objects to implement the \gls{nfw} density profile
    \citep{navarro_universal_1997}, normalized such that the total mass of the \gls{node} is enclosed with the virial radius and
    with the scale length $r_\mathrm{s} = r_\mathrm{virial}/c$ where $c$ is the halo concentration (see
    \refPhysics{darkMatterProfileConcentration}).
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMONFW
     !!{
     A dark matter halo profile class implementing \cite{navarro_universal_1997} dark matter halos.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_                 => null()
     logical                                    :: velocityDispersionUseSeriesExpansion
   contains
     final     ::        nfwDestructor
     procedure :: get => nfwGet
  end type darkMatterProfileDMONFW

  interface darkMatterProfileDMONFW
     !!{
     Constructors for the \refClass{darkMatterProfileDMONFW} dark matter halo profile class.
     !!}
     module procedure nfwConstructorParameters
     module procedure nfwConstructorInternal
  end interface darkMatterProfileDMONFW

contains

  function nfwConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMONFW} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (darkMatterProfileDMONFW )                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    class  (darkMatterHaloScaleClass), pointer       :: darkMatterHaloScale_
    logical                                          :: velocityDispersionUseSeriesExpansion

    !![
    <inputParameter>
      <name>velocityDispersionUseSeriesExpansion</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, radial velocity dispersion is computed using series expansion.</description>
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
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMONFW} dark matter halo profile class.
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
    !!{
    Destructor for the \refClass{darkMatterProfileDMONFW} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMONFW), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine nfwDestructor

  function nfwGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic   , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo, massTypeDark                   , weightByMass
    use :: Mass_Distributions        , only : massDistributionNFW  , kinematicsDistributionNFW
    implicit none
    class           (massDistributionClass         ), pointer                 :: massDistribution_
    type            (kinematicsDistributionNFW     ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMONFW       ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout)           :: node
    type            (enumerationWeightByType       ), intent(in   ), optional :: weightBy
    integer                                         , intent(in   ), optional :: weightIndex
    class           (nodeComponentBasic            ), pointer                 :: basic
    class           (nodeComponentDarkMatterProfile), pointer                 :: darkMatterProfile
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionNFW :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionNFW)
       basic             => node%basic            ()
       darkMatterProfile => node%darkMatterProfile()
       !![
       <referenceConstruct object="massDistribution_">
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
    end select
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
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function nfwGet

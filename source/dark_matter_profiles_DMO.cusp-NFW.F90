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
  An implementation of cusp-NFW \citep{delos_cusp-halo_2025} dark matter halo profiles.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOCuspNFW">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionCuspNFW} objects to implement the cusp-NFW density profile
    \citep{delos_cusp-halo_2025}, normalized such that the total mass of the \gls{node} is enclosed with the virial radius and
    with the scale length $r_\mathrm{s}$.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOCuspNFW
     !!{
     A dark matter halo profile class implementing cusp-CuspNFW \citep{delos_cusp-halo_2025} dark matter halos.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_                 => null()
     logical                                    :: velocityDispersionUseSeriesExpansion
     integer                                    :: promptCuspNFWYID                              , promptCuspNFWScaleID                      , &
          &                                        promptCuspNFWDensityID
     double precision                           :: toleranceRelativeVelocityDispersion           , toleranceRelativeVelocityDispersionMaximum
   contains
     final     ::        cuspNFWDestructor
     procedure :: get => cuspNFWGet
  end type darkMatterProfileDMOCuspNFW

  interface darkMatterProfileDMOCuspNFW
     !!{
     Constructors for the \refClass{darkMatterProfileDMOCuspNFW} dark matter halo profile class.
     !!}
     module procedure cuspNFWConstructorParameters
     module procedure cuspNFWConstructorInternal
  end interface darkMatterProfileDMOCuspNFW

contains

  function cuspNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOCuspNFW} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOCuspNFW)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass   ), pointer       :: darkMatterHaloScale_
    logical                                                      :: velocityDispersionUseSeriesExpansion
    double precision                                             :: toleranceRelativeVelocityDispersion , toleranceRelativeVelocityDispersionMaximum

    !![
    <inputParameter>
      <name>toleranceRelativeVelocityDispersion</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the velocity dispersion.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersionMaximum</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum relative tolerance to use in numerical solutions for the velocity dispersion.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityDispersionUseSeriesExpansion</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, radial velocity dispersion is computed using series expansion (but only for the case of $y=0$, i.e. an NFW profile).</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOCuspNFW(velocityDispersionUseSeriesExpansion,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function cuspNFWConstructorParameters

  function cuspNFWConstructorInternal(velocityDispersionUseSeriesExpansion,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOCuspNFW} dark matter halo profile class.
    !!}
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none
    type   (darkMatterProfileDMOCuspNFW)                        :: self
    class  (darkMatterHaloScaleClass   ), intent(in   ), target :: darkMatterHaloScale_
    logical                             , intent(in   )         :: velocityDispersionUseSeriesExpansion
    double precision                    , intent(in   )         :: toleranceRelativeVelocityDispersion , toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="velocityDispersionUseSeriesExpansion, toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, *darkMatterHaloScale_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWY"       id="self%promptCuspNFWYID"       isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWScale"   id="self%promptCuspNFWScaleID"   isEvolvable="no" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="promptCuspNFWDensity" id="self%promptCuspNFWDensityID" isEvolvable="no" isCreator="no"/>
    !!]
    ! Ensure that the dark matter profile component supports a "scale" property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                           &
         & call Error_Report                                                                                                &
         &      (                                                                                                           &
         &       'cuspNFW dark matter profile requires a dark matter profile component with a gettable "scale" property.'// &
         &       Component_List(                                                                                            &
         &                      'darkMatterProfile'                                                                       , &
         &                      defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)               &
         &                     )                                                                                         // &
         &      {introspection:location}                                                                                    &
         &      )
    return
  end function cuspNFWConstructorInternal

  subroutine cuspNFWDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOCuspNFW} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOCuspNFW), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine cuspNFWDestructor

  function cuspNFWGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic         , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo      , massTypeDark                  , weightByMass
    use :: Mass_Distributions        , only : massDistributionCuspNFW    , kinematicsDistributionCuspNFW , massDistributionNFW, kinematicsDistributionNFW, &
         &                                    kinematicsDistributionClass
    implicit none
    class           (massDistributionClass         ), pointer                 :: massDistribution_
    class           (kinematicsDistributionClass   ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOCuspNFW   ), intent(inout)           :: self
    type            (treeNode                      ), intent(inout)           :: node
    type            (enumerationWeightByType       ), intent(in   ), optional :: weightBy
    integer                                         , intent(in   ), optional :: weightIndex
    class           (nodeComponentBasic            ), pointer                 :: basic
    class           (nodeComponentDarkMatterProfile), pointer                 :: darkMatterProfile
    double precision                                                          :: amplitudeCuspScaleFree, radiusScaleCuspNFW, &
         &                                                                       densityScaleCuspNFW
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    basic                  => node             %basic                    (                           )
    darkMatterProfile      => node             %darkMatterProfile        (                           )
    amplitudeCuspScaleFree =  darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspNFWYID      )
    radiusScaleCuspNFW     =  darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspNFWScaleID  )
    densityScaleCuspNFW    =  darkMatterProfile%floatRank0MetaPropertyGet(self%promptCuspNFWDensityID)
    if (amplitudeCuspScaleFree <= 0.0d0) then
       ! For halos with no cusp, use an NFW mass distribution.
       allocate(massDistributionNFW           :: massDistribution_      )
       allocate(kinematicsDistributionNFW     :: kinematicsDistribution_)
    else
       ! For halos with a cusp, use the cusp-NFW distribution.
       allocate(massDistributionCuspNFW       :: massDistribution_      )
       allocate(kinematicsDistributionCuspNFW :: kinematicsDistribution_)
    end if
    select type(massDistribution_)
    type is (massDistributionCuspNFW)
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionCuspNFW(                                             &amp;
           &amp;                   densityNormalization=densityScaleCuspNFW   , &amp;
           &amp;                   radiusScale         =radiusScaleCuspNFW    , &amp;
           &amp;                   y                   =amplitudeCuspScaleFree, &amp;
           &amp;                   componentType       =componentTypeDarkHalo , &amp;
           &amp;                   massType            =massTypeDark            &amp;
           &amp;                  )
	 </constructor>
       </referenceConstruct>
       !!]
       select type (kinematicsDistribution_)
       type is (kinematicsDistributionCuspNFW)
          !![
	  <referenceConstruct object="kinematicsDistribution_">
	    <constructor>
              kinematicsDistributionCuspNFW(                                                                                            &amp;
               &amp;                        toleranceRelativeVelocityDispersion       =self%toleranceRelativeVelocityDispersion       , &amp; 
               &amp;                        toleranceRelativeVelocityDispersionMaximum=self%toleranceRelativeVelocityDispersionMaximum  &amp; 
	       &amp;                       )
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
    type is (massDistributionNFW    )
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionNFW(                                                                                  &amp;
            &amp;              mass         =basic            %mass                                      (    ), &amp;
            &amp;              virialRadius =self             %darkMatterHaloScale_%radiusVirial         (node), &amp;
            &amp;              scaleLength  =darkMatterProfile%scale                                     (    ), &amp;
            &amp;              componentType=                                       componentTypeDarkHalo      , &amp;
            &amp;              massType     =                                       massTypeDark                 &amp;
            &amp;             )
	 </constructor>
       </referenceConstruct>
       !!]
       select type (kinematicsDistribution_)
       type is (kinematicsDistributionNFW)
          !![
	  <referenceConstruct object="kinematicsDistribution_">
	    <constructor>
              kinematicsDistributionNFW(                                                                 &amp;
	       &amp;                    useSeriesApproximation=self%velocityDispersionUseSeriesExpansion &amp;
	       &amp;                   )
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
    end select
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function cuspNFWGet

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
An implementation of the hot halo mass distribution class which uses the model of \cite{patej_simple_2015}.
!!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <hotHaloMassDistribution name="hotHaloMassDistributionPatejLoeb2015">
   <description>Provides an implementation of the hot halo mass distribution class which uses the model of \cite{patej_simple_2015}.</description>
  </hotHaloMassDistribution>
  !!]
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionPatejLoeb2015
     !!{
     An implementation of the hot halo mass distribution class which uses the model of \cite{patej_simple_2015}.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     double precision                                     :: gamma                          , radiusShock
   contains
     final     ::        patejLoeb2015Destructor
     procedure :: get => patejLoeb2015Get
  end type hotHaloMassDistributionPatejLoeb2015

  interface hotHaloMassDistributionPatejLoeb2015
     !!{
     Constructors for the \refClass{hotHaloMassDistributionPatejLoeb2015} hot halo mass distribution class.
     !!}
     module procedure patejLoeb2015ConstructorParameters
     module procedure patejLoeb2015ConstructorInternal
  end interface hotHaloMassDistributionPatejLoeb2015

contains

  function patejLoeb2015ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily patejLoeb2015} hot halo mass distribution class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloMassDistributionPatejLoeb2015)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass           ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    double precision                                                      :: gamma               , radiusShock

    !![
    <inputParameter>
      <name>gamma</name>
      <defaultValue>1.15d0</defaultValue>
      <description>The parameter $\Gamma$ in the \cite{patej_simple_2015} hot halo gas mass distribution model.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusShock</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The shock radius, $s$, (in units of the halo virial radius) in the \cite{patej_simple_2015} hot halo gas mass distribution model.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=hotHaloMassDistributionPatejLoeb2015(gamma,radiusShock,darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function patejLoeb2015ConstructorParameters

  function patejLoeb2015ConstructorInternal(gamma,radiusShock,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{hotHaloMassDistributionPatejLoeb2015} hot halo mass distribution class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent, defaultHotHaloComponent
    implicit none
    type            (hotHaloMassDistributionPatejLoeb2015)                        :: self
    double precision                                      , intent(in   )         :: gamma                         , radiusShock
    class           (darkMatterProfileDMOClass           ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass            ), intent(in   ), target :: darkMatterHaloScale_
    logical                                               , save                  :: initialized           =.false.
    !![
    <constructorAssign variables="gamma, radiusShock, *darkMatterProfileDMO_, *darkMatterHaloScale_"/>
    !!]

    ! Check that required properties are gettable.
    if (.not.initialized) then
       !$omp critical(patejLoeb2015Initialized)
       if (.not.initialized) then
          if     (                                                                                                      &
               &  .not.(                                                                                                &
               &         defaultHotHaloComponent          %       massIsGettable()                                      &
               &        .and.                                                                                           &
               &         defaultHotHaloComponent          %outerRadiusIsGettable()                                      &
               &       )                                                                                                &
               & ) call Error_Report                                                                                    &
               & (                                                                                                      &
               &  'This method requires that the "mass" and "outerRadius" properties of the hot halo are gettable.'//   &
               &  Component_List(                                                                                       &
               &                 'hotHalo'                                                                           ,  &
               &                  defaultHotHaloComponent          %       massAttributeMatch(requireGettable=.true.)   &
               &                 .intersection.                                                                         &
               &                  defaultHotHaloComponent          %outerRadiusAttributeMatch(requireGettable=.true.)   &
               &                )                                                                                    // &
               &  {introspection:location}                                                                              &
               & )
          if     (                                                                                                      &
               &  .not.(                                                                                                &
               &         defaultDarkMatterProfileComponent%      scaleIsGettable()                                      &
               &       )                                                                                                &
               & ) call Error_Report                                                                                    &
               & (                                                                                                      &
               &  'This method requires that the "scale" property of the dark matter profile is gettable.'//            &
               &  Component_List(                                                                                       &
               &                 'darkMatterProfile'                                                                 ,  &
               &                  defaultDarkMatterProfileComponent%      scaleAttributeMatch(requireGettable=.true.)   &
               &                )                                                                                    // &
               &  {introspection:location}                                                                              &
               & )
          initialized=.true.
       end if
       !$omp end critical(patejLoeb2015Initialized)
    end if
    return
  end function patejLoeb2015ConstructorInternal

  subroutine patejLoeb2015Destructor(self)
    !!{
    Destructor for the \refClass{hotHaloMassDistributionPatejLoeb2015} hot halo mass distribution class.
    !!}
    implicit none
    type(hotHaloMassDistributionPatejLoeb2015), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine patejLoeb2015Destructor

  function patejLoeb2015Get(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the \cite{patej_simple_2015} hot halo mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Mass_Distributions        , only : massDistributionPatejLoeb2015
    use :: Galacticus_Nodes          , only : nodeComponentHotHalo, treeNode
    use :: Galactic_Structure_Options, only : componentTypeHotHalo, massTypeGaseous, weightByMass
    implicit none
    class           (massDistributionClass               ), pointer                 :: massDistribution_
    class           (hotHaloMassDistributionPatejLoeb2015), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    type            (enumerationWeightByType             ), intent(in   ), optional :: weightBy
    integer                                               , intent(in   ), optional :: weightIndex
    class           (nodeComponentHotHalo                ), pointer                 :: hotHalo
    class           (massDistributionClass               ), pointer                 :: massDistributionDMO
    double precision                                                                :: radiusOuter        , mass, &
         &                                                                             radiusShock
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Get properties of the hot halo.
    hotHalo     => node   %hotHalo    ()
    radiusOuter =  hotHalo%outerRadius()
    mass        =  hotHalo%mass       ()
    ! If outer radius is non-positive return a null profile.
    if (radiusOuter <= 0.0d0 .or. mass <= 0.0d0) return
    radiusShock         =  +self                      %radiusShock        &
         &                 *self%darkMatterHaloScale_ %radiusVirial(node)
    massDistributionDMO =>  self%darkMatterProfileDMO_%get         (node)
     ! Create the mass distribution.
    allocate(massDistributionPatejLoeb2015 :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionPatejLoeb2015)
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionPatejLoeb2015(                                             &amp;
             &amp;                       mass             =     mass                , &amp;
             &amp;                       gamma            =self%gamma               , &amp;
	     &amp;                       radiusShock      =     radiusShock         , &amp;
             &amp;                       radiusOuter      =     radiusOuter         , &amp;
             &amp;                       massDistribution_=     massDistributionDMO , &amp;
             &amp;                       componentType    =     componentTypeHotHalo, &amp;
             &amp;                       massType         =     massTypeGaseous       &amp;
             &amp;                      )
	 </constructor>
       </referenceConstruct>
       <objectDestructor name="massDistributionDMO"/>
       !!]
    end select
    return
  end function patejLoeb2015Get

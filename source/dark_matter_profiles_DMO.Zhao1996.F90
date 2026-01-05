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
  An implementation of \cite{zhao_analytical_1996} dark matter halo profiles.
  !!}

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOZhao1996">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionZhao1996} objects.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOZhao1996
     !!{
     A dark matter halo profile class implementing \cite{zhao_analytical_1996} dark matter halos.
     !!}
     private
     class           (darkMatterHaloScaleClass  ), pointer :: darkMatterHaloScale_ => null()
     double precision                                      :: alpha                         , beta, &
          &                                                   gamma
   contains   
     !![
     <methods>
       <method method="exponents"     description="Compute the exponents for the density profile."   />
       <method method="scaleRadius"   description="Compute the scale radius for the density profile."/>
       <method method="normalization" description="Compute the normalization of the density profile."/>
     </methods>
     !!]
     final     ::                  zhao1996Destructor
     procedure :: get           => zhao1996Get
     procedure :: exponents     => zhao1996Exponents
     procedure :: scaleRadius   => zhao1996ScaleRadius
     procedure :: normalization => zhao1996Normalization
  end type darkMatterProfileDMOZhao1996

  interface darkMatterProfileDMOZhao1996
     !!{
     Constructors for the \refClass{darkMatterProfileDMOZhao1996} dark matter halo profile class.
     !!}
     module procedure zhao1996ConstructorParameters
     module procedure zhao1996ConstructorInternal
  end interface darkMatterProfileDMOZhao1996

contains

  function zhao1996ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily zhao1996} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOZhao1996)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_
    double precision                                              :: alpha               , beta, &
         &                                                           gamma
    
    !![
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <description>The parameter $\alpha$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <description>The parameter $\beta$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <description>The parameter $\gamma$ of the \cite{zhao_analytical_1996} dark matter density profile.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOZhao1996(alpha,beta,gamma,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function zhao1996ConstructorParameters

  function zhao1996ConstructorInternal(alpha,beta,gamma,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOZhao1996} dark matter halo profile class.
    !!}
    implicit none
    type            (darkMatterProfileDMOZhao1996)                        :: self
    class           (darkMatterHaloScaleClass    ), intent(in   ), target :: darkMatterHaloScale_
    double precision                              , intent(in   )         :: alpha               , beta, &
         &                                                                   gamma
    !![
    <constructorAssign variables="alpha, beta, gamma, *darkMatterHaloScale_"/>
    !!]

    return
  end function zhao1996ConstructorInternal

  function zhao1996Get(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo   , massTypeDark                  , weightByMass
    use :: Mass_Distributions        , only : massDistributionZhao1996, kinematicsDistributionZhao1996
    implicit none
    class           (massDistributionClass          ), pointer                 :: massDistribution_
    type            (kinematicsDistributionZhao1996 ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOZhao1996   ), intent(inout)           :: self
    type            (treeNode                       ), intent(inout)           :: node
    type            (enumerationWeightByType        ), intent(in   ), optional :: weightBy
    integer                                          , intent(in   ), optional :: weightIndex
    double precision                                                           :: alpha                  , beta       , &
         &                                                                        gamma                  , scaleRadius, &
         &                                                                        mass
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionZhao1996 :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionZhao1996)
       mass       =self%normalization(node)
       scaleRadius=self%scaleRadius  (node)
       call self%exponents(node,alpha,beta,gamma)
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionZhao1996(                                                                       &amp;
           &amp;                    mass         =       mass                                            , &amp;
           &amp;                    radiusOuter  =self  %darkMatterHaloScale_%radiusVirial         (node), &amp;
           &amp;                    scaleLength  =       scaleRadius                                     , &amp;
           &amp;                    alpha        =       alpha                                           , &amp;
           &amp;                    beta         =       beta                                            , &amp;
           &amp;                    gamma        =       gamma                                           , &amp;
           &amp;                    componentType=                            componentTypeDarkHalo      , &amp;
           &amp;                    massType     =                            massTypeDark                 &amp;
           &amp;                   )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionZhao1996( &amp;
	  &amp;                       )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function zhao1996Get

  subroutine zhao1996Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOZhao1996} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOZhao1996), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine zhao1996Destructor

  subroutine zhao1996Exponents(self,node,alpha,beta,gamma)
    !!{
    Compute the exponents of the {\normalfont \ttfamily zhao1996} dark matter halo profile.
    !!}
    implicit none
    class           (darkMatterProfileDMOZhao1996  ), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(  out) :: alpha, beta, &
         &                                                             gamma

    alpha=self%alpha
    beta =self%beta
    gamma=self%gamma
    return
  end subroutine zhao1996Exponents
  
  double precision function zhao1996ScaleRadius(self,node)
    !!{
    Compute the scale radius of the {\normalfont \ttfamily zhao1996} dark matter halo profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(darkMatterProfileDMOZhao1996  ), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    darkMatterProfile   => node             %darkMatterProfile()
    zhao1996ScaleRadius =  darkMatterProfile%scale            ()
    return
  end function zhao1996ScaleRadius

  double precision function zhao1996Normalization(self,node)
    !!{
    Compute the mass normalization of the {\normalfont \ttfamily zhao1996} dark matter halo profile.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(darkMatterProfileDMOZhao1996), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    class(nodeComponentBasic          ), pointer       :: basic

    basic                 => node %basic()
    zhao1996Normalization =  basic%mass ()
    return
  end function zhao1996Normalization

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  An implementation of \cite{burkert_structure_1995} dark matter halo profiles.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOBurkert">
   <description>
    A dark matter profile DMO class which implements the \citep{burkert_structure_1995} density profile is used
    \begin{equation}
      \rho_\mathrm{dark matter}(r) \propto \left(1+{r\over r_\mathrm{s}}\right)^{-1} \left(1+[{r\over
      r_\mathrm{s}}]^2\right)^{-1},
    \end{equation}
    normalized such that the total mass of the \gls{node} is enclosed with the virial radius and with the scale length
    $r_\mathrm{s} = r_\mathrm{virial}/c$ where $c$ is the halo concentration (see
    \refPhysics{darkMatterProfileConcentration}). The mass enclosed within radius $r$ is given by
    \begin{equation}
    M(&lt;r) = M_\mathrm{virial} {2 \log(1 + R) + \log(1 + R^2) -2 \tan^{-1}(R) \over 2 \log(1 + c) + \log(1 + c^2) -2 \tan^{-1}(c)},
    \end{equation}
    where $R=r/r_\mathrm{s}$. The associated gravitational potential is
    \begin{equation}
    \Phi(r) = -\mathrm{G} \left(1+{1 \over R}\right) { 2 \tan^{-1}(R) - 2 \log(1 + R) + \log(1 + R^2) \over -2 \tan^{-1}(c) + 2
    \log(1 + c) + \log(1 + c^2) }.
    \end{equation}
    The peak of the rotation curve occurs at $R=3.2446257246042642$ (found by numerical solution), and the Fourier transform of
    the profile, $F(k) = \int_0^c 4 \pi r^2 \exp(-i k r) \rho(r) \mathrm{d} r / k r$ (needed in calculations of clustering
    using the halo model) is given by
    \begin{eqnarray}
      F(k) &amp;=&amp; \left\{2 \exp(-i k) \mathrm{C}_\mathrm{i}(k) - 2 \exp(-i k) \mathrm{C}_\mathrm{i}(k[1 + c]) + (1 + i)
      \left[-i \exp(k) \pi - \exp(k) \mathrm{E}_\mathrm{i}(-k) \right. \right. \nonumber \\
        &amp; &amp; +i \exp(-k) \mathrm{E}_\mathrm{i}(k) + \exp(k) \mathrm{E}_\mathrm{i}(i [i + c] k) - i \exp(-k)
        \mathrm{E}_\mathrm{i}(k [1 + i c ]) + (1 + i) \exp(-i k) \mathrm{S}_\mathrm{i}(k) \nonumber \\
        &amp; &amp; \left. \left. - (1 + i) \exp(-i k) \mathrm{S}_\mathrm{i}(k[1 + c])\right]\right\}/\left[k \left\{-2
        \tan^{-1}(c) + 2 \log(1 + c) + \log(1 + c^2)\right\}\right].
    \end{eqnarray}
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOBurkert
     !!{
     A dark matter halo profile class implementing \cite{burkert_structure_1995} dark matter halos.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::        burkertDestructor
     procedure :: get => burkertGet
  end type darkMatterProfileDMOBurkert

  interface darkMatterProfileDMOBurkert
     !!{
     Constructors for the {\normalfont \ttfamily burkert} dark matter halo profile class.
     !!}
     module procedure burkertConstructorParameters
     module procedure burkertConstructorInternal
  end interface darkMatterProfileDMOBurkert

contains

  function burkertConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMOBurkert)                :: self
    type (inputParameters            ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass   ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOBurkert(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function burkertConstructorParameters

  function burkertConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    !!}
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none
    type (darkMatterProfileDMOBurkert)                        :: self
    class(darkMatterHaloScaleClass   ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    ! Ensure that the dark matter profile component supports a "scale" property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                            &
         & call Error_Report                                                                                                 &
         &      (                                                                                                            &
         &       'Burkert dark matter profile requires a dark matter profile component with a gettable "scale" property.'//  &
         &       Component_List(                                                                                             &
         &                      'darkMatterProfile'                                                                        , &
         &                      defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)                &
         &                     )                                                                                         //  &
         &       {introspection:location}                                                                                    &
         &      )
    return
  end function burkertConstructorInternal

  subroutine burkertDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily burkert} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOBurkert), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine burkertDestructor

  function burkertGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic     , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo  , massTypeDark                  , weightByMass
    use :: Mass_Distributions        , only : massDistributionBurkert, kinematicsDistributionBurkert
    implicit none
    class  (massDistributionClass         ), pointer                 :: massDistribution_
    type   (kinematicsDistributionBurkert ), pointer                 :: kinematicsDistribution_
    class  (darkMatterProfileDMOBurkert   ), intent(inout)           :: self
    type   (treeNode                      ), intent(inout)           :: node
    type   (enumerationWeightByType       ), intent(in   ), optional :: weightBy
    integer                                , intent(in   ), optional :: weightIndex
    class  (nodeComponentBasic            ), pointer                 :: basic
    class  (nodeComponentDarkMatterProfile), pointer                 :: darkMatterProfile
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionBurkert :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionBurkert)
       basic             => node%basic            ()
       darkMatterProfile => node%darkMatterProfile()
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionBurkert(                                                                                  &amp;
           &amp;                   mass         =basic            %mass                                      (    ), &amp;
           &amp;                   radiusOuter  =self             %darkMatterHaloScale_%radiusVirial         (node), &amp;
           &amp;                   scaleLength  =darkMatterProfile%scale                                     (    ), &amp;
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
        kinematicsDistributionBurkert( &amp;
	 &amp;                       )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function burkertGet

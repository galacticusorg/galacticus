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
  An implementation of a dark matter density profile which includes the accretion flow surrounding the halo.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOAccretionFlowDiemerKravtsov2014">
    <description>
       An accretion flow class which models the accretion flow using the fitting function of
       \cite{diemer_dependence_2014}. Specifically, \refClass{massDistributionDiemerKravtsov2014} objects are built with
       parameters chosen using fits to the redshift and $\nu$ dependencies of the fitting parameters $b_\mathrm{e}$ and
       $s_\mathrm{e}$ chosen to match the results of their figure~18.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOAccretionFlowDiemerKravtsov2014
     !!{
     A dark matter halo profile class which implements a dark matter density profile which includes the accretion flow using the
     fitting function of \cite{diemer_dependence_2014}.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (darkMatterProfileDMOClass    ), pointer :: darkMatterProfileDMO_     => null()
     double precision                                         :: b0                                 , s0 , &
          &                                                      bz                                 , sz , &
          &                                                      bnu                                , snu     
   contains
     final     ::        accretionFlowDiemerKravtsov2014Destructor
     procedure :: get => accretionFlowDiemerKravtsov2014Get
  end type darkMatterProfileDMOAccretionFlowDiemerKravtsov2014

  interface darkMatterProfileDMOAccretionFlowDiemerKravtsov2014
     !!{
     Constructors for the \refClass{darkMatterProfileDMOAccretionFlowDiemerKravtsov2014} dark matter halo profile class.
     !!}
     module procedure accretionFlowDiemerKravtsov2014ConstructorParameters
     module procedure accretionFlowDiemerKravtsov2014ConstructorInternal
  end interface darkMatterProfileDMOAccretionFlowDiemerKravtsov2014

contains

  function accretionFlowDiemerKravtsov2014ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOAccretionFlowDiemerKravtsov2014} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOAccretionFlowDiemerKravtsov2014)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass                          ), pointer       :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass                            ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass                           ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                      ), pointer       :: cosmologicalMassVariance_
    double precision                                                                     :: b0                       , s0 , &
         &                                                                                  bz                       , sz , &
         &                                                                                  bnu                      , snu

    !![
    <inputParameter>
      <name>b0</name>
      <source>parameters</source>
      <defaultValue>+1.1250d0</defaultValue>
      <defaultSource>Derived by Andrew Benson by constructing simple functional forms which fit the plots in figure 18 of \cite{diemer_dependence_2014}.</defaultSource>
      <description>The parameter $b_0$ in the fitting function $b(\nu,z)=b_0 (1+z)^{b_z} \nu^{b_\nu}$ for the parameter $b(\nu,z)$ appearing in equation (4) of \cite{diemer_dependence_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>bz</name>
      <source>parameters</source>
      <defaultValue>+0.625d0</defaultValue>
      <defaultSource>Derived by Andrew Benson by constructing simple functional forms which fit the plots in figure 18 of \cite{diemer_dependence_2014}.</defaultSource>
      <description>The parameter $b_z$ in the fitting function $b(\nu,z)=b_0 (1+z)^{b_z} \nu^{b_\nu}$ for the parameter $b(\nu,z)$ appearing in equation (4) of \cite{diemer_dependence_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>bnu</name>
      <source>parameters</source>
      <defaultValue>-0.2250d0</defaultValue>
      <defaultSource>Derived by Andrew Benson by constructing simple functional forms which fit the plots in figure 18 of \cite{diemer_dependence_2014}.</defaultSource>
      <description>The parameter $b_\nu$ in the fitting function $b(\nu,z)=b_0 (1+z)^{b_z} \nu^{b_\nu}$ for the parameter $b(\nu,z)$ appearing in equation (4) of \cite{diemer_dependence_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>s0</name>
      <source>parameters</source>
      <defaultValue>+1.3925d0</defaultValue>
      <defaultSource>Derived by Andrew Benson by constructing simple functional forms which fit the plots in figure 18 of \cite{diemer_dependence_2014}.</defaultSource>
      <description>The parameter $s_0$ in the fitting function $s(\nu,z)=s_0 (1+z)^{s_z} \nu^{s_\nu}$ for the parameter $s(\nu,z)$ appearing in equation (4) of \cite{diemer_dependence_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>sz</name>
      <source>parameters</source>
      <defaultValue>-0.199d0</defaultValue>
      <defaultSource>Derived by Andrew Benson by constructing simple functional forms which fit the plots in figure 18 of \cite{diemer_dependence_2014}.</defaultSource>
      <description>The parameter $s_z$ in the fitting function $s(\nu,z)=s_0 (1+z)^{s_z} \nu^{s_\nu}$ for the parameter $s(\nu,z)$ appearing in equation (4) of \cite{diemer_dependence_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>snu</name>
      <source>parameters</source>
      <defaultValue>+0.0875d0</defaultValue>
      <defaultSource>Derived by Andrew Benson by constructing simple functional forms which fit the plots in figure 18 of \cite{diemer_dependence_2014}.</defaultSource>
      <description>The parameter $s_\nu$ in the fitting function $s(\nu,z)=s_0 (1+z)^{s_z} \nu^{s_\nu}$ for the parameter $s(\nu,z)$ appearing in equation (4) of \cite{diemer_dependence_2014}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    !!]
    self=darkMatterProfileDMOAccretionFlowDiemerKravtsov2014(b0,bz,bnu,s0,sz,snu,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="darkMatterProfileDMO_"    />
    !!]
    return
  end function accretionFlowDiemerKravtsov2014ConstructorParameters

  function accretionFlowDiemerKravtsov2014ConstructorInternal(b0,bz,bnu,s0,sz,snu,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOAccretionFlowDiemerKravtsov2014} dark matter profile class.
    !!}
    implicit none
    type            (darkMatterProfileDMOAccretionFlowDiemerKravtsov2014)                        :: self
    class           (darkMatterProfileDMOClass                          ), intent(in   ), target :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass                            ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass                           ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                      ), intent(in   ), target :: cosmologicalMassVariance_
    double precision                                                     , intent(in   )         :: b0                       , s0 , &
         &                                                                                          bz                       , sz , &
         &                                                                                          bnu                      , snu
    !![
    <constructorAssign variables="b0, bz, bnu, s0, sz, snu, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function accretionFlowDiemerKravtsov2014ConstructorInternal

  subroutine accretionFlowDiemerKravtsov2014Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOAccretionFlowDiemerKravtsov2014} dark matter profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOAccretionFlowDiemerKravtsov2014), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%darkMatterProfileDMO_"    />
    !!]
    return
  end subroutine accretionFlowDiemerKravtsov2014Destructor

  function accretionFlowDiemerKravtsov2014Get(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo      , massTypeDark                          , weightByMass
    use :: Mass_Distributions        , only : massDistributionSpherical  , massDistributionSphericalAccretionFlow, massDistributionDiemerKravtsov2014, kinematicsDistributionCollisionless, &
         &                                    nonAnalyticSolversNumerical
    implicit none
    class           (massDistributionClass                              ), pointer                 :: massDistribution_
    class           (massDistributionClass                              ), pointer                 :: massDistributionAccretionFlow_, massDistributionVirialized_
    type            (kinematicsDistributionCollisionless                ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOAccretionFlowDiemerKravtsov2014), intent(inout)           :: self
    type            (treeNode                                           ), intent(inout)           :: node
    type            (enumerationWeightByType                            ), intent(in   ), optional :: weightBy
    integer                                                              , intent(in   ), optional :: weightIndex
    class           (nodeComponentBasic                                 ), pointer                 :: basic
    double precision                                                                               :: time                          , mass                       , &
         &                                                                                            radius200Mean                 , densityMean                , &
         &                                                                                            nu                            , redshift                   , &
         &                                                                                            b                             , s                          , &
         &                                                                                            peakHeight                    , radiusTransition
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]
    
    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Get the virialized mass distribution.
    massDistributionVirialized_ => self%darkMatterProfileDMO_%get(node)
    select type (massDistributionVirialized_)
    class is (massDistributionSpherical)
       ! Create the accretion flow mass distribution.
       allocate(massDistributionDiemerKravtsov2014 :: massDistributionAccretionFlow_)
       select type(massDistributionAccretionFlow_)
       type is (massDistributionDiemerKravtsov2014)
          ! Extract basic quantities for the halo.
          basic => node %basic()
          time  =  basic%time ()
          mass  =  basic%mass ()
          ! Evaluate the control parameters.
          redshift=+self%cosmologyFunctions_      %redshiftFromExpansionFactor(                                &
               &    self%cosmologyFunctions_%expansionFactor                   (time=time                    ) &
               &                                                              )
          nu      =+self%criticalOverdensity_     %value                       (time=time,mass=mass,node=node) &
               &   /self%cosmologicalMassVariance_%rootVariance                (time=time,mass=mass          )
          ! Evaluate the parameters of the fitting function. These fits were derived by Andrew Benson by constructing simple functional
          ! forms which fit the plots in figure 18 of Diemer & Kravtsov (2014). There is no guarantee that these fits will perform
          ! sensibly outside the range of that plot (and, of course, they are only approximate even within the range of that plot).
          b=+self%b0*(1.0+redshift)**self%bz*nu**self%bnu
          s=+self%s0*(1.0+redshift)**self%sz*nu**self%snu
          ! Find the radius enclosing 200 times the mean density.
          densityMean       =  self                       %cosmologyFunctions_  %matterDensityEpochal  (         time       )
          radius200Mean     =  massDistributionVirialized_                      %radiusEnclosingDensity(+200.0d0*densityMean)
          ! Construct the accretion flow mass distribution. Note that we do not include the background density of the universe
          ! here, as (being uniform) it should have no effect on halo dynamics.
          !![
          <referenceConstruct object="massDistributionAccretionFlow_">
            <constructor>
              massDistributionDiemerKravtsov2014(                                     &amp;
               &amp;                             densityMean   =densityMean         , &amp;
               &amp;                             radius200Mean=radius200Mean        , &amp;
               &amp;                             includeMean  =.false.              , &amp;
               &amp;                             b            =b                    , &amp;
               &amp;                             s            =s                    , &amp;
               &amp;                             componentType=componentTypeDarkHalo, &amp;
               &amp;                             massType     =massTypeDark           &amp;
               &amp;                            )
            </constructor>
          </referenceConstruct>
          !!]
          ! Combine the virialized and accretion flow mass distributions.
          allocate(massDistributionSphericalAccretionFlow :: massDistribution_)
          select type(massDistribution_)
          type is (massDistributionSphericalAccretionFlow)
             ! Compute the transition radius following Diemer & Kravtsov (2014; equation 6).
             peakHeight      =+self%criticalOverdensity_     %value       (time=time,mass=mass) &
                  &           /self%cosmologicalMassVariance_%rootVariance(time=time,mass=mass)
             radiusTransition=+(             &
                  &             +1.90d0      &
                  &             -0.18d0      &
                  &             *peakHeight  &
                  &            )             &
                  &           *radius200Mean
             !![
             <referenceConstruct object="massDistribution_">
             <constructor>
             massDistributionSphericalAccretionFlow(                                                               &amp;
               &amp;                                radiusTransition              =radiusTransition              , &amp;
               &amp;                                nonAnalyticSolver             =nonAnalyticSolversNumerical   , &amp;
               &amp;                                massDistribution_             =massDistributionVirialized_   , &amp;
               &amp;                                massDistributionAccretionFlow_=massDistributionAccretionFlow_, &amp;
               &amp;                                componentType                 =componentTypeDarkHalo         , &amp;
               &amp;                                massType                      =massTypeDark                    &amp;
               &amp;                               )
             </constructor>
             </referenceConstruct>
             !!]
          end select
       end select
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
    <objectDestructor name="massDistributionAccretionFlow_"/>
    <objectDestructor name="massDistributionVirialized_"   />
    <objectDestructor name="kinematicsDistribution_"       />
    !!]
    return
  end function accretionFlowDiemerKravtsov2014Get

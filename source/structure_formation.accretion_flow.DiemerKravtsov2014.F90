!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  An accretion flow class which models the accretion flow using the fitting function of \cite{diemer_dependence_2014}.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMOClass
  
  !![
  <accretionFlows name="accretionFlowsDiemerKravtsov2014">
   <description>
    An accretion flow class which models the accretion flow using the fitting function of
    \cite{diemer_dependence_2014}. Specifically, the density profile of the accretion flow is modeled using their equation~(4),
    along with fits to the redshift and $\nu$ dependencies of the fitting parameters $b_\mathrm{e}$ and $s_\mathrm{e}$ chosen to
    match the results of their figure~18.
   </description>
  </accretionFlows>
  !!]
  type, extends(accretionFlowsClass) :: accretionFlowsDiemerKravtsov2014
     !!{
     An accretion flow class which models the accretion flow using the fitting function of \cite{diemer_dependence_2014}.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (darkMatterProfileDMOClass    ), pointer :: darkMatterProfileDMO_     => null()
     double precision                                         :: b0                                 , s0 , &
          &                                                      bz                                 , sz , &
          &                                                      bnu                                , snu
   contains
     final     ::             diemerKravtsov2014Destructor
     procedure :: density  => diemerKravtsov2014Density
     procedure :: velocity => diemerKravtsov2014Velocity
  end type accretionFlowsDiemerKravtsov2014

  interface accretionFlowsDiemerKravtsov2014
     !!{
     Constructors for the {\normalfont \ttfamily diemerKravtsov2014} accretion flows class.
     !!}
     module procedure diemerKravtsov2014ConstructorParameters
     module procedure diemerKravtsov2014ConstructorInternal
  end interface accretionFlowsDiemerKravtsov2014
  
contains
  
  function diemerKravtsov2014ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily diemerKravtsov2014} accretion flow class that takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (accretionFlowsDiemerKravtsov2014)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass   ), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass        ), pointer       :: criticalOverdensity_
    class           (darkMatterProfileDMOClass       ), pointer       :: darkMatterProfileDMO_
    double precision                                                  :: b0                       , s0 , &
         &                                                               bz                       , sz , &
         &                                                               bnu                      , snu
    
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
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    !!]
    self=accretionFlowsDiemerKravtsov2014(b0,bz,bnu,s0,sz,snu,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="darkMatterProfileDMO_"    />
    !!]
    return
  end function diemerKravtsov2014ConstructorParameters

  function diemerKravtsov2014ConstructorInternal(b0,bz,bnu,s0,sz,snu,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily diemerKravtsov2014} accretion flows class.
    !!}
    implicit none
    type            (accretionFlowsDiemerKravtsov2014)                        :: self
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass   ), intent(in   ), target :: cosmologicalMassVariance_
    class           (criticalOverdensityClass        ), intent(in   ), target :: criticalOverdensity_
    class           (darkMatterProfileDMOClass       ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                  , intent(in   )         :: b0                       , s0 , &
         &                                                                       bz                       , sz , &
         &                                                                       bnu                      , snu
    !![
    <constructorAssign variables="b0, bz, bnu, s0, sz, snu, *cosmologyFunctions_, *cosmologicalMassVariance_, *criticalOverdensity_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function diemerKravtsov2014ConstructorInternal

  subroutine diemerKravtsov2014Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily diemerKravtsov2014} accretion flows class.
    !!}
    implicit none
    type(accretionFlowsDiemerKravtsov2014), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%darkMatterProfileDMO_"    />
    !!]
    return
  end subroutine diemerKravtsov2014Destructor
  
  double precision function diemerKravtsov2014Density(self,node,radius)
    !!{
    Compute the density of the accretion flow at the given radius.
    !!}
    use :: Mass_Distributions, only : massDistributionClass
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    implicit none
    class           (accretionFlowsDiemerKravtsov2014), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    double precision                                  , intent(in   ) :: radius
    class           (nodeComponentBasic              ), pointer       :: basic
    class           (massDistributionClass           ), pointer       :: massDistribution_
    double precision                                                  :: time             , mass       , &
         &                                                               radius200Mean    , densityMean, &
         &                                                               nu               , redshift   , &
         &                                                               b                , s
    
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
    massDistribution_ => self             %darkMatterProfileDMO_%get                   (node                )
    densityMean       =  self             %cosmologyFunctions_  %matterDensityEpochal  (time                )
    radius200Mean     =  massDistribution_                      %radiusEnclosingDensity(+200.0d0*densityMean)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Evaluate equation (4) from Diemer & Kravtsov (2014).
    diemerKravtsov2014Density=+densityMean       &
         &                    *(                 &
         &                      +1.0d0           &
         &                      +b               &
         &                      /(               &
         &                        +radius        &
         &                        /5.0d0         &
         &                        /radius200Mean &
         &                       )**s            &
         &                     )
    return
  end function diemerKravtsov2014Density

  double precision function diemerKravtsov2014Velocity(self,node,radius)
    !!{
    Compute the velocity of the accretion flow at the given radius.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (accretionFlowsDiemerKravtsov2014), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    double precision                                  , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius
    
    diemerKravtsov2014Velocity=0.0d0
    call Error_Report('velocity is unsupported'//{introspection:location})
    return
  end function diemerKravtsov2014Velocity

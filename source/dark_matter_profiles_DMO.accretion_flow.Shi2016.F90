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
  An implementation of a dark matter density profile which includes the accretion flow surrounding the halo.
  !!}

  use :: Cosmology_Functions                      , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field               , only : cosmologicalMassVarianceClass          , criticalOverdensityClass
  use :: Spherical_Collapse_Solvers               , only : sphericalCollapseSolverClass
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass
  use :: Dark_Matter_Halo_Scales                  , only : darkMatterHaloScaleClass  
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOAccretionFlowShi2016">
    <description>
       A dark matter profile class which builds \refClass{massDistributionShi2016} objects to model accretion flows using the
       model of \cite{shi_outer_2016}.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOAccretionFlowShi2016
     !!{
     A dark matter halo profile class which implements a dark matter density profile which includes the accretion flow using the
     2-halo correlation function.
     !!}
     private
     class           (darkMatterHaloScaleClass               ), pointer :: darkMatterHaloScale_                => null()
     class           (darkMatterProfileDMOClass              ), pointer :: darkMatterProfileDMO_               => null()
     class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_                 => null()
     class           (criticalOverdensityClass               ), pointer :: criticalOverdensity_                => null()
     class           (cosmologicalMassVarianceClass          ), pointer :: cosmologicalMassVariance_           => null()
     class           (darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     class           (sphericalCollapseSolverClass           ), pointer :: sphericalCollapseSolver_            => null()
     double precision                                                   :: scaleFactorVelocity
    contains
     final     ::        accretionFlowShi2016Destructor
     procedure :: get => accretionFlowShi2016Get
  end type darkMatterProfileDMOAccretionFlowShi2016

  interface darkMatterProfileDMOAccretionFlowShi2016
     !!{
     Constructors for the \refClass{darkMatterProfileDMOAccretionFlowShi2016} dark matter halo profile class.
     !!}
     module procedure accretionFlowShi2016ConstructorParameters
     module procedure accretionFlowShi2016ConstructorInternal
  end interface darkMatterProfileDMOAccretionFlowShi2016

contains

  function accretionFlowShi2016ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOAccretionFlowShi2016} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOAccretionFlowShi2016)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass               ), pointer       :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass                 ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass                ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass           ), pointer       :: cosmologicalMassVariance_
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    class           (sphericalCollapseSolverClass            ), pointer       :: sphericalCollapseSolver_
    class           (darkMatterHaloMassAccretionHistoryClass ), pointer       :: darkMatterHaloMassAccretionHistory_
    double precision                                                          :: scaleFactorVelocity
    
    !![
    <inputParameter>
      <name>scaleFactorVelocity</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>A scale factor to be applied to inflow velocities.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"/>
    <objectBuilder class="criticalOverdensity"                name="criticalOverdensity_"                source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"           name="cosmologicalMassVariance_"           source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"               name="darkMatterProfileDMO_"               source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"                name="darkMatterHaloScale_"                source="parameters"/>
    <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    <objectBuilder class="sphericalCollapseSolver"            name="sphericalCollapseSolver_"            source="parameters"/>
    !!]
    self=darkMatterProfileDMOAccretionFlowShi2016(scaleFactorVelocity,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterHaloMassAccretionHistory_,sphericalCollapseSolver_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="criticalOverdensity_"               />
    <objectDestructor name="cosmologicalMassVariance_"          />
    <objectDestructor name="darkMatterProfileDMO_"              />
    <objectDestructor name="darkMatterHaloScale_"               />
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="sphericalCollapseSolver_"           />
    !!]
    return
  end function accretionFlowShi2016ConstructorParameters

  function accretionFlowShi2016ConstructorInternal(scaleFactorVelocity,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterHaloMassAccretionHistory_,sphericalCollapseSolver_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOAccretionFlowShi2016} dark matter profile class.
    !!}
    implicit none
    type            (darkMatterProfileDMOAccretionFlowShi2016)                        :: self
    class           (darkMatterProfileDMOClass               ), intent(in   ), target :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass                 ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass                ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass           ), intent(in   ), target :: cosmologicalMassVariance_
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterHaloMassAccretionHistoryClass ), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
    class           (sphericalCollapseSolverClass            ), intent(in   ), target :: sphericalCollapseSolver_
    double precision                                          , intent(in   )         :: scaleFactorVelocity
    !![
    <constructorAssign variables="scaleFactorVelocity, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *darkMatterProfileDMO_, *darkMatterHaloScale_, *sphericalCollapseSolver_, *darkMatterHaloMassAccretionHistory_"/>
    !!]

    return
  end function accretionFlowShi2016ConstructorInternal

  subroutine accretionFlowShi2016Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOAccretionFlowShi2016} dark matter profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOAccretionFlowShi2016), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%cosmologicalMassVariance_"          />
    <objectDestructor name="self%criticalOverdensity_"               />
    <objectDestructor name="self%darkMatterProfileDMO_"              />
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="self%sphericalCollapseSolver_"           />
    !!]
    return
  end subroutine accretionFlowShi2016Destructor

  function accretionFlowShi2016Get(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo      , massTypeDark                          , weightByMass
    use :: Mass_Distributions        , only : massDistributionSpherical  , massDistributionSphericalAccretionFlow, massDistributionShi2016, kinematicsDistributionShi2016, &
         &                                    nonAnalyticSolversNumerical
    use :: Numerical_Ranges          , only : Make_Range                 , rangeTypeLogarithmic
    use :: Tables                    , only : table1D
    implicit none
    class           (massDistributionClass                   ), pointer                 :: massDistribution_
    class           (massDistributionClass                   ), pointer                 :: massDistributionAccretionFlow_  , massDistributionVirialized_
    type            (kinematicsDistributionShi2016           ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOAccretionFlowShi2016), intent(inout)           :: self
    type            (treeNode                                ), intent(inout)           :: node
    type            (enumerationWeightByType                 ), intent(in   ), optional :: weightBy
    integer                                                   , intent(in   ), optional :: weightIndex
    class           (nodeComponentBasic                      ), pointer                 :: basic
    class           (table1D                                 ), allocatable             :: ratioRadiusTurnaroundVirialTable
    double precision                                                                    :: time                            , mass                       , &
         &                                                                                 densityMean                     , radius200Mean              , &
         &                                                                                 peakHeight                      , radiusTransition           , &
         &                                                                                 massAccretionRate               , radiusVirial               , &
         &                                                                                 ratioRadiusTurnaroundVirial
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]
    
    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Combine the virialized and accretion flow mass distributions.
    allocate(massDistributionSphericalAccretionFlow :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalAccretionFlow)
       ! Get the virialized mass distribution.
       massDistributionVirialized_ => self%darkMatterProfileDMO_%get(node)
       select type (massDistributionVirialized_)
       class is (massDistributionSpherical)
          ! Find the radius enclosing 200 times the mean density.
          densityMean       =  self                       %cosmologyFunctions_  %matterDensityEpochal  (         time       )
          radius200Mean     =  massDistributionVirialized_                      %radiusEnclosingDensity(+200.0d0*densityMean)
          ! Compute the transition radius following Diemer & Kravtsov (2014; equation 6).
          peakHeight      =+self%criticalOverdensity_     %value       (time=time,mass=mass) &
               &           /self%cosmologicalMassVariance_%rootVariance(time=time,mass=mass)
          radiusTransition=+(             &
               &             +1.90d0      &
               &             -0.18d0      &
               &             *peakHeight  &
               &            )             &
               &           *radius200Mean
          ! Create the accretion flow mass distribution.
          allocate(massDistributionShi2016 :: massDistributionAccretionFlow_)
          select type(massDistributionAccretionFlow_)
          type is (massDistributionShi2016)
             ! Extract basic quantities for the halo.
             basic             => node                                     %basic            (         )
             time              =  basic                                    %time             (         )
             mass              =  basic                                    %mass             (         )
             radiusVirial      =  self %darkMatterHaloScale_               %radiusVirial     (node     )
             massAccretionRate =  self %darkMatterHaloMassAccretionHistory_%massAccretionRate(node,time)
             call self%sphericalCollapseSolver_%radiusTurnaround(time,tableStore=.false.,radiusTurnaround_=ratioRadiusTurnaroundVirialTable)
             ratioRadiusTurnaroundVirial=ratioRadiusTurnaroundVirialTable%interpolate(time)
             !![
             <referenceConstruct object="massDistributionAccretionFlow_">
               <constructor>
		 massDistributionShi2016(                                                         &amp;
                   &amp;                 mass                       =mass                       , &amp;
                   &amp;                 massAccretionRate          =massAccretionRate          , &amp;
                   &amp;                 radiusVirial               =radiusVirial               , &amp;
                   &amp;                 ratioRadiusTurnaroundVirial=ratioRadiusTurnaroundVirial, &amp;
                   &amp;                 time                       =time                       , &amp;
                   &amp;                 scaleFactorVelocity        =self%scaleFactorVelocity   , &amp;
                   &amp;                 cosmologyFunctions_        =self%cosmologyFunctions_   , &amp;
		   &amp;                 componentType              =componentTypeDarkHalo      , &amp;
 		   &amp;                 massType                   =massTypeDark                 &amp;
	 	   &amp;                )
               </constructor>
             </referenceConstruct>
             !!]
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
    <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionShi2016()"/>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="massDistributionAccretionFlow_"/>
    <objectDestructor name="massDistributionVirialized_"   />
    <objectDestructor name="kinematicsDistribution_"       />
    !!]
    return
  end function accretionFlowShi2016Get

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

  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field     , only : cosmologicalMassVarianceClass   , criticalOverdensityClass
  use :: Correlation_Functions_Two_Point, only : correlationFunctionTwoPointClass
  use :: Dark_Matter_Halo_Biases        , only : darkMatterHaloBiasClass
  use :: Dark_Matter_Halo_Scales        , only : darkMatterHaloScaleClass
  use :: Linear_Growth                  , only : linearGrowthClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOAccretionFlowCorrelationFunction">
    <description>
       An accretion flow class which models the accretion flow using the 2-halo correlation function by building
       \refClass{massDistributionCorrelationFunction} objects.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOAccretionFlowCorrelationFunction
     !!{
     A dark matter halo profile class which implements a dark matter density profile which includes the accretion flow using the
     2-halo correlation function.
     !!}
     private
     class           (darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_         => null()
     class           (darkMatterProfileDMOClass       ), pointer :: darkMatterProfileDMO_        => null()
     class           (cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_          => null()
     class           (criticalOverdensityClass        ), pointer :: criticalOverdensity_         => null()
     class           (cosmologicalMassVarianceClass   ), pointer :: cosmologicalMassVariance_    => null()
     class           (correlationFunctionTwoPointClass), pointer :: correlationFunctionTwoPoint_ => null()
     class           (darkMatterHaloBiasClass         ), pointer :: darkMatterHaloBias_          => null()
     class           (linearGrowthClass               ), pointer :: linearGrowth_                => null()
     double precision                                            :: scaleFactorVelocity
    contains
     final     ::        accretionFlowCorrelationFunctionDestructor
     procedure :: get => accretionFlowCorrelationFunctionGet
  end type darkMatterProfileDMOAccretionFlowCorrelationFunction

  interface darkMatterProfileDMOAccretionFlowCorrelationFunction
     !!{
     Constructors for the \refClass{darkMatterProfileDMOAccretionFlowCorrelationFunction} dark matter halo profile class.
     !!}
     module procedure accretionFlowCorrelationFunctionConstructorParameters
     module procedure accretionFlowCorrelationFunctionConstructorInternal
  end interface darkMatterProfileDMOAccretionFlowCorrelationFunction

contains

  function accretionFlowCorrelationFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOAccretionFlowCorrelationFunction} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOAccretionFlowCorrelationFunction)                :: self
    type            (inputParameters                                     ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass                           ), pointer       :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass                             ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass                            ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                       ), pointer       :: cosmologicalMassVariance_
    class           (correlationFunctionTwoPointClass                    ), pointer       :: correlationFunctionTwoPoint_
    class           (darkMatterHaloBiasClass                             ), pointer       :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass                            ), pointer       :: darkMatterHaloScale_
    class           (linearGrowthClass                                   ), pointer       :: linearGrowth_
    double precision :: scaleFactorVelocity
    
    !![
    <inputParameter>
      <name>scaleFactorVelocity</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>A scale factor to be applied to inflow velocities.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters"/>
    <objectBuilder class="criticalOverdensity"         name="criticalOverdensity_"         source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"    name="cosmologicalMassVariance_"    source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"        name="darkMatterProfileDMO_"        source="parameters"/>
    <objectBuilder class="correlationFunctionTwoPoint" name="correlationFunctionTwoPoint_" source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"          name="darkMatterHaloBias_"          source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"         name="darkMatterHaloScale_"         source="parameters"/>
    <objectBuilder class="linearGrowth"                name="linearGrowth_"                source="parameters"/>
    !!]
    self=darkMatterProfileDMOAccretionFlowCorrelationFunction(scaleFactorVelocity,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,darkMatterHaloBias_,correlationFunctionTwoPoint_,darkMatterProfileDMO_,darkMatterHaloScale_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"         />
    <objectDestructor name="criticalOverdensity_"        />
    <objectDestructor name="cosmologicalMassVariance_"   />
    <objectDestructor name="correlationFunctionTwoPoint_"/>
    <objectDestructor name="darkMatterHaloBias_"         />
    <objectDestructor name="darkMatterProfileDMO_"       />
    <objectDestructor name="darkMatterHaloScale_"        />
    <objectDestructor name="linearGrowth_"               />
    !!]
    return
  end function accretionFlowCorrelationFunctionConstructorParameters

  function accretionFlowCorrelationFunctionConstructorInternal(scaleFactorVelocity,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,darkMatterHaloBias_,correlationFunctionTwoPoint_,darkMatterProfileDMO_,darkMatterHaloScale_,linearGrowth_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOAccretionFlowCorrelationFunction} dark matter profile class.
    !!}
    implicit none
    type            (darkMatterProfileDMOAccretionFlowCorrelationFunction)                        :: self
    class           (darkMatterProfileDMOClass                           ), intent(in   ), target :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass                             ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass                            ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                       ), intent(in   ), target :: cosmologicalMassVariance_
    class           (correlationFunctionTwoPointClass                    ), intent(in   ), target :: correlationFunctionTwoPoint_
    class           (darkMatterHaloBiasClass                             ), intent(in   ), target :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass                            ), intent(in   ), target :: darkMatterHaloScale_
    class           (linearGrowthClass                                   ), intent(in   ), target :: linearGrowth_
    double precision, intent(in   ) :: scaleFactorVelocity
    !![
    <constructorAssign variables="scaleFactorVelocity, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *darkMatterHaloBias_, *correlationFunctionTwoPoint_, *darkMatterProfileDMO_, *darkMatterHaloScale_, *linearGrowth_"/>
    !!]

    return
  end function accretionFlowCorrelationFunctionConstructorInternal

  subroutine accretionFlowCorrelationFunctionDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOAccretionFlowCorrelationFunction} dark matter profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOAccretionFlowCorrelationFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%cosmologicalMassVariance_"   />
    <objectDestructor name="self%criticalOverdensity_"        />
    <objectDestructor name="self%darkMatterProfileDMO_"       />
    <objectDestructor name="self%correlationFunctionTwoPoint_"/>
    <objectDestructor name="self%darkMatterHaloBias_"         />
    <objectDestructor name="self%linearGrowth_"               />
    !!]
    return
  end subroutine accretionFlowCorrelationFunctionDestructor

  function accretionFlowCorrelationFunctionGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo      , massTypeDark                          , weightByMass
    use :: Mass_Distributions        , only : massDistributionSpherical  , massDistributionSphericalAccretionFlow, massDistributionCorrelationFunction, kinematicsDistributionLam2013, &
         &                                    nonAnalyticSolversNumerical
    use :: Numerical_Ranges          , only : Make_Range                 , rangeTypeLogarithmic
    implicit none
    class           (massDistributionClass                               ), pointer                     :: massDistribution_
    class           (massDistributionClass                               ), pointer                     :: massDistributionAccretionFlow_          , massDistributionVirialized_
    type            (kinematicsDistributionLam2013                       ), pointer                     :: kinematicsDistribution_
    class           (darkMatterProfileDMOAccretionFlowCorrelationFunction), intent(inout)               :: self
    type            (treeNode                                            ), intent(inout)               :: node
    type            (enumerationWeightByType                             ), intent(in   ), optional     :: weightBy
    integer                                                               , intent(in   ), optional     :: weightIndex
    class           (nodeComponentBasic                                  ), pointer                     :: basic
    double precision                                                      , allocatable  , dimension(:) :: radius                                  , correlationFunction               , &
         &                                                                                                 correlationFunctionVolumeAveraged
    integer                                                               , parameter                   :: countRadiiPerDecade              =10
    double precision                                                      , parameter                   :: factorRadiusMinimum              =10.0d0, factorRadiusMaximum        =10.0d0
    integer                                                                                             :: countRadii                              , i
    double precision                                                                                    :: time                                    , mass                              , &
         &                                                                                                 radiusMinimum                           , radiusMaximum                     , &
         &                                                                                                 densityMean                             , radius200Mean                     , &
         &                                                                                                 peakHeight                              , radiusTransition
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
          allocate(massDistributionCorrelationFunction :: massDistributionAccretionFlow_)
          select type(massDistributionAccretionFlow_)
          type is (massDistributionCorrelationFunction)
             ! Extract basic quantities for the halo.
             basic => node %basic()
             time  =  basic%time ()
             mass  =  basic%mass ()
             ! Build a correlation function.
             radiusMinimum=radius200Mean/factorRadiusMinimum
             radiusMaximum=radius200Mean*factorRadiusMaximum
             countRadii   =int(log10(radiusMaximum/radiusMinimum)*countRadiiPerDecade)+1
             allocate(radius                           (countRadii))
             allocate(correlationFunction              (countRadii))
             allocate(correlationFunctionVolumeAveraged(countRadii))
             radius=Make_Range(radiusMinimum,radiusMaximum,countRadii,rangeTypeLogarithmic)
             do i=1,countRadii
                correlationFunction              (i)=+self%correlationFunctionTwoPoint_%correlation              (     radius(i),time) &
                     &                               *self%darkMatterHaloBias_         %bias                     (node,radius(i)     )
                correlationFunctionVolumeAveraged(i)=+self%correlationFunctionTwoPoint_%correlationVolumeAveraged(     radius(i),time) &
                     &                               *self%darkMatterHaloBias_         %bias                     (node,radius(i)     )
             end do
             !![
             <referenceConstruct object="massDistributionAccretionFlow_">
               <constructor>
		 massDistributionCorrelationFunction(                                              &amp;
		   &amp;                             mass               =mass                    , &amp;
		   &amp;                             time               =time                    , &amp;
		   &amp;                             radius             =radius                  , &amp;
		   &amp;                             correlationFunction=correlationFunction     , &amp;
		   &amp;                             cosmologyFunctions_=self%cosmologyFunctions_, &amp;
		   &amp;                             componentType      =componentTypeDarkHalo   , &amp;
 		   &amp;                             massType           =massTypeDark              &amp;
	 	   &amp;                            )
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
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionLam2013( &amp;
	  &amp;                       massVirial                       =                          mass                                                     , &amp;
	  &amp;                       radiusVirial                     =self%darkMatterHaloScale_%radiusVirial                        (node=node          ), &amp;
	  &amp;                       time                             =                          time                                                     , &amp;
	  &amp;                       overdensityCritical              =self%criticalOverdensity_%value                               (time=time,mass=mass), &amp;
	  &amp;                       rateLinearGrowth                 =self%linearGrowth_       %logarithmicDerivativeExpansionFactor(time=time          ), &amp;
	  &amp;                       scaleFactorVelocity              =self%scaleFactorVelocity                                                           , &amp;
	  &amp;                       radius                           =                          radius                                                   , &amp;
	  &amp;                       correlationFunctionVolumeAveraged=                          correlationFunctionVolumeAveraged                        , &amp;
	  &amp;                       cosmologyFunctions_              =self%cosmologyFunctions_                                                             &amp;
          &amp;                      )
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
  end function accretionFlowCorrelationFunctionGet

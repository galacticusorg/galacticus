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
  An implementation of virial orbits using a loss cone model.
  !!}

  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters           , only : cosmologyParametersClass
  use :: Cosmological_Velocity_Field    , only : cosmologicalVelocityFieldClass
  use :: Cosmological_Density_Field     , only : cosmologicalMassVarianceClass      , cosmologicalMassVariancePeakBackgroundSplit, criticalOverdensityClass, criticalOverdensityPeakBackgroundSplit, &
          &                                      haloEnvironmentNormal
  use :: Halo_Mass_Functions            , only : haloMassFunctionShethTormen
  use :: Dark_Matter_Halo_Biases        , only : darkMatterHaloBiasClass
  use :: Dark_Matter_Halo_Scales        , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO       , only : darkMatterProfileDMOClass
  use :: Linear_Growth                  , only : linearGrowthClass
  use :: Virial_Density_Contrast        , only : virialDensityContrastClass
  use :: Numerical_Interpolation        , only : interpolator
  use :: Correlation_Functions_Two_Point, only : correlationFunctionTwoPointClass
  use :: Galacticus_Nodes               , only : treeNode
  use :: Merger_Tree_Branching          , only : mergerTreeBranchingProbabilityClass

  !![
  <virialOrbit name="virialOrbitLossCone">
   <description>
    A virial orbits class using a loss cone model.
   </description>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitLossCone
     !!{
     A virial orbit class using a loss cone model.
     !!}
     private
     class           (cosmologyFunctionsClass            ), pointer                         :: cosmologyFunctions_             => null()
     class           (cosmologyParametersClass           ), pointer                         :: cosmologyParameters_            => null()
     class           (cosmologicalVelocityFieldClass     ), pointer                         :: cosmologicalVelocityField_      => null()     
     class           (darkMatterHaloBiasClass            ), pointer                         :: darkMatterHaloBias_             => null()
     class           (darkMatterHaloScaleClass           ), pointer                         :: darkMatterHaloScale_            => null()
     class           (darkMatterProfileDMOClass          ), pointer                         :: darkMatterProfileDMO_           => null()
     class           (cosmologicalMassVarianceClass      ), pointer                         :: cosmologicalMassVariance_       => null()
     class           (criticalOverdensityClass           ), pointer                         :: criticalOverdensity_            => null()
     class           (linearGrowthClass                  ), pointer                         :: linearGrowth_                   => null()
     class           (virialDensityContrastClass         ), pointer                         :: virialDensityContrast_          => null()
     class           (correlationFunctionTwoPointClass   ), pointer                         :: correlationFunctionTwoPoint_    => null()
     class           (mergerTreeBranchingProbabilityClass), pointer                         :: mergerTreeBranchingProbability_ => null()
     double precision                                                                       :: velocityMinimum                          , velocityMaximum                     , &
          &                                                                                    massMinimum                              , massMaximum                         , &
          &                                                                                    time                                     , velocityDispersionMultiplier
     integer                                                                                :: countMassesPerDecade                     , countVelocitiesPerUnit
     logical                                                                                :: includeInFlightGrowth
     double precision                                     , allocatable, dimension(:      ) :: mass                                     , velocity
     double precision                                     , allocatable, dimension(:,:    ) :: velocityRadialMeanVirial                 , velocityRadialDispersionVirial      , &
          &                                                                                    velocityTangentialMeanVirial             , velocityTangentialDispersionVirial  , &
          &                                                                                    velocityDistributionPeak                 , velocityTotalRMS
     double precision                                     , allocatable, dimension(:,:,:  ) :: velocityRadialDistributionOrbits         , velocityTangentialDistributionOrbits
     double precision                                     , allocatable, dimension(:,:,:,:) :: velocityDistributionOrbits
     type            (interpolator                       ), allocatable                     :: interpolatorMass                         , interpolatorVelocity
     type            (varying_string                     )                                  :: fileName
     logical                                                                                :: fileRead
     double precision                                                                       :: haloMassFunctionA                        , haloMassFunctionP                   , &
          &                                                                                    haloMassFunctionNormalization
   contains
     !![
     <methods>
      <method description="Tabulate the orbital velocity distribution."                                method="tabulate"    />
      <method description="Compute interpolating factors in the orbital velocity distribution tables." method="interpolants"/>
      <method description="Restore a tabulated solution from file."                                    method="restoreTable"/>
      <method description="Store a tabulated solution to file."                                        method="storeTable"  />
     </methods>
     !!]
     final     ::                                    lossConeDestructor
     procedure :: orbit                           => lossConeOrbit
     procedure :: velocityDistributionFunction    => lossConeVelocityDistributionFunction
     procedure :: densityContrastDefinition       => lossConeDensityContrastDefinition
     procedure :: velocityTangentialMagnitudeMean => lossConeVelocityTangentialMagnitudeMean
     procedure :: velocityTangentialVectorMean    => lossConeVelocityTangentialVectorMean
     procedure :: angularMomentumMagnitudeMean    => lossConeAngularMomentumMagnitudeMean
     procedure :: angularMomentumVectorMean       => lossConeAngularMomentumVectorMean
     procedure :: velocityTotalRootMeanSquared    => lossConeVelocityTotalRootMeanSquared
     procedure :: energyMean                      => lossConeEnergyMean
     procedure :: tabulate                        => lossConeTabulate
     procedure :: interpolants                    => lossConeInterpolants
     procedure :: storeTable                      => lossConeStoreTable
     procedure :: restoreTable                    => lossConeRestoreTable
  end type virialOrbitLossCone

  interface virialOrbitLossCone
     !!{
     Constructors for the \refClass{virialOrbitLossCone} virial orbit class.
     !!}
     module procedure lossConeConstructorParameters
     module procedure lossConeConstructorInternal
  end interface virialOrbitLossCone

  ! Submodule-scope objects used for OpenMP parallelism.
  class(cosmologyFunctionsClass                    ), pointer :: cosmologyFunctions_
  class(cosmologyParametersClass                   ), pointer :: cosmologyParameters_
  class(darkMatterHaloBiasClass                    ), pointer :: darkMatterHaloBias_
  class(darkMatterHaloScaleClass                   ), pointer :: darkMatterHaloScale_
  class(linearGrowthClass                          ), pointer :: linearGrowth_
  class(cosmologicalMassVarianceClass              ), pointer :: cosmologicalMassVariance_
  class(criticalOverdensityClass                   ), pointer :: criticalOverdensity_
  class(correlationFunctionTwoPointClass           ), pointer :: correlationFunctionTwoPoint_
  class(cosmologicalVelocityFieldClass             ), pointer :: cosmologicalVelocityField_
  class(mergerTreeBranchingProbabilityClass        ), pointer :: mergerTreeBranchingProbability_
  type (criticalOverdensityPeakBackgroundSplit     ), pointer :: criticalOverdensityEnvironmental_
  type (cosmologicalMassVariancePeakBackgroundSplit), pointer :: cosmologicalMassVarianceEnvironmental_
  type (haloEnvironmentNormal                      ), pointer :: haloEnvironment_
  type (haloMassFunctionShethTormen                ), pointer :: haloMassFunctionEnvironmental_
  !$omp threadprivate(cosmologyFunctions_,cosmologyParameters_,darkMatterHaloBias_,darkMatterHaloScale_,linearGrowth_,correlationFunctionTwoPoint_,cosmologicalVelocityField_,cosmologicalMassVariance_,criticalOverdensity_,cosmologicalMassVarianceEnvironmental_,criticalOverdensityEnvironmental_,haloEnvironment_,haloMassFunctionEnvironmental_,mergerTreeBranchingProbability_)

  ! Submodule-scope variables used in integrands.
  double precision                    :: massHost_, timeEvaluate_
  type            (treeNode), pointer :: nodeHost_
  !$omp threadprivate(massHost_,nodeHost_,timeEvaluate_)

contains

  function lossConeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialOrbitLossCone} virial orbits class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (virialOrbitLossCone                )                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass           ), pointer       :: cosmologyParameters_
    class           (cosmologicalVelocityFieldClass     ), pointer       :: cosmologicalVelocityField_
    class           (cosmologicalMassVarianceClass      ), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass           ), pointer       :: criticalOverdensity_
    class           (darkMatterHaloBiasClass            ), pointer       :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass          ), pointer       :: darkMatterProfileDMO_
    class           (linearGrowthClass                  ), pointer       :: linearGrowth_
    class           (virialDensityContrastClass         ), pointer       :: virialDensityContrast_
    class           (correlationFunctionTwoPointClass   ), pointer       :: correlationFunctionTwoPoint_
    class           (mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    double precision                                                     :: velocityMinimum                , velocityMaximum             , &
         &                                                                  haloMassFunctionA              , haloMassFunctionP           , &
         &                                                                  haloMassFunctionNormalization  , velocityDispersionMultiplier
    integer                                                              :: countMassesPerDecade           , countVelocitiesPerUnit
    logical                                                              :: includeInFlightGrowth

    !![
    <inputParameter>
      <name>velocityMinimum</name>
      <source>parameters</source>
      <description>The minimum velocity (in units of the host virial velocity) for which to compute velocity distributions.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityMaximum</name>
      <source>parameters</source>
      <description>The maximum velocity (in units of the host virial velocity) for which to compute velocity distributions.</description>
    </inputParameter>
    <inputParameter>
      <name>countVelocitiesPerUnit</name>
      <source>parameters</source>
      <description>The number of points per unit of velocity (in units of the host virial velocity) for which to compute velocity distributions.</description>
    </inputParameter>
    <inputParameter>
      <name>countMassesPerDecade</name>
      <source>parameters</source>
      <description>The number of points per decade of mass for which to compute infall properties.</description>
    </inputParameter>
    <inputParameter>
      <name>includeInFlightGrowth</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, linear growth of the velocity field during the flight of the secondary halo to the primary is included.</description>
    </inputParameter>
    <inputParameter>
      <name>haloMassFunctionA</name>
      <source>parameters</source>
      <defaultValue>0.707d0</defaultValue>
      <description>The parameter $a$ of the \cite{sheth_ellipsoidal_2001} halo mass function used in averaging over environment.</description>
    </inputParameter>
    <inputParameter>
      <name>haloMassFunctionP</name>
      <source>parameters</source>
      <defaultValue>0.300d0</defaultValue>
      <description>The parameter $p$ of the \cite{sheth_ellipsoidal_2001} halo mass function used in averaging over environment.</description>
    </inputParameter>
    <inputParameter>
      <name>haloMassFunctionNormalization</name>
      <source>parameters</source>
      <defaultValue>0.322d0</defaultValue>
      <description>The normalization parameter $A$ of the \cite{sheth_ellipsoidal_2001} halo mass function used in averaging over environment.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityDispersionMultiplier</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>A multiplier applied to the dispersion of the cosmological velocity field.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    <objectBuilder class="cosmologyParameters"            name="cosmologyParameters_"            source="parameters"/>
    <objectBuilder class="cosmologicalVelocityField"      name="cosmologicalVelocityField_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"            source="parameters"/>
    <objectBuilder class="linearGrowth"                   name="linearGrowth_"                   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"             name="darkMatterHaloBias_"             source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"           name="darkMatterProfileDMO_"           source="parameters"/>
    <objectBuilder class="virialDensityContrast"          name="virialDensityContrast_"          source="parameters"/>
    <objectBuilder class="correlationFunctionTwoPoint"    name="correlationFunctionTwoPoint_"    source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    !!]
    self=virialOrbitLossCone(velocityMinimum,velocityMaximum,countVelocitiesPerUnit,countMassesPerDecade,includeInFlightGrowth,haloMassFunctionA,haloMassFunctionP,haloMassFunctionNormalization,velocityDispersionMultiplier,cosmologyFunctions_,cosmologyParameters_,cosmologicalVelocityField_,linearGrowth_,darkMatterHaloBias_,darkMatterHaloScale_,virialDensityContrast_,correlationFunctionTwoPoint_,cosmologicalMassVariance_,criticalOverdensity_,mergerTreeBranchingProbability_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="cosmologyParameters_"           />
    <objectDestructor name="cosmologicalVelocityField_"     />
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="linearGrowth_"                  />
    <objectDestructor name="darkMatterHaloScale_"           />
    <objectDestructor name="darkMatterHaloBias_"            />
    <objectDestructor name="darkMatterProfileDMO_"          />
    <objectDestructor name="virialDensityContrast_"         />
    <objectDestructor name="correlationFunctionTwoPoint_"   />
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    return
  end function lossConeConstructorParameters

  function lossConeConstructorInternal(velocityMinimum,velocityMaximum,countVelocitiesPerUnit,countMassesPerDecade,includeInFlightGrowth,haloMassFunctionA,haloMassFunctionP,haloMassFunctionNormalization,velocityDispersionMultiplier,cosmologyFunctions_,cosmologyParameters_,cosmologicalVelocityField_,linearGrowth_,darkMatterHaloBias_,darkMatterHaloScale_,virialDensityContrast_,correlationFunctionTwoPoint_,cosmologicalMassVariance_,criticalOverdensity_,mergerTreeBranchingProbability_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{virialOrbitLossCone} virial orbits class.
    !!}
    use :: Input_Paths       , only : inputPath   , pathTypeDataDynamic
    use :: ISO_Varying_String, only : operator(//)
    use :: Numerical_Ranges  , only : Make_Range  , rangeTypeLinear
    implicit none
    type            (virialOrbitLossCone                )                        :: self
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass           ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologicalVelocityFieldClass     ), intent(in   ), target :: cosmologicalVelocityField_
    class           (cosmologicalMassVarianceClass      ), intent(in   ), target :: cosmologicalMassVariance_
    class           (criticalOverdensityClass           ), intent(in   ), target :: criticalOverdensity_
    class           (darkMatterHaloBiasClass            ), intent(in   ), target :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass          ), intent(in   ), target :: darkMatterProfileDMO_
    class           (linearGrowthClass                  ), intent(in   ), target :: linearGrowth_
    class           (virialDensityContrastClass         ), intent(in   ), target :: virialDensityContrast_
    class           (correlationFunctionTwoPointClass   ), intent(in   ), target :: correlationFunctionTwoPoint_
    class           (mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    double precision                                     , intent(in   )         :: velocityMinimum                , velocityMaximum             , &
         &                                                                          haloMassFunctionA              , haloMassFunctionP           , &
         &                                                                          haloMassFunctionNormalization  , velocityDispersionMultiplier
    integer                                              , intent(in   )         :: countMassesPerDecade           , countVelocitiesPerUnit
    logical                                              , intent(in   )         :: includeInFlightGrowth
    integer                                                                      :: countVelocities
    !![
    <constructorAssign variables="velocityMinimum, velocityMaximum, countVelocitiesPerUnit, countMassesPerDecade, includeInFlightGrowth, haloMassFunctionA, haloMassFunctionP, haloMassFunctionNormalization, velocityDispersionMultiplier, *cosmologyFunctions_, *cosmologyParameters_, *cosmologicalVelocityField_, *linearGrowth_, *darkMatterHaloBias_, *darkMatterHaloScale_, *virialDensityContrast_, *correlationFunctionTwoPoint_, *cosmologicalMassVariance_, *criticalOverdensity_, *mergerTreeBranchingProbability_, *darkMatterProfileDMO_"/>
    !!]

    ! Set an initial mass range, along with an unphysical initial time (so that retabulation will be forced on the first call).
    self%time       =-huge(0.0d0)
    self%massMinimum=1.0d06
    self%massMaximum=1.0d15
    ! Build the velocity array and interpolator.
    countVelocities=int((self%velocityMaximum-self%velocityMinimum)*dble(self%countVelocitiesPerUnit))+1    
    allocate(self%velocity            (countVelocities))
    allocate(self%interpolatorVelocity                 )
    self%velocity            =Make_Range(self%velocityMinimum,self%velocityMaximum,countVelocities,rangeTypeLinear)
    self%interpolatorVelocity=interpolator(self%velocity)
    ! Build a file name for storing the tabulated solution.
    self%fileName=inputPath(pathTypeDataDynamic)                                                       // &
         &        'darkMatterHalos/'                                                                   // &
         &        self%objectType      (                                                              )// &
         &        '_'                                                                                  // &
         &        self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &        '.hdf5'
    self%fileRead=.false.
    return
  end function lossConeConstructorInternal

  subroutine lossConeDestructor(self)
    !!{
    Destructor for the \refClass{virialOrbitLossCone} virial orbits class.
    !!}
    implicit none
    type(virialOrbitLossCone), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%cosmologyParameters_"           />
    <objectDestructor name="self%cosmologicalVelocityField_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"      />
    <objectDestructor name="self%criticalOverdensity_"           />
    <objectDestructor name="self%darkMatterHaloBias_"            />
    <objectDestructor name="self%darkMatterHaloScale_"           />
    <objectDestructor name="self%linearGrowth_"                  />
    <objectDestructor name="self%virialDensityContrast_"         />
    <objectDestructor name="self%correlationFunctionTwoPoint_"   />
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    <objectDestructor name="self%darkMatterProfileDMO_"          />
    !!]
    return
  end subroutine lossConeDestructor

  function lossConeOrbit(self,node,host,acceptUnboundOrbits) result(orbit)
    !!{
    Return lossCone orbital parameters for a satellite.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    type            (keplerOrbit        )                        :: orbit
    class           (virialOrbitLossCone), intent(inout), target :: self
    type            (treeNode           ), intent(inout)         :: host                         , node
    logical                              , intent(in   )         :: acceptUnboundOrbits
    class           (nodeComponentBasic ), pointer               :: basicSatellite               , basicHost
    double precision                     , parameter             :: boundTolerance        =1.0d-4 ! Tolerence to ensure that orbits are sufficiently bound.
    integer         (c_size_t           ), dimension(0:1)        :: iSatellite                   , iHost                     , &
         &                                                          iRadial                      , iTangential
    double precision                     , dimension(0:1)        :: hSatellite                   , hHost                     , &
         &                                                          hRadial                      , hTangential
    double precision                                             :: velocityHost                 , distributionMaximum       , &
         &                                                          radiusHost                   , radiusHostSelf            , &
         &                                                          massSatellite                , massHost                  , &
         &                                                          velocityRadialInternal       , velocityTangentialInternal, &
         &                                                          distributionFunction         , energyInternal            , &
         &                                                          uniformRandom
    integer                                                      :: jSatellite                   , jHost                     , &
         &                                                          jRadial                      , jTangential
    logical                                                      :: foundOrbit

    call self%tabulate    (node,host                                                                                     )
    call self%interpolants(node,host,massSatellite,massHost,velocityHost,radiusHostSelf,iSatellite,iHost,hSatellite,hHost)
    ! Find the peak of our distribution for use in rejection sampling.
    distributionMaximum=0.0d0
    do jSatellite=0,1
       do jHost  =0,1
          distributionMaximum=+distributionMaximum                                                &
               &              +self%velocityDistributionPeak(iHost(jHost),iSatellite(jSatellite)) &
               &              *                              hHost(jHost)                         &
               &              *                                           hSatellite(jSatellite)
       end do
    end do
    ! Perform rejection sampling to find the orbital velocities.
    foundOrbit=.false.
    do while(.not.foundOrbit)
       ! Reset the orbit.
       call orbit%reset    (                            )
       ! Set basic properties of the orbit.
       call orbit%massesSet(massSatellite,massHost      )
       call orbit%radiusSet(              radiusHostSelf)
       ! Select potential radial and tangential velocities.
       velocityRadialInternal    =node%hostTree%randomNumberGenerator_%uniformSample()*self%velocityMaximum
       velocityTangentialInternal=node%hostTree%randomNumberGenerator_%uniformSample()*self%velocityMaximum
       ! Evaluate distribution function for these parameters.
       call self%interpolatorVelocity%linearFactors(velocityRadialInternal    ,iRadial    (0),hRadial    )
       call self%interpolatorVelocity%linearFactors(velocityTangentialInternal,iTangential(0),hTangential)
       iRadial    (1)=iRadial    (0)+1_c_size_t
       iTangential(1)=iTangential(0)+1_c_size_t
       distributionFunction=0.0d0
       do jSatellite          =0,1
          do jHost            =0,1
             do jRadial       =0,1
                do jTangential=0,1
                   distributionFunction=+distributionFunction                                       &
                        &                +self%velocityDistributionOrbits(                          &
                        &                                                 iHost      (jHost      ), &
                        &                                                 iSatellite (jSatellite ), &
                        &                                                 iRadial    (jRadial    ), &
                        &                                                 iTangential(jTangential)  &
                        &                                                )                          &
                        &                *                                hHost      (jHost      )  &
                        &                *                                hSatellite (jSatellite )  &
                        &                *                                hRadial    (jRadial    )  &
                        &                *                                hTangential(jTangential)
                end do
             end do
          end do
       end do
       ! Perform rejection sampling.
       uniformRandom=distributionMaximum*node%hostTree%randomNumberGenerator_%uniformSample()
       if (uniformRandom <= distributionFunction) then
          foundOrbit=.true.
          ! If requested, check that the orbit is bound. We require it to have E<-boundTolerance to ensure that it is sufficiently
          ! bound that later rounding errors will not make it appear unbound.
          if (.not.acceptUnboundOrbits) then
             energyInternal=-1.0d0+0.5d0*(velocityRadialInternal**2+velocityTangentialInternal**2)*orbit%specificReducedMass()
             foundOrbit=(energyInternal < -boundTolerance)
          end if
       end if
       if (.not.foundOrbit) cycle
       call orbit%velocityRadialSet    (velocityRadialInternal    *velocityHost)
       call orbit%velocityTangentialSet(velocityTangentialInternal*velocityHost)
       ! Propagate the orbit to the virial radius under the default density contrast definition.
       radiusHost=self%darkMatterHaloScale_%radiusVirial(host)
       if (orbit%radiusApocenter() >= radiusHost .and. orbit%radiusPericenter() <= radiusHost) then
          foundOrbit     =  .true.
          basicHost      => host%basic()
          basicSatellite => node%basic()
          call orbit%propagate(radiusHost           ,infalling=.true.)
          call orbit%massesSet(basicSatellite%mass(),basicHost%mass())
       end if
    end do
    return
  end function lossConeOrbit

  double precision function lossConeVelocityDistributionFunction(self,node,host,velocityRadial,velocityTangential)
    !!{
    Return the orbital velocity distribution function.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (virialOrbitLossCone), intent(inout)  :: self
    type            (treeNode           ), intent(inout)  :: host          , node
    double precision                     , intent(in   )  :: velocityRadial, velocityTangential    
    integer         (c_size_t           ), dimension(0:1) :: iSatellite    , iHost             , &
         &                                                   iRadial       , iTangential
    double precision                     , dimension(0:1) :: hSatellite    , hHost             , &
         &                                                   hRadial       , hTangential
    double precision                                      :: velocityHost  , radiusHost        , &
         &                                                   massSatellite , massHost
    integer                                               :: jSatellite    , jHost             , &
         &                                                   jRadial       , jTangential

    ! Tabulate the distribution and get all interpolating factors.
    call self%tabulate                          (node,host                                                                                 )
    call self%interpolants                      (node,host,massSatellite,massHost,velocityHost,radiusHost,iSatellite,iHost,hSatellite,hHost)
    call self%interpolatorVelocity%linearFactors(velocityRadial    /velocityHost,iRadial    (0),hRadial    )
    call self%interpolatorVelocity%linearFactors(velocityTangential/velocityHost,iTangential(0),hTangential)
    iRadial    (1)=iRadial    (0)+1_c_size_t
    iTangential(1)=iTangential(0)+1_c_size_t
    ! Perform the interpolation.
    lossConeVelocityDistributionFunction=0.0d0
    do jSatellite          =0,1
       do jHost            =0,1
          do jRadial       =0,1
             do jTangential=0,1
                lossConeVelocityDistributionFunction=+lossConeVelocityDistributionFunction                      &
                     &                               +self%velocityDistributionOrbits(                          &
                     &                                                                iHost      (jHost      ), &
                     &                                                                iSatellite (jSatellite ), &
                     &                                                                iRadial    (jRadial    ), &
                     &                                                                iTangential(jTangential)  &
                     &                                                               )                          &
                     &                               *                                hHost      (jHost      )  &
                     &                               *                                hSatellite (jSatellite )  &
                     &                               *                                hRadial    (jRadial    )  &
                     &                               *                                hTangential(jTangential)
             end do
          end do
       end do
    end do
    return
  end function lossConeVelocityDistributionFunction

  subroutine lossConeInterpolants(self,node,host,massSatellite,massHost,velocityHost,radiusHost,iSatellite,iHost,hSatellite,hHost)
    !!{
    Compute interpolating factors in the orbital parameter tables.
    !!}
    use, intrinsic :: ISO_C_Binding                       , only : c_size_t
    use            :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use            :: Galacticus_Nodes                    , only : nodeComponentBasic
    implicit none
    class           (virialOrbitLossCone), intent(inout)                 :: self
    type            (treeNode           ), intent(inout)                 :: host          , node
    double precision                     , intent(  out)                 :: velocityHost  , radiusHost, &
         &                                                                  massSatellite , massHost
    integer         (c_size_t           ), intent(  out), dimension(0:1) :: iSatellite    , iHost
    double precision                     , intent(  out), dimension(0:1) :: hSatellite    , hHost
    class           (nodeComponentBasic ), pointer                       :: basicSatellite, basicHost

    ! Evaluate halo masses under our mass definition. Also gives us the host velocity scale.
    basicSatellite => node%basic()
    basicHost      => host%basic()
    massSatellite  =  Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                       node                                                                                                , &
         &                                                                       self%virialDensityContrast_%densityContrast(basicSatellite%mass(),basicSatellite%timeLastIsolated()), &
         &                                                cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                                cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                                virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                                darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &                                               )
    massHost       =  Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                       host                                                                                                , &
         &                                                                       self%virialDensityContrast_%densityContrast(basicHost     %mass(),basicHost     %timeLastIsolated()), &
         &                                                radius                =radiusHost                                                                                          , &
         &                                                velocity              =velocityHost                                                                                        , &
         &                                                cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                                cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                                virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                                darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &                                               )
    massSatellite  =  min(massSatellite,massHost)
    ! Compute interpolating factors.
    call self%interpolatorMass%linearFactors(log(massSatellite),iSatellite(0),hSatellite)
    call self%interpolatorMass%linearFactors(log(massHost     ),iHost     (0),hHost     )
    iSatellite(1)=iSatellite(0)+1_c_size_t
    iHost     (1)=iHost     (0)+1_c_size_t
    return
  end subroutine lossConeInterpolants
  
  function lossConeDensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the calculation of orbital parameters.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: lossConeDensityContrastDefinition
    class(virialOrbitLossCone       ), intent(inout) :: self

    lossConeDensityContrastDefinition => self%virialDensityContrast_
    return
  end function lossConeDensityContrastDefinition

  double precision function lossConeVelocityTangentialMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    implicit none
    class           (virialOrbitLossCone), intent(inout)  :: self
    type            (treeNode           ), intent(inout)  :: node         , host
    integer         (c_size_t           ), dimension(0:1) :: iSatellite   , iHost
    double precision                     , dimension(0:1) :: hSatellite   , hHost
    double precision                                      :: velocityHost , radiusHost, &
         &                                                   massSatellite, massHost
    integer                                               :: jSatellite   , jHost

    call self%tabulate    (node,host                                                                                 )
    call self%interpolants(node,host,massSatellite,massHost,velocityHost,radiusHost,iSatellite,iHost,hSatellite,hHost)
    lossConeVelocityTangentialMagnitudeMean=0.0d0
    do jSatellite=0,1
       do jHost  =0,1
          lossConeVelocityTangentialMagnitudeMean=+lossConeVelocityTangentialMagnitudeMean                    &
               &                                  +self%velocityTangentialMeanVirial(                         &
               &                                                                     iHost     (jHost      ), &
               &                                                                     iSatellite(jSatellite )  &
               &                                                                    )                         &
               &                                  *                                  hHost     (jHost      )  &
               &                                  *                                  hSatellite(jSatellite )
       end do
    end do
    return
  end function lossConeVelocityTangentialMagnitudeMean

  function lossConeVelocityTangentialVectorMean(self,node,host)
    !!{
    Return the mean of the vector tangential velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                     , dimension(3)  :: lossConeVelocityTangentialVectorMean
    class           (virialOrbitLossCone), intent(inout) :: self
    type            (treeNode           ), intent(inout) :: node                                , host
    !$GLC attributes unused :: self, node, host

    lossConeVelocityTangentialVectorMean=0.0d0
    call Error_Report('vector velocity is not defined for this class'//{introspection:location})
    return
  end function lossConeVelocityTangentialVectorMean

  double precision function lossConeAngularMomentumMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the angular momentum.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    implicit none
    class           (virialOrbitLossCone), intent(inout) :: self
    type            (treeNode           ), intent(inout) :: node        , host
    class           (nodeComponentBasic ), pointer       :: basic       , basicHost
    double precision                                     :: massHost    , radiusHost, &
         &                                                  velocityHost

    basic                                =>  node%basic()
    basicHost                            =>  host%basic()
    massHost                             =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()),radiusHost,velocityHost,self%cosmologyParameters_,self%cosmologyFunctions_,self%virialDensityContrast_,self%darkMatterProfileDMO_)
    lossConeAngularMomentumMagnitudeMean =  +self%velocityTangentialMagnitudeMean(node,host) &
         &                                  *radiusHost                                      &
         &                                  /(                                               & ! Account for reduced mass.
         &                                    +1.0d0                                         &
         &                                    +basic    %mass()                              &
         &                                    /basicHost%mass()                              &
         &                                   )
    return
  end function lossConeAngularMomentumMagnitudeMean

  function lossConeAngularMomentumVectorMean(self,node,host)
    !!{
    Return the mean of the vector angular momentum.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                     , dimension(3)  :: lossConeAngularMomentumVectorMean
    class           (virialOrbitLossCone), intent(inout) :: self
    type            (treeNode           ), intent(inout) :: node                             , host
    !$GLC attributes unused :: self, node, host

    lossConeAngularMomentumVectorMean=0.0d0
    call Error_Report('vector angular momentum is not defined for this class'//{introspection:location})
    return
  end function lossConeAngularMomentumVectorMean

  double precision function lossConeVelocityTotalRootMeanSquared(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    implicit none
    class           (virialOrbitLossCone), intent(inout)  :: self
    type            (treeNode           ), intent(inout)  :: node         , host
    integer         (c_size_t           ), dimension(0:1) :: iSatellite   , iHost
    double precision                     , dimension(0:1) :: hSatellite   , hHost
    double precision                                      :: velocityHost , radiusHost, &
         &                                                   massSatellite, massHost
    integer                                               :: jSatellite   , jHost

    call self%tabulate    (node,host                                                                                 )
    call self%interpolants(node,host,massSatellite,massHost,velocityHost,radiusHost,iSatellite,iHost,hSatellite,hHost)
    lossConeVelocityTotalRootMeanSquared=0.0d0
    do jSatellite=0,1
       do jHost  =0,1
          lossConeVelocityTotalRootMeanSquared=+lossConeVelocityTotalRootMeanSquared           &
               &                               +self%velocityTotalRMS(                         &
               &                                                      iHost     (jHost      ), &
               &                                                      iSatellite(jSatellite )  &
               &                                                     )                         &
               &                               *                      hHost     (jHost      )  &
               &                               *                      hSatellite(jSatellite )
       end do
    end do
    return
  end function lossConeVelocityTotalRootMeanSquared

  double precision function lossConeEnergyMean(self,node,host)
    !!{
    Return the mean energy of the orbits.
    !!}
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Galacticus_Nodes                    , only : nodeComponentBasic                 , treeNode
    use :: Numerical_Constants_Astronomical    , only : gravitationalConstant_internal
    implicit none
    class           (virialOrbitLossCone), intent(inout) :: self
    type            (treeNode           ), intent(inout) :: node        , host
    class           (nodeComponentBasic ), pointer       :: basic       , basicHost
    double precision                                     :: massHost    , radiusHost, &
         &                                                  velocityHost

    basic              =>  node%basic()
    basicHost          =>  host%basic()
    massHost           =   Dark_Matter_Profile_Mass_Definition(host,self%virialDensityContrast_%densityContrast(basicHost%mass(),basicHost%timeLastIsolated()),radiusHost,velocityHost,self%cosmologyParameters_,self%cosmologyFunctions_,self%virialDensityContrast_,self%darkMatterProfileDMO_)
    lossConeEnergyMean =  +0.5d0                                           &
         &                *self%velocityTotalRootMeanSquared(node,host)**2 &
         &                /(                                               & ! Account for reduced mass.
         &                  +1.0d0                                         &
         &                  +basic    %mass()                              &
         &                  /basicHost%mass()                              &
         &                 )                                               &
         &                -gravitationalConstant_internal                  &
         &                *massHost                                        &
         &                /radiusHost
    return
  end function lossConeEnergyMean

  subroutine lossConeTabulate(self,nodeSatelliteTarget,nodeHostTarget)
    !!{
    Compute properties of infalling halos.
    !!}
    use :: Cosmological_Density_Field          , only : haloEnvironmentNormal
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Display                             , only : displayCounter                     , displayCounterClear , displayIndent, displayUnindent, &
          &                                             verbosityLevelWorking
    use :: Error                               , only : errorStatusSuccess
    use :: Galacticus_Nodes                    , only : mergerTree                         , nodeComponentBasic  , treeNode
    use :: Calculations_Resets                 , only : Calculations_Reset
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Ranges                    , only : Make_Range                         , rangeTypeLogarithmic
    use :: Table_Labels                        , only : extrapolationTypeFix
    use :: Numerical_Integration               , only : integrator
    use :: OMP_Lib                             , only : OMP_Get_Thread_Num
    implicit none
    class           (virialOrbitLossCone           ), intent(inout)     , target      :: self
    type            (treeNode                      ), intent(inout)                   :: nodeHostTarget                                 , nodeSatelliteTarget
    type            (treeNode                      )                    , pointer     :: nodeHost                                       , nodeSatellite
    class           (nodeComponentBasic            )                    , pointer     :: basicHost                                      , basicSatellite
    type            (mergerTree                    )                    , pointer     :: tree
    double precision                                , dimension(:      ), allocatable :: radius                                         , velocityDispersionLinearTable
    double precision                                , dimension(:,:    ), allocatable :: velocityRadialMeanVirial                       , velocityRadialDispersionVirial            , &
         &                                                                               velocityTangentialMeanVirial                   , velocityTangentialDispersionVirial        , &
         &                                                                               velocityDistributionPeak                       , velocityTotalRMS
    double precision                                , dimension(:,:,:  ), allocatable :: velocityRadialDistributionOrbits               , velocityTangentialDistributionOrbits
    double precision                                , dimension(:,:,:,:), allocatable :: velocityDistributionOrbits
    double precision                                , dimension(2      )              :: radiiEvaluation
    double precision                                , parameter                       :: velocityRadialInfallMaximum           =6.0000d0
    double precision                                , parameter                       :: velocityRadialInfallStep              =0.0001d0
    double precision                                , parameter                       :: radiusTableMinimum                    =1.0d-2
    double precision                                , parameter                       :: radiusTableMaximum                    =1.0d+3
    integer                                         , parameter                       :: radiusTablePointsPerDecade            =20
    integer                                         , parameter                       :: overdensityLimitLower                 =-10.0d0
    integer                                                                           :: iHost                                          , iSatellite                                , &
         &                                                                               countProgress                                  , countTotal                                , &
         &                                                                               countMasses                                    , countVelocities                           , &
         &                                                                               indexVelocityRadial                            , indexVelocityTangential                   , &
         &                                                                               indexVelocityRadialInfall                      , countVelocityRadialInfall                 , &
         &                                                                               iEvaluate                                      , countRadii
    double precision                                                                  :: massHost                                       , massSatellite                             , &
         &                                                                               radiusPericenterVirial                         , velocityDispersionLinear                  , &
         &                                                                               velocityRadialVirial                           , velocityTangentialVirial                  , &
         &                                                                               velocityRadialInfall                           , velocityTangentialInfall                  , &
         &                                                                               distributionFunction                           , jacobianDeterminant                       , &
         &                                                                               velocityDispersionRadialEvaluateVirial         , velocityDispersionTangentialEvaluateVirial, &
         &                                                                               velocityMeanRadialEvaluateVirial               , distributionFunctionCumulated             , &
         &                                                                               radiusApocenterVirial                          , timeOfFlightVirial                        , &
         &                                                                               timeOfFlight                                   , timeEvaluate                              , &
         &                                                                               velocityDispersionEvaluate                     , velocityRadialMeanEvaluate                , &
         &                                                                               radiusEvaluate                                 , radiusEvaluateComoving                    , &
         &                                                                               radiusEvaluateVirial                           , velocityVirialHost                        , &
         &                                                                               radiusVirialHost                               , radiusInfallTerm1                         , &
         &                                                                               radiusInfallTerm2                              , factorEnvironmental                       , &
         &                                                                               massEnvironment                                , radiusEnvironment                         , &
         &                                                                               jacobianFactor                                 , jacobianSign                              , &
         &                                                                               radiusEvaluateLagrangian
    type            (interpolator                  )                                  :: interpolatorVelocityDispersionLinear
    type            (integrator                    )                                  :: integratorEnvironment                          , integratorEnvironmentNormalizer
    
    ! Read in any existing tabulation from file.
    call self%restoreTable()
    ! Check if we need to retabulate.
    basicSatellite => nodeSatelliteTarget%basic()
    basicHost      => nodeHostTarget     %basic()
    massSatellite  =  Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                       nodeSatelliteTarget                                                                                 , &
         &                                                                       self%virialDensityContrast_%densityContrast(basicSatellite%mass(),basicSatellite%timeLastIsolated()), &
         &                                                cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                                cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                                virialDensityContrast_=self%virialDensityContrast_                                                                         , &
         &                                                darkMatterProfileDMO_ =self%darkMatterProfileDMO_                                                                            &
         &  )
    massHost       =  Dark_Matter_Profile_Mass_Definition(                                                                                                                             &
         &                                                                       nodeHostTarget                                                                                      , &
         &                                                                       self%virialDensityContrast_%densityContrast(basicHost     %mass(),basicHost     %timeLastIsolated()), &
         &                                                cosmologyParameters_  =self%cosmologyParameters_                                                                           , &
         &                                                cosmologyFunctions_   =self%cosmologyFunctions_                                                                            , &
         &                                                virialDensityContrast_=self%virialDensityContrast_                                                                           &
         & )
    if     (                                           &
         &   basicHost     %time() == self%time        &
         &  .and.                                      &
         &   basicHost     %mass() >= self%massMinimum &
         &  .and.                                      &
         &   basicHost     %mass() <= self%massMaximum &
         &  .and.                                      &
         &   basicSatellite%mass() >= self%massMinimum &
         &  .and.                                      &
         &   basicSatellite%mass() <= self%massMaximum &
         & ) return
    ! Free existing tabulations.
    if (allocated(self%mass)) then
       deallocate(self%mass                                )
       deallocate(self%velocity                            )
       deallocate(self%velocityRadialMeanVirial            )
       deallocate(self%velocityRadialDispersionVirial      )
       deallocate(self%velocityTangentialMeanVirial        )
       deallocate(self%velocityTangentialDispersionVirial  )
       deallocate(self%velocityRadialDistributionOrbits    )
       deallocate(self%velocityTangentialDistributionOrbits)
       deallocate(self%velocityDistributionOrbits          )
       deallocate(self%velocityTotalRMS                    )
       deallocate(self%velocityDistributionPeak            )
       deallocate(self%interpolatorMass                    )
       deallocate(self%interpolatorVelocity                )
    end if
    ! Set time for this tabulation.
    self%time=basicHost%time()
    ! Get number of velocities to tabulate.
    countVelocities=size(self%velocity)
    ! Build range of masses.
    self%massMinimum=min(self%massMinimum,0.5d0*min(basicHost%mass(),basicSatellite%mass()))
    self%massMaximum=max(self%massMaximum,2.0d0*max(basicHost%mass(),basicSatellite%mass()))
    countMasses     =int(log10(self%massMaximum/self%massMinimum)*dble(self%countMassesPerDecade))+1    
    allocate(self%mass            (countMasses))
    allocate(self%interpolatorMass             )
    self%mass            =Make_Range(self%massMinimum,self%massMaximum,countMasses,rangeTypeLogarithmic)
    self%interpolatorMass=interpolator(log(self%mass))
    ! Allocate arrays for results and initialize.
    allocate(velocityRadialMeanVirial            (countMasses,countMasses                                ))
    allocate(velocityRadialDispersionVirial      (countMasses,countMasses                                ))
    allocate(velocityTangentialMeanVirial        (countMasses,countMasses                                ))
    allocate(velocityTangentialDispersionVirial  (countMasses,countMasses                                ))
    allocate(velocityRadialDistributionOrbits    (countMasses,countMasses,countVelocities                ))
    allocate(velocityTangentialDistributionOrbits(countMasses,countMasses,countVelocities                ))
    allocate(velocityDistributionOrbits          (countMasses,countMasses,countVelocities,countVelocities))
    allocate(velocityTotalRMS                    (countMasses,countMasses                                ))
    allocate(velocityDistributionPeak            (countMasses,countMasses                                ))
    velocityRadialMeanVirial                     =0.0d0
    velocityRadialDispersionVirial               =0.0d0
    velocityTangentialMeanVirial                 =0.0d0
    velocityTangentialDispersionVirial           =0.0d0
    velocityRadialDistributionOrbits             =0.0d0
    velocityTangentialDistributionOrbits         =0.0d0
    velocityDistributionOrbits                   =0.0d0
    velocityTotalRMS                             =0.0d0
    velocityDistributionPeak                     =0.0d0
    ! Iterate over host masses.
    countTotal   =(countMasses*(countMasses+1))/2
    countProgress=0
    !$omp parallel private(iHost,iSatellite,massHost,massSatellite,tree,nodeHost,nodeSatellite,basicHost,basicSatellite,radiusVirialHost,velocityVirialHost,velocityDispersionLinear,indexVelocityRadial,indexVelocityTangential,velocityRadialVirial,velocityTangentialVirial,countVelocityRadialInfall,indexVelocityRadialInfall,velocityRadialInfall,radiusInfallTerm1,radiusInfallTerm2,radiiEvaluation,iEvaluate,radiusEvaluateVirial,timeEvaluate,radiusApocenterVirial,radiusPericenterVirial,timeOfFlightVirial,timeOfFlight,radiusEvaluate,radiusEvaluateComoving,velocityDispersionEvaluate,velocityRadialMeanEvaluate,velocityDispersionRadialEvaluateVirial,velocityDispersionTangentialEvaluateVirial,velocityMeanRadialEvaluateVirial,velocityTangentialInfall,jacobianFactor,jacobianDeterminant,distributionFunction,interpolatorVelocityDispersionLinear,integratorEnvironment,integratorEnvironmentNormalizer,factorEnvironmental,massEnvironment,radiusEnvironment)
    allocate(tree                                                                          )
    allocate(     cosmologyFunctions_            ,mold=self%cosmologyFunctions_            )
    allocate(     cosmologyParameters_           ,mold=self%cosmologyParameters_           )
    allocate(     darkMatterHaloBias_            ,mold=self%darkMatterHaloBias_            )
    allocate(     darkMatterHaloScale_           ,mold=self%darkMatterHaloScale_           )
    allocate(     linearGrowth_                  ,mold=self%linearGrowth_                  )
    allocate(     correlationFunctionTwoPoint_   ,mold=self%correlationFunctionTwoPoint_   )
    allocate(     cosmologicalVelocityField_     ,mold=self%cosmologicalVelocityField_     )
    allocate(     criticalOverdensity_           ,mold=self%criticalOverdensity_           )
    allocate(     cosmologicalMassVariance_      ,mold=self%cosmologicalMassVariance_      )
    allocate(     mergerTreeBranchingProbability_,mold=self%mergerTreeBranchingProbability_)
    call tree%properties%initialize()
    !$omp critical(virialOrbitLossConeDeepCopy)
    !![
    <deepCopyReset variables="self%cosmologyFunctions_ self%cosmologyParameters_ self%darkMatterHaloBias_ self%darkMatterHaloScale_ self%linearGrowth_ self%correlationFunctionTwoPoint_ self%cosmologicalVelocityField_ self%cosmologicalMassVariance_ self%criticalOverdensity_ self%mergerTreeBranchingProbability_"/>
    <deepCopy source="self%cosmologyFunctions_"             destination="cosmologyFunctions_"            />
    <deepCopy source="self%cosmologyParameters_"            destination="cosmologyParameters_"           />
    <deepCopy source="self%darkMatterHaloBias_"             destination="darkMatterHaloBias_"            />
    <deepCopy source="self%darkMatterHaloScale_"            destination="darkMatterHaloScale_"           />
    <deepCopy source="self%linearGrowth_"                   destination="linearGrowth_"                  />
    <deepCopy source="self%correlationFunctionTwoPoint_"    destination="correlationFunctionTwoPoint_"   />
    <deepCopy source="self%cosmologicalVelocityField_"      destination="cosmologicalVelocityField_"     />
    <deepCopy source="self%cosmologicalMassVariance_"       destination="cosmologicalMassVariance_"      />
    <deepCopy source="self%mergerTreeBranchingProbability_" destination="mergerTreeBranchingProbability_"/>
    <deepCopy source="self%criticalOverdensity_"            destination="criticalOverdensity_"           />
    <deepCopyFinalize variables="cosmologyFunctions_ cosmologyParameters_ darkMatterHaloBias_ darkMatterHaloScale_ linearGrowth_ correlationFunctionTwoPoint_ cosmologicalVelocityField_ cosmologicalMassVariance_ criticalOverdensity_ mergerTreeBranchingProbability_"/>
    !!]
    !$omp end critical(virialOrbitLossConeDeepCopy)
    ! Construct a halo environment.
    radiusEnvironment=+20.0d0
    massEnvironment  =+ 4.0d0                                              &
         &            / 3.0d0                                              &
         &            *Pi                                                  &
         &            *cosmologyFunctions_%matterDensityEpochal(self%time) &
         &            *radiusEnvironment**3
    allocate(haloEnvironment_                      )
    allocate(cosmologicalMassVarianceEnvironmental_)
    allocate(criticalOverdensityEnvironmental_     )
    allocate(haloMassFunctionEnvironmental_        )
    !![
    <referenceConstruct object="haloEnvironment_"                      >
     <constructor>
      haloEnvironmentNormal                      (                                                                       &amp;
       &amp;                                      radiusEnvironment        =radiusEnvironment                          , &amp;
       &amp;                                      time                     =cosmologyFunctions_      %cosmicTime(1.0d0), &amp;
       &amp;                                      cosmologyParameters_     =cosmologyParameters_                       , &amp;
       &amp;                                      cosmologyFunctions_      =cosmologyFunctions_                        , &amp;
       &amp;                                      cosmologicalMassVariance_=cosmologicalMassVariance_                  , &amp;
       &amp;                                      linearGrowth_            =linearGrowth_                              , &amp;
       &amp;                                      criticalOverdensity_     =criticalOverdensity_                         &amp;
       &amp;                                      )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="cosmologicalMassVarianceEnvironmental_">
     <constructor>
      cosmologicalMassVariancePeakBackgroundSplit(                                                                       &amp;
       &amp;                                      haloEnvironment_         =haloEnvironment_                           , &amp;
       &amp;                                      cosmologicalMassVariance_=cosmologicalMassVariance_                  , &amp;
       &amp;                                      cosmologyParameters_     =cosmologyParameters_                       , &amp;
       &amp;                                      cosmologyFunctions_      =cosmologyFunctions_                          &amp;
       &amp;                                     )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="criticalOverdensityEnvironmental_"     >
     <constructor>
      criticalOverdensityPeakBackgroundSplit     (                                                                       &amp;
       &amp;                                      criticalOverdensity_     =criticalOverdensity_                       , &amp;
       &amp;                                      haloEnvironment_         =haloEnvironment_                           , &amp;
       &amp;                                      cosmologyFunctions_      =cosmologyFunctions_                        , &amp;
       &amp;                                      cosmologicalMassVariance_=cosmologicalMassVariance_                  , &amp;
       &amp;                                      linearGrowth_            =linearGrowth_                                &amp;
       &amp;                                     )
     </constructor>
    </referenceConstruct>
    <referenceConstruct object="haloMassFunctionEnvironmental_"        >
     <constructor>
      haloMassFunctionShethTormen                (                                                                       &amp;
       &amp;                                      cosmologyParameters_     =cosmologyParameters_                       , &amp;
       &amp;                                      cosmologicalMassVariance_=cosmologicalMassVarianceEnvironmental_     , &amp;
       &amp;                                      criticalOverdensity_     =criticalOverdensityEnvironmental_          , &amp;
       &amp;                                      a                        =self%haloMassFunctionA                     , &amp;
       &amp;                                      p                        =self%haloMassFunctionP                     , &amp;
       &amp;                                      normalization            =self%haloMassFunctionNormalization           &amp;
       &amp;                                     )
     </constructor>
    </referenceConstruct>
    !!]
    integratorEnvironment          =integrator(integrand=integrandEnvironment             ,toleranceRelative=1.0d-3)
    integratorEnvironmentNormalizer=integrator(integrand=integrandEnvironmentNormalization,toleranceRelative=1.0d-3)
    do iHost=1,countMasses
       ! Build host node.
       massHost           =  self%mass     (            iHost)
       massHost_          =  massHost
       nodeHost           => treeNode      (                 )
       nodeHost_          => nodeHost
       basicHost          => nodeHost%basic(autoCreate=.true.)
       nodeHost %hostTree => tree
       tree     %nodeBase => nodeHost
       call basicHost%massSet            (     massHost)
       call basicHost%timeSet            (self%time    )
       call basicHost%timeLastIsolatedSet(self%time    )
       call Calculations_Reset           (     nodeHost)
       ! Compute host virial properties.
       radiusVirialHost  =darkMatterHaloScale_%radiusVirial  (nodeHost)
       velocityVirialHost=darkMatterHaloScale_%velocityVirial(nodeHost)
       ! Compute the environmental boost factor for velocity dispersion.
       timeEvaluate_      =+                                                                                                                                 self%time
       factorEnvironmental=+integratorEnvironment          %integrate(overdensityLimitLower*cosmologicalMassVariance_%rootVariance(mass=massEnvironment,time=self%time),haloEnvironment_%overdensityLinearMaximum()) &
            &              /integratorEnvironmentNormalizer%integrate(overdensityLimitLower*cosmologicalMassVariance_%rootVariance(mass=massEnvironment,time=self%time),haloEnvironment_%overdensityLinearMaximum())
       ! Iterate over satellite masses.
       do iSatellite=1,countMasses          
          ! Only consider satellites less (or equally) massive than their host.
          if (iSatellite > iHost) cycle
          if (OMP_Get_Thread_Num() == 0) then
             call displayCounter(int(100.0d0*dble(countProgress)/dble(countTotal)),isNew=countProgress==0,verbosity=verbosityLevelWorking)
             countProgress=countProgress+1             
          end if
          ! Build satellite node.
          massSatellite           =  self%mass          (           iSatellite)
          nodeSatellite           => treeNode           (                     )
          basicSatellite          => nodeSatellite%basic(autoCreate=.true.    )
          nodeSatellite %hostTree => tree
          call basicSatellite%massSet            (     massSatellite)
          call basicSatellite%timeSet            (self%time         )
          call basicSatellite%timeLastIsolatedSet(self%time         )
          call Calculations_Reset                (     nodeSatellite)
          ! Tabulate the 1-D linear theory relative velocity of host-satellite halo pairs as a function of radius.
          !$omp single
          if (allocated(radius)) then
             deallocate(radius                       )
             deallocate(velocityDispersionLinearTable)
          end if
          countRadii=int(log10(radiusTableMaximum/radiusTableMinimum)*dble(radiusTablePointsPerDecade)+1)
          allocate(radius                       (countRadii))
          allocate(velocityDispersionLinearTable(countRadii))
          radius=Make_Range(radiusTableMinimum,radiusTableMaximum,countRadii,rangeTypeLogarithmic)
          !$omp end single
          !$omp do schedule(dynamic)
          do iEvaluate=1,countRadii
             velocityDispersionLinearTable(iEvaluate)=+                           factorEnvironmental                                                                  &
                  &                                   *self                      %velocityDispersionMultiplier                                                         &
                  &                                   *cosmologicalVelocityField_%velocityDispersion1DHaloPairwise(massHost,massSatellite,radius(iEvaluate),self%time)             
          end do
          !$omp end do
          interpolatorVelocityDispersionLinear=interpolator(radius,velocityDispersionLinearTable,extrapolationType=extrapolationTypeFix)          
          ! Initialize accumulators.
          !$omp single
          distributionFunctionCumulated=0.0d0
          !$omp end single
          !! Iterate over velocities.
          !$omp do reduction(+:velocityRadialDistributionOrbits,velocityTangentialDistributionOrbits,velocityDistributionOrbits,velocityRadialMeanVirial,velocityRadialDispersionVirial,velocityTangentialMeanVirial,velocityTangentialDispersionVirial,distributionFunctionCumulated) schedule(dynamic)
          do indexVelocityRadial       =1,countVelocities
             do indexVelocityTangential=1,countVelocities
                ! Extract velocities at the virial radius (in virial units).
                velocityRadialVirial    =self%velocity(indexVelocityRadial    )
                velocityTangentialVirial=self%velocity(indexVelocityTangential)
                ! Integrate over radial velocities at infall point.
                countVelocityRadialInfall=int(velocityRadialInfallMaximum/velocityRadialInfallStep)
                do indexVelocityRadialInfall=0,countVelocityRadialInfall
                   velocityRadialInfall=-dble(indexVelocityRadialInfall) &
                        &               *     velocityRadialInfallStep                   
                   ! Evaluate the radii at which this radial velocity is achieved. There will be (in general) two such
                   ! radii. Solve the quadratic equation for these radii.                   
                   !! Check that roots of the quadratic are real and finite.
                   radiusInfallTerm1=+  4.0d0                       &
                        &            +  4.0d0                       &
                        &            *  velocityTangentialVirial**2 &
                        &            *(                             &
                        &              -2.0d0                       &
                        &              +velocityRadialVirial    **2 &
                        &              +velocityTangentialVirial**2 &
                        &              -velocityRadialInfall    **2 &
                        &             )
                   radiusInfallTerm2=-  2.0d0                       &
                        &            +  velocityRadialVirial    **2 &
                        &            +  velocityTangentialVirial**2 &
                        &            -  velocityRadialInfall    **2
                    if     (                            &
                        &   radiusInfallTerm1 <  0.0d0 &
                        &  .or.                        &
                        &   radiusInfallTerm2 == 0.0d0 &
                        & ) cycle
                   !! Evaluate both roots of the equation to give the radii at which the current radial velocity is achieved.
                   radiiEvaluation(1)=(-2.0d0-sqrt(radiusInfallTerm1))/2.0d0/radiusInfallTerm2
                   radiiEvaluation(2)=(-2.0d0+sqrt(radiusInfallTerm1))/2.0d0/radiusInfallTerm2
                   !! Evaluate both possible radii.                    
                   do iEvaluate=1,2
                      radiusEvaluateVirial=radiiEvaluation(iEvaluate)
                      ! Evaluate only if the evaluation radius lies outside of the virial radius, and within a plausible range of influence.
                      if (radiusEvaluateVirial <= 1.0d0) cycle
                      ! Find the flight time from the evaluation radius to the virial radius, and evaluate the linear velocity
                      ! field at that lookback time.                         
                      !! Compute the apocenter of the orbit. For unbound orbits the apocentric radius will be negative - this
                      !! is acceptable and is handled correctly by the function which evaluates the time of flight.
                      radiusApocenterVirial =+(                                                                                                            &
                           &                   -2.0d0                                                                                                      &
                           &                   -sqrt(4.0d0+4.0d0*velocityTangentialVirial**2*(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2)) &
                           &                  )                                                                                                            &
                           &                 /  2.0d0                                                                                                      &
                           &                 /                                               (-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2)
                      radiusPericenterVirial=+(                                                                                                            &
                           &                   -2.0d0                                                                                                      &
                           &                   +sqrt(4.0d0+4.0d0*velocityTangentialVirial**2*(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2)) &
                           &                  )                                                                                                            &
                           &                 /  2.0d0                                                                                                      &
                           &                 /                                               (-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2)
                      timeOfFlightVirial   =+abs(                                                                                                            &
                           &                     +timeAlongOrbit(radiusEvaluateVirial,radiusApocenterVirial,radiusPericenterVirial,velocityTangentialVirial) &
                           &                     -timeAlongOrbit(1.0d0               ,radiusApocenterVirial,radiusPericenterVirial,velocityTangentialVirial) &
                           &                    )
                      timeOfFlight=+                     timeOfFlightVirial           &
                           &       *darkMatterHaloScale_%timescaleDynamical(nodeHost)
                      timeEvaluate=+self%time                                         &
                           &       -     timeOfFlight
                      ! Skip orbits where the evaluation time is prior to the Big Bang.
                      if (timeEvaluate <= 0.0d0) cycle
                      ! If growth of the linear velocity field during the flight of the secondary is to be included. Reset the time at which to evaluate to the present time.
                      if (self%includeInFlightGrowth) timeEvaluate=self%time
                      ! Compute the velocity field at the evaluation radius.
                      !! Get the correlation in velocities at the Lagrangian radius corresponding to the evaluation distance.
                      radiusEvaluate            =+radiusEvaluateVirial                              &
                           &                     *radiusVirialHost
                      radiusEvaluateComoving    =+radiusEvaluate                                    &
                           &                     /cosmologyFunctions_%expansionFactor(timeEvaluate)                      
                      radiusEvaluateLagrangian  =(                                                                                              &
                           &                      + darkMatterHaloScale_        %densityMean              (                      nodeHost     ) &
                           &                      / cosmologyFunctions_         %matterDensityEpochal     (                      self    %time) &
                           &                      *radiusVirialHost**3                                                                          &
                           &                      +(                                                                                            &
                           &                       +1.0d0                                                                                       &
                           &                       +correlationFunctionTwoPoint_%correlationVolumeAveraged(radiusEvaluateComoving,self   %time) &
                           &                      )                                                                                             &
                           &                      *radiusEvaluate  **3                                                                          &
                           &                     )**(1.0d0/3.0d0)
                      velocityDispersionEvaluate=+interpolatorVelocityDispersionLinear%interpolate(radiusEvaluateLagrangian)    &
                           &                     *(                                                                             &
                           &                       +cosmologyFunctions_%hubbleParameterEpochal              (     timeEvaluate) &
                           &                       *cosmologyFunctions_%expansionFactor                     (     timeEvaluate) &
                           &                       *linearGrowth_      %logarithmicDerivativeExpansionFactor(     timeEvaluate) &
                           &                       /cosmologyFunctions_%hubbleParameterEpochal              (self%time        ) &
                           &                       /cosmologyFunctions_%expansionFactor                     (self%time        ) &
                           &                       /linearGrowth_      %logarithmicDerivativeExpansionFactor(self%time        ) &
                           &                      )
                      ! Compute radial mean velocity at this separation and time. Note that we include the Hubble expansion term
                      ! here.
                      velocityRadialMeanEvaluate=+cosmologicalVelocityField_%velocityRadialMeanPairwise(radiusEvaluate,timeEvaluate,includeHubbleFlow=.true.)
                      ! Compute parameters of the velocity distribution function in virial units.
                      velocityDispersionRadialEvaluateVirial    =+velocityDispersionEvaluate/velocityVirialHost
                      velocityDispersionTangentialEvaluateVirial=+velocityDispersionEvaluate/velocityVirialHost
                      velocityMeanRadialEvaluateVirial          =+velocityRadialMeanEvaluate/velocityVirialHost
                      ! Compute the tangential velocity at the evaluation radius. Simply apply conservation of angular momentum.
                      velocityTangentialInfall=+velocityTangentialVirial &
                           &                   /radiusEvaluateVirial
                      ! Evaluate components of the Jacobian of the transformation from coordinates (r,vt') to (vr,vt).
                      !! First evaluate the factor which appears in the denominator of the Jacobian determinant. This
                      !! passes through zero at a certain radial infall velocity, leading to a caustic-like structure. We
                      !! exclude a small region around this caustic to avoid numerical divergences.
                      jacobianFactor=-2.0d0+velocityRadialVirial**2-velocityRadialInfall**2+velocityTangentialVirial**2
                      if (abs(jacobianFactor) < 1.0d-2) then
                         ! Close to the caustic - set the determinant of the Jacobian to zero.
                         if (iEvaluate == 1) then
                            ! "Negative" radius solution.
                            jacobianDeterminant=+velocityRadialVirial                                                                                            &
                                 &              *( velocityTangentialVirial**2/2.0d0-1.0d0/(2.0d0*(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2))) &
                                 &              +  velocityRadialVirial                                                                                          &
                                 &              /                       sqrt(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2)                         &
                                 &              /(+velocityRadialInfall-sqrt(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2))                        &
                                 &              +  velocityRadialVirial                                                                                          &
                                 &              *(                                                                                                               &
                                 &                + 1.0d0                                                                                                        &
                                 &                +12.0d0                        *velocityTangentialVirial**4                                                    &
                                 &                -12.0d0*velocityRadialVirial**2*velocityTangentialVirial**4                                                    &
                                 &                + 3.0d0*velocityRadialVirial**4*velocityTangentialVirial**4                                                    &
                                 &                -12.0d0                        *velocityTangentialVirial**6                                                    &
                                 &                + 6.0d0*velocityRadialVirial**2*velocityTangentialVirial**6                                                    &
                                 &                + 3.0d0                        *velocityTangentialVirial**8                                                    &
                                 &               )                                                                                                               &
                                 &              *(velocityRadialInfall-sqrt(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2))                         &
                                 &              /(4.0d0*(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2)**1.5d0)
                         else
                            ! "Positive" radius solution.
                            if (jacobianFactor < 0.0d0) then
                               ! Use the limiting value at this point.
                               jacobianDeterminant=-0.50d0                           *velocityRadialVirial   *velocityTangentialVirial**2
                            else
                               ! Use a series expansion.
                               jacobianDeterminant=-0.50d0                           *velocityRadialVirial   *velocityTangentialVirial**2                           &
                                    &              -0.75d0                           *velocityRadialVirial   *velocityTangentialVirial**4                           &
                                    &              *                      sqrt(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2)                          &
                                    &              *(velocityRadialInfall-sqrt(-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2))
                            end if
                         end if
                      else
                         ! Far from the caustic - evaluate the determinant of the Jacobian.
                         if      (iEvaluate == 1) then
                            jacobianSign=+1.0d0
                         else if (iEvaluate == 2) then
                            jacobianSign=-1.0d0
                         end if
                         jacobianDeterminant=+jacobianSign                                                                                                                                             &
                              &              *velocityRadialVirial                                                                                                                                     &
                              &              *(+1.0d0*jacobianSign+1.0d0/sqrt(1.0d0+(-2.0d0+velocityRadialVirial**2-velocityRadialInfall**2)*velocityTangentialVirial**2+velocityTangentialVirial**4)) &
                              &              /(+2.0d0-velocityRadialVirial**2+velocityRadialInfall**2-velocityTangentialVirial**2)
                      end if
                      ! Compute the distribution function at the virial radius, and multiply by the Jacobian of the transformation
                      ! to our velocity coordinates at the virial radius.
                      distributionFunction=+exp(                                                                                          & ! Gaussian distribution (offset to the mean radial velocity) for the radial component.
                           &                    -0.5d0                                                                                    &
                           &                    *(                                                                                        &
                           &                      +(                                                                                      &
                           &                        +velocityRadialInfall                                                                 &
                           &                        -velocityMeanRadialEvaluateVirial                                                     &
                           &                       )                                                                                      &
                           &                      /  velocityDispersionRadialEvaluateVirial                                               &
                           &                     )**2                                                                                     &
                           &                   )                                                                                          &
                           &               /sqrt(2.0d0*Pi)                                                                                &
                           &               /         velocityDispersionRadialEvaluateVirial                                               &
                           &               *         velocityTangentialInfall                                                             & ! Rayleigh distribution for the tangential component.
                           &               *exp(                                                                                          &
                           &                    -0.5d0                                                                                    &
                           &                    *(                                                                                        &
                           &                      +  velocityTangentialInfall                                                             &
                           &                      /  velocityDispersionTangentialEvaluateVirial                                           &
                           &                     )**2                                                                                     &
                           &                   )                                                                                          &
                           &               /         velocityDispersionTangentialEvaluateVirial**2                                        &
                           &               *4.0d0                                                                                         & ! Radial volume element.
                           &               *Pi                                                                                            &
                           &               *radiusEvaluateVirial**2                                                                       &
                           &               *(                                                                                             & ! Radial correlation term.
                           &                 +1.0d0                                                                                       &
                           &                 +darkMatterHaloBias_         %bias       (nodeHost     ,radiusEvaluate                     ) &
                           &                 *darkMatterHaloBias_         %bias       (nodeSatellite,radiusEvaluate                     ) &
                           &                 *correlationFunctionTwoPoint_%correlation(              radiusEvaluateComoving,timeEvaluate) &
                           &               )                                                                                              &
                           &              *abs(jacobianDeterminant) ! Coordinate transformation.
                      ! Weight the distribution function by the radial velocity to account for the fact that higher radial velocity halos merge more frequently.
                      distributionFunction=+distributionFunction &
                           &               *velocityRadialVirial
                      ! Accumulate the distribution, marginal distributions, and summary statistics.
                      velocityRadialDistributionOrbits    (iHost,iSatellite,indexVelocityRadial                        )=+velocityRadialDistributionOrbits    (iHost,iSatellite,indexVelocityRadial                        ) &
                           &                                                                                             +distributionFunction
                      velocityTangentialDistributionOrbits(iHost,iSatellite,                    indexVelocityTangential)=+velocityTangentialDistributionOrbits(iHost,iSatellite,                    indexVelocityTangential) &
                           &                                                                                             +distributionFunction
                      velocityDistributionOrbits          (iHost,iSatellite,indexVelocityRadial,indexVelocityTangential)=+velocityDistributionOrbits          (iHost,iSatellite,indexVelocityRadial,indexVelocityTangential) &
                           &                                                                                             +distributionFunction
                      velocityRadialMeanVirial            (iHost,iSatellite                                            )=+velocityRadialMeanVirial            (iHost,iSatellite                                            ) &
                           &                                                                                             +distributionFunction                                                                               &
                           &                                                                                             *velocityRadialVirial
                      velocityRadialDispersionVirial      (iHost,iSatellite                                            )=+velocityRadialDispersionVirial      (iHost,iSatellite                                            ) &
                           &                                                                                             +distributionFunction                                                                               &
                           &                                                                                             *velocityRadialVirial    **2
                      velocityTangentialMeanVirial        (iHost,iSatellite                                            )=+velocityTangentialMeanVirial        (iHost,iSatellite                                            ) &
                           &                                                                                             +distributionFunction                                                                               &
                           &                                                                                             *velocityTangentialVirial
                      velocityTangentialDispersionVirial  (iHost,iSatellite                                            )=+velocityTangentialDispersionVirial  (iHost,iSatellite                                            ) &
                           &                                                                                             +distributionFunction                                                                               &
                           &                                                                                             *velocityTangentialVirial**2
                      distributionFunctionCumulated                                                                     =+distributionFunctionCumulated                                                                      &
                           &                                                                                             +distributionFunction
                    end do
                end do
             end do
          end do
          !$omp end do
          !$omp single
          ! Compute the summary statistics.
          if (distributionFunctionCumulated > 0.0d0) then
             velocityRadialMeanVirial          (iHost,iSatellite)=     velocityRadialMeanVirial          (iHost,iSatellite)/distributionFunctionCumulated
             velocityTangentialMeanVirial      (iHost,iSatellite)=     velocityTangentialMeanVirial      (iHost,iSatellite)/distributionFunctionCumulated
             velocityRadialDispersionVirial    (iHost,iSatellite)=sqrt(velocityRadialDispersionVirial    (iHost,iSatellite)/distributionFunctionCumulated-velocityRadialMeanVirial    (iHost,iSatellite)**2)
             velocityTangentialDispersionVirial(iHost,iSatellite)=sqrt(velocityTangentialDispersionVirial(iHost,iSatellite)/distributionFunctionCumulated-velocityTangentialMeanVirial(iHost,iSatellite)**2)
          end if
          ! Normalize the distribution functions.
          if     (sum(velocityRadialDistributionOrbits    (iHost,iSatellite,:  )) > 0.0d0) &
               &      velocityRadialDistributionOrbits    (iHost,iSatellite,:  )=          &
               & +    velocityRadialDistributionOrbits    (iHost,iSatellite,:  )           &
               & /sum(velocityRadialDistributionOrbits    (iHost,iSatellite,:  ))          &
               & /(self%velocity(2)-self%velocity(1))
          if     (sum(velocityTangentialDistributionOrbits(iHost,iSatellite,  :)) > 0.0d0) &
               &      velocityTangentialDistributionOrbits(iHost,iSatellite,  :)=          &
               & +    velocityTangentialDistributionOrbits(iHost,iSatellite,  :)           &
               & /sum(velocityTangentialDistributionOrbits(iHost,iSatellite,  :))          &
               & /(self%velocity(2)-self%velocity(1))
          if (sum(    velocityDistributionOrbits          (iHost,iSatellite,:,:)) > 0.0d0) &
               &      velocityDistributionOrbits          (iHost,iSatellite,:,:)=          &
               & +    velocityDistributionOrbits          (iHost,iSatellite,:,:)           &
               & /sum(velocityDistributionOrbits          (iHost,iSatellite,:,:))          &
               & /(self%velocity(2)-self%velocity(1))**2
          ! Find the peak value in the distribution function.
          velocityDistributionPeak(iHost,iSatellite)=maxval(velocityDistributionOrbits(iHost,iSatellite,:,:))
          ! Find the rms total velocity.
          velocityTotalRMS(iHost,iSatellite)=0.0d0
          do indexVelocityRadial=1,countVelocities
             do indexVelocityTangential=1,countVelocities
                velocityTotalRMS(iHost,iSatellite)=+       velocityTotalRMS          (iHost,iSatellite                                            )    &
                     &                             +       velocityDistributionOrbits(iHost,iSatellite,indexVelocityRadial,indexVelocityTangential)    &
                     &                             *(                                                                                                  &
                     &                               +self%velocity                  (                 indexVelocityRadial                        )**2 &
                     &                               +self%velocity                  (                                     indexVelocityTangential)**2 &
                     &                              )
             end do
             velocityTotalRMS(iHost,iSatellite)=sqrt(                                                       &
                  &                                  +    velocityTotalRMS          (iHost,iSatellite    )  &
                  &                                  /sum(velocityDistributionOrbits(iHost,iSatellite,:,:)) &
                  &                                 )
          end do
          !$omp end single
          ! Clean up.
          deallocate(nodeSatellite)
       end do
       deallocate(nodeHost)
    end do
    ! Clean up.
    !![
    <objectDestructor name="cosmologyFunctions_"                   />
    <objectDestructor name="cosmologyParameters_"                  />
    <objectDestructor name="darkMatterHaloBias_"                   />
    <objectDestructor name="darkMatterHaloScale_"                  />
    <objectDestructor name="linearGrowth_"                         />
    <objectDestructor name="correlationFunctionTwoPoint_"          />
    <objectDestructor name="cosmologicalVelocityField_"            />
    <objectDestructor name="cosmologicalMassVariance_"             />
    <objectDestructor name="mergerTreeBranchingProbability_"       />
    <objectDestructor name="haloEnvironment_"                      />
    <objectDestructor name="cosmologicalMassVarianceEnvironmental_"/>
    <objectDestructor name="criticalOverdensityEnvironmental_"     />
    <objectDestructor name="haloMassFunctionEnvironmental_"        />
    !!]
    deallocate(tree)
    !$omp end parallel
    call displayCounterClear(verbosityLevelWorking)
    ! Transfer tabulated results to self.
    call move_alloc(velocityRadialMeanVirial            ,self%velocityRadialMeanVirial            )
    call move_alloc(velocityRadialDispersionVirial      ,self%velocityRadialDispersionVirial      )
    call move_alloc(velocityTangentialMeanVirial        ,self%velocityTangentialMeanVirial        )
    call move_alloc(velocityTangentialDispersionVirial  ,self%velocityTangentialDispersionVirial  )
    call move_alloc(velocityRadialDistributionOrbits    ,self%velocityRadialDistributionOrbits    )
    call move_alloc(velocityTangentialDistributionOrbits,self%velocityTangentialDistributionOrbits)
    call move_alloc(velocityDistributionOrbits          ,self%velocityDistributionOrbits          )
    call move_alloc(velocityTotalRMS                    ,self%velocityTotalRMS                    )
    call move_alloc(velocityDistributionPeak            ,self%velocityDistributionPeak            )
    ! Store data to file.
    call self%storeTable()
    return

  contains

    double precision function integrandEnvironment(overdensity)
      use :: Cosmology_Parameters, only : hubbleUnitsLittleH
      use :: Calculations_Resets , only : Calculations_Reset
      implicit none
      double precision, intent(in   ) :: overdensity
      double precision, parameter     :: radiusReferenceExternal=10.0d0 ! Reference radius (in units of Mpc/h) in the Sheth & Diaferio (2001) model for the environmental dependence of velocity dispersion.
      double precision                :: mu                            , massRadius     , &
           &                             massReference                 , radiusReference, &
           &                             massFunction                  , mergerRate     , &
           &                             overdensityNonlinear
      
      call haloEnvironment_%overdensityLinearSet(nodeHost_,overdensity)
      call Calculations_Reset                   (nodeHost_            )
      radiusReference     =+radiusReferenceExternal                                       &
           &               /cosmologyParameters_%hubbleConstant      (hubbleUnitsLittleH)
      massRadius          =+haloEnvironment_    %environmentMass     (                  )
      massReference       =+4.0d0                                                         &
           &               *Pi                                                            &
           &               /3.0d0                                                         &
           &               *cosmologyFunctions_ %matterDensityEpochal(timeEvaluate_     ) &
           &               *radiusReference**3
      mu                  =+0.6d0                                                                       &
           &               *(                                                                           &
           &                 +cosmologicalMassVariance_%rootVariance(mass=massRadius   ,time=self%time) &
           &                 /cosmologicalMassVariance_%rootVariance(mass=massReference,time=self%time) &
           &                )**2
      overdensityNonlinear=+haloEnvironment_               %overdensityNonlinear(                                  node=nodeHost_                                                                                                   )
      massFunction        =+haloMassFunctionEnvironmental_ %differential        (time=timeEvaluate_,mass=massHost_,node=nodeHost_                                                                                                   )
      mergerRate          =+mergerTreeBranchingProbability_%rate                (time=timeEvaluate_,mass=massHost_,node=nodeHost_,deltaCritical=criticalOverdensity_%value(time=self%time,mass=massHost_),massBranch=0.5d0*massHost_)
      integrandEnvironment=+(                                                                                                                                                                                                         &
           &                 +1.0d0                                                                                                                                                                                                   &
           &                 +overdensityNonlinear                                                                                                                                                                                    &
           &                )**mu                                                                                                                                                                                                     &
           &               *massFunction                                                                                                                                                                                              &
           &               *mergerRate                                                                                                                                                                                                &
           &               *haloEnvironment_               %pdf                 (overdensity                                                                                                                                        )
      return
    end function integrandEnvironment
    
    double precision function integrandEnvironmentNormalization(overdensity)
      !!{
      The normalizing integrand for the environmental dependence of velocity dispersion.
      !!}
      use :: Calculations_Resets, only : Calculations_Reset
      implicit none
      double precision, intent(in   ) :: overdensity
      double precision                :: mergerRate          , massFunction, &
           &                             overdensityNonlinear

      call haloEnvironment_%overdensityLinearSet(nodeHost_,overdensity)
      call Calculations_Reset                   (nodeHost_            )
      overdensityNonlinear             =+haloEnvironment_               %overdensityNonlinear(                                  node=nodeHost_                                                                                                   )
      massFunction                     =+haloMassFunctionEnvironmental_ %differential        (time=timeEvaluate_,mass=massHost_,node=nodeHost_                                                                                                   )
      mergerRate                       =+mergerTreeBranchingProbability_%rate                (time=timeEvaluate_,mass=massHost_,node=nodeHost_,deltaCritical=criticalOverdensity_%value(time=self%time,mass=massHost_),massBranch=0.5d0*massHost_)
      integrandEnvironmentNormalization=+massFunction                                                                                                                                                                                              &
           &                            *mergerRate                                                                                                                                                                                                &
           &                            *haloEnvironment_               %pdf                 (overdensity                                                                                                                                        )
      return
    end function integrandEnvironmentNormalization
    
  end subroutine lossConeTabulate

  double precision function timeAlongOrbit(radius,radiusApocenter,radiusPericenter,velocityTangentialVirial)
    !!{
    Compute the time taken along the orbit specified by the pericenter radius, {\normalfont \ttfamily radiusPericenter}, and the
    tangential velocity at the virial radius, {\normalfont \ttfamily velocityTangentialVirial}, to travel from the pericenter to
    the given radius, {\normalfont \ttfamily radius}. All quantities are in virial units. Writing
    \begin{equation}
     v^\prime_\mathrm{r}(r) = \left( -{2 \over r_\mathrm{p}} + {2 \over r} + {v_\theta^2 \over r_\mathrm{p}^2} - {v_\theta^2 \over r^2} \right)^{1/2}
    \end{equation}
    where $v^\prime_\mathrm{r}(r)$ is the radial velocity at radius $r$, $r_\mathrm{p}$ is the apocentric radius, and $v_\theta$ is
    the tangential velocity at the virial radius, the time of flight along the orbit is
    \begin{equation}
     t(r) = \int_{r_\mathrm{p}}^r {\mathrm{d}r \over v^\prime_\mathrm{r}(r)}.
    \end{equation}
    This was evaluated using Mathematica to give:
    \begin{equation}
     t(r) = \frac{\sqrt{2 r_\mathrm{p}-v_\theta^2} \left(r_\mathrm{p}^2 \left(v_\theta^2-2 r\right)+2 r_\mathrm{p} r^2-r^2 v_\theta^2\right)
            -2 r_\mathrm{p}^2 \sqrt{r_\mathrm{p}-r} \sqrt{-r_\mathrm{p} \left(r_\mathrm{p}-v_\theta^2\right)} \sqrt{-\frac{-2 r_\mathrm{p} r+r_\mathrm{p}
            v_\theta^2+r v_\theta^2}{r_\mathrm{p}^2-r_\mathrm{p} v_\theta^2}} \sinh ^{-1}\left(\frac{\sqrt{r_\mathrm{p}-r} \sqrt{2 r_\mathrm{p}-v_\theta^2}}{\sqrt{2}
            \sqrt{-r_\mathrm{p} \left(r_\mathrm{p}-v_\theta^2\right)}}\right)}{r \left(2 r_\mathrm{p}-v_\theta^2\right)^{3/2} \sqrt{\frac{v_\theta^2}{r_\mathrm{p}^2}
            -\frac{2}{r_\mathrm{p}}+\frac{2 r-v_\theta^2}{r^2}}}.
    \end{equation}
    Note that this expression works for unbound orbits.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radiusApocenter          , velocityTangentialVirial , &
         &                             radiusPericenter         , radius
    double precision, parameter     :: timeInfinite     =1.0d100
    double complex                  :: radiusPericenter_        , velocityTangentialVirial_, &
         &                             timeAlongOrbit_          , radius_
    
    if (radius == radiusPericenter) then
       ! The time at the pericenter is zero by construction.
       timeAlongOrbit=0.0d0
    else if (radiusApocenter > 0.0d0 .and. radius >= radiusApocenter) then
       ! The orbit is bound. The time to the apocenter is half the orbital period, which we can determine from Kepler's law
       ! using the semi-major axis.
       timeAlongOrbit=Pi*(0.5d0*(radiusPericenter+radiusApocenter))**1.5d0       
    else
       ! Calculations are performed using complex variables (as some of the square roots give imaginary results), and then the real
       ! component is taken at the end.
       radius_                  =radius
       radiusPericenter_        =radiusPericenter
       velocityTangentialVirial_=velocityTangentialVirial
       timeAlongOrbit_          =+(                                                                                                                                                                                                                &
            &                               +sqrt(+2.0d0*radiusPericenter_           -           velocityTangentialVirial_**2                                                                  )                                                   &
            &                      *             (+2.0d0*radiusPericenter_*radius_**2-radius_**2*velocityTangentialVirial_**2+radiusPericenter_**2*(-2.0d0*radius_+velocityTangentialVirial_**2))                                                  &
            &                      -               2.0d0*radiusPericenter_        **2                                                                                                                                                              &
            &                      *         sqrt(radiusPericenter_-radius_)                                                                                                                                                                       &
            &                      *         sqrt(-(+radiusPericenter_*(radiusPericenter_-velocityTangentialVirial_**2)))                                                                                                                          &
            &                      *         sqrt(-((-2.0d0*radiusPericenter_*radius_+radiusPericenter_*velocityTangentialVirial_**2+radius_*velocityTangentialVirial_**2)/(radiusPericenter_**2-radiusPericenter_*velocityTangentialVirial_**2))) &
            &                      *asinh(                                                                                                                                                                                                         &
            &                             +(                                                                                                                                                                                                       &
            &                               +sqrt(+      radiusPericenter_-radius_                     )                                                                                                                                           &
            &                               *sqrt(+2.0d0*radiusPericenter_-velocityTangentialVirial_**2)                                                                                                                                           &
            &                              )                                                                                                                                                                                                       &
            &                             /(                                                                                                                                                                                                       &
            &                               +sqrt(+2.0d0                                                             )                                                                                                                             &
            &                               *sqrt(-(radiusPericenter_*(radiusPericenter_-velocityTangentialVirial_**2)))                                                                                                                           &
            &                              )                                                                                                                                                                                                       &
            &                            )                                                                                                                                                                                                         &
            &                     )                                                                                                                                                                                                                &
            &                    /(                                                                                                                                                                                                                &
            &                      +                     radius_                                                                                                                                                                                   &
            &                      *             (+2.0d0*radiusPericenter_-velocityTangentialVirial_**2)**1.5d0                                                                                                                                    &
            &                      *         sqrt(-2.0d0/radiusPericenter_+velocityTangentialVirial_**2/radiusPericenter_**2+(2.0d0*radius_-velocityTangentialVirial_**2)/radius_**2)                                                              &
            &                     )
       timeAlongOrbit          =-real(timeAlongOrbit_)
    end if
    return
  end function timeAlongOrbit

  
  subroutine lossConeRestoreTable(self)
    !!{
    Attempt to restore a table from file.
    !!}
    use :: File_Utilities    , only : File_Exists, File_Lock, File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char
    implicit none
    class           (virialOrbitLossCone)             , intent(inout) :: self
    type            (hdf5Object                                      )                             :: file
    type            (lockDescriptor                                  )                             :: fileLock
    !$GLC attributes unused :: self

    if (.not.self%fileRead.and.File_Exists(self%fileName)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
       if (allocated(self%mass)) then
          deallocate(self%mass                                )
          deallocate(self%velocity                            )
          deallocate(self%velocityRadialMeanVirial            )
          deallocate(self%velocityRadialDispersionVirial      )
          deallocate(self%velocityTangentialMeanVirial        )
          deallocate(self%velocityTangentialDispersionVirial  )
          deallocate(self%velocityRadialDistributionOrbits    )
          deallocate(self%velocityTangentialDistributionOrbits)
          deallocate(self%velocityDistributionOrbits          )
          deallocate(self%interpolatorMass                    )
       end if
       !$ call hdf5Access%set()
       file=hdf5Object(char(self%fileName))
       call file%readAttribute('time'                                ,     self%time                                 )
       call file%readAttribute('massMinimum'                         ,     self%massMinimum                          )
       call file%readAttribute('massMaximum'                         ,     self%massMaximum                          )
       call file%readDataset  ('mass'                                ,     self%mass                                 )
       call file%readDataset  ('velocityRadialMeanVirial'            ,     self%velocityRadialMeanVirial             )
       call file%readDataset  ('velocityRadialDispersionVirial'      ,     self%velocityRadialDispersionVirial       )
       call file%readDataset  ('velocityTangentialMeanVirial'        ,     self%velocityTangentialMeanVirial         )
       call file%readDataset  ('velocityTangentialDispersionVirial'  ,     self%velocityTangentialDispersionVirial   )
       call file%readDataset  ('velocityRadialDistributionOrbits'    ,     self%velocityRadialDistributionOrbits     )
       call file%readDataset  ('velocityTangentialDistributionOrbits',     self%velocityTangentialDistributionOrbits )
       call file%readDataset  ('velocityDistributionOrbits'          ,     self%velocityDistributionOrbits           )
       call file%readDataset  ('velocityTotalRMS'                    ,     self%velocityTotalRMS                     )
       call file%readDataset  ('velocityDistributionPeak'            ,     self%velocityDistributionPeak             )
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%fileRead=.true.
       ! Rebuild interpolator.
       self%interpolatorMass=interpolator(log(self%mass))
    end if
    return
  end subroutine lossConeRestoreTable

  subroutine lossConeStoreTable(self)
    !!{
    Store the tabulated solution to file.
    !!}
    use :: File_Utilities    , only : Directory_Make, File_Lock, File_Path, File_Unlock, &
          &                           lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char
    implicit none
    class(virialOrbitLossCone), intent(inout) :: self
    type (hdf5Object         )                :: file
    type (lockDescriptor     )                :: fileLock

    call Directory_Make(File_Path(self%fileName))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock     (self%fileName,fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    file=hdf5Object(char(self%fileName),overWrite=.true.,readOnly=.false.)
    call file%writeAttribute(self%time                                ,'time'                                )
    call file%writeAttribute(self%massMinimum                         ,'massMinimum'                         )
    call file%writeAttribute(self%massMaximum                         ,'massMaximum'                         )
    call file%writeDataset  (self%mass                                ,'mass'                                )
    call file%writeDataset  (self%velocityRadialMeanVirial            ,'velocityRadialMeanVirial'            )
    call file%writeDataset  (self%velocityRadialDispersionVirial      ,'velocityRadialDispersionVirial'      )
    call file%writeDataset  (self%velocityTangentialMeanVirial        ,'velocityTangentialMeanVirial'        )
    call file%writeDataset  (self%velocityTangentialDispersionVirial  ,'velocityTangentialDispersionVirial'  )
    call file%writeDataset  (self%velocityRadialDistributionOrbits    ,'velocityRadialDistributionOrbits'    )
    call file%writeDataset  (self%velocityTangentialDistributionOrbits,'velocityTangentialDistributionOrbits')
    call file%writeDataset  (self%velocityDistributionOrbits          ,'velocityDistributionOrbits'          )
    call file%writeDataset  (self%velocityTotalRMS                    ,'velocityTotalRMS'                    )
    call file%writeDataset  (self%velocityDistributionPeak            ,'velocityDistributionPeak'            )
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    ! Mark the file as read, to avoid re-reading it later.
    self%fileRead=.true.
    return
  end subroutine lossConeStoreTable

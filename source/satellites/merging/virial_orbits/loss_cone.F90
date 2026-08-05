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

!+    Contributions to this file made by: Andrew Benson. The fixes for issue #1321, and the pinning of the mass tabulation to an
!+    absolute lattice for issue #1317, were diagnosed and drafted with assistance from Claude, and reviewed and verified by
!+    Andrew Benson.

  !!{RST
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
  use :: Numerical_Ranges               , only : rangeLattice
  use :: Correlation_Functions_Two_Point, only : correlationFunctionTwoPointClass
  use :: Galacticus_Nodes               , only : treeNode
  use :: Merger_Tree_Branching          , only : mergerTreeBranchingProbabilityClass

  !![
  <virialOrbit name="virialOrbitLossCone" docformat="rst">
   <description>
   A virial orbit class that draws satellite infall orbital parameters using a loss cone model, accounting for the depletion of nearly radial orbits due to merging. The velocity range and resolution of the orbital distribution grid are controlled by the ``[velocityMinimum]``, ``[velocityMaximum]``, ``[velocitiesPerUnit]``, and ``[massesPerDecade]`` parameters.
   </description>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitLossCone
     !!{RST
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
          &                                                                                    time                                     , velocityDispersionMultiplier
     integer                                                                                :: countMassesPerDecade                     , countVelocitiesPerUnit
     type            (rangeLattice                       )                                  :: latticeMass
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
     <methods docformat="rst">
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
     !!{RST
     Constructors for the :galacticus-class:`virialOrbitLossCone` virial orbits class.
     !!}
     module procedure lossConeConstructorParameters
     module procedure lossConeConstructorInternal
  end interface virialOrbitLossCone

  ! Default range of halo masses over which to tabulate. This seeds the pinned range: any tabulation, whatever masses prompted
  ! it, spans at least this range, which guarantees that any two tabulations overlap and so that one can always be extended to
  ! cover the other. These are whole decades, and so are themselves anchor points of the lattice for any density of points per
  ! decade.
  double precision, parameter :: massTableMinimum=1.0d06, massTableMaximum=1.0d15

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
    !!{RST
    Constructor for the :galacticus-class:`virialOrbitLossCone` virial orbits class which takes a parameter set as input.
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
    <inputParameter docformat="rst">
      <name>velocityMinimum</name>
      <source>parameters</source>
      <description>
      The minimum velocity (in units of the host virial velocity) for which to compute velocity distributions.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>velocityMaximum</name>
      <source>parameters</source>
      <description>
      The maximum velocity (in units of the host virial velocity) for which to compute velocity distributions.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>countVelocitiesPerUnit</name>
      <source>parameters</source>
      <description>
      The number of points per unit of velocity (in units of the host virial velocity) for which to compute velocity distributions.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>countMassesPerDecade</name>
      <source>parameters</source>
      <description>
      The number of points per decade of mass for which to compute infall properties.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>includeInFlightGrowth</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, linear growth of the velocity field during the flight of the secondary halo to the primary is included.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>haloMassFunctionA</name>
      <source>parameters</source>
      <defaultValue>0.707d0</defaultValue>
      <description>
      The parameter :math:`a` of the :cite:t:`sheth_ellipsoidal_2001` halo mass function used in averaging over environment.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>haloMassFunctionP</name>
      <source>parameters</source>
      <defaultValue>0.300d0</defaultValue>
      <description>
      The parameter :math:`p` of the :cite:t:`sheth_ellipsoidal_2001` halo mass function used in averaging over environment.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>haloMassFunctionNormalization</name>
      <source>parameters</source>
      <defaultValue>0.322d0</defaultValue>
      <description>
      The normalization parameter :math:`A` of the :cite:t:`sheth_ellipsoidal_2001` halo mass function used in averaging over environment.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>velocityDispersionMultiplier</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>
      A multiplier applied to the dispersion of the cosmological velocity field.
      </description>
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
    !!{RST
    Internal constructor for the :galacticus-class:`virialOrbitLossCone` virial orbits class.
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

    ! Set an unphysical initial time, so that tabulation is forced on the first call. The range of masses is left undefined here:
    ! it is established on that first call, seeded with the default range of "massTableMinimum" to "massTableMaximum".
    self%time=-huge(0.0d0)
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
    !!{RST
    Destructor for the :galacticus-class:`virialOrbitLossCone` virial orbits class.
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
    !!{RST
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
    !!{RST
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
    !!{RST
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
    !!{RST
    Return a virial density contrast object defining that used in the calculation of orbital parameters.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: lossConeDensityContrastDefinition
    class(virialOrbitLossCone       ), intent(inout) :: self

    lossConeDensityContrastDefinition => self%virialDensityContrast_
    return
  end function lossConeDensityContrastDefinition

  double precision function lossConeVelocityTangentialMagnitudeMean(self,node,host)
    !!{RST
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
    !!{RST
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
    !!{RST
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
    !!{RST
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
    !!{RST
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
    !!{RST
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
    !!{RST
    Compute properties of infalling halos.
    !!}
    use :: Cosmological_Density_Field          , only : haloEnvironmentNormal
    use :: Dark_Matter_Profile_Mass_Definitions, only : Dark_Matter_Profile_Mass_Definition
    use :: Display                             , only : displayCounter                     , displayCounterClear , displayIndent, displayUnindent     , &
          &                                             verbosityLevelWorking
    use :: Error                               , only : errorStatusSuccess
    use :: Galacticus_Nodes                    , only : mergerTree                         , nodeComponentBasic  , treeNode
    use :: Calculations_Resets                 , only : Calculations_Reset
    use :: Numerical_Constants_Math            , only : Pi
    use :: Numerical_Ranges                    , only : Make_Range                         , rangeTypeLogarithmic, Range_Pinned , Range_Lattice_Offset, &
          &                                             gridSchemePerDecade
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
    logical                                         , dimension(:,:    ), allocatable :: isComputed
    type            (rangeLattice                  )                                  :: latticeMass
    integer                                                                           :: offsetMass                                     , countMassesPrevious
    logical                                                                           :: reuseSolutions
    double precision                                , dimension(2      )              :: radiiEvaluation
    double precision                                , parameter                       :: velocityRadialInfallMaximum           =6.0000d0
    double precision                                , parameter                       :: velocityRadialInfallStep              =0.0001d0
    double precision                                , parameter                       :: radiusTableMinimum                    =1.0d-2
    double precision                                , parameter                       :: radiusTableMaximum                    =1.0d+3
    integer                                         , parameter                       :: radiusTablePointsPerDecade            =20
    double precision                                , parameter                       :: overdensityLimitLower                 =-10.0d0
    ! The lower limit of the environmental overdensity integrals below is nominally "overdensityLimitLower" standard deviations
    ! of the environmental density field. That can reach beyond the range over which the halo environment class is able to map
    ! linear to nonlinear overdensity (which is tabulated by its spherical collapse solver over a fixed range of linear
    ! overdensity), so the limit is clamped to the minimum of that range. Since the integrands are weighted by the (Gaussian)
    ! environmental PDF, the truncated tail contributes negligibly, and the same factor appears in the normalizing integral.
    double precision                                , parameter                       :: overdensityLinearMinimumMap           = -5.0d0
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
         &                                                                               radiusEvaluateLagrangian                       , overdensityEnvironmentMinimum             , &
         &                                                                               energyOrbitDoubled
    type            (interpolator                  )                    , allocatable :: interpolatorVelocityDispersionLinear
    type            (integrator                    )                    , allocatable :: integratorEnvironment                          , integratorEnvironmentNormalizer
    
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
    ! Find the range of masses required, pinned to an absolute lattice so that the tabulation---and hence every value
    ! interpolated from it---is independent of the masses at which it happened to be first requested, and so that it can be
    ! extended without recomputing any solution already found. Only the mass axis is pinned: the velocity axis is fixed by the
    ! parameters of this object and never changes, while the epoch is a key rather than an axis---the tabulation is discarded
    ! wholesale when it changes, since every tabulated quantity depends on it.
    !
    ! Both the raw node masses and those same masses under our own definition of the virial density contrast are covered. The
    ! latter are the coordinates at which the tables are actually interpolated (see "lossConeInterpolants"), while the former are
    ! what the range of the tabulation was historically tested against here; the two differ by the mass-definition conversion, so
    ! covering both leaves the tabulation no narrower than it was before being pinned. The interpolator built below aborts on
    ! extrapolation, so a mass falling outside the tabulated range is fatal rather than merely inaccurate.
    latticeMass=Range_Pinned(                                                                                 &
         &                                   [basicHost%mass(),basicSatellite%mass(),massHost,massSatellite], &
         &                                    self%countMassesPerDecade                                     , &
         &                                    gridSchemePerDecade                                           , &
         &                    rangeCurrent  =[massTableMinimum,massTableMaximum]                            , &
         &                    latticeCurrent=self%latticeMass                                                 &
         &                   )
    if     (                                                   &
         &   basicHost       %time  (           ) == self%time &
         &  .and.                                              &
         &   self%latticeMass%covers(latticeMass)              &
         & ) return
    ! Determine whether the solutions already found can be carried over into the new tabulation. They can only if the epoch is
    ! unchanged---every tabulated quantity depends on it---and if there is an existing tabulation, on a known lattice, to carry
    ! over. Where they can, they occupy a block of the new arrays offset by the difference in the lattices.
    reuseSolutions=  basicHost       %time     () == self%time &
         &         .and.                                       &
         &           self%latticeMass%isDefined()              &
         &         .and.                                       &
         &           allocated(self%velocityRadialMeanVirial)
    countMassesPrevious=0
    offsetMass         =0
    if (reuseSolutions) then
       countMassesPrevious=self%latticeMass%count
       offsetMass         =Range_Lattice_Offset(self%latticeMass,latticeMass)
    end if
    ! Set time for this tabulation.
    self%time=basicHost%time()
    ! Get number of velocities to tabulate. Note that "velocity" and "interpolatorVelocity" are built once, by the constructor,
    ! from parameters that do not change between tabulations, and are deliberately left alone here---they are used below, and by
    ! the interpolants, throughout this and subsequent tabulations.
    countVelocities=size(self%velocity)
    ! Build the range of masses. The abscissae are taken from the lattice, rather than by subdividing the range, so that they are
    ! bit-identical to those of any other tabulation built on the same lattice.
    countMasses=latticeMass%count
    if (allocated(self%mass            )) deallocate(self%mass            )
    if (allocated(self%interpolatorMass)) deallocate(self%interpolatorMass)
    allocate(self%mass            (countMasses))
    allocate(self%interpolatorMass             )
    self%mass            =latticeMass%values()
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
    allocate(isComputed                          (countMasses,countMasses                                ))
    velocityRadialMeanVirial                     =0.0d0
    velocityRadialDispersionVirial               =0.0d0
    velocityTangentialMeanVirial                 =0.0d0
    velocityTangentialDispersionVirial           =0.0d0
    velocityRadialDistributionOrbits             =0.0d0
    velocityTangentialDistributionOrbits         =0.0d0
    velocityDistributionOrbits                   =0.0d0
    velocityTotalRMS                             =0.0d0
    velocityDistributionPeak                     =0.0d0
    isComputed                                   =.false.
    ! Carry over the solutions already found. Both mass dimensions of every array are indexed on the same lattice, so the
    ! surviving solutions form a square block offset equally along both.
    if (reuseSolutions) then
       velocityRadialMeanVirial            (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious    )=self%velocityRadialMeanVirial
       velocityRadialDispersionVirial      (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious    )=self%velocityRadialDispersionVirial
       velocityTangentialMeanVirial        (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious    )=self%velocityTangentialMeanVirial
       velocityTangentialDispersionVirial  (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious    )=self%velocityTangentialDispersionVirial
       velocityRadialDistributionOrbits    (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious,:  )=self%velocityRadialDistributionOrbits
       velocityTangentialDistributionOrbits(offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious,:  )=self%velocityTangentialDistributionOrbits
       velocityDistributionOrbits          (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious,:,:)=self%velocityDistributionOrbits
       velocityTotalRMS                    (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious    )=self%velocityTotalRMS
       velocityDistributionPeak            (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious    )=self%velocityDistributionPeak
       isComputed                          (offsetMass+1:offsetMass+countMassesPrevious,offsetMass+1:offsetMass+countMassesPrevious    )=.true.
    end if
    ! Record the lattice on which this tabulation is built.
    self%latticeMass=latticeMass
    ! Iterate over host masses. Only those solutions not carried over need be found, and then only for satellites no more massive
    ! than their host.
    countTotal=0
    do iHost=1,countMasses
       countTotal=countTotal+count(.not.isComputed(iHost,1:iHost))
    end do
    countProgress=0
    !$omp parallel private(iHost,iSatellite,massHost,massSatellite,tree,nodeHost,nodeSatellite,basicHost,basicSatellite,radiusVirialHost,velocityVirialHost,velocityDispersionLinear,indexVelocityRadial,indexVelocityTangential,velocityRadialVirial,velocityTangentialVirial,countVelocityRadialInfall,indexVelocityRadialInfall,velocityRadialInfall,radiusInfallTerm1,radiusInfallTerm2,radiiEvaluation,iEvaluate,radiusEvaluateVirial,timeEvaluate,radiusApocenterVirial,radiusPericenterVirial,timeOfFlightVirial,timeOfFlight,radiusEvaluate,radiusEvaluateComoving,velocityDispersionEvaluate,velocityRadialMeanEvaluate,velocityDispersionRadialEvaluateVirial,velocityDispersionTangentialEvaluateVirial,velocityMeanRadialEvaluateVirial,velocityTangentialInfall,jacobianFactor,jacobianDeterminant,distributionFunction,interpolatorVelocityDispersionLinear,integratorEnvironment,integratorEnvironmentNormalizer,factorEnvironmental,massEnvironment,radiusEnvironment,overdensityEnvironmentMinimum,energyOrbitDoubled)
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
    !$omp barrier
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
       &amp;                                      factorMassEnvironment    =1.0d0                                      , &amp;
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
    ! Allocate the per-thread integrators and interpolator. These are declared "allocatable" and allocated here---rather than
    ! being plain (non-allocatable) variables in the "private" clause above---because gfortran does not apply default
    ! initialization to the private copy of a derived type, contrary to OpenMP, which requires the new list item to be
    ! initialized as if it had been locally declared without an initializer. Their reference-counting "resourceManager"
    ! components would therefore start filled with stack garbage, and the first assignment to them (which finalizes the
    ! left-hand side) would release that garbage: either blocking forever in OMP_Set_Lock() on a bogus lock pointer, or
    ! writing through a bogus counter pointer (issue #1321). OpenMP *does* guarantee that private copies of allocatables begin
    ! unallocated, and ALLOCATE applies default initialization, so this route is safe.
    allocate(integratorEnvironment               )
    allocate(integratorEnvironmentNormalizer     )
    allocate(interpolatorVelocityDispersionLinear)
    integratorEnvironment          =integrator(integrand=integrandEnvironment             ,toleranceRelative=1.0d-3)
    integratorEnvironmentNormalizer=integrator(integrand=integrandEnvironmentNormalization,toleranceRelative=1.0d-3)
    ! Find the lower limit for the integrals over environmental overdensity, clamped to the range over which the halo environment
    ! is able to map linear to nonlinear overdensity (see "overdensityLinearMinimumMap" above).
    overdensityEnvironmentMinimum=max(                                                                              &
         &                            +overdensityLimitLower                                                        &
         &                            *cosmologicalMassVariance_%rootVariance(mass=massEnvironment,time=self%time), &
         &                            +overdensityLinearMinimumMap                                                  &
         &                           )
    do iHost=1,countMasses
       ! Skip this host mass entirely if the solution for every satellite mass paired with it was carried over---this avoids
       ! rebuilding the host node and re-evaluating its environmental boost factor, as well as the orbit distributions
       ! themselves. Note that "isComputed" is shared, and is not modified anywhere within this parallel region, so every thread
       ! takes this branch identically and the work-sharing constructs below remain matched across the team.
       if (all(isComputed(iHost,1:iHost))) cycle
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
       timeEvaluate_      =+self                           %time
       factorEnvironmental=+integratorEnvironment          %integrate(overdensityEnvironmentMinimum,haloEnvironment_%overdensityLinearMaximum()) &
            &              /integratorEnvironmentNormalizer%integrate(overdensityEnvironmentMinimum,haloEnvironment_%overdensityLinearMaximum())
       ! Iterate over satellite masses.
       do iSatellite=1,countMasses          
          ! Only consider satellites less (or equally) massive than their host.
          if (iSatellite > iHost         ) cycle
          ! Skip solutions carried over from a previous tabulation.
          if (isComputed(iHost,iSatellite)) cycle
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
                   !! Skip cases where the roots are complex, or where the quadratic is degenerate. Note that
                   !! "radiusInfallTerm1" must be excluded when it is exactly zero (not merely when it is negative): that is the
                   !! tangency at which the two roots coincide, and the determinant of the Jacobian evaluated below diverges
                   !! there as 1/√(radiusInfallTerm1). This is an integrable singularity of measure zero in the integral over
                   !! infall radial velocity, so neighboring steps of that integral capture its contribution.
                    if     (                            &
                        &    radiusInfallTerm1 <= 0.0d0 &
                        &   .or.                        &
                        &    radiusInfallTerm2 == 0.0d0 &
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
                      !! Twice the specific energy of the orbit. The turning points are the roots of "2 E r² + 2 r - v_θ² = 0",
                      !! so this is also the coefficient which vanishes for a marginally bound (parabolic) orbit, for which the
                      !! quadratic degenerates to a linear equation with the single root r = v_θ²/2 and no apocenter.
                      energyOrbitDoubled=-2.0d0+velocityRadialVirial**2+velocityTangentialVirial**2
                      if (energyOrbitDoubled == 0.0d0) then
                         radiusPericenterVirial=+0.5d0*velocityTangentialVirial**2
                         !! There is no apocenter; retain the convention that a negative apocentric radius denotes an unbound
                         !! orbit. The value is never used, as timeAlongOrbit() branches on the energy before reaching it.
                         radiusApocenterVirial =-huge(0.0d0)
                      else
                         radiusApocenterVirial =+(                                                                  &
                              &                   -2.0d0                                                            &
                              &                   -sqrt(4.0d0+4.0d0*velocityTangentialVirial**2*energyOrbitDoubled) &
                              &                  )                                                                  &
                              &                 /  2.0d0                                                            &
                              &                 /                                               energyOrbitDoubled
                         radiusPericenterVirial=+(                                                                  &
                              &                   -2.0d0                                                            &
                              &                   +sqrt(4.0d0+4.0d0*velocityTangentialVirial**2*energyOrbitDoubled) &
                              &                  )                                                                  &
                              &                 /  2.0d0                                                            &
                              &                 /                                               energyOrbitDoubled
                      end if
                      timeOfFlightVirial   =+abs(                                                                                                                               &
                           &                     +timeAlongOrbit(radiusEvaluateVirial,radiusApocenterVirial,radiusPericenterVirial,velocityTangentialVirial,energyOrbitDoubled) &
                           &                     -timeAlongOrbit(1.0d0               ,radiusApocenterVirial,radiusPericenterVirial,velocityTangentialVirial,energyOrbitDoubled) &
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
      !!{RST
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

  double precision function timeAlongOrbit(radius,radiusApocenter,radiusPericenter,velocityTangentialVirial,energyDoubled)
    !!{RST
    Compute the time taken along the orbit specified by the pericenter radius, ``radiusPericenter``, and the tangential velocity at the virial radius, ``velocityTangentialVirial``, to travel from the pericenter to the given radius, ``radius``. All quantities are in virial units. The argument ``energyDoubled`` is twice the specific orbital energy, :math:`2 E = -2 + v_\mathrm{r}^2 + v_\theta^2`, which selects between the bound, marginally bound, and unbound forms below.

    The general expression given below is singular in two limits which are, in general, sampled by the tabulation grid, and which are therefore treated separately using their (elementary) closed forms:

    * **Radial orbits** (:math:`v_\theta=0`). The pericentric radius is then zero---the orbit plunges to the origin---and the general expression diverges since it divides by :math:`r_\mathrm{p}`. The radial Kepler problem instead gives, for :math:`a` the semi-major axis,

      .. math::

         r = a (1-\cos \eta),  \,\,\, t = a^{3/2} (\eta - \sin \eta)

      when bound (:math:`a = -1/2E`), and

      .. math::

         r = a (\cosh \eta - 1),  \,\,\, t = a^{3/2} (\sinh \eta - \eta)

      when unbound (:math:`a = +1/2E`).

    * **Marginally bound (parabolic) orbits** (:math:`E=0`). Here :math:`2 r_\mathrm{p} - v_\theta^2 = -2 E r_\mathrm{p}^2 = 0`, so the general expression divides by zero. Barker's equation instead gives

      .. math::

         t(r) = \sqrt{2 r_\mathrm{p}^3} \left( D + D^3/3 \right), \,\,\, D = \sqrt{r/r_\mathrm{p} - 1},

      which reduces to the radial, marginally bound result :math:`t = \sqrt{2} r^{3/2}/3` as :math:`r_\mathrm{p} \rightarrow 0`.

    Writing

    .. math::

       v^\prime_\mathrm{r}(r) = \left( -{2 \over r_\mathrm{p}} + {2 \over r} + {v_\theta^2 \over r_\mathrm{p}^2} - {v_\theta^2 \over r^2} \right)^{1/2}

    where :math:`v^\prime_\mathrm{r}(r)` is the radial velocity at radius :math:`r`, :math:`r_\mathrm{p}` is the apocentric radius, and :math:`v_\theta` is the tangential velocity at the virial radius, the time of flight along the orbit is

    .. math::

       t(r) = \int_{r_\mathrm{p}}^r {\mathrm{d}r \over v^\prime_\mathrm{r}(r)}.

    This was evaluated using Mathematica to give:

    .. math::

       t(r) = \frac{\sqrt{2 r_\mathrm{p}-v_\theta^2} \left(r_\mathrm{p}^2 \left(v_\theta^2-2 r\right)+2 r_\mathrm{p} r^2-r^2 v_\theta^2\right)
              -2 r_\mathrm{p}^2 \sqrt{r_\mathrm{p}-r} \sqrt{-r_\mathrm{p} \left(r_\mathrm{p}-v_\theta^2\right)} \sqrt{-\frac{-2 r_\mathrm{p} r+r_\mathrm{p}
              v_\theta^2+r v_\theta^2}{r_\mathrm{p}^2-r_\mathrm{p} v_\theta^2}} \sinh ^{-1}\left(\frac{\sqrt{r_\mathrm{p}-r} \sqrt{2 r_\mathrm{p}-v_\theta^2}}{\sqrt{2}
              \sqrt{-r_\mathrm{p} \left(r_\mathrm{p}-v_\theta^2\right)}}\right)}{r \left(2 r_\mathrm{p}-v_\theta^2\right)^{3/2} \sqrt{\frac{v_\theta^2}{r_\mathrm{p}^2}
              -\frac{2}{r_\mathrm{p}}+\frac{2 r-v_\theta^2}{r^2}}}.

    Note that this expression works for unbound orbits.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: radiusApocenter          , velocityTangentialVirial , &
         &                             radiusPericenter         , radius                   , &
         &                             energyDoubled
    double precision, parameter     :: timeInfinite     =1.0d100
    double precision                :: semiMajorAxis            , anomaly
    double complex                  :: radiusPericenter_        , velocityTangentialVirial_, &
         &                             timeAlongOrbit_          , radius_

    if      (velocityTangentialVirial == 0.0d0) then
       ! A purely radial orbit. The pericentric radius is zero, so the general expression below---which divides by it---is
       ! singular. Use instead the elementary solution of the radial Kepler problem.
       if      (energyDoubled <  0.0d0) then
          ! Bound.
          semiMajorAxis =-1.0d0/energyDoubled
          !! Clamp the argument of the arccosine, which can exceed unity in magnitude by a rounding error when the radius
          !! reaches the apocenter.
          anomaly       =acos(max(-1.0d0,min(+1.0d0,1.0d0-radius/semiMajorAxis)))
          timeAlongOrbit=semiMajorAxis**1.5d0*(anomaly-sin(anomaly))
       else if (energyDoubled         == 0.0d0) then
          ! Marginally bound.
          timeAlongOrbit=sqrt(2.0d0)*radius**1.5d0/3.0d0
       else
          ! Unbound.
          semiMajorAxis =+1.0d0/energyDoubled
          anomaly       =acosh(1.0d0+radius/semiMajorAxis)
          timeAlongOrbit=semiMajorAxis**1.5d0*(sinh(anomaly)-anomaly)
       end if
    else if (energyDoubled            == 0.0d0) then
       ! A marginally bound (parabolic) orbit with non-zero angular momentum. Here 2 rₚ - v_θ² = -2 E rₚ² vanishes, so the general
       ! expression below divides by zero. Use Barker's equation instead.
       anomaly       =sqrt(max(0.0d0,radius/radiusPericenter-1.0d0))
       timeAlongOrbit=sqrt(2.0d0*radiusPericenter**3)*(anomaly+anomaly**3/3.0d0)
    else if (radius == radiusPericenter) then
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
    !!{RST
    Attempt to restore a table from file.
    !!}
    use :: File_Utilities    , only : File_Exists              , File_Lock          , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: ISO_Varying_String, only : char
    use :: Numerical_Ranges  , only : enumerationGridSchemeType, gridSchemePerDecade
    implicit none
    class  (virialOrbitLossCone), intent(inout) :: self
    type   (hdf5File           )                :: file
    type   (lockDescriptor     )                :: fileLock
    type   (rangeLattice       )                :: latticeCached
    integer                                     :: schemeCached      , pointsPerCached, &
         &                                         indexMinimumCached, countCached

    if (.not.self%fileRead.and.File_Exists(self%fileName)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
       ! As in tabulate(), "velocity" and "interpolatorVelocity" are constructor-built and are not restored from file, so they
       ! must not be freed here.
       if (allocated(self%mass)) then
          deallocate(self%mass                                )
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
       end if
       !$ call hdf5Access%set()
       ! Open read-only: opening read-write would have HDF5 take an exclusive lock on the file, so that another process reading
       ! it concurrently---which the shared file lock taken above explicitly permits---would fail to open it at all.
       file=hdf5File(self%fileName,overWrite=.false.,readOnly=.true.)
       call file%readAttribute('time'                                ,     self%time                                )
       call file%readAttribute('massGridScheme'                      ,     schemeCached                             )
       call file%readAttribute('massPointsPer'                       ,     pointsPerCached                          )
       call file%readAttribute('massIndexMinimum'                    ,     indexMinimumCached                       )
       call file%readAttribute('massCount'                           ,     countCached                              )
       call file%readDataset  ('mass'                                ,     self%mass                                )
       call file%readDataset  ('velocityRadialMeanVirial'            ,     self%velocityRadialMeanVirial            )
       call file%readDataset  ('velocityRadialDispersionVirial'      ,     self%velocityRadialDispersionVirial      )
       call file%readDataset  ('velocityTangentialMeanVirial'        ,     self%velocityTangentialMeanVirial        )
       call file%readDataset  ('velocityTangentialDispersionVirial'  ,     self%velocityTangentialDispersionVirial  )
       call file%readDataset  ('velocityRadialDistributionOrbits'    ,     self%velocityRadialDistributionOrbits    )
       call file%readDataset  ('velocityTangentialDistributionOrbits',     self%velocityTangentialDistributionOrbits)
       call file%readDataset  ('velocityDistributionOrbits'          ,     self%velocityDistributionOrbits          )
       call file%readDataset  ('velocityTotalRMS'                    ,     self%velocityTotalRMS                    )
       call file%readDataset  ('velocityDistributionPeak'            ,     self%velocityDistributionPeak            )
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       self%fileRead=.true.
       ! Place the restored tabulation on its lattice. If the file is not self-consistent---the datasets not matching the lattice
       ! recorded alongside them, or that lattice not one which this object could have built---then discard everything read from
       ! it and retabulate from scratch, rather than leave a partially-restored tabulation behind.
       latticeCached=rangeLattice(enumerationGridSchemeType(schemeCached),pointsPerCached,indexMinimumCached,countCached)
       if     (                                                        &
            &   latticeCached%isDefined  (                 )           &
            &  .and.                                                   &
            &   latticeCached%scheme      ==      gridSchemePerDecade  &
            &  .and.                                                   &
            &   pointsPerCached           == self%countMassesPerDecade &
            &  .and.                                                   &
            &   size(self%mass                      ) == countCached   &
            &  .and.                                                   &
            &   size(self%velocityRadialMeanVirial,1) == countCached   &
            &  .and.                                                   &
            &   size(self%velocityRadialMeanVirial,2) == countCached   &
            & ) then
          self%latticeMass=latticeCached
          ! Take the abscissae from the lattice rather than from the file, so that they are bit-identical to those of any other
          ! tabulation built on the same lattice.
          self%mass       =latticeCached%values()
          ! Rebuild interpolator. It must be allocated first: "interpolator" has a defined assignment, so assigning to it does
          ! not cause automatic allocation, and its assignment finalizes the left-hand side---which, if unallocated, means
          ! finalizing garbage. On this path the interpolator is always unallocated: it is either freed above, or was never
          ! built, since it is constructed only by tabulate().
          if (.not.allocated(self%interpolatorMass)) allocate(self%interpolatorMass)
          self%interpolatorMass=interpolator(log(self%mass))
       else
          self%time       =-huge(0.0d0)
          self%latticeMass= rangeLattice()
          deallocate(self%mass                                )
          deallocate(self%velocityRadialMeanVirial            )
          deallocate(self%velocityRadialDispersionVirial      )
          deallocate(self%velocityTangentialMeanVirial        )
          deallocate(self%velocityTangentialDispersionVirial  )
          deallocate(self%velocityRadialDistributionOrbits    )
          deallocate(self%velocityTangentialDistributionOrbits)
          deallocate(self%velocityDistributionOrbits          )
          deallocate(self%velocityTotalRMS                    )
          deallocate(self%velocityDistributionPeak            )
       end if
    end if
    return
  end subroutine lossConeRestoreTable

  subroutine lossConeStoreTable(self)
    !!{RST
    Store the tabulated solution to file.
    !!}
    use :: File_Utilities    , only : Directory_Make, File_Lock, File_Path, File_Unlock, &
          &                           lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File
    use :: ISO_Varying_String, only : char
    implicit none
    class(virialOrbitLossCone), intent(inout) :: self
    type (hdf5File           )                :: file
    type (lockDescriptor     )                :: fileLock

    call Directory_Make(File_Path(self%fileName))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock     (self%fileName,fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    file=hdf5File(self%fileName,overWrite=.true.,readOnly=.false.)
    call file%writeAttribute(self%time                                ,'time'                                )
    ! Record the lattice, not the bounds: the bounds follow from the lattice, but the converse does not, and it is the lattice
    ! which determines whether a tabulation read back later can be extended to cover a new range without recomputation.
    call file%writeAttribute(self%latticeMass%scheme%ID               ,'massGridScheme'                      )
    call file%writeAttribute(self%latticeMass%pointsPer               ,'massPointsPer'                       )
    call file%writeAttribute(self%latticeMass%indexMinimum            ,'massIndexMinimum'                    )
    call file%writeAttribute(self%latticeMass%count                   ,'massCount'                           )
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

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  use :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScaleClass
  use :: Virial_Orbits                     , only : virialOrbitClass
  use :: Numerical_Random_Numbers          , only : randomNumberGeneratorClass
  use :: Spherical_Collapse_Accretion_Flows, only : accretionFlowsClass
  use :: Merger_Tree_Branching             , only : mergerTreeBranchingProbabilityClass
  use :: Cosmological_Density_Field        , only : criticalOverdensityClass           , cosmologicalMassVarianceClass
  use :: Halo_Mass_Functions               , only : haloMassFunctionClass
  use :: Cosmology_Functions               , only : cosmologyFunctionsClass

  !![
  <task name="taskMergingHaloOrbitDistribution">
   <description>A task which tabulates the orbital parameter distribution for merging halos.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskMergingHaloOrbitDistribution
     !!{
     A task which tabulates the orbital parameter distribution for merging halos.
     !!}
     private
     class           (darkMatterHaloScaleClass           ), pointer :: darkMatterHaloScale_            => null()
     class           (virialOrbitClass                   ), pointer :: virialOrbit_                    => null()
     class           (randomNumberGeneratorClass         ), pointer :: randomNumberGenerator_          => null()
     class           (accretionFlowsClass                ), pointer :: accretionFlows_                 => null()
     class           (mergerTreeBranchingProbabilityClass), pointer :: mergerTreeBranchingProbability_ => null()
     class           (criticalOverdensityClass           ), pointer :: criticalOverdensity_            => null()
     class           (cosmologicalMassVarianceClass      ), pointer :: cosmologicalMassVariance_       => null()
     class           (haloMassFunctionClass              ), pointer :: haloMassFunction_               => null()
     class           (cosmologyFunctionsClass            ), pointer :: cosmologyFunctions_             => null()
     double precision                                               :: velocityMinimum                          , velocityMaximum       , &
          &                                                            massMinimum                              , massMaximum           , &
          &                                                            time
     integer                                                        :: countMassesPerDecade                     , countVelocitiesPerUnit
   contains
     final     ::            mergingHaloOrbitDistributionDestructor
     procedure :: perform => mergingHaloOrbitDistributionPerform
  end type taskMergingHaloOrbitDistribution

  interface taskMergingHaloOrbitDistribution
     !!{
     Constructors for the {\normalfont \ttfamily mergingHaloOrbitDistribution} task.
     !!}
     module procedure mergingHaloOrbitDistributionConstructorParameters
     module procedure mergingHaloOrbitDistributionConstructorInternal
  end interface taskMergingHaloOrbitDistribution

contains
  
  function mergingHaloOrbitDistributionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily mergingHaloOrbitDistribution} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes   , only : nodeClassHierarchyInitialize
    use :: Input_Parameters   , only : inputParameter              , inputParameters
    use :: Node_Components    , only : Node_Components_Initialize  , Node_Components_Thread_Initialize
    implicit none
    type            (taskMergingHaloOrbitDistribution   )                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (virialOrbitClass                   ), pointer       :: virialOrbit_
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class           (randomNumberGeneratorClass         ), pointer       :: randomNumberGenerator_
    class           (accretionFlowsClass                ), pointer       :: accretionFlows_
    class           (mergerTreeBranchingProbabilityClass), pointer       :: mergerTreeBranchingProbability_
    class           (criticalOverdensityClass           ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), pointer       :: cosmologicalMassVariance_
    class           (haloMassFunctionClass              ), pointer       :: haloMassFunction_
    type            (inputParameters                    ), pointer       :: parametersRoot
    double precision                                                     :: massMinimum                    , massMaximum           , &
         &                                                                  velocityMinimum                , velocityMaximum       , &
         &                                                                  redshift                       , time
    integer                                                              :: countMassesPerDecade           , countVelocitiesPerUnit

    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize     (parametersRoot)
       call Node_Components_Initialize       (parametersRoot)
       call Node_Components_Thread_Initialize(parametersRoot)
    else
       parametersRoot => null()
       call nodeClassHierarchyInitialize     (parameters    )
       call Node_Components_Initialize       (parameters    )
       call Node_Components_Thread_Initialize(parameters    )
    end if
    !![
    <inputParameter>
      <name>velocityMinimum</name>
      <source>parameters</source>
      <description>The minimum velocity (in units of the host virial velocity) for which to compute velocity distributions.</description>
      <type>real</type>
      <cardinality>0..1</cardinality>
    </inputParameter>
    <inputParameter>
      <name>velocityMaximum</name>
      <source>parameters</source>
      <description>The maximum velocity (in units of the host virial velocity) for which to compute velocity distributions.</description>
      <type>real</type>
      <cardinality>0..1</cardinality>
    </inputParameter>
    <inputParameter>
      <name>countVelocitiesPerUnit</name>
      <source>parameters</source>
      <description>The number of points per unit of velocity (in units of the host virial velocity) for which to compute velocity distributions.</description>
      <type>real</type>
      <cardinality>0..1</cardinality>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <description>The minimum mass halo for which to compute mergingHaloOrbitDistribution properties.</description>
      <type>real</type>
      <cardinality>0..1</cardinality>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <source>parameters</source>
      <description>The maximum mass halo for which to compute mergingHaloOrbitDistribution properties.</description>
      <type>real</type>
      <cardinality>0..1</cardinality>
    </inputParameter>
    <inputParameter>
      <name>countMassesPerDecade</name>
      <source>parameters</source>
      <description>The number of points per decade of mass for which to compute mergingHaloOrbitDistribution properties.</description>
      <type>real</type>
      <cardinality>0..1</cardinality>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <description>The redshift.</description>
      <type>real</type>
      <cardinality>0..1</cardinality>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"/>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    <objectBuilder class="virialOrbit"                    name="virialOrbit_"                    source="parameters"/>
    <objectBuilder class="randomNumberGenerator"          name="randomNumberGenerator_"          source="parameters"/>
    <objectBuilder class="accretionFlows"                 name="accretionFlows_"                 source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"            source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"       source="parameters"/>
    <objectBuilder class="haloMassFunction"               name="haloMassFunction_"               source="parameters"/>
    !!]
    time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    self=taskMergingHaloOrbitDistribution(time,velocityMinimum,velocityMaximum,countVelocitiesPerUnit,massMinimum,massMaximum,countMassesPerDecade,virialOrbit_,cosmologyFunctions_,darkMatterHaloScale_,accretionFlows_,mergerTreeBranchingProbability_,criticalOverdensity_,cosmologicalMassVariance_,haloMassFunction_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="virialOrbit_"                   />
    <objectDestructor name="randomNumberGenerator_"         />
    <objectDestructor name="accretionFlows_"                />
    <objectDestructor name="darkMatterHaloScale_"           />
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="haloMassFunction_"              />
    !!]
    return
  end function mergingHaloOrbitDistributionConstructorParameters

  function mergingHaloOrbitDistributionConstructorInternal(time,velocityMinimum,velocityMaximum,countVelocitiesPerUnit,massMinimum,massMaximum,countMassesPerDecade,virialOrbit_,cosmologyFunctions_,darkMatterHaloScale_,accretionFlows_,mergerTreeBranchingProbability_,criticalOverdensity_,cosmologicalMassVariance_,haloMassFunction_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily mergingHaloOrbitDistribution} task class.
    !!}
    implicit none
    type            (taskMergingHaloOrbitDistribution   )                        :: self
    class           (virialOrbitClass                   ), intent(in   ), target :: virialOrbit_
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (randomNumberGeneratorClass         ), intent(in   ), target :: randomNumberGenerator_
    class           (accretionFlowsClass                ), intent(in   ), target :: accretionFlows_
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class           (mergerTreeBranchingProbabilityClass), intent(in   ), target :: mergerTreeBranchingProbability_
    class           (criticalOverdensityClass           ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), intent(in   ), target :: cosmologicalMassVariance_
    class           (haloMassFunctionClass              ), intent(in   ), target :: haloMassFunction_
    double precision                                     , intent(in   )         :: massMinimum                    , massMaximum           , &
         &                                                                          velocityMinimum                , velocityMaximum       , &
         &                                                                          time
    integer                                              , intent(in   )         :: countMassesPerDecade           , countVelocitiesPerUnit
    !![
    <constructorAssign variables="time, velocityMinimum, velocityMaximum, countVelocitiesPerUnit, massMinimum, massMaximum, countMassesPerDecade, *virialOrbit_, *cosmologyFunctions_, *darkMatterHaloScale_, *accretionFlows_, *mergerTreeBranchingProbability_, *criticalOverdensity_, *cosmologicalMassVariance_, *haloMassFunction_, *randomNumberGenerator_"/>
    !!]

    return
  end function mergingHaloOrbitDistributionConstructorInternal

  subroutine mergingHaloOrbitDistributionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily mergingHaloOrbitDistribution} task class.
    !!}
    use :: Node_Components, only : Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
    implicit none
    type(taskMergingHaloOrbitDistribution), intent(inout) :: self
    
    !![
    <objectDestructor name="self%virialOrbit_"                   />
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%darkMatterHaloScale_"           />
    <objectDestructor name="self%randomNumberGenerator_"         />
    <objectDestructor name="self%accretionFlows_"                />
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    <objectDestructor name="self%criticalOverdensity_"           />
    <objectDestructor name="self%cosmologicalMassVariance_"      />
    <objectDestructor name="self%haloMassFunction_"              />
    !!]
    call Node_Components_Uninitialize       ()
    call Node_Components_Thread_Uninitialize()
    return
  end subroutine mergingHaloOrbitDistributionDestructor

  subroutine mergingHaloOrbitDistributionPerform(self,status)
    !!{
    Compute properties of mergingHaloOrbitDistributioning halos.
    !!}
    use :: Display                       , only : displayIndent                , displayUnindent     , displayCounter, displayCounterClear, &
         &                                        verbosityLevelWorking
    use :: Galacticus_Nodes              , only : treeNode                     , nodeComponentBasic  , mergerTree
    use :: Galacticus_HDF5               , only : galacticusOutputFile
    use :: Galacticus_Calculations_Resets, only : Galacticus_Calculations_Reset
    use :: IO_HDF5                       , only : hdf5Access                   , hdf5Object
    use :: Numerical_Ranges              , only : Make_Range                   , rangeTypeLogarithmic, rangeTypeLinear
    implicit none
    class           (taskMergingHaloOrbitDistribution), intent(inout)     , target      :: self
    integer                                           , intent(  out)     , optional    :: status
    type            (treeNode                        )                    , pointer     :: nodeHost                        , nodeSatellite
    class           (nodeComponentBasic              )                    , pointer     :: basicHost                       , basicSatellite
    type            (mergerTree                      )                    , pointer     :: tree
    double precision                                  , dimension(:      ), allocatable :: mass                            , velocity                            , &
         &                                                                                 separation                      , density                             , &
         &                                                                                 haloMassFunction
    double precision                                  , dimension(:,:    ), allocatable :: velocityRadialMeanVirial        , velocityRadialDispersionVirial      , &
         &                                                                                 velocityTangentialMeanVirial    , velocityTangentialDispersionVirial  , &
         &                                                                                 rateMerging
    double precision                                  , dimension(:,:,:  ), allocatable :: velocityRadialDistributionOrbits, velocityTangentialDistributionOrbits 
    double precision                                  , dimension(:,:,:,:), allocatable :: velocityDistributionOrbits
    type            (hdf5Object                      )                                  :: output
    integer                                                                             :: iHost                           , iSatellite                          , &
         &                                                                                 iVelocityRadial                 , iVelocityTangential                 , &
         &                                                                                 countMasses                     , countVelocities
    double precision                                                                    :: massHost                        , massSatellite                       , &
         &                                                                                 velocityRadialVirial            , velocityTangentialVirial            , &
         &                                                                                 distributionFunction            , distributionFunctionSum             , &
         &                                                                                 velocityWidthBin                , overdensityCriticalGrowthRate       , &
         &                                                                                 overdensityCritical             , rootVarianceLogarithmicGrowthRate   , &
         &                                                                                 barrierEffectiveGrowthRate
    
    call displayIndent('Begin task: merging halo orbit distribution')
    ! Build range of velocities.
    countVelocities=int(     (self%velocityMaximum-self%velocityMinimum)*dble(self%countVelocitiesPerUnit))+1    
    allocate(velocity        (countVelocities))
    velocity=Make_Range(self%velocityMinimum,self%velocityMaximum,countVelocities,rangeTypeLinear     )
    ! Build range of masses.
    countMasses    =int(log10(self%massMaximum    /self%massMinimum    )*dble(self%countMassesPerDecade  ))+1    
    allocate(mass            (countMasses    ))
    allocate(separation      (countMasses    ))
    allocate(density         (countMasses    ))
    allocate(haloMassFunction(countMasses    ))
    mass    =Make_Range(self%massMinimum    ,self%massMaximum    ,countMasses    ,rangeTypeLogarithmic)
    ! Allocate arrays for results and initialize.
    allocate(rateMerging                         (countMasses,countMasses                                ))
    allocate(velocityRadialMeanVirial            (countMasses,countMasses                                ))
    allocate(velocityRadialDispersionVirial      (countMasses,countMasses                                ))
    allocate(velocityTangentialMeanVirial        (countMasses,countMasses                                ))
    allocate(velocityTangentialDispersionVirial  (countMasses,countMasses                                ))
    allocate(velocityRadialDistributionOrbits    (countMasses,countMasses,countVelocities                ))
    allocate(velocityTangentialDistributionOrbits(countMasses,countMasses,countVelocities                ))
    allocate(velocityDistributionOrbits          (countMasses,countMasses,countVelocities,countVelocities))
    rateMerging                         =0.0d0
    velocityRadialMeanVirial            =0.0d0
    velocityRadialDispersionVirial      =0.0d0
    velocityTangentialMeanVirial        =0.0d0
    velocityTangentialDispersionVirial  =0.0d0
    velocityRadialDistributionOrbits    =0.0d0
    velocityTangentialDistributionOrbits=0.0d0
    velocityDistributionOrbits          =0.0d0
    ! Iterate over host masses.
    allocate(tree                                                        )
    allocate(tree%randomNumberGenerator_,mold=self%randomNumberGenerator_) 
    call tree%properties%initialize()
    !$omp critical(taskMergingHaloOrbitDistributionDeepCopy)
    !![
    <deepCopyReset variables="self%randomNumberGenerator_"/>
    <deepCopy source="self%randomNumberGenerator_" destination="tree%randomNumberGenerator_"/>
    <deepCopyFinalize variables="tree%randomNumberGenerator_"/>
    !!]
    !$omp end critical(taskMergingHaloOrbitDistributionDeepCopy)
    do iHost=1,countMasses
       ! Build host node.
       massHost           =  mass          (            iHost)
       nodeHost           => treeNode      (                 )
       basicHost          => nodeHost%basic(autoCreate=.true.)
       nodeHost %hostTree => tree
       tree     %baseNode => nodeHost
       call basicHost%massSet            (     massHost)
       call basicHost%timeSet            (self%time    )
       call basicHost%timeLastIsolatedSet(self%time    )
       call Galacticus_Calculations_Reset(     nodeHost)
       ! Compute the separation of the merging pairs (i.e. the virial radius of the primary halo), the density of the
       ! accretion flow at that separation, and the halo mass function.
       separation      (iHost)=self%darkMatterHaloScale_%virialRadius(nodeHost                   )
       density         (iHost)=self%accretionFlows_     %density     (nodeHost ,separation(iHost))
       haloMassFunction(iHost)=self%haloMassFunction_   %differential(self%time,massHost,nodeHost)
       ! Iterate over satellite masses.
       do iSatellite=1,countMasses
          ! Only consider satellites less (or equally) massive than their host.
          if (iSatellite > iHost) cycle
          ! Build satellite node.
          massSatellite           =  mass               (           iSatellite)
          nodeSatellite           => treeNode           (                     )
          basicSatellite          => nodeSatellite%basic(autoCreate=.true.    )
          nodeSatellite %hostTree => tree
          call basicSatellite%massSet            (     massSatellite)
          call basicSatellite%timeSet            (self%time         )
          call basicSatellite%timeLastIsolatedSet(self%time         )
          call Galacticus_Calculations_Reset     (     nodeSatellite)
          ! Compute critical overdensity and its growth rate for the host halo.
          overdensityCritical                 =+self%criticalOverdensity_     %value                              (time=basicHost%time(),mass=basicHost%mass(),node=nodeHost)
          overdensityCriticalGrowthRate       =+self%criticalOverdensity_     %gradientTime                       (time=basicHost%time(),mass=basicHost%mass(),node=nodeHost)
          ! Compute the growth rate of the mass variance. 
          rootVarianceLogarithmicGrowthRate   =+self%cosmologicalMassVariance_%rootVarianceLogarithmicGradientTime(time=basicHost%time(),mass=basicHost%mass()              )
          ! Compute absolute value of the rate of change of the effective barrier.
          barrierEffectiveGrowthRate          =+abs(                                               &
               &                                    +          overdensityCriticalGrowthRate       &
               &                                    -          overdensityCritical                 &
               &                                    *          rootVarianceLogarithmicGrowthRate   &
               &                                    /basicHost%time                             () &
               &                                   )
          ! Compute the merger rate.
          rateMerging(iHost,iSatellite)=+self%mergerTreeBranchingProbability_%rate                      (massHost+massSatellite,self%criticalOverdensity_%value(self%time),self%time,massSatellite,nodeHost) &
               &                        *                                     barrierEffectiveGrowthRate
          ! Iterate over velocities, and compute the distribution functions.
          do iVelocityRadial       =1,countVelocities
             do iVelocityTangential=1,countVelocities
                ! Extract velocities at the virial radius (in virial units).
                velocityRadialVirial    =velocity(iVelocityRadial    )
                velocityTangentialVirial=velocity(iVelocityTangential)
                ! Evaluate the distribution function.
                distributionFunction=self%virialOrbit_%velocityDistributionFunction(                                                                                  &
                     &                                                                                                                                nodeSatellite , &
                     &                                                                                                                                nodeHost      , &
                     &                                                              velocityRadialVirial    *self%darkMatterHaloScale_%virialVelocity(nodeHost     ), &
                     &                                                              velocityTangentialVirial*self%darkMatterHaloScale_%virialVelocity(nodeHost     )  &
                     &                                                             )
                ! Accumulate the distribution function, marginal distribution, and moments.
                velocityDistributionOrbits                  (iHost,iSatellite,iVelocityRadial,iVelocityTangential)= &
                     & +velocityDistributionOrbits          (iHost,iSatellite,iVelocityRadial,iVelocityTangential)  &
                     & +distributionFunction
                velocityRadialDistributionOrbits            (iHost,iSatellite,iVelocityRadial                    )= &
                     & +velocityRadialDistributionOrbits    (iHost,iSatellite,iVelocityRadial                    )  &
                     & +distributionFunction
                velocityTangentialDistributionOrbits        (iHost,iSatellite,                iVelocityTangential)= &
                     & +velocityTangentialDistributionOrbits(iHost,iSatellite,                iVelocityTangential)  &
                     & +distributionFunction
                velocityRadialMeanVirial                    (iHost,iSatellite                                    )= &
                     & +velocityRadialMeanVirial            (iHost,iSatellite                                    )  &
                     & +distributionFunction                                                                        &
                     & *velocityRadialVirial
                velocityRadialDispersionVirial              (iHost,iSatellite                                    )= &
                     & +velocityRadialDispersionVirial      (iHost,iSatellite                                    )  &
                     & +distributionFunction                                                                        &
                     & *velocityRadialVirial    **2
                velocityTangentialMeanVirial                (iHost,iSatellite                                    )= &
                     & +velocityTangentialMeanVirial        (iHost,iSatellite                                    )  &
                     & +distributionFunction                                                                        &
                     & *velocityTangentialVirial
                velocityTangentialDispersionVirial          (iHost,iSatellite                                    )= &
                     & +velocityTangentialDispersionVirial  (iHost,iSatellite                                    )  &
                     & +distributionFunction                                                                        &
                     & *velocityTangentialVirial**2
             end do
          end do
          ! Normalize the distributions.
          distributionFunctionSum=sum(velocityDistributionOrbits(iHost,iSatellite,:,:))
          velocityWidthBin       =+velocity(2) &
               &                  -velocity(1)
          velocityDistributionOrbits                        (iHost,iSatellite,:,:)=   &
               &       +velocityDistributionOrbits          (iHost,iSatellite,:,:)    &
               &       /distributionFunctionSum                                       &
               &       /velocityWidthBin                                          **2
          velocityRadialDistributionOrbits                  (iHost,iSatellite,:  )=   &
               &       +velocityRadialDistributionOrbits    (iHost,iSatellite,:  )    &
               &       /distributionFunctionSum                                       &
               &       /velocityWidthBin
          velocityTangentialDistributionOrbits              (iHost,iSatellite,  :)=   &
               &       +velocityTangentialDistributionOrbits(iHost,iSatellite,  :)    &
               &       /distributionFunctionSum                                       &
               &       /velocityWidthBin
          velocityRadialMeanVirial                          (iHost,iSatellite    )=   &
               &       +velocityRadialMeanVirial            (iHost,iSatellite    )    &
               &       /distributionFunctionSum
          velocityRadialDispersionVirial                    (iHost,iSatellite    )=   &
               & +sqrt(                                                               &
               &       +velocityRadialDispersionVirial      (iHost,iSatellite    )    &
               &       /distributionFunctionSum                                       &
               &       -velocityRadialMeanVirial            (iHost,iSatellite    )**2 &
               &      )
          velocityTangentialMeanVirial                      (iHost,iSatellite    )=   &
               & +velocityTangentialMeanVirial              (iHost,iSatellite    )    &
               & /distributionFunctionSum
          velocityTangentialDispersionVirial                (iHost,iSatellite    )=   &
               & +sqrt(                                                               &
               &       +velocityTangentialDispersionVirial  (iHost,iSatellite    )    &
               &       /distributionFunctionSum                                       &
               &       -velocityTangentialMeanVirial        (iHost,iSatellite    )**2 &
               &      )
          ! Clean up.
          deallocate(nodeSatellite)
       end do
       deallocate(nodeHost)
    end do
    deallocate(tree)
    ! Write output.
    !$ call hdf5Access%set()
    output=galacticusOutputFile%openGroup('mergingHaloOrbitDistribution')
    call output%writeAttribute(self%cosmologyFunctions_%matterDensityEpochal                (self%time),'densityMean'                         )
    call output%writeDataset  (                         velocity                                       ,'velocity'                            )
    call output%writeDataset  (                         mass                                           ,'mass'                                )
    call output%writeDataset  (                         separation                                     ,'separation'                          )
    call output%writeDataset  (                         density                                        ,'density'                             )
    call output%writeDataset  (                         haloMassFunction                               ,'haloMassFunction'                    )
    call output%writeDataset  (                         rateMerging                                    ,'rateMerging'                         )
    call output%writeDataset  (                         velocityRadialMeanVirial                       ,'velocityRadialMeanVirial'            )
    call output%writeDataset  (                         velocityRadialDispersionVirial                 ,'velocityRadialDispersionVirial'      )
    call output%writeDataset  (                         velocityTangentialMeanVirial                   ,'velocityTangentialMeanVirial'        )
    call output%writeDataset  (                         velocityTangentialDispersionVirial             ,'velocityTangentialDispersionVirial'  )
    call output%writeDataset  (                         velocityRadialDistributionOrbits               ,'velocityRadialDistributionOrbits'    )
    call output%writeDataset  (                         velocityTangentialDistributionOrbits           ,'velocityTangentialDistributionOrbits')
    call output%writeDataset  (                         velocityDistributionOrbits                     ,'velocityDistributionOrbits'          )
    call output%close         (                                                                                                               )
    !$ call hdf5Access%unset()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: merging halo orbit distributions')
    return
  end subroutine mergingHaloOrbitDistributionPerform

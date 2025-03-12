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
Implements a normally-distributed halo environment.
!!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Kind_Numbers              , only : kind_int8
  use :: Linear_Growth             , only : linearGrowthClass
  use :: Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt
  use :: Statistics_Distributions  , only : distributionFunction1DNormal                    , distributionFunction1DPeakBackground
  use :: Tables                    , only : table2DLinLinLin

  !![
  <haloEnvironment name="haloEnvironmentNormal">
   <description>Implements a normally-distributed halo environment.</description>
   <deepCopy>
    <functionClass variables="sphericalCollapseSolver_, distributionOverdensity, distributionOverdensityMassive"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="sphericalCollapseSolver_, distributionOverdensity, distributionOverdensityMassive"/>
   </stateStorable>
  </haloEnvironment>
  !!]
  type, extends(haloEnvironmentClass) :: haloEnvironmentNormal
     !!{
     A normal halo environment class.
     !!}
     private
     class           (cosmologyParametersClass                         ), pointer :: cosmologyParameters_            => null()
     class           (cosmologyFunctionsClass                          ), pointer :: cosmologyFunctions_             => null()
     class           (cosmologicalMassVarianceClass                    ), pointer :: cosmologicalMassVariance_       => null()
     class           (linearGrowthClass                                ), pointer :: linearGrowth_                   => null()
     class           (criticalOverdensityClass                         ), pointer :: criticalOverdensity_            => null()
     type            (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt ), pointer :: sphericalCollapseSolver_        => null()
     type            (distributionFunction1DPeakBackground             ), pointer :: distributionOverdensity         => null()
     type            (distributionFunction1DNormal                     ), pointer :: distributionOverdensityMassive  => null()
     type            (table2DLinLinLin                                 )          :: linearToNonLinear
     double precision                                                             :: radiusEnvironment                        , variance           , &
          &                                                                          environmentalOverdensityMaximum          , overdensityPrevious, &
          &                                                                          includedVolumeFraction                   , redshift           , &
          &                                                                          time                                     , massEnvironment
     integer         (kind_int8                                        )          :: uniqueIDPrevious
     logical                                                                      :: linearToNonLinearInitialized
     type            (varying_string                                   )          :: propertyName
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                                  normalDestructor
     procedure :: overdensityLinear             => normalOverdensityLinear
     procedure :: overdensityLinearGradientTime => normalOverdensityLinearGradientTime
     procedure :: overdensityNonLinear          => normalOverdensityNonLinear
     procedure :: environmentRadius             => normalEnvironmentRadius
     procedure :: environmentMass               => normalEnvironmentMass
     procedure :: overdensityLinearMaximum      => normalOverdensityLinearMaximum
     procedure :: pdf                           => normalPDF
     procedure :: cdf                           => normalCDF
     procedure :: overdensityLinearSet          => normalOverdensityLinearSet
     procedure :: volumeFractionOccupied        => normalVolumeFractionOccupied
     procedure :: autoHook                      => normalAutoHook
     procedure :: calculationReset              => normalCalculationReset
  end type haloEnvironmentNormal

  interface haloEnvironmentNormal
     !!{
     Constructors for the {\normalfont \ttfamily normal} halo environment class.
     !!}
     module procedure normalConstructorParameters
     module procedure normalConstructorInternal
  end interface haloEnvironmentNormal

contains

  function normalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily normal} halo environment class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloEnvironmentNormal        )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass            ), pointer       :: linearGrowth_
    class           (criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    double precision                                               :: radiusEnvironment        , redshift, &
         &                                                            massEnvironment          , time

    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <inputParameter>
      <name>massEnvironment</name>
      <source>parameters</source>
      <defaultValue>1.0d15</defaultValue>
      <description>The mass within the sphere sphere used to determine the variance in the environmental density.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusEnvironment</name>
      <source>parameters</source>
      <defaultValue>7.0d0</defaultValue>
      <description>The radius of the sphere used to determine the variance in the environmental density.</description>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift at which the environment is defined.</description>
    </inputParameter>
    !!]
    time=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    !![
    <conditionalCall>
     <call>self=haloEnvironmentNormal(time,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,criticalOverdensity_{conditions})</call>
     <argument name="massEnvironment"   value="massEnvironment"   parameterPresent="parameters"/>
     <argument name="radiusEnvironment" value="radiusEnvironment" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="criticalOverdensity_"     />
    !!]
    return
  end function normalConstructorParameters

  function normalConstructorInternal(time,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,criticalOverdensity_,radiusEnvironment,massEnvironment) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily normal} halo mass function class.
    !!}
    use :: Error_Functions         , only : Error_Function
    use :: Numerical_Constants_Math, only : Pi
    use :: Error                   , only : Error_Report
    use :: ISO_Varying_String      , only : assignment(=)
    implicit none
    type            (haloEnvironmentNormal        )                           :: self
    class           (cosmologyParametersClass     ), target   , intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), target   , intent(in   ) :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass), target   , intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass            ), target   , intent(in   ) :: linearGrowth_
    class           (criticalOverdensityClass     ), target   , intent(in   ) :: criticalOverdensity_
    double precision                                          , intent(in   ) :: time
    double precision                               , optional , intent(in   ) :: radiusEnvironment               , massEnvironment
    double precision                                          , parameter     :: overdensityMean          =0.0d+0
    double precision                                          , parameter     :: limitUpperBuffer         =1.0d-4
    double precision                                                          :: overdensityVariance
    !![
    <constructorAssign variables="time, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_, *linearGrowth_, *criticalOverdensity_, radiusEnvironment, massEnvironment" />
    !!]

    ! Compute environmental radius and mass.
    if (present(radiusEnvironment).and.present(massEnvironment)) call Error_Report('only one of radiusEnvironment and massEnvironment may be specified'//{introspection:location})
    if (present(radiusEnvironment)) then
       self%massEnvironment=+4.0d0                                                                   &
            &               *Pi                                                                      &
            &               *self%radiusEnvironment                                              **3 &
            &               *self%cosmologyFunctions_%matterDensityEpochal(expansionFactor=1.0d0)    &
            &               /3.0d0
    else if (present(massEnvironment)) then
       self%radiusEnvironment=+(                                                                      &
            &                   +3.0d0                                                                &
            &                   *self%massEnvironment                                                 &
            &                   /self%cosmologyFunctions_%matterDensityEpochal(expansionFactor=1.0d0) &
            &                   /4.0d0                                                                &
            &                   /Pi                                                                   &
            &                  )**(1.0d0/3.0d0)
    else
       call Error_Report('one of radiusEnvironment and massEnvironment must be specified'//{introspection:location})
    end if
    ! Set the redshift.
    self%redshift                       =self%cosmologyFunctions_      %redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%time))
    ! Find the root-variance in the linear density field on the given scale.
    self%variance                       =self%cosmologicalMassVariance_%rootVariance               (time=time,mass=self%massEnvironment)**2
    ! Build the distribution function.
    overdensityVariance                 =+self%variance                                                                                     &
         &                               *self%linearGrowth_           %value                      (time=time                          )**2
    ! Construct the distribution for δ. This assumes a normal distribution for the densities, but conditioned on the fact that the
    ! region has not collapsed on any larger scale. The resulting distribution is given by eqn. (9) of Mo & White (1996; MNRAS;
    ! 282; 347). We include some small buffer to the collapse threshold to avoid rounding errors. 
    self%environmentalOverdensityMaximum=+self%criticalOverdensity_    %value                      (time=time,mass=self%massEnvironment)    &
         &                               *(                                                                                                 &
         &                                 +1.0d0                                                                                           &
         &                                 -limitUpperBuffer                                                                                &
         &                                )
    allocate(self%distributionOverdensity)
    !![
    <referenceConstruct owner="self" isResult="yes" object="distributionOverdensity" constructor="distributionFunction1DPeakBackground(overdensityVariance,self%environmentalOverdensityMaximum)"/>
    !!]
    ! Construct a standard normal distribution function which will be used for assigning the overdensity for trees which exceed
    ! the mass of the background.
    allocate(self%distributionOverdensityMassive)
    !![
    <referenceConstruct owner="self" isResult="yes" object="distributionOverdensityMassive" constructor="distributionFunction1DNormal(mean=0.0d0,variance=1.0d0)"/>
    !!]
    ! Find the fraction of cosmological volume which is included in regions below the collapse threshold. This is used to scale
    ! the PDF such that when the mass function is averaged over the PDF we get the correct mass function.
    self%includedVolumeFraction         =Error_Function(self%environmentalOverdensityMaximum/sqrt(2.0d0)/sqrt(overdensityVariance))
    ! Initialize optimizer.
    self%uniqueIDPrevious=-1_kind_int8
    ! Set initialization states.
    self%linearToNonLinearInitialized=.false.
    ! Set name of property to use for environment.
    self%propertyName='haloEnvironmentOverdensity'
    ! Construct a spherical collapse solver.
    allocate(self%sphericalCollapseSolver_)
    !![
    <referenceConstruct owner="self" isResult="yes" object="sphericalCollapseSolver_" constructor="sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt(self%cosmologyFunctions_,self%linearGrowth_)"/>
    !!]
    return
  end function normalConstructorInternal

  subroutine normalAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self

    call calculationResetEvent%attach(self,normalCalculationReset,openMPThreadBindingAllLevels,label='haloEnvironmentNormal')
    return
  end subroutine normalAutoHook

  subroutine normalDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily normal} halo mass function class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(haloEnvironmentNormal), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"          />
    <objectDestructor name="self%cosmologicalMassVariance_"     />
    <objectDestructor name="self%cosmologyFunctions_"           />
    <objectDestructor name="self%criticalOverdensity_"          />
    <objectDestructor name="self%linearGrowth_"                 />
    <objectDestructor name="self%sphericalCollapseSolver_"      />
    <objectDestructor name="self%distributionOverdensity"       />
    <objectDestructor name="self%distributionOverdensityMassive"/>
    !!]
    if (calculationResetEvent%isAttached(self,normalCalculationReset)) call calculationResetEvent%detach(self,normalCalculationReset)
    return
  end subroutine normalDestructor

  subroutine normalCalculationReset(self,node,uniqueID)
    !!{
    Reset the normal halo environment calculation.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    class  (haloEnvironmentNormal), intent(inout) :: self
    type   (treeNode             ), intent(inout) :: node
    integer(kind_int8            ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node, uniqueID

    self%overdensityPrevious=-huge(0.0d0)
    self%uniqueIDPrevious   =-1_kind_int8
    return
  end subroutine normalCalculationReset

  double precision function normalOverdensityLinear(self,node,presentDay)
    !!{
    Return the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (haloEnvironmentNormal               ), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    logical                                               , intent(in   ), optional :: presentDay
    class           (nodeComponentBasic                  ), pointer                 :: basic
    double precision                                                                :: variance
    !![
    <optionalArgument name="presentDay" defaultsTo=".false." />
    !!]

    if (node%hostTree%nodeBase%uniqueID() /= self%uniqueIDPrevious) then
       self%uniqueIDPrevious=node%hostTree%nodeBase%uniqueID()
       if (node%hostTree%properties%exists(self%propertyName)) then
          self%overdensityPrevious=node%hostTree%properties%value(self%propertyName)
       else
          ! Choose an overdensity.
          basic    => node%hostTree%nodeBase        %basic       (                      )
          variance =  self%cosmologicalMassVariance_%rootVariance(basic%mass(),self%time)**2
          if (variance > self%variance) then
             ! The variance on the mass scale of the tree exceeds that of the environment. Therefore, the overdensity is
             ! drawn from the distribution expected for the background scale given that it has not collapsed to become a halo on
             ! any larger scale.
             self%overdensityPrevious=+self%distributionOverdensity       %sample   (                                                                &
                  &                                                                  randomNumberGenerator_=node %hostTree%randomNumberGenerator_    &
                  &                                                                 )
          else
             ! The variance on the mass scale of the tree is less than that of the background. Given that the base halo of the tree
             ! collapsed into a halo of mass Mₜ>Mₑ, the distribution of overdensities on the scale of Mₑ is just a Gaussian, with
             ! mean of δ_c, and variance equal to the difference in variance between the two scales.
             self%overdensityPrevious=+self%distributionOverdensityMassive%sample   (                                                                &
                  &                                                                  randomNumberGenerator_=node %hostTree%randomNumberGenerator_    &
                  &                                                                 )                                                                &
                  &                   *sqrt(                                                                                                         &
                  &                         +self                          %variance                                                                 &
                  &                         -                               variance                                                                 &
                  &                        )                                                                                                         &
                  &                   +      self%criticalOverdensity_      %value  (time                 =basic         %time                 ())
          end if
          call node%hostTree%properties%set(self%propertyName,self%overdensityPrevious)
       end if
    end if
    normalOverdensityLinear=self%overdensityPrevious
    if (.not.presentDay_) then
       basic                   =>  node                                 %basic(                                               )
       normalOverdensityLinear =  +normalOverdensityLinear                                                                      &
            &                     *self                   %linearGrowth_%value(time=basic%time                         (     )) &
            &                     /self                   %linearGrowth_%value(time=self %time                                )
    else
       normalOverdensityLinear =  +normalOverdensityLinear                                                                      &
            &                     *self                   %linearGrowth_%value(time=self%cosmologyFunctions_%cosmicTime(1.0d0)) &
            &                     /self                   %linearGrowth_%value(time=self%time                                 )
    end if
    return
  end function normalOverdensityLinear

  double precision function normalOverdensityLinearGradientTime(self,node)
    !!{
    Return the time gradient of the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentBasic   ), pointer       :: basic

    basic                               =>  node%basic()
    normalOverdensityLinearGradientTime =  +self%overdensityLinear(node)                                                      &
         &                                 *self%linearGrowth_      %logarithmicDerivativeExpansionFactor( time=basic%time()) &
         &                                 *self%cosmologyFunctions_%expansionRate                       (                    &
         &                                  self%cosmologyFunctions_%expansionFactor                      (     basic%time()) &
         &                                                                                               )
    return
  end function normalOverdensityLinearGradientTime

  double precision function normalOverdensityNonLinear(self,node)
    !!{
    Return the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentBasic   ), pointer       :: basic

    ! Get a table of linear vs. nonlinear density.
    if (.not.self%linearToNonLinearInitialized) then
       call self%sphericalCollapseSolver_%linearNonlinearMap(self%time,self%linearToNonLinear)
       self%linearToNonLinearInitialized=.true.
    end if
    ! Find the nonlinear overdensity.
    basic                      =>            node                  %basic      (                                         )
    normalOverdensityNonLinear =  max(-1.0d0,self%linearToNonLinear%interpolate(self%overdensityLinear(node),basic%time()))
    return
  end function normalOverdensityNonLinear

  double precision function normalEnvironmentRadius(self)
    !!{
    Return the radius of the environment.
    !!}
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self

    normalEnvironmentRadius=self%radiusEnvironment
    return
  end function normalEnvironmentRadius

  double precision function normalEnvironmentMass(self)
    !!{
    Return the mass of the environment.
    !!}
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self

    normalEnvironmentMass=self%massEnvironment
    return
  end function normalEnvironmentMass

  double precision function normalOverdensityLinearMaximum(self)
    !!{
    Return the maximum overdensity for which the \gls{pdf} is non-zero.
    !!}
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self

    normalOverdensityLinearMaximum=self%environmentalOverdensityMaximum
    return
  end function normalOverdensityLinearMaximum

  double precision function normalPDF(self,overdensity)
    !!{
    Return the PDF of the environmental overdensity.
    !!}
    implicit none
    class           (haloEnvironmentNormal), intent(inout) :: self
    double precision                       , intent(in   ) :: overdensity

    normalPDF=+self%distributionOverdensity%density(overdensity)
    return
  end function normalPDF

  double precision function normalCDF(self,overdensity)
    !!{
    Return the CDF of the environmental overdensity.
    !!}
    implicit none
    class           (haloEnvironmentNormal), intent(inout) :: self
    double precision                       , intent(in   ) :: overdensity

    normalCDF=self%distributionOverdensity%cumulative(overdensity)
    return
  end function normalCDF

  subroutine normalOverdensityLinearSet(self,node,overdensity)
    !!{
    Set the environmental linear overdensity in the given {\normalfont \ttfamily node}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (haloEnvironmentNormal), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: overdensity

    if (overdensity > self%environmentalOverdensityMaximum) call Error_Report('δ≥δ_c is inconsistent with normal (peak-background) density field'//{introspection:location})
    call node%hostTree%properties%set(self%propertyName,overdensity)
    self%uniqueIDPrevious   =-1_kind_int8
    self%overdensityPrevious=overdensity
    return
  end subroutine normalOverdensityLinearSet

  double precision function normalVolumeFractionOccupied(self)
    !!{
    Return the fraction of the volume occupied by regions described by this environment.
    !!}
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self

    normalVolumeFractionOccupied=self%includedVolumeFraction
    return
  end function normalVolumeFractionOccupied

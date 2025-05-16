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
Implements a log-normal halo environment.
!!}

  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Linear_Growth           , only : linearGrowthClass
  use :: Statistics_Distributions, only : distributionFunction1DLogNormal

  !![
  <haloEnvironment name="haloEnvironmentLogNormal">
   <description>Implements a log-normal halo environment.</description>
  </haloEnvironment>
  !!]
  type, extends(haloEnvironmentClass) :: haloEnvironmentLogNormal
     !!{
     A logNormal halo environment class.
     !!}
     private
     class           (cosmologyParametersClass       ), pointer :: cosmologyParameters_        => null()
     class           (cosmologyFunctionsClass        ), pointer :: cosmologyFunctions_         => null()
     class           (cosmologicalMassVarianceClass  ), pointer :: cosmologicalMassVariance_   => null()
     class           (linearGrowthClass              ), pointer :: linearGrowth_               => null()
     class           (criticalOverdensityClass       ), pointer :: criticalOverdensity_        => null()
     type            (distributionFunction1DLogNormal), pointer :: distributionDensityContrast => null()
     double precision                                           :: radiusEnvironment                    , variance
   contains
     final     ::                                  logNormalDestructor
     procedure :: overdensityLinear             => logNormalOverdensityLinear
     procedure :: overdensityLinearGradientTime => logNormalOverdensityLinearGradientTime
     procedure :: overdensityNonLinear          => logNormalOverdensityNonLinear
     procedure :: environmentRadius             => logNormalEnvironmentRadius
     procedure :: environmentMass               => logNormalEnvironmentMass
     procedure :: overdensityLinearMinimum      => logNormalOverdensityLinearMinimum
     procedure :: pdf                           => logNormalPDF
     procedure :: cdf                           => logNormalCDF
     procedure :: overdensityLinearSet          => logNormalOverdensityLinearSet
  end type haloEnvironmentLogNormal

  interface haloEnvironmentLogNormal
     !!{
     Constructors for the \refClass{haloEnvironmentLogNormal} halo environment class.
     !!}
     module procedure logNormalConstructorParameters
     module procedure logNormalConstructorInternal
  end interface haloEnvironmentLogNormal

contains

  function logNormalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloEnvironmentLogNormal} halo environment class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloEnvironmentLogNormal     )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass            ), pointer       :: linearGrowth_
    class           (criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    double precision                                               :: radiusEnvironment

    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <inputParameter>
      <name>radiusEnvironment</name>
      <source>parameters</source>
      <variable>radiusEnvironment</variable>
      <defaultValue>7.0d0</defaultValue>
      <description>The radius of the sphere used to determine the variance in the environmental density.</description>
    </inputParameter>
    !!]
    self=haloEnvironmentLogNormal(radiusEnvironment,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,criticalOverdensity_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="criticalOverdensity_"     />
    !!]
    return
  end function logNormalConstructorParameters

  function logNormalConstructorInternal(radiusEnvironment,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,criticalOverdensity_) result(self)
    !!{
    Internal constructor for the \refClass{haloEnvironmentLogNormal} halo mass function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (haloEnvironmentLogNormal     )                           :: self
    class           (cosmologyParametersClass     ), target   , intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), target   , intent(in   ) :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass), target   , intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass            ), target   , intent(in   ) :: linearGrowth_
    class           (criticalOverdensityClass     ), target   , intent(in   ) :: criticalOverdensity_
    double precision                                          , intent(in   ) :: radiusEnvironment
    double precision                               , parameter                :: densityContrastMean        =1.0d0
    !![
    <constructorAssign variables="radiusEnvironment, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_, *linearGrowth_, *criticalOverdensity_" />
    !!]

    self%variance                   =self%cosmologicalMassVariance_%rootVariance(                                                        &
         &                                                                       +4.0d0                                                  &
         &                                                                       /3.0d0                                                  &
         &                                                                       *Pi                                                     &
         &                                                                       *self%cosmologyParameters_%OmegaMatter      (     )     &
         &                                                                       *self%cosmologyParameters_%densityCritical  (     )     &
         &                                                                       *self                     %radiusEnvironment       **3, &
         &                                                                        self%cosmologyFunctions_ %cosmicTime       (1.0d0)     &
         &                                                                      )                                                   **2
    ! Construct a log-normal distribution, for 1+δ.
    allocate(self%distributionDensityContrast)
    !![
    <referenceConstruct owner="self" isResult="yes" object="distributionDensityContrast">
     <constructor>
      distributionFunction1DLogNormal(                                                                    &amp;
        &amp;                         mean      =+densityContrastMean                                   , &amp;
        &amp;                         variance  =+self%variance                                         , &amp;
        &amp;                         limitUpper=+1.0d0                                                   &amp;
        &amp;                                    +self%criticalOverdensity_%value(expansionFactor=1.0d0)  &amp;
        &amp;                        )
     </constructor>
    </referenceConstruct>
    !!]
    return
  end function logNormalConstructorInternal

  subroutine logNormalDestructor(self)
    !!{
    Destructor for the \refClass{haloEnvironmentLogNormal} halo mass function class.
    !!}
    implicit none
    type(haloEnvironmentLogNormal), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"      />
    <objectDestructor name="self%cosmologicalMassVariance_" />
    <objectDestructor name="self%linearGrowth_"             />
    <objectDestructor name="self%criticalOverdensity_"      />
    <objectDestructor name="self%cosmologyFunctions_"       />
    !!]
    return
  end subroutine logNormalDestructor

  double precision function logNormalOverdensityLinear(self,node,presentDay)
    !!{
    Return the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Kind_Numbers    , only : kind_int8
    implicit none
    class           (haloEnvironmentLogNormal       ), intent(inout)           :: self
    type            (treeNode                       ), intent(inout)           :: node
    logical                                          , intent(in   ), optional :: presentDay
    class           (nodeComponentBasic             ), pointer                 :: basic
    integer         (kind_int8                      ), save                    :: uniqueIDPrevious           =-1_kind_int8
    double precision                                 , save                    :: overdensityPrevious
    !$omp threadprivate(uniqueIDPrevious,overdensityPrevious)
    !![
    <optionalArgument name="presentDay" defaultsTo=".false." />
    !!]

    if (node%hostTree%nodeBase%uniqueID() /= uniqueIDPrevious) then
       uniqueIDPrevious=node%hostTree%nodeBase%uniqueID()
       if (node%hostTree%properties%exists('haloEnvironmentOverdensity')) then
          overdensityPrevious=node%hostTree%properties%value('haloEnvironmentOverdensity')
       else
           ! Choose an overdensity.
          overdensityPrevious=+self%distributionDensityContrast%sample(                                                            &
               &                                                       randomNumberGenerator_=node%hostTree%randomNumberGenerator_ &
               &                                                      )                                                            &
               &              -1.0d0
          call node%hostTree%properties%set('haloEnvironmentOverdensity',overdensityPrevious)
       end if
    end if
    logNormalOverdensityLinear=overdensityPrevious
    if (.not.presentDay_) then
      basic                      =>  node                                    %basic(                 )
      logNormalOverdensityLinear =  +logNormalOverdensityLinear                                        &
           &                        *self                      %linearGrowth_%value(time=basic%time())
    end if
    return
  end function logNormalOverdensityLinear

  double precision function logNormalOverdensityLinearGradientTime(self,node)
    !!{
    Return the time gradient of the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(haloEnvironmentLogNormal), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node
    class(nodeComponentBasic      ), pointer       :: basic

    basic                                  =>  node%basic()
    logNormalOverdensityLinearGradientTime =  +self%overdensityLinear(node)                                                      &
         &                                    *self%linearGrowth_      %logarithmicDerivativeExpansionFactor( time=basic%time()) &
         &                                    *self%cosmologyFunctions_%expansionRate                       (                    &
         &                                     self%cosmologyFunctions_%expansionFactor                      (     basic%time()) &
         &                                                                                                  )
    return
  end function logNormalOverdensityLinearGradientTime

  double precision function logNormalOverdensityNonLinear(self,node)
    !!{
    Return the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(haloEnvironmentLogNormal), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node

    ! In this model the nonlinear and linear density fields are identified.
    logNormalOverdensityNonLinear=self%overdensityLinear(node)
    return
  end function logNormalOverdensityNonLinear

  double precision function logNormalEnvironmentRadius(self)
    !!{
    Return the radius of the environment.
    !!}
    implicit none
    class(haloEnvironmentLogNormal), intent(inout) :: self

    logNormalEnvironmentRadius=self%radiusEnvironment
    return
  end function logNormalEnvironmentRadius

  double precision function logNormalEnvironmentMass(self)
    !!{
    Return the mass of the environment.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(haloEnvironmentLogNormal), intent(inout) :: self

    logNormalEnvironmentMass=+4.0d0                                                                   &
         &                   *Pi                                                                      &
         &                   *self%radiusEnvironment                                              **3 &
         &                   *self%cosmologyFunctions_%matterDensityEpochal(expansionFactor=1.0d0)    &
         &                   /3.0d0
    return
  end function logNormalEnvironmentMass

  double precision function logNormalOverdensityLinearMinimum(self)
    !!{
    Return the minimum overdensity for which the \gls{pdf} is non-zero.
    !!}
    implicit none
    class(haloEnvironmentLogNormal), intent(inout) :: self
    !$GLC attributes unused :: self

    logNormalOverdensityLinearMinimum=-1.0d0
    return
  end function logNormalOverdensityLinearMinimum

  double precision function logNormalPDF(self,overdensity)
    !!{
    Return the PDF of the environmental overdensity.
    !!}
    implicit none
    class           (haloEnvironmentLogNormal), intent(inout) :: self
    double precision                          , intent(in   ) :: overdensity

    logNormalPDF=self%distributionDensityContrast%density(1.0d0+overdensity)
    return
  end function logNormalPDF

  double precision function logNormalCDF(self,overdensity)
    !!{
    Return the CDF of the environmental overdensity.
    !!}
    implicit none
    class           (haloEnvironmentLogNormal), intent(inout) :: self
    double precision                          , intent(in   ) :: overdensity

    logNormalCDF=self%distributionDensityContrast%cumulative(1.0d0+overdensity)
    return
  end function logNormalCDF

  subroutine logNormalOverdensityLinearSet(self,node,overdensity)
    !!{
    Return the CDF of the environmental overdensity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (haloEnvironmentLogNormal), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    double precision                          , intent(in   ) :: overdensity
    !$GLC attributes unused :: self

    if (overdensity <= -1.0d0) call Error_Report('δ≤-1 is inconsistent with log-normal density field'//{introspection:location})
    call node%hostTree%properties%set('haloEnvironmentOverdensity',overdensity)
    return
  end subroutine logNormalOverdensityLinearSet


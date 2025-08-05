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
  Implementation of a merger tree masses class which samples masses from a distribution.
  !!}

  use :: Merger_Trees_Build_Masses_Distributions, only : mergerTreeBuildMassDistributionClass

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesSampledDistribution" abstract="yes">
   <description>A merger tree masses class which samples masses from a distribution.</description>
  </mergerTreeBuildMasses>
  !!]
  type, extends(mergerTreeBuildMassesClass) :: mergerTreeBuildMassesSampledDistribution
     !!{
     Implementation of a merger tree masses class which samples masses from a distribution.
     !!}
     private
     class           (mergerTreeBuildMassDistributionClass), pointer :: mergerTreeBuildMassDistribution_ => null()
     double precision                                                :: massTreeMinimum                           , massTreeMaximum, &
          &                                                             treesPerDecade
   contains
     !![
     <methods>
       <method description="Handles construction of the abstract parent class." method="construct" />
       <method description="Return a set of values {\normalfont \ttfamily sampleCount} in the interval 0--1, corresponding to values of the cumulative mass distribution." method="sampleCMF" />
     </methods>
     !!]
     final     ::              sampledDistributionDestructor
     procedure :: construct => sampledDistributionConstruct
     procedure :: sampleCMF => sampledDistributionCMF
  end type mergerTreeBuildMassesSampledDistribution

  interface mergerTreeBuildMassesSampledDistribution
     module procedure sampledDistributionConstructorParameters
     module procedure sampledDistributionConstructorInternal
  end interface mergerTreeBuildMassesSampledDistribution

contains

  function sampledDistributionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassesSampledDistribution} merger tree masses class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildMassesSampledDistribution)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (mergerTreeBuildMassDistributionClass    ), pointer       :: mergerTreeBuildMassDistribution_
    double precision                                                          :: massTreeMinimum                 , massTreeMaximum, &
         &                                                                       treesPerDecade
    
    !![
    <inputParameter>
      <name>massTreeMinimum</name>
      <defaultValue>1.0d10</defaultValue>
      <description>The minimum mass of merger tree base halos to consider when sampled masses from a distribution, in units of $M_\odot$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massTreeMaximum</name>
      <defaultValue>1.0d15</defaultValue>
      <description>The maximum mass of merger tree base halos to consider when sampled masses from a distribution, in units of $M_\odot$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>treesPerDecade</name>
      <defaultValue>10.0d0</defaultValue>
      <description>The number of merger trees masses to sample per decade of base halo mass.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerTreeBuildMassDistribution" name="mergerTreeBuildMassDistribution_" source="parameters"/>
    !!]
    self=mergerTreeBuildMassesSampledDistribution(massTreeMinimum,massTreeMaximum,treesPerDecade,mergerTreeBuildMassDistribution_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerTreeBuildMassDistribution_"/>
    !!]
    return
  end function sampledDistributionConstructorParameters

  function sampledDistributionConstructorInternal(massTreeMinimum,massTreeMaximum,treesPerDecade,mergerTreeBuildMassDistribution_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassesSampledDistribution} merger tree masses class.
    !!}
    use :: Display, only : displayMessage, verbosityLevelWarn
    use :: Error  , only : Error_Report
    implicit none
    type            (mergerTreeBuildMassesSampledDistribution)                        :: self
    class           (mergerTreeBuildMassDistributionClass    ), intent(in   ), target :: mergerTreeBuildMassDistribution_
    double precision                                          , intent(in   )         :: massTreeMinimum                 , massTreeMaximum, &
         &                                                                               treesPerDecade
    !![
    <constructorAssign variables="massTreeMinimum, massTreeMaximum, treesPerDecade, *mergerTreeBuildMassDistribution_"/>
    !!]
    
    if (self%massTreeMaximum >= 1.0d16              )                                               &
         & call displayMessage(                                                                     &
         &                     '[massHaloMaximum] > 10¹⁶ - this seems very large and may lead '//   &
         &                     'to failures in merger tree construction'                         ,  &
         &                     verbosityLevelWarn                                                   &
         &                    )
    if (self%massTreeMaximum <= self%massTreeMinimum)                                               &
         & call Error_Report  (                                                                     &
         &                     '[massHaloMaximum] > [massHaloMinimum] is required'             //   &
         &                     {introspection:location}                                             &
         &                    )
    return
  end function sampledDistributionConstructorInternal

  subroutine sampledDistributionDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildMassesSampledDistribution} merger tree masses class.
    !!}
    implicit none
    type(mergerTreeBuildMassesSampledDistribution), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBuildMassDistribution_"/>
    !!]
    return
  end subroutine sampledDistributionDestructor

  subroutine sampledDistributionConstruct(self,time,mass,massMinimum,massMaximum,weight)
    !!{
    Construct a set of merger tree masses by sampling from a distribution.
    !!}
    use            :: Error                  , only : Error_Report
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Numerical_Integration  , only : integrator
    use            :: Numerical_Interpolation, only : interpolator
    use            :: Numerical_Ranges       , only : Make_Range          , rangeTypeLinear
    use            :: Sorting                , only : sort
    use            :: Table_Labels           , only : extrapolationTypeFix
    implicit none
    class           (mergerTreeBuildMassesSampledDistribution), intent(inout)                            :: self
    double precision                                          , intent(in   )                            :: time
    double precision                                          , intent(  out), allocatable, dimension(:) :: mass                                 , weight                                   , &
         &                                                                                                  massMinimum                          , massMaximum
    double precision                                          , parameter                                :: toleranceAbsolute            =1.0d-12, toleranceRelative                 =1.0d-3
    integer                                                   , parameter                                :: massFunctionSamplePerDecade  =100
    double precision                                                         , allocatable, dimension(:) :: massFunctionSampleLogMass            , massFunctionSampleLogMassMonotonic       , &
         &                                                                                                  massFunctionSampleProbability
    integer         (c_size_t                                )                                           :: treeCount                            , iTree
    integer                                                                                              :: iSample                              , jSample                                  , &
         &                                                                                                  massFunctionSampleCount
    double precision                                                                                     :: probability                          , massFunctionSampleLogPrevious
    type            (integrator                              )                                           :: integrator_
    type            (interpolator                            )                                           :: interpolator_
    logical                                                                                              :: monotonize
    !$GLC attributes unused :: weight

    ! Generate a randomly sampled set of halo masses.
    treeCount=max(2_c_size_t,int(log10(self%massTreeMaximum/self%massTreeMinimum)*self%treesPerDecade,kind=c_size_t))
    allocate(mass       (treeCount))
    allocate(massMinimum(treeCount))
    allocate(massMaximum(treeCount))
    ! Create a cumulative probability for sampling halo masses.
    massFunctionSampleCount=max(2,int(log10(self%massTreeMaximum/self%massTreeMinimum)*massFunctionSamplePerDecade))
    allocate(massFunctionSampleLogMass         (massFunctionSampleCount))
    allocate(massFunctionSampleLogMassMonotonic(massFunctionSampleCount))
    allocate(massFunctionSampleProbability     (massFunctionSampleCount))
    massFunctionSampleLogMass    =Make_Range(log10(self%massTreeMinimum),log10(self%massTreeMaximum),massFunctionSampleCount,rangeType=rangeTypeLinear)
    massFunctionSampleLogPrevious=           log10(self%massTreeMinimum)
    jSample=0
    integrator_=integrator(distributionIntegrand,toleranceAbsolute=toleranceAbsolute,toleranceRelative=toleranceRelative)
    do iSample=1,massFunctionSampleCount
       if (massFunctionSampleLogMass(iSample) > massFunctionSampleLogPrevious) then
          probability=integrator_%integrate(                                        &
               &                            massFunctionSampleLogPrevious         , &
               &                            massFunctionSampleLogMass    (iSample)  &
               &                           )
       else
          probability=0.0d0
       end if
       if (iSample == 1) then
          monotonize=.true.
       else if (jSample > 0) then
          monotonize= massFunctionSampleProbability(jSample)+probability &
               &     >                                                   &
               &      massFunctionSampleProbability(jSample)
       else
          monotonize=.false.
       end if
       if (monotonize) then
          jSample=jSample+1
          massFunctionSampleProbability     (jSample)=probability
          massFunctionSampleLogMassMonotonic(jSample)=massFunctionSampleLogMass(iSample)
          if (jSample > 1)                                                                        &
               & massFunctionSampleProbability(jSample)=+massFunctionSampleProbability(jSample  ) &
               &                                        +massFunctionSampleProbability(jSample-1)
       end if
       massFunctionSampleLogPrevious=massFunctionSampleLogMass(iSample)
    end do
    massFunctionSampleCount=jSample
    if (massFunctionSampleCount < 2) call Error_Report('tabulated mass function sampling density has fewer than 2 non-zero points'//{introspection:location})
    ! Normalize the cumulative probability distribution.
    massFunctionSampleProbability(1:massFunctionSampleCount)=+massFunctionSampleProbability(1:massFunctionSampleCount) &
         &                                                   /massFunctionSampleProbability(  massFunctionSampleCount)
    ! Generate a set of points in the cumulative distribution, sort them, and find the mass ranges which they occupy.
    call self%sampleCMF(mass)
    call sort(mass)
    do iTree=1,treeCount
        if (iTree == 1       ) then
          massMinimum(iTree)=+0.0d0
       else
          massMinimum(iTree)=+0.5d0*(mass(iTree)+mass(iTree-1))
       end if
       if (iTree == treeCount) then
          massMaximum(iTree)=+1.0d0
        else
          massMaximum(iTree)=+0.5d0*(mass(iTree)+mass(iTree+1))
       end if
    end do
    ! Compute the corresponding halo masses by interpolation in the cumulative probability distribution function.
    interpolator_=interpolator(                                                                                 &
            &                                    massFunctionSampleProbability     (1:massFunctionSampleCount), &
            &                                    massFunctionSampleLogMassMonotonic(1:massFunctionSampleCount), &
            &                  extrapolationType=extrapolationTypeFix                                           &
            &                 )
    do iTree=1,treeCount
       mass       (iTree)=interpolator_%interpolate(mass       (iTree ))
       massMinimum(iTree)=interpolator_%interpolate(massMinimum(iTree ))
       massMaximum(iTree)=interpolator_%interpolate(massMaximum(iTree ))
    end do
    mass       =10.0d0**mass
    massMinimum=10.0d0**massMinimum
    massMaximum=10.0d0**massMaximum
    deallocate(massFunctionSampleLogMass         )
    deallocate(massFunctionSampleProbability     )
    deallocate(massFunctionSampleLogMassMonotonic)
    return

  contains

    double precision function distributionIntegrand(logMass)
      !!{
      The integrand over the mass function sampling density function.
      !!}
      implicit none
      double precision, intent(in   ) :: logMass

      distributionIntegrand=self%mergerTreeBuildMassDistribution_%sample(10.0d0**logMass,time,self%massTreeMinimum,self%massTreeMaximum)
      return
    end function distributionIntegrand

  end subroutine sampledDistributionConstruct

  subroutine sampledDistributionCMF(self,x)
    !!{
    Stub function for cumulative mass function.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (mergerTreeBuildMassesSampledDistribution), intent(inout)               :: self
    double precision                                          , intent(  out), dimension(:) :: x
    !$GLC attributes unused :: self, x

    call Error_Report('attempt to call function in abstract type'//{introspection:location})
    return
  end subroutine sampledDistributionCMF


!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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

!+    Contributions to this file made by: Andrew Benson, Ethan Nadler.

!!{
Contains a module which implements a excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012},
but using a midpoint method to perform the integrations \citep{du_substructure_2017}, and with a Brownian bridge constraint.
!!}
  
  use :: Cosmological_Density_Field, only : criticalOverdensityClass
  use :: Linear_Growth             , only : linearGrowthClass 

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahiMidpointBrownianBridge">
    <description>
      An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method
      to perform the integrations \citep{du_substructure_2017}, and with a Brownian bridge constraint.
    </description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingFarahiMidpoint) :: excursionSetFirstCrossingFarahiMidpointBrownianBridge
     !!{
     An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, but using a midpoint method
     to perform the integrations \citep{du_substructure_2017}, and with a Brownian bridge constraint.
     !!}
     private
     class           (excursionSetFirstCrossingClass), pointer :: excursionSetFirstCrossing_     => null()
     class           (linearGrowthClass             ), pointer :: linearGrowth_                  => null()
     class           (criticalOverdensityClass      ), pointer :: criticalOverdensity_           => null()
     double precision                                          :: criticalOverdensityConstrained          , varianceConstrained, &
          &                                                       timeConstrained                         , massConstrained
   contains
     final     ::                     farahiMidpointBrownianBridgeDestructor
     procedure :: rate             => farahiMidpointBrownianBridgeRate
     procedure :: rateNonCrossing  => farahiMidpointBrownianBridgeRateNonCrossing
     procedure :: varianceLimit    => farahiMidpointBrownianBridgeVarianceLimit
     procedure :: varianceResidual => farahiMidpointBrownianBridgeVarianceResidual
     procedure :: offsetEffective  => farahiMidpointBrownianBridgeOffsetEffective
     procedure :: fileWrite        => farahiMidpointBrownianBridgeFileWrite
  end type excursionSetFirstCrossingFarahiMidpointBrownianBridge

  interface excursionSetFirstCrossingFarahiMidpointBrownianBridge
     !!{
     Constructors for the Farahi-midpoint Brownian bridge excursion set barrier class.
     !!}
     module procedure farahiMidpointBrownianBridgeConstructorParameters
     module procedure farahiMidpointBrownianBridgeConstructorInternal
  end interface excursionSetFirstCrossingFarahiMidpointBrownianBridge

contains

  function farahiMidpointBrownianBridgeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the Farahi-midpoint excursion set class first crossing class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (excursionSetFirstCrossingFarahiMidpointBrownianBridge)                :: self
    type            (inputParameters                                      ), intent(inout) :: parameters
    double precision                                                                       :: criticalOverdensityConstrained, varianceConstrained, &
         &                                                                                    timeConstrained               , massConstrained    , &
         &                                                                                    timePresent

    self%excursionSetFirstCrossingFarahiMidpoint=excursionSetFirstCrossingFarahiMidpoint(parameters)
    !![
    <objectBuilder class="linearGrowth"              name="self%linearGrowth_"              source="parameters"/>
    <objectBuilder class="criticalOverdensity"       name="self%criticalOverdensity_"       source="parameters"/>
    <objectBuilder class="excursionSetFirstCrossing" name="self%excursionSetFirstCrossing_" source="parameters"/>
    !!]
    timePresent=self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    if      (parameters%isPresent('criticalOverdensityConstrained')) then
       if     (                                                                                                                                                                   &
            &  .not.parameters%isPresent('varianceConstrained')                                                                                                                   &
            & ) call Error_Report('both "criticalOverdensityConstrained" and "varianceConstrained" must be provided'                                  //{introspection:location})
       if     (                                                                                                                                                                   &
            &       parameters%isPresent('timeConstrained'               )                                                                                                        &
            &  .or.                                                                                                                                                               &
            &       parameters%isPresent('massConstrained'               )                                                                                                        &
            & ) call Error_Report('can not mix "criticalOverdensityConstrained/varianceConstrained" and "timeConstrained/massConstrained" constraints'//{introspection:location})
       !![
       <inputParameter>
         <name>criticalOverdensityConstrained</name>
         <source>parameters</source>
         <description>The critical overdensity at the end of the Brownian bridge.</description>
       </inputParameter>
       <inputParameter>
         <name>varianceConstrained</name>
         <source>parameters</source>
         <description>The variance at the end of the Brownian bridge.</description>
       </inputParameter>
       !!]
       massConstrained=self%cosmologicalMassVariance_%mass          (time               =timePresent                   ,rootVariance=sqrt(varianceConstrained))
       timeConstrained=self%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=criticalOverdensityConstrained,mass        =     massConstrained     )
    else if (parameters%isPresent('timeConstrained               ')) then
       if     (                                                                                                                                                                   &
            &  .not.parameters%isPresent('massConstrained'    )                                                                                                                   &
            & ) call Error_Report('both "timeConstrained" and "massConstrained" must be provided'                                                     //{introspection:location})
       if     (                                                                                                                                                                   &
            &       parameters%isPresent('criticalOverdensityConstrained')                                                                                                        &
            &  .or.                                                                                                                                                               &
            &       parameters%isPresent('varianceConstrained'           )                                                                                                        &
            & ) call Error_Report('can not mix "criticalOverdensityConstrained/varianceConstrained" and "timeConstrained/massConstrained" constraints'//{introspection:location})
       !![
       <inputParameter>
         <name>timeConstrained</name>
         <source>parameters</source>
         <description>The time at the end of the Brownian bridge.</description>
       </inputParameter>
       <inputParameter>
         <name>massConstrained</name>
         <source>parameters</source>
         <description>The halo mass at the end of the Brownian bridge.</description>
       </inputParameter>
       !!]
       criticalOverdensityConstrained=+self%criticalOverdensity_     %value       (time=timeConstrained,mass=massConstrained)    &
            &                         /self%linearGrowth_            %value       (time=timeConstrained                     )
       varianceConstrained           =+self%cosmologicalMassVariance_%rootVariance(time=timePresent    ,mass=massConstrained)**2
    else
       criticalOverdensityConstrained=0.0d0
       varianceConstrained           =0.0d0
       timeConstrained               =0.0d0
       massConstrained               =0.0d0
       call Error_Report('must provide either [criticalOverdensityConstrained] and [varianceConstrained], or [timeConstrained] and [massConstrained]')
    end if
    self%criticalOverdensityConstrained=criticalOverdensityConstrained
    self%varianceConstrained           =varianceConstrained
    self%timeConstrained               =timeConstrained
    self%massConstrained               =massConstrained
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function farahiMidpointBrownianBridgeConstructorParameters

  function farahiMidpointBrownianBridgeConstructorInternal(varianceConstrained,criticalOverdensityConstrained,timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,varianceIsUnlimited,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_,criticalOverdensity_,linearGrowth_,excursionSetFirstCrossing_) result(self)
    !!{
    Internal constructor for the Farahi-midpoint excursion set class first crossing class.
    !!}
    implicit none
    type            (excursionSetFirstCrossingFarahiMidpointBrownianBridge)                        :: self
    double precision                                                       , intent(in   )         :: varianceConstrained             , criticalOverdensityConstrained, &
         &                                                                                            timeStepFractional
    integer                                                                , intent(in   )         :: varianceNumberPerUnitProbability, varianceNumberPerUnit         , &
         &                                                                                            timeNumberPerDecade             , varianceNumberPerDecade
    logical                                                                , intent(in   )         :: varianceIsUnlimited
    type            (varying_string                                       ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass                              ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetBarrierClass                             ), intent(in   ), target :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass                        ), intent(in   ), target :: cosmologicalMassVariance_
    class           (linearGrowthClass                                    ), intent(in   ), target :: linearGrowth_
    class           (criticalOverdensityClass                             ), intent(in   ), target :: criticalOverdensity_
    class           (excursionSetFirstCrossingClass                       ), intent(in   ), target :: excursionSetFirstCrossing_
    double precision                                                                               :: timePresent
    !![
    <constructorAssign variables="varianceConstrained, criticalOverdensityConstrained, *criticalOverdensity_, *linearGrowth_, *excursionSetFirstCrossing_"/>
    !!]
    
    self%excursionSetFirstCrossingFarahiMidpoint=excursionSetFirstCrossingFarahiMidpoint(timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,varianceIsUnlimited,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_)
    ! Find mass and time corresponding to the constraint point.
    timePresent         =self%cosmologyFunctions_      %cosmicTime    (expansionFactor    =1.0d0                                                                          )
    self%massConstrained=self%cosmologicalMassVariance_%mass          (time               =timePresent                        ,rootVariance=sqrt(self%varianceConstrained))
    self%timeConstrained=self%criticalOverdensity_     %timeOfCollapse(criticalOverdensity=self%criticalOverdensityConstrained,mass        =     self%massConstrained     )
    return
  end function farahiMidpointBrownianBridgeConstructorInternal

  subroutine farahiMidpointBrownianBridgeDestructor(self)
    !!{
    Destructor for the piecewise Farahi excursion Brownian bridge set first crossing class.
    !!}
    implicit none
    type(excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self

    !![
    <objectDestructor name="self%linearGrowth_"             />
    <objectDestructor name="self%criticalOverdensity_"      />
    <objectDestructor name="self%excursionSetFirstCrossing_"/>
    !!]
    return
  end subroutine farahiMidpointBrownianBridgeDestructor

  double precision function farahiMidpointBrownianBridgeRate(self,variance,varianceProgenitor,time,node)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    double precision                                                       , intent(in   ) :: variance                      , varianceProgenitor , &
         &                                                                                    time
    type            (treeNode                                             ), intent(inout) :: node
    double precision                                                                       :: rootVarianceConstrained       , varianceConstrained, &
         &                                                                                    criticalOverdensityConstrained

    ! For progenitor variances less than or equal to the original variance, return zero.
    if (varianceProgenitor <= variance) then
       farahiMidpointBrownianBridgeRate=0.0d0
       return
    end if
    ! Determine effective constraint point at this epoch.
    rootVarianceConstrained       =+self%cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained           =+rootVarianceConstrained**2
    criticalOverdensityConstrained=+self%criticalOverdensityConstrained             &
         &                         *sqrt(                                           &
         &                               +varianceConstrained                       &
         &                               /self%varianceConstrained                  &
         &                              )
    ! Determine whether to use the conditioned or unconditioned solutions.
    if (.not.node%isOnMainBranch() .or. self%excursionSetBarrier_%barrier(varianceProgenitor,time,node,rateCompute=.true.) > criticalOverdensityConstrained) then
       ! Node is either not on the main branch, or it is on the main branch, but the time corresponds to a barrier above the
       ! constrained point. In either case we want the unconstrained solution.
       farahiMidpointBrownianBridgeRate=self%excursionSetFirstCrossing_%rate           (variance,varianceProgenitor,time,node)
    else if (varianceProgenitor >= varianceConstrained) then
       ! For progenitor variances in excess of the constrained variance the first crossing rate must be zero.
       farahiMidpointBrownianBridgeRate=0.0d0
    else
       ! Use the constrained solution.
       farahiMidpointBrownianBridgeRate=self                           %rateInterpolate(variance,varianceProgenitor,time,node)
    end if
    return
  end function farahiMidpointBrownianBridgeRate

  double precision function farahiMidpointBrownianBridgeRateNonCrossing(self,variance,time,node)
    !!{
    Return the rate for excursion set non-crossing.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    double precision                                                       , intent(in   ) :: time , variance
    type            (treeNode                                             ), intent(inout) :: node
    double precision                                                                       :: rootVarianceConstrained       , varianceConstrained, &
         &                                                                                    criticalOverdensityConstrained

    ! Determine effective constraint point at this epoch.
    rootVarianceConstrained       =+self%cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained           =+rootVarianceConstrained**2
    criticalOverdensityConstrained=+self%criticalOverdensityConstrained             &
         &                         *sqrt(                                           &
         &                               +varianceConstrained                       &
         &                               /self%varianceConstrained                  &
         &                              )
    ! Determine whether to use the conditioned or unconditioned solutions.
    if (.not.node%isOnMainBranch() .or. self%excursionSetBarrier_%barrier(variance,time,node,rateCompute=.true.) > criticalOverdensityConstrained) then
       ! Node is either not on the main branch, or it is on the main branch, but the time corresponds to a barrier above the
       ! constrained point. In either case we want the unconstrained solution.
       farahiMidpointBrownianBridgeRateNonCrossing=self%excursionSetFirstCrossing_%rateNonCrossing           (variance,time,node)
    else if (variance >= varianceConstrained) then
       ! Fo progenitor variances in excess of the constrained variance the non-crossing rate must be zero.
       farahiMidpointBrownianBridgeRateNonCrossing=0.0d0
    else
       ! Use the constrained solution.
       farahiMidpointBrownianBridgeRateNonCrossing=self                           %rateNonCrossingInterpolate(variance,time,node)
    end if
    return
  end function farahiMidpointBrownianBridgeRateNonCrossing
  
  double precision function farahiMidpointBrownianBridgeVarianceLimit(self,varianceProgenitor)
    !!{
    Return the maximum variance to which to tabulate.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    double precision                                                       , intent(in   ) :: varianceProgenitor

    farahiMidpointBrownianBridgeVarianceLimit=min(                                                                                      &
         &                                        self%excursionSetFirstCrossingFarahiMidpoint%varianceLimit      (varianceProgenitor), &
         &                                        self                                        %varianceConstrained                      &
         &                                       )
    return
  end function farahiMidpointBrownianBridgeVarianceLimit

  function farahiMidpointBrownianBridgeVarianceResidual(self,time,variance0,variance1,variance2,cosmologicalMassVariance_) result(varianceResidual)
    !!{
    Return the residual variance between two points.
    !!}
    use :: Kind_Numbers, only : kind_quad
    implicit none
    real (kind_quad                                            )                :: varianceResidual
    class(excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    real (kind_quad                                            ), intent(in   ) :: variance0                 , variance1           , &
         &                                                                         variance2
    double precision                                            , intent(in   ) :: time
    class           (cosmologicalMassVarianceClass             ), intent(inout) :: cosmologicalMassVariance_
    real (kind_quad                                            )                :: rootVarianceConstrained  , varianceConstrained

    rootVarianceConstrained=+cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained    =+rootVarianceConstrained**2
    varianceResidual       =+(varianceConstrained-variance1-variance0) &
         &                  *(variance1          -variance2          ) &
         &                  /(varianceConstrained          -variance0)
    return
  end function farahiMidpointBrownianBridgeVarianceResidual

  function farahiMidpointBrownianBridgeOffsetEffective(self,time,variance0,variance1,variance2,delta0,delta1,delta2,cosmologicalMassVariance_) result(offsetEffective)
    !!{
    Return the residual variance between two points.
    !!}
    use :: Kind_Numbers, only : kind_quad
    implicit none
    real (kind_quad                                            )                :: offsetEffective
    class(excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout) :: self
    real (kind_quad                                            ), intent(in   ) :: delta0                        , delta1             , &
         &                                                                         delta2                        , variance0          , &
         &                                                                         variance1                     , variance2
    double precision                                            , intent(in   ) :: time
    class           (cosmologicalMassVarianceClass             ), intent(inout) :: cosmologicalMassVariance_
    real (kind_quad                                            )                :: rootVarianceConstrained       , varianceConstrained, &
         &                                                                         criticalOverdensityConstrained

    rootVarianceConstrained       =+cosmologicalMassVariance_%rootVariance(self%massConstrained,time)
    varianceConstrained           =+rootVarianceConstrained**2
    criticalOverdensityConstrained=+self%criticalOverdensityConstrained             &
         &                         *sqrt(                                           &
         &                               +varianceConstrained                       &
         &                               /self%varianceConstrained                  &
         &                              )
    offsetEffective              =+(+delta1                        -   delta2) &
         &                        -(+criticalOverdensityConstrained-   delta0) &
         &                        *(+           variance1          -variance2) &
         &                        /(+           varianceConstrained-variance0)    
    return
  end function farahiMidpointBrownianBridgeOffsetEffective
 
  subroutine farahiMidpointBrownianBridgeFileWrite(self)
    !!{
    Write tabulated data on excursion set first crossing probabilities to file.
    !!}
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class           (excursionSetFirstCrossingFarahiMidpointBrownianBridge), intent(inout)               :: self
    type            (hdf5Object                                           )                              :: dataFile          , dataGroup
    double precision                                                       , allocatable  , dimension(:) :: linearGrowthFactor
    integer                                                                                              :: i
    
    ! Write the primary data.
    call self%excursionSetFirstCrossingFarahiMidpoint%fileWrite()
    ! Don't write anything if neither table is initialized.
    if (.not.self%tableInitializedRate) return
    ! Open the data file.
    !$ call hdf5Access%set()
    call dataFile%openFile(char(self%fileName),overWrite=.false.)
    ! Check if the rate table is populated.
    if (self%tableInitializedRate) then
       allocate(linearGrowthFactor(size(self%timeTableRate)))
       do i=1,size(self%timeTableRate)
          linearGrowthFactor(i)=self%linearGrowth_%value(time=self%timeTableRate(i))
       end do
       dataGroup=dataFile%openGroup("rate")
       call dataGroup%writeDataset(linearGrowthFactor,'linearGrowthFactor','The linear growth factors at the times at which results are tabulated.')
       call dataGroup%close()
    end if
    ! Close the data file.
    call dataFile%close()
    !$ call hdf5Access%unset()
    return
  end subroutine farahiMidpointBrownianBridgeFileWrite

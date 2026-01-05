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
  Implements a transfer function envelope class which computes a transfer function that is a monotonically-decreasing (as a
  function of wavenumber) envelope to another transfer function.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: Tables              , only : table1DMonotoneCSpline

  !![
  <transferFunction name="transferFunctionEnvelope">
    <description>
      A transfer function envelope class which computes a transfer function that is a monotonically-decreasing (as a function of
      wavenumber) envelope to another transfer function. There is no unique monotonic envelope to an arbitrary function. The
      approach taken here largely follows the algorithm described by \cite{mcclain_algorithm_1991}. Briefly, inflection points in
      the transfer function are identified, and the midpoints of the corresponding concave upward intervals are used as initial
      guesses for the tangent points of the envelope and original transfer function. These guesses for the tangent points are
      iteratively updated to remove regions when the envelope function fails (i.e. is below the original function).

      If {\normalfont \ttfamily [transferFunctionReference]} is supplied then half-, quarter-, and fraction-mode masses relative
      to that reference transfer function can be computed using the envelope function. If {\normalfont \ttfamily
      [envelopeModeMassesOnly]} is set to true, then the envelope transfer function is used \emph{only} for calculation of these
      mode masses---the original (non-envelope) transfer function is returned in all other cases. If {\normalfont \ttfamily
      [envelopeModeMassesOnly]} is set to false then the enveloped transfer function is used for \emph{all} calculations.
    </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionEnvelope
     !!{
     A transfer function envelope class which computes a transfer function that is a monotonically-decreasing (as a function of
     wavenumber) envelope to another transfer function.
     !!}
     private
     type            (table1DMonotoneCSpline  )          :: transferTable
     class           (transferFunctionClass   ), pointer :: transferFunction_            => null() , transferFunctionReference    => null()
     double precision                                    :: wavenumberMinimum                      , wavenumberMaximum
     logical                                             :: envelopeModeMassesOnly                 , envelopeRatio
     double precision                                    :: wavenumberMinimumLogarithmic           , wavenumberMaximumLogarithmic           , &
          &                                                 normalization                          , normalizationReference
     integer                                             :: tablePointsPerDecade
     logical                                             :: tableInitialized             =  .false., modeMassSolving              =  .false.
   contains
     final     ::                          envelopeDestructor
     procedure :: value                 => envelopeValue
     procedure :: logarithmicDerivative => envelopeLogarithmicDerivative
     procedure :: halfModeMass          => envelopeHalfModeMass
     procedure :: quarterModeMass       => envelopeQuarterModeMass
     procedure :: fractionModeMass      => envelopeFractionModeMass
     procedure :: epochTime             => envelopeEpochTime
  end type transferFunctionEnvelope

  interface transferFunctionEnvelope
     !!{
     Constructors for the envelope transfer function class.
     !!}
     module procedure envelopeConstructorParameters
     module procedure envelopeConstructorInternal
  end interface transferFunctionEnvelope

contains

  function envelopeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the envelope transfer function class which takes a parameter set as input.
    !!}
    implicit none
    type            (transferFunctionEnvelope)                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (transferFunctionClass   ), pointer       :: transferFunction_     , transferFunctionReference
    integer                                                   :: tablePointsPerDecade
    double precision                                          :: wavenumberMinimum     , wavenumberMaximum
    logical                                                   :: envelopeModeMassesOnly, envelopeRatio

    !![
    <inputParameter>
      <name>envelopeRatio</name>
      <source>parameters</source>
      <defaultValue>parameters%isPresent('transferFunctionReference')</defaultValue>
      <description>If true, the envelope is computed on the ratio of the transfer function to the reference transfer function, otherwise it is computed directly on the transfer function.</description>
    </inputParameter>
    <inputParameter>
      <name>envelopeModeMassesOnly</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, the envelope is used \emph{only} when computing fractional mode masses, and does not affect the transfer function itself.</description>
    </inputParameter>
    <inputParameter>
      <name>tablePointsPerDecade</name>
      <source>parameters</source>
      <defaultValue>100</defaultValue>
      <description>The number of points per decade of wavenumber at which to tabulate the transfer function.</description>
    </inputParameter>
    <inputParameter>
      <name>wavenumberMinimum</name>
      <source>parameters</source>
      <defaultValue>1.0d-6</defaultValue>
      <description>The minimum wavenumber to which the envelope should be computed initially (the envelope range will be extended to small wavenumbers as needed).</description>
    </inputParameter>
    <inputParameter>
      <name>wavenumberMaximum</name>
      <source>parameters</source>
      <defaultValue>1.0d+6</defaultValue>
      <description>The maximum wavenumber to which the envelope should be computed.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunction_"    source="parameters"/>
    !!]
    if (parameters%isPresent('transferFunctionReference')) then
       !![
       <objectBuilder class="transferFunction" name="transferFunctionReference" source="parameters" parameterName="transferFunctionReference"/>
       !!]
    else
       transferFunctionReference => null()
    end if
    !![
    <conditionalCall>
      <call>self=transferFunctionEnvelope(tablePointsPerDecade,wavenumberMinimum,wavenumberMaximum,envelopeRatio,envelopeModeMassesOnly,cosmologyParameters_,transferFunction_{conditions})</call>
      <argument name="transferFunctionReference" value="transferFunctionReference" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="transferFunction_"   />
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    if (parameters%isPresent('transferFunctionReference')) then
       !![
       <objectDestructor name="transferFunctionReference"/>
       !!]
    end if
    return
  end function envelopeConstructorParameters

  function envelopeConstructorInternal(tablePointsPerDecade,wavenumberMinimum,wavenumberMaximum,envelopeRatio,envelopeModeMassesOnly,cosmologyParameters_,transferFunction_,transferFunctionReference) result(self)
    !!{
    Internal constructor for the envelope transfer function class.
    !!}
    implicit none
    type            (transferFunctionEnvelope)                                  :: self
    class           (cosmologyParametersClass), intent(in   ), target           :: cosmologyParameters_
    class           (transferFunctionClass   ), intent(in   ), target           :: transferFunction_
    class           (transferFunctionClass   ), intent(in   ), target, optional :: transferFunctionReference
    integer                                   , intent(in   )                   :: tablePointsPerDecade
    double precision                          , intent(in   )                   :: wavenumberMinimum               , wavenumberMaximum
    logical                                   , intent(in   )                   :: envelopeModeMassesOnly          , envelopeRatio
    double precision                          , parameter                       :: wavenumberTiny           =1.0d-4
    
    !![
    <constructorAssign variables="*cosmologyParameters_, *transferFunction_, *transferFunctionReference, tablePointsPerDecade, wavenumberMinimum, wavenumberMaximum, envelopeRatio, envelopeModeMassesOnly"/>
    !!]
    
    self%modeMassSolving             =.false.
    self%tableInitialized            =.false.
    self%wavenumberMinimumLogarithmic=log(wavenumberMinimum)
    self%wavenumberMaximumLogarithmic=log(wavenumberMaximum)
    self%normalization               =self%transferFunction_        %value(wavenumber=wavenumberTiny)
    if (present(transferFunctionReference)) then
       self%normalizationReference   =self%transferFunctionReference%value(wavenumber=wavenumberTiny)
    else
       self%transferFunctionReference => null()
       if (envelopeRatio) call Error_Report('can not compute envelope on ratio with a reference transfer function'//{introspection:location})
    end if
    return
  end function envelopeConstructorInternal

  subroutine envelopeDestructor(self)
    !!{
    Destructor for the envelope transfer function class.
    !!}
    implicit none
    type(transferFunctionEnvelope), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%transferFunction_"   />
    !!]
    if (associated(self%transferFunctionReference)) then
       !![
       <objectDestructor name="self%transferFunctionReference"/>
       !!]
    end if
    if (self%tableInitialized) call self%transferTable%destroy()
    return
  end subroutine envelopeDestructor

  double precision function envelopeValue(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionEnvelope), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    double precision                                          :: wavenumberLogarithmic

    if (self%envelopeModeMassesOnly .and. .not.self%modeMassSolving) then
       envelopeValue=+self%transferFunction_%value(wavenumber) &
            &        /self%normalization
    else
       wavenumberLogarithmic=log(wavenumber)
       call envelopeTabulate(self,wavenumberLogarithmic)
       envelopeValue=exp(self%transferTable%interpolate(wavenumberLogarithmic))
       if (self%envelopeRatio)                                                        &
            & envelopeValue=+                               envelopeValue             &
            &               *self%transferFunctionReference%        value(wavenumber) &
            &               /self%normalizationReference
    end if
    return
  end function envelopeValue

  double precision function envelopeLogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionEnvelope), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber
    double precision                                          :: wavenumberLogarithmic

    if (self%envelopeModeMassesOnly) then
       envelopeLogarithmicDerivative=self%transferFunction_%logarithmicDerivative(wavenumber)
    else
       wavenumberLogarithmic=log(wavenumber)
       call envelopeTabulate(self,wavenumberLogarithmic)
       envelopeLogarithmicDerivative=+self%transferTable%interpolateGradient(wavenumberLogarithmic)
       if (self%envelopeRatio)                                                                                        &
            & envelopeLogarithmicDerivative=+                               envelopeLogarithmicDerivative             &
            &                               +self%transferFunctionReference%        logarithmicDerivative(wavenumber)
    end if
    return
  end function envelopeLogarithmicDerivative

  double precision function envelopeEpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionEnvelope), intent(inout) :: self

    envelopeEpochTime=self%transferFunction_%epochTime()
    return
  end function envelopeEpochTime

  double precision function envelopeHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of two relative to a \gls{cdm} transfer function
    !!}
    implicit none
    class  (transferFunctionEnvelope), intent(inout), target   :: self
    integer                          , intent(  out), optional :: status
    
    envelopeHalfModeMass=self%fractionModeMass(0.50d0,status)
    return
  end function envelopeHalfModeMass

  double precision function envelopeQuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of four relative to a \gls{cdm} transfer function
    !!}
    implicit none
    class  (transferFunctionEnvelope), intent(inout), target   :: self
    integer                          , intent(  out), optional :: status

    envelopeQuarterModeMass=self%fractionModeMass(0.25d0,status)
    return
  end function envelopeQuarterModeMass

  double precision function envelopeFractionModeMass(self,fraction,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a given fraction relative to a \gls{cdm} transfer function
    !!}
    use :: Error                   , only : Error_Report             , errorStatusSuccess           , errorStatusFail
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (transferFunctionEnvelope), intent(inout), target   :: self
    double precision                          , intent(in   )           :: fraction
    integer                                   , intent(  out), optional :: status
    type            (rootFinder              )                          :: finder
    double precision                                                    :: wavenumberFractionMode, matterDensity

    if (associated(self%transferFunctionReference)) then
       ! There is no analytic solution for the fraction-mode mass so we resort to numerical root finding.
       finder                  = rootFinder(                                                             &
            &                               rootFunction                 =modeSolver                   , &
            &                               toleranceRelative            =1.000d-3                     , &
            &                               rangeExpandUpward            =2.000d+0                     , &
            &                               rangeExpandDownward          =0.500d+0                     , &
            &                               rangeExpandType              =rangeExpandMultiplicative    , &
            &                               rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                               rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
            &                              )
       call envelopeTabulate(self,wavenumberLogarithmic=0.0d0)
       self%modeMassSolving    =.true.
       wavenumberFractionMode  = finder%find(rootGuess=1.0d0)
       matterDensity           =+self%cosmologyParameters_%OmegaMatter    () &
            &                   *self%cosmologyParameters_%densityCritical()
       ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
       ! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].
       envelopeFractionModeMass=+4.0d0                    &
            &                   *Pi                       &
            &                   /3.0d0                    &
            &                   *matterDensity            &
            &                   *(                        &
            &                     +Pi                     &
            &                     /wavenumberFractionMode &
            &                    )**3
       self%modeMassSolving    =.false.
       if (present(status)) status=errorStatusSuccess
    else
       envelopeFractionModeMass=0.0d0
       if (present(status)) then
          status=errorStatusFail
       else
          call Error_Report('fractional mode mass can not be computed as no reference transfer function was supplied'//{introspection:location})
       end if
    end if
    return

  contains
    
    double precision function modeSolver(wavenumber)
      !!{
      Function used in solving for fraction-mode masses in the enveloped transfer function.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      modeSolver=+self                          %value   (wavenumber) &
           &     /self%transferFunctionReference%value   (wavenumber) &
           &     *self%normalizationReference                         &
           &     -                               fraction
      return
    end function modeSolver
    
  end function envelopeFractionModeMass
  
  subroutine envelopeTabulate(self,wavenumberLogarithmic)
    !!{
    Tabulate the envelope to the transfer function.
    !!}
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLinear
    implicit none
    class           (transferFunctionEnvelope), intent(inout)               :: self
    double precision                          , intent(in   )               :: wavenumberLogarithmic
    integer                                   , parameter                   :: iterationMaximum             =20
    double precision                          , dimension(:  ), allocatable :: wavenumberLogarithmicDense      , wavenumberLogarithmicSparse, &
         &                                                                     transferFunctionDense           , transferFunctionSparse     , &
         &                                                                     transferFunctionGradientDense
    integer                                   , dimension(:  ), allocatable :: iPivot                          , iPivotNew                  , &
         &                                                                     iInflection
    integer                                   , dimension(:,:), allocatable :: iFail
    double precision                                                        :: transferFunction
    logical                                                                 :: makeTable                       , failing
    integer                                                                 :: pointCountDense                 , pointCountSparse           , &
         &                                                                     pointCountSparseNew             , iLastNonZero               , &
         &                                                                     countInflections                , countFail                  , &
         &                                                                     i                               , j                          , &
         &                                                                     k                               , iteration

    makeTable=.not.self%tableInitialized
    if (.not.makeTable) &
         & makeTable=wavenumberLogarithmic < self%wavenumberMinimumLogarithmic
    if (wavenumberLogarithmic > self%wavenumberMaximumLogarithmic .and. self%tableInitialized) &
         & call Error_Report('wavenumber exceeds the maximum for which the envelope was computed'//{introspection:location})
    if (makeTable) then
       self%wavenumberMinimumLogarithmic=min(self%wavenumberMinimumLogarithmic,wavenumberLogarithmic-1.0d0)
       self%wavenumberMaximumLogarithmic=max(self%wavenumberMaximumLogarithmic,wavenumberLogarithmic      )
       pointCountDense=int((self%wavenumberMaximumLogarithmic-self%wavenumberMinimumLogarithmic)*dble(self%tablePointsPerDecade)/log(10.0d0))+1
       allocate(wavenumberLogarithmicDense   (pointCountDense  ))
       allocate(transferFunctionDense        (pointCountDense  ))
       allocate(transferFunctionGradientDense(pointCountDense  ))
       allocate(iInflection                  (pointCountDense  ))
       allocate(iFail                        (pointCountDense,2))
       allocate(iPivot                       (pointCountDense  ))
       allocate(iPivotNew                    (pointCountDense  ))
       wavenumberLogarithmicDense=Make_Range(self%wavenumberMinimumLogarithmic,self%wavenumberMaximumLogarithmic,pointCountDense,rangeType=rangeTypeLinear)
       ! Tabulate the transfer function, identifying the last non-zero point. Note that we tabulate the absolute value here - some
       ! oscillatory transfer functions can be negative [since all that is actually used in calculations is T²(k)].
       iLastNonZero=0
       do i=1,pointCountDense
          transferFunctionDense       (i)=+abs(self%transferFunction_        %value(exp(wavenumberLogarithmicDense(i)))/self%normalization         )
          if (self%envelopeRatio) &
               & transferFunctionDense(i)=+transferFunctionDense(i)                                                                                  &
               &                          /abs(self%transferFunctionReference%value(exp(wavenumberLogarithmicDense(i)))/self%normalizationReference)
          if (transferFunctionDense(i) > 0.0d0) iLastNonZero=i
       end do
       ! Compute estimates of the gradient of the transfer function, ∂T (defined as dT(k)/dlogk)), using a symmetric finite-difference approach.
       do i=2,pointCountDense-1
          transferFunctionGradientDense(i)=+(     transferFunctionDense(i+1)-     transferFunctionDense(i-1)) &
               &                           /(wavenumberLogarithmicDense(i+1)-wavenumberLogarithmicDense(i-1))
       end do
       transferFunctionGradientDense(              1)=transferFunctionGradientDense(                2)
       transferFunctionGradientDense(pointCountDense)=transferFunctionGradientDense(pointCountDense-1)
       ! Find inflection points, corresponding to ∂²T=0. These correspond to minima and maxima of the gradient.
       countInflections=0
       do i=2,pointCountDense-1
          if     (                                                                       &
               &   transferFunctionGradientDense(i) > transferFunctionGradientDense(i-1) &
               &  .and.                                                                  &
               &   transferFunctionGradientDense(i) > transferFunctionGradientDense(i+1) &
               & ) then
             ! Inflection point found (maximum of ∂T).
             countInflections                     =countInflections+1
             iInflection        (countInflections)=i
          end if
          if     (                                                                       &
               &   transferFunctionGradientDense(i) < transferFunctionGradientDense(i-1) &
               &  .and.                                                                  &
               &   transferFunctionGradientDense(i) < transferFunctionGradientDense(i+1) &
               & ) then
             ! Inflection point found (minimum of ∂T).
             if (countInflections == 0) then
                ! If no prior inflection point has been found, add one corresponding to the first point in our table. We will
                ! later identify our initial guesses for the tangent points with the midpoints between of concave upward regions -
                ! these corresponding to the midpoint between inflection points corresponding to successive maxima and minima. If
                ! our first point is a minimum, we add this extra point to ensure that a tangent point is introduced before this
                ! first minimum.
                countInflections                  =countInflections+1
                iInflection     (countInflections)=1
             end if
             countInflections                     =countInflections+1
             iInflection        (countInflections)=i
           end if
       end do
       ! Construct set of initial guesses for our pivot points for interpolation. These are set to the midpoints of concave upward
       ! regions of the transfer function, plus the boundary points at the start and end of our tabulation.
       pointCountSparse=countInflections/2+2
       allocate(wavenumberLogarithmicSparse(pointCountSparse))
       allocate(transferFunctionSparse     (pointCountSparse))
       iPivot(               1)=1
       iPivot(pointCountSparse)=iLastNonZero
       do i=1,pointCountSparse-2
          iPivot(i+1)=+(                        &
               &        +iInflection(2*(i-1)+1) &
               &        +iInflection(2*(i-1)+2) &
               &       )                        &
               &      /2
       end do
       wavenumberLogarithmicSparse=wavenumberLogarithmicDense(iPivot(1:pointCountSparse))
       transferFunctionSparse     =     transferFunctionDense(iPivot(1:pointCountSparse))
       ! Remove non-monotonic points, as we want a monotonic envelope function.
       i=pointCountSparse-1
       do while (i > 1)
          if (transferFunctionSparse(i) <= transferFunctionSparse(i+1)) then
             iPivot          (i:pointCountSparse-1)=iPivot          (i+1:pointCountSparse)
             pointCountSparse                      =pointCountSparse                      -1
          end if
          i=i-1
       end do
       ! Build the interpolation table.
       wavenumberLogarithmicSparse=wavenumberLogarithmicDense(iPivot(1:pointCountSparse))
       transferFunctionSparse     =     transferFunctionDense(iPivot(1:pointCountSparse))
       call self%transferTable%create  (    wavenumberLogarithmicSparse(1:pointCountSparse) )
       call self%transferTable%populate(log(     transferFunctionSparse(1:pointCountSparse)))
       ! Begin iteratively improving our envelope function.
       countFail=huge(0)
       iteration=     0
       do while (countFail > 0 .and. iteration < iterationMaximum)
          iteration=iteration+1
          ! Check for regions where the interpolation fails to envelope.
          countFail=0
          failing  =.false.
          do i=2,pointCountDense
             ! Ignore wavenumbers beyond the non-zero range.
             if (wavenumberLogarithmicDense(i) > wavenumberLogarithmicSparse(pointCountSparse)) exit
             ! Check if the interpolated envelope function lies below the original transfer function at this point.
             transferFunction=exp(self%transferTable%interpolate(wavenumberLogarithmicDense(i)))
             if (transferFunction < transferFunctionDense(i).and..not.failing) then
                ! The interpolated envelope function does lie below the original transfer function, and we are not currently in a
                ! failing region. Therefore, begin a new failing region.
                countFail             =countFail+1
                iFail    (countFail,1)=i        -1
                failing               =.true.
             else if (transferFunction >= transferFunctionDense(i).and.failing) then
                ! The interpolated envelope function lies above the original transfer function, and we are currently in a failing
                ! region. Therefore, we have reached the end of the failing region.
                iFail    (countFail,2)=i
                failing               =.false.
             end if
          end do
          ! If the final failing region extends to the end of the tabulated region, end that region now.
          if (failing) iFail(countFail,2)=i
          ! If there are failing regions, perform an update of our pivot points.
          if (countFail > 0) then
             ! First, for each existing pivot point, if it lies within a failing region, move it to the center of that failing
             ! region. This gives an improved estimate of the tangent point.
             pointCountSparseNew=pointCountSparse
             do i=1,pointCountSparse
                iPivotNew(i)=iPivot(i)
                do j=1,countFail
                   if (iPivot(i) >= iFail(j,1) .and. iPivot(i) < iFail(j,2)) then
                      iPivotNew(i)=sum(iFail(j,:))/2
                      exit
                   end if
                end do
             end do
             ! If the first pivot was moved by the above, add a new pivot to anchor the start of the table.
             if (iPivotNew(1) /= 1) then
                do j=pointCountSparse,1,-1
                   iPivotNew(j+1)=iPivotNew(j)
                end do
                iPivotNew          (1)=                    1
                pointCountSparseNew   =pointCountSparseNew+1
             end if
             ! Look for failed regions containing no pivot. We want to add a new pivot point inside these regions - we add this at
             ! the region midpoint.
             do i=1,countFail
                if (.not.any(iPivotNew(1:pointCountSparseNew) >= iFail(i,1) .and. iPivotNew(1:pointCountSparseNew) < iFail(i,2))) then
                   j=1
                   do while (j <= pointCountSparseNew)
                      if (iPivotNew(j) > iFail(i,1)) then
                         do k=pointCountSparseNew,j,-1
                            iPivotNew(k+1)=iPivotNew(k)
                         end do
                         iPivotNew          (j)=sum(iFail(i,:))/2
                         pointCountSparseNew   =pointCountSparseNew+1
                         exit
                      end if
                      j=j+1
                   end do
                end if
             end do
             ! Remove any duplicated pivot points.
             i=1
             do while (i < pointCountSparseNew)
                if (iPivotNew(i+1)==iPivotNew(i)) then
                   do j=i,pointCountSparseNew-1
                      iPivotNew(j)=iPivotNew(j+1)
                   end do
                   pointCountSparseNew=pointCountSparseNew-1
                else
                   i=i+1
                end if
             end do
             ! Construct the updated interpolation.           
             iPivot                     =          iPivotNew
             pointCountSparse           =pointCountSparseNew
             wavenumberLogarithmicSparse=wavenumberLogarithmicDense(iPivot(1:pointCountSparse))
             transferFunctionSparse     =     transferFunctionDense(iPivot(1:pointCountSparse))
             call self%transferTable%destroy (                                )
             call self%transferTable%create  (    wavenumberLogarithmicSparse )
             call self%transferTable%populate(log(     transferFunctionSparse))
          end if
       end do
       ! Clean up.
       deallocate(wavenumberLogarithmicDense   )
       deallocate(transferFunctionDense        )
       deallocate(transferFunctionGradientDense)
       deallocate(iInflection                  )
       deallocate(iFail                        )
       deallocate(iPivot                       )
       deallocate(iPivotNew                    )
       deallocate(wavenumberLogarithmicSparse  )
       deallocate(transferFunctionSparse       )
       ! Mark the table as constructed
       self%tableInitialized=.true.
    end if
    return
  end subroutine envelopeTabulate


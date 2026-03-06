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

  ! Include an explicit dependency on the low-level multivariate distribution integration code.
  !: $(BUILDPATH)/Genz/mvndst.o

  !!{
  Implementation of a multivariate normal distribution function.
  !!}

  use :: Linear_Algebra, only : matrix

  !![
  <distributionFunctionMultivariate name="distributionFunctionMultivariateNormal">
   <description>
    A multivariate normal distribution.
   </description>
  </distributionFunctionMultivariate>
  !!]
  type, extends(distributionFunctionMultivariateClass) :: distributionFunctionMultivariateNormal
     !!{
     Implementation of a normal 1D distribution function.
     !!}
     private
     double precision                                      :: errorAbsolute        , errorRelative
     integer                                               :: countTrialsMaximum
     double precision        , dimension(:  ), allocatable :: mean                 , rootVariance            , &
          &                                                   correlation
     double precision        , dimension(:,:), allocatable :: covariance
     type            (matrix)                              :: covariance_
     double precision                                      :: normalization        , logNormalization
     logical                                               :: normalizationComputed, logNormalizationComputed
   contains
     !![
     <methods>
       <method method="cumulativeMonteCarlo" description="Compute the cumulative distribution function used Monte Carlo methods."/>
     </methods>
     !!]
     procedure :: density              => multivariateNormalDensity
     procedure :: cumulative           => multivariateNormalCumulative
     procedure :: cumulativeMonteCarlo => multivariateNormalCumulativeMonteCarlo
  end type distributionFunctionMultivariateNormal

  interface distributionFunctionMultivariateNormal
     !!{
     Constructors for the {\normalfont \ttfamily normal} 1D distribution function class.
     !!}
     module procedure multivariateNormalConstructorParameters
     module procedure multivariateNormalConstructorInternal
  end interface distributionFunctionMultivariateNormal

  ! Variables used for debugging purposes.
  double precision, dimension(:), allocatable :: xlow__                          , xhigh__
  logical         , dimension(:), allocatable :: useTransform_
  double precision                            :: probabilityLogarithmicReference_
  integer                                     :: signalDirect
  !$omp threadprivate(xLow__,xHigh__,usetransform_,probabilityLogarithmicReference_,signalDirect)
  
contains

  function multivariateNormalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily multivariateNormal} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (distributionFunctionMultivariateNormal)                              :: self
    type            (inputParameters                       ), intent(inout)               :: parameters
    class           (randomNumberGeneratorClass            ), pointer                     :: randomNumberGenerator_
    double precision                                        , dimension(:  ), allocatable :: mean
    double precision                                        , dimension(:,:), allocatable :: covariance
    double precision                                                                      :: errorAbsolute         , errorRelative
    integer                                                                               :: countTrialsMaximum

    !![
    <inputParameter>
      <name>mean</name>
      <description>The mean of the multivariate normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>covariance</name>
      <description>The covariance of the multivariate normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>errorAbsolute</name>
      <description>The absolute error tolerance in determining the cumulative distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>errorRelative</name>
      <description>The relative error tolerance in determining the cumulative distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countTrialsMaximum</name>
      <description>The maximum number of trials allowed in Monte Carlo evaluation of the cumulative distribution function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=distributionFunctionMultivariateNormal(mean,covariance,errorAbsolute,errorRelative,countTrialsMaximum,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function multivariateNormalConstructorParameters

  function multivariateNormalConstructorInternal(mean,covariance,errorAbsolute,errorRelative,countTrialsMaximum,randomNumberGenerator_) result(self)
    !!{
    Constructor for {\normalfont \ttfamily multivariateNormal} multivariate distribution function class.
    !!}
    implicit none
    type            (distributionFunctionMultivariateNormal)                                        :: self
    double precision                                        , intent(in   ), dimension(:  )         :: mean
    double precision                                        , intent(in   ), dimension(:,:)         :: covariance
    class           (randomNumberGeneratorClass            ), intent(in   ), optional      , target :: randomNumberGenerator_
    double precision                                        , intent(in   )                         :: errorAbsolute         , errorRelative
    integer                                                 , intent(in   )                         :: countTrialsMaximum
    integer                                                                                         :: i                     , j
    !![
    <constructorAssign variables="mean, covariance ,errorAbsolute, errorRelative, countTrialsMaximum, *randomNumberGenerator_"/>
    !!]

    if     (                                      &
         &   size(covariance,dim=1) /= size(mean) &
         &  .or.                                  &
         &   size(covariance,dim=2) /= size(mean) &
         & ) call Error_Report('`covariance` shape is invalid'//{introspection:location})
    self%covariance_             =matrix(covariance)
    self%normalizationComputed   =.false.
    self%logNormalizationComputed=.false.
    allocate(self%rootVariance(size(mean)                 ))
    allocate(self%correlation (size(mean)*(size(mean)-1)/2))
    do i=1,size(mean)
       if (self%covariance(i,i) <= 0.0d0) call Error_Report('non-positive variance detected'//{introspection:location})
       self   %rootVariance(i                )=+sqrt(                      &
            &                                        +self%covariance(i,i) &
            &                                       )
       if (i == 1) cycle
       do j=1,i-1
          self%correlation (j+((i-2)*(i-1))/2)=+      self%covariance(i,j) &
               &                               /sqrt(                      &
               &                                     +self%covariance(i,i) &
               &                                     *self%covariance(j,j) &
               &                                    )
       end do
    end do
    if (any(abs(self%correlation) > 1.0d0)) call Error_Report('invalid correlation matrix'//{introspection:location})
    return
  end function multivariateNormalConstructorInternal

  double precision function multivariateNormalDensity(self,x,logarithmic,status) result(density)
    !!{
    Return the density of a multivariate normal distribution.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Linear_Algebra          , only : assignment(=), vector
    use :: Error                   , only : Error_Report
    use :: Interface_GSL           , only : GSL_Success
    use :: Sorting                 , only : sort
    use :: String_Handling         , only : operator(//)
    use :: ISO_Varying_String      , only : operator(//) , var_str
    implicit none
    class           (distributionFunctionMultivariateNormal), intent(inout)               :: self
    double precision                                        , intent(in   ), dimension(:) :: x
    logical                                                 , intent(in   ), optional     :: logarithmic
    integer                                                 , intent(  out), optional     :: status
    type            (vector                                )                              :: offset
    double precision                                                                      :: argumentExponential
    integer                                                                               :: status_
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    density=0.0d0
    if (present(status)) status=GSL_Success
    if (size(x) /= size(self%mean)) call Error_Report('vector `x` has incorrect size'//{introspection:location})
    offset             =+     x    &
         &              -self%mean
    argumentExponential=self%covariance_%covarianceProduct(offset,status_)
    if (status_ /= GSL_Success) then
       if (present(status)) then
          status=status_
          return
       else
          call Error_Report(var_str('covariance product failed  (GSL error ')//status_//')'//{introspection:location})
       end if
    end if
    if (logarithmic_) then
       if (.not.self%logNormalizationComputed) then
          self%logNormalization        =-0.5d0*size(self%mean)*log(2.0d0*Pi)             &
               &                        -0.5d0*self%covariance_%logarithmicDeterminant()
          self%logNormalizationComputed=.true.
       end if
       density=+self%logNormalization &
            &  -0.5d0                 &
            &  *argumentExponential
    else
       if (.not.self%normalizationComputed) then
          self%normalization        =+1.0d0                                &
               &                     /sqrt(                                &
               &                           +(+2.0d0*Pi)**size(self%mean)   &
               &                           *self%covariance_%determinant() &
               &                          )
          self%normalizationComputed=.true.
       end if
       density=+     self%normalization  &
            &  *exp(                     &
            &       -0.5d0               &
            &       *argumentExponential &
            &      )
    end if
    return
  end function multivariateNormalDensity

  double precision function multivariateNormalCumulative(self,xLow,xHigh,logarithmic,status) result(probability)
    !!{
    Return the cumulative probability of a multivariate normal distribution.
    !!}
    use :: Models_Likelihoods_Constants, only : logImprobable
    use :: Interface_GSL               , only : GSL_ERange           , GSL_ETol               , GSL_Success
    use :: Error                       , only : signalHandlerRegister, signalHandlerDeregister, signalHandlerInterface
    implicit none
    class           (distributionFunctionMultivariateNormal), intent(inout)                             :: self
    double precision                                        , intent(in   ), dimension(         :     ) :: xLow         , xHigh
    logical                                                 , intent(in   ), optional                   :: logarithmic
    integer                                                 , intent(  out), optional                   :: status
    double precision                                                       , dimension(size(self%mean)) :: offsetLow    , offsetHigh
    integer                                                                , dimension(size(self%mean)) :: infinite
    procedure       (signalHandlerInterface                ), pointer                                   :: handler
    integer                                                                                             :: maximumValues, status_
    double precision                                                                                    :: error
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    probability=0.0d0
    if (present(status)) status=GSL_Success
    if (size(xLow         ) /= size(self%mean)) call Error_Report('vector `xLow` has incorrect size' //{introspection:location})
    if (size(        xHigh) /= size(self%mean)) call Error_Report('vector `xHigh` has incorrect size'//{introspection:location})
    if (any (xLow >  xHigh)                   ) call Error_Report('`xLow` > `xHigh` is not allowed'  //{introspection:location})
    if (any (xLow == xHigh)                   ) then
       if (logarithmic_) then
          probability=logImprobable
       else
          probability=0.0d0
       end if
       return
    end if
    where (xLow  > -huge(0.0d0))
       offsetLow    =(+xLow -self%mean)/self%rootVariance
    end where
    where (xHigh < +huge(0.0d0))
       offsetHigh   =(+xHigh-self%mean)/self%rootVariance
    end where
    infinite     =-1
    where (xLow  > -huge(0.0d0))
       infinite=infinite+2
    end where
    where (xHigh < +huge(0.0d0))
       infinite=infinite+1
    end where
    maximumValues=1000*size(self%mean)
    ! Register an error handler so that we can diagnose any floating point errors.
    handler      => fpeHandlerDirect
    signalDirect =  0
    call signalHandlerRegister(handler)
    call mvndst(size(self%mean),offsetLow,offsetHigh,infinite,self%correlation,maximumValues,self%errorAbsolute,self%errorRelative,error,probability,status_)
    ! Deregister our error handler.
    call signalHandlerDeregister(handler)
    ! Check for an error in direct calculation.
    if (signalDirect == 0) then
       select case (status_)
       case (1)
          if (present(status)) then
             status=GSL_ETol
          else
             call Error_Report('requested tolerance not obtained'//{introspection:location})
          end if
       case (2)
          if (present(status)) then
             status=GSL_ERange
          else
             if (size(self%mean) <   1) call Error_Report('too few dimensions' //{introspection:location})
             if (size(self%mean) > 500) call Error_Report('too many dimensions'//{introspection:location})
          end if
       end select
       if (logarithmic_) then
          if (probability > 0.0d0) then
             probability=log(probability)
          else
             ! Integration failed to give a non-zero answer. But, logarithmic probability was requested. Attempt to provide an
             ! approximate answer using Monte Carlo integration. Use a relatively small maximum number of trials here - we only need
             ! an approximate answer.
             if (present(status)) status=GSL_Success
             probability=self%cumulativeMonteCarlo(xLow,xHigh,logarithmic,status)
          end if
       end if
    else
       ! Direct calculation failed (likely a floating point error). Try Monte Carlo approach.
       probability=self%cumulativeMonteCarlo(xLow,xHigh,logarithmic,status)
    end if
    return
  end function multivariateNormalCumulative

  double precision function multivariateNormalCumulativeMonteCarlo(self,xLow,xHigh,logarithmic,status) result(probability)
    !!{
    Return the cumulative probability of a multivariate normal distribution computed using Monte Carlo methods.
    !!}
    use :: Interface_GSL               , only : GSL_Success  , GSL_EMaxIter
    use :: Models_Likelihoods_Constants, only : logImprobable
    use :: Error                       , only : Error_Report , signalHandlerRegister, signalHandlerDeregister, signalHandlerInterface
    implicit none
    class           (distributionFunctionMultivariateNormal), intent(inout), target         :: self
    double precision                                        , intent(in   ), dimension(:  ) :: xLow                                           , xHigh
    logical                                                 , intent(in   ), optional       :: logarithmic
    integer                                                 , intent(  out), optional       :: status
    integer                                                 , parameter                     :: countMonteCarlo                      = 1000
    integer                                                 , parameter                     :: probabilityRelativeLogarithmicMaximum=  300.0d0
    class           (distributionFunctionMultivariateNormal), pointer                       :: self_
    procedure       (signalHandlerInterface                ), pointer                       :: handler
    double precision                                        , allocatable  , dimension(:  ) :: yMinimum                                       , yMaximum                       , &
         &                                                                                     xLow_                                          , xHigh_                         , &
         &                                                                                     x                                              , y                              , &
         &                                                                                     mean
    double precision                                        , allocatable  , dimension(:,:) :: covariance
    logical                                                 , allocatable  , dimension(:  ) :: useTransform
    logical                                                                                 :: failed
    integer                                                                                 :: i                                               , j                             , &
         &                                                                                     iSubset                                         , jSubset                       , &
         &                                                                                     iSample                                         , iTrial                        , &
         &                                                                                     countSubset
    double precision                                                                        :: probabilityRelativeLogarithmic                 , probabilityLogarithmicReference, &
         &                                                                                     integralCumulative                             , integralSquareCumulative       , &
         &                                                                                     errorEstimate                                  , probability_                   , &
         &                                                                                     integralEstimate
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    if (logarithmic_) then
       probability=logImprobable
    else
       probability=0.0d0
    end if
    if (present(status)) status=GSL_Success
    if (size(xLow         ) /= size(self%mean)) call Error_Report('vector `xLow` has incorrect size' //{introspection:location})
    if (size(        xHigh) /= size(self%mean)) call Error_Report('vector `xHigh` has incorrect size'//{introspection:location})
    if (any (xLow >  xHigh)                   ) call Error_Report('`xLow` > `xHigh` is not allowed'  //{introspection:location})
    if (any (xLow == xHigh)                   ) then
       if (logarithmic_) then
          probability=logImprobable
       else
          probability=0.0d0
       end if
       return
    end if
    ! Find any dimensions which are fully integrated out. We can marginalize over these trivially by dropping the relevant entries
    ! in the multivariate Gaussian.
    countSubset=count(                      &
         &             xLow  > -huge(0.0d0) &
         &            .or.                  &
         &             xHigh < +huge(0.0d0) &
         &           )
    if (countSubset == size(self%mean)) then
       ! No dimensions are fully integrated out - simply set a pointer to self since we can use it unmodified.
       self_ => self
       allocate(xLow_ ,source=xLow )
       allocate(xHigh_,source=xHigh)
    else
       ! One or more dimensions can be integrated out. Construct a new multivariate normal distribution with the marginalized
       ! dimensions removed.
       allocate(xLow_     (countSubset            ))
       allocate(xHigh_    (countSubset            ))
       allocate(mean      (countSubset            ))
       allocate(covariance(countSubset,countSubset))
       !! Iterate over rows.
       iSubset=0
       do       i=1,size(self%mean)
          if       (xLow(i) > -huge(0.0d0) .or. xHigh(i) < +huge(0.0d0)) then
             iSubset         =iSubset+1
             xLow_  (iSubset)=     xLow (i)
             xHigh_ (iSubset)=     xHigh(i)
             mean   (iSubset)=self%mean (i)
             !! Iterate over columns.
             jSubset=0
             do j=1,size(self%mean)
                if (xLow(j) > -huge(0.0d0) .or. xHigh(j) < +huge(0.0d0)) then
                   jSubset                    =jSubset+1
                   covariance(iSubset,jSubset)=self%covariance(i,j)
                end if
             end do
          end if
       end do
       allocate(distributionFunctionMultivariateNormal :: self_)
       select type (self_)
       type is (distributionFunctionMultivariateNormal)
          self_=distributionFunctionMultivariateNormal(mean,covariance,self%errorAbsolute,self%errorRelative,self%countTrialsMaximum)
       end select
    end if
    ! Allocate work arrays.
    allocate(x           (size(xLow_)))
    allocate(y           (size(xLow_)))
    allocate(yMinimum    (size(xLow_)))
    allocate(yMaximum    (size(xLow_)))
    allocate(useTransform(size(xLow_)))
    ! Find a suitable point at which to compute a reference probability.
    where (sign(1.0d0,xLow_)*sign(1.0d0,xHigh_) <= 0.0d0)
       ! Zero is spanned.
       x=0.0d0
    elsewhere
       x=min(abs(xLow_),abs(xHigh_))
    end where
    probabilityLogarithmicReference=self_%density(x,logarithmic=.true.,status=status)
    if (present(status) .and. status /= GSL_Success) then
       call cleanUp()
       return
    end if
    ! Use a tanh⁻¹ transform to allow integration over infinite intervals.
    useTransform=xLow_ == -huge(0.0d0) .or. xHigh_ == +huge(0.0d0)
    where (useTransform)
       where (xLow_  > -huge(0.0d0))
          yMinimum=+tanh(xLow_ )
       elsewhere
          yMinimum=-1.0d0
       end where
       where (xHigh_ < +huge(0.0d0))
          yMaximum=+tanh(xHigh_)
       elsewhere
          yMaximum=+1.0d0
       end where
    elsewhere
       yMinimum   =+xLow_
       yMaximum   =+xHigh_
    end where
    if (any(yMinimum >= yMaximum)) then
       if (logarithmic_) then
          probability=logImprobable
       else
          probability=0.0d0
       end if
       call cleanUp()
       return
    end if
    ! Register an error handler so that we can diagnose any floating point errors that occur during update of the reference
    ! probability. Keep copies of the work arrays to use in error reporting.
    allocate(xLow__       ,source=xLow_       )
    allocate(xHigh__      ,source=xHigh_      )
    allocate(useTransform_,source=useTransform)
    probabilityLogarithmicReference_ =  probabilityLogarithmicReference
    handler                          => fpeHandlerMonteCarlo
    call signalHandlerRegister(handler)
    !! Adjust the offset probability by any tanh transformation terms.    
    probabilityLogarithmicReference=+probabilityLogarithmicReference                                   &
         &                          +2.0d0*sum(logCosh(min(abs(xLow_),abs(xHigh_))),mask=useTransform)
    ! Deregister our handler and free copied work arrays.
    call signalHandlerDeregister(handler)
    deallocate(xLow__       )
    deallocate(xHigh__      )
    deallocate(useTransform_)    
    !! Repeatedly try to evaluate the integral. This is necessary as our first guess at a suitable probability density offset may
    !! be incorrect (i.e. some points in the integration volume may have hugely higher probability density), forcing us to update
    !! our guess and try again.
    failed=.true.
    do while (failed)
       failed=.false.
       !! Evaluate until sufficient precision is reached.
       iTrial                  =      0
       errorEstimate           =+huge(0.0d0)
       integralEstimate        =-     1.0d0
       integralCumulative      =+     0.0d0
       integralSquareCumulative=+     0.0d0
       probability             =+     0.0d0
       do while (errorEstimate > self%errorRelative*integralEstimate)
          iTrial=iTrial+1
          if (iTrial > self%countTrialsMaximum) then
             if (present(status)) then
                status=GSL_EMaxIter
                exit
             else
                call Error_Report('integral failed to converge'//{introspection:location})
             end if
          end if
          probability_=0.0d0
          do iSample=1,countMonteCarlo
             do i=1,size(self_%mean)
                y(i)=yMinimum(i)+(yMaximum(i)-yMinimum(i))*self%randomNumberGenerator_%uniformSample()
             end do
             if (any(useTransform .and. abs(y) >= 1.0d0)) then
                ! Point at infinity, zero density.
             else
                where (useTransform)
                   x=atanh(y)
                elsewhere
                   x=      y
                end where
                probabilityRelativeLogarithmic=+      self_%density                        (        x ,logarithmic=.true.      ,status=status) &
                     &                         +2.0d0*      sum                            (logCosh(x),mask       =useTransform              ) &
                     &                         -            probabilityLogarithmicReference
                if (present(status) .and. status /= GSL_Success) then
                   if (logarithmic_) then
                      probability=logImprobable
                   else
                      probability=0.0d0
                   end if
                   call cleanUp()
                   return
                end if
                if (probabilityRelativeLogarithmic > probabilityRelativeLogarithmicMaximum) then
                   ! The relative probability is too large, which means that our initial offset was too small. Update the offset
                   ! probability, and start again.
                   probabilityLogarithmicReference=probabilityRelativeLogarithmic+probabilityLogarithmicReference
                   failed=.true.
                   exit
                else
                   probability_=+    probability_                    &
                        &       +exp(probabilityRelativeLogarithmic)
                end if
             end if
          end do
          ! Accumulate the probability.
          probability      =+probability  &
               &            +probability_
          ! Estimate the error in the current estimate of the probability.
          integralCumulative      =+integralCumulative      +probability_
          integralSquareCumulative=+integralSquareCumulative+probability_**2
          integralEstimate        =+integralCumulative                       &
               &                   /dble(iTrial)
          if (iTrial > 1) then
             errorEstimate        =+sqrt(                                               &
                  &                      +(                                             &
                  &                        + integralSquareCumulative/dble(iTrial )     &
                  &                        -(integralCumulative      /dble(iTrial ))**2 &
                  &                       )                                             &
                  &                      /                            dble(iTrial-1)    &
                  &                     )
          else
             errorEstimate        =+huge(0.0d0)
          end if
       end do
    end do
    if (probability > 0.0d0) then
       probability=+    log(     probability              )  & ! Take the log of the accumulated Monte Carlo probability.
            &      -    log(dble(countMonteCarlo         ))  & ! Find the mean by "dividing" by the number of Monte Carlo samples per trial...
            &      -    log(dble(iTrial                  ))  & ! ...and by the number of Monte Carlo trials.
            &      +sum(log(     yMaximum       -yMinimum )) & ! "Multiply" by the volume of the region integrated over.
            &      +probabilityLogarithmicReference            ! "Multiply" back in the reference probability.
       ! Convert from logarithmic form if requested.
       if (.not.logarithmic_) probability=exp(probability  )
    else
       if (     logarithmic_) probability=    logImprobable
    end if
    call cleanUp()
    return

  contains
    
    subroutine cleanUp()
      !!{
      Perform clean up before exiting Monte Carlo integration of the multivariate normal distribution.
      !!}
      implicit none

      ! If a subset of dimensions was used, free the temporary distribution object that was created.
      if (countSubset < size(self%mean)) deallocate(self_)
      return
    end subroutine cleanUp
    
  end function multivariateNormalCumulativeMonteCarlo
  
  pure elemental double precision function logCosh(x)
    !!{
    Evaluate the logarithm of $\cosh x$, handling large values.
    !!}
    implicit none
    double precision, intent(in   ) :: x
    double precision, parameter     :: xMaximum=300.0d0
    
    if (abs(x) < xMaximum) then
       logCosh=log(cosh(x))
    else
       logCosh=abs(x)-log(2.0d0)
    end if
    return
  end function logCosh

  subroutine fpeHandlerMonteCarlo(signal)
    !!{
    Report useful information if a floating point error occurs.
    !!}
    use :: Display           , only : displayIndent , displayUnindent, displayMessage, verbosityLevelSilent
    implicit none
    integer          , intent(in   ) :: signal
    integer                          :: i
    character(len=48)                :: message
    !$GLC attributes unused :: signal
    
    call displayIndent("multivariate normal CDF MCMC evaluation - reference probability update error occured",verbosityLevelSilent)
    write (message,'(a,e14.6)') "probabilityLogarithmicReference = ",probabilityLogarithmicReference_
    call displayMessage(message,verbosityLevelSilent)
    do i=1,size(xLow__)
       if (useTransform_(i)) then
          write (message,'(i3," ",e14.6," ",e14.6," ",e14.6)') i,xLow__(i),xHigh__(i),logCosh(min(abs(xlow__(i)),abs(xHigh__(i))))
       else
          write (message,'(i3," ",e14.6," ",e14.6," ",a3)'   ) i,xLow__(i),xHigh__(i),"n/a"
       end if
       call displayMessage(message,verbosityLevelSilent)
    end do
    call displayUnindent("done",verbosityLevelSilent)
    return
  end subroutine fpeHandlerMonteCarlo

  subroutine fpeHandlerDirect(signal)
    !!{
    Handle floating point exceptions when evaluating the cumulative distribution function.
    !!}
    implicit none
    integer, intent(in   ) :: signal

    signalDirect=signal
    return
  end subroutine fpeHandlerDirect

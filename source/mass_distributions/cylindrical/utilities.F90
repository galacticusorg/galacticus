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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides utilities shared by cylindrical mass distributions.
!!}

module Mass_Distributions_Cylindrical_Utilities
  !!{RST
  Provides utilities shared by cylindrical mass distributions: the vertical kernel, :math:`K(k,z)`, of a disk with a
  :math:`\mathrm{sech}^2` vertical density profile, and the evaluation of Hankel transforms (integrals over wavenumber of
  Bessel functions of the first kind) such as those which give the potential and acceleration of a thick disk
  :cite:p:`kuijken_mass_1989`.
  !!}
  use :: Numerical_Integration, only : integrator
  implicit none
  private
  public :: kernelVerticalSechSquared             , Hankel_Transform_Integrate, Kernel_Vertical_Integrators_Allocate, &
       &    Kernel_Vertical_Integrators_Deallocate

  type :: kernelVerticalSechSquared
     !!{RST
     The vertical kernel, :math:`K(k,z)`, and its first and second :math:`z` derivatives, for a disk with a
     :math:`\mathrm{sech}^2` vertical density profile, tabulated as a function of wavenumber at a set of heights. All lengths
     are in units of the radial scale length of the disk.
     !!}
     private
     double precision                                :: scaleHeight                  , wavenumberLogarithmicMinimum, &
          &                                             wavenumberLogarithmicInterval
     integer                                         :: countWavenumbers
     double precision, allocatable, dimension(:    ) :: heights
     double precision, allocatable, dimension(:,:,:) :: table
   contains
     !![
     <methods docformat="rst">
       <method description="Return the vertical kernel (``indexKernel=1``), or its first (``indexKernel=2``) or second (``indexKernel=3``) :math:`z` derivative, at the given wavenumber and the height with index ``indexHeight`` in the table." method="value"/>
     </methods>
     !!]
     procedure :: value => kernelVerticalSechSquaredValue
  end type kernelVerticalSechSquared

  interface kernelVerticalSechSquared
     !!{RST
     Constructors for the :math:`\mathrm{sech}^2` vertical kernel class.
     !!}
     module procedure kernelVerticalSechSquaredConstructor
  end interface kernelVerticalSechSquared

  ! Range and density of the tabulation of the vertical kernel in wavenumber.
  double precision                , parameter                :: wavenumberKernelMinimum  = 1.0d-3, wavenumberKernelMaximum =1.0d7
  double precision                , parameter                :: wavenumberKernelPerDecade=32.0d+0

  ! Integrators, and the scaled wavenumber and height, used in direct evaluation of the vertical kernel.
  type            (integrator    ), allocatable              :: integratorKernel                 , integratorKernelGradient     , &
       &                                                        integratorKernelGradient2
  double precision                                           :: wavenumberScaled                 , heightScaled
  !$omp threadprivate(integratorKernel,integratorKernelGradient,integratorKernelGradient2,wavenumberScaled,heightScaled)

contains

  function kernelVerticalSechSquaredConstructor(scaleHeight,heights) result(self)
    !!{RST
    Tabulate the vertical kernel, :math:`K(k,z)`, and its derivatives, as a function of wavenumber at each of the given
    ``heights``, for a disk with dimensionless scale height ``scaleHeight``. The kernel depends only on wavenumber and height, so
    tabulating it once avoids recomputing it (each evaluation of which requires a numerical integration) for every radius, and for
    every quantity computed from it.
    !!}
    implicit none
    type            (kernelVerticalSechSquared)                              :: self
    double precision                           , intent(in   )               :: scaleHeight
    double precision                           , intent(in   ), dimension(:) :: heights
    integer                                                                  :: iHeight    , iWavenumber
    double precision                                                         :: wavenumber

    self%scaleHeight                  =scaleHeight
    self%heights                      =heights
    self%countWavenumbers             =int(log10(wavenumberKernelMaximum/wavenumberKernelMinimum)*wavenumberKernelPerDecade)+1
    self%wavenumberLogarithmicMinimum =log(wavenumberKernelMinimum)
    self%wavenumberLogarithmicInterval=log(wavenumberKernelMaximum/wavenumberKernelMinimum)/dble(self%countWavenumbers-1)
    allocate(self%table(self%countWavenumbers,size(heights),3))
    !$omp parallel private(iWavenumber,wavenumber)
    call Kernel_Vertical_Integrators_Allocate()
    !$omp do schedule(dynamic)
    do iHeight=1,size(heights)
       do iWavenumber=1,self%countWavenumbers
          wavenumber=exp(self%wavenumberLogarithmicMinimum+dble(iWavenumber-1)*self%wavenumberLogarithmicInterval)
          self%table(iWavenumber,iHeight,1)=kernelVerticalDirect         (self,wavenumber,heights(iHeight))
          self%table(iWavenumber,iHeight,2)=kernelVerticalGradientDirect (self,wavenumber,heights(iHeight))
          self%table(iWavenumber,iHeight,3)=kernelVerticalGradient2Direct(self,wavenumber,heights(iHeight))
       end do
    end do
    !$omp end do
    call Kernel_Vertical_Integrators_Deallocate()
    !$omp end parallel
    return
  end function kernelVerticalSechSquaredConstructor

  subroutine Kernel_Vertical_Integrators_Allocate()
    !!{RST
    Allocate and construct the (thread private) integrators used in direct evaluation of the vertical kernel. These are needed
    whenever the kernel is evaluated outside of its tabulated range of wavenumber, so must be allocated by each thread which
    evaluates the kernel.
    !!}
    implicit none

    allocate(integratorKernel         )
    allocate(integratorKernelGradient )
    allocate(integratorKernelGradient2)
    integratorKernel         =integrator(kernelIntegrand         ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
    integratorKernelGradient =integrator(kernelGradientIntegrand ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
    integratorKernelGradient2=integrator(kernelGradient2Integrand,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
    return
  end subroutine Kernel_Vertical_Integrators_Allocate

  subroutine Kernel_Vertical_Integrators_Deallocate()
    !!{RST
    Destroy the (thread private) integrators used in direct evaluation of the vertical kernel.
    !!}
    implicit none

    deallocate(integratorKernel         )
    deallocate(integratorKernelGradient )
    deallocate(integratorKernelGradient2)
    return
  end subroutine Kernel_Vertical_Integrators_Deallocate

  double precision function kernelVerticalSechSquaredValue(self,wavenumber,indexHeight,indexKernel) result(kernel)
    !!{RST
    Interpolate the vertical kernel (``indexKernel=1``), or its first (``indexKernel=2``) or second (``indexKernel=3``)
    :math:`z` derivative, at the height with index ``indexHeight`` in the table. Four-point Lagrange interpolation in
    :math:`\log k` is used, as this requires no state (and so is safe for concurrent use by multiple OpenMP threads). Outside of
    the tabulated range of wavenumber the kernel is evaluated directly, which requires that the calling thread has allocated the
    integrators (see ``Kernel_Vertical_Integrators_Allocate``).
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (kernelVerticalSechSquared), intent(inout) :: self
    double precision                           , intent(in   ) :: wavenumber
    integer                                    , intent(in   ) :: indexHeight, indexKernel
    double precision                           , dimension(4)  :: weight
    double precision                                           :: position
    integer                                                    :: i

    if (wavenumber < wavenumberKernelMinimum .or. wavenumber > wavenumberKernelMaximum) then
       select case (indexKernel)
       case (1)
          kernel=kernelVerticalDirect         (self,wavenumber,self%heights(indexHeight))
       case (2)
          kernel=kernelVerticalGradientDirect (self,wavenumber,self%heights(indexHeight))
       case (3)
          kernel=kernelVerticalGradient2Direct(self,wavenumber,self%heights(indexHeight))
       case default
          kernel=0.0d0
          call Error_Report('incorrect kernel index'//{introspection:location})
       end select
    else
       ! Find the position in the table, and the index of the second of the four points used in the interpolation.
       position=+(                                         &
            &     +log(wavenumber)                         &
            &     -    self%wavenumberLogarithmicMinimum   &
            &    )                                         &
            &   /      self%wavenumberLogarithmicInterval
       i       =min(max(int(position)+1,2),self%countWavenumbers-2)
       position=position-dble(i-1)
       ! Evaluate the Lagrange interpolation weights.
       weight  =[                                                           &
            &    -(position      )*(position-1.0d0)*(position-2.0d0)/6.0d0, &
            &    +(position+1.0d0)*(position-1.0d0)*(position-2.0d0)/2.0d0, &
            &    -(position+1.0d0)*(position      )*(position-2.0d0)/2.0d0, &
            &    +(position+1.0d0)*(position      )*(position-1.0d0)/6.0d0  &
            &   ]
       kernel=sum(weight*self%table(i-1:i+2,indexHeight,indexKernel))
    end if
    return
  end function kernelVerticalSechSquaredValue

  double precision function kernelVerticalDirect(self,wavenumber,height)
    !!{RST
    The vertical kernel, :math:`K(k,z)`, appearing in the expression for the potential of the disk,
    :math:`\Phi(R,z) = -2\pi G \int_0^\infty \mathrm{d}k \, J_0(kR) \tilde{\Sigma}(k) K(k,z)`, where :math:`\tilde{\Sigma}(k)`
    is the Hankel transform of the surface density. :math:`K(k,z)` is the convolution of the normalized vertical density
    profile, :math:`\rho_z(z)=\mathrm{sech}^2(z/h)/2h`, with the razor-thin kernel, :math:`\exp(-k|z|)`. Writing
    :math:`u=z/h` and :math:`\kappa=kh` this can be expressed as a single Laplace-type integral:

    .. math::

       K(k,z) = {1 \over 2} \int_0^\infty \mathrm{d}x \, \mathrm{e}^{-\kappa x} \left[ \mathrm{sech}^2(u+x) + \mathrm{sech}^2(u-x) \right],

    which is normalized such that :math:`K(0,z)=1`. This function evaluates the kernel directly, by numerical integration.
    !!}
    implicit none
    class           (kernelVerticalSechSquared), intent(inout) :: self
    double precision                           , intent(in   ) :: wavenumber, height

    kernelVerticalDirect=+0.5d0                                                    &
         &               *kernelIntegrate(self,wavenumber,height,integratorKernel)
    return
  end function kernelVerticalDirect

  double precision function kernelVerticalGradientDirect(self,wavenumber,height)
    !!{RST
    The :math:`z` derivative of the vertical kernel, :math:`K(k,z)`:

    .. math::

       {\partial K \over \partial z} = {k \over 2} \int_0^\infty \mathrm{d}x \, \mathrm{e}^{-\kappa x} \left[ \mathrm{sech}^2(u+x) - \mathrm{sech}^2(u-x) \right].
    !!}
    implicit none
    class           (kernelVerticalSechSquared), intent(inout) :: self
    double precision                           , intent(in   ) :: wavenumber, height

    kernelVerticalGradientDirect=+0.5d0                                                            &
         &                       *wavenumber                                                       &
         &                       *kernelIntegrate(self,wavenumber,height,integratorKernelGradient)
    return
  end function kernelVerticalGradientDirect

  double precision function kernelVerticalGradient2Direct(self,wavenumber,height)
    !!{RST
    The second :math:`z` derivative of the vertical kernel, :math:`K(k,z)`, found by differentiating the integrand of
    :math:`\partial K/\partial z` with respect to :math:`u`:

    .. math::

       {\partial^2 K \over \partial z^2} = {k \over 2h} \int_0^\infty \mathrm{d}x \, \mathrm{e}^{-\kappa x} \left[ S^\prime(u+x) - S^\prime(u-x) \right],

    where :math:`S(y)=\mathrm{sech}^2 y`. This is equivalent to :math:`k^2 K - 2 k \rho_z(z)` (which follows from
    :math:`(\partial^2/\partial z^2 - k^2)K = -2k\rho_z`), but avoids the catastrophic cancellation between those two terms
    which occurs for :math:`\kappa \gg 1`.
    !!}
    implicit none
    class           (kernelVerticalSechSquared), intent(inout) :: self
    double precision                           , intent(in   ) :: wavenumber, height

    kernelVerticalGradient2Direct=+0.5d0                                                             &
         &                        *wavenumber                                                        &
         &                        /self%scaleHeight                                                  &
         &                        *kernelIntegrate(self,wavenumber,height,integratorKernelGradient2)
    return
  end function kernelVerticalGradient2Direct

  double precision function kernelIntegrate(self,wavenumber,height,integrator_)
    !!{RST
    Evaluate the integral over :math:`x` appearing in the vertical kernel, :math:`K(k,z)`, or its derivatives, using the
    provided integrator.

    The integrand has structure on three scales: near :math:`x=0` (where it decays as :math:`\mathrm{e}^{-(\kappa+2)x}`), and
    near :math:`x=u` (where :math:`\mathrm{sech}^2(u-x)` peaks, with width of order unity). For large :math:`\kappa` the
    structure near :math:`x=0` (and the onset of the peak at :math:`x=u`) has width :math:`1/\kappa`. Since an adaptive
    integrator can entirely miss a narrow feature lying between its initial sample points (and so report a confident, but
    incorrect, result) the range of integration is broken at each of these scales, so that every feature lies close to the
    boundary of some sub-interval.
    !!}
    implicit none
    class           (kernelVerticalSechSquared), intent(inout) :: self
    double precision                           , intent(in   ) :: wavenumber             , height
    type            (integrator               ), intent(inout) :: integrator_
    ! Range (in x) beyond which sech²(u±x) has decayed by a factor e⁻⁵⁰ from its peak.
    double precision                           , parameter     :: rangeSech       =25.0d0
    ! Range (in κx) beyond which exp(-κx) has decayed by a factor e⁻⁵⁰.
    double precision                           , parameter     :: rangeExponential=50.0d0
    double precision                           , dimension(7)  :: limits
    double precision                                           :: rangeDecay
    integer                                                    :: i

    wavenumberScaled=+wavenumber       &
         &           *self%scaleHeight
    heightScaled    =+height           &
         &           /self%scaleHeight
    if (wavenumberScaled*rangeSech > rangeExponential) then
       rangeDecay=+rangeExponential &
            &     /wavenumberScaled
    else
       rangeDecay=+rangeSech
    end if
    limits=[                                                         &
         &  0.0d0                                                  , &
         &  min(heightScaled,rangeSech,rangeDecay)                 , &
         &  min(heightScaled,rangeSech           )                 , &
         &  max(min(heightScaled,rangeSech),heightScaled-rangeSech), &
         &  heightScaled                                           , &
         &  heightScaled+rangeDecay                                , &
         &  heightScaled+rangeSech                                   &
         & ]
    kernelIntegrate=0.0d0
    do i=1,size(limits)-1
       if (limits(i+1) > limits(i)) kernelIntegrate=+kernelIntegrate                              &
            &                                       +integrator_%integrate(limits(i),limits(i+1))
    end do
    return
  end function kernelIntegrate

  double precision function kernelIntegrand(x)
    !!{RST
    Integrand for the vertical kernel, :math:`K(k,z)`.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    kernelIntegrand=+exp(-wavenumberScaled*x)      &
         &          *(                             &
         &            +sechSquared(heightScaled+x) &
         &            +sechSquared(heightScaled-x) &
         &           )
    return
  end function kernelIntegrand

  double precision function kernelGradientIntegrand(x)
    !!{RST
    Integrand for the :math:`z` derivative of the vertical kernel, :math:`K(k,z)`.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    kernelGradientIntegrand=+exp(-wavenumberScaled*x)      &
         &                  *(                             &
         &                    +sechSquared(heightScaled+x) &
         &                    -sechSquared(heightScaled-x) &
         &                   )
    return
  end function kernelGradientIntegrand

  double precision function kernelGradient2Integrand(x)
    !!{RST
    Integrand for the second :math:`z` derivative of the vertical kernel, :math:`K(k,z)`.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    kernelGradient2Integrand=+exp(-wavenumberScaled*x)                &
         &                   *(                                       &
         &                     +sechSquaredDerivative(heightScaled+x) &
         &                     -sechSquaredDerivative(heightScaled-x) &
         &                    )
    return
  end function kernelGradient2Integrand

  double precision function sechSquared(y)
    !!{RST
    Evaluate :math:`\mathrm{sech}^2 y`, written in terms of :math:`\mathrm{e}^{-2|y|}` so that it can not overflow at large
    :math:`|y|`.
    !!}
    implicit none
    double precision, intent(in   ) :: y
    double precision                :: expTerm

    expTerm    =+exp(-2.0d0*abs(y))
    sechSquared=+4.0d0     &
         &      *expTerm   &
         &      /(         &
         &        +1.0d0   &
         &        +expTerm &
         &       )**2
    return
  end function sechSquared

  double precision function sechSquaredDerivative(y)
    !!{RST
    Evaluate :math:`\mathrm{d}\,\mathrm{sech}^2 y/\mathrm{d}y = -2\,\mathrm{sech}^2 y \tanh y`, written in terms of
    :math:`\mathrm{e}^{-2|y|}` so that it can not overflow at large :math:`|y|`.
    !!}
    implicit none
    double precision, intent(in   ) :: y
    double precision                :: expTerm

    expTerm              =+exp(-2.0d0*abs(y))
    sechSquaredDerivative=-2.0d0          &
         &                *sechSquared(y) &
         &                *sign(1.0d0,y)  &
         &                *(              &
         &                  +1.0d0        &
         &                  -expTerm      &
         &                 )              &
         &                /(              &
         &                  +1.0d0        &
         &                  +expTerm      &
         &                 )
    return
  end function sechSquaredDerivative

  double precision function Hankel_Transform_Integrate(integrator_,order,radius) result(integral)
    !!{RST
    Evaluate the integral over wavenumber, from zero to infinity, of the integrand of the given integrator, which is assumed to
    contain a Bessel function of the first kind of the given ``order`` with argument :math:`kR`, where :math:`R=`\ ``radius``.
    The integral is accumulated over the intervals between successive zeros of the Bessel function. The partial sums so obtained
    form a (typically alternating) sequence which converges slowly - for the components of the tidal tensor of a disk, whose
    integrands carry additional powers of :math:`k`, many thousands of intervals would be needed for the contribution from an
    interval to fall below :math:`10^{-6}` of the total. The limit of the sequence is therefore estimated using the
    :math:`\epsilon`-algorithm of :cite:t:`wynn_device_1956` applied to the most recent partial sums, and accumulation stops
    once successive estimates agree to a relative precision of :math:`10^{-8}`. Accumulation also stops (and the partial sum is
    returned) if the contribution from an interval falls below :math:`10^{-6}` of the total, which occurs first if the integrand
    decays rapidly.
    !!}
    use :: Bessel_Functions, only : Bessel_Function_J0_Zero, Bessel_Function_J1_Zero, Bessel_Function_Jn_Zero
    use :: Error           , only : Error_Report
    implicit none
    type            (integrator), intent(inout)               :: integrator_
    integer                     , intent(in   )               :: order
    double precision            , intent(in   )               :: radius
    integer                     , parameter                   :: countPartialSums=16     , countIntervalsMinimum=8
    double precision            , parameter                   :: toleranceDirect = 1.0d-6, toleranceAccelerated =1.0d-8
    double precision            , dimension(countPartialSums) :: partialSums
    double precision            , dimension(3               ) :: estimates
    double precision                                          :: wavenumberLow           , wavenumberHigh              , &
         &                                                       integralDelta
    integer                                                   :: iBesselZero             , countStored
    logical                                                   :: converged

    integral      =0.0d0
    wavenumberHigh=0.0d0
    iBesselZero   =0
    countStored   =0
    estimates     =huge(0.0d0)
    converged     =.false.
    do while (.not.converged)
       iBesselZero  =+iBesselZero    &
            &        +1
       wavenumberLow=+wavenumberHigh
       select case (order)
       case (0)
          wavenumberHigh=+Bessel_Function_J0_Zero(      iBesselZero) &
               &         /radius
       case (1)
          wavenumberHigh=+Bessel_Function_J1_Zero(      iBesselZero) &
               &         /radius
       case (2)
          wavenumberHigh=+Bessel_Function_Jn_Zero(2.0d0,iBesselZero) &
               &         /radius
       case default
          wavenumberHigh=+0.0d0
          call Error_Report('incorrect Bessel function order'//{introspection:location})
       end select
       integralDelta=+integrator_%integrate(wavenumberLow,wavenumberHigh)
       converged    = abs(integralDelta) < toleranceDirect*abs(integral)
       integral     =+integral      &
            &        +integralDelta
       ! Store the partial sum, retaining only the most recent.
       if (countStored < countPartialSums) then
          countStored             =countStored+1
       else
          partialSums(1:countStored-1)=partialSums(2:countStored)
       end if
       partialSums(countStored)=integral
       ! Estimate the limit of the sequence of partial sums, and test for convergence of successive estimates.
       if (.not.converged .and. iBesselZero >= countIntervalsMinimum) then
          estimates=[estimates(2:3),wynnEpsilon(partialSums(1:countStored))]
          if     (                                                                         &
               &   abs(estimates(3)-estimates(2)) < toleranceAccelerated*abs(estimates(3)) &
               &  .and.                                                                    &
               &   abs(estimates(3)-estimates(1)) < toleranceAccelerated*abs(estimates(3)) &
               & ) then
             integral =estimates(3)
             converged=.true.
          end if
       end if
    end do
    return
  end function Hankel_Transform_Integrate

  double precision function wynnEpsilon(sequence) result(estimate)
    !!{RST
    Estimate the limit of the given sequence using the :math:`\epsilon`-algorithm of :cite:t:`wynn_device_1956`. The estimate is
    the final element of the highest even-order column of the :math:`\epsilon` table which can be constructed. Construction
    stops early if two successive elements of a column are equal (to within rounding), as the sequence has then converged.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:             ) :: sequence
    double precision               , dimension(size(sequence)) :: columnPrevious            , columnCurrent, &
         &                                                        columnNext
    double precision, parameter                                :: differenceMinimum=1.0d-280
    double precision                                           :: difference
    integer                                                    :: countSequence             , iColumn      , &
         &                                                        i
    logical                                                    :: degenerate

    countSequence =size(sequence)
    estimate      =sequence(countSequence)
    columnPrevious=0.0d0
    columnCurrent =sequence
    degenerate    =.false.
    do iColumn=1,countSequence-1
       do i=1,countSequence-iColumn
          difference=columnCurrent(i+1)-columnCurrent(i)
          if (abs(difference) <= max(epsilon(0.0d0)*max(abs(columnCurrent(i+1)),abs(columnCurrent(i))),differenceMinimum)) then
             degenerate=.true.
             exit
          end if
          columnNext(i)=columnPrevious(i+1)+1.0d0/difference
       end do
       if (degenerate) then
          ! Two successive elements of this column are equal. If this is an even-order column, this is the limit.
          if (mod(iColumn-1,2) == 0) estimate=columnCurrent(i+1)
          exit
       end if
       if (mod(iColumn,2) == 0) estimate=columnNext(countSequence-iColumn)
       columnPrevious(1:countSequence-iColumn+1)=columnCurrent(1:countSequence-iColumn+1)
       columnCurrent (1:countSequence-iColumn  )=columnNext   (1:countSequence-iColumn  )
    end do
    return
  end function wynnEpsilon

end module Mass_Distributions_Cylindrical_Utilities

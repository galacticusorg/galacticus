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
Contains a module which implements adaptive multi-dimensional integration (cubature).
!!}

module Numerical_Integration_nD
  !!{RST
  Implements adaptive multi-dimensional integration (cubature) over hyper-rectangular regions.

  The alternative to a genuine cubature is to nest one-dimensional integrators, one per dimension---which is what the
  computational domain volume integrators do. That is expensive for a reason no amount of tuning removes: a nest of a
  :math:`K`-point rule in :math:`d` dimensions cannot cost fewer than :math:`K^d` integrand evaluations, however smooth the
  integrand, because each level must apply its rule at least once. With the 61-point Gauss-Kronrod rule that is 3,721
  evaluations in two dimensions and 226,981 in three, before any subdivision at all.

  This module instead uses the degree-7 rule of :cite:t:`genz_remarks_1980`, which requires only :math:`2^d+2d^2+2d+1` points
  per region---17 in two dimensions and 33 in three. The rule carries an embedded degree-5 rule sharing all of its points, so an
  error estimate costs nothing extra.

  Regions are refined adaptively: the region of largest error estimate is repeatedly bisected, as the one-dimensional
  integrators of :galacticus-mod:`Numerical_Integration2` bisect intervals, with two differences.

  #. A region is split along a single coordinate rather than into :math:`2^d` children. Splitting into :math:`2^d` children
     multiplies the region count by eight in three dimensions at every refinement, and spends effort in directions along which
     the integrand is already smooth. The coordinate chosen is that of largest fourth difference, which the rule computes for
     free as part of its own error estimate.
  #. Regions are held in a binary heap rather than in the sorted linked list used in one dimension, whose insertion sort costs
     :math:`\mathcal{O}(N)` per insertion. That is tolerable in one dimension but not here, where region counts can reach
     millions for a difficult integrand.

  A degree-7 rule is the right tool at loose tolerance. At tight tolerance it loses to a nest of one-dimensional Gauss-Kronrod
  rules, which pay their :math:`K^d` floor once and then exploit the smoothness of the integrand with a much higher order rule.
  It should not be assumed to be faster in all circumstances.
  !!}
  use :: Numerical_Integration2, only : integrator2
  implicit none
  private
  public :: integratorGenzMalikND, integrandND

  abstract interface
     double precision function integrandND(x)
       !!{RST
       Interface for multi-dimensional integrands. The argument ``x`` holds the coordinates of the point at which the integrand is
       to be evaluated, and has extent equal to the dimensionality of the integral.
       !!}
       double precision, intent(in   ), dimension(:) :: x
     end function integrandND
  end interface

  type :: region
     !!{RST
     A hyper-rectangular region of the integration domain, with the estimates of the integral over it and of the error in that
     estimate, and the coordinate along which it is to be split if it is refined.
     !!}
     private
     double precision, allocatable, dimension(:) :: lower         , upper
     double precision                            :: integral      , error
     integer                                     :: dimensionSplit
  end type region

  type, extends(integrator2) :: integratorGenzMalikND
     !!{RST
     A multi-dimensional integrator using the adaptive degree-7 cubature rule of :cite:t:`genz_remarks_1980`.
     !!}
     private
     procedure       (integrandND), pointer, nopass :: integrand        => null()
     integer                                        :: countDimensions           , countEvaluationsMaximum, &
          &                                            countEvaluations
     ! Weights of the degree-7 rule and of the embedded degree-5 rule, for, in order: the center point; points with a single
     ! non-zero coordinate at ±λ₂; points with a single non-zero coordinate at ±λ₃; points with two non-zero coordinates at ±λ₄;
     ! and the 2ᵈ points with every coordinate at ±λ₅. The degree-5 rule does not use the last of these sets.
     double precision                               :: weight7Center             , weight7Lambda2         , &
          &                                            weight7Lambda3            , weight7Lambda4         , &
          &                                            weight7Lambda5            , weight5Center          , &
          &                                            weight5Lambda2            , weight5Lambda3         , &
          &                                            weight5Lambda4
   contains
     !![
     <methods docformat="rst">
       <method description="Set the integrand function to be integrated."                                 method="integrandSet"     />
       <method description="Evaluate the integral over a hyper-rectangular region."                       method="evaluate"         />
       <method description="Return the number of integrand evaluations used by the most recent integral." method="countEvaluated"   />
       <method description="Apply the cubature rules over a single region."                               method="regionEvaluate"   />
       <method description="Evaluate the integrand, counting the evaluation."                             method="integrandEvaluate"/>
       <method description="Test whether estimates of the integral and its error meet the tolerances."    method="hasConverged"     />
     </methods>
     !!]
     procedure :: integrandSet      => genzMalikNDIntegrandSet
     procedure :: evaluate          => genzMalikNDEvaluate
     procedure :: countEvaluated    => genzMalikNDCountEvaluated
     procedure :: regionEvaluate    => genzMalikNDRegionEvaluate
     procedure :: integrandEvaluate => genzMalikNDIntegrandEvaluate
     procedure :: hasConverged      => genzMalikNDHasConverged
  end type integratorGenzMalikND

  interface integratorGenzMalikND
     !!{RST
     Constructors for the :galacticus-class:`integratorGenzMalikND` class.
     !!}
     module procedure genzMalikNDConstructor
  end interface integratorGenzMalikND

  ! Squares of the generators of the rule---these are how the generators are specified, and λ₃²/λ₂² is used directly in forming
  ! the fourth differences which select the coordinate along which a region is split.
  double precision, parameter :: lambda2Squared    =9.0d0/70.0d0                 , lambda3Squared=9.0d0/10.0d0, &
       &                         lambda4Squared    =9.0d0/10.0d0                 , lambda5Squared=9.0d0/19.0d0
  double precision, parameter :: ratioLambdaSquared=lambda3Squared/lambda2Squared

contains

  function genzMalikNDConstructor(countDimensions,integrand,toleranceAbsolute,toleranceRelative,countEvaluationsMaximum) result(self)
    !!{RST
    Constructor for the :galacticus-class:`integratorGenzMalikND` class. ``countDimensions`` is the dimensionality of the integral,
    which must be at least two---in one dimension the integrators of :galacticus-mod:`Numerical_Integration2` should be used
    instead. At least one of ``toleranceAbsolute`` and ``toleranceRelative`` must be given. ``countEvaluationsMaximum`` limits the number
    of integrand evaluations performed before the integration is abandoned.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (integratorGenzMalikND)                          :: self
    integer                                , intent(in   )           :: countDimensions
    procedure       (integrandND          )                          :: integrand
    double precision                       , intent(in   ), optional :: toleranceAbsolute      , toleranceRelative
    integer                                , intent(in   ), optional :: countEvaluationsMaximum
    double precision                                                 :: countDimensions_
    !![
    <optionalArgument name="countEvaluationsMaximum" defaultsTo="100000000"/>
    !!]

    if (countDimensions < 2) call Error_Report('a multi-dimensional integral must have at least two dimensions'//{introspection:location})
    self%countDimensions        =  countDimensions
    self%countEvaluationsMaximum=  countEvaluationsMaximum_
    self%countEvaluations       =  0
    self%integrand              => integrand
    call self%toleranceSet(toleranceAbsolute=toleranceAbsolute,toleranceRelative=toleranceRelative)
    ! Construct the weights of the two rules. These are the closed forms in the dimensionality given by
    ! :cite:t:`genz_remarks_1980`, normalized such that the integral is the volume of the region multiplied by the weighted sum
    ! of the integrand over the points of the rule.
    countDimensions_   =dble(countDimensions)
    self%weight7Center =+(12824.0d0-9120.0d0*countDimensions_+400.0d0*countDimensions_**2)/19683.0d0
    self%weight7Lambda2=+   980.0d0                                                       / 6561.0d0
    self%weight7Lambda3=+( 1820.0d0- 400.0d0*countDimensions_                            )/19683.0d0
    self%weight7Lambda4=+   200.0d0                                                       /19683.0d0
    self%weight7Lambda5=+  6859.0d0                                                       /19683.0d0/2.0d0**countDimensions
    self%weight5Center =+(  729.0d0- 950.0d0*countDimensions_+ 50.0d0*countDimensions_**2)/  729.0d0
    self%weight5Lambda2=+   245.0d0                                                       /  486.0d0
    self%weight5Lambda3=+(  265.0d0- 100.0d0*countDimensions_                            )/ 1458.0d0
    self%weight5Lambda4=+    25.0d0                                                       /  729.0d0
    return
  end function genzMalikNDConstructor

  subroutine genzMalikNDIntegrandSet(self,integrand)
    !!{RST
    Set the integrand for a :galacticus-class:`integratorGenzMalikND` integrator.
    !!}
    implicit none
    class    (integratorGenzMalikND), intent(inout) :: self
    procedure(integrandND          )                :: integrand

    self%integrand => integrand
    return
  end subroutine genzMalikNDIntegrandSet

  integer function genzMalikNDCountEvaluated(self)
    !!{RST
    Return the number of integrand evaluations used by the most recent call to  evaluate.
    !!}
    implicit none
    class(integratorGenzMalikND), intent(in   ) :: self

    genzMalikNDCountEvaluated=self%countEvaluations
    return
  end function genzMalikNDCountEvaluated

  logical function genzMalikNDHasConverged(self,integral,error)
    !!{RST
    Return true if the given estimates of the integral and of its error meet the requested tolerances.
    !!}
    implicit none
    class           (integratorGenzMalikND), intent(in   ) :: self
    double precision                       , intent(in   ) :: integral, error

    genzMalikNDHasConverged=  abs(error) <= self%toleranceAbsolute               &
         &                  .or.                                                 &
         &                    abs(error) <= self%toleranceRelative*abs(integral)
    return
  end function genzMalikNDHasConverged

  double precision function genzMalikNDEvaluate(self,lower,upper,status) result(integral)
    !!{RST
    Evaluate a multi-dimensional integral over the hyper-rectangular region with the given  lower and  upper bounds.

    The region of largest error estimate is repeatedly split in two along a single coordinate until the total error estimate
    meets the requested tolerances. Regions are held in a binary heap ordered on their error estimates, so the region to refine
    is always at its root.
    !!}
    use :: Error, only : Error_Report            , errorStatusSuccess
    use :: Error, only : errorStatusMaxIterations
    implicit none
    class           (integratorGenzMalikND), intent(inout)               :: self
    double precision                       , intent(in   ), dimension(:) :: lower         , upper
    integer                                , intent(  out), optional     :: status
    type            (region               ), allocatable  , dimension(:) :: heap          , heapTemporary
    type            (region               )                              :: regionParent  , regionChild1 , &
         &                                                                  regionChild2
    integer                                                              :: countRegions  , sizeHeap     , &
         &                                                                  dimensionSplit
    double precision                                                     :: error         , boundarySplit
    logical                                                              :: converged

    if     (                                                                                                                    &
         &   size(lower) /= self%countDimensions                                                                                &
         &  .or.                                                                                                                &
         &   size(upper) /= self%countDimensions                                                                                &
         & ) call Error_Report('bounds must have extent equal to the dimensionality of the integral'//{introspection:location})
    if (present(status)) status=errorStatusSuccess
    self%countEvaluations=0
    ! Begin with a single region spanning the whole domain.
    sizeHeap    =16
    countRegions= 1
    allocate(heap(sizeHeap))
    heap(1)%lower=lower
    heap(1)%upper=upper
    call self%regionEvaluate(heap(1))
    integral =heap(1)%integral
    error    =heap(1)%error
    converged=self%hasConverged(integral,error)
    do while (.not.converged)
       ! Abandon the integration if the budget of integrand evaluations is exhausted.
       if (self%countEvaluations >= self%countEvaluationsMaximum) then
          if (present(status)) then
             status=errorStatusMaxIterations
             exit
          else
             call Error_Report('maximum number of integrand evaluations exceeded'//{introspection:location})
          end if
       end if
       ! Pop the region of largest error estimate from the root of the heap.
       regionParent=heap(1)
       if (countRegions > 1) heap(1)=heap(countRegions)
       countRegions=countRegions-1
       call heapSiftDown(heap,countRegions,1)
       ! Split that region in two along the coordinate of largest fourth difference.
       dimensionSplit                    =regionParent%dimensionSplit
       boundarySplit                     =+0.5d0                                &
            &                             *(                                    &
            &                               +regionParent%lower(dimensionSplit) &
            &                               +regionParent%upper(dimensionSplit) &
            &                              )
       regionChild1%lower                =   regionParent%lower
       regionChild1%upper                =   regionParent%upper
       regionChild2%lower                =   regionParent%lower
       regionChild2%upper                =   regionParent%upper
       regionChild1%upper(dimensionSplit)=boundarySplit
       regionChild2%lower(dimensionSplit)=boundarySplit
       call self%regionEvaluate(regionChild1)
       call self%regionEvaluate(regionChild2)
       ! Update the running estimates of the integral and of its error, replacing the contribution of the parent region with
       ! those of its two children.
       integral =+integral              &
            &    -regionParent%integral &
            &    +regionChild1%integral &
            &    +regionChild2%integral
       error    =+error                 &
            &    -regionParent%error    &
            &    +regionChild1%error    &
            &    +regionChild2%error
       converged=self%hasConverged(integral,error)
       ! Grow the heap if necessary, then push both children onto it.
       if (countRegions+2 > sizeHeap) then
          call move_alloc(heap,heapTemporary)
          sizeHeap=2*sizeHeap
          allocate(heap(sizeHeap))
          heap(1:countRegions)=heapTemporary(1:countRegions)
          deallocate(heapTemporary)
       end if
       countRegions      =countRegions+1
       heap(countRegions)=regionChild1
       call heapSiftUp  (heap,countRegions)
       countRegions      =countRegions+1
       heap(countRegions)=regionChild2
       call heapSiftUp  (heap,countRegions)
    end do
    return
  end function genzMalikNDEvaluate

  subroutine genzMalikNDRegionEvaluate(self,region_)
    !!{RST
    Apply the degree-7 rule, and the embedded degree-5 rule, over the given region. The estimate of the integral is that of the
    degree-7 rule, the estimate of the error is the difference between the two rules, and the coordinate along which the region
    is to be split if refined is that of largest fourth difference.
    !!}
    implicit none
    class           (integratorGenzMalikND), intent(inout)                   :: self
    type            (region               ), intent(inout)                   :: region_
    double precision                       , dimension(self%countDimensions) :: center        , halfWidth       , &
         &                                                                      x             , differenceFourth
    double precision                                                         :: sum7          , sum5            , &
         &                                                                      volume        , integrand0      , &
         &                                                                      integrandPlus2, integrandMinus2 , &
         &                                                                      integrandPlus3, integrandMinus3 , &
         &                                                                      sumLambda2    , sumLambda3      , &
         &                                                                      sumLambda4    , sumLambda5      , &
         &                                                                      lambda2       , lambda3         , &
         &                                                                      lambda4       , lambda5
    integer                                                                  :: i             , j               , &
         &                                                                      k             , signs

    lambda2  =sqrt(lambda2Squared)
    lambda3  =sqrt(lambda3Squared)
    lambda4  =sqrt(lambda4Squared)
    lambda5  =sqrt(lambda5Squared)
    center   =0.5d0*(region_%upper+region_%lower)
    halfWidth=0.5d0*(region_%upper-region_%lower)
    volume   =product(region_%upper-region_%lower)
    ! The center point.
    integrand0      =self%integrandEvaluate(center)
    ! Points with a single non-zero coordinate, at ±λ₂ and ±λ₃. These also provide the fourth differences used to select the
    ! coordinate along which the region is to be split.
    sumLambda2      =0.0d0
    sumLambda3      =0.0d0
    differenceFourth=0.0d0
    do i=1,self%countDimensions
       x                  =center
       x               (i)=center(i)+halfWidth(i)*lambda2
       integrandPlus2     =self%integrandEvaluate(x)
       x               (i)=center(i)-halfWidth(i)*lambda2
       integrandMinus2    =self%integrandEvaluate(x)
       x               (i)=center(i)+halfWidth(i)*lambda3
       integrandPlus3     =self%integrandEvaluate(x)
       x               (i)=center(i)-halfWidth(i)*lambda3
       integrandMinus3    =self%integrandEvaluate(x)
       sumLambda2         =+sumLambda2                     &
            &              +integrandPlus2+integrandMinus2
       sumLambda3         =+sumLambda3                     &
            &              +integrandPlus3+integrandMinus3
       differenceFourth(i)=+abs(                                                                      &
            &                   +                    integrandPlus3+integrandMinus3-2.0d0*integrand0  &
            &                   -ratioLambdaSquared*(integrandPlus2+integrandMinus2-2.0d0*integrand0) &
            &                  )
    end do
    ! Points with two non-zero coordinates, at ±λ₄ in each.
    sumLambda4=0.0d0
    do i=1,self%countDimensions-1
       do j=i+1,self%countDimensions
          x         =center
          x      (i)=center(i)+halfWidth(i)*lambda4
          x      (j)=center(j)+halfWidth(j)*lambda4
          sumLambda4=sumLambda4+self%integrandEvaluate(x)
          x      (j)=center(j)-halfWidth(j)*lambda4
          sumLambda4=sumLambda4+self%integrandEvaluate(x)
          x      (i)=center(i)-halfWidth(i)*lambda4
          sumLambda4=sumLambda4+self%integrandEvaluate(x)
          x      (j)=center(j)+halfWidth(j)*lambda4
          sumLambda4=sumLambda4+self%integrandEvaluate(x)
       end do
    end do
    ! The 2ᵈ points with every coordinate at ±λ₅. The sign of each coordinate is taken from the bits of a counter.
    sumLambda5=0.0d0
    do signs=0,2**self%countDimensions-1
       do k=1,self%countDimensions
          if (btest(signs,k-1)) then
             x(k)=center(k)+halfWidth(k)*lambda5
          else
             x(k)=center(k)-halfWidth(k)*lambda5
          end if
       end do
       sumLambda5=sumLambda5+self%integrandEvaluate(x)
    end do
    ! Form the two rules.
    sum7   =+self%weight7Center *integrand0 &
         &  +self%weight7Lambda2*sumLambda2 &
         &  +self%weight7Lambda3*sumLambda3 &
         &  +self%weight7Lambda4*sumLambda4 &
         &  +self%weight7Lambda5*sumLambda5
    sum5   =+self%weight5Center *integrand0 &
         &  +self%weight5Lambda2*sumLambda2 &
         &  +self%weight5Lambda3*sumLambda3 &
         &  +self%weight5Lambda4*sumLambda4
    region_%integral=volume*    sum7
    region_%error   =volume*abs(sum7-sum5)
    ! Select the coordinate along which the region is to be split if it is refined: that of largest fourth difference, with ties
    ! broken in favor of the longest edge. Tie-breaking matters: for a discontinuous integrand the fourth differences take only a
    ! few discrete values and so tie frequently, and always resolving a tie the same way would split the region along a single
    ! coordinate forever, slivering it instead of refining it.
    region_%dimensionSplit=1
    do i=2,self%countDimensions
       if     (                                                                      &
            &        differenceFourth(i) >  differenceFourth(region_%dimensionSplit) &
            &  .or.                                                                  &
            &   (                                                                    &
            &        differenceFourth(i) == differenceFourth(region_%dimensionSplit) &
            &    .and.                                                               &
            &        halfWidth       (i) >  halfWidth       (region_%dimensionSplit) &
            &   )                                                                    &
            & ) region_%dimensionSplit=i
    end do
    return
  end subroutine genzMalikNDRegionEvaluate

  double precision function genzMalikNDIntegrandEvaluate(self,x)
    !!{RST
    Evaluate the integrand, counting the evaluation.
    !!}
    implicit none
    class           (integratorGenzMalikND), intent(inout)               :: self
    double precision                       , intent(in   ), dimension(:) :: x

    self%countEvaluations       =self%countEvaluations+1
    genzMalikNDIntegrandEvaluate=self%integrand(x)
    return
  end function genzMalikNDIntegrandEvaluate

  subroutine heapSiftUp(heap,i)
    !!{RST
    Restore the heap property by sifting the element at position ``i`` towards the root.
    !!}
    implicit none
    type   (region), intent(inout), dimension(:) :: heap
    integer        , intent(in   )               :: i
    type   (region)                              :: swap
    integer                                      :: child, parent

    child=i
    do while (child > 1)
       parent=child/2
       if (heap(parent)%error >= heap(child)%error) exit
       swap        =heap(parent)
       heap(parent)=heap(child )
       heap(child )=swap
       child       =parent
    end do
    return
  end subroutine heapSiftUp

  subroutine heapSiftDown(heap,countRegions,i)
    !!{RST
    Restore the heap property by sifting the element at position ``i`` away from the root.
    !!}
    implicit none
    type   (region), intent(inout), dimension(:) :: heap
    integer        , intent(in   )               :: countRegions, i
    type   (region)                              :: swap
    integer                                      :: parent      , child, &
         &                                          largest

    parent=i
    do
       child  =2*parent
       largest=  parent
       if (child   <= countRegions .and. heap(child  )%error > heap(largest)%error) largest=child
       if (child+1 <= countRegions .and. heap(child+1)%error > heap(largest)%error) largest=child+1
       if (largest == parent) exit
       swap         =heap(parent )
       heap(parent )=heap(largest)
       heap(largest)=swap
       parent       =largest
    end do
    return
  end subroutine heapSiftDown

end module Numerical_Integration_nD

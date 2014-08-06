!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements Poisson binomial distributions.

module Math_Distributions_Poisson_Binomial
  !% Implements Poisson binomial distributions.
  private
  public :: Poisson_Binomial_Distribution           , Poisson_Binomial_Distribution_Mean               , &
       &    Poisson_Binomial_Distribution_Mean_Pairs, Poisson_Binomial_Distribution_Mean_Pairs_Jacobian

contains

  double precision function Poisson_Binomial_Distribution(k,p)
    !% Computes the Poisson binomial distribution with event probabilities {\tt p} at argument
    !% {\tt k}. Uses the discrete Fourier transform method proposed by \cite{fernandez_closed-form_2010}.
    use Numerical_Constants_Math
    implicit none
    integer         , intent(in   )               :: k
    double precision, intent(in   ), dimension(:) :: p
    double complex                                :: Cl, product, probability
    integer                                       :: l , m

    Poisson_Binomial_Distribution=0.0d0
    if (k < 0 .or. k > size(p)) return
    probability=dcmplx(0.0d0,0.0d0)
    do l=0,size(p)
       product=dcmplx(1.0d0,0.0d0)
       Cl     =exp(dble(l)*2.0d0*Pi*dcmplx(0.0d0,1.0d0)/dble(1+size(p)))
       do m=1,size(p)
          product=product*(1.0d0+(Cl-1.0d0)*p(m))
       end do
       probability=probability+product*exp(-dble(l*k)*2.0d0*Pi*dcmplx(0.0d0,1.0d0)/dble(1+size(p)))
    end do
    Poisson_Binomial_Distribution=max(dreal(probability)/dble(1+size(p)),0.0d0)
    return
  end function Poisson_Binomial_Distribution

  function Poisson_Binomial_Distribution_Jacobian(k,p)
    !% Computes the Jacobian of the Poisson binomial distribution with event probabilities {\tt p} at argument
    !% {\tt k}. Uses the discrete Fourier transform method proposed by \cite{fernandez_closed-form_2010}.
    use Numerical_Constants_Math
    implicit none
    integer         , intent(in   )                     :: k
    double precision, intent(in   ), dimension(     : ) :: p
    double precision               , dimension(size(p)) :: Poisson_Binomial_Distribution_Jacobian
    double complex                 , dimension(size(p)) :: probability
    double complex                                      :: Cl, product
    integer                                             :: l , m

    Poisson_Binomial_Distribution_Jacobian=0.0d0
    if (k < 0 .or. k > size(p)) return
    probability=dcmplx(0.0d0,0.0d0)
    do l=0,size(p)
       product=dcmplx(1.0d0,0.0d0)
       Cl     =exp(dble(l)*2.0d0*Pi*dcmplx(0.0d0,1.0d0)/dble(1+size(p)))
       do m=1,size(p)
          product=product*(1.0d0+(Cl-1.0d0)*p(m))
       end do
       probability=probability+product*exp(-dble(l*k)*2.0d0*Pi*dcmplx(0.0d0,1.0d0)/dble(1+size(p)))*(Cl-1.0d0)/(1.0d0+(Cl-1.0d0)*p)
    end do
    Poisson_Binomial_Distribution_Jacobian=max(dreal(probability)/dble(1+size(p)),0.0d0)
    return
  end function Poisson_Binomial_Distribution_Jacobian

  double precision function Poisson_Binomial_Distribution_Mean(p)
    !% Computes the mean of a Poisson binomial distribution.
    implicit none
    double precision, intent(in   ), dimension(:) :: p

    Poisson_Binomial_Distribution_Mean=sum(p)
    return
  end function Poisson_Binomial_Distribution_Mean

  double precision function Poisson_Binomial_Distribution_Mean_Pairs(p)
    !% Computes the mean number of pairs expected from a Poisson binomial distribution with
    !% event probabilities {\tt p}. Assumes that pair order is significant, i.e. both $AB$ and
    !% $BA$ are counted.
    implicit none
    double precision, intent(in   ), dimension(:) :: p
    integer                                       :: k

    Poisson_Binomial_Distribution_Mean_Pairs=0.0d0
    if (size(p) < 2) return
    do k=2,size(p)
       Poisson_Binomial_Distribution_Mean_Pairs=Poisson_Binomial_Distribution_Mean_Pairs+dble(k*(k-1))*Poisson_Binomial_Distribution(k,p)
    end do
    return
  end function Poisson_Binomial_Distribution_Mean_Pairs

  function Poisson_Binomial_Distribution_Mean_Pairs_Jacobian(p)
    !% Computes the Jacobian of the mean number of pairs expected from a Poisson binomial distribution with
    !% event probabilities {\tt p}. Assumes that pair order is significant, i.e. both $AB$ and
    !% $BA$ are counted.
    implicit none
    double precision, intent(in   ), dimension(     : ) :: p
    double precision               , dimension(size(p)) :: Poisson_Binomial_Distribution_Mean_Pairs_Jacobian
    integer                                             :: k

    Poisson_Binomial_Distribution_Mean_Pairs_Jacobian=0.0d0
    if (size(p) < 2) return
    do k=2,size(p)
       Poisson_Binomial_Distribution_Mean_Pairs_Jacobian=+              Poisson_Binomial_Distribution_Mean_Pairs_Jacobian      &
            &                                            +dble(k*(k-1))*Poisson_Binomial_Distribution_Jacobian           (k,p)
    end do
    return
  end function Poisson_Binomial_Distribution_Mean_Pairs_Jacobian

end module Math_Distributions_Poisson_Binomial

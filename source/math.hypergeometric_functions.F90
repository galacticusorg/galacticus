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

!: ./work/build/pFq/pfq.new.o

!% Contains a module which implements hypergeometric functions.

module Hypergeometric_Functions
  !% Implements hypergeometric functions.
  use FGSL
  implicit none
  private
  public :: Hypergeometric_1F1, Hypergeometric_2F1, Hypergeometric_pFq

  interface Hypergeometric_pFq
     module procedure :: Hypergeometric_pFq_Real
     module procedure :: Hypergeometric_pFq_Complex
  end interface Hypergeometric_pFq
  
contains

  double precision function Hypergeometric_1F1(a,b,x)
    !% Evaluate the $_1F_1(a_1;b_1;x)$ hypergeometric function.
    implicit none
    double precision, intent(in   ) :: a(1), b(1), x

    Hypergeometric_1F1=FGSL_SF_Hyperg_1F1(a(1),b(1),x)
    return
  end function Hypergeometric_1F1

  double precision function Hypergeometric_2F1(a,b,x)
    !% Evaluate the $_2F_1(a_1,a_2;b_1;x)$ hypergeometric function.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: a(2), b(1), x

    ! GSL only evaluates this function for |x|<1.
    if (abs(x) <= 1.0d0) then
       ! |x|<1 so simply call the GSL function to compute the function.
       Hypergeometric_2F1=FGSL_SF_Hyperg_2F1(a(1),a(2),b(1),x)
    else if (x < -1.0d0) then
       ! x<-1 so use a Pfaff transformation to evaluate in terms of a hypergeometric function with |x|<1.
       Hypergeometric_2F1=FGSL_SF_Hyperg_2F1(a(2),b(1)-a(1),b(1),x/(x-1.0d0))/(1.0d0-x)**a(2)
    else
       call Galacticus_Error_Report('Hypergeometric_2F1','function cannot be evaluated for x>1')
    end if
    return
  end function Hypergeometric_2F1

  double complex function Hypergeometric_pFq_Complex(a,b,x)
    !% Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x), using the algorithm of
    !% \cite{perger_numerical_1993}.
    implicit none
    double complex, intent(in   ), dimension(:) :: a    , b
    double complex, intent(in   )               :: x
    double complex                              :: PFQ
    integer                                     :: LNPFQ, IX, NSIGFIG

    LNPFQ  = 0
    IX     = 0
    NSIGFIG=10
    Hypergeometric_pFq_Complex=PFQ(a,b,size(a),size(b),x,LNPFQ,IX,NSIGFIG)
    return
  end function Hypergeometric_pFq_Complex

  double precision function Hypergeometric_pFq_Real(a,b,x)
    !% Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x) for real arguments.
    implicit none
    double precision, intent(in   ), dimension(:) :: a  , b
    double precision, intent(in   )               :: x
    
    Hypergeometric_pFq_Real=real(Hypergeometric_pFq_Complex(dcmplx(a),dcmplx(b),dcmplx(x)))
    return
  end function Hypergeometric_pFq_Real

end module Hypergeometric_Functions

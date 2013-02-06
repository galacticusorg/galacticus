!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program to test mathematical special functions.

program Test_Math_Special_Functions
  !% Tests of mathematical special functions.
  use Unit_Tests
  use Bessel_Functions
  use Exponential_Integrals
  use Factorials
  use Gamma_Functions
  use Hypergeometric_Functions
  implicit none
  double precision, dimension(10) :: argument=[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0]
  double precision, dimension(10) :: P=[       &  
       &0.1353352832d0,&
       &                         0.4060058496d0,&
       &                         0.6766764160d0,&
       &                         0.8571234602d0,&
       &                         0.9473469824d0,&
       &                         0.9834363913d0,&
       &                         0.9954661943d0,&
       &                         0.9989032808d0,&
       &                         0.9997625524d0,&
       &                         0.9999535016d0]
  double precision, dimension(10) :: Q=[ 0.8646647168d0,&
       &                         0.5939941504d0,&
       &                         0.3233235840d0,&
       &                         0.1428765398d0,&
       &                         0.0526530176d0,&
       &                         0.0165636087d0,&
       &                         0.0045338057d0,&
       &                         0.0010967192d0,&
       &                         0.0002374473d0,&
       &                         0.0000464981d0]
  double precision, dimension(10) :: BesselK0,BesselK1,BesselI0,BesselI1,sineIntegral,cosineIntegral,factorials,doubleFactorial &
       &,gammaFunction,logGammaFunction,incompleteGammaFunction,incompleteComplementaryGammaFunction&
       &,inverseGammaFunctionIncompleteComplementary,inverseGammaFunctionIncomplete,hypergeometric1F1,hypergeometric2F1
  integer                         :: i

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Math: special functions")

  ! Evaluate functions.
  do i=1,10
     BesselK0                                   (i)=Bessel_Function_K0                             (                             argument(i)       )
     BesselK1                                   (i)=Bessel_Function_K1                             (                             argument(i)       )
     BesselI0                                   (i)=Bessel_Function_I0                             (                             argument(i)       )
     BesselI1                                   (i)=Bessel_Function_I1                             (                             argument(i)       )
     sineIntegral                               (i)=Sine_Integral                                  (                             argument(i)       )
     cosineIntegral                             (i)=Cosine_Integral                                (                             argument(i)       )
     factorials                                 (i)=Factorial                                      (                                      i        )
     doubleFactorial                            (i)=Logarithmic_Double_Factorial                   (                                      i        )
     gammaFunction                              (i)=Gamma_Function                                 (                             argument(i)       ) 
     logGammaFunction                           (i)=Gamma_Function_Logarithmic                     (                             argument(i)       )
     incompleteGammaFunction                    (i)=Gamma_Function_Incomplete                      (                             argument(i),2.0d0 )
     incompleteComplementaryGammaFunction       (i)=Gamma_Function_Incomplete_Complementary        (                             argument(i),2.0d0 ) 
     inverseGammaFunctionIncomplete             (i)=Inverse_Gamma_Function_Incomplete              (                             argument(i),P(i)  )
     inverseGammaFunctionIncompleteComplementary(i)=Inverse_Gamma_Function_Incomplete_Complementary(                             argument(i),Q(i)  )
     hypergeometric1F1                          (i)=Hypergeometric_1F1                             ([1.0d0      ],[2.0d0],       argument(i)       )
     hypergeometric2F1                          (i)=Hypergeometric_2F1                             ([1.5d0,0.5d0],[1.5d0],1.0d0/(argument(i)+1.0d0))
  end do

  ! Test Bessel function results.
  call Assert("Bessel K₀(x)",       &
       &       BesselK0,            &
       &       [                    &
       &        0.4210244382d0    , &
       &        0.1138938727d0    , &
       &        0.03473950439d0   , &
       &        0.01115967609d0   , &
       &        0.003691098334d0  , &
       &        0.001243994328d0  , &
       &        0.0004247957419d0 , &
       &        0.0001464707052d0 , &
       &        0.00005088131296d0, &
       &        0.00001778006232d0  &
       &       ],                   &
       &       relTol=1.0d-6        &
       &     )
  call Assert("Bessel K₁(x)",       &
       &       BesselK1,            &
       &       [                    &
       &        0.6019072302d0,     &
       &        0.1398658818d0,     &
       &        0.04015643113d0,    &
       &        0.01248349889d0,    &
       &        0.004044613445d0,   &
       &        0.001343919718d0,   &
       &        0.0004541824869d0,  &
       &        0.0001553692118d0,  &
       &        0.00005363701638d0, &
       &        0.00001864877345d0  &
       &       ],                   &
       &       relTol=1.0d-6        &
       &     )
  call Assert("Bessel I₀(x)",       &
       &       BesselI0,            &
       &       [                    &
       &           1.266065878d0,   &
       &           2.279585302d0,   &
       &           4.880792586d0,   &
       &          11.30192195d0,    &
       &          27.23987182d0,    &
       &          67.23440698d0,    &
       &         168.5939085d0,     &
       &         427.5641157d0,     &
       &        1093.588355d0,      &
       &        2815.716628d0       &
       &       ],                   &
       &       relTol=1.0d-6        &
       &     )
  call Assert("Bessel I₁(x)",       &
       &       BesselI1,            &
       &       [                    &
       &           0.5651591040d0,  &
       &           1.590636855d0,   &
       &           3.953370217d0,   &
       &           9.759465154d0,   &
       &          24.33564214d0,    &
       &          61.34193678d0,    &
       &         156.0390929d0,     &
       &         399.8731368d0,     &
       &        1030.914723d0,      &
       &        2670.988304d0       &
       &       ],                   &
       &       relTol=1.0d-6        &
       &     )

  ! Test exponential integrals.
  call Assert("sine integral, Si(x)",   &
       &       sineIntegral,            &
       &       [                        &
       &        0.9460830704d0,         &
       &        1.605412977d0,          &
       &        1.848652528d0,          &
       &        1.758203139d0,          &
       &        1.549931245d0,          &
       &        1.424687551d0,          &
       &        1.454596614d0,          &
       &        1.574186822d0,          &
       &        1.665040076d0,          &
       &        1.658347594d0           &
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )
  call Assert("cosine integral, Ci(x)", &
       &       cosineIntegral,          &
       &       [                        &
       &         0.3374039229d0,        &
       &         0.4229808288d0,        &
       &         0.1196297860d0,        &
       &        -0.1409816979d0,        &
       &        -0.1900297497d0,        &
       &        -0.06805724389d0,       &
       &         0.07669527848d0,       &
       &         0.1224338825d0,        &
       &         0.05534753133d0,       &
       &        -0.04545643300d0        &
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )

  ! Test factorials.
  call Assert("factorial, x!",          &
       &       factorials,              &
       &       [                        &
       &            1.0d0,              &
       &            2.0d0,              &
       &            6.0d0,              &
       &           24.0d0,              &
       &          120.0d0,              &
       &          720.0d0,              &
       &         5040.0d0,              &
       &        40320.0d0,              &
       &            3.62880d5,          &
       &            3.628800d6          &
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )
  call Assert("double factorial, x!!",  &
       &       doubleFactorial,         &
       &       [                        &
       &        0.0d0,                  &
       &        0.6931471806d0,         &
       &        1.098612289d0,          &
       &        2.079441542d0,          &
       &        2.708050201d0,          &
       &        3.871201011d0,          &
       &        4.653960350d0,          &
       &        5.950642553d0,          &
       &        6.851184927d0,          &
       &        8.253227646d0           &
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )

  ! Test gamma functions.
  call Assert("gamma function, Γ(x)",   &
       &       gammaFunction,           &
       &       [                        &
       &            1.0d0,              &
       &            1.0d0,              &
       &            2.0d0,              &
       &            6.0d0,              &
       &           24.0d0,              &
       &          120.0d0,              &
       &          720.0d0,              &
       &         5040.0d0,              &
       &        40320.0d0,              &
       &            3.62880d5           &
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )
  call Assert("log of gamma functions, ln(Γ(x))", &
       &       logGammaFunction,                  &
       &       [                                  &
       &         0.0d0,                           &
       &         0.0d0,                           &
       &         0.6931471806d0,                  &
       &         1.791759469d0,                   &
       &         3.178053830d0,                   &
       &         4.787491743d0,                   &
       &         6.579251212d0,                   &
       &         8.525161361d0,                   &
       &        10.60460290d0,                    &
       &        12.80182748d0                     &
       &       ],                                 &
       &       relTol=1.0d-6                      &
       &     )
  call Assert("incomplete gamma function, Γ(x,2)", &
       &       incompleteGammaFunction,            &
       &       [                                   &
       &        0.1353352832d0,                    &
       &        0.4060058496d0,&
       &        0.6766764160d0,&
       &        0.8571234602d0,&
       &        0.9473469824d0,&
       &        0.9834363913d0,&
       &        0.9954661943d0,&
       &        0.9989032808d0,&
       &        0.9997625524d0,&
       &        0.9999535016d0&
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )
  call Assert("complementary gamma function, 1-Γ(x,2)",   &
       &       incompleteComplementaryGammaFunction,           &
       &       [                        &
       &                         0.8646647168d0,&
       &                         0.5939941504d0,&
       &                         0.3233235840d0,&
       &                         0.1428765398d0,&
       &                         0.0526530176d0,&
       &                         0.0165636087d0,&
       &                         0.0045338057d0,&
       &                         0.0010967192d0,&
       &                         0.0002374473d0,&
       &                         0.0000464981d0 &
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )
  call Assert("inverse incomplete gamma function, Γ⁻¹(x,2)", &
       &       inverseGammaFunctionIncomplete,               &
       &       [                                             &
       &        2.0d0,                                       &
       &        2.0d0,                                       &
       &        2.0d0,                                       &
       &        2.0d0,                                       &
       &        2.0d0,                                       &
       &        2.0d0,                                       &
       &        2.0d0,                                       &
       &        2.0d0,                                       &
       &        2.0d0,                                       &
       &        2.0d0                                        &
       &       ],                                            &
       &       relTol=1.0d-6                                 &
       &     )
  call Assert("inverse complementary gamma function, {1-Γ(x,2)}⁻¹", &
       &       inverseGammaFunctionIncompleteComplementary,         &
       &       [                                                    &
       &        2.0d0,                                              &
       &        2.0d0,                                              &
       &        2.0d0,                                              &
       &        2.0d0,                                              &
       &        2.0d0,                                              &
       &        2.0d0,                                              &
       &        2.0d0,                                              &
       &        2.0d0,                                              &
       &        2.0d0,                                              &
       &        2.0d0                                               &
       &       ],                                                   &
       &       relTol=1.0d-6                                        &
       &     )

  ! Test hypergeometric functions.
  call Assert("hypergeometric, ₁F₁([1],[2],x)", &
       &       hypergeometric1F1,               &
       &       [                                &
       &           1.718281828d0,               &
       &           3.194528049d0,               &
       &           6.361845641d0,               &
       &          13.39953751d0,                &
       &          29.48263182d0,                &
       &          67.07146558d0,                &
       &         156.5190226d0,                 &
       &         372.4947484d0,                 &
       &         900.2315475d0,                 &
       &        2202.546579d0                   &
       &       ],                               &
       &       relTol=1.0d-6                    &
       &     )
  call Assert("hypergeometric, ₂F₁([3/2,1/2],[3/2],x)", &
       &       hypergeometric2F1,                       &
       &       [                                        &
       &        1.414213562d0,                          &
       &        1.224744871d0,                          &
       &        1.154700538d0,                          &
       &        1.118033989d0,                          &
       &        1.095445115d0,                          &
       &        1.080123450d0,                          &
       &        1.069044968d0,                          &
       &        1.060660172d0,                          &
       &        1.054092553d0,                          &
       &        1.048808848d0                           &
       &       ],                                       &
       &       relTol=1.0d-6                            &
       &     )

  ! Test hypergeometric 2F1 function transitions for |x|>1.
  hypergeometric2F1(1)=Hypergeometric_2F1([1.0d0,1.0d0],[2.0d0],-2.0d0)
  call Assert("hypergeometric, ₂F₁([1,1],[2],-2)",hypergeometric2F1(1),log(3.0d0)/2.0d0,relTol=1.0d-6 )

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Math_Special_Functions

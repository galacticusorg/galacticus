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
Contains a program to test mathematical special functions.
!!}

program Test_Math_Special_Functions
  !!{
  Tests of mathematical special functions.
  !!}
  use :: Bessel_Functions        , only : Bessel_Function_I0               , Bessel_Function_I1                             , Bessel_Function_J0                     , Bessel_Function_J0_Zero       , &
          &                               Bessel_Function_J1               , Bessel_Function_J1_Zero                        , Bessel_Function_Jn                     , Bessel_Function_Jn_Zero       , &
          &                               Bessel_Function_K0               , Bessel_Function_K1                             , Bessel_Function_In
  use :: Binomial_Coefficients   , only : Binomial_Coefficient
  use :: Dilogarithms            , only : Dilogarithm
  use :: Display                 , only : displayVerbositySet              , verbosityLevelStandard
  use :: Error_Functions         , only : Error_Function
  use :: Exponential_Integrals   , only : Cosine_Integral                  , Sine_Integral
  use :: Factorials              , only : Factorial                        , Logarithmic_Double_Factorial
  use :: Beta_Functions          , only : Beta_Function                    , Beta_Function_Incomplete_Normalized
  use :: Gamma_Functions         , only : Gamma_Function                   , Gamma_Function_Incomplete                      , Gamma_Function_Incomplete_Complementary, Gamma_Function_Logarithmic    , &
          &                               Inverse_Gamma_Function_Incomplete, Inverse_Gamma_Function_Incomplete_Complementary, Digamma_Function
  use :: Hypergeometric_Functions, only : Hypergeometric_1F1               , Hypergeometric_2F1                             , Hypergeometric_pFq                     , Hypergeometric_pFq_Regularized
  use :: Polylogarithms          , only : Polylogarithm_2                  , Polylogarithm_3
  use :: Numerical_Constants_Math, only : Pi
  use :: Error                   , only : Error_Handler_Register
  use :: Unit_Tests              , only : Assert                           , Unit_Tests_Begin_Group                         , Unit_Tests_End_Group                   , Unit_Tests_Finish
  implicit none
  double precision, dimension(10) :: argument                            =[1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0]
  double precision, dimension(10) :: P                                   =[0.1353352832d0,0.4060058496d0,0.6766764160d0,0.8571234602d0,0.9473469824d0,0.9834363913d0,0.9954661943d0,0.9989032808d0,0.9997625524d0,0.9999535016d0]
  double precision, dimension(10) :: Q                                   =[0.8646647168d0,0.5939941504d0,0.3233235840d0,0.1428765398d0,0.0526530176d0,0.0165636087d0,0.0045338057d0,0.0010967192d0,0.0002374473d0,0.0000464981d0]
  double precision, dimension(10) :: BesselI0                                                                                                                                                                                    , BesselI1                                   , &
       &                             BesselK0                                                                                                                                                                                    , BesselK1                                   , &
       &                             BesselJ0                                                                                                                                                                                    , BesselJ1                                   , &
       &                             BesselJ2                                                                                                                                                                                    , BesselJ0Zero                               , &
       &                             BesselJ1Zero                                                                                                                                                                                , BesselJ2Zero                               , &
       &                             cosineIntegral                                                                                                                                                                              , doubleFactorial                            , &
       &                             factorials                                                                                                                                                                                  , gammaFunction                              , &
       &                             hypergeometric1F1                                                                                                                                                                           , hypergeometric2F1                          , &
       &                             incompleteComplementaryGammaFunction                                                                                                                                                        , incompleteGammaFunction                    , &
       &                             incompleteComplementaryGammaFunction2                                                                                                                                                       , incompleteGammaFunction2                   , &
       &                             inverseGammaFunctionIncomplete                                                                                                                                                              , inverseGammaFunctionIncompleteComplementary, &
       &                             logGammaFunction                                                                                                                                                                            , sineIntegral                               , &
       &                             hypergeometric1F2                                                                                                                                                                           , hypergeometric2F1approx                    , &
       &                             hypergeometric2F1NegativeArgument                                                                                                                                                           , hypergeometric2F1approxNegativeArgument    , &
       &                             hypergeometric3F2NegativeArgument                                                                                                                                                           , hypergeometric3F2Accelerated               , &
       &                             polylogarithm2                                                                                                                                                                              , polylogarithm3                             , &
       &                             hypergeometric1F2Regularized                                                                                                                                                                , BesselI2                                   , &
       &                             BesselIHalf                                                                                                                                                                                 , digammaFunction                            , &
       &                             dilogarithm_
  double complex  , dimension(17) :: errorFunctionComplex
  integer                         :: i


  ! Establish error handlers.
  call Error_Handler_Register()
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Math: special functions")

  ! Evaluate functions.
  do i=1,10
     BesselK0                                   (i)=Bessel_Function_K0                             (                                          argument(i)                                      )
     BesselK1                                   (i)=Bessel_Function_K1                             (                                          argument(i)                                      )
     BesselJ0                                   (i)=Bessel_Function_J0                             (                                          argument(i)                                      )
     BesselJ1                                   (i)=Bessel_Function_J1                             (                                          argument(i)                                      )
     BesselJ2                                   (i)=Bessel_Function_Jn                             (2                                ,        argument(i)                                      )
     BesselJ0Zero                               (i)=Bessel_Function_J0_Zero                        (                                                   i                                       )
     BesselJ1Zero                               (i)=Bessel_Function_J1_Zero                        (                                                   i                                       )
     BesselJ2Zero                               (i)=Bessel_Function_Jn_Zero                        (2.0d0                            ,                 i                                       )
     BesselI0                                   (i)=Bessel_Function_I0                             (                                          argument(i)                                      )
     BesselI1                                   (i)=Bessel_Function_I1                             (                                          argument(i)                                      )
     BesselI2                                   (i)=Bessel_Function_In                             (2                                ,        argument(i)                                      )
     BesselIHalf                                (i)=Bessel_Function_In                             (0.5d0                            ,        argument(i)                                      )
     sineIntegral                               (i)=Sine_Integral                                  (                                          argument(i)                                      )
     cosineIntegral                             (i)=Cosine_Integral                                (                                          argument(i)                                      )
     factorials                                 (i)=Factorial                                      (                                                   i                                       )
     doubleFactorial                            (i)=Logarithmic_Double_Factorial                   (                                                   i                                       )
     gammaFunction                              (i)=Gamma_Function                                 (                                          argument(i)                                      )
     digammaFunction                            (i)=Digamma_Function                               (                                          argument(i)                                      )
     logGammaFunction                           (i)=Gamma_Function_Logarithmic                     (                                          argument(i)                                      )
     incompleteGammaFunction                    (i)=Gamma_Function_Incomplete                      (                                          argument(i), 2.0d0                               )
     incompleteComplementaryGammaFunction       (i)=Gamma_Function_Incomplete_Complementary        (                                          argument(i), 2.0d0                               )
     incompleteGammaFunction2                   (i)=Gamma_Function_Incomplete                      (                                          argument(i),16.67d0                              )
     incompleteComplementaryGammaFunction2      (i)=Gamma_Function_Incomplete_Complementary        (                                          argument(i),16.67d0                              )
     inverseGammaFunctionIncomplete             (i)=Inverse_Gamma_Function_Incomplete              (                                          argument(i),P(i)                                 )
     inverseGammaFunctionIncompleteComplementary(i)=Inverse_Gamma_Function_Incomplete_Complementary(                                          argument(i),Q(i)                                 )
     hypergeometric1F1                          (i)=Hypergeometric_1F1                             ([1.0d0            ],[2.0d0      ],        argument(i)                                      )
     hypergeometric2F1                          (i)=Hypergeometric_2F1                             ([1.5d0,0.5d0      ],[1.5d0      ],1.0d0/( argument(i)+1.0d0)                               )
     hypergeometric2F1approx                    (i)=Hypergeometric_2F1                             ([1.5d0,0.5d0      ],[1.5d0      ],1.0d0/( argument(i)+1.0d0)      ,toleranceRelative=1.0d-6)
     hypergeometric1F2                          (i)=Hypergeometric_pFq                             ([1.5d0            ],[1.5d0,0.5d0],        argument(i)                                      )
     hypergeometric2F1NegativeArgument          (i)=Hypergeometric_2F1                             ([1.5d0,0.7d0      ],[2.0d0      ],  -exp( argument(i)      )                               )
     hypergeometric2F1approxNegativeArgument    (i)=Hypergeometric_2F1                             ([1.5d0,0.7d0      ],[2.0d0      ],  -exp( argument(i)      )      ,toleranceRelative=1.0d-6)
     hypergeometric3F2NegativeArgument          (i)=Hypergeometric_pFq                             ([1.5d0,0.7d0,1.1d0],[2.0d0,0.2d0],  -exp( argument(i)      )                               )
     hypergeometric3F2Accelerated               (i)=Hypergeometric_pFq                             ([1.5d0,0.7d0,1.1d0],[2.0d0,0.2d0],      ( argument(i)-6.0d0)/5.0d0,useAcceleration=.true.  )
     hypergeometric1F2Regularized               (i)=Hypergeometric_pFq_Regularized                 ([1.5d0            ],[1.5d0,0.5d0],        argument(i)                                      )
     polylogarithm2                             (i)=Polylogarithm_2                                (                                  -1.0d0/ argument(i)                                      )
     polylogarithm3                             (i)=Polylogarithm_3                                (                                  -1.0d0/ argument(i)                                      )
     dilogarithm_                               (i)=Dilogarithm                                    (                                  -1.0d0/ argument(i)                                      )
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
  call Assert("Bessel J₀(x)",            &
       &       BesselJ0,                 &
       &       [                         &
       &        +0.765197686557966600d0, &
       &        +0.223890779141235700d0, &
       &        -0.260051954901933400d0, &
       &        -0.397149809863847400d0, &
       &        -0.177596771314338300d0, &
       &        +0.150645257250996900d0, &
       &        +0.300079270519555600d0, &
       &        +0.171650807137553900d0, &
       &        -0.090333611182876120d0, &
       &        -0.245935764451348000d0  &
       &       ],                        &
       &       relTol=1.0d-6             &
       &     )
  call Assert("Bessel J₁(x)",            &
       &       BesselJ1,                 &
       &       [                         &
       &        +0.440050585744933500d0, &
       &        +0.576724807756873400d0, &
       &        +0.339058958525936500d0, &
       &        -0.066043328023549140d0, &
       &        -0.327579137591465200d0, &
       &        -0.276683858127565600d0, & 
       &        -0.004682823482345833d0, &
       &        +0.234636346853914600d0, &
       &        +0.245311786573325300d0, &
       &        +0.04347274616886136d0   &
       &       ],                        &
       &       relTol=1.0d-6             &
       &     )
  call Assert("Bessel J₂(x)",            &
       &       BesselJ2,                 &
       &       [                         &
       &        +0.114903484931900500d0, &
       &        +0.352834028615637700d0, &
       &        +0.486091260585891100d0, &
       &        +0.364128145852072800d0, &
       &        +0.046565116277752220d0, &
       &        -0.242873209960185500d0, & 
       &        -0.301417220085940100d0, &
       &        -0.112991720424075200d0, &
       &        +0.144847341532504000d0, &
       &        +0.254630313685120300d0  &
       &       ],                        &
       &       relTol=1.0d-6             &
       &     )
  call Assert("Bessel J₀(xᵢ)=0",          &
       &       BesselJ0Zero,             &
       &       [                         &
       &         2.404825557695773d0,    &
       &         5.520078110286311d0,    &
       &         8.653727912911012d0,    & 
       &        11.791534439014280d0,    &
       &        14.930917708487790d0,    &
       &        18.071063967910920d0,    &
       &        21.211636629879260d0,    &
       &        24.352471530749300d0,    &
       &        27.493479132040250d0,    & 
       &        30.634606468431980d0     &
       &       ],                        &
       &       relTol=1.0d-6             &
       &     )
  call Assert("Bessel J₁(xᵢ)=0",          &
       &       BesselJ1Zero,             &
       &       [                         &
       &         3.831705970207515d0,    &
       &         7.015586669815619d0,    &
       &        10.173468135062720d0,    &
       &        13.323691936314220d0,    &
       &        16.470630050877630d0,    &
       &        19.615858510468240d0,    &
       &        22.760084380592770d0,    &
       &        25.903672087618380d0,    &
       &        29.046828534916860d0,    & 
       &        32.189679910974400d0     &
       &       ],                        &
       &       relTol=1.0d-6             &
       &     )
  call Assert("Bessel J₂(xᵢ)=0",          &
       &       BesselJ2Zero,             &
       &       [                         &
       &         5.135622301840683d0,    &
       &         8.417244140399857d0,    &
       &        11.619841172149060d0,    &
       &        14.795951782351260d0,    &
       &        17.959819494987830d0,    &
       &        21.116997053021850d0,    &
       &        24.270112313573100d0,    &
       &        27.420573549984560d0,    &
       &        30.569204495516400d0,    & 
       &        33.716519509222700d0     &
       &       ],                        &
       &       relTol=1.0d-6             &
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
  call Assert("Bessel I₂(x)",          &
       &       BesselI2,               &
       &       [                       &
       &        +1.357476697670383d-1, &
       &        +6.889484476987382d-1, &
       &        +2.245212440929951d+0, &
       &        +6.422189375284105d+0, &
       &        +1.750561496662423d+1, &
       &        +4.678709471726457d+1, &
       &        +1.240113105474453d+2, &
       &        +3.275958315261646d+2, &
       &        +8.644961939520510d+2, &
       &        +2.281518967726002d+3  &
       &       ],                      &
       &       relTol=1.0d-6           &
       &     )
  call Assert("Bessel I½(x)",          &
       &       BesselIHalf,            &
       &       [                       &
       &        +9.376748882454880d-1, &
       &        +2.046236863089055d+0, &
       &        +4.614822903407601d+0, &
       &        +1.088710179858842d+1, &
       &        +2.647754749755906d+1, &
       &        +6.570503691665827d+1, &
       &        +1.653567995485437d+2, &
       &        +4.204563140044776d+2, &
       &        +1.077554243705911d+3, &
       &        +2.778784603874571d+3  &
       &       ],                      &
       &       relTol=1.0d-6           &
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

  ! Test Gamma functions.
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
  call Assert("digamma function, ψ(x)", &
       &       digammaFunction,         &
       &       [                        &
       &        -5.772156649015328d-1 , &
       &        +4.227843350984672d-1 , &
       &        +9.227843350984670d-1 , &
       &        +1.256117668431800d+0 , &
       &        +1.506117668431801d+0 , &
       &        +1.706117668431800d+0 , &
       &        +1.872784335098467d+0 , &
       &        +2.015641477955610d+0 , &
       &        +2.140641477955610d+0 , &
       &        +2.251752589066721d+0   &
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
       &        0.1353352832d0                   , &
       &        0.4060058496d0                   , &
       &        0.6766764160d0                   , &
       &        0.8571234602d0                   , &
       &        0.9473469824d0                   , &
       &        0.9834363913d0                   , &
       &        0.9954661943d0                   , &
       &        0.9989032808d0                   , &
       &        0.9997625524d0                   , &
       &        0.9999535016d0                     &
       &       ],                                  &
       &       relTol=1.0d-6                       &
       &     )
  call Assert("complementary gamma function, 1-Γ(x,2)", &
       &       incompleteComplementaryGammaFunction   , &
       &       [                                        &
       &        0.8646647168d0                        , &
       &        0.5939941504d0                        , &
       &        0.3233235840d0                        , &
       &        0.1428765398d0                        , &
       &        0.0526530176d0                        , &
       &        0.0165636087d0                        , &
       &        0.0045338057d0                        , &
       &        0.0010967192d0                        , &
       &        0.0002374473d0                        , &
       &        0.0000464981d0                          &
       &       ],                                       &
       &       relTol=1.0d-6                            &
       &     )
  call Assert("incomplete gamma function, Γ(x,16.67)"     , &
       &       incompleteGammaFunction2,                    &
       &       [                                            &
       &        5.758521420655205d-8,                       &
       &        1.017530735029775d-6,                       &
       &        9.018676651091340d-6,                       &
       &        5.347837745800675d-5,                       &
       &        2.387641805708268d-4,                       &
       &        8.565070481489690d-4,                       &
       &        2.572802648570240d-3,                       &
       &        6.660038028430611d-3,                       &
       &        1.517681475121466d-2,                       &
       &        3.095177785886023d-2                        &
       &       ],                                           &
       &       relTol=1.0d-6                                &
       &     )
  call Assert("complementary gamma function, 1-Γ(x,16.67)", &
       &       incompleteComplementaryGammaFunction2      , &
       &       [                                            &
       &        9.99999942414786d-1,                        &
       &        9.99998982469265d-1,                        &
       &        9.99990981323349d-1,                        &
       &        9.99946521622542d-1,                        &
       &        9.99761235819429d-1,                        &
       &        9.99143492951851d-1,                        &
       &        9.97427197351430d-1,                        &
       &        9.93339961971569d-1,                        &
       &        9.84823185248785d-1,                        &
       &        9.69048222141140d-1                         &
       &       ],                                           &
       &                                                    &
       &       relTol=1.0d-6                                &
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
  call Assert("hypergeometric (approximate), ₂F₁([3/2,1/2],[3/2],x)", &
       &       hypergeometric2F1approx,                               &
       &       [                                                      &
       &        1.414213562d0,                                        &
       &        1.224744871d0,                                        &
       &        1.154700538d0,                                        &
       &        1.118033989d0,                                        &
       &        1.095445115d0,                                        &
       &        1.080123450d0,                                        &
       &        1.069044968d0,                                        &
       &        1.060660172d0,                                        &
       &        1.054092553d0,                                        &
       &        1.048808848d0                                         &
       &       ],                                                     &
       &       relTol=1.0d-6                                          &
       &     )
  call Assert("hypergeometric, ₁F₂([3/2],[3/2,1/2],x)", &
       &       hypergeometric1F2,                       &
       &       [                                        &
       &        3.762195691d0,                          &
       &        8.488967212d0,                          &
       &        15.98952330d0,                          &
       &        27.30823283d0,                          &
       &        43.77746767d0,                          &
       &        67.08012955d0,                          &
       &        99.32335942d0,                          &
       &        143.1251286d0,                          &
       &        201.7156361d0,                          &
       &        279.0556851d0                           &
       &       ],                                       &
       &       relTol=1.0d-6                            &
       &     )
  call Assert("hypergeometric regularized, ₁F₂([3/2],[3/2,1/2],x)/Γ(3/2)Γ(1/2)", &
       &       hypergeometric1F2Regularized,                                     &
       &       [                                                                 &
       &          2.395088164459957d0,                                           &
       &          5.404244374495762d0,                                           &
       &         10.179246689601360d0,                                           &
       &         17.384960971825730d0,                                           &
       &         27.869601505963740d0,                                           &
       &         42.704536809925250d0,                                           &
       &         63.231214467615680d0,                                           &
       &         91.116286835144320d0,                                           &
       &        128.416162351259700d0,                                           &
       &        177.652366745316700d0                                            &
       &       ],                                                                &
       &       relTol=1.0d-6                                                     &
       &     )

  ! Test hypergeometric 2F1 function transitions for |x|>1.
  hypergeometric2F1(1)=Hypergeometric_2F1([1.0d0,1.0d0],[2.0d0],-2.0d0)
  call Assert("hypergeometric, ₂F₁([1,1],[2],-2)",hypergeometric2F1(1),log(3.0d0)/2.0d0,relTol=1.0d-6 )

  ! Test bugs in the GSL implementation of the hypergeometric function 2F1(a,b,c,x) at x>0.5 when c-a-b is an integer
  ! and at least one of a, b, or c is negative.
  hypergeometric2F1(1)=Hypergeometric_2F1([1.5d0,-0.7d0],[2.8d0],0.7d0)
  call Assert("hypergeometric, ₂F₁([3/2,-7/10],[14/5],7/10)",hypergeometric2F1(1),0.7132641626d0,relTol=1.0d-6)
  hypergeometric2F1approx(1)=Hypergeometric_2F1([1.5d0,-0.7d0],[2.8d0],0.7d0,toleranceRelative=1.0d-6)
  call Assert("hypergeometric (approximate), ₂F₁([3/2,-7/10],[14/5],7/10)",hypergeometric2F1approx(1),0.7132641626d0,relTol=1.0d-6)

  ! Test hypergeometric 2F1 function for x<-1.
  call Assert("hypergeometric, ₂F₁([3/2,7/10],[2],x) (x<-1)", &
       &       hypergeometric2F1NegativeArgument,             &
       &       [                                              &
       &        0.4789885350d0,                               &
       &        0.2913537719d0,                               &
       &        0.1612432319d0,                               &
       &        0.08457261386d0,                              &
       &        0.04313811262d0,                              &
       &        0.02169967246d0,                              &
       &        0.010841876494d0,                             &
       &        0.005399413105d0,                             &
       &        0.002684860664d0,                             &
       &        0.0013340880235d0                             &
       &       ],                                             &
       &       relTol=1.0d-6                                  &
       &     )
  call Assert("hypergeometric (approximate), ₂F₁([3/2,7/10],[2],x) (x<-1)", &
       &       hypergeometric2F1approxNegativeArgument,                     &
       &       [                                                            &
       &        0.4789885350d0,                                             &
       &        0.2913537719d0,                                             &
       &        0.1612432319d0,                                             &
       &        0.08457261386d0,                                            &
       &        0.04313811262d0,                                            &
       &        0.02169967246d0,                                            &
       &        0.010841876494d0,                                           &
       &        0.005399413105d0,                                           &
       &        0.002684860664d0,                                           &
       &        0.0013340880235d0                                           &
       &       ],                                                           &
       &       relTol=1.0d-6                                                &
       &     )

  ! Test hypergeometric 3F2 function for x<-1.
  call Assert("hypergeometric, ₃F₂([3/2,7/10,11/10],[2,1/5],x) (x<-1)", &
       &       hypergeometric3F2NegativeArgument,                       &
       &       [                                                        &
       &        -0.5235443366d0,                                        &
       &        -0.5125266122d0,                                        &
       &        -0.3583570728d0,                                        &
       &        -0.2128977422d0,                                        &
       &        -0.1165025820d0,                                        &
       &        -0.06108436160d0,                                       &
       &        -0.03130060634d0,                                       &
       &        -0.01583617207d0,                                       &
       &        -0.007954088823d0,                                      &
       &        -0.003978068717d0                                       &
       &       ],                                                       &
       &       relTol=1.0d-6                                            &
       &     )

  ! Test hypergeometric 3F2 function with acceleration algorithm.
  call Assert("hypergeometric (accelerated), ₃F₂([3/2,7/10,11/10],[2,1/5],x)", &
       &       hypergeometric3F2accelerated,                                   &
       &       [                                                               &
       &        -0.2171553593d0,                                               &
       &        -0.1108443669d0,                                               &
       &         0.03478627832d0,                                              &
       &         0.2398652189d0,                                               &
       &         0.5395008837d0,                                               &
       &         1.000000000d0,                                                &
       &         1.761012496d0,                                                &
       &         3.167565670d0,                                                &
       &         6.326044401d0,                                                &
       &         17.27894316d0                                                 &
       &       ],                                                              &
       &       relTol=1.0d-6                                                   &
       &     )
  
  ! Test polylogarithm functions.
  call Assert("polylogarithm, Li₂(x)" , &
       &       polylogarithm2         , &
       &       [                        &
       &        -0.82246703342411320d0, &
       &        -0.44841420692364620d0, &
       &        -0.30903312648780850d0, &
       &        -0.23590029768626350d0, &
       &        -0.19080013777753560d0, &
       &        -0.16019301354439550d0, &
       &        -0.13805517651807560d0, &
       &        -0.12129662872272650d0, &
       &        -0.10816821022923270d0, &
       &        -0.09760523522932158d0  &
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )
  call Assert("polylogarithm, Li₃(x)"  , &
       &       polylogarithm3          , &
       &       [                         &
       &        -0.90154267736969570d0 , &
       &        -0.47259784465889690d0 , &
       &        -0.32065094800515400d0 , &
       &        -0.24271200333891630d0 , &
       &        -0.19527359293105430d0 , &
       &        -0.16335479483218130d0 , &
       &        -0.14040803431674340d0 , &
       &        -0.12311562602883610d0 , &
       &        -0.10961645233796000d0 , &
       &        -0.09878555018070007d0   &
       &       ],                        &
       &       relTol=1.0d-6             &
       &     )
    
  ! Test dilogarithm functions.
  call Assert("dilogarithm, Li₂(x)"   , &
       &       dilogarithm_           , &
       &       [                        &
       &        -0.82246703342411320d0, &
       &        -0.44841420692364620d0, &
       &        -0.30903312648780850d0, &
       &        -0.23590029768626350d0, &
       &        -0.19080013777753560d0, &
       &        -0.16019301354439550d0, &
       &        -0.13805517651807560d0, &
       &        -0.12129662872272650d0, &
       &        -0.10816821022923270d0, &
       &        -0.09760523522932158d0  &
       &       ],                       &
       &       relTol=1.0d-6            &
       &     )
  call Assert("dilogarithm, Li₂(2) (complex)" , &
       &       [real(Dilogarithm(dcmplx(2.0d0,0.0d0))),imag(Dilogarithm(dcmplx(2.0d0,0.0d0)))], &
       &       [Pi**2/4.0d0                           ,-Pi*log(2.0d0)]                        , &
       &       relTol=1.0d-6                                                                    &
       &     )
  
  ! Test error function with complex argument.
  errorFunctionComplex=Error_Function(                             &
       &                              [                            &
       &                               dcmplx(10.000d0, 10.000d0), &
       &                               dcmplx(10.000d0,  5.000d0), &
       &                               dcmplx( 5.000d0,  5.000d0), &
       &                               dcmplx( 5.000d0,  1.000d0), &
       &                               dcmplx( 1.000d0,  1.000d0), &
       &                               dcmplx( 1.000d0,  0.500d0), &
       &                               dcmplx( 0.500d0,  0.500d0), &
       &                               dcmplx( 0.500d0,  0.100d0), &
       &                               dcmplx( 0.100d0,  0.100d0), &
       &                               dcmplx( 0.100d0,  0.050d0), &
       &                               dcmplx( 0.050d0,  0.050d0), &
       &                               dcmplx( 0.050d0,  0.010d0), &
       &                               dcmplx( 0.010d0,  0.010d0), &
       &                               dcmplx( 0.010d0,  0.005d0), &
       &                               dcmplx( 0.005d0,  0.005d0), &
       &                               dcmplx( 0.005d0,  0.001d0), &
       &                               dcmplx( 0.001d0,  0.001d0)  &
       &                              ]                            &
       &                             )
  call Assert("error function of complex argument, erf(z)"           , &
       &      errorFunctionComplex                                   , &
       &      [                                                        &
       &       dcmplx(+9.616493742724747d-01, -1.098768460819404d-02), &
       &       dcmplx(+1.000000000000000d+00, -9.495949264558077d-36), &
       &       dcmplx(+9.303796037430947d-01, +3.893619089512146d-02), &
       &       dcmplx(+1.000000000002960d+00, -2.846018382085604d-12), &
       &       dcmplx(+1.316151281697949d+00, +1.904534692378354d-01), &
       &       dcmplx(+9.507097283189570d-01, +1.879734672233839d-01), &
       &       dcmplx(+6.426129148548198d-01, +4.578813944351928d-01), &
       &       dcmplx(+5.249121488205361d-01, +8.802479434588868d-02), &
       &       dcmplx(+1.135856345618654d-01, +1.120811719910652d-01), &
       &       dcmplx(+1.127425509896926d-01, +5.590323090214489d-02), &
       &       dcmplx(+5.651284873688534d-02, +5.632478587819852d-02), &
       &       dcmplx(+5.637760588665819d-02, +1.125599074671486d-02), &
       &       dcmplx(+1.128454387859423d-02, +1.128303937304405d-02), &
       &       dcmplx(+1.128369762595771d-02, +5.641378676150062d-03), &
       &       dcmplx(+5.641989865663000d-03, +5.641801802469853d-03), &
       &       dcmplx(+5.641854461787554d-03, +1.128351334067245d-03), &
       &       dcmplx(+1.128379919345890d-03, +1.128378414842284d-03)  &
       &      ]                                                      , &
       &      relTol=dcmplx(1.0d-6,1.0d-6)&
       &     )

  ! Test binomial coefficients.
  call Assert(                                                                                                                                                                               &
       &             "binomial coefficient, ₂C₁, ₁₀C₃, ₋₂C₀, ₋₂C₁₀, ₋₂C₂₀, ₋₂C₂₀₀"                                                                                                         , &
       &             [Binomial_Coefficient(2,1),Binomial_Coefficient(10,3),Binomial_Coefficient(-2,0),Binomial_Coefficient(-2,10),Binomial_Coefficient(-2,20),Binomial_Coefficient(-2,200)], &
       &             [2.0d0                    ,120.0d0                   ,1.0d0                     ,11.0d0                     ,21.0d0                     ,201.0d0                     ], &
       &      relTol=1.0d-9)

  ! Test beta functions.
  call Assert(                                                                                                &
       &             "beta function, B₀.₅(4000,4000), B₀.₅(4000,3990), B₀.₅(4000,3500), B₀.₅(4000,5000)"    , &
       &             [                                                                                        &
       &              Beta_Function_Incomplete_Normalized(4000.0d0,4.00000000000000d3,0.500000000000000d-00), &
       &              Beta_Function_Incomplete_Normalized(4000.0d0,3.99000000000000d3,0.500000000000000d-00), &
       &              Beta_Function_Incomplete_Normalized(4000.0d0,3.50000000000000d3,0.500000000000000d-00), &
       &              Beta_Function_Incomplete_Normalized(4000.0d0,5.00000000000000d3,0.500000000000000d-00)  &
       &             ]                                                                                      , &
       &             [                                                                                        &
       &              0.5000000000000000d0                                                                  , &
       &              0.4554595983119333d0                                                                  , &
       &              0.0000000000000000d0                                                                  , &
       &              1.0000000000000000d0                                                                    &
       &             ]                                                                                      , &
       &      relTol=1.0d-9,absTol=5.0d-9)
  
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Math_Special_Functions

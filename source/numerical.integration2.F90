!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a variety of numerical integrators.

module Numerical_Integration2
  !% Implements a variety of numerical integrators.
  private
  public :: integrand1D                                , integrandVectorized1D                     , &
       &    integratorCompositeTrapezoidal1D           , integratorVectorizedCompositeTrapezoidal1D, & 
       &    integratorAdaptiveCompositeTrapezoidal1D   , integratorCompositeGaussKronrod1D         , &
       &    integratorVectorizedCompositeGaussKronrod1D, integrator                                , &
       &    integrator1D                               , integratorVectorized1D

  ! Interval type.
  type :: interval
     private
     integer                             :: depth
     double precision                    :: a    , b       , &
          &                                 fa   , fb      , &
          &                                 error, integral
     type            (interval), pointer :: next
  end type interval

  ! Generic integrator.
  type :: integrator
     !% Generic numerical integrator class.
     double precision :: toleranceAbsolute, toleranceRelative
   contains
     !@ <objectMethods>
     !@   <object>integrator</object>
     !@   <objectMethod>
     !@     <method>toleranceSet</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ [toleranceAbsolute]\argin, \doublezero\ [toleranceRelative]\argin</arguments>
     !@     <description>Set tolerances to use in this integrator.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: toleranceSet => toleranceSetGeneric
  end type integrator

  ! Generic one-dimensional integrator.
  type, abstract, extends(integrator) :: integrator1D
     !% Generic one-dimensional numerical integrator class.
     private
     procedure       (integrand1D), pointer, nopass :: integrand
   contains
     !@ <objectMethods>
     !@   <object>integrator1D</object>
     !@   <objectMethod>
     !@     <method>integrandSet</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless procedure(\doublezero\ (\doublezero\ x\argin))\textgreater}\ integrand\argin</arguments>
     !@     <description>Set the integrand function to be integrated.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>evaluate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ a\argin, \doublezero\ b,\argin</arguments>
     !@     <description>Evaluate the integral.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                       :: integrandSet => integrandSet1D
     procedure(evaluate1D), deferred :: evaluate
  end type integrator1D
  abstract interface
     double precision function evaluate1D(self,a,b)
       import integrator1D
       class           (integrator1D), intent(inout) :: self
       double precision              , intent(in   ) :: a,b
     end function evaluate1D
  end interface
  abstract interface
     double precision function integrand1D(x)
       double precision, intent(in   ) :: x
     end function integrand1D
  end interface

  ! Composite trapezoidal 1D integrator.
  type, extends(integrator1D) :: integratorCompositeTrapezoidal1D
     !% One-dimensional numerical integrator class using a composite trapezoidal rule.
     private
     integer :: iterationsMaximum
   contains
     !@ <objectMethods>
     !@   <object>integratorCompositeTrapezoidal1D</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ iterationsMaximum\argin</arguments>
     !@     <description>Set the maximum number of iterations allowed in the integrator.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: initialize => compositeTrapezoidalInitialize1D
     procedure :: evaluate   => compositeTrapezoidalEvaluate1D
  end type integratorCompositeTrapezoidal1D

  ! Adaptive-composite trapezoidal 1D integrator.
  type, extends(integratorCompositeTrapezoidal1D) :: integratorAdaptiveCompositeTrapezoidal1D
     !% One-dimensional numerical integrator class using an adaptive composite trapezoidal rule.
     private
   contains
     procedure :: evaluate => adaptiveCompositeTrapezoidalEvaluate1D
  end type integratorAdaptiveCompositeTrapezoidal1D

  ! Composite Gauss-Kronrod 1D integrator.
  type, extends(integrator1D) :: integratorCompositeGaussKronrod1D
     !% One-dimensional numerical integrator class using a composite Gauss-Kronrod rule.
     private
     integer                                     :: iterationsMaximum
     double precision, allocatable, dimension(:) :: xKronrod         , wGauss, wKronrod
   contains
     !@ <objectMethods>
     !@   <object>integratorCompositeGaussKronrod1D</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ iterationsMaximum\argin, \intzero\ order\argin</arguments>
     !@     <description>Initialize the integrator.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>evaluateInterval</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ a\argin, \doublezero\ b\argin, \doublezero\ integralKronrod\argout, \doublezero\ error\argout</arguments>
     !@     <description>Evaluate the integral over an interval and also return the error on the integral.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: initialize       => compositeGaussKronrod1DInitialize
     procedure :: evaluate         => compositeGaussKronrod1DEvaluate
     procedure :: evaluateInterval => compositeGaussKronrod1DEvaluateInterval
  end type integratorCompositeGaussKronrod1D

  ! Generic one-dimensional vectorized integrator.
  type, abstract, extends(integrator) :: integratorVectorized1D
     !% Generic one-dimensional vectorized numerical integrator class.
     private
     procedure       (integrandVectorized1D), pointer, nopass :: integrand
   contains
     !@ <objectMethods>
     !@   <object>integratorVectorized1D</object>
     !@   <objectMethod>
     !@     <method>integrandSet</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless \doubleone\ function(doubleone\ x\argin)\textgreater}\ integrand\argin</arguments>
     !@     <description>Set the integrand function to be integrated.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure                                 :: integrandSet => integrandVectorizedSet1D
     procedure(evaluateVectorized1D), deferred :: evaluate
  end type integratorVectorized1D
  abstract interface
     double precision function evaluateVectorized1D(self,a,b)
       import integratorVectorized1D
       class           (integratorVectorized1D), intent(inout) :: self
       double precision                        , intent(in   ) :: a,b
     end function evaluateVectorized1D
  end interface
  abstract interface
     function integrandVectorized1D(x)
       double precision, intent(in   ), dimension(     : ) :: x
       double precision               , dimension(size(x)) :: integrandVectorized1D
     end function integrandVectorized1D
  end interface

  ! Vectorized composite trapezoidal 1D integrator.
  type, extends(integratorVectorized1D) :: integratorVectorizedCompositeTrapezoidal1D
     !% One-dimensional numerical integrator class using a vectorized composite trapezoidal rule.
     private
     integer                                     :: iterationsMaximum
     double precision, allocatable, dimension(:) :: d
   contains
     !@ <objectMethods>
     !@   <object>integratorVectorizedCompositeTrapezoidal1D</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ iterationsMaximum\argin</arguments>
     !@     <description>Set the maximum number of iterations allowed in the integrator.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>evaluate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ a\argin, \doublezero\ b,\argin</arguments>
     !@     <description>Evaluate the integral.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: initialize => vectorizedCompositeTrapezoidalInitialize1D
     procedure :: evaluate   => vectorizedCompositeTrapezoidalEvaluate1D
  end type integratorVectorizedCompositeTrapezoidal1D

  ! Vectorized composite Gauss-Kronrod 1D integrator.
  type, extends(integratorVectorized1D) :: integratorVectorizedCompositeGaussKronrod1D
     !% One-dimensional numerical integrator class using a vectorized composite Gauss-Kronrod rule.
     private
     integer                                     :: iterationsMaximum
     double precision, allocatable, dimension(:) :: xKronrod         , wGauss, wKronrod
   contains
     !@ <objectMethods>
     !@   <object>integratorVectorizedCompositeGaussKronrod1D</object>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ iterationsMaximum\argin, \intzero\ order\argin</arguments>
     !@     <description>Set the maximum number of iterations allowed, and the order of the integrator.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>evaluate</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ a\argin, \doublezero\ b,\argin</arguments>
     !@     <description>Evaluate the integral.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>evaluateInterval</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ a\argin, \doublezero\ b\argin, \doublezero\ integralKronrod\argout, \doublezero\ error\argout</arguments>
     !@     <description>Evaluate the integral over an interval and also return the error on the integral.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: initialize       => vectorizedCompositeGaussKronrod1DInitialize
     procedure :: evaluate         => vectorizedCompositeGaussKronrod1DEvaluate
     procedure :: evaluateInterval => vectorizedCompositeGaussKronrod1DEvaluateInterval
  end type integratorVectorizedCompositeGaussKronrod1D

contains

  ! Generic functions.
  subroutine toleranceSetGeneric(self,toleranceAbsolute,toleranceRelative)
    !% Initialize the tolerances for numerical integrators.
    use Galacticus_Error
    implicit none
    class           (integrator), intent(inout)           :: self
    double precision            , intent(in   ), optional :: toleranceAbsolute,toleranceRelative

    if (.not.(present(toleranceAbsolute).or.present(toleranceRelative)))                                &
         &  call Galacticus_Error_Report(                                                               &
         &                               'toleranceSetGeneric'                                        , &
         &                               'at least one of absolute and relative tolerance must be set'  &
         &                              )
    self%toleranceAbsolute=0.0d0
    self%toleranceRelative=0.0d0
    if (present(toleranceAbsolute)) self%toleranceAbsolute=toleranceAbsolute
    if (present(toleranceRelative)) self%toleranceRelative=toleranceRelative
    return
  end subroutine toleranceSetGeneric

  ! Generic one-dimensional integrator.
  subroutine integrandSet1D(self,integrand)
    !% Initialize the integrand for one-dimensional numerical integrators.
    implicit none
    class    (integrator1D), intent(inout) :: self
    procedure(integrand1D )                :: integrand

    self%integrand => integrand
    return
  end subroutine integrandSet1D

  ! Composite trapezoidal 1D integrator.
  subroutine compositeTrapezoidalInitialize1D(self,iterationsMaximum)
    !% Initialize a one-dimensional, composite trapezoidal numerical integrator.
    use Galacticus_Error
    implicit none
    class  (integratorCompositeTrapezoidal1D), intent(inout) :: self
    integer                                  , intent(in   ) :: iterationsMaximum

    if (iterationsMaximum < 3)                                                &
         & call Galacticus_Error_Report(                                      &
         &                              'compositeTrapezoidalInitialize1D'  , &
         &                              'at least 3 iterations are required'  &
         &                             )
    self%iterationsMaximum=iterationsMaximum
    return
  end subroutine compositeTrapezoidalInitialize1D

  double precision function compositeTrapezoidalEvaluate1D(self,a,b)
    !% Evaluate a one-dimension integral using a numerical composite trapezoidal rule.
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    class           (integratorCompositeTrapezoidal1D), intent(inout) :: self
    double precision                                  , intent(in   ) :: a,b
    integer                                                           :: iteration        , j        , &
         &                                                               n
    logical                                                           :: converged
    double precision                                                  :: cumulant         , integrand, &
         &                                                               integrandPrevious, x        , &
         &                                                               step

    cumulant =0.5d0*(self%integrand(b)+self%integrand(a))
    iteration=1
    converged=.false.
    do while (.not.converged)
       iteration=iteration+1
       if (iteration > self%iterationsMaximum)                               &
            & call Galacticus_Error_Report(                                  &
            &                              'compositeTrapezoidalEvaluate1D', &
            &                              'maximum iterations exceeded'     &
            &                             )
       n   =2**iteration
       step=(b-a)/dble(n)
       ! Evaluate only every second point, as the others were computed in previous iterations.
       do j=1,n-1,2
          x=a+dble(j)*step
          cumulant=cumulant+self%integrand(x)
       end do
       integrand=step*cumulant
       if (iteration > 2)                                    &
            & converged=Values_Agree(                        &
            &                        integrand             , &
            &                        integrandPrevious     , &
            &                        self%toleranceAbsolute, &
            &                        self%toleranceRelative  &
            &                       )
      integrandPrevious=integrand
   end do
    compositeTrapezoidalEvaluate1D=integrand
    return
  end function compositeTrapezoidalEvaluate1D

  ! Composite Gauss-Kronrod 1D integrator.

  subroutine compositeGaussKronrod1DInitialize(self,iterationsMaximum,order)
    !% Initialize a one-dimensional, composite Gauss-Kronrod numerical integrator. Evaluation points and weights are taken from
    !% those used in the \gls{gsl}.
    use Galacticus_Error
    use Memory_Management
    implicit none
    class  (integratorCompositeGaussKronrod1D), intent(inout) :: self
    integer                                   , intent(in   ) :: iterationsMaximum, order

    if (iterationsMaximum < 3)                                                 &
         & call Galacticus_Error_Report(                                       &
         &                              'compositeGaussKronrod1DInitialize'  , &
         &                              'at least 3 iterations are required'   &
         &                             )
    self%iterationsMaximum=iterationsMaximum
    ! Choose order.
    select case (order)
    case (15) ! 15-point Kronrod rule.
       call Alloc_Array(self%xKronrod,[8])
       call Alloc_Array(self%wKronrod,[8])
       call Alloc_Array(self%wGauss  ,[4])
       self%xKronrod =[                                       &
            &          0.991455371120812639206854697526329d0, &
            &          0.949107912342758524526189684047851d0, &
            &          0.864864423359769072789712788640926d0, &
            &          0.741531185599394439863864773280788d0, &
            &          0.586087235467691130294144838258730d0, &
            &          0.405845151377397166906606412076961d0, &
            &          0.207784955007898467600689403773245d0, &
            &          0.000000000000000000000000000000000d0  &
            &        ]
       self%wGauss  =[                                        &
            &          0.129484966168869693270611432679082d0, &
            &          0.279705391489276667901467771423780d0, &
            &          0.381830050505118944950369775488975d0, &
            &          0.417959183673469387755102040816327d0  &
            &        ]
       self%wKronrod=[                                        &
            &          0.022935322010529224963732008058970d0, &
            &          0.063092092629978553290700663189204d0, &
            &          0.104790010322250183839876322541518d0, &
            &          0.140653259715525918745189590510238d0, &
            &          0.169004726639267902826583426598550d0, &
            &          0.190350578064785409913256402421014d0, &
            &          0.204432940075298892414161999234649d0, &
            &          0.209482141084727828012999174891714d0  &
            &        ]
    case (61) ! 61-point Kronrod rule.
       call Alloc_Array(self%xKronrod,[31])
       call Alloc_Array(self%wKronrod,[31])
       call Alloc_Array(self%wGauss  ,[15])
       self%xKronrod =[                                       &
            &          0.999484410050490637571325895705811d0, &
            &          0.996893484074649540271630050918695d0, &
            &          0.991630996870404594858628366109486d0, &
            &          0.983668123279747209970032581605663d0, &
            &          0.973116322501126268374693868423707d0, &
            &          0.960021864968307512216871025581798d0, &
            &          0.944374444748559979415831324037439d0, &
            &          0.926200047429274325879324277080474d0, &
            &          0.905573307699907798546522558925958d0, &
            &          0.882560535792052681543116462530226d0, &
            &          0.857205233546061098958658510658944d0, &
            &          0.829565762382768397442898119732502d0, &
            &          0.799727835821839083013668942322683d0, &
            &          0.767777432104826194917977340974503d0, &
            &          0.733790062453226804726171131369528d0, &
            &          0.697850494793315796932292388026640d0, &
            &          0.660061064126626961370053668149271d0, &
            &          0.620526182989242861140477556431189d0, &
            &          0.579345235826361691756024932172540d0, &
            &          0.536624148142019899264169793311073d0, &
            &          0.492480467861778574993693061207709d0, &
            &          0.447033769538089176780609900322854d0, &
            &          0.400401254830394392535476211542661d0, &
            &          0.352704725530878113471037207089374d0, &
            &          0.304073202273625077372677107199257d0, &
            &          0.254636926167889846439805129817805d0, &
            &          0.204525116682309891438957671002025d0, &
            &          0.153869913608583546963794672743256d0, &
            &          0.102806937966737030147096751318001d0, &
            &          0.051471842555317695833025213166723d0, &
            &          0.000000000000000000000000000000000d0  &
            &         ]
       self%wGauss   =[                                       &
            &          0.007968192496166605615465883474674d0, &
            &          0.018466468311090959142302131912047d0, &
            &          0.028784707883323369349719179611292d0, &
            &          0.038799192569627049596801936446348d0, &
            &          0.048402672830594052902938140422808d0, &
            &          0.057493156217619066481721689402056d0, &
            &          0.065974229882180495128128515115962d0, &
            &          0.073755974737705206268243850022191d0, &
            &          0.080755895229420215354694938460530d0, &
            &          0.086899787201082979802387530715126d0, &
            &          0.092122522237786128717632707087619d0, &
            &          0.096368737174644259639468626351810d0, &
            &          0.099593420586795267062780282103569d0, &
            &          0.101762389748405504596428952168554d0, &
            &          0.102852652893558840341285636705415d0  &
            &         ]
       self%wKronrod =[                                       &
            &          0.001389013698677007624551591226760d0, &
            &          0.003890461127099884051267201844516d0, &
            &          0.006630703915931292173319826369750d0, &
            &          0.009273279659517763428441146892024d0, &
            &          0.011823015253496341742232898853251d0, &
            &          0.014369729507045804812451432443580d0, &
            &          0.016920889189053272627572289420322d0, &
            &          0.019414141193942381173408951050128d0, &
            &          0.021828035821609192297167485738339d0, &
            &          0.024191162078080601365686370725232d0, &
            &          0.026509954882333101610601709335075d0, &
            &          0.028754048765041292843978785354334d0, &
            &          0.030907257562387762472884252943092d0, &
            &          0.032981447057483726031814191016854d0, &
            &          0.034979338028060024137499670731468d0, &
            &          0.036882364651821229223911065617136d0, &
            &          0.038678945624727592950348651532281d0, &
            &          0.040374538951535959111995279752468d0, &
            &          0.041969810215164246147147541285970d0, &
            &          0.043452539701356069316831728117073d0, &
            &          0.044814800133162663192355551616723d0, &
            &          0.046059238271006988116271735559374d0, &
            &          0.047185546569299153945261478181099d0, &
            &          0.048185861757087129140779492298305d0, &
            &          0.049055434555029778887528165367238d0, &
            &          0.049795683427074206357811569379942d0, &
            &          0.050405921402782346840893085653585d0, &
            &          0.050881795898749606492297473049805d0, &
            &          0.051221547849258772170656282604944d0, &
            &          0.051426128537459025933862879215781d0, &
            &          0.051494729429451567558340433647099d0  &
            &         ]
    case default
       call Galacticus_Error_Report('compositeGaussKronrod1DInitialize','unknown order')
    end select
    return
  end subroutine compositeGaussKronrod1DInitialize

  double precision function compositeGaussKronrod1DEvaluate(self,a,b)
    !% Evaluate a one-dimension integral using a numerical composite Gauss-Kronrod rule.
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    class           (integratorCompositeGaussKronrod1D), intent(inout) :: self
    double precision                                   , intent(in   ) :: a           , b
    integer                                                            :: iInterval
    logical                                                            :: converged
    double precision                                                   :: error       , midpoint
    type            (interval                         ), pointer       :: head        , newInterval1, &
         &                                                                newInterval2, current     , &
         &                                                                previous    , newInterval

    ! Create our first estimate of the integral, using a single interval.
    allocate(head)
    head%a        =  a
    head%b        =  b
    head%next     => null()
    call self%evaluateInterval(head%a,head%b,head%integral,head%error)
    ! Initialize current integral and error estimates, and flag that we're not converged.
    compositeGaussKronrod1DEvaluate=head%integral
    error                          =head%error
    converged                      =                                                 &
         &  abs(error) < self%toleranceAbsolute                                      &
         & .or.                                                                      &
         &  abs(error) < self%toleranceRelative*abs(compositeGaussKronrod1DEvaluate)
    ! Iterate until convergence is reached.
    do while (.not.converged)
       ! Bisect the head interval. By construction, this will always be the interval with the largest absolute error.
       current  => head      ! Pop the head from the list.
       head     => head%next
       midpoint =  (current%b+current%a)/2.0d0
       ! Create two new subintervals.
       allocate(newInterval1)
       allocate(newInterval2)
       newInterval1%next => null()
       newInterval2%next => null()
       newInterval1%a    =  current%a
       newInterval1%b    =  midpoint
       newInterval2%a    =  midpoint
       newInterval2%b    =  current%b
       ! Compute the integral and error estimate in each new subinterval.
       call self%evaluateInterval(newInterval1%a,newInterval1%b,newInterval1%integral,newInterval1%error)
       call self%evaluateInterval(newInterval2%a,newInterval2%b,newInterval2%integral,newInterval2%error)
       ! Update the total integral and error estimate for the new subintervals.
       compositeGaussKronrod1DEvaluate         &
            & =compositeGaussKronrod1DEvaluate &
            & -current     %integral           &
            & +newInterval1%integral           &
            & +newInterval2%integral
       error                                   &
            & =error                           &
            & -current     %error              &
            & +newInterval1%error              &
            & +newInterval2%error
       ! Test for convergence.
       converged=                                                                       &
            &  abs(error) < self%toleranceAbsolute                                      &
            & .or.                                                                      &
            &  abs(error) < self%toleranceRelative*abs(compositeGaussKronrod1DEvaluate)
       ! Destroy the old interval.
       deallocate(current)
       ! Insert the new intervals into our stack.
       do iInterval=1,2
          if (iInterval == 1) then
             newInterval => newInterval1
          else
             newInterval => newInterval2
          end if
          ! Check if the stack is empty.
          if (associated(head)) then
             ! Stack is not empty. Perform an insertion sort to insert our new subinterval into the sorted stack.
             current  => head
             previous => head 
             do while (associated(current%next).and.abs(current%error) > abs(newInterval%error))
                previous => current
                current  => current%next
             end do
             ! Check if we reached the end of the stack.
             if (.not.associated(current%next)) then
                ! End of stack reached. Check if new interval goes before or after the final stack entry.
                if (abs(current%error) > abs(newInterval%error)) then
                   ! New interval goes at the end of the stack.
                   current%next => newInterval
                else
                   ! New interval goes before the final interval.
                   newInterval%next => previous   %next
                   previous   %next => newInterval
                end if
             else
                ! End of stack not reached. Insert new interval after the previous interval.
                newInterval%next => previous   %next
                previous   %next => newInterval
             end if
          else
             ! Stack is empty - simply point the head to our new subinterval.
             head => newInterval
          end if
       end do
    end do
    ! Destroy the stack.
    current => head
    do while (associated(current))
       previous => current
       current  => current%next
       deallocate(previous)
    end do
    return
  end function compositeGaussKronrod1DEvaluate


! This illustrates how we might vectorize this. On each loop, you want to split the stack in two such that the sum of the errors below the split is less than the required error on the total integral. This means that we know that we need to bisect *every* interval above the split since in the best-case scenario the errors after bisecting these go to zero and we have then just sneaked in under the required error budget. Then, build an array of all midpoints in the intervals to be bisected - do a vector eval of the function at those points - then create the new intervals and populate them with the required function values. Could then look at optimizing the creation of the new intervals:
!==> suppose we're making N new intervals
!==> allocate a pointer array of intervals of size N instead of doing these one at a time
!==> efficient means to sort them? use quickSort for example instead of insertion?


  subroutine compositeGaussKronrod1DEvaluateInterval(self,a,b,integralKronrod,error)
    !% Evaluate the integral over an interval using the Gauss-Kronrod method (also estimates the error). Specific implementation
    !% is based on that in the \gls{gsl}.
    implicit none
    class           (integratorCompositeGaussKronrod1D), intent(inout)                    :: self
    double precision                                   , intent(in   )                    :: a                 , b
    double precision                                   , intent(  out)                    :: integralKronrod   , error
    double precision                                   , dimension(size(self%xKronrod)-1) :: fValue1           , fValue2
    double precision                                                                      :: integralGauss     , fCenter         , &
         &                                                                                   xCenter           , halfLength      , &
         &                                                                                   fSum              , x               , &
         &                                                                                   halfLengthAbsolute, integralAbsolute, &
         &                                                                                   mean              , integralAsc     , &
         &                                                                                   scale
    integer                                                                               :: pointCountKronrod , pointCountGauss , &
         &                                                                                   jtw               , jtwm1           , &
         &                                                                                   j
    
    pointCountKronrod =size(self%wKronrod)
    pointCountGauss   =size(self%wGauss  )
    xCenter           =0.5d0*(b+a)
    halfLength        =0.5d0*(b-a)
    halfLengthAbsolute=abs(halfLength)
    fCenter           =self%integrand(xCenter)
    integralKronrod   =fCenter*self%wKronrod(pointCountKronrod)
    integralGauss     =0.0d0
    integralAbsolute  =abs(integralKronrod)
    if (mod(pointCountKronrod,2) == 0) integralGauss=fCenter*self%wGauss(pointCountKronrod/2)
    do j=0,(pointCountKronrod-1)/2-1
       jtw             =j*2+1
       x               =halfLength*self%xKronrod(jtw+1)
       fValue1(jtw+1)  =self%integrand(xCenter-x)
       fValue2(jtw+1)  =self%integrand(xCenter+x)
       fSum            =fValue1(jtw+1)+fValue2(jtw+1)
       integralGauss   =integralGauss   +self%wGauss  (j  +1)*fSum
       integralKronrod =integralKronrod +self%wKronrod(jtw+1)*fSum
       integralAbsolute=integralAbsolute+self%wKronrod(jtw+1)*(abs(fValue1(jtw+1))+abs(fValue2(jtw+1)))
    end do
    do j=0,pointCountKronrod/2-1
       jtwm1           =j*2
       x               =halfLength*self%xKronrod(jtwm1+1)
       fValue1(jtwm1+1)=self%integrand(xCenter-x)
       fValue2(jtwm1+1)=self%integrand(xCenter+x)
       integralKronrod =integralKronrod +self%wKronrod(jtwm1+1)*(    fValue1(jtwm1+1) +    fValue2(jtwm1+1))
       integralAbsolute=integralAbsolute+self%wKronrod(jtwm1+1)*(abs(fValue1(jtwm1+1))+abs(fValue2(jtwm1+1)))
    end do
    mean        =integralKronrod*0.5d0
    integralAsc=self%wKronrod(pointCountKronrod)*abs(fCenter-mean)
    do j=0,pointCountKronrod-1
       integralAsc=integralAsc+self%wKronrod(j+1)*(abs(fValue1(j+1)-mean)+abs(fValue2(j+1)-mean))
    end do
    error           =abs((integralKronrod-integralGauss)*halfLength)
    integralKronrod =integralKronrod *halfLength
    integralAbsolute=integralAbsolute*halfLengthAbsolute
    integralAsc     =integralAsc     *halfLengthAbsolute
    ! Compute error.
    if (integralAsc /= 0.0d0 .and. error /= 0.0d0) then
       scale=(200.0d0*error/integralAsc)**1.5d0
       if (scale < 1.0d0) then
          error=integralAsc*scale
       else
          error=integralAsc
       end if
    end if
    if (integralAbsolute > tiny(1.0d0)/(50.0d0*epsilon(1.0d0))) error=max(error,50.0d0*epsilon(1.0d0)*integralAbsolute)
    return
  end subroutine compositeGaussKronrod1DEvaluateInterval
  
  ! Vectorized composite Gauss-Kronrod 1D integrator.

  subroutine vectorizedCompositeGaussKronrod1DInitialize(self,iterationsMaximum,order)
    !% Initialize a one-dimensional, vectorized composite Gauss-Kronrod numerical integrator. Evaluation points and weights are
    !% taken from those used in the \gls{gsl}.
    use Galacticus_Error
    use Memory_Management
    implicit none
    class  (integratorVectorizedCompositeGaussKronrod1D), intent(inout) :: self
    integer                                             , intent(in   ) :: iterationsMaximum, order

    if (iterationsMaximum < 3)                                                         &
         & call Galacticus_Error_Report(                                               &
         &                              'vectorizedCompositeGaussKronrod1DInitialize', &
         &                              'at least 3 iterations are required'           &
         &                             )
    self%iterationsMaximum=iterationsMaximum
    ! Choose order.
    select case (order)
    case (15) ! 15-point Kronrod rule.
       call Alloc_Array(self%xKronrod,[8])
       call Alloc_Array(self%wKronrod,[8])
       call Alloc_Array(self%wGauss  ,[4])
       self%xKronrod =[                                       &
            &          0.991455371120812639206854697526329d0, &
            &          0.949107912342758524526189684047851d0, &
            &          0.864864423359769072789712788640926d0, &
            &          0.741531185599394439863864773280788d0, &
            &          0.586087235467691130294144838258730d0, &
            &          0.405845151377397166906606412076961d0, &
            &          0.207784955007898467600689403773245d0, &
            &          0.000000000000000000000000000000000d0  &
            &        ]
       self%wGauss  =[                                        &
            &          0.129484966168869693270611432679082d0, &
            &          0.279705391489276667901467771423780d0, &
            &          0.381830050505118944950369775488975d0, &
            &          0.417959183673469387755102040816327d0  &
            &        ]
       self%wKronrod=[                                        &
            &          0.022935322010529224963732008058970d0, &
            &          0.063092092629978553290700663189204d0, &
            &          0.104790010322250183839876322541518d0, &
            &          0.140653259715525918745189590510238d0, &
            &          0.169004726639267902826583426598550d0, &
            &          0.190350578064785409913256402421014d0, &
            &          0.204432940075298892414161999234649d0, &
            &          0.209482141084727828012999174891714d0  &
            &        ]
    case (61) ! 61-point Kronrod rule.
       call Alloc_Array(self%xKronrod,[31])
       call Alloc_Array(self%wKronrod,[31])
       call Alloc_Array(self%wGauss  ,[15])
       self%xKronrod =[                                       &
            &          0.999484410050490637571325895705811d0, &
            &          0.996893484074649540271630050918695d0, &
            &          0.991630996870404594858628366109486d0, &
            &          0.983668123279747209970032581605663d0, &
            &          0.973116322501126268374693868423707d0, &
            &          0.960021864968307512216871025581798d0, &
            &          0.944374444748559979415831324037439d0, &
            &          0.926200047429274325879324277080474d0, &
            &          0.905573307699907798546522558925958d0, &
            &          0.882560535792052681543116462530226d0, &
            &          0.857205233546061098958658510658944d0, &
            &          0.829565762382768397442898119732502d0, &
            &          0.799727835821839083013668942322683d0, &
            &          0.767777432104826194917977340974503d0, &
            &          0.733790062453226804726171131369528d0, &
            &          0.697850494793315796932292388026640d0, &
            &          0.660061064126626961370053668149271d0, &
            &          0.620526182989242861140477556431189d0, &
            &          0.579345235826361691756024932172540d0, &
            &          0.536624148142019899264169793311073d0, &
            &          0.492480467861778574993693061207709d0, &
            &          0.447033769538089176780609900322854d0, &
            &          0.400401254830394392535476211542661d0, &
            &          0.352704725530878113471037207089374d0, &
            &          0.304073202273625077372677107199257d0, &
            &          0.254636926167889846439805129817805d0, &
            &          0.204525116682309891438957671002025d0, &
            &          0.153869913608583546963794672743256d0, &
            &          0.102806937966737030147096751318001d0, &
            &          0.051471842555317695833025213166723d0, &
            &          0.000000000000000000000000000000000d0  &
            &         ]
       self%wGauss   =[                                       &
            &          0.007968192496166605615465883474674d0, &
            &          0.018466468311090959142302131912047d0, &
            &          0.028784707883323369349719179611292d0, &
            &          0.038799192569627049596801936446348d0, &
            &          0.048402672830594052902938140422808d0, &
            &          0.057493156217619066481721689402056d0, &
            &          0.065974229882180495128128515115962d0, &
            &          0.073755974737705206268243850022191d0, &
            &          0.080755895229420215354694938460530d0, &
            &          0.086899787201082979802387530715126d0, &
            &          0.092122522237786128717632707087619d0, &
            &          0.096368737174644259639468626351810d0, &
            &          0.099593420586795267062780282103569d0, &
            &          0.101762389748405504596428952168554d0, &
            &          0.102852652893558840341285636705415d0  &
            &         ]
       self%wKronrod =[                                       &
            &          0.001389013698677007624551591226760d0, &
            &          0.003890461127099884051267201844516d0, &
            &          0.006630703915931292173319826369750d0, &
            &          0.009273279659517763428441146892024d0, &
            &          0.011823015253496341742232898853251d0, &
            &          0.014369729507045804812451432443580d0, &
            &          0.016920889189053272627572289420322d0, &
            &          0.019414141193942381173408951050128d0, &
            &          0.021828035821609192297167485738339d0, &
            &          0.024191162078080601365686370725232d0, &
            &          0.026509954882333101610601709335075d0, &
            &          0.028754048765041292843978785354334d0, &
            &          0.030907257562387762472884252943092d0, &
            &          0.032981447057483726031814191016854d0, &
            &          0.034979338028060024137499670731468d0, &
            &          0.036882364651821229223911065617136d0, &
            &          0.038678945624727592950348651532281d0, &
            &          0.040374538951535959111995279752468d0, &
            &          0.041969810215164246147147541285970d0, &
            &          0.043452539701356069316831728117073d0, &
            &          0.044814800133162663192355551616723d0, &
            &          0.046059238271006988116271735559374d0, &
            &          0.047185546569299153945261478181099d0, &
            &          0.048185861757087129140779492298305d0, &
            &          0.049055434555029778887528165367238d0, &
            &          0.049795683427074206357811569379942d0, &
            &          0.050405921402782346840893085653585d0, &
            &          0.050881795898749606492297473049805d0, &
            &          0.051221547849258772170656282604944d0, &
            &          0.051426128537459025933862879215781d0, &
            &          0.051494729429451567558340433647099d0  &
            &         ]
    case default
       call Galacticus_Error_Report('vectorizedCompositeGaussKronrod1DInitialize','unknown order')
    end select
    return
  end subroutine vectorizedCompositeGaussKronrod1DInitialize

  double precision function vectorizedCompositeGaussKronrod1DEvaluate(self,a,b)
    !% Evaluate a one-dimension integral using a numerical composite Gauss-Kronrod rule.
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    class           (integratorVectorizedCompositeGaussKronrod1D), intent(inout) :: self
    double precision                                   , intent(in   ) :: a           , b
    integer                                                            :: iInterval
    logical                                                            :: converged
    double precision                                                   :: error       , midpoint
    type            (interval                         ), pointer       :: head        , newInterval1, &
         &                                                                newInterval2, current     , &
         &                                                                previous    , newInterval

    ! Create our first estimate of the integral, using a single interval.
    allocate(head)
    head%a        =  a
    head%b        =  b
    head%next     => null()
    call self%evaluateInterval(head%a,head%b,head%integral,head%error)
    ! Initialize current integral and error estimates, and flag that we're not converged.
    vectorizedCompositeGaussKronrod1DEvaluate=head%integral
    error                          =head%error
    converged                      =                                                 &
         &  abs(error) < self%toleranceAbsolute                                      &
         & .or.                                                                      &
         &  abs(error) < self%toleranceRelative*abs(vectorizedCompositeGaussKronrod1DEvaluate)
    ! Iterate until convergence is reached.
    do while (.not.converged)
       ! Bisect the head interval. By construction, this will always be the interval with the largest absolute error.
       current  => head      ! Pop the head from the list.
       head     => head%next
       midpoint =  (current%b+current%a)/2.0d0
       ! Create two new subintervals.
       allocate(newInterval1)
       allocate(newInterval2)
       newInterval1%next => null()
       newInterval2%next => null()
       newInterval1%a    =  current%a
       newInterval1%b    =  midpoint
       newInterval2%a    =  midpoint
       newInterval2%b    =  current%b
       ! Compute the integral and error estimate in each new subinterval.
       call self%evaluateInterval(newInterval1%a,newInterval1%b,newInterval1%integral,newInterval1%error)
       call self%evaluateInterval(newInterval2%a,newInterval2%b,newInterval2%integral,newInterval2%error)
       ! Update the total integral and error estimate for the new subintervals.
       vectorizedCompositeGaussKronrod1DEvaluate         &
            & =vectorizedCompositeGaussKronrod1DEvaluate &
            & -current     %integral                     &
            & +newInterval1%integral                     &
            & +newInterval2%integral
       error                                             &
            & =error                                     &
            & -current     %error                        &
            & +newInterval1%error                        &
            & +newInterval2%error
       ! Test for convergence.
       converged=                                                                       &
            &  abs(error) < self%toleranceAbsolute                                      &
            & .or.                                                                      &
            &  abs(error) < self%toleranceRelative*abs(vectorizedCompositeGaussKronrod1DEvaluate)
       ! Destroy the old interval.
       deallocate(current)
       ! Insert the new intervals into our stack.
       do iInterval=1,2
          if (iInterval == 1) then
             newInterval => newInterval1
          else
             newInterval => newInterval2
          end if
          ! Check if the stack is empty.
          if (associated(head)) then
             ! Stack is not empty. Perform an insertion sort to insert our new subinterval into the sorted stack.
             current  => head
             previous => head 
             do while (associated(current%next).and.abs(current%error) > abs(newInterval%error))
                previous => current
                current  => current%next
             end do
             ! Check if we reached the end of the stack.
             if (.not.associated(current%next)) then
                ! End of stack reached. Check if new interval goes before or after the final stack entry.
                if (abs(current%error) > abs(newInterval%error)) then
                   ! New interval goes at the end of the stack.
                   current%next => newInterval
                else
                   ! New interval goes before the final interval.
                   newInterval%next => previous   %next
                   previous   %next => newInterval
                end if
             else
                ! End of stack not reached. Insert new interval after the previous interval.
                newInterval%next => previous   %next
                previous   %next => newInterval
             end if
          else
             ! Stack is empty - simply point the head to our new subinterval.
             head => newInterval
          end if
       end do
    end do
    ! Destroy the stack.
    current => head
    do while (associated(current))
       previous => current
       current  => current%next
       deallocate(previous)
    end do
    return
  end function vectorizedCompositeGaussKronrod1DEvaluate

  subroutine vectorizedCompositeGaussKronrod1DEvaluateInterval(self,a,b,integralKronrod,error)
    !% Evaluate the integral over an interval using the Gauss-Kronrod method (also estimates the error). Specific implementation
    !% is based on that in the \gls{gsl}.
    implicit none
    class           (integratorVectorizedCompositeGaussKronrod1D), intent(inout)                                     :: self
    double precision                                             , intent(in   )                                     :: a                 , b
    double precision                                             , intent(  out)                                     :: integralKronrod   , error
    double precision                                             , pointer      , dimension(:                      ) :: fValue1           , fValue2
    double precision                                             , target       , dimension(2*size(self%xKronrod)-1) :: fUnion
    double precision                                                            , dimension(2*size(self%xKronrod)-1) :: xUnion
    double precision                                                                                                 :: integralGauss     , scale           , &
         &                                                                                                              xCenter           , halfLength      , &
         &                                                                                                              halfLengthAbsolute, integralAbsolute, &
         &                                                                                                              mean              , integralAsc
    integer                                                                                                          :: pointCountKronrod , pointCountGauss

    ! Establish point counts and interval.
    pointCountKronrod =size(self%wKronrod)
    pointCountGauss   =size(self%wGauss  )
    if (mod(pointCountKronrod,2) == 0) pointCountGauss=pointCountGauss-1
    xCenter           =0.5d0*(b+a)
    halfLength        =0.5d0*(b-a)
    halfLengthAbsolute=abs(halfLength)
    ! Evaluate function at all required points.
    xUnion(1                                            )=  xCenter
    xUnion(2                    :2+  pointCountKronrod-2)=  xCenter-self%xKronrod(1:pointCountKronrod-1)*halfLength
    xUnion(2+pointCountKronrod-1:2+2*pointCountKronrod-3)=  xCenter+self%xKronrod(1:pointCountKronrod-1)*halfLength
    fUnion                                               =  self%integrand(xUnion)
    fValue1                                              => fUnion(2                    :2+  pointCountKronrod-2)
    fValue2                                              => fUnion(2+pointCountKronrod-1:2+2*pointCountKronrod-3)
    ! Evaluate integrals.
    integralGauss=sum(                                    &
         &            +self%wGauss(1:  pointCountGauss  ) &
         &            *(                                  &
         &              +fValue1  (2:2*pointCountGauss:2) &
         &              +fValue2  (2:2*pointCountGauss:2) &
         &             )                                  &
         &           )
    if (mod(pointCountKronrod,2) == 0) integralGauss=integralGauss+fUnion(1)*self%wGauss(pointCountKronrod/2)
    integralKronrod=+    self%wKronrod  (  pointCountKronrod  )  &
         &               *       fUnion (1                    )  &
         &          +sum(                                        &
         &               +self%wKronrod (1:pointCountKronrod-1)  &
         &               *(                                      &
         &                 +     fvalue1(1:pointCountKronrod-1)  &
         &                 +     fvalue2(1:pointCountKronrod-1)  &
         &                )                                      &
         &              )
    integralAbsolute=+self%wKronrod     (  pointCountKronrod  )  &
         &                *  abs(fUnion (1                    )) &
         &           +sum(                                       &
         &                +self%wKronrod(1:pointCountKronrod-1)  &
         &                *(                                     &
         &                  +abs(fvalue1(1:pointCountKronrod-1)) &
         &                  +abs(fvalue2(1:pointCountKronrod-1)) &
         &                 )                                     &
         &               ) 
    mean       =integralKronrod*0.5d0
    integralAsc=+sum(self%wKronrod                   *(abs(fValue1  -mean)+abs(fValue2-mean))) &
         &      +    self%wKronrod(pointCountKronrod)* abs(fUnion(1)-mean)
    ! Evaluate error.
    error           =abs((integralKronrod-integralGauss)*halfLength)
    integralKronrod =integralKronrod *halfLength
    integralAbsolute=integralAbsolute*halfLengthAbsolute
    integralAsc     =integralAsc     *halfLengthAbsolute
    ! Compute error.
    if (integralAsc /= 0.0d0 .and. error /= 0.0d0) then
       scale=(200.0d0*error/integralAsc)**1.5d0
       if (scale < 1.0d0) then
          error=integralAsc*scale
       else
          error=integralAsc
       end if
    end if
    if (integralAbsolute > tiny(1.0d0)/(50.0d0*epsilon(1.0d0))) error=max(error,50.0d0*epsilon(1.0d0)*integralAbsolute)
    return
  end subroutine vectorizedCompositeGaussKronrod1DEvaluateInterval
  
  ! Adaptive composite trapezoidal integrator.

  double precision function adaptiveCompositeTrapezoidalEvaluate1D(self,a,b)
    !% Evaluate a one-dimension integral using a numerical composite trapezoidal rule.
    use Numerical_Comparison
    use Galacticus_Error
    implicit none
    class           (integratorAdaptiveCompositeTrapezoidal1D), intent(inout) :: self
    double precision                                          , intent(in   ) :: a,b
    logical                                                                   :: converged
    double precision                                                          :: error,width,midpoint,fa,fb,fm
    type            (interval                                ), pointer       :: head,newInterval1,newInterval2,current,previous

    ! Create our first estimate of the integral, using a single trapezoid.
    allocate(head)
    head%a        =                  a
    head%b        =                  b
    head%fa       =   self%integrand(a)
    head%fb       =   self%integrand(b)
    head%depth    =   1
    head%next     =>  null()
    head%integral =   0.5d0             &
         &           *(head%fb+head%fa) &
         &           *(head% b-head% a)
    head%error    =   head%integral
    ! Initialize current integral and error estimates, and flag that we're not converged.
    adaptiveCompositeTrapezoidalEvaluate1D=head%integral
    error                                 =head%integral
    converged                             =.false.
    ! Iterate until convergence is reached.
    do while (.not.converged)
       ! Bisect the head interval. By construction, this will always be the interval with the largest absolute error.
       current  => head      ! Pop the head from the list.
       head     => head%next
       width    =  (current%b-current%a)/2.0d0
       midpoint =  (current%b+current%a)/2.0d0
       ! Evaluate the function at the endpoints and midpoints of the current interval. For endpoints, reuse stored values.
       fa       =  current%fa
       fb       =  current%fb
       fm       =  self%integrand(midpoint)
       ! Create two new subintervals.
       allocate(newInterval1)
       allocate(newInterval2)
       newInterval1%next     => null()
       newInterval2%next     => null()
       newInterval1%a        =  current%a
       newInterval1%b        =  midpoint
       newInterval2%a        =  midpoint
       newInterval2%b        =  current%b
       newInterval1%fa       =  fa
       newInterval1%fb       =  fm
       newInterval2%fa       =  fm
       newInterval2%fb       =  fb
       newInterval1%depth    =  current     %depth+1
       newInterval2%depth    =  newInterval1%depth
       ! Compute the integral and error estimate in each new subinterval.
       newInterval1%integral =  (fm+fa)*0.5d0*width
       newInterval2%integral =  (fm+fb)*0.5d0*width
       newInterval1%error    =  abs(0.5d0*width*fm-0.25d0*current%integral)
       newInterval2%error    =  newInterval1%error
       ! Update the total integral and error estimate for the new subintervals.
       adaptiveCompositeTrapezoidalEvaluate1D         &
            & =adaptiveCompositeTrapezoidalEvaluate1D &
            & -current     %integral                  &
            & +newInterval1%integral                  &
            & +newInterval2%integral
       error                                          &
            & =error                                  &
            & -current     %error                     &
            & +newInterval1%error                     &
            & +newInterval2%error
       ! Test for convergence.
       converged=                                                                              &
            &  abs(error) < self%toleranceAbsolute                                             &
            & .or.                                                                             &
            &  abs(error) < self%toleranceRelative*abs(adaptiveCompositeTrapezoidalEvaluate1D)
       ! Destroy the old interval.
       deallocate(current)
       ! Insert the new intervals into our stack.
       ! First, link the two subintervals - they will always be consecutive in the list as they have the same error estimate.
       newInterval1%next => newInterval2
       ! Check if the stack is empty.
       if (associated(head)) then
          ! Stack is not empty. Perform an insertion sort to insert our new subintervals into the sorted stack.
          current  => head
          previous => head 
          do while (associated(current%next).and.abs(current%error) > abs(newInterval1%error))
             previous => current
             current  => current%next
          end do
          ! Check if we reached the end of the stack.
          if (.not.associated(current%next)) then
             ! End of stack reached. Check if new intervals go before or after the final stack entry.
             if (abs(current%error) > abs(newInterval1%error)) then
                ! New intervals go at the end of the stack.
                current%next => newInterval1
             else
                ! New intervals go before the final interval.
                newInterval2%next => previous    %next
                previous    %next => newInterval1
             end if
          else
             ! End of stack not reached. Insert new intervals after the previous interval.
             newInterval2%next => previous    %next
             previous    %next => newInterval1
          end if
       else
          ! Stack is empty - simply point the head to the first of our new subintervals.
          head => newInterval1
       end if
    end do
    ! Destroy the stack.
    current => head
    do while (associated(current))
       previous => current
       current  => current%next
       deallocate(previous)
    end do
    return
  end function adaptiveCompositeTrapezoidalEvaluate1D

  ! Generic one-dimensional vectorized integrator.
  subroutine integrandVectorizedSet1D(self,integrand)
    !% Initialize the integrand for one-dimensional numerical integrators.
    implicit none
    class    (integratorVectorized1D), intent(inout) :: self
    procedure(integrandVectorized1D )                :: integrand

    self%integrand => integrand
    return
  end subroutine integrandVectorizedSet1D

  ! Vectorized composite trapezoidal 1D integrator.
  subroutine vectorizedCompositeTrapezoidalInitialize1D(self,iterationsMaximum)
    !% Initialize a one-dimensional, vectorized composite trapezoidal numerical integrator.
    use Galacticus_Error
    use Memory_Management
    implicit none
    class           (integratorVectorizedCompositeTrapezoidal1D), intent(inout) :: self
    integer                                                     , intent(in   ) :: iterationsMaximum
    integer                                                                     :: workspaceSize    , i, &
         &                                                                         j                , n
    double precision                                                            :: step

    if (iterationsMaximum < 3)                                                        &
         & call Galacticus_Error_Report(                                              &
         &                              'vectorizedCompositeTrapezoidalInitialize1D', &
         &                              'at least 3 iterations are required'          &
         &                             )
    self%iterationsMaximum=iterationsMaximum
    workspaceSize=2**iterationsMaximum
    call Alloc_Array(self%d,[workspaceSize])
    self%d(1:2)=[0.0d0,1.0d0]
    do i=2,iterationsMaximum
       n=2**i
       step=1.0d0/dble(n)
       forall(j=1:n-1:2)
          self%d(n/2+1+j/2)=dble(j)*step
       end forall
    end do
    return
  end subroutine vectorizedCompositeTrapezoidalInitialize1D

  double precision function vectorizedCompositeTrapezoidalEvaluate1D(self,a,b)
    !% Evaluate a one-dimension integral using a numerical vectorized composite trapezoidal rule.
    use Numerical_Comparison
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
#ifdef YEPPP
    use yepCore
    use yepMath
#endif
    implicit none
    class           (integratorVectorizedCompositeTrapezoidal1D), intent(inout)                        :: self
    double precision                                            , intent(in   )                        :: a,b
    integer                                                                                            :: iteration        , n
    logical                                                                                            :: converged
    double precision                                                                                   :: cumulant         , integrand, &
         &                                                                                                integrandPrevious, step
#ifdef YEPPP
    integer         (c_int                                     )                                       :: s
    double precision                                                                                   :: summand
    double precision                                            , dimension(2**self%iterationsMaximum) :: functionValues
#endif

    cumulant =0.5d0*sum(self%integrand([b,a]))
    iteration=1
    converged=.false.
    do while (.not.converged)
       iteration=iteration+1
       if (iteration > self%iterationsMaximum)                                         &
            & call Galacticus_Error_Report(                                            &
            &                              'vectorizedCompositeTrapezoidalEvaluate1D', &
            &                              'maximum iterations exceeded'               &
            &                             )
        n   =2**iteration
        step=(b-a)/dble(n)
        ! Evaluate only every second point, as the others were computed in previous iterations.
#ifdef YEPPP
        functionValues(1:n/2)=self%integrand(self%d(n/2+1:n)*(b-a)+a)
        s                    =yepCore_Sum_V64f_S64f(functionValues,summand,int(n/2,kind=c_size_t))
        cumulant             =cumulant+summand
#else
        cumulant             =cumulant+sum(self%integrand(self%d(n/2+1:n)*(b-a)+a))
#endif
        integrand=step*cumulant
        if (iteration > 2)                                    &
             & converged=Values_Agree(                        &
             &                        integrand             , &
             &                        integrandPrevious     , &
             &                        self%toleranceAbsolute, &
             &                        self%toleranceRelative  &
             &                       )
       integrandPrevious=integrand
    end do
    vectorizedCompositeTrapezoidalEvaluate1D=integrand
    return
  end function vectorizedCompositeTrapezoidalEvaluate1D

end module Numerical_Integration2

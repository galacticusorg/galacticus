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

!!{
Contains a module which implements quasi-random sequences.
!!}

! Specify an explicit dependence on the interface.GSL.C.odeiv2.o object file.
!: $(BUILDPATH)/interface.GSL.C.odeiv2.o
! Specify an explicit dependence on  the Multi-Step Backward Differentiation Formula with Active/Inactive
! variables, and LU decomposition with Active/Inactive variables object files.
!: $(BUILDPATH)/gslODEInitVal2/driver2.o
!: $(BUILDPATH)/gslODEInitVal2/cscal2.o
!: $(BUILDPATH)/gslODEInitVal2/msbdfactive.o
!: $(BUILDPATH)/gslODEInitVal2/lu.o

! Add dependency on GSL library.
!; gsl

module Numerical_ODE_Solvers
  !!{
  Implements an ODE solver class.
  !!}
  use, intrinsic :: ISO_C_Binding         , only : c_double                   , c_funptr  , c_int, c_ptr, &
          &                                        c_size_t                   , c_null_ptr
  use            :: Numerical_Integration2, only : integratorMultiVectorized1D
  use            :: Resource_Manager      , only : resourceManager
  implicit none
  private
  public :: odeSolver

  ! Stepper types.
  integer, public, parameter :: gsl_odeiv2_step_rk2        = 1
  integer, public, parameter :: gsl_odeiv2_step_rk4        = 2
  integer, public, parameter :: gsl_odeiv2_step_rkf45      = 3
  integer, public, parameter :: gsl_odeiv2_step_rkck       = 4
  integer, public, parameter :: gsl_odeiv2_step_rk8pd      = 5
  integer, public, parameter :: gsl_odeiv2_step_rk1imp     = 6
  integer, public, parameter :: gsl_odeiv2_step_rk2imp     = 7
  integer, public, parameter :: gsl_odeiv2_step_rk4imp     = 8
  integer, public, parameter :: gsl_odeiv2_step_bsimp      = 9
  integer, public, parameter :: gsl_odeiv2_step_msadams    =10
  integer, public, parameter :: gsl_odeiv2_step_msbdf      =11
  integer, public, parameter :: gsl_odeiv2_step_msbdfactive=12

  interface
     function gsl_odeiv2_step_type_get(i) bind(c,name='gsl_odeiv2_step_type_get')
       !!{
       Template for GSL interface ODEIV2 stepper type function.
       !!}
       import c_ptr, c_int
       type   (c_ptr)                       :: gsl_odeiv2_step_type_get
       integer(c_int), intent(in   ), value :: i
     end function gsl_odeiv2_step_type_get
     function gsl_odeiv2_system_init(dim,func,jacobian) bind(c,name='gsl_odeiv2_system_init')
       !!{
       Template for GSL interface ODEIV2 system initialization function.
       !!}
       import c_ptr, c_funptr, c_size_t
       type   (c_ptr   )                       :: gsl_odeiv2_system_init
       integer(c_size_t), intent(in   ), value :: dim
       type   (c_funptr), intent(in   ), value :: func                  , jacobian
     end function gsl_odeiv2_system_init
     subroutine gsl_odeiv2_system_free(system) bind(c,name='gsl_odeiv2_system_free')
       !!{
       Template for GSL interface ODEIV2 system free function.
       !!}
       import c_ptr
       type(c_ptr), intent(in   ), value :: system
     end subroutine gsl_odeiv2_system_free
     function gsl_odeiv2_driver_alloc_y_new(sys,T,hstart,epsabs,epsrel) bind(c,name='gsl_odeiv2_driver_alloc_y_new')
       !!{
       Template for GSL interface ODEIV2 step allocation function.
       !!}
       import c_ptr, c_double
       type(c_ptr   )                       :: gsl_odeiv2_driver_alloc_y_new
       type(c_ptr   ), intent(in   ), value :: sys                          , T
       real(c_double), intent(in   ), value :: epsabs                       , epsrel, &
            &                                  hstart
     end function gsl_odeiv2_driver_alloc_y_new
     function gsl_odeiv2_driver_alloc_scaled2_new(sys,T,hstart,epsabs,epsrel,a_y,a_dydt,scale_abs,is_non_negative) bind(c,name='gsl_odeiv2_driver_alloc_scaled2_new')
       !!{
       Template for GSL interface ODEIV2 step allocation function.
       !!}
       import c_ptr, c_double, c_int
       type   (c_ptr   )                              :: gsl_odeiv2_driver_alloc_scaled2_new
       type   (c_ptr   ), intent(in   ), value        :: sys                                , T
       real   (c_double), intent(in   ), value        :: epsabs                             , epsrel, &
            &                                            a_y                                , a_dydt, &
            &                                            hstart
       real   (c_double), intent(in   ), dimension(*) :: scale_abs
       integer(c_int   ), intent(in   ), dimension(*) :: is_non_negative
     end function gsl_odeiv2_driver_alloc_scaled2_new
     subroutine gsl_odeiv2_driver_free(d) bind(c,name='gsl_odeiv2_driver_free')
       !!{
       Template for GSL interface ODEIV2 driver free function.
       !!}
       import c_ptr
       type(c_ptr), intent(in   ), value :: d
     end subroutine gsl_odeiv2_driver_free     
     function gsl_odeiv2_driver_reset(d) bind(c,name='gsl_odeiv2_driver_reset')
       !!{
       Template for GSL interface ODEIV2 driver reset function.
       !!}
       import c_ptr, c_int
       integer(c_int)                       :: gsl_odeiv2_driver_reset
       type   (c_ptr), intent(in   ), value :: d
     end function gsl_odeiv2_driver_reset
     function gsl_odeiv2_driver_reset_hstart(d,hstart) bind(c,name='gsl_odeiv2_driver_reset_hstart')
       !!{
       Template for GSL interface ODEIV2 driver timestep reset function.
       !!}
       import c_ptr, c_int, c_double
       integer(c_int   )                       :: gsl_odeiv2_driver_reset_hstart
       type   (c_ptr   ), intent(in   ), value :: d
       real   (c_double), intent(in   ), value :: hstart
     end function gsl_odeiv2_driver_reset_hstart
     function gsl_odeiv2_driver2_apply(d,t,t1,y,ps,li,sa) bind(c,name='gsl_odeiv2_driver2_apply')
       !!{
       Interface for GSL interface ODEIV2 driver apply function.
       !!}
       import c_ptr, c_funptr, c_double, c_int
       integer(c_int   )                              :: gsl_odeiv2_driver2_apply
       type   (c_ptr   ), intent(in   ), value        :: d
       real   (c_double), intent(inout)               :: t
       real   (c_double), intent(in   ), value        :: t1
       real   (c_double), intent(inout), dimension(*) :: y
       type   (c_funptr), intent(in   ), value        :: ps
       type   (c_funptr), intent(in   ), value        :: li
       type   (c_funptr), intent(in   ), value        :: sa
     end function gsl_odeiv2_driver2_apply
     function gsl_odeiv2_driver_h(d) bind(c,name='gsl_odeiv2_driver_h')
       !!{
       Template for GSL interface to ODEIV2 driver timestep get function.
       !!}
       import c_ptr, c_double
       real(c_double)                       :: gsl_odeiv2_driver_h
       type(c_ptr   ), intent(in   ), value :: d
     end function gsl_odeiv2_driver_h     
     subroutine msbdfactive_context(d,dim,t,y,dydt) bind(c)
       !!{
       Template for MSBDF active stepper context.
       !!}
       import c_ptr, c_size_t, c_double
       type   (c_ptr   ), intent(in   ), value        :: d
       integer(c_size_t), intent(in   ), value        :: dim
       real   (c_double), intent(in   ), value        :: t
       real   (c_double), intent(inout), dimension(*) :: y  , dydt
     end subroutine msbdfactive_context
     subroutine msbdfactive_state(d,dim,y) bind(c,name='msbdfactive_state')
       !!{
       Template for MSBDF active ODE stepper state function.
       !!}
       import c_ptr, c_size_t, c_double
       type   (c_ptr   ), intent(in   ), value        :: d
       integer(c_size_t), intent(in   ), value        :: dim
       real   (c_double), intent(inout), dimension(*) :: y
     end subroutine msbdfactive_state
     subroutine gsl_odeiv2_driver_errors(d,yerr) bind(c,name='gsl_odeiv2_driver_errors')
       !!{
       Template for GSL ODE driver errors.
       !!}
       import c_ptr, c_double
       type   (c_ptr   ), intent(in   ), value        :: d
       real   (c_double), intent(  out), dimension(*) :: yerr
     end subroutine gsl_odeiv2_driver_errors
     subroutine gsl_odeiv2_driver_init_errors(d) bind(c,name='gsl_odeiv2_driver_init_errors')
       !!{
       Template for GSL ODE driver error initialization.
       !!}
       import c_ptr
       type(c_ptr), intent(in   ), value :: d
     end subroutine gsl_odeiv2_driver_init_errors
  end interface
  
  type :: gslODEDriverWrapper
     !!{
     Wrapper class for managing GSL ODE drivers.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: gslODEDriverWrapperDestructor
  end type gslODEDriverWrapper
  
  type :: gslODESystemWrapper
     !!{
     Wrapper class for managing GSL ODE drivers.
     !!}
     type(c_ptr) :: gsl=c_null_ptr
   contains
     final :: gslODESystemWrapperDestructor
  end type gslODESystemWrapper
  
  type :: odeSolver
     !!{
     Type providing ODE solving.
     !!}
     private
     procedure(derivativesTemplate        ), pointer    , nopass :: derivatives
     procedure(jacobianTemplate           ), pointer    , nopass :: jacobian
     procedure(integrandTemplate          ), pointer    , nopass :: integrands              => null()
     procedure(finalStateTemplate         ), pointer    , nopass :: finalState              => null()
     procedure(postStepTemplate           ), pointer    , nopass :: postStep                => null()
     procedure(errorAnalyzerTemplate      ), pointer    , nopass :: errorAnalyzer           => null()
     procedure(errorHandlerTemplate       ), pointer    , nopass :: errorHandler            => null()
     class    (integratorMultiVectorized1D), pointer             :: integrator              => null()
     type     (gslODEDriverWrapper        ), pointer             :: driver                  => null()
     type     (gslODESystemWrapper        ), pointer             :: system                  => null()
     type     (c_ptr                      ), allocatable         :: gsl_odeiv2_step_type
     type     (resourceManager            )                      :: driverManager                    , systemManager
     integer                                                     :: stepperType
     integer  (c_size_t                   )                      :: dim
     logical                                                     :: integratorErrorTolerant
   contains
     !![
     <methods>
       <method description="Solve the ODE system."                            method="solve" />
       <method description="Return estimates of the errors in ODE variables." method="errors"/>
     </methods>
     !!]
     final     ::           odeSolverDestructor
     procedure :: solve  => odeSolverSolve
     procedure :: errors => odeSolverErrors
  end type odeSolver
  
  interface odeSolver
     !!{
     Constructor for the {\normalfont \ttfamily odeSolver} class.
     !!}
     module procedure odeSolverConstructor
  end interface odeSolver
  
  ! Interface for the derivatives function.
  abstract interface
     integer function derivativesTemplate(x,y,dydx)
       double precision, intent(in   )               :: x
       double precision, intent(in   ), dimension(:) :: y
       double precision, intent(  out), dimension(:) :: dydx
     end function derivativesTemplate
  end interface

  ! Interface for the Jacobian function
  abstract interface
     integer function jacobianTemplate(x,y,dfdy,dfdx)
       double precision, intent(in   )               :: x
       double precision, intent(in   ), dimension(:) :: y
       double precision, intent(  out), dimension(:) :: dfdy
       double precision, intent(  out), dimension(:) :: dfdx
     end function jacobianTemplate
  end interface

  ! Integrand interface.
  abstract interface
     subroutine integrandTemplate(x,y,dydx,z0,e,dzdx)
       double precision, intent(in   ), dimension(  :) :: x
       double precision, intent(in   ), dimension(:,:) :: y   , dydx
       double precision, intent(in   ), dimension(:  ) :: z0
       logical         , intent(inout), dimension(:  ) :: e
       double precision, intent(  out), dimension(:,:) :: dzdx
     end subroutine integrandTemplate
  end interface

  ! Final state interface.
  abstract interface
     subroutine finalStateTemplate(x,y)
       double precision, intent(in   )               :: x
       double precision, intent(in   ), dimension(:) :: y
     end subroutine finalStateTemplate
  end interface
  
  ! Post-step interface.
  abstract interface
     subroutine postStepTemplate(x,y,status)
       import c_int, c_double
       real   (c_double), intent(in   ), value        :: x
       real   (c_double), intent(inout), dimension(*) :: y
       integer(c_int   ), intent(inout)               :: status
     end subroutine postStepTemplate
  end interface

  ! Error analyzer interface.
  abstract interface
     subroutine errorAnalyzerTemplate(x,x1,y,yerr,xStep,status)
       import c_int, c_double
       real   (c_double), intent(in   ), dimension(*) :: y     , yerr
       real   (c_double), intent(in   ), value        :: x     , x1  , &
            &                                            xStep
       integer(c_int   ), intent(in   ), value        :: status
     end subroutine errorAnalyzerTemplate
  end interface

  ! Error handler interface.
  abstract interface
     subroutine errorHandlerTemplate(status,x,xStep,y)
       import c_int, c_double
       integer(c_int   ), intent(in   )               :: status
       real   (c_double), intent(in   )               :: x     , xStep
       real   (c_double), intent(in   ), dimension(:) :: y
     end subroutine errorHandlerTemplate
  end interface

  type :: solverList
     !!{
     Type used to maintain a list of ODE solvers when ODE solving is performed recursively.
     !!}
     type(odeSolver), pointer :: solver => null()
  end type solverList

  ! List of currently active ODE solvers.
  integer                                        :: active =0
  type   (solverList), allocatable, dimension(:) :: solvers
  !$omp threadprivate(solvers,active)

contains

  function odeSolverConstructor(dim,derivatives,jacobian,integrator,integrands,integratorErrorTolerant,stepperType,toleranceAbsolute,toleranceRelative,hStart,dydtScale,yScale,scale,finalState,postStep,errorAnalyzer,errorHandler,isNonNegative) result(self)
    !!{
    Constructor for {\normalfont \ttfamily odeSolver} objects.
    !!}
    use            :: Error        , only : Error_Report
    use, intrinsic :: ISO_C_Binding, only : c_funloc    , c_null_funptr
    implicit none
    type            (odeSolver                  )                                        :: self
    integer         (c_size_t                   ), intent(in   )                         :: dim
    procedure       (derivativesTemplate        )                                        :: derivatives
    procedure       (jacobianTemplate           )               , optional               :: jacobian
    procedure       (integrandTemplate          )               , optional               :: integrands
    class           (integratorMultiVectorized1D), intent(in   ), optional, target       :: integrator
    logical                                      , intent(in   ), optional               :: integratorErrorTolerant
    procedure       (finalStateTemplate         )               , optional               :: finalState
    procedure       (postStepTemplate           )               , optional               :: postStep
    procedure       (errorAnalyzerTemplate      )               , optional               :: errorAnalyzer
    procedure       (errorHandlerTemplate       )               , optional               :: errorHandler
    integer                                      , intent(in   ), optional               :: stepperType
    double precision                             , intent(in   ), optional               :: toleranceAbsolute      , toleranceRelative, &
         &                                                                                  yScale                 , dydtScale        , &
         &                                                                                  hStart
    double precision                             , intent(in   ), optional, dimension(:) :: scale
    class           (*                          )                         , pointer      :: dummyPointer_
    logical                                      , intent(in   ), optional, dimension(:) :: isNonNegative
    integer         (c_int                      ), allocatable            , dimension(:) :: is_non_negative
    !![
    <optionalArgument name="stepperType"             defaultsTo="gsl_odeiv2_step_rkck"/>
    <optionalArgument name="toleranceAbsolute"       defaultsTo="0.0d0"               />
    <optionalArgument name="toleranceRelative"       defaultsTo="0.0d0"               />
    <optionalArgument name="yScale"                  defaultsTo="1.0d0"               />
    <optionalArgument name="dydtScale"               defaultsTo="0.0d0"               />
    <optionalArgument name="hStart"                  defaultsTo="1.0d0"               />
    <optionalArgument name="integratorErrorTolerant" defaultsTo=".false."             />
    <constructorAssign variables="dim, *derivatives, *jacobian, *integrator, *integrands, *finalState, *postStep, *errorAnalyzer, *errorHandler"/>
    !!]

    ! Validate.
    if (toleranceAbsolute_ <= 0.0d0 .and. toleranceRelative_ <= 0.0d0) &
         & call Error_Report('at least one of absolute and relative tolerance must be greater than zero'//{introspection:location})
    ! Get the stepper type.
    self   %stepperType=stepperType_
    allocate(self%gsl_odeiv2_step_type)
    self   %gsl_odeiv2_step_type=gsl_odeiv2_step_type_get         (stepperType_)
    ! Allocate and initialize the system object.
    allocate(self%system)
    if (present(jacobian)) then
       self%system%gsl=gsl_odeiv2_system_init(self%dim,c_funloc(derivativesWrapper),c_funloc     (jacobianWrapper))
    else    
       self%system%gsl=gsl_odeiv2_system_init(self%dim,c_funloc(derivativesWrapper),c_null_funptr                 )
    end if
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%system
    self%systemManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    ! Allocate and initialize the driver object.
    allocate(self%driver)
    if (present(scale)) then
       allocate(is_non_negative(size(scale)))
       is_non_negative=0
       if (present(isNonNegative)) then
          where (isNonNegative)
             is_non_negative=1
          end where
       end if
       self%driver%gsl=gsl_odeiv2_driver_alloc_scaled2_new(self%system%gsl,self%gsl_odeiv2_step_type,hStart_,toleranceAbsolute_,toleranceRelative_,yScale_,dydtScale_,scale,is_non_negative)
       call gsl_odeiv2_driver_init_errors(self%driver%gsl)
    else
       self%driver%gsl=gsl_odeiv2_driver_alloc_y_new      (self%system%gsl,self%gsl_odeiv2_step_type,hStart_,toleranceAbsolute_,toleranceRelative_                                         )
    end if
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
      <description>ICE when passing a derived type component to a class(*) function argument.</description>
    !!]
    dummyPointer_      => self%driver
    self%driverManager =  resourceManager(dummyPointer_)
    !![
    </workaround>
    !!]
    ! Set integrator error tolerance behavior.
    self%integratorErrorTolerant=integratorErrorTolerant_
    return
  end function odeSolverConstructor

  subroutine gslODEDriverWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily gslODEDriverWrapper} object.
    !!}
    implicit none
    type(gslODEDriverWrapper), intent(inout) :: self

    call gsl_odeiv2_driver_free(self%gsl)
    return
  end subroutine gslODEDriverWrapperDestructor

  subroutine gslODESystemWrapperDestructor(self)
    !!{
    Destroy a {\normalfont \ttfamily gslODESystemWrapper} object.
    !!}
    implicit none
    type(gslODESystemWrapper), intent(inout) :: self

    call gsl_odeiv2_system_free(self%gsl)
    return
  end subroutine gslODESystemWrapperDestructor

  subroutine odeSolverDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily odeSolver} objects.
    !!}
    implicit none
    type(odeSolver), intent(inout) :: self

    nullify(self%integrator)
    return
  end subroutine odeSolverDestructor

  subroutine odeSolverSolve(self,x0,x1,y,z,xStep,status)
    !!{
    Solve the ODE system.
    !!}
    use            :: Error                 , only : Error_Report
    use, intrinsic :: ISO_C_Binding         , only : c_funloc      , c_null_funptr     , c_null_ptr
    use            :: ISO_Varying_String    , only : operator(//)  , var_str
    use            :: Interface_GSL         , only : GSL_Failure   , GSL_Success
    use            :: ODE_Solver_Error_Codes, only : interruptedAtX, odeSolverInterrupt
    use            :: String_Handling       , only : operator(//)
    implicit none
    class           (odeSolver ), intent(inout), target                        :: self
    double precision            , intent(inout)                                :: x0
    double precision            , intent(in   )                                :: x1
    double precision            , intent(inout), dimension(self%dim)           :: y
    double precision            , intent(inout), dimension(:       ), optional :: z
    double precision                           , dimension(self%dim)           :: y0
    double precision            , intent(inout)                     , optional :: xStep
    integer                     , intent(  out)                     , optional :: status
    integer                     , parameter                                    :: solversIncrement =3
    type            (solverList), allocatable  , dimension(:       )           :: solversTmp
    double precision            , allocatable  , dimension(:       )           :: z0
    integer                                                                    :: zCount
    type            (c_funptr  )                                               :: latentIntegrator_  , errorAnalyzer_, &
         &                                                                        postStep_
    double precision                                                           :: xStep_             , x             , &
         &                                                                        x1_                , x0Step
    logical                                                                    :: evolveForward
    integer         (c_int     )                                               :: status_

    ! Add the current solver to the list of solvers. This allows us to track back to the previously used solver if this function
    ! is called recursively.
    active=active+1
    if (allocated(solvers)) then
       if (size(solvers) < active) then
          call move_alloc(solvers,solversTmp)
          allocate(solvers(size(solversTmp)+solversIncrement))
          solvers(1:size(solversTmp))=solversTmp
          deallocate(solversTmp)
       end if
    else
       allocate(solvers(solversIncrement))
    end if
    solvers(active)%solver => self
    ! Make initial guess for timestep.
    xStep_=(x1-x0)
    if (present(xStep)) then
       if (xStep > 0.0d0) xStep_=min(xStep,xStep_)
    end if
    ! Keep a local copy of the initial y values so that we can repeat the step if necessary.
    y0 =y    
    ! Keep a local copy of the end point as we may reset it.
    x1_=x1
    ! Set initial value of x variable.
    x  =x0
    ! Determine if we want forward or backward evolution.
    evolveForward=x1 > x0
    ! Reset the driver.
    status_   =GSL_ODEIV2_Driver_Reset       (self%driver%gsl       )
    if    (status_ /= GSL_Success) call Error_Report('failed to reset ODE driver'    //{introspection:location})
    if (xStep_ /= 0.0d0) then
       status_=GSL_ODEIV2_Driver_Reset_hStart(self%driver%gsl,xStep_)
       if (status_ /= GSL_Success) call Error_Report('failed to reset ODE step size'//{introspection:location})
    end if
    ! Initialize integrator.
    if (present(z)) then
       zCount=size(z)
       allocate(z0(zCount))
       z0=z
       call self%integrator%integrandSet(zCount,integrandsWrapper)
       latentIntegrator_=C_FunLoc(latentIntegrator)
    else
       allocate(z0(0))
       latentIntegrator_=C_NULL_FUNPTR 
    end if
    ! Initialize error analyzer.
    if (associated(self%errorAnalyzer)) then
       errorAnalyzer_=c_funloc(self%errorAnalyzer)
    else
       errorAnalyzer_=c_null_funptr
    end if
    ! Initialize post-step processor.
    if (associated(self%postStep)) then
       postStep_=c_funloc(self%postStep)
    else
       postStep_=c_null_funptr
    end if
    ! Evolve the system until the final time is reached.
    do while (                                  &
         &          evolveForward .and. x < x1_ &
         &    .or.                              &
         &     .not.evolveForward .and. x > x1_ &
         &   )
       ! Store current time.
       x0Step=x
       ! Apply the ODE solver.
       status_=GSL_ODEIV2_Driver2_Apply(self%driver%gsl,x,x1_,y,postStep_,latentIntegrator_,errorAnalyzer_)
       select case (status_)
       case (GSL_Success)
          ! Successful completion of the step - do nothing except store the step-size used.
          if (present(xStep)) xStep=GSL_ODEIV2_Driver_h(self%driver%gsl)
       case (GSL_Failure)
          ! Generic failure - most likely a stepsize underflow.
          if (associated(self%errorHandler)) then
             xStep_=GSL_ODEIV2_Driver_h(self%gsl_odeiv2_driver)
             call self%errorHandler(status_,x,xStep_,y)
          end if
          ! If ODE status was requested, then return it instead of aborting.
          if (present(status)) then
             x0    =x
             status=status_
             ! Restore state.
             active=active-1
             return
          end if
          call Error_Report(var_str('ODE integration failed with status ')//status_//' [generic failure] => most likely a stepsize underflow'//{introspection:location})
       case (odeSolverInterrupt)
          ! The evolution was interrupted. Reset the end time of the evolution and continue.
          x1_=interruptedAtX
          if (x > x1_) then
             ! The timestep exceeded the time at which an interrupt occurred. To maintain accuracy we need to repeat the step.
             y=y0
             x=x0
             if (present(z)) z=z0
             status_=GSL_ODEIV2_Driver_Reset(self%driver%gsl)
             if (status_ /= GSL_Success) call Error_Report('failed to reset ODE driver'//{introspection:location})
          end if
       case default
          ! If ODE status was requested, then return it instead of aborting.
          if (present(status)) then
             x0    =x
             status=status_
             ! Restore state.
             active=active-1
             return
          end if
          call Error_Report(var_str('ODE integration failed with status ')//status//{introspection:location})
       end select
    end do
    ! Return the new value of x.
    x0=x
    if (present(status)) status=status_
    ! Move back to the previous active solver.
    active=active-1
    return

  contains

    subroutine latentIntegrator(x)
      !!{
      Wrapper function which performs integration of latent variables.
      !!}
      use :: Display, only : displayMessage
      use :: Error  , only : Error_Report  , errorStatusSuccess
      implicit none
      double precision, intent(in   )       :: x
      double precision, dimension(self%dim) :: y
      integer                               :: status
      
      ! Call with the final state.
      if (associated(solvers(active)%solver%finalState)) then
         call MSBDFActive_State(self%driver%gsl,solvers(active)%solver%dim,y)
         call solvers(active)%solver%finalState(x,y)
      end if
      ! Evaluate the integrals, and update the stored time ready for the next step.
      z      =+z                                         &
           &  +self%integrator%evaluate(x0Step,x,status)
      x0Step =+                                x
      if (status /= errorStatusSuccess) then
         if (self%integratorErrorTolerant) then
            call displayMessage('integration of latent variables failed - ignoring'                          )
         else
            call Error_Report  ('integration of latent variables failed'           //{introspection:location})
         end if
      end if
      return
    end subroutine latentIntegrator
    
    subroutine integrandsWrapper(nz,x,e,dzdx)
      !!{
      Wrapper function which calls the integrands functions.
      !!}
      implicit none
      integer         , intent(in   )                              :: nz
      double precision, intent(in   ), dimension(              : ) :: x
      logical         , intent(inout), dimension(              : ) :: e
      double precision, intent(  out), dimension(nz      ,size(x)) :: dzdx
      double precision               , dimension(self%dim,size(x)) :: dydx, y
      integer                                                      :: i

      ! Evaluate the active parameters.
      do i=1,size(x)
         call msbdfactive_context(self%driver%gsl,self%dim,x(i),y(:,i),dydx(:,i))
      end do
      ! Call the integrand function.
      call self%integrands(x,y,dydx,z,e,dzdx)
      return
    end subroutine integrandsWrapper

  end subroutine odeSolverSolve
  
  function derivativesWrapper(x,y,dydx,parameterPointer) bind(c)
    !!{
    Wrapper function used for \gls{gsl} ODEIV2 derivative functions.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_double, c_int, c_ptr
    implicit none
    integer(c_int   )                              :: derivativesWrapper
    real   (c_double), value                       :: x
    real   (c_double), dimension(*), intent(in   ) :: y
    real   (c_double), dimension(*)                :: dydx
    type   (c_ptr   ), value                       :: parameterPointer
    !$GLC attributes unused :: parameterPointer
    
    derivativesWrapper=solvers(active)%solver%derivatives(                                       &
         &                                                x                                    , &
         &                                                y   (1:solvers(active)%solver%dim   ), &
         &                                                dydx(1:solvers(active)%solver%dim   )  &
         &                                               )
    return
  end function derivativesWrapper

  function jacobianWrapper(x,y,dfdy,dfdx,parameterPointer) bind(c)
    !!{
    Wrapper function used for \gls{gsl} ODEIV2 Jacobian functions.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_double, c_int, c_ptr
    implicit none
    integer(c_int   )                              :: jacobianWrapper
    real   (c_double), value                       :: x
    real   (c_double), dimension(*), intent(in   ) :: y
    real   (c_double), dimension(*)                :: dfdy              , dfdx
    type   (c_ptr   ), value                       :: parameterPointer
    !$GLC attributes unused :: parameterPointer
    
    jacobianWrapper=solvers(active)%solver%jacobian      (                                       &
         &                                                x                                    , &
         &                                                y   (1:solvers(active)%solver%dim   ), &
         &                                                dfdy(1:solvers(active)%solver%dim**2), &
         &                                                dfdx(1:solvers(active)%solver%dim   )  &
         &                                               )
    return
  end function jacobianWrapper

  subroutine odeSolverErrors(self,yError)
    !!{
    Return estimates of the errors in ODE properties.
    !!}
    implicit none
    class           (odeSolver), intent(inout)               :: self
    double precision           , intent(  out), dimension(:) :: yError

    call gsl_odeiv2_driver_errors(self%driver%gsl,yError)
    return
  end subroutine odeSolverErrors
  
end module Numerical_ODE_Solvers

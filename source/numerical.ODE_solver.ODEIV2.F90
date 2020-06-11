!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which provides an interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library} \href{http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html}{ODEIV2} differential equation solvers.

module ODEIV2_Solver
  !% Contains an interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library} \href{http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html}{ODEIV2} differential equation solvers.
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  private
  public :: ODEIV2_Solve, ODEIV2_Solver_Free

  ! ODE interface.
  abstract interface
     integer function odesTemplate(x,y,dydx)
       double precision, intent(in   )               :: x
       double precision, intent(in   ), dimension(:) :: y
       double precision, intent(  out), dimension(:) :: dydx
     end function odesTemplate
  end interface

  ! Jacobian interface.
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
     subroutine integrandTemplate(ny,nz,x,y,dydx,z0,e,dzdx)
       integer         , intent(in   )                        :: ny  , nz
       double precision, intent(in   ), dimension(        : ) :: x
       double precision, intent(in   ), dimension(ny,size(x)) :: y   , dydx
       double precision, intent(in   ), dimension(nz        ) :: z0
       logical         , intent(inout), dimension(        : ) :: e
       double precision, intent(  out), dimension(nz,size(x)) :: dzdx
     end subroutine integrandTemplate
  end interface

  ! Final state interface.
  abstract interface
     subroutine finalStateTemplate(ny,y)
       integer         , intent(in   )                :: ny
       double precision, intent(in   ), dimension(ny) :: y
     end subroutine finalStateTemplate
  end interface

  type :: odeiv2ODEsList
     !% Type used to maintain a list of ODEs when ODE solving is performed recursively.
     procedure(      odesTemplate), pointer, nopass :: ODEs
     procedure(  jacobianTemplate), pointer, nopass :: jacobian
     procedure( integrandTemplate), pointer, nopass :: integrands
     procedure(finalStateTemplate), pointer, nopass :: finalState
     integer  (c_size_t          )                  :: ODENumber , integrandsNumber
  end type odeiv2ODEsList

  ! List of currently active root ODE systems.
  integer                                            :: currentODEsIndex=0
  type   (odeiv2ODEsList), allocatable, dimension(:) :: currentODEs
  !$omp threadprivate(currentODEs,currentODEsIndex)

contains

  subroutine ODEIV2_Solve(                                                                             &
       &                  odeDriver,odeSystem,x0,x1,yCount,y,odes,toleranceAbsolute,toleranceRelative, &
       &                  postStep,Error_Analyzer                                                             , &
       &                  yScale,errorHandler,algorithm,reset,odeStatus,stepSize,jacobian,zCount,z,integrands,finalState,integrator_,integratorErrorTolerate  &
       &                 )
    !% Interface to the \href{http://www.gnu.org/software/gsl/}{GNU Scientific Library} \href{http://www.gnu.org/software/gsl/manual/html_node/Ordinary-Differential-Equations.html}{ODEIV2} differential equation solvers.
    use            :: Interface_GSL         , only : GSL_Failure                    , GSL_Success
    use            :: FODEIV2               , only : FODEIV2_Driver_Alloc_Scaled_New, FODEIV2_Driver_Alloc_y_New, FODEIV2_Driver_Apply, FODEIV2_Driver_Reset, &
          &                                          FODEIV2_Driver_Status          , FODEIV2_Driver_h          , FODEIV2_System_Init , Fodeiv2_Step_RKCK   , &
          &                                          fodeiv2_driver                 , fodeiv2_step_type         , fodeiv2_system
    use             :: Galacticus_Error      ,only : Galacticus_Error_Report
    use            :: ISO_Varying_String    , only : assignment(=)                  , operator(//)              , varying_string
    use, intrinsic :: ISO_C_Binding         , only : C_Null_FunPtr                  , C_FunPtr                  , C_FunLoc            , C_Ptr
    use            :: Numerical_Integration2, only : integratorMultiVectorized1D
    use            :: ODE_Solver_Error_Codes, only : interruptedAtX                 , odeSolverInterrupt
    use            :: String_Handling       , only : operator(//)
    implicit none
    double precision                              , intent(in   )                         :: toleranceAbsolute        , toleranceRelative        , x1
    integer                                       , intent(in   )                         :: yCount
    double precision                              , intent(inout)                         :: x0                       , y                (yCount)
    double precision                              , intent(in   ), optional               :: yScale           (yCount)
    type            (fodeiv2_driver              ), intent(inout)                         :: odeDriver
    type            (fodeiv2_system              ), intent(inout)                         :: odeSystem
    logical                                       , intent(inout), optional               :: reset
    procedure       (                            )               , optional, pointer      :: errorHandler
    type            (fodeiv2_step_type           ), intent(in   ), optional               :: algorithm
    integer                                       , intent(  out), optional               :: odeStatus
    double precision                              , intent(inout), optional               :: stepSize
    type            (c_funptr                    ), intent(in   ), optional               :: postStep
    type            (c_funptr                    ), intent(in   ), optional               :: Error_Analyzer
    procedure       (odesTemplate                )                                        :: odes
    procedure       (jacobianTemplate            ), optional                              :: jacobian
    integer                                       , parameter                             :: genericFailureCountMaximum=10
    procedure       (integrandTemplate           ), optional                              :: integrands
    procedure       (finalStateTemplate          ), optional                              :: finalState
    integer                                       , intent(in   ), optional               :: zCount
    double precision                              , intent(inout), optional, dimension(:) :: z
    class            (integratorMultiVectorized1D), intent(inout), optional               :: integrator_
    logical                                       , intent(in   ), optional               :: integratorErrorTolerate
    double precision                              , dimension(:) , allocatable            :: z0
    type            (odeiv2ODEsList              ), dimension(:), allocatable             :: currentODEsTmp
    double precision                                                                      :: y0(yCount)
    double precision                              , parameter                             :: dydtScaleUniform          =0.0d0, yScaleUniform=1.0d0
    integer                                       , parameter                             :: odesIncrement             =3
    integer                                                                               :: status
    double precision                                                                      :: h                               , x                       , &
         &                                                                                   x1Internal                      , xStepBegin
    logical                                                                               :: forwardEvolve                   , resetActual             , &
         &                                                                                   doReset
    type            (fodeiv2_step_type           )                                        :: algorithmActual
    type            (varying_string              )                                        :: message
    type            (c_ptr                       )                                        :: parameterPointer
    type            (c_funptr                    )                                        :: latentIntegrator_               , postStep_               , &
         &                                                                                   Error_Analyzer_
    !# <optionalArgument name="integratorErrorTolerate" defaultsTo=".false." />
    !# <optionalArgument name="zCount"                  defaultsTo="0"       />

    ! Add the current finder to the list of finders. This allows us to track back to the previously used finder if this function is called recursively.
    currentODEsIndex=currentODEsIndex+1
    if (allocated(currentODEs)) then
       if (size(currentODEs) < currentODEsIndex) then
          call move_alloc(currentODEs,currentODEsTmp)
          allocate(currentODEs(size(currentODEsTmp)+odesIncrement))
          currentODEs(1:size(currentODEsTmp))=currentODEsTmp
          deallocate(currentODEsTmp)
       end if
    else
       allocate(currentODEs(odesIncrement))
    end if
    currentODEs(currentODEsIndex)%ODEs      => ODEs
    currentODEs(currentODEsIndex)%ODENumber =  yCount
    if (present(jacobian)) then
       currentODEs(currentODEsIndex)%jacobian => jacobian
    else
       currentODEs(currentODEsIndex)%jacobian => null()
    end if
    if (present(integrands)) then
       currentODEs(currentODEsIndex)%integrands       => integrands
       currentODEs(currentODEsIndex)%integrandsNumber =  zCount
    else
       currentODEs(currentODEsIndex)%integrands       => null()
       currentODEs(currentODEsIndex)%integrandsNumber =  0
    end if
    if (present(finalState)) then
       currentODEs(currentODEsIndex)%finalState       => finalState
    else
       currentODEs(currentODEsIndex)%finalState       => null()
    end if
    ! Set initial values of integrands.
    if (present(integrands)) then
       allocate(z0(zCount))
       z0=z
    end if
    ! Initialize integrator if required.
    if (zCount_ > 0) call integrator_%integrandSet (zCount,integrandsWrapper)
    ! Decide whether to reset.
    resetActual=.false.
    if (present(reset)) then
       resetActual=reset
       reset=.false.
    end if
    if (resetActual) then
       doReset=.true.
    else if (.not.FODEIV2_Driver_Status(odeDriver)) then
       doReset=.true.
    else
       doReset=.false.
    end if
    if (doReset) then
       ! Make initial guess for timestep.
       h=(x1-x0)
       if (present(stepSize)) then
          if (stepSize > 0.0d0) h=min(stepSize,h)
       end if
       if (present(jacobian)) then
          odeSystem=FODEIV2_System_Init(odesWrapperIV2,currentODEs(currentODEsIndex)%ODENumber,parameterPointer,jacobianWrapperIV2)
       else
          odeSystem=FODEIV2_System_Init(odesWrapperIV2,currentODEs(currentODEsIndex)%ODENumber,parameterPointer                   )
       end if
       ! Select the algorithm to use.
       if (present(algorithm)) then
          algorithmActual=algorithm
       else
          algorithmActual=Fodeiv2_Step_RKCK
       end if
       if (present(yScale)) then
          ! Scales for the absolute tolerance have been given, so use them.
          odeDriver=FODEIV2_Driver_Alloc_Scaled_New(odeSystem,algorithmActual,h,toleranceAbsolute,toleranceRelative&
               &,yScaleUniform,dydtScaleUniform,yScale)
       else
          ! No scales given, assume they are all unity.
          odeDriver=FODEIV2_Driver_Alloc_y_New     (odeSystem,algorithmActual,h,toleranceAbsolute,toleranceRelative)
       end if
    end if
    ! Keep a local copy of the initial y values so that we can repeat the step if necessary.
    y0=y
    ! Keep a local copy of the end point as we may reset it.
    x1Internal=x1
    ! Set initial value of x variable.
    x=x0
    ! Determine if we want forward or backward evolution.
    forwardEvolve=x1>x0
    ! Reset the driver.
    status=FODEIV2_Driver_Reset(odeDriver)
    ! Get a C-pointer to our latent integrator.
    if (zCount_ > 0) then
       latentIntegrator_=C_FunLoc(latentIntegrator)
    else
       latentIntegrator_=C_NULL_FUNPTR
    end if
    ! Evolve the system until the final time is reached.
    do while ((forwardEvolve.and.x<x1Internal).or.(.not.forwardEvolve.and.x>x1Internal))
       ! Store current time.
       xStepBegin=x
       if (present(Error_Analyzer)) then
          Error_Analyzer_=Error_Analyzer
       else
          Error_Analyzer_=C_NULL_FUNPTR
       end if
       if (present(postStep)) then
          postStep_=postStep
       else
          postStep_=C_NULL_FUNPTR
       end if
       status=FODEIV2_Driver_Apply(odeDriver,x,x1Internal,y,postStep_,latentIntegrator_,Error_Analyzer_)
       select case (status)
       case (GSL_Success)
          ! Successful completion of the step - do nothing except store the step-size used.
          if (present(stepSize)) stepSize=FODEIV2_Driver_h(odeDriver)
       case (GSL_Failure)
          ! Generic failure - most likely a stepsize underflow.
          if (present(errorHandler)) call errorHandler(status,x,y)
          ! If ODE status was requested, then return it instead of aborting.
          if (present(odeStatus)) then
             x0=x
             odeStatus=status
             ! Restore state.
             currentODEsIndex=currentODEsIndex-1
             return
          end if
          message='ODE integration failed with status '
          message=message//status//' [generic failure]'//char(10)
          message=message//' => most likely a stepsize underflow'
          call Galacticus_Error_Report(message//{introspection:location})
       case (odeSolverInterrupt)
          ! The evolution was interrupted. Reset the end time of the evolution and continue.
          x1Internal=interruptedAtX
          if (x > x1Internal) then
             ! The timestep exceeded the time at which an interrupt occured. To maintain accuracy we need to repeat the step.
             y=y0
             x=x0
             if (present(z)) z=z0
             status=FODEIV2_Driver_Reset(odeDriver)
          end if
       case default
          ! Some other error condition.
          if (present(errorHandler)) call errorHandler(status,x,y)
          ! If ODE status was requested, then return it instead of aborting.
          if (present(odeStatus)) then
             x0=x
             odeStatus=status
             ! Restore state.
             currentODEsIndex=currentODEsIndex-1
             return
          end if
          message='ODE integration failed with status '
          message=message//status
          call Galacticus_Error_Report(message//{introspection:location})
       end select
    end do
    ! Return the new value of x.
    x0=x
    if (present(odeStatus)) odeStatus=status
    ! Restore state.
    currentODEsIndex=currentODEsIndex-1
    return

  contains

    subroutine latentIntegrator(x)
      !% Wrapper function which performs integration of latent variables.
      use :: FODEIV2           , only : FODEIV2_Driver_MSBDFActive_State
      use :: Galacticus_Display, only : Galacticus_Display_Message
      use :: Galacticus_Error  , only : Galacticus_Error_Report         , errorStatusSuccess
      implicit none
      double precision, intent(in   )     :: x
      double precision, dimension(yCount) :: y
      integer                             :: status

      ! Call with the final state.
      if (associated(currentODEs(currentODEsIndex)%finalState)) then
         call FODEIV2_Driver_MSBDFActive_State(odeDriver,currentODEs(currentODEsIndex)%ODENumber,y)
         call currentODEs(currentODEsIndex)%finalState(yCount,y)
      end if
      ! Evaluate the integrals, and update the stored time ready for the next step.
      z         =+z                                         &
           &     +integrator_%evaluate(xStepBegin,x,status)
      xStepBegin=+                                x
      if (status /= errorStatusSuccess) then
         if (integratorErrorTolerate_) then
            call Galacticus_Display_Message('integration of latent variables failed - ignoring'                          )
         else
            call Galacticus_Error_Report   ('integration of latent variables failed'           //{introspection:location})
         end if
      end if
      return
    end subroutine latentIntegrator

    subroutine integrandsWrapper(nz,x,e,dzdx)
      !% Wrapper function which calls the integrands functions.
      use :: FODEIV2, only : FODEIV2_Driver_MSBDFActive_Context
      implicit none
      integer         , intent(in   )                            :: nz
      double precision, intent(in   ), dimension(            : ) :: x
      logical         , intent(inout), dimension(            : ) :: e
      double precision, intent(  out), dimension(nz    ,size(x)) :: dzdx
      double precision               , dimension(yCount,size(x)) :: dydx, y
      integer                                                    :: i

      ! Evaluate the active parameters.
      do i=1,size(x)
         call FODEIV2_Driver_MSBDFActive_Context(odeDriver,currentODEs(currentODEsIndex)%ODENumber,x(i),y(:,i),dydx(:,i))
      end do
      ! Call the integrand function.
      call currentODEs(currentODEsIndex)%integrands(yCount,nz,x,y,dydx,z,e,dzdx)
      return
    end subroutine integrandsWrapper

  end subroutine ODEIV2_Solve

  function odesWrapperIV2(x,y,dydx,parameterPointer) bind(c)
    !% Wrapper function used for \gls{gsl} ODEIV2 functions.
    use, intrinsic :: ISO_C_Binding, only : c_double, c_int, c_ptr
    implicit none
    integer(kind=c_int   )                              :: odesWrapperIV2
    real   (kind=c_double), value                       :: x
    real   (kind=c_double), dimension(*), intent(in   ) :: y
    real   (kind=c_double), dimension(*)                :: dydx
    type   (     c_ptr   ), value                       :: parameterPointer
    !$GLC attributes unused :: parameterPointer

    odesWrapperIV2=currentODEs(currentODEsIndex)%ODes(x,y(1:currentODEs(currentODEsIndex)%ODENumber),dydx(1:currentODEs(currentODEsIndex)%ODENumber))
    return
  end function odesWrapperIV2

  function jacobianWrapperIV2(x,y,dfdy,dfdx,parameterPointer) bind(c)
    !% Wrapper function used for \gls{gsl} ODEIV2 Jacobian functions.
    use, intrinsic :: ISO_C_Binding, only : c_double, c_int, c_ptr
    implicit none
    integer(kind=c_int   )                              :: jacobianWrapperIV2
    real   (kind=c_double), value                       :: x
    real   (kind=c_double), dimension(*), intent(in   ) :: y
    real   (kind=c_double), dimension(*)                :: dfdy              , dfdx
    type   (     c_ptr   ), value                       :: parameterPointer
    !$GLC attributes unused :: parameterPointer

    jacobianWrapperIV2=currentODEs(currentODEsIndex)%jacobian(x,y(1:currentODEs(currentODEsIndex)%ODENumber),dfdy(1:currentODEs(currentODEsIndex)%ODENumber**2),dfdx(1:currentODEs(currentODEsIndex)%ODENumber))
    return
  end function jacobianWrapperIV2

  subroutine ODEIV2_Solver_Free(odeDriver,odeSystem)
    !% Free up workspace allocated to ODE solving.
    use :: FODEIV2, only : Fodeiv2_Driver_Free, Fodeiv2_System_Free, fodeiv2_driver, fodeiv2_system
    implicit none
    type(fodeiv2_driver), intent(inout) :: odeDriver
    type(fodeiv2_system), intent(inout) :: odeSystem

    call Fodeiv2_Driver_Free(odeDriver)
    call Fodeiv2_System_Free(odeSystem)
    return
  end subroutine ODEIV2_Solver_Free

end module ODEIV2_Solver

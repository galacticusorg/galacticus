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

!% Contains a module which does root finding.

module Root_Finder
  !% Implements root finding.
  use, intrinsic :: ISO_C_Binding
  use               FGSL
  implicit none
  private
  public :: rootFinder

  ! Enumeration of range expansion types.
  !@ <enumeration>
  !@  <name>rangeExpand</name>
  !@  <description>Used to specify the way in which the bracketing range should be expanded when searching for roots using a {\tt rootFinder} object.</description>
  !@  <entry label="rangeExpandNull"           />
  !@  <entry label="rangeExpandAdditive"       />
  !@  <entry label="rangeExpandMultiplicative" />
  !@ </enumeration>
  integer, parameter, public :: rangeExpandNull              =0
  integer, parameter, public :: rangeExpandAdditive          =1
  integer, parameter, public :: rangeExpandMultiplicative    =2

  ! Enumeration of sign expectations.
  !@ <enumeration>
  !@  <name>rangeExpandSignExpect</name>
  !@  <description>Used to specify the expected sign of the root function when searching for roots using a {\tt rootFinder} object.</description>
  !@  <entry label="rangeExpandSignExpectNegative" />
  !@  <entry label="rangeExpandSignExpectNone"     />
  !@  <entry label="rangeExpandSignExpectPositive" />
  !@ </enumeration>
  integer, parameter, public :: rangeExpandSignExpectNegative=-1
  integer, parameter, public :: rangeExpandSignExpectNone    =0
  integer, parameter, public :: rangeExpandSignExpectPositive=+1

  type :: rootFinder
     !% Type containing all objects required when calling the FGSL root solver function.
     private
     type            (fgsl_function                 )                  :: fgslFunction
     type            (fgsl_function_fdf             )                  :: fgslFunctionDerivative
     type            (fgsl_root_fsolver             )                  :: solver
     type            (fgsl_root_fdfsolver           )                  :: solverDerivative
     type            (fgsl_root_fsolver_type        )                  :: solverType                   =FGSL_Root_fSolver_Brent
     type            (fgsl_root_fdfsolver_type      )                  :: solverDerivativeType         =FGSL_Root_fdfSolver_Steffenson
     double precision                                                  :: toleranceAbsolute            =1.0d-10
     double precision                                                  :: toleranceRelative            =1.0d-10
     logical                                                           :: initialized                  =.false.
     logical                                                           :: resetRequired                =.false.
     logical                                                           :: useDerivative
     integer                                                           :: rangeExpandType              =rangeExpandNull
     double precision                                                  :: rangeExpandUpward            =1.0d0
     double precision                                                  :: rangeExpandDownward          =1.0d0
     double precision                                                  :: rangeUpwardLimit
     double precision                                                  :: rangeDownwardLimit
     logical                                                           :: rangeUpwardLimitSet          =.false.
     logical                                                           :: rangeDownwardLimitSet        =.false.
     integer                                                           :: rangeExpandDownwardSignExpect=rangeExpandSignExpectNone
     integer                                                           :: rangeExpandUpwardSignExpect  =rangeExpandSignExpectNone
     procedure       (rootFunctionTemplate          ), nopass, pointer :: finderFunction
     procedure       (rootFunctionDerivativeTemplate), nopass, pointer :: finderFunctionDerivative
     procedure       (rootFunctionBothTemplate      ), nopass, pointer :: finderFunctionBoth
   contains
     !@ <objectMethods>
     !@   <object>rootFinder</object>
     !@   <objectMethod>
     !@     <method>rootFunction</method>
     !@     <description>Set the function that evaluates $f(x)$ to use in a {\tt rootFinder} object.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless function(\textless double\textgreater} x\argin\textcolor{red}{)\textgreater} rootFunction</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>type</method>
     !@     <description>Set the type of algorithm to use in a {\tt rootFinder} object.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(fgsl\_root\_fsolver\_type)\textgreater} solverType\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>tolerance</method>
     !@     <description>Set the tolerance to use in a {\tt rootFinder} object.</description>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ [toleranceAbsolute]\argin, \doublezero\ [toleranceRelative]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>rangeExpand</method>
     !@     <description>Specify how the initial range will be expanded in a {\tt rootFinder} object to bracket the root.</description>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ [rangeExpandUpward]\argin,\doublezero\ [rangeExpandDownward]\argin, \enumRangeExpand\ [rangeExpandType]\argin, \doublezero\ [rangeUpwardLimit]\argin, \doublezero\ [rangeDownwardLimit]\argin, \enumRangeExpandSignExpect\ [rangeExpandDownwardSignExpect]\argin, \enumRangeExpandSignExpect\ [rangeExpandUpwardSignExpect]\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>find</method>
     !@     <description>Find the root of the function given an initial guess or range.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ [rootGuess]|\textcolor{red}{\textless double(2)\textgreater} [rootRange]</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isInitialized</method>
     !@     <description>Return the initialization state of a {\tt rootFinder} object.</description>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                            Root_Finder_Destroy
     procedure :: rootFunction            => Root_Finder_Root_Function
     procedure :: rootFunctionDerivative  => Root_Finder_Root_Function_Derivative
     procedure :: type                    => Root_Finder_Type
     procedure :: typeDerivative          => Root_Finder_Derivative_Type
     procedure :: tolerance               => Root_Finder_Tolerance
     procedure :: rangeExpand             => Root_Finder_Range_Expand
     procedure :: find                    => Root_Finder_Find
     procedure :: isInitialized           => Root_Finder_Is_Initialized
  end type rootFinder

  abstract interface
     double precision function rootFunctionTemplate(x)
       double precision, intent(in   ) :: x
     end function rootFunctionTemplate
  end interface

  abstract interface
     double precision function rootFunctionDerivativeTemplate(x)
       double precision, intent(in   ) :: x
     end function rootFunctionDerivativeTemplate
  end interface

  abstract interface
     subroutine rootFunctionBothTemplate(x,f,df)
       double precision, intent(in   ) :: x
       double precision, intent(  out) :: f, df
     end subroutine rootFunctionBothTemplate
  end interface

  class(rootFinder), pointer :: currentFinder
  !$omp threadprivate(currentFinder)

contains

  subroutine Root_Finder_Destroy(self)
    !% Destroy a root finder object.
    implicit none
    type(rootFinder), intent(inout) :: self
    
    if (FGSL_Well_Defined(self%solver                )) then
       call FGSL_Root_FSolver_Free  (self%solver                )
       call FGSL_Function_Free      (self%fgslFunction          )
    end if
    if (FGSL_Well_Defined(self%solverDerivative      )) then
       call FGSL_Root_FdFSolver_Free(self%solverDerivative      )
       call FGSL_Function_FdF_Free  (self%fgslFunctionDerivative)
    end if
    return
  end subroutine Root_Finder_Destroy

  logical function Root_Finder_Is_Initialized(self)
    !% Return whether a {\tt rootFinder} object is initalized.
    implicit none
    class(rootFinder), intent(in   ) :: self

    Root_Finder_Is_Initialized=self%initialized
    return
  end function Root_Finder_Is_Initialized

  recursive double precision function Root_Finder_Find(self,rootGuess,rootRange,status)
    !% Finds the root of the supplied {\tt root} function.
    use Galacticus_Error
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class           (rootFinder    )              , intent(inout), target   :: self
    real            (kind=c_double )              , intent(in   ), optional :: rootGuess
    real            (kind=c_double ), dimension(2), intent(in   ), optional :: rootRange
    integer                                       , intent(  out), optional :: status
    class           (rootFinder    ), pointer                               :: previousFinder
    integer                         , parameter                             :: iterationMaximum=1000
    logical                                                                 :: rangeChanged         , rangeLowerAsExpected, rangeUpperAsExpected
    integer                                                                 :: iteration            , statusActual
    double precision                                                        :: xHigh                , xLow                , xRoot               , &
         &                                                                     xRootPrevious
    type            (c_ptr         )                                        :: parameterPointer
    type            (varying_string)                                        :: message
    character       (len= 30       )                                        :: label

    ! Store a pointer to the previous rootFinder object. This is necessary as this function can be called recursively, so we must
    ! be able to return state to its original form before exiting the function.
    previousFinder => currentFinder
    currentFinder  => self
    ! Initialize the root finder variables if necessary.
    if (self%useDerivative) then
       if (.not.FGSL_Well_Defined(self%solverDerivative).or.self%resetRequired) then
          self%fgslFunctionDerivative=FGSL_Function_fdf_Init   (                                         &
               &                                                Root_Finder_Wrapper_Function           , &
               &                                                Root_Finder_Wrapper_Function_Derivative, &
               &                                                Root_Finder_Wrapper_Function_Both      , &
               &                                                parameterPointer                         &
               &                                               )
          self%solverDerivative      =FGSL_Root_fdfSolver_Alloc(self%solverDerivativeType)
          self%resetRequired         =.false.
       end if
    else
       if (.not.FGSL_Well_Defined(self%solver).or.self%resetRequired) then
          self%fgslFunction          =FGSL_Function_Init       (Root_Finder_Wrapper_Function,parameterPointer)
          self%solver                =FGSL_Root_fSolver_Alloc  (self%solverType                              )
          self%resetRequired         =.false.
       end if
    end if
    ! Initialize range.
    if      (present(rootRange)) then
       xLow =rootRange(1)
       xHigh=rootRange(2)
    else if (present(rootGuess)) then
       xLow =rootGuess
       xHigh=rootGuess
    else
       call Galacticus_Error_Report('Root_Finder_Find','either "rootGuess" or "rootRange" must be specified')
    end if
    ! Expand the range as necessary.
    if (self%useDerivative) then
       xRoot       =0.5d0*(xLow+xHigh)
       statusActual=FGSL_Root_fdfSolver_Set(self%solverDerivative,self%fgslFunctionDerivative,xRoot)
    else
       do while (self%finderFunction(xLow)*self%finderFunction(xHigh) > 0.0d0)
          rangeChanged=.false.
          select case (self%rangeExpandDownwardSignExpect)
          case (rangeExpandSignExpectNegative)
             rangeLowerAsExpected=(self%finderFunction(xLow ) < 0.0d0)
          case (rangeExpandSignExpectPositive)
             rangeLowerAsExpected=(self%finderFunction(xLow ) > 0.0d0)
          case default
             rangeLowerAsExpected=.false.
          end select
          select case (self%rangeExpandUpwardSignExpect  )
          case (rangeExpandSignExpectNegative)
             rangeUpperAsExpected=(self%finderFunction(xHigh) < 0.0d0)
          case (rangeExpandSignExpectPositive)
             rangeUpperAsExpected=(self%finderFunction(xHigh) > 0.0d0)
          case default
             rangeUpperAsExpected=.false.
          end select
          select case (self%rangeExpandType)
          case (rangeExpandAdditive      )
             if     (                                  &
                  &   self%rangeExpandUpward   > 1.0d0 &
                  &  .and.                             &
                  &  .not.rangeUpperAsExpected         &
                  &  .and.                             &
                  &  (                                 &
                  &   xHigh < self%rangeUpwardLimit    &
                  &   .or.                             &
                  &   .not.self%rangeUpwardLimitSet    &
                  &  )                                 &
                  & ) then
                xHigh=xHigh+self%rangeExpandUpward
                if (self%rangeUpwardLimitSet  ) xHigh=min(xHigh,self%rangeUpwardLimit  )
                rangeChanged=.true.
             end if
             if     (                                  &
                  &   self%rangeExpandDownward < 1.0d0 &
                  &  .and.                             &
                  &  .not.rangeLowerAsExpected         &
                  &  .and.                             &
                  &  (                                 &
                  &   xLow  > self%rangeDownwardLimit  &
                  &   .or.                             &
                  &   .not.self%rangeDownwardLimitSet  &
                  &  )                                 &
                  & ) then
                xLow =xLow +self%rangeExpandDownward
                if (self%rangeDownwardLimitSet) xLow =max(xLow ,self%rangeDownwardLimit)
                rangeChanged=.true.
             end if
          case (rangeExpandMultiplicative)
             if     (                                    &
                  &  (                                   &
                  &   (                                  &
                  &     self%rangeExpandUpward   > 1.0d0 &
                  &    .and.                             &
                  &     xHigh                    > 0.0d0 &
                  &   )                                  &
                  &   .or.                               &
                  &   (                                  &
                  &     self%rangeExpandUpward   < 1.0d0 &
                  &    .and.                             &
                  &     xHigh                    < 0.0d0 &
                  &   )                                  &
                  &  )                                   &
                  &  .and.                               &
                  &  .not.rangeUpperAsExpected           &
                  &  .and.                               &
                  &  (                                   &
                  &   xHigh < self%rangeUpwardLimit      &
                  &   .or.                               &
                  &   .not.self%rangeUpwardLimitSet      &
                  &  )                                   &
                  & ) then
                xHigh=xHigh*self%rangeExpandUpward
                if (self%rangeUpwardLimitSet  ) xHigh=min(xHigh,self%rangeUpwardLimit  )
                rangeChanged=.true.
             end if
             if     (                                    &
                  &  (                                   &
                  &   (                                  &
                  &     self%rangeExpandDownward < 1.0d0 &
                  &    .and.                             &
                  &     xLow                     > 0.0d0 &
                  &   )                                  &
                  &   .or.                               &
                  &   (                                  &
                  &     self%rangeExpandDownward > 1.0d0 &
                  &    .and.                             &
                  &     xLow                     < 0.0d0 &
                  &   )                                  &
                  &  )                                   &
                  &  .and.                               &
                  &  .not.rangeLowerAsExpected           &
                  &  .and.                               &
                  &  (                                   &
                  &   xLow  > self%rangeDownwardLimit    &
                  &   .or.                               &
                  &   .not.self%rangeDownwardLimitSet    &
                  &  )                                   &
                  & ) then
                xLow =xLow *self%rangeExpandDownward
                if (self%rangeDownwardLimitSet) xLow =max(xLow ,self%rangeDownwardLimit)
                rangeChanged=.true.
             end if
          end select
          if (.not.rangeChanged) then
             message='unable to expand range to bracket root'
             write (label,'(e12.6,a1,e12.6)') xLow ,":",self%finderFunction(xLow )
             message=message//char(10)//'xLow :f(xLow )='//trim(label)
             write (label,'(e12.6,a1,e12.6)') xHigh,":",self%finderFunction(xHigh)
             message=message//char(10)//'xHigh:f(xHigh)='//trim(label)
             if (self%rangeExpandDownwardSignExpect /= rangeExpandSignExpectNone) then
                if (rangeLowerAsExpected) then
                   message=message//char(10)//"f(xLow ) has expected sign"
                else
                   message=message//char(10)//"f(xLow ) does not have expected sign"
                end if
             end if
             if (self%rangeExpandUpwardSignExpect   /= rangeExpandSignExpectNone) then
                if (rangeUpperAsExpected) then
                   message=message//char(10)//"f(xHigh) has expected sign"
                else
                   message=message//char(10)//"f(xHigh) does not have expected sign"
                end if
             end if
             if (self%rangeDownwardLimitSet) then
                write (label,'(e12.6)') self%rangeDownwardLimit
                message=message//char(10)//"xLow  > "//trim(label)//" being enforced"
             end if
             if (self%rangeUpwardLimitSet  ) then
                write (label,'(e12.6)') self%rangeUpwardLimit
                message=message//char(10)//"xHigh < "//trim(label)//" being enforced"
             end if
             if (present(status)) then
                call Galacticus_Display_Message(message,verbosityWarn)
                status=errorStatusOutOfRange
                return
             else
                call Galacticus_Error_Report('Root_Finder_Find',message)
             end if
          end if
       end do
       statusActual=FGSL_Root_fSolver_Set(self%solver,self%fgslFunction,xLow,xHigh)
    end if
    ! Find the root.
    if (statusActual /= FGSL_Success) then
       if (present(status)) then
          status=statusActual
          return
       else
          call Galacticus_Error_Report('Root_Finder_Find','failed to initialize solver')
       end if
    end if
    iteration=0
    do
       iteration=iteration+1
       if (self%useDerivative) then
          statusActual=FGSL_Root_fdfSolver_Iterate(self%solverDerivative)
          if (statusActual /= FGSL_Success .or. iteration > iterationMaximum) exit
          xRootPrevious=xRoot
          xRoot        =FGSL_Root_fdfSolver_Root(self%solverDerivative)
          statusActual =FGSL_Root_Test_Delta(xRoot,xRootPrevious,self%toleranceAbsolute,self%toleranceRelative)
       else
          statusActual=FGSL_Root_fSolver_Iterate  (self%solver          )
          if (statusActual /= FGSL_Success .or. iteration > iterationMaximum) exit
          xRoot =FGSL_Root_fSolver_Root  (self%solver)
          xLow  =FGSL_Root_fSolver_x_Lower(self%solver)
          xHigh =FGSL_Root_fSolver_x_Upper(self%solver)
          statusActual=FGSL_Root_Test_Interval(xLow,xHigh,self%toleranceAbsolute,self%toleranceRelative)
       end if 
       if (statusActual == FGSL_Success) exit
    end do
    if (statusActual /= FGSL_Success) then
       if (present(status)) then 
          status=statusActual
          return
       else
          call Galacticus_Error_Report('Root_Finder_Find','failed to find root')
       end if
    end if
    if (present(status)) status=FGSL_Success
    Root_Finder_Find=xRoot
    ! Restore state.
    currentFinder => previousFinder
    return
  end function Root_Finder_Find

  subroutine Root_Finder_Root_Function(self,rootFunction)
    !% Sets the function to use in a {\tt rootFinder} object.
    implicit none
    class    (rootFinder          ), intent(inout) :: self
    procedure(rootFunctionTemplate)                :: rootFunction

    self%finderFunction => rootFunction
    self%initialized    =  .true.
    self%useDerivative  =  .false.
    return
  end subroutine Root_Finder_Root_Function

  subroutine Root_Finder_Root_Function_Derivative(self,rootFunction,rootFunctionDerivative,rootFunctionBoth)
    !% Sets the function to use in a {\tt rootFinder} object.
    implicit none
    class    (rootFinder                    ), intent(inout) :: self
    procedure(rootFunctionTemplate          )                :: rootFunction
    procedure(rootFunctionDerivativeTemplate)                :: rootFunctionDerivative
    procedure(rootFunctionBothTemplate      )                :: rootFunctionBoth

    self%finderFunction           => rootFunction
    self%finderFunctionDerivative => rootFunctionDerivative
    self%finderFunctionBoth       => rootFunctionBoth
    self%initialized              =  .true.
    self%useDerivative            =  .true.
    return
  end subroutine Root_Finder_Root_Function_Derivative

  subroutine Root_Finder_Type(self,solverType)
    !% Sets the type to use in a {\tt rootFinder} object.
    implicit none
    class(rootFinder            ), intent(inout) :: self
    type (fgsl_root_fsolver_type), intent(in   ) :: solverType

    ! Set the solver type and indicate that a reset will be required to update the internal FGSL objects.
    self%solverType   =solverType
    self%resetRequired=.true.
    return
  end subroutine Root_Finder_Type

  subroutine Root_Finder_Derivative_Type(self,solverDerivativeType)
    !% Sets the type to use in a {\tt rootFinder} object.
    implicit none
    class(rootFinder              ), intent(inout) :: self
    type (fgsl_root_fdfsolver_type), intent(in   ) :: solverDerivativeType

    ! Set the solver type and indicate that a reset will be required to update the internal FGSL objects.
    self%solverDerivativeType=solverDerivativeType
    self%resetRequired       =.true.
    return
  end subroutine Root_Finder_Derivative_Type

  subroutine Root_Finder_Tolerance(self,toleranceAbsolute,toleranceRelative)
    !% Sets the tolerances to use in a {\tt rootFinder} object.
    implicit none
    class           (rootFinder), intent(inout)           :: self
    double precision            , intent(in   ), optional :: toleranceAbsolute, toleranceRelative

    if (present(toleranceAbsolute)) self%toleranceAbsolute=toleranceAbsolute
    if (present(toleranceRelative)) self%toleranceRelative=toleranceRelative
    return
  end subroutine Root_Finder_Tolerance

  subroutine Root_Finder_Range_Expand(self,rangeExpandUpward,rangeExpandDownward,rangeExpandType,rangeUpwardLimit,rangeDownwardLimit,rangeExpandDownwardSignExpect,rangeExpandUpwardSignExpect)
    !% Sets the rules for range expansion to use in a {\tt rootFinder} object.
    use Galacticus_Error
    implicit none
    class           (rootFinder), intent(inout)           :: self
    integer                     , intent(in   ), optional :: rangeExpandDownwardSignExpect, rangeExpandType    , &
         &                                                   rangeExpandUpwardSignExpect
    double precision            , intent(in   ), optional :: rangeDownwardLimit           , rangeExpandDownward, &
         &                                                   rangeExpandUpward            , rangeUpwardLimit

    if (present(rangeExpandUpward            )) self%rangeExpandUpward  =rangeExpandUpward
    if (present(rangeExpandDownward          )) self%rangeExpandDownward=rangeExpandDownward
    if (present(rangeExpandType              )) then
       self%rangeExpandType    =rangeExpandType
    else
       self%rangeExpandType    =rangeExpandNull
    end if
    select case (self%rangeExpandType)
    case (rangeExpandAdditive      )
       if (.not.present(rangeExpandUpward  )) self%rangeExpandUpward=0.0d0
       if (.not.present(rangeExpandDownward)) self%rangeExpandUpward=0.0d0
    case (rangeExpandMultiplicative)
       if (.not.present(rangeExpandUpward  )) self%rangeExpandUpward=1.0d0
       if (.not.present(rangeExpandDownward)) self%rangeExpandUpward=1.0d0
    end select
    if (present(rangeUpwardLimit             )) then
       self%rangeUpwardLimit             =rangeUpwardLimit
       self%rangeUpwardLimitSet          =.true.
    else
       self%rangeUpwardLimit             =0.0d0
       self%rangeUpwardLimitSet          =.false.
    end if
    if (present(rangeDownwardLimit           )) then
       self%rangeDownwardLimit           =rangeDownwardLimit
       self%rangeDownwardLimitSet        =.true.
    else
       self%rangeDownwardLimit           =0.0d0
       self%rangeDownwardLimitSet        =.false.
    end if
    if (present(rangeExpandDownwardSignExpect)) then
       self%rangeExpandDownwardSignExpect=rangeExpandDownwardSignExpect
    else
       self%rangeExpandDownwardSignExpect=rangeExpandSignExpectNone
    end if
    if (present(rangeExpandUpwardSignExpect  )) then
       self%rangeExpandUpwardSignExpect  =rangeExpandUpwardSignExpect
    else
       self%rangeExpandUpwardSignExpect  =rangeExpandSignExpectNone
    end if
    return
  end subroutine Root_Finder_Range_Expand

  function Root_Finder_Wrapper_Function(x,parameterPointer) bind(c)
    !% Wrapper function callable by {\tt FGSL} used in root finding.
    implicit none
    real(kind=c_double), value :: x
    type(c_ptr        ), value :: parameterPointer
    real(kind=c_double)        :: Root_Finder_Wrapper_Function

    Root_Finder_Wrapper_Function=currentFinder%finderFunction(x)
    return
  end function Root_Finder_Wrapper_Function

  double precision function Root_Finder_Wrapper_Function_Derivative(x,parameterPointer) bind(c)
    !% Wrapper function callable by {\tt FGSL} used in root finding.
    implicit none
    real(kind=c_double), value :: x
    type(c_ptr        ), value :: parameterPointer

    Root_Finder_Wrapper_Function_Derivative=currentFinder%finderFunctionDerivative(x)
    return
  end function Root_Finder_Wrapper_Function_Derivative

  subroutine Root_Finder_Wrapper_Function_Both(x,parameterPointer,f,df) bind(c)
    !% Wrapper function callable by {\tt FGSL} used in root finding.
    implicit none
    real(kind=c_double), value         :: x
    type(c_ptr        ), value         :: parameterPointer
    real(kind=c_double), intent(  out) :: f               , df

    call currentFinder%finderFunctionBoth(x,f,df)
    return
  end subroutine Root_Finder_Wrapper_Function_Both

end module Root_Finder

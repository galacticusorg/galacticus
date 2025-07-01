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

! Specify an explicit dependence on the interface.GSL.C.root_finding.o object file.
!: $(BUILDPATH)/interface.GSL.C.root_finding.o

! Add dependency on GSL library.
!; gsl

!!{
Contains a module which does root finding.
!!}
module Root_Finder
  !!{
  Implements root finding.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double, c_int, c_null_ptr, c_ptr
  implicit none
  private
  public :: rootFinder

  ! Enumeration of range expansion types.
  !![
  <enumeration>
   <name>rangeExpand</name>
   <description>Used to specify the way in which the bracketing range should be expanded when searching for roots using a {\normalfont \ttfamily rootFinder} object.</description>
   <visibility>public</visibility>
   <entry label="null"           />
   <entry label="additive"       />
   <entry label="multiplicative" />
  </enumeration>
  !!]

  ! Enumeration of sign expectations.
  !![
  <enumeration>
   <name>rangeExpandSignExpect</name>
   <description>Used to specify the expected sign of the root function when searching for roots using a {\normalfont \ttfamily rootFinder} object.</description>
   <entry label="negative" />
   <entry label="none"     />
   <entry label="positive" />
  </enumeration>
  !!]
  
  ! Enumeration of stopping criteria.
  !![
  <enumeration>
   <name>stoppingCriterion</name>
   <description>Used to specify the stopping criterion to use when searching for roots using a {\normalfont \ttfamily rootFinder} object.</description>
   <visibility>public</visibility>
   <entry label="delta"    />
   <entry label="interval" />
  </enumeration>
  !!]
  
  !![
  <deepCopyActions class="rootFinder">
   <rootFinder>
    <setTo variables="functionInitialized" state=".false."/>
   </rootFinder>
  </deepCopyActions>
  !!]
  
  ! Solver types.
  integer, public, parameter :: gsl_root_fsolver_bisection   =1
  integer, public, parameter :: gsl_root_fsolver_brent       =2
  integer, public, parameter :: gsl_root_fsolver_falsepos    =3
  integer, public, parameter :: gsl_root_fdfsolver_newton    =4
  integer, public, parameter :: gsl_root_fdfsolver_secant    =5
  integer, public, parameter :: gsl_root_fdfsolver_steffenson=6

  type :: rootFinder
     !!{
     Type containing all objects required when calling the GSL root solver function.
     !!}
     private
     type            (c_ptr                               )                  :: gslFunction                  =c_null_ptr
     type            (c_ptr                               )                  :: solver                       =c_null_ptr
     type            (c_ptr                               )                  :: solverType                   =c_null_ptr
     integer                                                                 :: solverTypeID
     double precision                                                        :: toleranceAbsolute
     double precision                                                        :: toleranceRelative
     logical                                                                 :: initialized
     logical                                                                 :: functionInitialized          =.false.
     logical                                                                 :: resetRequired
     logical                                                                 :: useDerivative
     logical                                                                 :: testLimits
     type            (enumerationStoppingCriterionType    )                  :: stoppingCriterion
     type            (enumerationRangeExpandType          )                  :: rangeExpandType
     double precision                                                        :: rangeExpandUpward
     double precision                                                        :: rangeExpandDownward
     double precision                                                        :: rangeUpwardLimit
     double precision                                                        :: rangeDownwardLimit
     logical                                                                 :: rangeUpwardLimitSet
     logical                                                                 :: rangeDownwardLimitSet
     type            (enumerationRangeExpandSignExpectType)                  :: rangeExpandUpwardSignExpect  
     type            (enumerationRangeExpandSignExpectType)                  :: rangeExpandDownwardSignExpect
     procedure       (rootFunctionTemplate                ), nopass, pointer :: finderFunction
     procedure       (rootFunctionDerivativeTemplate      ), nopass, pointer :: finderFunctionDerivative
     procedure       (rootFunctionBothTemplate            ), nopass, pointer :: finderFunctionBoth
   contains
     !![
     <methods>
       <method description="Set the function that evaluates $f(x)$ to use in a {\normalfont \ttfamily rootFinder} object."                                          method="rootFunction"          />
       <method description="Set the functions that evaluate $f(x)$ and derivatives to use in a {\normalfont \ttfamily rootFinder} object."                          method="rootFunctionDerivative"/>
       <method description="Set the type of algorithm to use in a {\normalfont \ttfamily rootFinder} object."                                                       method="type"                  />
       <method description="Set the tolerance to use in a {\normalfont \ttfamily rootFinder} object."                                                               method="tolerance"             />
       <method description="Specify how the initial range will be expanded in a {\normalfont \ttfamily rootFinder} object to bracket the root."                     method="rangeExpand"           />
       <method description="Wrapper function to find the root of the function given an initial guess or range."                                                     method="find"                  />
       <method description="Wrapper function to find the root of the function given an initial guess or range plus the function value at the low end of the range." method="findWithFLower"        />
       <method description="Wrapper function to find the root of the function given an initial guess or range plus the function value at the low end of the range." method="findWithFUpper"        />
       <method description="Find the root of the function given an initial guess or range."                                                                         method="find_"                 />
       <method description="Return the initialization state of a {\normalfont \ttfamily rootFinder} object."                                                        method="isInitialized"         />
       <method description="Destroy the {\normalfont \ttfamily rootFinder} object."                                                                                 method="destroy"               />
       <method description="Return true if the solver type is valid."                                                                                               method="solverTypeIsValid"     />
     </methods>
     !!]
     final     ::                            rootFinderDestructor
     procedure :: destroy                 => rootFinderDestroy
     procedure :: rootFunction            => rootFinderRootFunction
     procedure :: rootFunctionDerivative  => rootFinderRootFunctionDerivative
     procedure :: type                    => rootFinderType
     procedure :: tolerance               => rootFinderTolerance
     procedure :: rangeExpand             => rootFinderRangeExpand
     procedure ::                            rootFinderFindGuess
     procedure ::                            rootFinderFindRange
     procedure ::                            rootFinderFindRangeValues
     generic   :: find                    => rootFinderFindGuess             , rootFinderFindRange, &
          &                                  rootFinderFindRangeValues
     procedure :: findWithFLower          => rootFinderFindRangeValueLow
     procedure :: findWithFUpper          => rootFinderFindRangeValueHigh
     procedure :: find_                   => rootFinderFind
     procedure :: isInitialized           => rootFinderIsInitialized
     procedure :: solverTypeIsValid       => rootFinderSolverTypeIsValid
  end type rootFinder

  interface rootFinder
     !!{
     Interface to constructors for root finders.
     !!}
     module procedure rootFinderConstructorInternal
  end interface rootFinder
  
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

  type :: rootFinderList
     !!{
     Type used to maintain a list of root finder objects when root finding is performed recursively.
     !!}
     class           (rootFinder), pointer :: finder         => null()
     logical                               :: lowInitialUsed          , highInitialUsed
     double precision                      :: xLowInitial             , xHighInitial   , &
          &                                   fLowInitial             , fHighInitial
  end type rootFinderList

  ! List of currently active root finders.
  integer                                            :: currentFinderIndex=0
  type   (rootFinderList), allocatable, dimension(:) :: currentFinders
  !$omp threadprivate(currentFinders,currentFinderIndex)

  interface
     function gsl_root_fsolver_alloc(T) bind(c,name='gsl_root_fsolver_alloc')
       !!{
       Template for the GSL root solver allocate function.
       !!}
       import
       type(c_ptr)        :: gsl_root_fsolver_alloc
       type(c_ptr), value :: T
     end function gsl_root_fsolver_alloc

     function gsl_root_fdfsolver_alloc(T) bind(c,name='gsl_root_fdfsolver_alloc')
       !!{
       Template for the GSL root solver allocate function.
       !!}
       import
       type(c_ptr)        :: gsl_root_fdfsolver_alloc
       type(c_ptr), value :: T
     end function gsl_root_fdfsolver_alloc

     subroutine gsl_root_fsolver_free(s) bind(c,name='gsl_root_fsolver_free')
       !!{
       Template for the GSL root solver free function.
       !!}
       import
       type(c_ptr), value :: s
     end subroutine gsl_root_fsolver_free

     subroutine gsl_root_fdfsolver_free(s) bind(c,name='gsl_root_fdfsolver_free')
       !!{
       Template for the GSL root solver free function.
       !!}
       import
       type(c_ptr), value :: s
     end subroutine gsl_root_fdfsolver_free
     
     integer(c_int) function gsl_root_fsolver_set(s,f,x_lower,x_upper) bind(c,name='gsl_root_fsolver_set')
       !!{
       Template for the GSL root solver set function.
       !!}
       import
       type(c_ptr   ), value :: s      , f
       real(c_double), value :: x_lower, x_upper
     end function gsl_root_fsolver_set

     integer(c_int) function gsl_root_fdfsolver_set(s,fdf,root) bind(c,name='gsl_root_fdfsolver_set')
       !!{
       Template for the GSL root solver set function.
       !!}
       import
       type(c_ptr   ), value :: s    , fdf
       real(c_double), value :: root
     end function gsl_root_fdfsolver_set

     integer(c_int) function gsl_root_fsolver_iterate(s) bind(c,name='gsl_root_fsolver_iterate')
       !!{
       Template for the GSL root solver iterate function.
       !!}
       import
       type(c_ptr), value :: s
     end function gsl_root_fsolver_iterate

     integer(c_int) function gsl_root_fdfsolver_iterate(s) bind(c,name='gsl_root_fdfsolver_iterate')
       !!{
       Template for the GSL root solver iterate function.
       !!}
       import
       type(c_ptr), value :: s
     end function gsl_root_fdfsolver_iterate

     real(c_double) function gsl_root_fsolver_root(s) bind(c,name='gsl_root_fsolver_root')
       !!{
       Template for the GSL root solver root function.
       !!}
       import
       type(c_ptr), value :: s
     end function gsl_root_fsolver_root

     real(c_double) function gsl_root_fdfsolver_root(s) bind(c,name='gsl_root_fdfsolver_root')
       !!{
       Template for the GSL root solver root function.
       !!}
       import
       type(c_ptr), value :: s
     end function gsl_root_fdfsolver_root

     integer(c_int) function gsl_root_test_delta(x1,x0,epsabs,epsrel) bind(c,name='gsl_root_test_delta')
       !!{
       Template for the GSL root solver test delta function.
       !!}
       import
       real(c_double), value :: x1    , x0    , &
            &                   epsabs, epsrel
     end function gsl_root_test_delta

     integer(c_int) function gsl_root_test_interval(x_lower,x_upper,epsabs,epsrel) bind(c,name='gsl_root_test_interval')
       !!{
       Template for the GSL root solver test delta function.
       !!}
       import
       real(c_double), value :: x_lower, x_upper, &
            &                   epsabs , epsrel
     end function gsl_root_test_interval

     real(c_double) function gsl_root_fsolver_x_lower(s) bind(c,name='gsl_root_fsolver_x_lower')
       !!{
       Template for the GSL root solver x-lower function.
       !!}
       import
       type(c_ptr), value :: s
     end function gsl_root_fsolver_x_lower

     real(c_double) function gsl_root_fsolver_x_upper(s) bind(c,name='gsl_root_fsolver_x_upper')
       !!{
       Template for the GSL root solver x-upper function.
       !!}
       import
       type(c_ptr), value :: s
     end function gsl_root_fsolver_x_upper

     function gsl_fsolver_type_get(i) bind(c,name='gsl_fsolver_type_get')
       !!{
       Template for GSL interface {\normalfont \ttfamily fsolver} type function.
       !!}
       import c_ptr, c_int
       type   (c_ptr)                       :: gsl_fsolver_type_get
       integer(c_int), intent(in   ), value :: i
     end function gsl_fsolver_type_get

     function gsl_fdfsolver_type_get(i) bind(c,name='gsl_fdfsolver_type_get')
       !!{
       Template for GSL interface {\normalfont \ttfamily fdfsolver} type function.
       !!}
       import c_ptr, c_int
       type   (c_ptr)                       :: gsl_fdfsolver_type_get
       integer(c_int), intent(in   ), value :: i
     end function gsl_fdfsolver_type_get
  end interface
  
contains
  
  function rootFinderConstructorInternal(rootFunction,rootFunctionDerivative,rootFunctionBoth,solverType,toleranceAbsolute,toleranceRelative,rangeExpandType,rangeExpandUpward,rangeExpandDownward,rangeUpwardLimit,rangeDownwardLimit,rangeExpandUpwardSignExpect,rangeExpandDownwardSignExpect,testLimits,stoppingCriterion) result(self)
    !!{
    Internal constructor for root finders.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (rootFinder                          )                          :: self
    double precision                                      , intent(in   ), optional :: toleranceAbsolute            , toleranceRelative
    integer                                               , intent(in   ), optional :: solverType
    type            (enumerationRangeExpandType          ), intent(in   ), optional :: rangeExpandType
    type            (enumerationRangeExpandSignExpectType), intent(in   ), optional :: rangeExpandDownwardSignExpect, rangeExpandUpwardSignExpect
    type            (enumerationStoppingCriterionType    ), intent(in   ), optional :: stoppingCriterion
    double precision                                      , intent(in   ), optional :: rangeDownwardLimit           , rangeExpandDownward        , &
         &                                                                             rangeExpandUpward            , rangeUpwardLimit
    logical                                               , intent(in   ), optional :: testLimits
    procedure       (rootFunctionTemplate                )               , optional :: rootFunction
    procedure       (rootFunctionDerivativeTemplate      )               , optional :: rootFunctionDerivative
    procedure       (rootFunctionBothTemplate            )               , optional :: rootFunctionBoth
    
    ! Initialize GSL objects to null pointers.
    self%gslFunction                  =c_null_ptr
    self%solver                       =c_null_ptr
    self%solverType                   =c_null_ptr
    ! Initialize to a null solver type.
    self%solverTypeID                 =0
    ! Initialize to tolerances at machine precision.
    self%toleranceAbsolute            =        0.0d0
    self%toleranceRelative            =epsilon(0.0d0)
    ! Initialize state.
    self%initialized                  =.false.
    self%functionInitialized          =.false.
    self%resetRequired                =.false.
    self%useDerivative                =.false.
    ! Initialize range expansion to no expansion.
    self%rangeExpandType              =rangeExpandNull
    self%rangeExpandUpward            =1.0d0
    self%rangeExpandDownward          =1.0d0
    self%rangeUpwardLimitSet          =.false.
    self%rangeDownwardLimitSet        =.false.
    self%rangeExpandDownwardSignExpect=rangeExpandSignExpectNone
    self%rangeExpandUpwardSignExpect  =rangeExpandSignExpectNone
    ! Initialize stopping criterion to an interval test.
    self%stoppingCriterion            =stoppingCriterionInterval
    ! If functions are provided, set them.
    if (present(rootFunction)) then
       if (present(rootFunctionDerivative).or.present(rootFunctionBoth)) then
          if      (.not.present(rootFunctionDerivative)) then
             call Error_Report('missing "rootFunctionDerivative"'//{introspection:location})
          else if (.not.present(rootFunctionBoth      )) then
             call Error_Report('missing "rootFunctionBoth"'      //{introspection:location})
          else
             call self%rootFunctionDerivative(rootFunction,rootFunctionDerivative,rootFunctionBoth)
          end if
       else
          call self%rootFunction(rootFunction)
       end if
    else if (present(rootFunctionDerivative).or.present(rootFunctionBoth)) then
       call Error_Report('missing "rootFunction"'//{introspection:location})
    end if
    ! If a stopping criterion is provided, set it.
    if (present(stoppingCriterion)) self%stoppingCriterion=stoppingCriterion
    ! Validate stopping criterion.
    if (self%useDerivative .and. self%stoppingCriterion == stoppingCriterionInterval) &
         & call Error_Report('"interval" stopping criteria is not valid when using a derivative-based method'//{introspection:location})
    ! If a solver type is provided, set that.
    if (present(solverType)) call self%type(solverType)
    ! If tolerances are provided, set them.
    call self%tolerance(toleranceAbsolute,toleranceRelative)
    ! If range expansion is defined, set it.
    call self%rangeExpand(rangeExpandUpward,rangeExpandDownward,rangeExpandType,rangeUpwardLimit,rangeDownwardLimit,rangeExpandDownwardSignExpect,rangeExpandUpwardSignExpect,testLimits)
    return
  end function rootFinderConstructorInternal
  
  subroutine rootFinderDestroy(self)
    !!{
    Destroy a root finder object.
    !!}
    use :: Interface_GSL, only : gslFunctionDestroy
    implicit none
    class(rootFinder), intent(inout) :: self

    if (self%functionInitialized) then
       if (self%useDerivative) then
          call GSL_Root_FdFSolver_Free(self%solver)
       else
          call GSL_Root_FSolver_Free  (self%solver)
       end if
       call gslFunctionDestroy(self%gslFunction)
       self%functionInitialized=.false.
    end if
    return
  end subroutine rootFinderDestroy

  subroutine rootFinderDestructor(self)
    !!{
    Finalize a root finder object.
    !!}
    implicit none
    type(rootFinder), intent(inout) :: self

    call self%destroy()
    return
  end subroutine rootFinderDestructor

  logical function rootFinderIsInitialized(self)
    !!{
    Return whether a {\normalfont \ttfamily rootFinder} object is initialized.
    !!}
    implicit none
    class(rootFinder), intent(in   ) :: self

    rootFinderIsInitialized=self%initialized
    return
  end function rootFinderIsInitialized

  recursive double precision function rootFinderFindGuess(self,rootGuess,report,status)
    !!{
    Wrapper function for root finder finding that accepts an initial guess for the root. Multiple, generic entry points are used
    here (instead of a single entry point with optional arguments) as it seems to result in faster performance, and this function
    is performance critical.
    !!}
    implicit none
    class           (rootFinder   ), intent(inout), target   :: self
    real            (kind=c_double), intent(in   )           :: rootGuess
    logical                        , intent(in   ), optional :: report
    integer                        , intent(  out), optional :: status
    double precision               , dimension(2)            :: rootRange

    rootRange=rootGuess
    rootFinderFindGuess=self%find(rootRange,report,status)
    return
  end function rootFinderFindGuess
  
  recursive double precision function rootFinderFindRange(self,rootRange,report,status)
    !!{
    Wrapper function for root finder finding that accepts an initial range for the root. Multiple, generic entry points are used
    here (instead of a single entry point with optional arguments) as it seems to result in faster performance, and this function
    is performance critical.
    !!}
    implicit none
    class           (rootFinder   )              , intent(inout), target   :: self
    real            (kind=c_double), dimension(2), intent(in   )           :: rootRange
    logical                                      , intent(in   ), optional :: report
    integer                                      , intent(  out), optional :: status
    double precision               , dimension(2)                          :: rootRangeValues

    if (self%useDerivative) then
       ! Values will not be used, so do not compute them.
       rootRangeValues=-huge(0.0d0)
    else
       rootRangeValues(1)=self%finderFunction(rootRange(1))
       rootRangeValues(2)=self%finderFunction(rootRange(2))
    end if
    rootFinderFindRange=self%find(rootRange,rootRangeValues,report,status)
    return
  end function rootFinderFindRange
  
  recursive double precision function rootFinderFindRangeValueLow(self,rootRange,rootRangeValueLow,report,status)
    !!{
    Wrapper function for root finder finding that accepts an initial range for the root and function value at the lower end of that
    range. Multiple, generic entry points are used here (instead of a single entry point with optional arguments) as it seems to
    result in faster performance, and this function is performance critical.    
    !!}
    implicit none
    class  (rootFinder   )              , intent(inout), target   :: self
    real   (kind=c_double), dimension(2), intent(in   )           :: rootRange
    real   (kind=c_double)              , intent(in   )           :: rootRangeValueLow
    logical                             , intent(in   ), optional :: report
    integer                             , intent(  out), optional :: status
    double precision      , dimension(2)                          :: rootRangeValues

    if (self%useDerivative) then
       ! Values will not be used, so do not compute them.
       rootRangeValues=-huge(0.0d0)
    else
       rootRangeValues(1)=                    rootRangeValueLow
       rootRangeValues(2)=self%finderFunction(rootRange        (2))
    end if
    rootFinderFindRangeValueLow=self%find_(rootRange,rootRangeValues,report,status)
    return
  end function rootFinderFindRangeValueLow

  recursive double precision function rootFinderFindRangeValueHigh(self,rootRange,rootRangeValueHigh,report,status)
    !!{
    Wrapper function for root finder finding that accepts an initial range for the root and function value at the upper end of that
    range. Multiple, generic entry points are used here (instead of a single entry point with optional arguments) as it seems to
    result in faster performance, and this function is performance critical.    
    !!}
    implicit none
    class  (rootFinder   )              , intent(inout), target   :: self
    real   (kind=c_double), dimension(2), intent(in   )           :: rootRange
    real   (kind=c_double)              , intent(in   )           :: rootRangeValueHigh
    logical                             , intent(in   ), optional :: report
    integer                             , intent(  out), optional :: status
    double precision      , dimension(2)                          :: rootRangeValues

    if (self%useDerivative) then
       ! Values will not be used, so do not compute them.
       rootRangeValues=-huge(0.0d0)
    else
       rootRangeValues(1)=self%finderFunction(rootRange         (1))
       rootRangeValues(2)=                    rootRangeValueHigh
    end if
    rootFinderFindRangeValueHigh=self%find_(rootRange,rootRangeValues,report,status)
    return
  end function rootFinderFindRangeValueHigh

  recursive double precision function rootFinderFindRangeValues(self,rootRange,rootRangeValues,report,status)
    !!{
    Wrapper function for root finder finding that accepts an initial range for the root and function values at the ends of that
    range. Multiple, generic entry points are used here (instead of a single entry point with optional arguments) as it seems to
    result in faster performance, and this function is performance critical.    
    !!}
    implicit none
    class  (rootFinder   )              , intent(inout), target   :: self
    real   (kind=c_double), dimension(2), intent(in   )           :: rootRange, rootRangeValues
    logical                             , intent(in   ), optional :: report
    integer                             , intent(  out), optional :: status
    
    rootFinderFindRangeValues=self%find_(rootRange,rootRangeValues,report,status)
    return
  end function rootFinderFindRangeValues

  recursive double precision function rootFinderFind(self,rootRange,rootRangeValues,report,status)
    !!{
    Finds the root of the supplied {\normalfont \ttfamily root} function.
    !!}
    use            :: Display           , only : displayMessage            , verbosityLevelWarn   , displayIndent     , displayUnindent
    use            :: Error             , only : Error_Report              , errorStatusOutOfRange, errorStatusSuccess, GSL_Error_Handler_Abort_Off, &
         &                                       GSL_Error_Handler_Abort_On
    use, intrinsic :: ISO_C_Binding     , only : c_funptr
    use            :: ISO_Varying_String, only : assignment(=)             , operator(//)         , varying_string
    use            :: Interface_GSL     , only : GSL_Success               , gslFunction          , gslFunctionFdF
    implicit none
    class           (rootFinder          )              , intent(inout), target   :: self
    real            (kind=c_double       ), dimension(2), intent(in   )           :: rootRange             , rootRangeValues
    logical                                             , intent(in   ), optional :: report
    integer                                             , intent(  out), optional :: status
    type            (rootFinderList      ), dimension(:), allocatable             :: currentFindersTmp
    integer                               , parameter                             :: iterationMaximum =1000
    integer                               , parameter                             :: findersIncrement =   3
    logical                                                                       :: rangeChanged          , rangeLowerAsExpected, rangeUpperAsExpected
    integer                                                                       :: iteration             , statusActual
    double precision                                                              :: xHigh                 , xLow                , xRoot               , &
         &                                                                           xRootPrevious         , fLow                , fHigh               , &
         &                                                                           fDownwardLimit        , fUpwardLimit
    type            (varying_string      ), save                                  :: message
    !$omp threadprivate(message)
    character       (len= 30             )                                        :: label
    !![
    <optionalArgument name="report" defaultsTo=".false."/>
    !!]

    ! Begin reporting.
    if (report_) call displayIndent("Root finder report:")
    ! Add the current finder to the list of finders. This allows us to track back to the previously used finder if this function is called recursively.
    currentFinderIndex=currentFinderIndex+1
    if (allocated(currentFinders)) then
       if (size(currentFinders) < currentFinderIndex) then
          call move_alloc(currentFinders,currentFindersTmp)
          allocate(currentFinders(size(currentFindersTmp)+findersIncrement))
          currentFinders(1:size(currentFindersTmp))=currentFindersTmp
          deallocate(currentFindersTmp)
       end if
    else
       allocate(currentFinders(findersIncrement))
    end if
    currentFinders(currentFinderIndex)%finder => self
    ! Initialize the root finder variables if necessary.
    if (self%useDerivative) then
       if (.not.self%functionInitialized.or.self%resetRequired) then
          if (     self%functionInitialized  ) call GSL_Root_fdfSolver_Free(self%solver)
          if (.not.self%solverTypeIsValid()) then
             self%solverTypeID    =gsl_root_fdfsolver_steffenson
             self%solverType      =gsl_fdfsolver_type_get  (self%solverTypeID)
          end if
          self%gslFunction        =gslFunctionFdF          (                               &
               &                                            rootFunctionWrapper          , &
               &                                            rootFunctionDerivativeWrapper, &
               &                                            rootFunctionBothWrapper        &    
               &                                           )
          self%solver             =GSL_Root_fdfSolver_Alloc(self%solverType)
          self%resetRequired      =.false.
          self%functionInitialized=.true.
       end if
    else
       if (.not.self%functionInitialized.or.self%resetRequired) then
          if (     self%functionInitialized  ) call GSL_Root_fSolver_Free(self%solver)
          if (.not.self%solverTypeIsValid()) then
             self%solverTypeID    =gsl_root_fsolver_brent
             self%solverType      =gsl_fsolver_type_get  (self%solverTypeID           )
          end if
          self%gslFunction        =gslFunction           (rootFunctionWrapper)
          self%solver             =GSL_Root_fSolver_Alloc(self%solverType             )
          self%resetRequired      =.false.
          self%functionInitialized=.true.
      end if
    end if
    ! Initialize range.
    xLow =rootRange(1)
    xHigh=rootRange(2)
    if (report_) then
       write (label,'(e12.6,a1,e12.6)') xLow,":",xHigh
       message="initial range = "//trim(label)
       call displayMessage(message)
    end if
    ! If we have a downward limit, and know the sign to expect at the downward limit we can check if a solution can be found at
    ! all.
    if (self%testLimits .and. self%rangeDownwardLimitSet .and. self%rangeExpandDownwardSignExpect%ID /= rangeExpandSignExpectNone%ID) then
       fDownwardLimit=self%finderFunction(self%rangeDownwardLimit) 
       select case (self%rangeExpandDownwardSignExpect%ID)
       case (rangeExpandSignExpectNegative%ID)
          rangeLowerAsExpected=(fDownwardLimit  <= 0.0d0)
       case (rangeExpandSignExpectPositive%ID)
          rangeLowerAsExpected=(fDownwardLimit  >= 0.0d0)
       case default
          rangeLowerAsExpected=.false.
          call Error_Report('inconsistent expectation'//{introspection:location})
       end select
       if (.not.rangeLowerAsExpected) then
          if (present(status)) then
             status            =errorStatusOutOfRange
             currentFinderIndex=currentFinderIndex-1
             rootFinderFind    =self%rangeDownwardLimit
             return
          else
             message='root function has incorrect sign at downward limit'
             write (label,'(e12.6,a1,e12.6)') self%rangeDownwardLimit ,":",fDownwardLimit
             message=message//char(10)//'xDownwardLimit :f(xDownwardLimit)='//trim(label)
             call Error_Report(message)
          end if
       else if (report) then
          write (label,'(e12.6,a1,e12.6)') self%rangeDownwardLimit ,":",fDownwardLimit
          message="function at downward limit: x, f(x) = "//trim(label)
          call displayMessage(message)
       end if
    end if
    ! If we have a upward limit, and know the sign to expect at the upward limit we can check if a solution can be found at
    ! all.
    if (self%testLimits .and. self%rangeUpwardLimitSet .and. self%rangeExpandUpwardSignExpect%ID /= rangeExpandSignExpectNone%ID) then
       fUpwardLimit=self%finderFunction(self%rangeUpwardLimit) 
       select case (self%rangeExpandUpwardSignExpect%ID)
       case (rangeExpandSignExpectNegative%ID)
          rangeUpperAsExpected=(fUpwardLimit  <= 0.0d0)
       case (rangeExpandSignExpectPositive%ID)
          rangeUpperAsExpected=(fUpwardLimit  >= 0.0d0)
       case default
          rangeUpperAsExpected=.false.
          call Error_Report('inconsistent expectation'//{introspection:location})
       end select
       if (.not.rangeUpperAsExpected) then
          if (present(status)) then
             status            =errorStatusOutOfRange
             currentFinderIndex=currentFinderIndex-1
             rootFinderFind    =self%rangeUpwardLimit
             return
          else
             message='root function has incorrect sign at upward limit'
             write (label,'(e12.6,a1,e12.6)') self%rangeUpwardLimit ,":",fUpwardLimit
             message=message//char(10)//'xUpwardLimit :f(xUpwardLimit)='//trim(label)
             call Error_Report(message)
          end if
       else if (report) then
          write (label,'(e12.6,a1,e12.6)') self%rangeUpwardLimit ,":",fUpwardLimit
          message="function at upward limit: x, f(x) = "//trim(label)
          call displayMessage(message)
       end if
    end if
    ! Expand the range as necessary.
    if (self%useDerivative) then
       xRoot       =0.5d0*(xLow+xHigh)
       statusActual=GSL_Root_fdfSolver_Set(self%solver,self%gslFunction,xRoot)
    else
       currentFinders(currentFinderIndex)%lowInitialUsed =.true.
       currentFinders(currentFinderIndex)%highInitialUsed=.true.
       fLow=rootRangeValues(1)
       if (xHigh == xLow) then
          ! If a root guess was used, the initial xHigh will equal xLow, so we can avoid re-evaluating the function here.
          fHigh=fLow
       else
          fHigh=rootRangeValues(2)
       end if
       if (report_) call displayIndent("expanding range")
       do while (sign(1.0d0,fLow)*sign(1.0d0,fHigh) > 0.0d0 .and. fLow /= 0.0d0 .and. fHigh /= 0.0d0)
          rangeChanged=.false.
          select case (self%rangeExpandDownwardSignExpect%ID)
          case (rangeExpandSignExpectNegative%ID)
             rangeLowerAsExpected=(fLow  < 0.0d0)
          case (rangeExpandSignExpectPositive%ID)
             rangeLowerAsExpected=(fLow  > 0.0d0)
          case default
             rangeLowerAsExpected=.false.
          end select
          select case (self%rangeExpandUpwardSignExpect  %ID)
          case (rangeExpandSignExpectNegative%ID)
             rangeUpperAsExpected=(fHigh < 0.0d0)
          case (rangeExpandSignExpectPositive%ID)
             rangeUpperAsExpected=(fHigh > 0.0d0)
          case default
             rangeUpperAsExpected=.false.
          end select
          call reportState(includeExpectations=.true.)
          select case (self%rangeExpandType%ID)
          case (rangeExpandAdditive      %ID)
             if     (                                  &
                  &   self%rangeExpandUpward   > 0.0d0 &
                  &  .and.                             &
                  &  .not.rangeUpperAsExpected         &
                  &  .and.                             &
                  &  (                                 &
                  &   xHigh < self%rangeUpwardLimit    &
                  &   .or.                             &
                  &   .not.self%rangeUpwardLimitSet    &
                  &  )                                 &
                  & ) then
                if (rangeLowerAsExpected) then
                   ! The lower end of the range has the expected sign. Therefore, we can shift the lower end of the range to the
                   ! current upper end before updating the upper end.
                   xLow=xHigh
                   fLow=fHigh
                end if
                xHigh=xHigh+self%rangeExpandUpward
                if (self%rangeUpwardLimitSet  ) xHigh=min(xHigh,self%rangeUpwardLimit  )
                fHigh=self%finderFunction(xHigh)
                rangeChanged=.true.
             end if
             if     (                                  &
                  &   self%rangeExpandDownward < 0.0d0 &
                  &  .and.                             &
                  &  .not.rangeLowerAsExpected         &
                  &  .and.                             &
                  &  (                                 &
                  &   xLow  > self%rangeDownwardLimit  &
                  &   .or.                             &
                  &   .not.self%rangeDownwardLimitSet  &
                  &  )                                 &
                  & ) then
                if (rangeUpperAsExpected) then
                   ! The upper end of the range has the expected sign. Therefore, we can shift the upper end of the range to the
                   ! current lower end before updating the lower end.
                   xHigh=xLow
                   fHigh=fLow
                end if
                xLow =xLow +self%rangeExpandDownward
                if (self%rangeDownwardLimitSet) xLow =max(xLow ,self%rangeDownwardLimit)
                fLow =self%finderFunction(xLow )
                rangeChanged=.true.
             end if
          case (rangeExpandMultiplicative%ID)
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
                if (rangeLowerAsExpected) then
                   ! The lower end of the range has the expected sign. Therefore, we can shift the lower end of the range to the
                   ! current upper end before updating the upper end.
                   xLow=xHigh
                   fLow=fHigh
                end if
                xHigh=xHigh*self%rangeExpandUpward
                if (self%rangeUpwardLimitSet  ) xHigh=min(xHigh,self%rangeUpwardLimit  )
                fHigh=self%finderFunction(xHigh)
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
                if (rangeUpperAsExpected) then
                   ! The upper end of the range has the expected sign. Therefore, we can shift the upper end of the range to the
                   ! current lower end before updating the lower end.
                   xHigh=xLow
                   fHigh=fLow
                end if
                xLow =xLow *self%rangeExpandDownward
                if (self%rangeDownwardLimitSet) xLow =max(xLow ,self%rangeDownwardLimit)
                fLow =self%finderFunction(xLow )
                rangeChanged=.true.
             end if
          end select
          if (.not.rangeChanged) then
             rootFinderFind=0.0d0
             if (report_) then
                call displayMessage ('failed to bracket root')
                call displayUnindent('done'                  )
             end if
             if (present(status)) then
                status=errorStatusOutOfRange
                currentFinderIndex=currentFinderIndex-1
                return
             else
                fLow =self%finderFunction(xLow )
                fHigh=self%finderFunction(xHigh)
                message='unable to expand range to bracket root'
                write (label,'(e12.6,a1,e12.6)') xLow ,":",fLow
                message=message//char(10)//'xLow :f(xLow )='//trim(label)
                write (label,'(e12.6,a1,e12.6)') xHigh,":",fHigh
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
                call Error_Report(message//{introspection:location})
             end if
          end if
       end do
       call reportState(includeExpectations=.false.)
       if (report_) call displayUnindent("done")
       ! Store the values of the function at the lower and upper extremes of the range.
       currentFinders(currentFinderIndex)%xLowInitial    =xLow
       currentFinders(currentFinderIndex)%xHighInitial   =xHigh
       currentFinders(currentFinderIndex)%fLowInitial    =fLow
       currentFinders(currentFinderIndex)%fHighInitial   =fHigh
       currentFinders(currentFinderIndex)%lowInitialUsed =.false.
       currentFinders(currentFinderIndex)%highInitialUsed=.false.
       ! Set the initial range for the solver.
       statusActual=GSL_Root_fSolver_Set(self%solver,self%gslFunction,xLow,xHigh)
    end if
    ! Set error handler if necessary.
    if (present(status)) then
       call GSL_Error_Handler_Abort_Off()
       statusActual=errorStatusSuccess
    end if
    ! Find the root.
    if (statusActual /= GSL_Success) then
       rootFinderFind=0.0d0
       if (present(status)) then
          status=statusActual
       else
          call Error_Report('failed to initialize solver'//{introspection:location})
       end if
    else
       iteration=0
       xRoot    =0.0d0
       if (report_) call displayIndent("seeking root")
       do
          iteration=iteration+1
          if (self%useDerivative) then
             statusActual=GSL_Root_fdfSolver_Iterate(self%solver)
          else
             statusActual=GSL_Root_fSolver_Iterate  (self%solver)
          end if
          if (statusActual /= GSL_Success .or. iteration > iterationMaximum) exit
          if (iteration > 1) then
             select case (self%stoppingCriterion%ID)
             case (stoppingCriterionDelta   %ID)
                xRootPrevious=xRoot
                xRoot        =GSL_Root_fdfSolver_Root(self%solver)
                if (report_) then
                   write (label,'(e12.6,a2,e12.6)') xRoot,", ",self%finderFunction(xRoot)
                   message="xRoot, fRoot  = "//trim(label)
                   call displayMessage(message)
                end if
                statusActual =GSL_Root_Test_Delta(xRoot,xRootPrevious,self%toleranceAbsolute,self%toleranceRelative)
             case (stoppingCriterionInterval%ID)
                xRoot =GSL_Root_fSolver_Root   (self%solver)
                xLow  =GSL_Root_fSolver_x_Lower(self%solver)
                xHigh =GSL_Root_fSolver_x_Upper(self%solver)
                if (report_) then
                   write (label,'(e12.6,a2,e12.6)') xRoot,", ",self%finderFunction(xRoot)
                   message="xRoot, fRoot  = "//trim(label)
                   call displayMessage(message)
                   fLow =self%finderFunction(xLow )
                   fHigh=self%finderFunction(xHigh)
                   call reportState(includeExpectations=.false.)
                end if
                statusActual=GSL_Root_Test_Interval(xLow,xHigh,self%toleranceAbsolute,self%toleranceRelative)
             case default
                call Error_Report('unknown stopping criterion'//{introspection:location})
             end select
             if (statusActual == GSL_Success) exit
          end if
       end do
       if (statusActual /= GSL_Success) then
          rootFinderFind=0.0d0
          if (present(status)) then
             status=statusActual
             call displayMessage("failed to find root")
          else
             if (report_) call displayUnindent("done")
             call Error_Report('failed to find root'//{introspection:location})
          end if
       else
          if (present(status)) status=GSL_Success
          rootFinderFind=xRoot
       end if
       if (report_) call displayUnindent("done")
    end if
    ! Reset error handler.
    if (present(status)) call GSL_Error_Handler_Abort_On()
    ! Restore state.
    currentFinderIndex=currentFinderIndex-1
    ! Finish reporting.
    if (report_) call displayUnindent("done")
    return

  contains

    subroutine reportState(includeExpectations)
      !!{
      Report on the state of the root finder.
      !!}
      implicit none
      logical, intent(in   ) :: includeExpectations
      
      if (.not.report_) return
      write (label,'(e12.6,a2,e12.6)') xLow,", ",fLow
      message="xLow , fLow  = "//trim(label)
      if (includeExpectations .and. self%rangeExpandDownwardSignExpect /= rangeExpandSignExpectNone) then
         if (rangeLowerAsExpected) then
            message=message//" [sign correct]"
         else
            message=message//" [sign incorrect]"
         end if
      end if
      if (self%rangeDownwardLimitSet) then
         if (xLow <= self%rangeDownwardLimit) then
            message=message//" {at limit}"
         else
            message=message//" {not at limit}"
         end if
      end if
      call displayMessage(message)      
      write (label,'(e12.6,a2,e12.6)') xHigh,", ",fHigh
      message="xHigh, fHigh = "//trim(label)
      if (includeExpectations .and. self%rangeExpandUpwardSignExpect /= rangeExpandSignExpectNone) then
         if (rangeUpperAsExpected) then
            message=message//" [sign correct]"
         else
            message=message//" [sign incorrect]"
         end if
      end if
      if (self%rangeUpwardLimitSet) then
         if (xHigh >= self%rangeUpwardLimit) then
            message=message//" {at limit}"
         else
            message=message//" {not at limit}"
         end if
      end if
      call displayMessage(message)      
      return
    end subroutine reportState
    
  end function rootFinderFind

  subroutine rootFinderRootFunction(self,rootFunction)
    !!{
    Sets the function to use in a {\normalfont \ttfamily rootFinder} object.
    !!}
    implicit none
    class    (rootFinder          ), intent(inout) :: self
    procedure(rootFunctionTemplate)                :: rootFunction

    call self%destroy()
    self%finderFunction => rootFunction
    self%initialized    =  .true.
    self%useDerivative  =  .false.
    self%resetRequired  =  .true.
    return
  end subroutine rootFinderRootFunction

  subroutine rootFinderRootFunctionDerivative(self,rootFunction,rootFunctionDerivative,rootFunctionBoth)
    !!{
    Sets the function to use in a {\normalfont \ttfamily rootFinder} object.
    !!}
    implicit none
    class    (rootFinder                    ), intent(inout) :: self
    procedure(rootFunctionTemplate          )                :: rootFunction
    procedure(rootFunctionDerivativeTemplate)                :: rootFunctionDerivative
    procedure(rootFunctionBothTemplate      )                :: rootFunctionBoth

    call self%destroy()
    self%finderFunction           => rootFunction
    self%finderFunctionDerivative => rootFunctionDerivative
    self%finderFunctionBoth       => rootFunctionBoth
    self%initialized              =  .true.
    self%useDerivative            =  .true.
    self%resetRequired            =  .true.
    return
  end subroutine rootFinderRootFunctionDerivative

  subroutine rootFinderType(self,solverType)
    !!{
    Sets the type to use in a {\normalfont \ttfamily rootFinder} object.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (rootFinder), intent(inout) :: self
    integer            , intent(in   ) :: solverType

    ! Set the solver type and indicate that a reset will be required to update the internal GSL objects.
    self%solverTypeID =solverType
    select case (solverType)
    case (gsl_root_fsolver_bisection,gsl_root_fsolver_brent   ,gsl_root_fsolver_falsepos    )
       self%solverType=gsl_fsolver_type_get  (self%solverTypeID)
    case (gsl_root_fdfsolver_newton ,gsl_root_fdfsolver_secant,gsl_root_fdfsolver_steffenson)
       self%solverType=gsl_fdfsolver_type_get(self%solverTypeID)
    case default
       call Error_Report('unknown solver type'//{introspection:location})
    end select
    self%resetRequired=.true.
    if (.not.self%solverTypeIsValid()) call Error_Report('invalid solver type'//{introspection:location})
    return
  end subroutine rootFinderType

  subroutine rootFinderTolerance(self,toleranceAbsolute,toleranceRelative)
    !!{
    Sets the tolerances to use in a {\normalfont \ttfamily rootFinder} object.
    !!}
    implicit none
    class           (rootFinder), intent(inout)           :: self
    double precision            , intent(in   ), optional :: toleranceAbsolute, toleranceRelative

    if (present(toleranceAbsolute)) self%toleranceAbsolute=toleranceAbsolute
    if (present(toleranceRelative)) self%toleranceRelative=toleranceRelative
    return
  end subroutine rootFinderTolerance

  subroutine rootFinderRangeExpand(self,rangeExpandUpward,rangeExpandDownward,rangeExpandType,rangeUpwardLimit,rangeDownwardLimit,rangeExpandDownwardSignExpect,rangeExpandUpwardSignExpect,testLimits)
    !!{
    Sets the rules for range expansion to use in a {\normalfont \ttfamily rootFinder} object.
    !!}
    implicit none
    class           (rootFinder                          ), intent(inout)           :: self
    type            (enumerationRangeExpandType          ), intent(in   ), optional :: rangeExpandType
    type            (enumerationRangeExpandSignExpectType), intent(in   ), optional :: rangeExpandDownwardSignExpect, rangeExpandUpwardSignExpect
    double precision                                      , intent(in   ), optional :: rangeDownwardLimit           , rangeExpandDownward        , &
         &                                                                             rangeExpandUpward            , rangeUpwardLimit
    logical                                               , intent(in   ), optional :: testLimits

    if (present(rangeExpandUpward            )) self%rangeExpandUpward  =rangeExpandUpward
    if (present(rangeExpandDownward          )) self%rangeExpandDownward=rangeExpandDownward
    if (present(rangeExpandType              )) then
       self%rangeExpandType    =rangeExpandType
    else
       self%rangeExpandType    =rangeExpandNull
    end if
    select case (self%rangeExpandType%ID)
    case (rangeExpandAdditive      %ID)
       if (.not.present(rangeExpandUpward  )) self%rangeExpandUpward  =0.0d0
       if (.not.present(rangeExpandDownward)) self%rangeExpandDownward=0.0d0
    case (rangeExpandMultiplicative%ID)
       if (.not.present(rangeExpandUpward  )) self%rangeExpandUpward  =1.0d0
       if (.not.present(rangeExpandDownward)) self%rangeExpandDownward=1.0d0
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
    if (present(testLimits)) then
       self%testLimits=testLimits
    else
       self%testLimits=.false.
    end if
    return
  end subroutine rootFinderRangeExpand
  
  logical function rootFinderSolverTypeIsValid(self)
    !!{
    Sets the tolerances to use in a {\normalfont \ttfamily rootFinder} object.
    !!}
    implicit none
    class(rootFinder), intent(inout) :: self

    rootFinderSolverTypeIsValid=self%solverTypeID /= 0
    if (.not.rootFinderSolverTypeIsValid) return
    if (self%useDerivative) then
       rootFinderSolverTypeIsValid= self%solverTypeID == gsl_root_fdfsolver_newton     &
            &                      .or.                                                &
            &                       self%solverTypeID == gsl_root_fdfsolver_secant     &
            &                      .or.                                                &
            &                       self%solverTypeID == gsl_root_fdfsolver_steffenson
    else
       rootFinderSolverTypeIsValid= self%solverTypeID == gsl_root_fsolver_bisection    &
            &                      .or.                                                &
            &                       self%solverTypeID == gsl_root_fsolver_brent        &
            &                      .or.                                                &
            &                       self%solverTypeID == gsl_root_fsolver_falsepos    
    end if
    return
  end function rootFinderSolverTypeIsValid

  recursive function rootFunctionWrapper(x) bind(c)
    !!{
    Wrapper function callable by {\normalfont \ttfamily GSL} used in root finding.
    !!}
    implicit none
    real(c_double), intent(in   ), value :: x
    real(c_double)                       :: rootFunctionWrapper

    ! Attempt to use previously computed solutions if possible.
    if      (.not.currentFinders(currentFinderIndex)%lowInitialUsed  .and. x == currentFinders(currentFinderIndex)%xLowInitial ) then
       rootFunctionWrapper=currentFinders(currentFinderIndex)%fLowInitial
       currentFinders(currentFinderIndex)%lowInitialUsed =.true.
    else if (.not.currentFinders(currentFinderIndex)%highInitialUsed .and. x == currentFinders(currentFinderIndex)%xHighInitial) then
       rootFunctionWrapper=currentFinders(currentFinderIndex)%fHighInitial
       currentFinders(currentFinderIndex)%highInitialUsed=.true.
    else
       ! No previously computed solution available - evaluate the function.
       rootFunctionWrapper=currentFinders(currentFinderIndex)%finder%finderFunction(x)
    end if
    return
  end function rootFunctionWrapper

  recursive function rootFunctionDerivativeWrapper(x) bind(c)
    !!{
    Wrapper function callable by {\normalfont \ttfamily GSL} used in root finding.
    !!}
    implicit none
    real(c_double)                       :: rootFunctionDerivativeWrapper
    real(c_double), intent(in   ), value :: x

    rootFunctionDerivativeWrapper=currentFinders(currentFinderIndex)%finder%finderFunctionDerivative(x)
    return
  end function rootFunctionDerivativeWrapper

  recursive subroutine rootFunctionBothWrapper(x,parameters,f,df) bind(c)
    !!{
    Wrapper function callable by {\normalfont \ttfamily GSL} used in root finding.
    !!}
    implicit none
    real(c_double), intent(in   ), value :: x
    real(c_double), intent(  out)        :: f         , df
    type(c_ptr   ), intent(in   ), value :: parameters
    !$GLC attributes unused :: parameters

    call currentFinders(currentFinderIndex)%finder%finderFunctionBoth(x,f,df)
    return
  end subroutine rootFunctionBothWrapper

end module Root_Finder

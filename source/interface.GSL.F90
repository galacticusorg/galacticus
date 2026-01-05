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
Contains a module which interfaces with low-level aspects of the GSL library.
!!}

! Specify an explicit dependence on the interface.GSL.C.o object file.
!: $(BUILDPATH)/interface.GSL.C.o

! Add dependency on GSL library.
!; gsl

module Interface_GSL
  !!{
  Interfaces with low-level aspects of the GSL library.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_funptr, c_ptr   , c_double, c_int, &
       &                                  c_char  , c_size_t
  private
  public :: gslFunction        , gslFunctionFdF        , gslFunctionDestroy, &
       &    gslFunctionTemplate, gslFunctionFdFTemplate, gslSetErrorHandler, &
       &    gslFileOpen        , gslFileClose          , gslErrorDecode

  abstract interface
     !!{
     Interface for {\normalfont \ttfamily gslFunction} type. We ignore the {\normalfont \ttfamily parameters} argument here as
     it is not used by \glc.
     !!}
     double precision function gslFunctionTemplate(x)
       import c_ptr
       double precision, intent(in   ), value :: x
     end function gslFunctionTemplate
  end interface

  abstract interface
     !!{
     Interface for {\normalfont \ttfamily gslFunctionFdF} type.
     !!}
     subroutine gslFunctionFdFTemplate(x,parameters,f,df)
       import c_ptr
       double precision       , intent(in   ), value :: x
       type            (c_ptr), intent(in   ), value :: parameters
       double precision       , intent(  out)        :: f, df
     end subroutine gslFunctionFdFTemplate
  end interface

  abstract interface
     !!{
     Interface for {\normalfont \ttfamily  gsl\_error\_handler\_t} type.
     !!}
     subroutine gslErrorHandlerTemplate(reason,file,line,errorNumber)
       import c_char, c_int
       character(c_char), dimension(*) :: file       , reason
       integer  (c_int ), value        :: errorNumber, line
     end subroutine gslErrorHandlerTemplate
  end interface

  interface
     !!{
     Interfaces to C functions.
     !!}
     function gslFunctionConstructor(f) bind(c,name="gslFunctionConstructor")
       !!{
       Interface to a C function which establishes a {\normalfont \ttfamily gslFunction} type.
       !!}
       import c_ptr, c_funptr
       type(c_ptr   )        :: gslFunctionConstructor
       type(c_funptr), value :: f
     end function gslFunctionConstructor

     function gslFunctionFdFConstructor(f,df,fdf) bind(c,name="gslFunctionFdFConstructor")
       !!{
       Interface to a C function which establishes a {\normalfont \ttfamily gslFunctionFdF} type.
       !!}
       import c_ptr, c_funptr
       type(c_ptr   )        :: gslFunctionFdFConstructor
       type(c_funptr), value :: f                        , df, &
            &                   fdf
     end function gslFunctionFdFConstructor

     subroutine gslFunctionDestructor(f) bind(c,name="gslFunctionDestructor")
       !!{
       Interface to a C function which destroys a {\normalfont \ttfamily gslFunction} type.
       !!}
       import c_funptr
       type(c_funptr), value :: f
     end subroutine gslFunctionDestructor

     function gsl_set_error_handler(new_handler) bind(c,name="gsl_set_error_handler")
       import c_funptr
       type(c_funptr)        :: gsl_set_error_handler
       type(c_funptr), value :: new_handler
     end function gsl_set_error_handler
     
     function gslFileOpenC(fileName,access) bind(c,name='gslFileOpenC')
       !!{
       Template for a C function that opens a file for GSL state output.
       !!}
       import c_ptr, c_char
       type     (c_ptr )               :: gslFileOpenC
       character(c_char), dimension(*) :: fileName    , access
     end function gslFileOpenC
     
     function gslFileCloseC(stream) bind(c,name='gslFileCloseC')
       !!{
       Template for a C function that opens a file for GSL state output.
       !!}
       import c_ptr, c_int
       integer(c_int)        :: gslFileCloseC
       type   (c_ptr), value :: stream
     end function gslFileCloseC
 
     subroutine gslErrorDecoder(gsl_errno,gsl_str,gsl_strlen) bind(c,name='gslErrorDecoder')
       !!{
       Template for a C function that returns a string describing a GSL error.
       !!}
       import c_int, c_char, c_size_t
       integer  (c_int   ), value        :: gsl_errno
       character(c_char  ), dimension(*) :: gsl_str
       integer  (c_size_t), value        :: gsl_strlen
     end subroutine gslErrorDecoder
  end interface

  interface gslSetErrorHandler
     !!{
     Generic interface to GSL error handler set functions.
     !!}
     module procedure gslSetErrorHandlerFunction
     module procedure gslSetErrorHandlerPointer
  end interface gslSetErrorHandler
  
  type, public, bind(c) :: gsl_sf_result
     !!{
     Type for GSL special function results.
     !!}
     real(c_double) :: val, err
  end type gsl_sf_result

  ! Error codes.
  !![
  <constant variable="GSL_Success"  gslSymbol="GSL_SUCCESS"  gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for success." group="GSL"/>
  <constant variable="GSL_Failure"  gslSymbol="GSL_FAILURE"  gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for generic failure." group="GSL"/>
  <constant variable="GSL_ESing"    gslSymbol="GSL_ESING"    gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for apparent singularity detected." group="GSL"/>
  <constant variable="GSL_EDom"     gslSymbol="GSL_EDOM"     gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for input domain error, e.g sqrt(-1)." group="GSL"/>
  <constant variable="GSL_ERange"   gslSymbol="GSL_ERANGE"   gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for output range error, e.g. exp(1e100)." group="GSL"/>
  <constant variable="GSL_EZeroDiv" gslSymbol="GSL_EZERODIV" gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for attempt to divide by zero." group="GSL"/>
  <constant variable="GSL_EUndrFlw" gslSymbol="GSL_EUNDRFLW" gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for underflow." group="GSL"/>
  <constant variable="GSL_ENoProg"  gslSymbol="GSL_ENOPROG"  gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for iteration not making progress towards solution." group="GSL"/>
  <constant variable="GSL_Continue" gslSymbol="GSL_CONTINUE" gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for iteration has not converged." group="GSL"/>
  <constant variable="GSL_EBadFunc" gslSymbol="GSL_EBADFUNC" gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for problem with user-supplied function." group="GSL"/>
  <constant variable="GSL_EBadTol"  gslSymbol="GSL_EBADTOL"  gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for user-specified invalid tolerance." group="GSL"/>
  <constant variable="GSL_ETol"     gslSymbol="GSL_ETOL"     gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for failure to reach the specified tolerance." group="GSL"/>
  <constant variable="GSL_ERound"   gslSymbol="GSL_EROUND"   gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for failure because of roundoff error." group="GSL"/>
  <constant variable="GSL_EMaxIter" gslSymbol="GSL_EMAXITER" gslHeader="gsl_errno" type="integer" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/err.html#error-codes" description="Error code for exceeding the maximum number of iterations." group="GSL"/>
  !!]

  ! Precision modes.
  !![
  <constant variable="GSL_Prec_Double"  gslSymbol="GSL_PREC_DOUBLE"  gslHeader="gsl_mode" type="integer" description="Specifies GSL double-precision mode." reference="Gnu Scientific Library" group="GSL"/>
  <constant variable="GSL_Prec_Single"  gslSymbol="GSL_PREC_SINGLE"  gslHeader="gsl_mode" type="integer" description="Specifies GSL single-precision mode." reference="Gnu Scientific Library" group="GSL"/>
  <constant variable="GSL_Prec_Approx"  gslSymbol="GSL_PREC_APPROX"  gslHeader="gsl_mode" type="integer" description="Specifies GSL approximate-precision mode." reference="Gnu Scientific Library" group="GSL"/>
  !!]
  
contains

  function gslFileOpen(fileName,access) result(stream)
    !!{
    Open a file for output of GSL state.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_null_char
    implicit none
    type     (c_ptr)                :: stream
    character(len=*), intent(in   ) :: fileName   , access

    stream=gslFileOpenC(trim(fileName)//c_null_char,trim(access)//c_null_char)
    return
  end function gslFileOpen

  subroutine gslFileClose(stream)
    !!{
    Close a file used for output of GSL state.
    !!}
    implicit none
    type   (c_ptr), intent(in   ) :: stream
    integer(c_int)                :: status

    status=gslFileCloseC(stream)
    return
  end subroutine gslFileClose
  
  function gslSetErrorHandlerFunction(newHandler)
    !!{
    Set the GSL error handler to the provided function.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_funloc
    implicit none
    type     (c_funptr               ) :: gslSetErrorHandlerFunction
    procedure(gslErrorHandlerTemplate) :: newHandler

    gslSetErrorHandlerFunction=gsl_set_error_handler(c_funloc(newHandler))
    return
  end function gslSetErrorHandlerFunction

  function gslSetErrorHandlerPointer(handler)
    !!{
    Set the GSL error handler to the provided pointer.
    !!}
    implicit none
    type(c_funptr) :: gslSetErrorHandlerPointer
    type(c_funptr) :: handler

    gslSetErrorHandlerPointer=gsl_set_error_handler(handler)
    return
  end function gslSetErrorHandlerPointer

  function gslFunction(f)
    !!{
    Return a {\normalfont \ttfamily c\_ptr} object for the given function {\normalfont \ttfamily f}.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_funloc
    implicit none
    type     (c_ptr              ) :: gslFunction
    procedure(gslFunctionTemplate) :: f

    gslFunction=gslFunctionConstructor(c_funloc(f))
    return
  end function gslFunction

  function gslFunctionFdF(f,df,fdf)
    !!{
    Return a {\normalfont \ttfamily c\_ptr} object for the given function {\normalfont \ttfamily f}.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_funloc
    implicit none
    type     (c_ptr                 ) :: gslFunctionFdF
    procedure(gslFunctionTemplate   ) :: f             , df
    procedure(gslFunctionFdFTemplate) :: fdf

    gslFunctionFdF=gslFunctionFdFConstructor(c_funloc(f),c_funloc(df),c_funloc(fdf))
    return
  end function gslFunctionFdF

  subroutine gslFunctionDestroy(f)
    !!{
    Destroy a {\normalfont \ttfamily c\_ptr} to a {\normalfont \ttfamily gslFunction} object.
    !!}
    implicit none
    type(c_ptr), intent(in   ) :: f

    call gslFunctionDestructor(f)
    return
  end subroutine gslFunctionDestroy

  function gslErrorDecode(errorNumber) result(description)
    !!{
    Decode a GSL error number.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_f_pointer
    use            :: ISO_Varying_String, only : varying_string     , assignment(=)
    use            :: String_Handling   , only : String_C_to_Fortran
    implicit none
    type     (varying_string)                               :: description
    integer                  , intent(in   )                :: errorNumber
    integer  (c_size_t      ), parameter                    :: descriptionLength=128

    character(c_char        ), dimension(descriptionLength) :: description_
    
    call gslErrorDecoder(errorNumber,description_,descriptionLength)
    description=String_C_to_Fortran(description_)
    return
  end function gslErrorDecode

end module Interface_GSL

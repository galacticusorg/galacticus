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
Contains a module which implements beta functions.
!!}

! Add dependency on GSL library.
!; gsl

module Beta_Functions
  !!{
  Implements beta functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double     , c_int
  use            :: Interface_GSL, only : gsl_sf_result
  implicit none
  private
  public :: Beta_Function, Beta_Function_Incomplete_Normalized

  interface
     function gsl_sf_beta_e(a,b,r) bind(c,name='gsl_sf_beta_e')
       !!{
       Template for the GSL beta C function.
       !!}
       import
       integer(c_int        )        :: gsl_sf_beta_e
       real   (c_double     ), value :: a            , b
       type   (gsl_sf_result)        :: r
     end function gsl_sf_beta_e
  end interface

  interface
     function gsl_sf_beta_inc_e(a,b,x,r) bind(c,name='gsl_sf_beta_inc_e')
       !!{
       Template for the GSL incomplete beta C function.
       !!}
       import
       integer(c_int        )        :: gsl_sf_beta_inc_e
       real   (c_double     ), value :: a                , b, &
            &                           x
       type   (gsl_sf_result)        :: r
     end function gsl_sf_beta_inc_e
  end interface

contains

  double precision function Beta_Function(a,b,status) result(r)
    !!{
    Evaluate the beta function, $B(a,b)$.
    !!}
    use :: Error             , only : Error_Report, GSL_Error_Handler_Abort_Off, GSL_Error_Handler_Abort_On
    use :: Interface_GSL     , only : GSL_Success , gslErrorDecode
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    double precision               , intent(in   )           :: a      , b
    integer                        , intent(  out), optional :: status
    integer         (c_int        )                          :: status_
    type            (gsl_sf_result)                          :: r_

    call GSL_Error_Handler_Abort_Off()
    status_=GSL_SF_Beta_E(a,b,r_    )
    r      =                  r_%val
    if (present(status)) status=status_
    select case (status_)
    case (GSL_Success )
       ! Success - nothing to do.
    case default
       ! Other error - abort or return status code.
       if (.not.present(status)) then
          block
            character(len=44) :: label
            write (label,'(e21.15,", ",e21.15)') a,b
            call Error_Report('incomplete beta function B(a,b) evaluation failed for (a,b)='//label//': '//gslErrorDecode(status_)//{introspection:location})
          end block
       end if
    end select
    call GSL_Error_Handler_Abort_On ()
    return
  end function Beta_Function

  double precision function Beta_Function_Incomplete_Normalized(a,b,x,status) result(r)
    !!{
    Evaluate the normalized incomplete beta function, $B_x(a,b)/B(a,b)$.
    !!}
    use :: Error             , only : Error_Report, GSL_Error_Handler_Abort_Off, GSL_Error_Handler_Abort_On
    use :: Interface_GSL     , only : GSL_Success , GSL_EUndrFlw               , gslErrorDecode
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    double precision               , intent(in   )           :: a                , b     , &
         &                                                      x
    integer                        , intent(  out), optional :: status
    double precision               , parameter               :: countWidths=5.0d0
    integer         (c_int        )                          :: status_
    type            (gsl_sf_result)                          :: r_
    logical                                                  :: reportError
    double precision                                         :: b50              , bWidth
    character       (len=67       )                          :: label

    call GSL_Error_Handler_Abort_Off()
    status_=GSL_SF_Beta_Inc_E(a,b,x,r_    )
    r      =                        r_%val
    reportError=.false.
    if (present(status)) status=status_
    select case (status_)
    case (GSL_Success )
       ! Success - nothing to do.
    case (GSL_EUndrFlw)
       ! Underflow. The normalized, incomplete beta function rises from 0 to 1 as b goes from 0 to ∞. It reaches a value of 0.5 at
       ! a(1-x)/x, and the increase from 0 to 1 occurs over a scale of around √[a(1-x)/x²] around this. Using these facts, if an
       ! underflow error has occured (typically when a and b are large) we can return the limiting behavior if we are sufficiently
       ! far below or above the midpoint.
       !! Compute the midpoint and scale over which the function increases.
       b50   =     a*(1.0d0-x)/x
       bWidth=sqrt(a*(1.0d0-x)/x**2)
       ! Test if we are in either limit.
       if      (b > b50+countWidths*bWidth) then
          r=1.0d0
       else if (b < b50-countWidths*bWidth) then
          r=0.0d0
       else
          ! We are not in a limiting regime, report this error.
          reportError=.true.
       end if
    case default
       ! Other error - abort.
       reportError=.true.
    end select
    if (reportError.and..not.present(status)) then
       write (label,'(e21.15,", ",e21.15,", ",e21.15)') x,a,b
       call Error_Report('incomplete beta function Bₓ(a,b) evaluation failed for (x,a,b)='//label//': '//gslErrorDecode(status_)//{introspection:location})
    end if
    call GSL_Error_Handler_Abort_On ()
    return
  end function Beta_Function_Incomplete_Normalized
  
end module Beta_Functions
